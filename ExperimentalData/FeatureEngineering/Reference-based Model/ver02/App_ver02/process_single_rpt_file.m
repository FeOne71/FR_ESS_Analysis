function ch_data = process_single_rpt_file(filepath, static_step, ocv_steps, win_size, crate_steps)
    % PROCESS_SINGLE_RPT_FILE Parses Lab RPT CSV and extracts V-Q data
    %
    % Inputs:
    %   filepath: path to CSV file
    %   static_step: step index for Static Capacity (default: 3)
    %   ocv_steps: [charge, discharge] step indices for OCV (default: [8, 10])
    %   win_size: moving average window (default: 20)
    %   crate_steps: struct with C-rate step info (optional)
    %     .labels = {'c01','c05','c1','c2','c3'}
    %     .charge = [28, 32, 36, 40, 44]
    %     .discharge = [48, 52, 56, 60, 64]
    %
    % Output:
    %   ch_data: struct with fields:
    %     Static, OCV_charge, OCV_discharge,
    %     c01_charge, c01_discharge, c05_charge, c05_discharge, ...
    
    if nargin < 2 || isempty(static_step), static_step = 3; end
    if nargin < 3 || isempty(ocv_steps), ocv_steps = [8, 10]; end
    if nargin < 4 || isempty(win_size), win_size = 20; end
    if nargin < 5 || isempty(crate_steps)
        crate_steps.labels    = {'c01','c05','c1','c2','c3'};
        crate_steps.charge    = [28, 32, 36, 40, 44];
        crate_steps.discharge = [48, 52, 56, 60, 64];
    end
    
    ch_data = struct();
    try
        % Ultra-fast textscan method for specific Lab RPT format
        fid = fopen(filepath, 'r');
        fgetl(fid);
        formatSpec = '%*f %f %*q %f %s %*s %f %f %f %*[^\n]';
        dataArray = textscan(fid, formatSpec, 'Delimiter', ',', 'EmptyValue', NaN, 'ReturnOnError', false);
        fclose(fid);
        
        steps = dataArray{1};
        subs = dataArray{2}; 
        t_strs = dataArray{3};
        I = dataArray{4};
        V = dataArray{5};
        Q = dataArray{6};
        
        % Convert time strings to seconds (vectorized char matrix method)
        try
            t = zeros(length(t_strs), 1);
            if ~isempty(t_strs)
                if length(t_strs{1}) >= 8
                    strMat = char(t_strs);
                    if size(strMat, 2) >= 8
                        H = (strMat(:,1)-'0')*10 + (strMat(:,2)-'0');
                        M = (strMat(:,4)-'0')*10 + (strMat(:,5)-'0');
                        S = (strMat(:,7)-'0')*10 + (strMat(:,8)-'0');
                        t = H*3600 + M*60 + S;
                    else
                        t = (1:length(V))';
                    end
                else
                    t = (1:length(V))';
                end
            end
        catch
            t = (1:length(V))';
        end
        
    catch ME
        fprintf('    [Error] Failed to read %s (%s)\n', filepath, ME.message);
        if exist('fid','var') && fid > 0, fclose(fid); end
        return;
    end
    
    if isempty(V)
        return;
    end
    
    % 1. Static Capacity
    mask_static = (steps == static_step);
    ch_data.Static = extract_and_interpolate_arrays(V, Q, t, mask_static, win_size);
    
    % 2. OCV Charge
    mask_chg = (steps == ocv_steps(1));
    ch_data.OCV_charge = extract_and_interpolate_arrays(V, Q, t, mask_chg, win_size);
    
    % 3. OCV Discharge (Requires flipping for monotonic V)
    mask_dch = (steps == ocv_steps(2));
    ch_data.OCV_discharge = extract_and_interpolate_flip_arrays(V, Q, t, mask_dch, win_size);
    
    % 4. C-rate data (charge/discharge for each C-rate)
    for r = 1:length(crate_steps.labels)
        label = crate_steps.labels{r};
        
        % Charge
        mask_cr_chg = (steps == crate_steps.charge(r));
        ch_data.([label '_charge']) = extract_and_interpolate_arrays(V, Q, t, mask_cr_chg, win_size);
        
        % Discharge (flip for monotonic V)
        mask_cr_dch = (steps == crate_steps.discharge(r));
        ch_data.([label '_discharge']) = extract_and_interpolate_flip_arrays(V, Q, t, mask_cr_dch, win_size);
    end
end

function seg = extract_and_interpolate_arrays(V_all, Q_all, t_all, mask, win_size)
    seg = struct();
    if sum(mask) < 2, return; end
    
    V = V_all(mask); 
    Q = Q_all(mask);
    t = t_all(mask);
    
    if numel(V) >= win_size
        V = movmean(V, win_size);
        Q = movmean(Q, win_size);
        t = movmean(t, win_size); 
    end
    
    dV = 0.001;
    [V_sorted, si] = sort(V);
    Q_sorted = Q(si);
    t_sorted = t(si);
    
    [V_unique, ~, ic] = unique(V_sorted);
    Q_unique = accumarray(ic, Q_sorted, [], @mean);
    if numel(V_unique) >= 2
        V_grid = (min(V_unique) : dV : max(V_unique))';
        seg.V_grid = V_grid;
        seg.Q = interp1(V_unique, Q_unique, V_grid, 'linear');
        
        t_unique = accumarray(ic, t_sorted, [], @mean);
        seg.t = interp1(V_unique, t_unique, V_grid, 'linear');
    end
end

function seg = extract_and_interpolate_flip_arrays(V_all, Q_all, t_all, mask, win_size)
    seg = struct();
    if sum(mask) < 2, return; end
    
    V = flipud(V_all(mask)); 
    Q = flipud(Q_all(mask));
    t = flipud(t_all(mask));
     
    if numel(V) >= win_size
        V = movmean(V, win_size);
        Q = movmean(Q, win_size);
        t = movmean(t, win_size);
    end
    
    dV = 0.001;
    [V_sorted, si] = sort(V);
    Q_sorted = Q(si);
    t_sorted = t(si);
    
    [V_unique, ~, ic] = unique(V_sorted);
    Q_unique = accumarray(ic, Q_sorted, [], @mean);
    if numel(V_unique) >= 2
        V_grid = (min(V_unique) : dV : max(V_unique))';
        seg.V_grid = V_grid;
        seg.Q = interp1(V_unique, Q_unique, V_grid, 'linear');
        
        t_unique = accumarray(ic, t_sorted, [], @mean);
        seg.t = interp1(V_unique, t_unique, V_grid, 'linear');
    end
end
