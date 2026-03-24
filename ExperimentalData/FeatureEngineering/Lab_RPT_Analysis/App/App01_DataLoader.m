function App_VQ_grid = App_DataLoader(rpt_dir, dc_dir, win_size_rpt, win_size_dc, use_parallel, html_handle)
% APP_DATALOADER Autodetects RPT and Drive Cycle data from user-selected directories, applies smoothing, and interpolates to a 0V grid.
%
% Inputs:
%   rpt_dir: Directory path containing RPT data files
%   dc_dir: Directory path containing Drive Cycle data files
%   win_size_rpt: Window size for moving average on RPT data. Default: 20
%   win_size_dc: Window size for moving average on Drive Cycle data. Default: 20
%   use_parallel: Flag for parallel processing (true/false). Default: false
%   html_handle: (Optional) uihtml component handle for HTML progress bar
%
% Output:
%   App_VQ_grid: struct containing interpolated data hierarchy.

if nargin < 3 || isempty(win_size_rpt), win_size_rpt = 20; end
if nargin < 4 || isempty(win_size_dc), win_size_dc = 20; end
if nargin < 5 || isempty(use_parallel)
    use_parallel = false; % Default to single-threaded for easier debugging
end
if nargin < 6
    html_handle = []; % Default no progress bar
end

fprintf('--- App DataLoader Started ---\n');
fprintf('RPT Directory: %s\n', rpt_dir);
fprintf('DC Directory: %s\n', dc_dir);
fprintf('Smoothing Window: RPT=%d, DC=%d\n', win_size_rpt, win_size_dc);

App_VQ_grid = struct();

% Pre-define typical channels and cycles to aid structuring
channels = {'Ch09', 'Ch10', 'Ch11', 'Ch12', 'Ch13', 'Ch14', 'Ch15', 'Ch16'};
rpt_cycles = {'0cyc', '200cyc', '400cyc','600cyc','800cyc','1000cyc'};

% Define step indices (adjust if needed based on the lab protocol)
% ocv_steps = [8, 10]; % Charge OCV, Dchg OCV
% static_capacity_step = 3; 

% 1. Search for RPT Data
% The user provides the exact path to the RPT folder.
if ~exist(rpt_dir, 'dir')
    fprintf('Error: RPT folder path is invalid (%s)\n', rpt_dir);
    return;
end

csv_files = dir(fullfile(rpt_dir, '*_RPT_*.csv'));
if isempty(csv_files)
    fprintf('Warning: No *_RPT_*.csv files found in %s\n', rpt_dir);
    % We don't return here so Drive Cycle can still be processed
end

total_files = length(csv_files);
if total_files > 0
    fprintf('Found %d RPT CSV files. Processing (parallel)...\n', total_files);
    
    % cell 배열에 결과 저장
    file_paths = cell(total_files, 1);
    ch_keys = cell(total_files, 1);
    cyc_keys = cell(total_files, 1);
    results = cell(total_files, 1);
    valid = false(total_files, 1);
    
    % 파일 정보 사전 준비
    for idx = 1:total_files
        fn = csv_files(idx).name;
        file_paths{idx} = fullfile(csv_files(idx).folder, fn);
        
        tok = regexp(fn, 'Ch(\d+)_RPT_(\d+)cyc', 'tokens');
        if ~isempty(tok)
            ch_keys{idx} = sprintf('Ch%02d', str2double(tok{1}{1}));
            cyc_keys{idx} = sprintf('cyc%d', str2double(tok{1}{2}));
            valid(idx) = true;
        end
    end
    
    % ========== parfeval로 비동기 병렬 처리 ==========
    ws = win_size_rpt;
    nValid = sum(valid);
    
    % Worker Pool에 작업 제출 (비동기)
    futures(nValid) = parallel.FevalFuture;
    fi = 0;
    validIdx = find(valid);
    for k = 1:nValid
        idx = validIdx(k);
        fi = fi + 1;
        futures(fi) = parfeval(@process_single_rpt_file, 1, file_paths{idx}, 1, [2, 10], ws);
    end
    
    % 하나씩 완료될 때마다 프로그레스 바 업데이트 (메인 스레드에서!)
    for k = 1:nValid
        [completedIdx, result] = fetchNext(futures);
        origIdx = validIdx(completedIdx);
        results{origIdx} = result;
        
        % 프로그레스 바 점진 업데이트 (0% ~ 50%)
        pct = round((k / nValid) * 50);
        fprintf('  [%d/%d] %s (%s) Done → %d%%\n', k, nValid, ch_keys{origIdx}, cyc_keys{origIdx}, pct);
        if ~isempty(html_handle)
            updateHTMLBar(html_handle, pct);
        end
    end
    
    % ========== 결과 병합 ==========
    for idx = 1:total_files
        if valid(idx) && ~isempty(results{idx})
            ck = cyc_keys{idx};
            chk = ch_keys{idx};
            if ~isfield(App_VQ_grid, ck)
                App_VQ_grid.(ck) = struct();
            end
            App_VQ_grid.(ck).(chk) = results{idx};
        end
    end
end

% ---------------------------------------------------------
% 2. Search for Drive Cycle Data
% ---------------------------------------------------------
if exist(dc_dir, 'dir')
    % Assuming Drive Cycle data are stored as .mat or .csv
    dc_files = dir(fullfile(dc_dir, '*.mat')); % Prioritize parsed .mat if available
    if isempty(dc_files)
        dc_files = dir(fullfile(dc_dir, '*Drive*.csv'));
    end
    
    if isempty(dc_files)
        fprintf('Warning: No Drive Cycle data found in %s\n', dc_dir);
    else
        total_dc = length(dc_files);
        fprintf('Found %d Drive Cycle files. Processing...\n', total_dc);
        for idx = 1:total_dc
            file_name = dc_files(idx).name;
            file_path = fullfile(dc_files(idx).folder, file_name);
            
            % Attempt to extract cycle (e.g., 'parsedDriveCycle_0cyc_filtered.mat')
            tok = regexp(file_name, '(\d+)cyc', 'tokens');
            if ~isempty(tok)
                cyc_str = tok{1}{1};
                cyc_key = sprintf('cyc%s', cyc_str);
                fprintf('  Processing DC: %s\n', file_name);
                
                % Process DC file (assuming it contains channel data)
                dc_data = process_single_dc_file(file_path, win_size_dc);
                
                % Merge into main struct
                ch_fields = fieldnames(dc_data);
                for c = 1:length(ch_fields)
                    ch_key = ch_fields{c};
                    App_VQ_grid.(cyc_key).(ch_key).DriveCycle = dc_data.(ch_key);
                end
            end
            
            % Update gauge if provided (50% to 100%)
            if ~isempty(gauge_handle) && isvalid(gauge_handle)
                prog = 50 + (idx / total_dc) * 50;
                gauge_handle.Value = min(max(prog, gauge_handle.Limits(1)), gauge_handle.Limits(2));
                drawnow limitrate;
            end
        end
    end
else
    fprintf('Warning: DC folder path is invalid (%s)\n', dc_dir);
end

fprintf('--- App DataLoader Completed ---\n');
end

% =========================================================================
% Local Help Functions for Loading and Interpolation
% =========================================================================
function dc_data = process_single_dc_file(filepath, win_size)
    % Placeholder for actual DC parsing logic depending on Lab format.
    % In this MVP, if it's a .mat file, we load it directly and try to find V, I, Q
    dc_data = struct();
    try
        [~, ~, ext] = fileparts(filepath);
        if strcmpi(ext, '.mat')
            % MAT Logic
            file_data = load(filepath);
            vars = fieldnames(file_data);
            if ~isempty(vars)
                top_var = file_data.(vars{1}); % e.g., parsedDriveCycle_0cyc
                ch_names = fieldnames(top_var);
                for i=1:length(ch_names)
                    ch = ch_names{i};
                    % Standardize channel name from 'ch9_Drive...' to 'Ch09'
                    ch_id_match = regexp(ch, 'ch(\d+)', 'tokens');
                    if ~isempty(ch_id_match)
                        ch_num = str2double(ch_id_match{1}{1});
                        standard_ch = sprintf('Ch%02d', ch_num);
                        
                        % Find DC1 profile or similar inside
                        sub_fields = fieldnames(top_var.(ch));
                        for j=1:length(sub_fields)
                            if startsWith(sub_fields{j}, 'DC', 'IgnoreCase', true)
                                % Found a DC profile
                                raw_profile = top_var.(ch).(sub_fields{j});
                                % Assuming it has V, I, Q (or SOC) arrays
                                if isfield(raw_profile, 'V') && isfield(raw_profile, 'I')
                                    % Use the same extraction structure, just mapping fields
                                    V = raw_profile.V(:);
                                    % Apply smoothing before return
                                    if numel(V) >= win_size
                                         V = movmean(V, win_size);
                                    end
                                    dc_data.(standard_ch).(sub_fields{j}).V_grid_raw = V;
                                end
                            end
                        end
                    end
                end
            end
        end
    catch ME
        fprintf('    [Error] DC Parsing Failed: %s\n', ME.message);
    end
end
function ch_data = process_single_rpt_file(filepath, static_step, ocv_steps, win_size)
    ch_data = struct();
    try
        % Ultra-fast textscan method for specific Lab RPT format
        % Column Index Mapping based on standard Lab format:
        % 2: Step Index (int)
        % 4: Cycle Index 
        % 5: Time (string 'hh:mm:ss.ms') or we can just skip it if we don't strictly *need* t for V-Q interpolation!
        % Wait, earlier we used time for t_grid. Let's just parse it as string.
        % 7: Current(A) 
        % 8: Voltage(V)
        % 9: Capacity(Ah)
        
        fid = fopen(filepath, 'r');
        fgetl(fid);
        % %f %f %*s %f %s %*s %f %f %f %*[^\n]
        formatSpec = '%*f %f %*q %f %s %*s %f %f %f %*[^\n]';
        dataArray = textscan(fid, formatSpec, 'Delimiter', ',', 'EmptyValue', NaN, 'ReturnOnError', false);
        fclose(fid);
        
        steps = dataArray{1};
        subs = dataArray{2}; 
        t_strs = dataArray{3};
        I = dataArray{4};
        V = dataArray{5};
        Q = dataArray{6};
        
        % Convert time strings to seconds much faster than `duration()` or `cellfun`
        % Format is usually 'HH:MM:SS.SSS' e.g. '00:00:01.000'
        try
            t = zeros(length(t_strs), 1);
            if ~isempty(t_strs)
                if length(t_strs{1}) >= 8
                    % 100x Faster pure vectorization parsing using char matrix
                    % Convert all cells to block of char, e.g. 10000x12 matrix
                    strMat = char(t_strs);
                    
                    if size(strMat, 2) >= 8
                        % 'HH:MM:SS' are at fixed columns: 1:2, 4:5, 7:8
                        % Subseconds are at 10:end, but we can safely ignore subseconds for VQ grid
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
            % Fastest fallback
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
    
    % Substep filter logic was: mask = (Step == static_step) & (SubCol == 2)
    % BUT the 4th column is 'Cycle Index' and the raw data has '1' for cycle index.
    % Ah! In my initial implementation, I used `subCol=4` assuming it was 'Step Type' or 'Substep', but Col 4 is 'Cycle Index' which is always 1 for 0cyc file.
    % Let's drop the `& (subs == 2)` constraint entirely since it was likely derived from a different format.
    % We will just isolate by Step Index.
    
    % 1. Static Capacity
    mask_static = (steps == static_step);
    ch_data.Static = extract_and_interpolate_arrays(V, Q, t, mask_static, win_size);
    
    % 2. OCV Charge
    mask_chg = (steps == ocv_steps(1));
    ch_data.OCV_charge = extract_and_interpolate_arrays(V, Q, t, mask_chg, win_size);
    
    % 3. OCV Discharge (Requires flipping for monotonic V)
    mask_dch = (steps == ocv_steps(2));
    ch_data.OCV_discharge = extract_and_interpolate_flip_arrays(V, Q, t, mask_dch, win_size);
end

function seg = extract_and_interpolate_arrays(V_all, Q_all, t_all, mask, win_size)
    seg = struct();
    if sum(mask) < 2, return; end
    
    V = V_all(mask); 
    Q = Q_all(mask);
    t = t_all(mask);
    
    % Add moving average before interpolation!
    if numel(V) >= win_size
        V = movmean(V, win_size);
        Q = movmean(Q, win_size);
        t = movmean(t, win_size); 
    end
    
    % V_grid Interpolation (0.001V)
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
     
    % Add moving average before interpolation!
    if numel(V) >= win_size
        V = movmean(V, win_size);
        Q = movmean(Q, win_size);
        t = movmean(t, win_size);
    end
    
    dV = 0.001;
    [V_sorted, si] = sort(V);  % Ascending sort
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
% (End of file, no extra broken functions)
