function FeatureTable = App_FeatureExtractor(App_VQ_grid, master_ruler_path)
% APP_FEATUREEXTRACTOR Extracts features and labels from interpolated V-Q data
%
% Matches logic from RPT_Feature_Extraction_Advanced.m exactly:
%   - Features extracted from 5 C-rates: c01, c05, c1, c2, c3
%   - Labels: SOH (Static capacity), LLI/LAM (OCV_charge vs Fresh)
%   - Normalization: C-rate-independent Z-score
%
% Features (14): [Chg_dQ(1-5), Dch_dQ(1-5), Chg_PkH, Chg_PkA, Dch_PkH, Dch_PkA]
% Labels (3): [SOH, LLI, LAM]
%
% Output:
%   FeatureTable: table with CellID, Cycle, CrateLabel, CrateNum,
%                 X_Features, Y_Labels, X_Normalized

fprintf('--- App FeatureExtractor Started ---\n');

if nargin < 2 || isempty(master_ruler_path)
    master_ruler_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\Dataset\MasterRulers.mat';
end

if ~exist(master_ruler_path, 'file')
    fprintf('Error: Master Ruler file not found at %s\n', master_ruler_path);
    FeatureTable = table();
    return;
end

% 1. Load Master Ruler
loaded = load(master_ruler_path);
MasterRulers = loaded.MasterRulers;
fprintf('Master Ruler loaded successfully.\n');

% Standard voltage windows
win_chg_min = 3.7; win_chg_max = 3.95;
win_dch_min = 3.75; win_dch_max = 3.88;

% C-rate configuration (same as RPT_Feature_Extraction_Advanced.m)
target_crates = {'c01', 'c05', 'c1', 'c2', 'c3'};
target_crates_val = [0.1, 0.5, 1, 2, 3];

% Arrays for table construction
data_CellID = {};
data_Cycle = [];
data_Crate = {};
data_CrateNum = [];
data_Features = [];
data_Labels = [];
cnt = 0;

% 2. Extract Features & Labels
cycles = fieldnames(App_VQ_grid);
for i_ch_global = 1:100  % iterate channels across all cycles
    % Find all unique channels
    break;
end

% Collect all channels
all_channels = {};
for c = 1:length(cycles)
    ch_list = fieldnames(App_VQ_grid.(cycles{c}));
    all_channels = union(all_channels, ch_list);
end
channels = sort(all_channels);

for i = 1:length(channels)
    ch = channels{i};
    
    if ~isfield(MasterRulers, ch)
        continue;
    end
    
    V_ruler_chg = MasterRulers.(ch).V_bounds_chg;
    V_ruler_dch = MasterRulers.(ch).V_bounds_dch;
    
    % Get Q_rated from cyc0 Static (same as existing script)
    Q_rated_ch = 55.6;  % default
    if isfield(App_VQ_grid, 'cyc0') && isfield(App_VQ_grid.cyc0, ch) && ...
       isfield(App_VQ_grid.cyc0.(ch), 'Static') && isfield(App_VQ_grid.cyc0.(ch).Static, 'Q')
        q_static = App_VQ_grid.cyc0.(ch).Static.Q;
        if ~isempty(q_static)
            Q_rated_ch = max(abs(q_static));
            if Q_rated_ch == 0, Q_rated_ch = 55.6; end
        end
    end
    
    for c = 1:length(cycles)
        cyc_key = cycles{c};
        cyc_num = sscanf(cyc_key, 'cyc%d');
        
        if ~isfield(App_VQ_grid.(cyc_key), ch)
            continue;
        end
        
        ch_data = App_VQ_grid.(cyc_key).(ch);
        
        % --- Iterate through C-rates (same as existing script) ---
        for r = 1:length(target_crates)
            crate_label = target_crates{r};
            crate_val = target_crates_val(r);
            f_chg = [crate_label '_charge'];
            f_dch = [crate_label '_discharge'];
            
            if ~isfield(ch_data, f_chg) || ~isfield(ch_data, f_dch)
                continue;
            end
            
            % Check data validity
            if ~isfield(ch_data.(f_chg), 'V_grid') || ~isfield(ch_data.(f_dch), 'V_grid')
                continue;
            end
            
            % Feature extraction
            [fea_chg_dQ, fea_chg_diff] = extract_features_half(ch_data.(f_chg), ...
                V_ruler_chg, 'charge', win_chg_min, win_chg_max);
            [fea_dch_dQ, fea_dch_diff] = extract_features_half(ch_data.(f_dch), ...
                V_ruler_dch, 'discharge', win_dch_min, win_dch_max);
            
            row_features = [fea_chg_dQ, fea_dch_dQ, fea_chg_diff, fea_dch_diff];
            
            % --- Label 1: SOH (Static Capacity) ---
            if isfield(ch_data, 'Static') && isfield(ch_data.Static, 'Q') && ~isempty(ch_data.Static.Q)
                lbl_soh = max(abs(ch_data.Static.Q));
            else
                lbl_soh = NaN;
            end
            
            % --- Label 2 & 3: LLI, LAM (OCV_charge vs Fresh_OCV_Charge) ---
            if isfield(ch_data, 'OCV_charge') && isfield(ch_data.OCV_charge, 'V_grid') && ...
               isfield(MasterRulers.(ch), 'Fresh_OCV_Charge') && ~isempty(MasterRulers.(ch).Fresh_OCV_Charge)
                [lbl_lli, lbl_lam] = analyze_ica_aging(ch_data.OCV_charge, ...
                    MasterRulers.(ch).Fresh_OCV_Charge, Q_rated_ch);
            else
                lbl_lli = NaN;
                lbl_lam = NaN;
            end
            
            % Store
            cnt = cnt + 1;
            data_CellID{cnt, 1} = ch;
            data_Cycle(cnt, 1) = cyc_num;
            data_Crate{cnt, 1} = crate_label;
            data_CrateNum(cnt, 1) = crate_val;
            data_Features(cnt, :) = row_features;
            data_Labels(cnt, :) = [lbl_soh, lbl_lli, lbl_lam];
        end
    end
end

% 3. Create Table
if isempty(data_Features)
    fprintf('Warning: No features were extracted. Check input data structure.\n');
    FeatureTable = table();
    return;
end

FeatureTable = table(data_CellID, data_Cycle, data_Crate, data_CrateNum, ...
    data_Features, data_Labels, 'VariableNames', ...
    {'CellID', 'Cycle', 'CrateLabel', 'CrateNum', 'X_Features', 'Y_Labels'});

% 4. C-rate-independent Z-score Normalization (same as existing script)
FeatureTable.X_Normalized = FeatureTable.X_Features;
for r = 1:length(target_crates)
    c_val = target_crates_val(r);
    idx = FeatureTable.CrateNum == c_val;
    X_sub = FeatureTable.X_Features(idx, :);
    
    if isempty(X_sub), continue; end
    
    mu = mean(X_sub, 1, 'omitnan');
    sigma = std(X_sub, 0, 1, 'omitnan');
    sigma(sigma==0) = 1;
    FeatureTable.X_Normalized(idx, :) = (X_sub - mu) ./ sigma;
end

fprintf('Extracted %d samples (5 C-rates): 14 features + 3 labels.\n', height(FeatureTable));
fprintf('--- App FeatureExtractor Completed ---\n');
end

% =========================================================================
% Local Helper Functions
% =========================================================================
function [features_dQ, features_diff] = extract_features_half(data_struct, V_ruler, mode, minV, maxV)
    V_grid = data_struct.V_grid;
    Q_val  = data_struct.Q;
    [V_u, uid] = unique(V_grid); 
    Q_u = Q_val(uid);
    
    if numel(V_u) < 2
         features_dQ = nan(1, length(V_ruler)-1);
         features_diff = [NaN, NaN];
         return;
    end
    
    % Capacity segments delta Q
    Q_ruler = interp1(V_u, Q_u, V_ruler, 'linear');
    features_dQ = abs(diff(Q_ruler));
    
    % dQ/dV characteristics
    dV = gradient(V_u);
    dQ = gradient(Q_u);
    dV(dV==0) = NaN;
    dQdV_raw = dQ ./ dV;
    
    window_size = 20;
    dQdV_filt = movmean(dQdV_raw, window_size);
    
    if strcmp(mode, 'charge')
        mask = V_u >= minV & V_u <= maxV;
    else
        mask = V_u <= maxV & V_u >= minV;
    end
    
    dQdV_win = dQdV_filt(mask);
    V_win = V_u(mask);
    
    if isempty(dQdV_win)
        features_diff = [NaN, NaN];
    else
        pk_height = max(dQdV_win);
        pk_area = trapz(V_win, dQdV_win);
        features_diff = [pk_height, abs(pk_area)];
    end
end

function [lli, lam] = analyze_ica_aging(curr_ocv_struct, fresh_ocv_struct, Q_rated)
    if nargin < 3, Q_rated = 55.6; end
    
    if isempty(curr_ocv_struct) || isempty(fresh_ocv_struct)
        lli = NaN;
        lam = NaN;
        return;
    end
    
    roi_min = 3.40;
    roi_max = 4.00;
    [peak_V_fresh, peak_H_fresh] = get_main_peak(fresh_ocv_struct, roi_min, roi_max);
    [peak_V_aged, peak_H_aged] = get_main_peak(curr_ocv_struct, roi_min, roi_max);
    
    if isnan(peak_V_fresh) || isnan(peak_V_aged)
        lli = NaN;
        lam = NaN;
    else
        % LLI (%) = |dV_peak| * H_peak_fresh / Q_rated * 100
        lli = (abs(peak_V_aged - peak_V_fresh) * peak_H_fresh / Q_rated) * 100;
        
        % LAM (%): Integrate dQ/dV area around peak [3.4~4.0V]
        mask_f = fresh_ocv_struct.V_grid >= 3.4 & fresh_ocv_struct.V_grid <= 4.0;
        mask_a = curr_ocv_struct.V_grid >= 3.4 & curr_ocv_struct.V_grid <= 4.0;
        
        dQdV_f = gradient(fresh_ocv_struct.Q)./gradient(fresh_ocv_struct.V_grid);
        dQdV_f(isinf(dQdV_f)) = NaN; dQdV_f = fillmissing(dQdV_f, 'linear');
        dQdV_a = gradient(curr_ocv_struct.Q)./gradient(curr_ocv_struct.V_grid);
        dQdV_a(isinf(dQdV_a)) = NaN; dQdV_a = fillmissing(dQdV_a, 'linear');
        
        area_peak_fresh = trapz(fresh_ocv_struct.V_grid(mask_f), abs(movmean(dQdV_f(mask_f), 21)));
        area_peak_aged  = trapz(curr_ocv_struct.V_grid(mask_a),  abs(movmean(dQdV_a(mask_a), 21)));
        
        lam = (1 - area_peak_aged / area_peak_fresh) * 100;
    end
end

function [pk_V, pk_H] = get_main_peak(data_struct, minV, maxV)
    V = data_struct.V_grid;
    Q = data_struct.Q;
    
    if numel(V) < 10
        pk_V = NaN;
        pk_H = NaN;
        return;
    end
    
    [V_u, uid] = unique(V);
    Q_u = Q(uid);
    dV = gradient(V_u);
    dQ = gradient(Q_u);
    dV(dV==0) = NaN;
    dQdV = dQ ./ dV;
    dQdV(isinf(dQdV)) = NaN;
    dQdV = fillmissing(dQdV, 'linear');
    
    window_size = 21;
    dQdV_filt = movmean(dQdV, window_size);
    
    mask = V_u >= minV & V_u <= maxV;
    V_roi = V_u(mask);
    dH_roi = dQdV_filt(mask);
    
    if isempty(V_roi)
        pk_V = NaN;
        pk_H = NaN;
        return;
    end
    
    [pks, locs] = findpeaks(dH_roi);
    
    if isempty(pks)
        [pk_H, idx] = max(dH_roi); 
        pk_V = V_roi(idx);
    else
        [pk_H, idx] = max(pks);
        pk_V = V_roi(locs(idx));
    end
end
