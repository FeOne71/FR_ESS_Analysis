% Phase2_Feature_Extraction.m
% Extracts 25+ features and 3 labels from the parsed field data events
% Labels: Energy Efficiency, Dynamic IR, SOH_BMS
% Features: V, I, T statistics, dynamic rates (dV/dt, dT/dt), and duration

clear; clc; close all;
warning('off', 'all');

%% Configuration
fieldProcessDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Field_Process';
parsedDir = fullfile(fieldProcessDir, 'Parsed_Events');
outputDir = fullfile(fieldProcessDir, 'Features');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Find all parsed event files
parsedFiles = dir(fullfile(parsedDir, 'Parsed_Events_*.mat'));

% Initialize storage for the Master Dataset
MasterData = struct();
master_idx = 1;

%% Feature Extraction Loop
for i = 1:length(parsedFiles)
    inFile = fullfile(parsedFiles(i).folder, parsedFiles(i).name);
    fprintf('Extracting features from %s ...\n', parsedFiles(i).name);
    
    % Extract date string
    token = regexp(parsedFiles(i).name, 'Parsed_Events_(\d{8})\.mat', 'tokens', 'once');
    date_str = token{1};
    
    load(inFile, 'events');
    
    %% 1. Process Charge Events
    if isfield(events, 'Charge')
        chg_names = fieldnames(events.Charge);
        for k = 1:length(chg_names)
            evt = events.Charge.(chg_names{k});
            
            % Skip if data is too short
            if length(evt.t) < 5
                continue;
            end
            
            % Compute Features
            feat = extract_event_features(evt, 'Charge');
            feat.Date = date_str;
            feat.EventType = 'Charge';
            feat.EventID = chg_names{k};
            
            % Store
            MasterData(master_idx).Features = feat;
            master_idx = master_idx + 1;
        end
    end
    
    %% 2. Process Discharge Events
    if isfield(events, 'Discharge')
        dch_names = fieldnames(events.Discharge);
        for k = 1:length(dch_names)
            evt = events.Discharge.(dch_names{k});
            if length(evt.t) < 5, continue; end
            feat = extract_event_features(evt, 'Discharge');
            feat.Date = date_str;
            feat.EventType = 'Discharge';
            feat.EventID = dch_names{k};
            MasterData(master_idx).Features = feat;
            master_idx = master_idx + 1;
        end
    end
end

%% Convert to Table and Save
if master_idx == 1
    fprintf('No valid events found to extract features.\n');
else
    % Flatten struct array into a table
    featureNames = fieldnames(MasterData(1).Features);
    T_data = table();
    for i = 1:length(featureNames)
        fn = featureNames{i};
        colData = cell(length(MasterData), 1);
        for j = 1:length(MasterData)
            val = MasterData(j).Features.(fn);
            if isnumeric(val)
                colData{j} = val;
            else
                colData{j} = {val}; % keep strings in cells
            end
        end
        
        if isnumeric(colData{1})
            T_data.(fn) = cell2mat(colData);
        else
            T_data.(fn) = [colData{:}].';
        end
    end

    save(fullfile(outputDir, 'Master_Feature_Table.mat'), 'T_data', '-v7.3');
    writetable(T_data, fullfile(outputDir, 'Master_Feature_Table.csv'));

    fprintf('\nPhase 2 Feature Extraction Complete. Total Events Processed: %d\n', length(MasterData));
    fprintf('Saved CSV and MAT to %s\n', outputDir);
end

%% Helper Function for Feature Extraction
function feat = extract_event_features(evt, type)
    feat = struct();
    
    % Extract arrays
    V = evt.V;
    I = abs(evt.I); % Use absolute current for statistical magnitude
    T_raw = evt.T;
    t_sec = evt.t_sec;
    dt_arr = diff(t_sec);
    
    valid_T = T_raw(T_raw > 0); % Basic noise filtering for T
    if isempty(valid_T)
        valid_T = [0];
    end
    
    % Time duration feature
    feat.duration_sec = t_sec(end) - t_sec(1);
    
    % ==========================================
    % 1. V Features (Voltage)
    % ==========================================
    feat.V_avg = mean(V, 'omitnan');
    feat.V_max = max(V);
    feat.V_min = min(V);
    feat.V_std = std(V, 'omitnan');
    feat.V_skewness = skewness(V);
    feat.V_kurtosis = kurtosis(V);
    feat.V_diff = V(end) - V(1); % V_end - V_start
    feat.dV_dt_avg = feat.V_diff / max(feat.duration_sec, 1e-6);
    
    % ==========================================
    % 2. I Features (Current & Capacity)
    % ==========================================
    feat.I_avg = mean(I, 'omitnan');
    feat.I_max = max(I);
    feat.I_min = min(I);
    feat.I_std = std(I, 'omitnan');
    
    % Ah Throughput (Integral of I dt)
    % I is in Amperes, t_sec is in seconds. Ah = A * (sec/3600)
    feat.Ah_throughput = trapz(t_sec, I) / 3600; 
    
    % Energy (Wh) for Efficiency Labeling
    % Energy = Integral of V * I dt
    feat.Energy_Wh = trapz(t_sec, V .* I) / 3600;
    
    % ==========================================
    % 3. T Features (Temperature)
    % ==========================================
    feat.T_avg = mean(valid_T, 'omitnan');
    feat.T_max = max(valid_T);
    feat.T_min = min(valid_T);
    feat.T_diff = valid_T(end) - valid_T(1);
    feat.dT_dt_avg = feat.T_diff / max(feat.duration_sec, 1e-6);
    
    % ==========================================
    % 4. Labels & Target Candidates
    % ==========================================
    
    % Label 1: BMS SOH Map
    feat.SOH_BMS_avg = mean(evt.SOH, 'omitnan');
    feat.SOC_start = evt.SOC(1);
    feat.SOC_end = evt.SOC(end);
    
    % Label 2: Dynamic IR (Short-term Local IR at 5 levels: 1s, 3s, 5s, 10s, 30s)
    % Assuming sampling rate dt is approximately 1 second
    % We compute dV/dI from the start of the event to each specific second
    ir_levels = [1, 3, 5, 10, 30];
    for lvl_idx = 1:length(ir_levels)
        sec = ir_levels(lvl_idx);
        ptr = sec + 1; % indices are 1-based, e.g., 1s delta -> index 2 vs index 1
        
        feat_name = sprintf('Local_IR_%ds_Ohm', sec);
        
        if length(I) >= ptr && abs(I(ptr) - I(1)) > 5 % Delta I > 5A logic
            feat.(feat_name) = abs((V(ptr) - V(1)) / (I(ptr) - I(1) + 1e-6));
        else
            feat.(feat_name) = NaN;
        end
    end
end
