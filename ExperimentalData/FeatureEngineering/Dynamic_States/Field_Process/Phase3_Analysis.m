% Phase3_Analysis.m
% Aggregates event-level features into daily records
% Calculates the 3 Target Labels: Energy Efficiency, Standard IR, SOH_BMS
% Performs Correlation Analysis and generates heatmaps

clear; clc; close all;
warning('off', 'all');

%% Configuration
fieldProcessDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Field_Process';
featuresDir = fullfile(fieldProcessDir, 'Features');
outputDir = fullfile(fieldProcessDir, 'Results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Load Master Feature Table
inFile = fullfile(featuresDir, 'Master_Feature_Table.mat');
if ~exist(inFile, 'file')
    error('Master_Feature_Table.mat not found. Please run Phase 2 first.');
end
load(inFile, 'T_data');

%% Daily Aggregation & Label Generation
% We will group by Date
unique_dates = unique(T_data.Date);
num_days = length(unique_dates);

% Preallocate arrays for the Daily Table
DailyFeatures = {};
df_idx = 1;

for i = 1:num_days
    d_str = unique_dates{i};
    day_data = T_data(strcmp(T_data.Date, d_str), :);
    
    chg_data = day_data(strcmp(day_data.EventType, 'Charge'), :);
    dch_data = day_data(strcmp(day_data.EventType, 'Discharge'), :);
    
    % Skip days with no charge or discharge (cannot compute efficiency)
    if isempty(chg_data) || isempty(dch_data)
        continue;
    end
    
    % --- Calculate Y Labels ---
    
    % 1. Energy Efficiency (%)
    total_chg_energy = sum(chg_data.Energy_Wh, 'omitnan');
    total_dch_energy = sum(dch_data.Energy_Wh, 'omitnan');
    if total_chg_energy > 0
        % Discharge energy is negative in some conventions, using abs just in case
        eff = (abs(total_dch_energy) / abs(total_chg_energy)) * 100; 
    else
        eff = NaN;
    end
    
    % 2. Standard IR (Macro Labels for 5 levels)
    ir_levels = [1, 3, 5, 10, 30];
    daily_ir_vals = struct();
    for lvl = ir_levels
        colName = sprintf('Local_IR_%ds_Ohm', lvl);
        if ismember(colName, day_data.Properties.VariableNames)
            valid_vals = day_data.(colName)(~isnan(day_data.(colName)));
            if ~isempty(valid_vals)
                daily_ir_vals.(sprintf('Label_Std_IR_%ds', lvl)) = mean(valid_vals);
            else
                daily_ir_vals.(sprintf('Label_Std_IR_%ds', lvl)) = NaN;
            end
        else
            daily_ir_vals.(sprintf('Label_Std_IR_%ds', lvl)) = NaN;
        end
    end
    
    % 3. SOH_BMS
    soh_bms = mean(day_data.SOH_BMS_avg, 'omitnan');
    
    % Filtering outliers
    if eff > 120 || eff < 40 % Extremely unrealistic efficiency
        eff = NaN;
    end
    
    % --- Calculate X Features (Daily Aggregates) ---
    feat = struct();
    feat.Date = d_str;
    
    % Add the stored IR values into feat
    feat_fields = fieldnames(daily_ir_vals);
    for f_idx = 1:length(feat_fields)
        f_name = feat_fields{f_idx};
        feat.(f_name) = daily_ir_vals.(f_name);
    end
    
    % Voltage
    feat.Mean_V_avg = mean(day_data.V_avg, 'omitnan');
    feat.Max_V_max = max(day_data.V_max);
    feat.Min_V_min = min(day_data.V_min);
    feat.Mean_dV_dt = mean(day_data.dV_dt_avg, 'omitnan');
    
    % Current
    feat.Total_Ah_Throughput = sum(abs(day_data.Ah_throughput), 'omitnan');
    feat.Max_I = max(day_data.I_max);
    
    % Temperature
    feat.Mean_T_avg = mean(day_data.T_avg, 'omitnan');
    feat.Max_T_max = max(day_data.T_max);
    feat.Mean_dT_dt = mean(day_data.dT_dt_avg, 'omitnan');
    
    % Duration
    feat.Total_Active_Time_sec = sum(day_data.duration_sec, 'omitnan');
    
    % Store Labels (others are generated iteratively above)
    feat.Label_Energy_Eff = eff;
    feat.Label_SOH_BMS = soh_bms;
    
    % Append
    DailyFeatures{df_idx} = feat;
    df_idx = df_idx + 1;
end

if df_idx == 1
    fprintf('No valid daily data generated.\n');
    return;
end

% Convert to Table
DailyTable = struct2table([DailyFeatures{:}]);

% Drop rows where key labels are NaN
DailyTable = rmmissing(DailyTable, 'DataVariables', {'Label_Energy_Eff', 'Label_Std_IR_1s', 'Label_SOH_BMS'});

fprintf('Generated %d daily aggregated records.\n', height(DailyTable));

% Save table
save(fullfile(outputDir, 'Daily_Features_Labels.mat'), 'DailyTable');
writetable(DailyTable, fullfile(outputDir, 'Daily_Features_Labels.csv'));

%% Correlation Analysis & Heatmap
if height(DailyTable) > 2
    % Extract numerical features and labels
    num_vars = DailyTable(:, 2:end); % Skip Date string
    var_names = num_vars.Properties.VariableNames;
    data_mat = table2array(num_vars);
    
    % Compute Pearson Correlation
    [R, P] = corrcoef(data_mat, 'Rows', 'pairwise');
    
    % Separate features (X) and labels (Y) for focused heatmap
    % We have 7 labels now: Eff, IR (1,3,5,10,30s), SOH
    % Let's find their indices dynamically to be safe
    label_vars = {'Label_Energy_Eff', 'Label_Std_IR_1s', 'Label_Std_IR_3s', ...
                  'Label_Std_IR_5s', 'Label_Std_IR_10s', 'Label_Std_IR_30s', 'Label_SOH_BMS'};
    
    y_idx = find(ismember(var_names, label_vars));
    x_idx = setdiff(1:length(var_names), y_idx);
    
    R_xy = R(x_idx, y_idx);
    
    % Plot Heatmap
    fig = figure('Name', 'Feature-Label Correlation Heatmap', 'Position', [100, 100, 900, 700], 'Color', 'w');
    
    h = heatmap(var_names(y_idx), var_names(x_idx), R_xy, ...
        'Colormap', parula, ...
        'ColorLimits', [-1 1], ...
        'Title', 'Pearson Correlation: Field Features vs. Target Labels');
    
    h.XLabel = 'Target Labels';
    h.YLabel = 'Extracted Features (Daily Aggregated)';
    
    % Save Plot
    saveas(fig, fullfile(outputDir, 'Correlation_Heatmap.fig'));
    fprintf('  Saved Correlation_Heatmap.fig\n');
else
    fprintf('Not enough daily data points to run correlation analysis.\n');
end

fprintf('\nPhase 3 Analysis Complete.\n');
