%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCIR Statistical Comparison (0cyc vs 200cyc)
% SOC90의 DC1-DC8에 대해 시간별 DCIR 값을 p-value로 비교
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Load Data
fprintf('Loading DCIR data...\n');
load('Lab_DC_DCIR_0cyc_Events.mat');           % 0cyc data
load('Lab_DC_DCIR_200cyc_Events.mat');    % 200cyc data

% Check what variables are loaded
fprintf('=== Loaded variables ===\n');
vars = who;
for i = 1:length(vars)
    fprintf('%s\n', vars{i});
end

% Check if variables exist and use correct names
if exist('Lab_DC_DCIR_0cyc', 'var')
    fprintf('\nUsing Lab_DC_DCIR_0cyc for 0cyc data\n');
    data_0cyc = Lab_DC_DCIR_0cyc;
else
    error('Lab_DC_DCIR_0cyc variable not found in 0cyc file!');
end

if exist('Lab_DC_DCIR_200cyc', 'var')
    fprintf('Using Lab_DC_DCIR_200cyc for 200cyc data\n');
    data_200cyc = Lab_DC_DCIR_200cyc;
else
    error('Lab_DC_DCIR_200cyc variable not found in 200cyc file!');
end

%% Settings
dt_list = [1, 3, 5, 10, 30, 50];  % DCIR time intervals
soc_level = 'SOC50';
dc_profiles = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

% Extract channel names
channels_0cyc = {};
fields_0cyc = fieldnames(data_0cyc);
for i = 1:length(fields_0cyc)
    if contains(fields_0cyc{i}, '_ChgEvent')
        channel_name = strrep(fields_0cyc{i}, '_Drive_0cyc_ChgEvent', '');
        channels_0cyc{end+1} = channel_name;
    end
end
channels_0cyc = unique(channels_0cyc);

channels_200cyc = {};
fields_200cyc = fieldnames(data_200cyc);
for i = 1:length(fields_200cyc)
    if contains(fields_200cyc{i}, '_ChgEvent')
        channel_name = strrep(fields_200cyc{i}, '_Drive_200cyc_ChgEvent', '');
        channels_200cyc{end+1} = channel_name;
    end
end
channels_200cyc = unique(channels_200cyc);

fprintf('0cyc channels: ');
for i = 1:length(channels_0cyc)
    fprintf('%s ', channels_0cyc{i});
end
fprintf('\n');

fprintf('200cyc channels: ');
for i = 1:length(channels_200cyc)
    fprintf('%s ', channels_200cyc{i});
end
fprintf('\n');

%% Collect DCIR data for each time interval
fprintf('\n=== Collecting DCIR data ===\n');

% Initialize results structure
comparison_results = struct();

% Get common channels between 0cyc and 200cyc
common_channels = intersect(channels_0cyc, channels_200cyc);
fprintf('Common channels for comparison: ');
for i = 1:length(common_channels)
    fprintf('%s ', common_channels{i});
end
fprintf('\n');

for dt_idx = 1:length(dt_list)
    dt_sec = dt_list(dt_idx);
    field_name = sprintf('DCIR_%ds', dt_sec);
    
    fprintf('\n--- %ds DCIR Analysis ---\n', dt_sec);
    
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        
        fprintf('  %s:\n', dc_profile);
        
        for ch_idx = 1:length(common_channels)
            channel = common_channels{ch_idx};
            
            % === CHARGING EVENTS ===
            % 0cyc charging data for this channel and profile
            dcir_0cyc_chg = [];
            chg_struct_name_0cyc = sprintf('%s_Drive_0cyc_ChgEvent', channel);
            
            % Collect 0cyc charging events
            if isfield(data_0cyc, chg_struct_name_0cyc) && isfield(data_0cyc.(chg_struct_name_0cyc), soc_level) && isfield(data_0cyc.(chg_struct_name_0cyc).(soc_level), dc_profile)
                events = fieldnames(data_0cyc.(chg_struct_name_0cyc).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = data_0cyc.(chg_struct_name_0cyc).(soc_level).(dc_profile).(evt_name);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_0cyc_chg = [dcir_0cyc_chg; evt_data.(field_name).val];
                    end
                end
            end
            
            % 200cyc charging data for this channel and profile
            dcir_200cyc_chg = [];
            chg_struct_name_200cyc = sprintf('%s_Drive_200cyc_ChgEvent', channel);
            
            % Collect 200cyc charging events
            if isfield(data_200cyc, chg_struct_name_200cyc) && isfield(data_200cyc.(chg_struct_name_200cyc), soc_level) && isfield(data_200cyc.(chg_struct_name_200cyc).(soc_level), dc_profile)
                events = fieldnames(data_200cyc.(chg_struct_name_200cyc).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = data_200cyc.(chg_struct_name_200cyc).(soc_level).(dc_profile).(evt_name);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_200cyc_chg = [dcir_200cyc_chg; evt_data.(field_name).val];
                    end
                end
            end
            
            % === DISCHARGING EVENTS ===
            % 0cyc discharging data for this channel and profile
            dcir_0cyc_dchg = [];
            dchg_struct_name_0cyc = sprintf('%s_Drive_0cyc_DchEvent', channel);
            
            % Collect 0cyc discharging events
            if isfield(data_0cyc, dchg_struct_name_0cyc) && isfield(data_0cyc.(dchg_struct_name_0cyc), soc_level) && isfield(data_0cyc.(dchg_struct_name_0cyc).(soc_level), dc_profile)
                events = fieldnames(data_0cyc.(dchg_struct_name_0cyc).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = data_0cyc.(dchg_struct_name_0cyc).(soc_level).(dc_profile).(evt_name);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_0cyc_dchg = [dcir_0cyc_dchg; evt_data.(field_name).val];
                    end
                end
            end
            
            % 200cyc discharging data for this channel and profile
            dcir_200cyc_dchg = [];
            dchg_struct_name_200cyc = sprintf('%s_Drive_200cyc_DchEvent', channel);
            
            % Collect 200cyc discharging events
            if isfield(data_200cyc, dchg_struct_name_200cyc) && isfield(data_200cyc.(dchg_struct_name_200cyc), soc_level) && isfield(data_200cyc.(dchg_struct_name_200cyc).(soc_level), dc_profile)
                events = fieldnames(data_200cyc.(dchg_struct_name_200cyc).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = data_200cyc.(dchg_struct_name_200cyc).(soc_level).(dc_profile).(evt_name);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_200cyc_dchg = [dcir_200cyc_dchg; evt_data.(field_name).val];
                    end
                end
            end
            
            % === STATISTICAL COMPARISON FOR CHARGING ===
            if ~isempty(dcir_0cyc_chg) && ~isempty(dcir_200cyc_chg)
                % Normality test (only if we have enough samples)
                h_0cyc_chg = 0; h_200cyc_chg = 0;  % Default to normal
                
                if length(dcir_0cyc_chg) >= 4
                    [h_0cyc_chg, ~] = lillietest(dcir_0cyc_chg);
                end
                
                if length(dcir_200cyc_chg) >= 4
                    [h_200cyc_chg, ~] = lillietest(dcir_200cyc_chg);
                end
                
                % Choose test based on normality
                if h_0cyc_chg == 0 && h_200cyc_chg == 0  % Both normal
                    [h, p_value_chg] = ttest2(dcir_0cyc_chg, dcir_200cyc_chg);
                    test_type_chg = 't-test';
                else  % Non-normal
                    [p_value_chg, h] = ranksum(dcir_0cyc_chg, dcir_200cyc_chg);
                    test_type_chg = 'Wilcoxon rank-sum';
                end
                
                % Calculate effect size (Cohen's d)
                pooled_std_chg = sqrt(((length(dcir_0cyc_chg)-1)*var(dcir_0cyc_chg) + (length(dcir_200cyc_chg)-1)*var(dcir_200cyc_chg)) / (length(dcir_0cyc_chg) + length(dcir_200cyc_chg) - 2));
                cohens_d_chg = (mean(dcir_200cyc_chg) - mean(dcir_0cyc_chg)) / pooled_std_chg;
                
                % Store results for charging
                comparison_results.(field_name).(dc_profile).(channel).charging.p_value = p_value_chg;
                comparison_results.(field_name).(dc_profile).(channel).charging.h = h;
                comparison_results.(field_name).(dc_profile).(channel).charging.test_type = test_type_chg;
                comparison_results.(field_name).(dc_profile).(channel).charging.cohens_d = cohens_d_chg;
                comparison_results.(field_name).(dc_profile).(channel).charging.mean_0cyc = mean(dcir_0cyc_chg);
                comparison_results.(field_name).(dc_profile).(channel).charging.mean_200cyc = mean(dcir_200cyc_chg);
                comparison_results.(field_name).(dc_profile).(channel).charging.std_0cyc = std(dcir_0cyc_chg);
                comparison_results.(field_name).(dc_profile).(channel).charging.std_200cyc = std(dcir_200cyc_chg);
                comparison_results.(field_name).(dc_profile).(channel).charging.n_0cyc = length(dcir_0cyc_chg);
                comparison_results.(field_name).(dc_profile).(channel).charging.n_200cyc = length(dcir_200cyc_chg);
                
                significance_chg = '';
                if p_value_chg < 0.001
                    significance_chg = '***';
                elseif p_value_chg < 0.01
                    significance_chg = '**';
                elseif p_value_chg < 0.05
                    significance_chg = '*';
                end
                
                fprintf('    %s (Charging): p=%.4f%s (%s), d=%.3f, 0cyc=%.2f±%.2f (n=%d), 200cyc=%.2f±%.2f (n=%d)\n', ...
                    channel, p_value_chg, significance_chg, test_type_chg, cohens_d_chg, ...
                    mean(dcir_0cyc_chg), std(dcir_0cyc_chg), length(dcir_0cyc_chg), ...
                    mean(dcir_200cyc_chg), std(dcir_200cyc_chg), length(dcir_200cyc_chg));
            else
                fprintf('    %s (Charging): No data available\n', channel);
            end
            
            % === STATISTICAL COMPARISON FOR DISCHARGING ===
            if ~isempty(dcir_0cyc_dchg) && ~isempty(dcir_200cyc_dchg)
                % Normality test (only if we have enough samples)
                h_0cyc_dchg = 0; h_200cyc_dchg = 0;  % Default to normal
                
                if length(dcir_0cyc_dchg) >= 4
                    [h_0cyc_dchg, ~] = lillietest(dcir_0cyc_dchg);
                end
                
                if length(dcir_200cyc_dchg) >= 4
                    [h_200cyc_dchg, ~] = lillietest(dcir_200cyc_dchg);
                end
                
                % Choose test based on normality
                if h_0cyc_dchg == 0 && h_200cyc_dchg == 0  % Both normal
                    [h, p_value_dchg] = ttest2(dcir_0cyc_dchg, dcir_200cyc_dchg);
                    test_type_dchg = 't-test';
                else  % Non-normal
                    [p_value_dchg, h] = ranksum(dcir_0cyc_dchg, dcir_200cyc_dchg);
                    test_type_dchg = 'Wilcoxon rank-sum';
                end
                
                % Calculate effect size (Cohen's d)
                pooled_std_dchg = sqrt(((length(dcir_0cyc_dchg)-1)*var(dcir_0cyc_dchg) + (length(dcir_200cyc_dchg)-1)*var(dcir_200cyc_dchg)) / (length(dcir_0cyc_dchg) + length(dcir_200cyc_dchg) - 2));
                cohens_d_dchg = (mean(dcir_200cyc_dchg) - mean(dcir_0cyc_dchg)) / pooled_std_dchg;
                
                % Store results for discharging
                comparison_results.(field_name).(dc_profile).(channel).discharging.p_value = p_value_dchg;
                comparison_results.(field_name).(dc_profile).(channel).discharging.h = h;
                comparison_results.(field_name).(dc_profile).(channel).discharging.test_type = test_type_dchg;
                comparison_results.(field_name).(dc_profile).(channel).discharging.cohens_d = cohens_d_dchg;
                comparison_results.(field_name).(dc_profile).(channel).discharging.mean_0cyc = mean(dcir_0cyc_dchg);
                comparison_results.(field_name).(dc_profile).(channel).discharging.mean_200cyc = mean(dcir_200cyc_dchg);
                comparison_results.(field_name).(dc_profile).(channel).discharging.std_0cyc = std(dcir_0cyc_dchg);
                comparison_results.(field_name).(dc_profile).(channel).discharging.std_200cyc = std(dcir_200cyc_dchg);
                comparison_results.(field_name).(dc_profile).(channel).discharging.n_0cyc = length(dcir_0cyc_dchg);
                comparison_results.(field_name).(dc_profile).(channel).discharging.n_200cyc = length(dcir_200cyc_dchg);
                
                significance_dchg = '';
                if p_value_dchg < 0.001
                    significance_dchg = '***';
                elseif p_value_dchg < 0.01
                    significance_dchg = '**';
                elseif p_value_dchg < 0.05
                    significance_dchg = '*';
                end
                
                fprintf('    %s (Discharging): p=%.4f%s (%s), d=%.3f, 0cyc=%.2f±%.2f (n=%d), 200cyc=%.2f±%.2f (n=%d)\n', ...
                    channel, p_value_dchg, significance_dchg, test_type_dchg, cohens_d_dchg, ...
                    mean(dcir_0cyc_dchg), std(dcir_0cyc_dchg), length(dcir_0cyc_dchg), ...
                    mean(dcir_200cyc_dchg), std(dcir_200cyc_dchg), length(dcir_200cyc_dchg));
            else
                fprintf('    %s (Discharging): No data available\n', channel);
            end
        end
    end
end

%% Create summary table
fprintf('\n=== Statistical Summary ===\n');
fprintf('Format: DC Profile | Channel | Time | p-value (Charging/Discharging) | Test | Effect Size | 0cyc (mean±std, n) | 200cyc (mean±std, n)\n');
fprintf('--------------------------------------------------------------------------------\n');

for dt_idx = 1:length(dt_list)
    dt_sec = dt_list(dt_idx);
    field_name = sprintf('DCIR_%ds', dt_sec);
    
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        
        for ch_idx = 1:length(common_channels)
            channel = common_channels{ch_idx};
            
            if isfield(comparison_results, field_name) && isfield(comparison_results.(field_name), dc_profile) && isfield(comparison_results.(field_name).(dc_profile), channel)
                result = comparison_results.(field_name).(dc_profile).(channel);
                
                % Charging results
                if isfield(result, 'charging')
                    significance_chg = '';
                    if result.charging.p_value < 0.001
                        significance_chg = '***';
                    elseif result.charging.p_value < 0.01
                        significance_chg = '**';
                    elseif result.charging.p_value < 0.05
                        significance_chg = '*';
                    end
                    
                    fprintf('%s | %s | %ds | Charging: %.4f%s | %s | %.3f | %.2f±%.2f (%d) | %.2f±%.2f (%d)\n', ...
                        dc_profile, channel, dt_sec, result.charging.p_value, significance_chg, ...
                        result.charging.test_type, result.charging.cohens_d, ...
                        result.charging.mean_0cyc, result.charging.std_0cyc, result.charging.n_0cyc, ...
                        result.charging.mean_200cyc, result.charging.std_200cyc, result.charging.n_200cyc);
                end
                
                % Discharging results
                if isfield(result, 'discharging')
                    significance_dchg = '';
                    if result.discharging.p_value < 0.001
                        significance_dchg = '***';
                    elseif result.discharging.p_value < 0.01
                        significance_dchg = '**';
                    elseif result.discharging.p_value < 0.05
                        significance_dchg = '*';
                    end
                    
                    fprintf('%s | %s | %ds | Discharging: %.4f%s | %s | %.3f | %.2f±%.2f (%d) | %.2f±%.2f (%d)\n', ...
                        dc_profile, channel, dt_sec, result.discharging.p_value, significance_dchg, ...
                        result.discharging.test_type, result.discharging.cohens_d, ...
                        result.discharging.mean_0cyc, result.discharging.std_0cyc, result.discharging.n_0cyc, ...
                        result.discharging.mean_200cyc, result.discharging.std_200cyc, result.discharging.n_200cyc);
                end
            end
        end
    end
end

%% Create visualization
% Create directory if it doesn't exist
if ~exist('figures/DC_Events', 'dir')
    mkdir('figures/DC_Events');
    fprintf('Created directory: figures/DC_Events\n');
end

% Create separate figures for each channel
for ch_idx = 1:length(common_channels)
    channel = common_channels{ch_idx};
    
    % Create separate figures for each time interval
    for dt_idx = 1:length(dt_list)
        dt_sec = dt_list(dt_idx);
        field_name = sprintf('DCIR_%ds', dt_sec);
        
        figure('Name', sprintf('DCIR Comparison: %s %ds (0cyc vs 200cyc)', channel, dt_sec), 'Position', [100, 100, 1200, 800]);
        
        for dc_idx = 1:length(dc_profiles)
            dc_profile = dc_profiles{dc_idx};
            
            subplot(2, 4, dc_idx);
            
            % Collect data for this channel, DC profile, and specific time interval
            dcir_0cyc_chg = [];
            dcir_200cyc_chg = [];
            dcir_0cyc_dchg = [];
            dcir_200cyc_dchg = [];
            
            % 0cyc charging data for this specific time interval
            chg_struct_name_0cyc = sprintf('%s_Drive_0cyc_ChgEvent', channel);
            if isfield(data_0cyc, chg_struct_name_0cyc) && isfield(data_0cyc.(chg_struct_name_0cyc), soc_level) && isfield(data_0cyc.(chg_struct_name_0cyc).(soc_level), dc_profile)
                events = fieldnames(data_0cyc.(chg_struct_name_0cyc).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = data_0cyc.(chg_struct_name_0cyc).(soc_level).(dc_profile).(evt_name);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_0cyc_chg = [dcir_0cyc_chg; evt_data.(field_name).val];
                    end
                end
            end
            
            % 200cyc charging data for this specific time interval
            chg_struct_name_200cyc = sprintf('%s_Drive_200cyc_ChgEvent', channel);
            if isfield(data_200cyc, chg_struct_name_200cyc) && isfield(data_200cyc.(chg_struct_name_200cyc), soc_level) && isfield(data_200cyc.(chg_struct_name_200cyc).(soc_level), dc_profile)
                events = fieldnames(data_200cyc.(chg_struct_name_200cyc).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = data_200cyc.(chg_struct_name_200cyc).(soc_level).(dc_profile).(evt_name);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_200cyc_chg = [dcir_200cyc_chg; evt_data.(field_name).val];
                    end
                end
            end
            
            % 0cyc discharging data for this specific time interval
            dchg_struct_name_0cyc = sprintf('%s_Drive_0cyc_DchEvent', channel);
            if isfield(data_0cyc, dchg_struct_name_0cyc) && isfield(data_0cyc.(dchg_struct_name_0cyc), soc_level) && isfield(data_0cyc.(dchg_struct_name_0cyc).(soc_level), dc_profile)
                events = fieldnames(data_0cyc.(dchg_struct_name_0cyc).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = data_0cyc.(dchg_struct_name_0cyc).(soc_level).(dc_profile).(evt_name);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_0cyc_dchg = [dcir_0cyc_dchg; evt_data.(field_name).val];
                    end
                end
            end
            
            % 200cyc discharging data for this specific time interval
            dchg_struct_name_200cyc = sprintf('%s_Drive_200cyc_DchEvent', channel);
            if isfield(data_200cyc, dchg_struct_name_200cyc) && isfield(data_200cyc.(dchg_struct_name_200cyc), soc_level) && isfield(data_200cyc.(dchg_struct_name_200cyc).(soc_level), dc_profile)
                events = fieldnames(data_200cyc.(dchg_struct_name_200cyc).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = data_200cyc.(dchg_struct_name_200cyc).(soc_level).(dc_profile).(evt_name);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_200cyc_dchg = [dcir_200cyc_dchg; evt_data.(field_name).val];
                    end
                end
            end
            
            % Create scatter plot for this specific time interval
            all_data = [];
            group_labels = {};
            
            % Add charging data
            if ~isempty(dcir_0cyc_chg)
                all_data = [all_data; dcir_0cyc_chg];
                group_labels = [group_labels; repmat({'0cyc-Charging'}, length(dcir_0cyc_chg), 1)];
            end
            
            if ~isempty(dcir_200cyc_chg)
                all_data = [all_data; dcir_200cyc_chg];
                group_labels = [group_labels; repmat({'200cyc-Charging'}, length(dcir_200cyc_chg), 1)];
            end
            
            % Add discharging data
            if ~isempty(dcir_0cyc_dchg)
                all_data = [all_data; dcir_0cyc_dchg];
                group_labels = [group_labels; repmat({'0cyc-Discharging'}, length(dcir_0cyc_dchg), 1)];
            end
            
            if ~isempty(dcir_200cyc_dchg)
                all_data = [all_data; dcir_200cyc_dchg];
                group_labels = [group_labels; repmat({'200cyc-Discharging'}, length(dcir_200cyc_dchg), 1)];
            end
            
            if ~isempty(all_data)
                % Create scatter plot with different groups on x-axis
                hold on;
                
                % Define x positions for each group
                x_positions = [1, 2, 3, 4];  % 0cyc-Charging, 200cyc-Charging, 0cyc-Discharging, 200cyc-Discharging
                group_names = {'0cyc-Charging', '200cyc-Charging', '0cyc-Discharging', '200cyc-Discharging'};
                colors = {'blue', 'red', 'green', 'magenta'};
                
                % Plot data points for each group
                if ~isempty(dcir_0cyc_chg)
                    scatter(ones(size(dcir_0cyc_chg)), dcir_0cyc_chg, 'filled', 'MarkerFaceColor', colors{1}, 'DisplayName', '0cyc-Charging');
                end
                
                if ~isempty(dcir_200cyc_chg)
                    scatter(2*ones(size(dcir_200cyc_chg)), dcir_200cyc_chg, 'filled', 'MarkerFaceColor', colors{2}, 'DisplayName', '200cyc-Charging');
                end
                
                if ~isempty(dcir_0cyc_dchg)
                    scatter(3*ones(size(dcir_0cyc_dchg)), dcir_0cyc_dchg, 'filled', 'MarkerFaceColor', colors{3}, 'DisplayName', '0cyc-Discharging');
                end
                
                if ~isempty(dcir_200cyc_dchg)
                    scatter(4*ones(size(dcir_200cyc_dchg)), dcir_200cyc_dchg, 'filled', 'MarkerFaceColor', colors{4}, 'DisplayName', '200cyc-Discharging');
                end
                
                hold off;
                
                title(sprintf('%s', dc_profile));
                ylabel('DCIR (mΩ)');
                xlabel('Groups');
                xticks(x_positions);
                xticklabels(group_names);
                xtickangle(45);
                grid on;
                
                % Create legend with statistical information
                legend_entries = {};
                
                % Charging comparison info
                if ~isempty(dcir_0cyc_chg) && ~isempty(dcir_200cyc_chg)
                    if isfield(comparison_results, field_name) && isfield(comparison_results.(field_name), dc_profile) && isfield(comparison_results.(field_name).(dc_profile), channel) && isfield(comparison_results.(field_name).(dc_profile).(channel), 'charging')
                        result = comparison_results.(field_name).(dc_profile).(channel).charging;
                        significance = '';
                        if result.p_value < 0.001
                            significance = '***';
                        elseif result.p_value < 0.01
                            significance = '**';
                        elseif result.p_value < 0.05
                            significance = '*';
                        end
                        legend_entries{end+1} = sprintf('Charging: p=%.4f%s, 0cyc=%.2f±%.2f, 200cyc=%.2f±%.2f', ...
                            result.p_value, significance, result.mean_0cyc, result.std_0cyc, result.mean_200cyc, result.std_200cyc);
                    end
                end
                
                % Discharging comparison info
                if ~isempty(dcir_0cyc_dchg) && ~isempty(dcir_200cyc_dchg)
                    if isfield(comparison_results, field_name) && isfield(comparison_results.(field_name), dc_profile) && isfield(comparison_results.(field_name).(dc_profile), channel) && isfield(comparison_results.(field_name).(dc_profile).(channel), 'discharging')
                        result = comparison_results.(field_name).(dc_profile).(channel).discharging;
                        significance = '';
                        if result.p_value < 0.001
                            significance = '***';
                        elseif result.p_value < 0.01
                            significance = '**';
                        elseif result.p_value < 0.05
                            significance = '*';
                        end
                        legend_entries{end+1} = sprintf('Discharging: p=%.4f%s, 0cyc=%.2f±%.2f, 200cyc=%.2f±%.2f', ...
                            result.p_value, significance, result.mean_0cyc, result.std_0cyc, result.mean_200cyc, result.std_200cyc);
                    end
                end
                
                % Add legend if we have statistical information
                if ~isempty(legend_entries)
                    legend(legend_entries, 'Location', 'best', 'FontSize', 8);
                end
            else
                title(sprintf('%s (No data)', dc_profile));
            end
        end
        
        sgtitle(sprintf('DCIR Comparison: %s %ds (SOC50)', channel, dt_sec));
        
        % Save figure
        saveas(gcf, sprintf('figures/DC_Events/%s_%ds_DCIR_Comparison.fig', channel, dt_sec));
        saveas(gcf, sprintf('figures/DC_Events/%s_%ds_DCIR_Comparison.png', channel, dt_sec));
        close(gcf);
    end
end

%% Save results
save('DCIR_Statistical_Comparison_Results.mat', 'comparison_results');
fprintf('\nResults saved to: DCIR_Statistical_Comparison_Results.mat\n');

fprintf('\nStatistical comparison completed!\n');
fprintf('Legend: * p<0.05, ** p<0.01, *** p<0.001\n'); 