clear; clc; close all;

%% Configuration
base_path = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Experiment_DriveCycle_Fitting\Results';
channels = {'ch9', 'ch10', 'ch11', 'ch12', 'ch13', 'ch14', 'ch15', 'ch16'};
cycle_conditions = {'cyc0', 'cyc200'};
soc_conditions = {'SOC50', 'SOC70', 'SOC90'};
dc_names = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

%% Load and Compare Results
fprintf('=== Channel Comparison Analysis ===\n');

for current_cycle = cycle_conditions
    current_cycle = current_cycle{1};
    fprintf('\n=== Processing %s Condition ===\n', current_cycle);
    
    for current_soc = soc_conditions
        current_soc = current_soc{1};
        fprintf('\n--- Processing %s ---\n', current_soc);
        
        for dc_idx = 1:length(dc_names)
            DC_name = dc_names{dc_idx};
            fprintf('Processing %s...\n', DC_name);
            
            % Initialize results table for this DC
            channel_results = table();
            
            % Load results from each channel
            for ch_idx = 1:length(channels)
                current_channel = channels{ch_idx};
                
                % Construct file path
                result_file = fullfile(base_path, current_cycle, current_channel, current_soc, DC_name, ...
                    sprintf('%s_%s_%s_result.mat', current_cycle, current_soc, DC_name));
                
                % Check if file exists
                if exist(result_file, 'file')
                    % Load the result
                    load(result_file);
                    
                    % Extract parameters
                    R0 = params_optimal(1) * 1000;  % Convert to mΩ
                    R1 = params_optimal(2) * 1000;
                    R2 = params_optimal(3) * 1000;
                    tau1 = params_optimal(4);
                    tau2 = params_optimal(5);
                    C1 = (tau1 / params_optimal(2)) / 1000;  % Convert to kF
                    C2 = (tau2 / params_optimal(3)) / 1000;
                    
                    % Calculate performance metrics from cost function
                    rmse = cost_optimal * 1000;  % Convert to mV
                    
                    % Add row to table
                    new_row = table({current_channel}, R0, R1, R2, tau1, tau2, C1, C2, rmse, ...
                        'VariableNames', {'Channel', 'R0_mOhm', 'R1_mOhm', 'R2_mOhm', ...
                        'tau1_s', 'tau2_s', 'C1_kF', 'C2_kF', 'RMSE_mV'});
                    
                    channel_results = [channel_results; new_row];
                    
                else
                    fprintf('  Warning: File not found for %s\n', result_file);
                end
            end
            
            % Save comparison table
            if height(channel_results) > 0
                % Create output directory
                output_dir = fullfile(base_path, current_cycle, sprintf('%s_%s_comparison', current_soc, DC_name));
                if ~exist(output_dir, 'dir')
                    mkdir(output_dir);
                end
                
                % Save table and statistics
                output_file = fullfile(output_dir, sprintf('%s_%s_%s_channel_comparison.mat', current_cycle, current_soc, DC_name));
                save(output_file, 'channel_results', 'statistics');
                
                % Calculate statistics
                param_names = {'R0_mOhm', 'R1_mOhm', 'R2_mOhm', 'tau1_s', 'tau2_s', 'C1_kF', 'C2_kF', 'RMSE_mV'};
                statistics = struct();
                
                for p = 1:length(param_names)
                    param = param_names{p};
                    values = channel_results.(param);
                    statistics.([param '_avg']) = mean(values);
                    statistics.([param '_min']) = min(values);
                    statistics.([param '_max']) = max(values);
                    statistics.([param '_std']) = std(values);
                end
                
                % Display summary
                fprintf('  Results for %d channels loaded\n', height(channel_results));
                fprintf('  Mean RMSE: %.2f ± %.2f mV\n', statistics.RMSE_mV_avg, statistics.RMSE_mV_std);
                fprintf('  RMSE range: %.2f - %.2f mV\n', statistics.RMSE_mV_min, statistics.RMSE_mV_max);
                
                % Display table
                fprintf('\n  Channel Comparison Table:\n');
                disp(channel_results);
                
                % Display statistics
                fprintf('\n  Parameter Statistics:\n');
                for p = 1:length(param_names)
                    param = param_names{p};
                    fprintf('  %s: %.3f ± %.3f (%.3f - %.3f)\n', param, ...
                        statistics.([param '_avg']), statistics.([param '_std']), ...
                        statistics.([param '_min']), statistics.([param '_max']));
                end
                
            else
                fprintf('  No results found for %s\n', DC_name);
            end
        end
    end
end

fprintf('\n=== Analysis Complete ===\n'); 