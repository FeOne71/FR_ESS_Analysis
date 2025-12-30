%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECM 2RC Fitting Script Based on Drive Cycle Data
%
% Purpose: Fit ECM 2RC model to Channel 9 drive cycle data
% Parameters: R0, R1, R2, tau1, tau2
% SOC Calculation: Docker method (two-point OCV)
%   - SOC1: Calculated from OCV at the end of first rest period
%   - SOC2: Calculated from OCV at the end of second rest period
%   - Intermediate period: Linear interpolation based on current integration ratio
% Cycle-specific OCV_integrated_{cycle_num} is used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. Configuration and Data Loading
% =========================================================================
fprintf('=== ECM 2RC Fitting Script Started ===\n');

% Data path configuration
ocv_data_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat';
drive_cycle_folder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';

% Load OCV data
load(ocv_data_path);
if ~exist('OCV_data', 'var'), error("'OCV_data' variable not found."); end

% Analysis parameters (test: SOC50, DC1 only)
TARGET_CHANNEL_PREFIX = 'ch9_Drive_';
TARGET_SOC_LEVELS = {'SOC50'};  % Test: SOC50 only
TARGET_DC_PROFILES = {'DC1'};   % Test: DC1 only
TARGET_CYCLES = [0];%, 200, 400, 600];

% Extract channel name and number from prefix (e.g., 'ch16_Drive_' -> 'ch16', channel_num = 16)
channel_match = regexp(TARGET_CHANNEL_PREFIX, 'ch(\d+)_', 'tokens');
if ~isempty(channel_match)
    channel_name = ['ch' channel_match{1}{1}];
    channel_num = str2double(channel_match{1}{1});  % Channel number for indexing individual_static_capacity
else
    channel_name = 'ch_unknown';
    channel_num = 1;  % Default to channel 1 if extraction fails
end

fprintf('OCV data loaded: %s\n', ocv_data_path);
fprintf('Drive cycle data folder: %s\n', drive_cycle_folder);

%% 2. ECM Fitting Configuration
% =========================================================================
fprintf('\n=== ECM Fitting Configuration ===\n');

% Initial parameter estimation [R0, R1, R2, tau1, tau2] (BSL Reference)
params_initial = [0.003, 0.005, 0.0005, 54, 3228];  % [Ω, Ω, Ω, s, s]

% Parameter bounds (BSL Reference)
lb = [0, 0, 0, 0.01, 0.01];   % Lower bound
ub = [0.05, 0.005, 0.03, 100, 5000];   % Upper bound

% Linear constraint: tau1 < tau2 (A_ineq * x <= b_ineq)
% [0, 0, 0, 1, -1] * [R0, R1, R2, tau1, tau2]' <= 0
% i.e., tau1 - tau2 <= 0  =>  tau1 <= tau2
A_ineq = [0, 0, 0, 1, -1];  % tau1 < tau2 constraint
b_ineq = 0;

% Optimization options (BSL Reference: optimset 방식)
options = optimset('display', 'off', ...
    'MaxIter', 1e3, ...
    'MaxFunEvals', 1e4, ...
    'TolFun', 1e-14, ...
    'TolX', 1e-15);

fprintf('Initial parameters: R0=%.6f, R1=%.6f, R2=%.6f, tau1=%.2f, tau2=%.2f\n', ...
    params_initial(1), params_initial(2), params_initial(3), params_initial(4), params_initial(5));

%% 3. Drive Cycle File List Generation
% =========================================================================
fprintf('\n=== Drive Cycle File List Generation ===\n');

drive_cycle_files = dir(fullfile(drive_cycle_folder, 'parsedDriveCycle_*_filtered.mat'));
fprintf('Found drive cycle files: %d\n', length(drive_cycle_files));

% Extract cycle numbers from files
file_cycles = zeros(length(drive_cycle_files), 1);
for i = 1:length(drive_cycle_files)
    filename = drive_cycle_files(i).name;
    match = regexp(filename, 'parsedDriveCycle_(\w+)_filtered.mat', 'tokens');
    cycle_str = match{1}{1};
    file_cycles(i) = str2double(regexp(cycle_str, '\d+', 'match'));
end

fprintf('Cycle numbers: %s\n', mat2str(unique(file_cycles)));

%% 4. ECM Fitting Results Storage Structure
% =========================================================================
ecm_results = struct();
fitting_summary = [];

%% 5. Drive Cycle Data Processing and ECM Fitting
% =========================================================================
fprintf('\n=== Drive Cycle Data Processing and ECM Fitting Started ===\n');

for soc_name_cell = TARGET_SOC_LEVELS
    soc_name = soc_name_cell{1};
    target_soc = str2double(soc_name(4:end)); % SOC90 -> 90
    
    fprintf('\n--- Processing %s... ---\n', soc_name);
    
    for dc_name_cell = TARGET_DC_PROFILES
        dc_name = dc_name_cell{1};
        
        for cycle_num = TARGET_CYCLES
            fprintf('Processing: %s-%s-Cycle%d\n', soc_name, dc_name, cycle_num);
            
            % Select cycle-specific OCV data
            ocv_field_name = sprintf('OCV_integrated_%d', cycle_num);
            if isfield(OCV_data, ocv_field_name)
                ocv_integrated_data = OCV_data.(ocv_field_name);
                fprintf('  Using %s OCV data for cycle %d\n', ocv_field_name, cycle_num);
                
                % Extract capacity value for this channel from individual_static_capacity
                % Note: individual_static_capacity contains 8 channels (ch9~ch16), indexed as 1~8
                % Therefore: channel_num 9 -> index 1, channel_num 10 -> index 2, ..., channel_num 16 -> index 8
                if isfield(ocv_integrated_data, 'individual_static_capacity')
                    capacity_array = ocv_integrated_data.individual_static_capacity;
                    capacity_length = length(capacity_array);
                    
                    % Convert channel number to array index (ch9~ch16 -> index 1~8)
                    if channel_num >= 9 && channel_num <= 16
                        capacity_index = channel_num - 8;  % ch9 -> 1, ch10 -> 2, ..., ch16 -> 8
                        if capacity_index >= 1 && capacity_index <= capacity_length
                            channel_capacity = capacity_array(capacity_index);
                            fprintf('  Channel capacity: %.4f Ah (ch%d -> index %d)\n', channel_capacity, channel_num, capacity_index);
                        else
                            channel_capacity = NaN;
                            fprintf('  Warning: Calculated index %d is out of range [1, %d]\n', capacity_index, capacity_length);
                        end
                    else
                        channel_capacity = NaN;
                        fprintf('  Warning: Channel number %d is not in expected range [9, 16]\n', channel_num);
                        fprintf('  Available capacity array (length %d): %s\n', capacity_length, mat2str(capacity_array));
                    end
                else
                    channel_capacity = NaN;
                    fprintf('  Warning: individual_static_capacity field not found in OCV_integrated_data\n');
                end
            else
                channel_capacity = NaN;
            end
            
            % Find cycle file
            cycle_file_idx = find(file_cycles == cycle_num);
            if isempty(cycle_file_idx)
                fprintf('  Cycle file not found\n');
                continue;
            end
            
            % Load file
            current_filepath = fullfile(drive_cycle_folder, drive_cycle_files(cycle_file_idx).name);
            loaded_data = load(current_filepath);
            
            % Generate variable name
            filename = drive_cycle_files(cycle_file_idx).name;
            match = regexp(filename, 'parsedDriveCycle_(\w+)_filtered.mat', 'tokens');
            cycle_str = match{1}{1};
            expected_var_name = ['parsedDriveCycle_', cycle_str];
            
            if ~isfield(loaded_data, expected_var_name)
                fprintf('  Expected variable not found\n');
                continue;
            end
            
            drive_cycle_data_struct = loaded_data.(expected_var_name);
            target_channel_name = [TARGET_CHANNEL_PREFIX, cycle_str];
            
            if ~isfield(drive_cycle_data_struct, target_channel_name)
                fprintf('  Target channel not found\n');
                continue;
            end
            
            dc_data_all = drive_cycle_data_struct.(target_channel_name);
            
            % Check SOC/DC data
            if ~isfield(dc_data_all, soc_name) || ~isfield(dc_data_all.(soc_name), dc_name)
                fprintf('  SOC/DC data not found\n');
                continue;
            end
            
            dc_data = dc_data_all.(soc_name).(dc_name);
            
            if isempty(dc_data.V) || isempty(dc_data.I)
                fprintf('  Empty voltage/current data\n');
                continue;
            end
            
            % Extract data (full dataset)
            V_full = dc_data.V;
            I_full = dc_data.I;
            t_full = seconds(dc_data.t);  % Convert duration to seconds
            
            % Check actual sampling interval from loaded data
            if length(t_full) > 1
                dt_sample = diff(t_full);
                mean_dt = mean(dt_sample);
                fprintf('  Actual sampling interval from data: mean=%.4f seconds\n', mean_dt);
                if abs(mean_dt - 0.1) < 0.01
                    fprintf('    -> 0.1 second sampling confirmed\n');
                elseif abs(mean_dt - 1.0) < 0.01
                    fprintf('    -> 1.0 second sampling confirmed\n');
                else
                    fprintf('    -> Non-standard sampling interval detected\n');
                end
            else
                mean_dt = 0.1;  % Default assumption if cannot determine
                fprintf('  Warning: Cannot determine sampling interval, assuming 0.1s\n');
            end
            
            % Find rest periods (from full data)
            rest_mask = abs(I_full) <= 3.2;  % 64*0.05 = 3.2A
            
            % Minimum rest time in seconds (e.g., 10 minutes = 600 seconds)
            min_rest_time_seconds = 590;  % 10 minutes in seconds
            % Convert to number of points based on actual sampling interval
            min_rest_time_points = round(min_rest_time_seconds / mean_dt);
            fprintf('  Minimum rest time: %.1f minutes (%d points at %.4f s interval)\n', ...
                min_rest_time_seconds/60, min_rest_time_points, mean_dt);
            rest_periods = [];
            in_rest = false;
            rest_start = 0;
            
            for i = 1:length(rest_mask)
                if rest_mask(i) && ~in_rest
                    rest_start = i;
                    in_rest = true;
                elseif ~rest_mask(i) && in_rest
                    rest_duration_points = (i-1) - rest_start + 1;
                    if rest_duration_points >= min_rest_time_points
                        rest_periods = [rest_periods; rest_start, i-1];
                    end
                    in_rest = false;
                end
            end
            
            % Check last period
            if in_rest
                rest_duration_points = length(rest_mask) - rest_start + 1;
                if rest_duration_points >= min_rest_time_points
                    rest_periods = [rest_periods; rest_start, length(rest_mask)];
                end
            end
            
            % Rest period debugging
            fprintf('  Rest period detection results:\n');
            fprintf('    Number of rest periods found: %d\n', size(rest_periods, 1));
            min_rest_time_minutes = min_rest_time_seconds / 60;
            fprintf('    Minimum rest time: %.1f minutes (%d points at %.4f s interval)\n', ...
                min_rest_time_minutes, min_rest_time_points, mean_dt);
            fprintf('    Rest mask ratio: %.1f%%\n', sum(rest_mask)/length(rest_mask)*100);
            
            for j = 1:size(rest_periods, 1)
                duration_points = rest_periods(j, 2) - rest_periods(j, 1) + 1;
                % Calculate time using actual time intervals
                if rest_periods(j, 2) <= length(t_full)
                    duration_seconds = t_full(rest_periods(j, 2)) - t_full(rest_periods(j, 1));
                    duration_minutes = duration_seconds / 60;
                else
                    % Use average interval if index is out of range
                    duration_minutes = duration_points * mean_dt / 60;
                end
                fprintf('    Rest period %d: indices %d~%d (%.1f minutes)\n', j, rest_periods(j, 1), rest_periods(j, 2), duration_minutes);
            end
            
            % --- Data preprocessing: Determine fitting range and trim data ---
            if size(rest_periods, 1) < 2
                fprintf('  WARNING: Could not find two or more rest periods. Skipping this cycle.\n');
                continue; % Cannot fit without second rest period
            end
            
            % Define end point of data for fitting as 'end of second rest period'
            end_idx = rest_periods(2, 2);
            
            % Trim all data based on this end point
            V_measured = V_full(1:end_idx);
            I_measured = I_full(1:end_idx);
            t_measured = t_full(1:end_idx);
            
            fprintf('  Data preprocessing: Using %d out of %d original data points (excluding final SOC recovery period)\n', end_idx, length(V_full));
            
            % Create rest period information for trimmed data (use only first two for fitting)
            % This information maintains the absolute indices from the original data
            rest_periods_fit = rest_periods(1:2, :);
            
            % --- End of data preprocessing ---
            
            % All calculations from now on use trimmed V_measured, I_measured, t_measured
            
            % SOC calculation (using trimmed data and rest period information)
            [~, ~, ~, SOC_full, ~] = ...
                calculate_soc_docker(V_measured, I_measured, t_measured, ocv_integrated_data, rest_periods_fit);
            fprintf('  SOC calculation completed\n');
            
            % Generate OCV vector from SOC profile using existing OCV_SOC_func
            fprintf('  Generating OCV vector from SOC profile...\n');
            if isfield(ocv_integrated_data, 'OCV_SOC_func')
                % Use existing OCV_SOC_func from OCV_integrated structure
                OCV_vec = ocv_integrated_data.OCV_SOC_func(SOC_full);
            else
                % Fallback: create interpolation function if OCV_SOC_func doesn't exist
                fprintf('  Warning: OCV_SOC_func not found. Creating interpolation function...\n');
                soc_grid_points = ocv_integrated_data.SOC_grid;
                ocv_values_at_soc_grid = ocv_integrated_data.V_avg_SOC;
                [ocv_values_sorted, uniqueIdx] = unique(ocv_values_at_soc_grid, 'stable');
                soc_grid_sorted = soc_grid_points(uniqueIdx);
                soc_to_ocv_interp_func = @(soc) interp1(soc_grid_sorted, ocv_values_sorted, soc, 'linear', 'extrap');
                OCV_vec = soc_to_ocv_interp_func(SOC_full);
            end
            fprintf('  OCV vector generated (length: %d)\n', length(OCV_vec));
            
            % Perform ECM fitting
            % Create cost function
            cost_handle = @(params) ecm_cost_function(params, I_measured, V_measured, t_measured, OCV_vec);
            
            % MultiStart optimization
            ms = MultiStart('Display', 'off', 'UseParallel', true);
            nStartPts = 30;
            startPts = RandomStartPointSet('NumStartPoints', nStartPts);
            
            % Define optimization problem (including tau1 < tau2 constraint)
            problem = createOptimProblem('fmincon', ...
                'objective', cost_handle, ...
                'x0', params_initial, ...
                'lb', lb, ...
                'ub', ub, ...
                'Aineq', A_ineq, ...
                'bineq', b_ineq, ...
                'options', options);
            
            % Execute optimization
            tic;
            [params_optimal, cost_optimal, exitflag, output, solutions] = run(ms, problem, startPts);
            optimization_time = toc;
            
            % Model validation
            V_model = ecm_2rc_model(params_optimal, I_measured, t_measured, OCV_vec);
            residual = V_measured - V_model;
            rmse = sqrt(mean(residual.^2));
            mae = mean(abs(residual));
            max_error = max(abs(residual));
            
            % Calculate capacitances
            C1 = params_optimal(4) / params_optimal(2);
            C2 = params_optimal(5) / params_optimal(3);
            
            % Save results
            result_key = sprintf('%s_%s_Cycle%d', soc_name, dc_name, cycle_num);
            ecm_results.(result_key) = struct();
            ecm_results.(result_key).SOC_target = target_soc;
            ecm_results.(result_key).DC_profile = dc_name;
            ecm_results.(result_key).Cycle = cycle_num;
            ecm_results.(result_key).R0 = params_optimal(1);
            ecm_results.(result_key).R1 = params_optimal(2);
            ecm_results.(result_key).R2 = params_optimal(3);
            ecm_results.(result_key).tau1 = params_optimal(4);
            ecm_results.(result_key).tau2 = params_optimal(5);
            ecm_results.(result_key).C1 = C1;
            ecm_results.(result_key).C2 = C2;
            ecm_results.(result_key).rmse = rmse;
            ecm_results.(result_key).mae = mae;
            ecm_results.(result_key).max_error = max_error;
            ecm_results.(result_key).optimization_time = optimization_time;
            ecm_results.(result_key).exitflag = exitflag;
            ecm_results.(result_key).V_model = V_model;
            ecm_results.(result_key).V_measured = V_measured;
            ecm_results.(result_key).I_measured = I_measured;
            ecm_results.(result_key).SOC = SOC_full;
            ecm_results.(result_key).t_measured = t_measured;
            ecm_results.(result_key).Capacity_Ah = channel_capacity;  % Store capacity value
            
            % Add to summary table (include capacity)
            summary_row = {target_soc, dc_name, cycle_num, channel_capacity, params_optimal(1), params_optimal(2), params_optimal(3), params_optimal(4), params_optimal(5), C1, C2, rmse, mae, max_error, optimization_time};
            fitting_summary = [fitting_summary; summary_row];
            
            fprintf('  ECM fitting completed: RMSE=%.4f V\n', rmse);
            
            % Fitting visualization
            fprintf('  Generating fitting visualization...\n');
            create_fitting_visualization(channel_name, soc_name, dc_name, cycle_num, V_measured, V_model, I_measured, SOC_full, t_measured, rmse, mae, max_error, channel_capacity);
        end
    end
end

%% 6. Results Summary and Saving
% =========================================================================

if ~isempty(fitting_summary)
    % Create table
    ecm_summary_table = cell2table(fitting_summary, 'VariableNames', ...
        {'SOC_target', 'DC_Profile', 'Cycle', 'Capacity_Ah', 'R0', 'R1', 'R2', 'tau1', 'tau2', 'C1', 'C2', 'RMSE', 'MAE', 'Max_Error', 'Time_s'});
    
    % Sort
    ecm_summary_table = sortrows(ecm_summary_table, {'SOC_target', 'DC_Profile', 'Cycle'});
    
    % Output results
    fprintf('\nECM 2RC Fitting Results Summary:\n');
    disp(ecm_summary_table);
    
    % Cycle-by-cycle statistics summary (table format)
    fprintf('\n=== Cycle-by-Cycle Statistics Summary ===\n');
    fprintf('Total processed combinations: %d\n', height(ecm_summary_table));
    
    % Cycle-by-cycle analysis
    unique_cycles = unique(ecm_summary_table.Cycle);
    
    % Table header (with capacity)
    fprintf('\nCycle | Capacity (Ah) | R0 (mΩ)     | R1 (mΩ)     | R2 (mΩ)     | tau1 (s)    | tau2 (s)    | RMSE (mV)   | Time (s)\n');
    fprintf('------|---------------|-------------|-------------|-------------|-------------|-------------|-------------|----------\n');
    
    for i = 1:length(unique_cycles)
        cycle_num = unique_cycles(i);
        cycle_mask = ecm_summary_table.Cycle == cycle_num;
        cycle_data = ecm_summary_table(cycle_mask, :);
        
        % Calculate values (convert to mΩ) - currently one result per cycle
        capacity_val = cycle_data.Capacity_Ah(1);
        r0_val = cycle_data.R0(1) * 1000;
        r1_val = cycle_data.R1(1) * 1000;
        r2_val = cycle_data.R2(1) * 1000;
        tau1_val = cycle_data.tau1(1);
        tau2_val = cycle_data.tau2(1);
        rmse_val = cycle_data.RMSE(1) * 1000;
        time_val = cycle_data.Time_s(1);
        
        if ~isnan(capacity_val)
            fprintf('%5d | %13.4f | %11.3f | %11.3f | %11.3f | %11.3f | %11.3f | %11.3f | %8.2f\n', ...
                cycle_num, capacity_val, r0_val, r1_val, r2_val, tau1_val, tau2_val, rmse_val, time_val);
        else
            fprintf('%5d | %13s | %11.3f | %11.3f | %11.3f | %11.3f | %11.3f | %11.3f | %8.2f\n', ...
                cycle_num, 'N/A', r0_val, r1_val, r2_val, tau1_val, tau2_val, rmse_val, time_val);
        end
    end
    
    % Save results (channel-specific folder)
    save_dir = fullfile('Results', 'ECM_Fitting', channel_name);
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    
    % Save MAT file
    mat_file = fullfile(save_dir, sprintf('ECM_2RC_Fitting_Results_%s.mat', channel_name));
    save(mat_file, 'ecm_results', 'ecm_summary_table', 'ocv_integrated_data', '-v7.3');
    fprintf('\nResults saved: %s\n', mat_file);
    
    % Generate overall error visualization
    fprintf('\n=== Generating Overall Error Visualization ===\n');
    create_overall_error_visualization(channel_name, ecm_results);
    
else
    fprintf('No processed data found.\n');
end

fprintf('\n=== ECM 2RC Fitting Completed ===\n');

%% Helper Functions
% =========================================================================

function create_fitting_visualization(channel_name, soc_name, dc_name, cycle_num, V_measured, V_model, I_measured, SOC_full, t_measured, rmse, mae, max_error, channel_capacity)
    % ECM fitting visualization function
    % channel_capacity: Capacity value in Ah (from individual_static_capacity)
    
    % Convert time to hours
    t_hours = t_measured / 3600;
    
    % Create new figure (vertical layout, larger size for better visibility)
    figure('Position', [100, 100, 1600, 1600]);
    
    % 1. Voltage comparison and SOC (same subplot)
    subplot(3, 1, 1);
    set(gca, 'Position', [0.1, 0.68, 0.85, 0.28]);  % [left, bottom, width, height] - larger subplot
    yyaxis left;
    plot(t_hours, V_measured, '-', 'Color', '#0073C2', 'LineWidth', 1.5, 'DisplayName', 'Measured Voltage');
    hold on;
    plot(t_hours, V_model, '-', 'Color', '#CD534C', 'LineWidth', 1.5, 'DisplayName', 'ECM Model');
    ylabel('Voltage (V)', 'FontWeight', 'bold');
    ylim([3.6, 4.2]);  % Fixed voltage range
    
    yyaxis right;
    plot(t_hours, SOC_full, 'g-', 'LineWidth', 2, 'DisplayName', 'SOC');
    ylabel('SOC (%)', 'FontWeight', 'bold');
    ylim([50, 100]);  % Fixed SOC range
    
    xlabel('Time (h)', 'FontWeight', 'bold');
    title(sprintf('Voltage & SOC\nRMSE: %.4f V', rmse), 'FontWeight', 'bold', 'FontSize', 10);
    legend('Measured Voltage', 'ECM Model', 'SOC', 'Location', 'best');
    grid on;
    
    % Set axis colors to black
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.YAxis(1).Color = 'k';  % Left y-axis
    ax.YAxis(2).Color = 'k';  % Right y-axis
    
    % 2. Model error analysis
    subplot(3, 1, 2);
    set(gca, 'Position', [0.1, 0.36, 0.85, 0.28]);  % Larger subplot
    residual = V_measured - V_model;
    plot(t_hours, residual * 1000, 'k-', 'LineWidth', 1);
    xlabel('Time (h)', 'FontWeight', 'bold');
    ylabel('Voltage Error (mV)', 'FontWeight', 'bold');
    title(sprintf('Model Error\nMAE: %.4f V, Max Error: %.4f V', mae, max_error), 'FontWeight', 'bold', 'FontSize', 10);
    grid on;
    
    % Set axis colors to black
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    
    % 3. Drive cycle current profile
    subplot(3, 1, 3);
    set(gca, 'Position', [0.1, 0.04, 0.85, 0.28]);  % Larger subplot
    plot(t_hours, I_measured, 'Color', '#E18727', 'LineWidth', 1);
    xlabel('Time (h)', 'FontWeight', 'bold');
    ylabel('Current (A)', 'FontWeight', 'bold');
    title('Current Profile', 'FontWeight', 'bold', 'FontSize', 10);
    grid on;
    
    % Set axis colors to black
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    
    % Overall title with capacity information
    if ~isnan(channel_capacity)
        title_str = sprintf('%s-%s-%s-Cycle%d ECM 2RC Fitting Results\nCapacity: %.4f Ah', ...
            channel_name, soc_name, dc_name, cycle_num, channel_capacity);
    else
        title_str = sprintf('%s-%s-%s-Cycle%d ECM 2RC Fitting Results', ...
            channel_name, soc_name, dc_name, cycle_num);
    end
    sgtitle(title_str, 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save visualization (channel-specific folder)
    save_dir = fullfile('Results', 'ECM_Fitting', channel_name);
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    
    fig_filename = sprintf('%s_%s_%s_Cycle%d_Fitting.fig', channel_name, soc_name, dc_name, cycle_num);
    fig_path = fullfile(save_dir, fig_filename);
    savefig(fig_path);
    
    fprintf('    Visualization saved: %s\n', fig_filename);
    
    % Close figure to save memory
    close(gcf);
end

function create_overall_error_visualization(channel_name, ecm_results)
    % Overall error visualization function
    
    result_keys = fieldnames(ecm_results);
    if isempty(result_keys)
        fprintf('No fitting results to visualize.\n');
        return;
    end
    
    % Collect RMSE, MAE, Max Error values
    rmse_values = [];
    mae_values = [];
    max_error_values = [];
    cycles = [];
    
    for i = 1:length(result_keys)
        key = result_keys{i};
        rmse_values = [rmse_values; ecm_results.(key).rmse];
        mae_values = [mae_values; ecm_results.(key).mae];
        max_error_values = [max_error_values; ecm_results.(key).max_error];
        cycles = [cycles; ecm_results.(key).Cycle];
    end
    
    % Create new figure
    figure('Position', [100, 100, 1200, 800]);
    
    % RMSE distribution
    subplot(2, 2, 1);
    histogram(rmse_values, 20, 'FaceColor', [0.2, 0.6, 0.8]);
    xlabel('RMSE (V)', 'FontWeight', 'bold');
    ylabel('Frequency', 'FontWeight', 'bold');
    title('RMSE Distribution', 'FontWeight', 'bold');
    grid on;
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    
    % MAE distribution
    subplot(2, 2, 2);
    histogram(mae_values, 20, 'FaceColor', [0.8, 0.4, 0.2]);
    xlabel('MAE (V)', 'FontWeight', 'bold');
    ylabel('Frequency', 'FontWeight', 'bold');
    title('MAE Distribution', 'FontWeight', 'bold');
    grid on;
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    
    % Max Error distribution
    subplot(2, 2, 3);
    histogram(max_error_values, 20, 'FaceColor', [0.6, 0.8, 0.2]);
    xlabel('Max Error (V)', 'FontWeight', 'bold');
    ylabel('Frequency', 'FontWeight', 'bold');
    title('Max Error Distribution', 'FontWeight', 'bold');
    grid on;
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    
    % Cycle-by-cycle RMSE boxplot
    subplot(2, 2, 4);
    unique_cycles = unique(cycles);
    cycle_rmse_data = [];
    cycle_labels = {};
    
    for i = 1:length(unique_cycles)
        cycle_mask = cycles == unique_cycles(i);
        if sum(cycle_mask) > 0
            cycle_rmse_data = [cycle_rmse_data; rmse_values(cycle_mask)'];
            cycle_labels{end+1} = sprintf('Cycle %d', unique_cycles(i));
        end
    end
    
    if ~isempty(cycle_rmse_data)
        boxplot(cycle_rmse_data, cycle_labels);
        ylabel('RMSE (V)', 'FontWeight', 'bold');
        title('RMSE by Cycle', 'FontWeight', 'bold');
        grid on;
        ax = gca;
        ax.XColor = 'k';
        ax.YColor = 'k';
    end
    
    % Overall title
    sgtitle(sprintf('%s: ECM 2RC Fitting Error Analysis', channel_name), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save visualization (channel-specific folder)
    save_dir = fullfile('Results', 'ECM_Fitting', channel_name);
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    
    fig_filename = sprintf('Overall_Error_Analysis_%s.fig', channel_name);
    fig_path = fullfile(save_dir, fig_filename);
    savefig(fig_path);
    
    fprintf('Overall error visualization saved: %s\n', fig_filename);
    
    % Close figure to save memory
    close(gcf);
end

function cost = ecm_cost_function(params, I, V_measured, t, OCV_vec)
    % ECM 2RC model cost function
    % Input: params, I (current), V_measured, t (time), OCV_vec (pre-calculated OCV vector)
    V_model = ecm_2rc_model(params, I, t, OCV_vec);
    residual = V_measured - V_model;
    cost = sqrt(mean(residual.^2));  % RMSE
end

function V_terminal = ecm_2rc_model(params, I, t, OCV_vec)
    % ECM 2RC model function
    % Input: params [R0, R1, R2, tau1, tau2], I (current), t (time), OCV_vec (pre-calculated OCV vector)
    
    % Extract parameters
    R0 = params(1);
    R1 = params(2);
    R2 = params(3);
    tau1 = params(4);
    tau2 = params(5);
    
    % Calculate time intervals
    % dt = diff(t) creates [t(2)-t(1), t(3)-t(2), ...] with length N-1
    % dt(k-1) represents the time interval between k-1 and k
    dt = diff(t);
    
    N = length(t);
    V_terminal = zeros(N, 1);
    
    % Initialize RC network state variables
    Vrc1 = 0;
    Vrc2 = 0;
    
    for k = 1:N
        % Use pre-calculated OCV value
        V_oc = OCV_vec(k);
        
        % R0 voltage drop
        IR0 = R0 * I(k);
        
        if k == 1
            % Initial state: Vrc1 and Vrc2 are zero
            Vrc1 = 0;
            Vrc2 = 0;
        else
            % Update RC1, RC2 using actual time intervals
            % dt(k-1) is the time interval between k-1 and k
            % I(k-1) is the current at time k-1
            alpha1 = exp(-dt(k-1)/tau1);
            alpha2 = exp(-dt(k-1)/tau2);
            Vrc1 = Vrc1*alpha1 + R1*(1 - alpha1)*I(k-1);
            Vrc2 = Vrc2*alpha2 + R2*(1 - alpha2)*I(k-1);
        end
        
        % Final voltage
        V_terminal(k) = V_oc + IR0 + Vrc1 + Vrc2;
    end
end

function [SOC_initial, validation_result, initial_rest_end, SOC_full, final_rest_end] = calculate_soc_docker(V_measured, I_measured, t_measured, ocv_integrated_data, rest_periods_fit)
    % Docker method SOC calculation function (two-point OCV method: SOC1, SOC2)
    % This function assumes it receives only trimmed 'rest-drive-rest' data
    % rest_periods_fit is a 2x2 array in [start, end] format (first and second rest periods only)
    % Rest periods are provided from the main script to avoid duplicate detection and ensure consistency
    
    % Create OCV inverse function (Voltage -> SOC)
    soc_grid_points = ocv_integrated_data.SOC_grid;
    ocv_values_at_soc_grid = ocv_integrated_data.V_avg_SOC;
    [ocv_values_sorted, uniqueIdx] = unique(ocv_values_at_soc_grid, 'stable');
    soc_grid_sorted = soc_grid_points(uniqueIdx);
    voltage_to_soc_interp_func = @(voltage) interp1(ocv_values_sorted, soc_grid_sorted, voltage, 'linear', 'extrap');
    
    % Validate and use provided rest period information
    if size(rest_periods_fit, 1) < 2
        error('Provided rest period information is invalid. Expected at least 2 rest periods.');
    end
    
    % Extract rest period end points
    initial_rest_end = rest_periods_fit(1, 2);
    if initial_rest_end > length(V_measured)
        error('First rest period end index exceeds data range.');
    end
    
    final_rest_end = rest_periods_fit(2, 2);
    if final_rest_end > length(V_measured)
        error('Second rest period end index exceeds data range.');
    end
    
    % Calculate SOC1 (end of first rest period)
    V_ocv_initial = V_measured(initial_rest_end);
    SOC_initial = voltage_to_soc_interp_func(V_ocv_initial);
    
    % Calculate SOC2 (end of second rest period)
    V_ocv_final = V_measured(final_rest_end);
    SOC_final = voltage_to_soc_interp_func(V_ocv_final);
    
    % Calculate full SOC profile (Docker method: two-point OCV based)
    % SOC(t) = SOC1 + (SOC2-SOC1) * (∫₀ᵗ I*dt) / (∫₀ᵉⁿᵈ I*dt)
    % where SOC1 = initial_rest_end, SOC2 = final_rest_end
    N = length(I_measured);
    
    % Calculate time intervals (using actual data time intervals)
    dt = diff(t_measured);  % Calculate actual time differences
    
    SOC_full = zeros(N, 1);
    
    % Calculate cumulative current integration (from initial_rest_end + 1 to final_rest_end)
    cumulative_current = zeros(N, 1);
    for i = (initial_rest_end + 1):final_rest_end
        % Use current at i-1 time point for integration
        cumulative_current(i) = cumulative_current(i-1) + I_measured(i-1) * dt(i-1);
    end
    
    % Total current integration (from initial_rest_end to final_rest_end)
    total_current_integration = cumulative_current(final_rest_end);
    
    % Prevent division by zero
    if total_current_integration == 0
        total_current_integration = eps;
    end
    
    % SOC calculation (Docker method)
    for i = 1:N
        if i <= initial_rest_end
            % Phase 1: Initial rest period - constant SOC1
            SOC_full(i) = SOC_initial;
        elseif i >= final_rest_end
            % Phase 3: Final rest period - constant SOC2
            SOC_full(i) = SOC_final;
        else
            % Phase 2: Active period - linear interpolation based on current integration ratio
            SOC_full(i) = SOC_initial + (SOC_final - SOC_initial) * (cumulative_current(i) / total_current_integration);
        end
    end
    
    % Validation results
    validation_result.has_validation = true;
    validation_result.initial_ocv = V_ocv_initial;
    validation_result.final_ocv = V_ocv_final;
    validation_result.initial_soc = SOC_initial;
    validation_result.final_soc = SOC_final;
    validation_result.total_current_integration = total_current_integration;
    validation_result.final_rest_end = final_rest_end;
    validation_result.soc_method = 'docker_method_two_point_ocv';
end
