%% ECM 2RC Parameter Fitting
% Equivalent Circuit Model (ECM) with 2RC network parameter optimization
% Uses MultiStart optimization

clear; clc; close all;

%% 1. Load Data
fprintf('=== Loading Data ===\n');

% Load actual load data and OCV data
drive_cycle_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\parsedDriveCycle\parsedDriveCycle_0cyc.mat';
load(drive_cycle_path);

ocv_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary\OCV_interpolation_functions.mat';

load(ocv_path);
FigPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\Experiment_DriveCycle_Fitting\Figures';
matPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\Experiment_DriveCycle_Fitting\MatPath';

% Create directories if they don't exist
if ~exist(FigPath, 'dir')
    mkdir(FigPath);
    fprintf('Created FigPath directory: %s\n', FigPath);
end
if ~exist(matPath, 'dir')
    mkdir(matPath);
    fprintf('Created matPath directory: %s\n', matPath);
end

% Extract channel 9 data (actual load data)
ch9_data = parsedDriveCycle_0cyc.ch9_Drive_0cyc;

% OCV function check
if isstruct(OCV_integrated_RPT0cyc) && isfield(OCV_integrated_RPT0cyc, 'SOC_OCV_func')
    ocv_func = OCV_integrated_RPT0cyc.SOC_OCV_func;
    fprintf('SOC_OCV function loaded successfully!\n');
else
    error('SOC_OCV function not found in the structure!');
end

fprintf('Data loaded successfully!\n');

%% 2. Optimization Setup
fprintf('\n=== Optimization Setup ===\n');

% Initial parameter guess [R0, R1, R2, tau1, tau2]
params_initial = [0.001, 0.003, 0.01, 30, 600];  % [Ω, Ω, Ω, s, s]

% Parameter bounds (realistic ECM ranges for lithium-ion battery)
% R0: Ohmic resistance 0.5-10 mΩ
% R1: Fast RC resistance 1-30 mΩ  
% R2: Slow RC resistance 5-100 mΩ
% tau1: Fast time constant 1-300 s
% tau2: Slow time constant 60-1800 s (1-30 minutes)
lb = [0.00005, 0.0001, 0.0005, 0.1, 1];  % lower bounds
ub = [0.05, 0.05, 0.1, 500, 5000];       % upper bounds (tau2 max 30 min)

% Optimization options using optimoptions (modern approach)
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 3000, ...
    'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-15, ...
    'FiniteDifferenceType', 'central', ...
    'FunctionTolerance', 1e-14);

% Result storage structure
results = struct();

%% 3. DC1-DC8 Processing and Parameter Optimization
fprintf('\n=== Processing DC1-DC8 ===\n');

for dc_idx = 1:2
    DC_name = sprintf('DC%d', dc_idx);
    
    % Check if DC data exists
    if ~isfield(ch9_data.SOC90, DC_name)
        fprintf('Skipping %s - not found\n', DC_name);
        continue;
    end
    
    fprintf('\n=== Processing %s ===\n', DC_name);
    
    % Extract current DC data
    dc_data = ch9_data.SOC90.(DC_name);
    V_measured = dc_data.V;
    I_measured = dc_data.I;
    t_measured = seconds(dc_data.t);  % Convert duration to seconds
    
    % Calculate time interval and data size
    dt = 0.1;  % Fixed time step as specified
    N = length(t_measured);
    
    fprintf('Data points: %d, Time step: %.3f s\n', N, dt);
    fprintf('Time range: %.2f to %.2f seconds\n', min(t_measured), max(t_measured));
    fprintf('Current range: %.3f to %.3f A\n', min(I_measured), max(I_measured));
    fprintf('Voltage range: %.3f to %.3f V\n', min(V_measured), max(V_measured));
    
    % Check if SOC data already exists in the dataset
    fprintf('Available fields in dc_data:\n');
    disp(fieldnames(dc_data));
    
    % Use existing SOC data if available
    if isfield(dc_data, 'SOC')
        SOC = dc_data.SOC;
        fprintf('Using existing SOC data from dataset\n');
        fprintf('SOC range: %.1f to %.1f%%\n', min(SOC), max(SOC));
    elseif isfield(dc_data, 'soc')
        SOC = dc_data.soc;
        fprintf('Using existing soc data from dataset\n');
        fprintf('SOC range: %.1f to %.1f%%\n', min(SOC), max(SOC));
    else
        % Calculate initial SOC using estimate_and_validate_soc function
        [SOC_initial, battery_capacity, validation_result, initial_rest_end] = estimate_and_validate_soc(V_measured, I_measured, t_measured, OCV_integrated_RPT0cyc);
        
        % Calculate SOC properly
        SOC = zeros(N, 1);
        
        % Phase 1: Rest period (SOC constant)
        SOC(1:initial_rest_end) = SOC_initial;
        
        % Phase 2: Current integration from rest period end
        for i = (initial_rest_end+1):N
            % Current sign: charge(+) → SOC increase, discharge(-) → SOC decrease
            SOC_change = (I_measured(i) * dt) / (battery_capacity * 3600) * 100;
            SOC(i) = SOC(i-1) + SOC_change;
        end
        fprintf('SOC range: %.1f to %.1f%%\n', min(SOC), max(SOC));
    end
    
    % Parameter Optimization
    fprintf('\n=== Parameter Optimization for %s ===\n', DC_name);
    
    % Create cost function handle
    cost_handle = @(params) cost_function(params, I_measured, SOC, V_measured, t_measured, ocv_func);
    
    % Initial state output
    fprintf('Initial parameters: R0=%.6f, R1=%.6f, R2=%.6f, tau1=%.2f, tau2=%.2f\n', ...
        params_initial(1), params_initial(2), params_initial(3), params_initial(4), params_initial(5));
    
    initial_cost = cost_handle(params_initial);
    fprintf('Initial cost (RMSE): %.6f V\n', initial_cost);
    
    % MultiStart optimization (HM_2RC.m method)
    fprintf('Starting MultiStart optimization...\n');
    
    % Check Parallel Computing Toolbox availability
    try
        % Check parallel pool
        if isempty(gcp('nocreate'))
            fprintf('Starting parallel pool...\n');
        end
        use_parallel = true;
    catch
        fprintf('Parallel Computing Toolbox not available. Using serial computation.\n');
        use_parallel = false;
    end
    
    % Create MultiStart object
    if use_parallel
        ms = MultiStart( ...
            "UseParallel", true, ...  
            "Display", "iter");
    else
        ms = MultiStart( ...
            "UseParallel", false, ...  
            "Display", "iter");
    end
    
    % Generate random starting points
    nStartPts = 20;  % Number of starting points
    startPts = RandomStartPointSet('NumStartPoints', nStartPts);
    
    % Define optimization problem
    problem = createOptimProblem('fmincon', ...
        'objective', cost_handle, ...
        'x0', params_initial, ...
        'lb', lb, ...
        'ub', ub, ...
        'options', options);
    
    % Start total time measurement
    total_tic = tic;
    
    % Run MultiStart
    [params_optimal, cost_optimal, exitflag, output, solutions] = run(ms, problem, startPts);
    
    total_time = toc(total_tic);
    
    fprintf('\n=== MultiStart optimization completed for %s ===\n', DC_name);
    fprintf('Total start points: %d\n', nStartPts);
    fprintf('Solutions found: %d\n', length(solutions));
    fprintf('Best solution exit flag: %d\n', exitflag);
    
    % Result validation
    if cost_optimal < initial_cost
        fprintf('SUCCESS: Optimization improved the cost function!\n');
        fprintf('Cost improvement: %.6f V -> %.6f V (%.2f%% reduction)\n', ...
            initial_cost, cost_optimal, (1 - cost_optimal/initial_cost)*100);
    else
        fprintf('WARNING: Optimization did not improve the initial guess.\n');
        fprintf('Initial cost: %.6f V, Final cost: %.6f V\n', initial_cost, cost_optimal);
    end
    
    % Extract optimal parameters
    R0 = params_optimal(1);
    R1 = params_optimal(2);
    R2 = params_optimal(3);
    tau1 = params_optimal(4);
    tau2 = params_optimal(5);
    
    % Calculate capacitors
    C1 = tau1 / R1;
    C2 = tau2 / R2;
    
    % Result output
    fprintf('MultiStart optimization completed in %.2f seconds!\n', total_time);
    fprintf('Best RMSE achieved: %.6f V\n', cost_optimal);
    
    % Solutions analysis
    if ~isempty(solutions)
        all_costs = [solutions.Fval];
        fprintf('Solution statistics:\n');
        fprintf('  - Best RMSE: %.6f V\n', min(all_costs));
        fprintf('  - Worst RMSE: %.6f V\n', max(all_costs));
        fprintf('  - Mean RMSE: %.6f V\n', mean(all_costs));
        fprintf('  - Std RMSE: %.6f V\n', std(all_costs));
        
        % Optimal solution iteration info
        best_idx = find([solutions.Fval] == cost_optimal, 1);
        if ~isempty(best_idx)
            fprintf('  - Best solution iterations: %d\n', solutions(best_idx).Output.iterations);
        end
    end
    
    fprintf('Optimal parameters: R0=%.6f, R1=%.6f, R2=%.6f, tau1=%.2f, tau2=%.2f\n', ...
        R0, R1, R2, tau1, tau2);
    fprintf('Calculated capacitors: C1=%.2f F, C2=%.2f F\n', C1, C2);
    fprintf('Time constants: τ1=%.2f s (fast), τ2=%.2f s (slow)\n', tau1, tau2);
    
    % Model Validation
    fprintf('\n=== Model Validation for %s ===\n', DC_name);
    
    % Generate model voltage with optimal parameters
    V_model = ECM_2RC_model(params_optimal, I_measured, SOC, t_measured, ocv_func);
    
    % Calculate error metrics
    residual = V_measured - V_model;
    rmse = sqrt(mean(residual.^2));
    mae = mean(abs(residual));
    max_error = max(abs(residual));
    
    % Store results in structure
    results.(DC_name) = struct();
    results.(DC_name).params = params_optimal;
    results.(DC_name).R0 = R0;
    results.(DC_name).R1 = R1;
    results.(DC_name).R2 = R2;
    results.(DC_name).tau1 = tau1;
    results.(DC_name).tau2 = tau2;
    results.(DC_name).C1 = C1;
    results.(DC_name).C2 = C2;
    results.(DC_name).rmse = rmse;
    results.(DC_name).mae = mae;
    results.(DC_name).max_error = max_error;
    results.(DC_name).V_model = V_model;
    results.(DC_name).V_measured = V_measured;
    results.(DC_name).I_measured = I_measured;
    results.(DC_name).SOC = SOC;
    results.(DC_name).t_measured = t_measured;
    results.(DC_name).optimization_time = total_time;
    results.(DC_name).solutions = solutions;
    results.(DC_name).initial_cost = initial_cost;
    results.(DC_name).final_cost = cost_optimal;
    results.(DC_name).exitflag = exitflag;
    
    fprintf('Final RMSE: %.6f V, MAE: %.6f V\n', rmse, mae);
    fprintf('Maximum error: %.6f V\n', max_error);
    
    % Voltage comparison visualization
    fprintf('\n=== Voltage Comparison Visualization for %s ===\n', DC_name);
    
    % Create comprehensive visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % 1. Voltage comparison
    subplot(2, 2, 1);
    plot(t_measured/3600, V_measured, 'Color', '#0073C2', 'LineWidth', 1.5, 'DisplayName', 'Measured');
    hold on;
    plot(t_measured/3600, V_model, 'Color', '#CD534C','LineWidth', 1.5, 'DisplayName', 'ECM Model');
    xlabel('Time [h]');
    ylabel('Terminal Voltage [V]');
    title('Voltage Comparison');
    legend('show');
    grid on;
    
    % 2. Current profile
    subplot(2, 2, 2);
    plot(t_measured/3600, I_measured, 'Color', '#0073C2', 'LineWidth', 1);
    xlabel('Time [h]');
    ylabel('Current [A]');
    title('Current Profile');
    grid on;
    
    % 3. Error analysis
    subplot(2, 2, 3);
    plot(t_measured/3600, residual*1000, 'Color', '#CD534C', 'LineWidth', 1);
    xlabel('Time [h]');
    ylabel('Voltage Error [mV]');
    title('Model Error');
    grid on;
    
    % 4. Error histogram
    subplot(2, 2, 4);
    plot(t_measured/3600, SOC, 'Color', '#20854E', 'LineWidth', 1);
    xlabel('Time [h]');
    ylabel('SOC');
    title('SOC');
    grid on;
    
    % Add overall title
    sgtitle(sprintf('SOC90-%s ECM 2RC Model Fitting Results (RMSE: %.3f mV)', DC_name, rmse*1000));
    
    % Save voltage comparison figure
    fig_file = fullfile(FigPath, sprintf('SOC90_%s_VoltageComparison.fig', DC_name));
    savefig(fig_file);
    fprintf('Voltage comparison figure saved: %s\n', fig_file);
    
    % Cost Surface Visualization
    fprintf('\n=== Cost Surface Visualization for %s ===\n', DC_name);
    
    % Parameter names for visualization
    param_names = {'R0', 'R1', 'R2', 'tau1', 'tau2'};
    
    % Create cost surface visualization
    visualize_cost_surface(cost_handle, params_optimal, param_names);
    
    % Set figure name and save
    set(gcf, 'Name', sprintf('SOC90-%s Cost Surface', DC_name));
    fig_file = fullfile(FigPath, sprintf('SOC90_%s_CostSurface.fig', DC_name));
    savefig(fig_file);
    fprintf('Cost surface figure saved: %s\n', fig_file);
    
    fprintf('Completed %s processing\n', DC_name);
end

%% 4. Overall Results Summary
fprintf('\n=== Summary of All DC Results ===\n');
fprintf('%-4s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s\n', ...
    'DC', 'R0(Ω)', 'R1(Ω)', 'R2(Ω)', 'τ1(s)', 'τ2(s)', 'C1(F)', 'C2(F)', 'RMSE(V)', 'Time(s)');
fprintf('%s\n', repmat('-', 1, 100));

for dc_idx = 1:2
    DC_name = sprintf('DC%d', dc_idx);
    if isfield(results, DC_name)
        r = results.(DC_name);
        fprintf('%-4s | %.6f | %.6f | %.6f | %8.2f | %8.2f | %8.2f | %8.2f | %.6f | %8.1f\n', ...
            DC_name, r.R0, r.R1, r.R2, r.tau1, r.tau2, r.C1, r.C2, r.rmse, r.optimization_time);
    end
end

%% 5. Save Results
% Save all results
mat_file = fullfile(matPath, 'ECM_2RC_fitting_results_SOC90_DC1to2.mat');
save(mat_file, 'results');
fprintf('\nAll results saved to: %s\n', mat_file);

fprintf('\n=== Analysis Complete ===\n');
fprintf('Total processing completed for SOC90 DC1-DC8\n');
fprintf('Results structure contains:\n');
fprintf('  - Individual DC results: results.DC1, results.DC2, ..., results.DC8\n');
fprintf('  - Voltage comparison figures: SOC90_DCX_VoltageComparison.fig\n');
fprintf('  - Cost surface figures: SOC90_DCX_CostSurface.fig\n'); 