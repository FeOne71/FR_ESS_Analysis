clear; clc; close all;

%% Dual Polarization ECM Model - Battery Parameter Identification
% This model distinguishes between:
% - Rohm: Ohmic resistance (immediate response)
% - Rct: Charge transfer resistance (fast electrochemical, τ1 = 1~10s)
% - Rp: Polarization resistance (slow diffusion, τ2 = 100~1000s)

%% File Directory
drive_cycle_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\parsed_data\parsedDriveCycle_0cyc_filtered.mat';
load(drive_cycle_path);
ocv_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary\OCV_interpolation_functions.mat';

load(ocv_path);

% Create output directories
FigPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\Experiment_DriveCycle_Fitting\Figures\DualPolarization';
matPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\Experiment_DriveCycle_Fitting\MatPath\DualPolarization';

if ~exist(FigPath, 'dir'), mkdir(FigPath); end
if ~exist(matPath, 'dir'), mkdir(matPath); end

% Extract channel 9 data
ch9_data = parsedDriveCycle_0cyc.ch9_Drive_0cyc;

% OCV function check
if isstruct(OCV_integrated_RPT0cyc) && isfield(OCV_integrated_RPT0cyc, 'OCV_SOC_func')
    ocv_func = OCV_integrated_RPT0cyc.OCV_SOC_func;
    fprintf('OCV_SOC function loaded successfully\n');
else
    error('OCV_SOC function not found in the structure');
end

fprintf('Dual Polarization ECM Model initialized\n');

%% Model Parameters
% Dual Polarization Model: V_terminal = V_OCV - I*Rohm - V1 - V2
% dV1/dt = -V1/(R1*C1) + I/C1  % Fast electrochemical (Rct)
% dV2/dt = -V2/(R2*C2) + I/C2  % Slow diffusion (Rp)

% Parameter definitions:
% params = [Rohm, Rct, Rp, tau1, tau2]
% where tau1 = Rct*C1 (fast), tau2 = Rp*C2 (slow)

% Initial parameter estimates
params_initial = [
    0.0001,  % Rohm (Ω) - Ohmic resistance
    0.0005,  % Rct (Ω) - Charge transfer resistance  
    0.0008,  % Rp (Ω) - Polarization resistance
    2.0,     % tau1 (s) - Fast time constant (electrochemical)
    120.0    % tau2 (s) - Slow time constant (diffusion)
];

% Multiple initial candidates for robust optimization
params_candidates = [
    params_initial;
    params_initial .* [0.8, 1.2, 0.9, 1.5, 0.8];
    params_initial .* [1.2, 0.8, 1.1, 0.7, 1.3];
    params_initial .* [0.9, 1.1, 1.2, 1.2, 0.9];
    params_initial .* [1.1, 0.9, 0.8, 0.9, 1.1];
];

% Parameter bounds
lb = [1e-5, 1e-5, 1e-5, 0.1, 10];      % Lower bounds
ub = [0.01, 0.01, 0.01, 20, 2000];     % Upper bounds

% Display parameter info
fprintf('\n=== Dual Polarization Model Parameters ===\n');
fprintf('Rohm: Ohmic resistance (immediate response)\n');
fprintf('Rct:  Charge transfer resistance (fast electrochemical, τ1 = 1~10s)\n');
fprintf('Rp:   Polarization resistance (slow diffusion, τ2 = 100~1000s)\n');
fprintf('\nParameter bounds:\n');
fprintf('  Rohm: %.1f - %.1f mΩ\n', lb(1)*1000, ub(1)*1000);
fprintf('  Rct:  %.1f - %.1f mΩ\n', lb(2)*1000, ub(2)*1000);
fprintf('  Rp:   %.1f - %.1f mΩ\n', lb(3)*1000, ub(3)*1000);
fprintf('  τ1:   %.1f - %.1f s\n', lb(4), ub(4));
fprintf('  τ2:   %.1f - %.1f s\n', lb(5), ub(5));

% Optimization options
options = optimset( ...
    'Display', 'iter', ...
    'MaxIter', 3000, ...
    'MaxFunEvals', 1e5, ...
    'TolFun', 1e-14, ...
    'TolX', 1e-15);

%% Process Drive Cycles
results = struct();

fprintf('\n=== Processing Drive Cycles with Dual Polarization Model ===\n');

for dc_idx = 1:2  % Process DC1 and DC2 first
    DC_name = sprintf('DC%d', dc_idx);
    
    if ~isfield(ch9_data.SOC90, DC_name)
        fprintf('Skipping %s - not found\n', DC_name);
        continue;
    end
    
    fprintf('\n=== Processing %s ===\n', DC_name);
    
    % Extract data
    dc_data = ch9_data.SOC90.(DC_name);
    V_measured = dc_data.V;
    I_measured = dc_data.I;
    t_measured = seconds(dc_data.t);
    
    dt = 0.1;  % Time step
    N = length(t_measured);
    
    fprintf('Data points: %d, Time step: %.3f s\n', N, dt);
    fprintf('Current range: %.3f to %.3f A\n', min(I_measured), max(I_measured));
    fprintf('Voltage range: %.3f to %.3f V\n', min(V_measured), max(V_measured));
    
    % Calculate SOC
    [SOC_initial, battery_capacity, validation_result, initial_rest_end] = estimate_and_validate_soc(V_measured, I_measured, t_measured, OCV_integrated_RPT0cyc);
    
    % Calculate SOC properly
    SOC = zeros(N, 1);
    
    % Phase 1: Rest period (SOC constant)
    SOC(1:initial_rest_end) = SOC_initial;
    
    % Phase 2: Current integration from rest period end
    for i = (initial_rest_end+1):N
        SOC_change = (I_measured(i) * dt) / (battery_capacity * 3600) * 100;
        SOC(i) = SOC(i-1) + SOC_change;
    end
    
    fprintf('SOC range: %.1f to %.1f%%\n', min(SOC), max(SOC));
    fprintf('Battery capacity: %.4f Ah\n', battery_capacity);
    
    %% Parameter Optimization
    fprintf('\n=== Parameter Optimization (%s) ===\n', DC_name);
    
    % Cost function for Dual Polarization model
    cost_function = @(params) dual_polarization_cost(params, I_measured, SOC, V_measured, t_measured, ocv_func, dt);
    
    % Test initial cost
    initial_cost = cost_function(params_initial);
    fprintf('Initial cost (RMSE): %.6f V\n', initial_cost);
    
    % Multi-start optimization
    num_runs = size(params_candidates, 1);
    all_results = [];
    
    fprintf('Starting optimization with %d initial points...\n', num_runs);
    
    for run_idx = 1:num_runs
        fprintf('\n--- Run %d/%d ---\n', run_idx, num_runs);
        p0 = params_candidates(run_idx, :);
        
        fprintf('Starting: Rohm=%.1f, Rct=%.1f, Rp=%.1f mΩ, τ1=%.1f, τ2=%.1f s\n', ...
            p0(1)*1000, p0(2)*1000, p0(3)*1000, p0(4), p0(5));
        
        [params_opt, cost_opt, exitflag, output] = fmincon( ...
            cost_function, p0, [],[],[],[], lb, ub, [], options);
        
        all_results(run_idx).params = params_opt;
        all_results(run_idx).cost = cost_opt;
        all_results(run_idx).exitflag = exitflag;
        all_results(run_idx).output = output;
        
        fprintf('Result: RMSE=%.6f V, exitflag=%d\n', cost_opt, exitflag);
    end
    
    % Select best result
    [best_cost, best_idx] = min([all_results.cost]);
    best_params = all_results(best_idx).params;
    
    fprintf('\n=== Best Result for %s ===\n', DC_name);
    fprintf('RMSE: %.6f V\n', best_cost);
    fprintf('Rohm: %.3f mΩ (Ohmic resistance)\n', best_params(1)*1000);
    fprintf('Rct:  %.3f mΩ (Charge transfer resistance)\n', best_params(2)*1000);
    fprintf('Rp:   %.3f mΩ (Polarization resistance)\n', best_params(3)*1000);
    fprintf('τ1:   %.2f s (Fast electrochemical)\n', best_params(4));
    fprintf('τ2:   %.2f s (Slow diffusion)\n', best_params(5));
    
    % Calculate capacitances
    C1 = best_params(4) / best_params(2);  % C1 = tau1 / Rct
    C2 = best_params(5) / best_params(3);  % C2 = tau2 / Rp
    
    fprintf('Calculated capacitances:\n');
    fprintf('C1: %.1f F (Electrochemical)\n', C1);
    fprintf('C2: %.1f F (Diffusion)\n', C2);
    
    %% Model Validation and Visualization
    fprintf('\n=== Model Validation ===\n');
    
    % Simulate with best parameters
    [V_sim, V1_sim, V2_sim] = dual_polarization_simulate(best_params, I_measured, SOC, ocv_func, dt);
    
    % Calculate errors
    V_error = V_measured - V_sim;
    rmse = sqrt(mean(V_error.^2));
    mae = mean(abs(V_error));
    max_error = max(abs(V_error));
    
    fprintf('Validation metrics:\n');
    fprintf('  RMSE: %.6f V\n', rmse);
    fprintf('  MAE:  %.6f V\n', mae);
    fprintf('  Max error: %.6f V\n', max_error);
    
    % Store results
    results.(DC_name) = struct();
    results.(DC_name).params = best_params;
    results.(DC_name).cost = best_cost;
    results.(DC_name).rmse = rmse;
    results.(DC_name).mae = mae;
    results.(DC_name).max_error = max_error;
    results.(DC_name).all_runs = all_results;
    results.(DC_name).capacitances = [C1, C2];
    
    %% Create Visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % Voltage comparison
    subplot(3,2,1);
    plot(t_measured/60, V_measured, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Measured');
    hold on;
    plot(t_measured/60, V_sim, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
    xlabel('Time (min)');
    ylabel('Voltage (V)');
    title(sprintf('%s - Voltage Comparison', DC_name));
    legend('Location', 'best');
    grid on;
    
    % Error plot
    subplot(3,2,2);
    plot(t_measured/60, V_error*1000, 'k-', 'LineWidth', 1);
    xlabel('Time (min)');
    ylabel('Error (mV)');
    title(sprintf('Voltage Error (RMSE: %.3f mV)', rmse*1000));
    grid on;
    
    % Current profile
    subplot(3,2,3);
    plot(t_measured/60, I_measured, 'Color', '#CD534C', 'LineWidth', 1);
    xlabel('Time (min)');
    ylabel('Current (A)');
    title('Current Profile');
    grid on;
    
    % SOC profile
    subplot(3,2,4);
    plot(t_measured/60, SOC, 'Color', '#20854E', 'LineWidth', 1.5);
    xlabel('Time (min)');
    ylabel('SOC (%)');
    title('SOC Profile');
    grid on;
    
    % Polarization voltages
    subplot(3,2,5);
    plot(t_measured/60, V1_sim*1000, 'g-', 'LineWidth', 1.5, 'DisplayName', 'V1 (Electrochemical)');
    hold on;
    plot(t_measured/60, V2_sim*1000, 'm-', 'LineWidth', 1.5, 'DisplayName', 'V2 (Diffusion)');
    xlabel('Time (min)');
    ylabel('Polarization Voltage (mV)');
    title('Polarization Components');
    legend('Location', 'best');
    grid on;
    
    % Parameter bar chart
    subplot(3,2,6);
    param_names = {'Rohm', 'Rct', 'Rp'};
    param_values = best_params(1:3) * 1000;  % Convert to mΩ
    bar(param_values);
    set(gca, 'XTickLabel', param_names);
    ylabel('Resistance (mΩ)');
    title('Identified Parameters');
    grid on;
    
    % Add text annotations for time constants
    text(0.5, max(param_values)*0.9, sprintf('τ1 = %.1f s', best_params(4)), 'FontSize', 10);
    text(0.5, max(param_values)*0.8, sprintf('τ2 = %.1f s', best_params(5)), 'FontSize', 10);
    
    sgtitle(sprintf('Dual Polarization ECM Results - %s', DC_name), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    fig_name = sprintf('%s_DualPolarization_Results.fig', DC_name);
    savefig(fullfile(FigPath, fig_name));
    fprintf('Figure saved: %s\n', fig_name);
    
    close(gcf);
end

%% Save Results
save(fullfile(matPath, 'DualPolarization_Results.mat'), 'results');
fprintf('\nResults saved to: DualPolarization_Results.mat\n');

%% Summary
fprintf('\n=== Dual Polarization Model Summary ===\n');
for dc_idx = 1:2
    DC_name = sprintf('DC%d', dc_idx);
    if isfield(results, DC_name)
        r = results.(DC_name);
        fprintf('%s: RMSE=%.3f mV, Rohm=%.2f mΩ, Rct=%.2f mΩ, Rp=%.2f mΩ\n', ...
            DC_name, r.rmse*1000, r.params(1)*1000, r.params(2)*1000, r.params(3)*1000);
    end
end

fprintf('\nDual Polarization ECM fitting completed!\n'); 