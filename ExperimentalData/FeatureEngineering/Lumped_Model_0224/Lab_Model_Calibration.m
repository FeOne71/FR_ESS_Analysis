% Lab_Model_Calibration.m
% Implements Phase 1 of the lumped thermal model workflow.
% Extracts intrinsic cell electrical (R0, R1, C1) and thermal (C_cell, R_housing_cell)
% parameters from laboratory C-rate data.

clear; clc; close all;

%% 1. Configuration & Data Loading
fprintf('--- Starting Phase 1: Lab Data Calibration ---\n');

% The path to the processed Lab RPT C-rate data
crate_mat_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Crate_integrated\Crate_integrated.mat';
ocv_mat_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat';

if ~exist(crate_mat_path, 'file')
    error('Crate_integrated.mat not found. Please verify the path.');
end
if ~exist(ocv_mat_path, 'file')
    error('OCV_integrated.mat not found. Please verify the path.');
end

fprintf('Loading Lab C-rate and OCV data...\n');
crate_load = load(crate_mat_path);
ocv_load = load(ocv_mat_path);

% We will use Channel 09, Cycle 0, 1C Charge profile for the baseline cell calibration
target_ch = 'Ch09';
target_cyc = 'cyc0';
target_rate = 'c1';

% Extract the dynamic time series data
cell_data = crate_load.Crate_data.(target_ch).(target_cyc).(target_rate).charge;
t_exp = cell_data.t;
if isduration(t_exp)
    t_exp = seconds(t_exp);
end
I_exp = cell_data.I;
V_exp = cell_data.V;
T_exp = cell_data.T1_raw;

% Extract OCV function for this cycle
ocv_func = ocv_load.OCV_data.OCV_integrated_0.OCV_SOC_func;
capacity_Ah = ocv_load.OCV_data.OCV_integrated_0.mean_capacity;

% Start from time 0
t_exp = t_exp - t_exp(1);

% Prepare SOC sequence (Assuming starts from ~0% since it's a charge step)
% We calculate continuous SOC by coulomb counting the current
SOC_exp = cumsum(abs(I_exp)) .* [0; diff(t_exp)] / 3600 / capacity_Ah * 100;
% Make sure SOC is within bounds
SOC_exp = min(max(SOC_exp, 0), 100);

% Derive Voc at each timestep
Voc_exp = ocv_func(SOC_exp);

% Downsample for optimization speed (e.g., take 1 out of 10 points if too large)
% 1C charge takes ~1 hour = 3600s. The data might have 1Hz -> 3600 points.
downsample_rate = 5;
t_fit = t_exp(1:downsample_rate:end);
I_fit = I_exp(1:downsample_rate:end);
V_fit = V_exp(1:downsample_rate:end);
T_fit = T_exp(1:downsample_rate:end);
Voc_fit = Voc_exp(1:downsample_rate:end);
T_amb = T_fit(1); % Chamber ambient roughly equals initial soaked temperature

%% 2. Electrical Model Parameter Extraction (R0, R1, C1)
fprintf('\nStarting Electrical Optimization (R0, R1, C1)...\n');
% Variables: x(1) = R0, x(2) = R1, x(3) = C1
lb_elec = [0.00001, 0.00001, 10]; 
ub_elec = [0.1,     0.1,     50000];

% Objective Function: Sum of Squared Errors between V_sim and V_fit
obj_elec = @(x) calc_electrical_error(x, t_fit, I_fit, V_fit, Voc_fit);

% Optimize using Particle Swarm (Avoid getting stuck in generic local minima)
options_ps = optimoptions('particleswarm', 'SwarmSize', 30, 'Display', 'iter', 'MaxIterations', 50);
[x_elec, fval_elec] = particleswarm(obj_elec, 3, lb_elec, ub_elec, options_ps);

R0_opt = x_elec(1);
R1_opt = x_elec(2);
C1_opt = x_elec(3);

fprintf('--> Electrical Calibration Results:\n');
fprintf('    R0 = %.5f Ohm\n', R0_opt);
fprintf('    R1 = %.5f Ohm\n', R1_opt);
fprintf('    C1 = %.1f F\n', C1_opt);

%% 3. Thermal Model Parameter Extraction (C_cell, R_housing_cell)
fprintf('\nStarting Thermal Optimization (C_cell, R_housing_cell)...\n');
% Variables: x(1) = C_cell (J/K), x(2) = R_housing_cell (K/W)
lb_therm = [100, 0.01]; 
ub_therm = [50000, 50];

obj_therm = @(x) calc_thermal_error(x, t_fit, I_fit, V_fit, Voc_fit, T_fit, T_amb);

[x_therm, fval_therm] = particleswarm(obj_therm, 2, lb_therm, ub_therm, options_ps);

C_cell_opt = x_therm(1);
R_housing_opt = x_therm(2);

fprintf('--> Thermal Calibration Results:\n');
fprintf('    C_cell = %.1f J/K\n', C_cell_opt);
fprintf('    R_housing_cell = %.4f K/W\n', R_housing_opt);

%% 4. Save Extracted Parameters
save_path = fullfile(pwd, 'Cell_Parameters_Lab.mat');
save(save_path, 'R0_opt', 'R1_opt', 'C1_opt', 'C_cell_opt', 'R_housing_opt');
fprintf('\nSuccessfully saved Cell Parameters to: %s\n', save_path);

%% 5. Visualization (Saved strictly as .fig)
fig_handle = figure('Name', 'Lab Model Calibration', 'Position', [100, 100, 1200, 500], 'Color', 'w');

% Top left: Electrical
[~, V_sim] = calc_electrical_error(x_elec, t_fit, I_fit, V_fit, Voc_fit);
subplot(1,2,1);
plot(t_fit/60, V_fit, 'k', 'LineWidth', 2); hold on;
plot(t_fit/60, V_sim, 'r--', 'LineWidth', 2);
title('Electrical Validation (1C Charge)');
xlabel('Time [min]'); ylabel('Voltage [V]');
legend('Experimental', 'Simulated', 'Location', 'best');
grid on;

% Top right: Thermal
[~, T_sim] = calc_thermal_error(x_therm, t_fit, I_fit, V_fit, Voc_fit, T_fit, T_amb);
subplot(1,2,2);
plot(t_fit/60, T_fit, 'k', 'LineWidth', 2); hold on;
plot(t_fit/60, T_sim, 'b--', 'LineWidth', 2);
title('Thermal Validation (1C Charge)');
xlabel('Time [min]'); ylabel('Temperature [^\circC]');
legend('Experimental', 'Simulated', 'Location', 'best');
grid on;

fig_save_path = fullfile(pwd, 'Lab_Model_Calibration_Results.fig');
savefig(fig_handle, fig_save_path);
close(fig_handle);
fprintf('Saved figure strictly as .fig to: %s\n', fig_save_path);


%% ========================================================================
% HELPER FUNCTIONS
function [sse, V_sim] = calc_electrical_error(x, t, I, V_exp, Voc)
    R0 = x(1);
    R1 = x(2);
    C1 = x(3);
    
    dt = [0; diff(t)];
    N = length(t);
    V_rc = zeros(N, 1);
    V_sim = zeros(N, 1);
    
    for k = 2:N
        % dV_rc/dt = -V_rc/(R1*C1) + I/C1
        dV_rc = (-V_rc(k-1)/(R1*C1) + I(k-1)/C1) * dt(k);
        V_rc(k) = V_rc(k-1) + dV_rc;
        
        V_sim(k) = Voc(k) + I(k)*R0 + V_rc(k); % Because it's charge, current is negative usually, check signs!
        % Wait, in KEPCO data, charge current is actually often negative or positive depending on convention.
        % If I is positive for charge, then V = Voc + I*R0 + Vrc.
        % We will use V_sim = Voc + abs(I)*R0 + V_rc if I is positive for charge
    end
    
    % Adjust sign dynamically based on current polarity
    sign_I = sign(mean(I));
    if sign_I > 0
        V_sim = Voc + I.*R0 + V_rc;
    else
        V_sim = Voc + I.*R0 - V_rc; % if I is negative for charge, V = Voc - I*R0 - Vrc, wait! I*R0 is already negative. So V = Voc - I*R0
        V_sim = Voc - I.*R0 - V_rc;
    end
    
    V_sim(1) = V_exp(1);
    
    % Sum of Squared Errors
    sse = sum((V_exp - V_sim).^2);
end

function [sse, T_sim] = calc_thermal_error(x, t, I, V, Voc, T_exp, T_amb)
    C_cell = x(1);
    R_th = x(2); 
    
    dt = [0; diff(t)];
    N = length(t);
    T_sim = zeros(N, 1);
    T_sim(1) = T_exp(1); 
    
    for k = 2:N
        % Joule heating + Reaction Heat (abs(I * (Voc - V)))
        Q_gen = abs(I(k-1) .* (Voc(k-1) - V(k-1))); 
        
        Q_diss = (T_amb - T_sim(k-1)) / R_th;
        
        dT = ((Q_gen + Q_diss) / C_cell) * dt(k);
        T_sim(k) = T_sim(k-1) + dT;
    end
    
    sse = sum((T_exp - T_sim).^2);
end
