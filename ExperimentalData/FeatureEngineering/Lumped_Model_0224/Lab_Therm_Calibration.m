% Lab_Therm_Calibration.m
% Implements Phase 1-B of the lumped thermal model workflow.
% Extracts intrinsic cell thermal parameters (C_cell, R_th = R_housing_cell) 
% from laboratory continuous C-rate data (1C charge) using RPT_VQ_grid.mat.

clear; clc; close all;

fprintf('--- Starting Phase 1-B: Thermal Model Calibration (1C Charge) ---\n');

%% 1. Configuration & Data Loading
% Load RPT_VQ_grid
vq_file = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat';
if ~exist(vq_file, 'file')
    error('RPT_VQ_grid.mat not found: %s', vq_file);
end

fprintf('Loading RPT_VQ_grid.mat...\n');
load(vq_file); % Loads 'RPT_VQ_grid' structure

% Target Channel and Cycle: cyc0, Ch09, 1C charge
target_struc = RPT_VQ_grid.cyc0.Ch09.c1_charge;

% The thermal model needs time, current, voltage, and temperature
% Use the _raw arrays which correspond to actual sampling times
t_duration = target_struc.t_raw;
I_exp = target_struc.I_raw;
V_exp = target_struc.V_raw;
T_exp = target_struc.T1_raw;

% Convert duration array to seconds (double)
t_exp = seconds(t_duration);
t_exp = t_exp - t_exp(1); % Start at 0

% Environmental parameters for Lab test
% The RPT was likely conducted in a 25C chamber
T_amb = 25; 

% OCV Assumption for Entropy/Heat Generation
% For a simplified Lumped Model per Eq 1-2, Q_gen = |I * (V_oc - V)|.
% In a continuous charge, Voc is technically changing with SOC. 
% For this demonstration, we'll approximate the heat generation 
% using the average resting OCV or simply rely on Joules heating + reaction.
% To be completely strictly adhering to the paper without a full OCV lookup table:
% We will use a dynamically approximated Voc derived from the voltage curve itself 
% minus the IR drop, or simply use V_exp directly with the extracted R0.

% Load R0 from 1-A to help estimate heat generation
elec_param_file = fullfile(pwd, 'Cell_Parameters_Elec_Lab.mat');
if exist(elec_param_file, 'file')
    load(elec_param_file, 'R0_opt');
else
    warning('Electrical parameters not found, using default R0=1mOhm for Q_gen estimation.');
    R0_opt = 0.001;
end

% Approximate Voc curve = V_exp - I_exp*R0 (Since it's charge, I > 0)
% Thus Q_gen = |I * (Voc - V)| = |I * (-I*R0)| = I^2 * R0
% Wait, reaction heat also includes entropy. A better approximation for total Q_gen 
% if we lack the exact OCV curve is simply the measured overpotential * I.
Voc_approx = V_exp - (I_exp .* R0_opt);

%% 2. Thermal Model Parameter Extraction (C_cell, R_th)
fprintf('\nStarting Thermal Optimization (C_cell, R_th)...\n');
fprintf('Duration: %.1fs, Avg Current: %.2fA, T_start: %.1fC, T_end: %.1fC\n', ...
    max(t_exp), mean(I_exp), T_exp(1), T_exp(end));

% Variables: x(1) = C_cell [J/K], x(2) = R_th [K/W]
lb_th = [10,  0.01]; % Minimum constraints
ub_th = [5000, 10];  % Maximum constraints

% Objective Function: Sum of Squared Errors between T_sim and T_exp
obj_th = @(x) calc_thermal_error(x, t_exp, I_exp, V_exp, Voc_approx, T_exp, T_amb);

% Optimize using Particle Swarm
options_ps = optimoptions('particleswarm', 'SwarmSize', 30, 'Display', 'iter', 'MaxIterations', 50);
[x_th, fval_th] = particleswarm(obj_th, 2, lb_th, ub_th, options_ps);

C_cell_opt = x_th(1);
R_housing_opt = x_th(2); % R_th maps to R_housing,cell in the Field paper

fprintf('--> Thermal Calibration Results:\n');
fprintf('    C_cell = %.2f J/K\n', C_cell_opt);
fprintf('    R_th (R_housing,cell) = %.4f K/W\n', R_housing_opt);
fprintf('    Final SSE = %.4f\n', fval_th);

%% 3. Save Extracted Parameters
% Combine both electrical and thermal into one final Lab file
save_path = fullfile(pwd, 'Cell_Parameters_Lab.mat');

if exist(elec_param_file, 'file')
    load(elec_param_file); % Reload R0, R1, C1
    save(save_path, 'R0_opt', 'R1_opt', 'C1_opt', 'C_cell_opt', 'R_housing_opt');
else
    save(save_path, 'C_cell_opt', 'R_housing_opt');
end

fprintf('\nSuccessfully saved Final Lab Parameters to: %s\n', save_path);

%% 4. Visualization (Saved strictly as .fig)
fig_handle = figure('Name', 'Lab Thermal Calibration', 'Position', [150, 150, 800, 500], 'Color', 'w');

[~, T_sim] = calc_thermal_error(x_th, t_exp, I_exp, V_exp, Voc_approx, T_exp, T_amb);
plot(t_exp/60, T_exp, 'k', 'LineWidth', 2); hold on;
plot(t_exp/60, T_sim, 'r--', 'LineWidth', 2);
title('Thermal Validation (1C Continuous Charge)');
xlabel('Time [Minutes]'); ylabel('Cell Temperature [^\circC]');
legend('Experimental T1', 'Simulated T_{cell}', 'Location', 'best');
grid on;

fig_save_path = fullfile(pwd, 'Lab_Therm_Calibration_Results.fig');
savefig(fig_handle, fig_save_path);
close(fig_handle);
fprintf('Saved figure strictly as .fig to: %s\n', fig_save_path);

%% ========================================================================
% HELPER FUNCTIONS
function [sse, T_sim] = calc_thermal_error(x, t, I, V, Voc, T_exp, T_amb)
    C_cell = x(1);
    R_th = x(2); 
    
    dt = [0; diff(t)];
    N = length(t);
    T_sim = zeros(N, 1);
    
    % Initial condition matching the start of the experiment
    T_sim(1) = T_exp(1); 
    
    for k = 2:N
        % Q_gen = |I * (Voc - V)|
        % For constant current charge, heat generation is dictated by overpotential
        Q_gen = abs(I(k-1) .* (Voc(k-1) - V(k-1))); 
        
        % Q_diss = (T_amb - T_cell) / R_th
        % (Assuming chamber is vast enough to stay exactly at T_amb)
        Q_diss = (T_amb - T_sim(k-1)) / R_th;
        
        dT = ((Q_gen + Q_diss) / C_cell) * dt(k);
        T_sim(k) = T_sim(k-1) + dT;
    end
    
    % Sum of Squared Errors
    sse = sum((T_exp - T_sim).^2);
end
