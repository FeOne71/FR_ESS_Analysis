% Lab_Elec_Calibration_LUT.m
% Phase 1-A: Extracting R0, R1, C1 Lookup Tables (LUT) across 0-100% SOC
% Uses multiple 1C discharge pulses from the RPT test.

clear; clc; close all;

fprintf('--- Starting Phase 1-A: Multi-Pulse Electrical Calibration (LUT) ---\n');

%% 1. Configuration & Data Loading
parsed_file = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed\RPT0_ch09_parsed.mat';
if ~exist(parsed_file, 'file')
    error('Parsed DCIR file not found: %s', parsed_file);
end

fprintf('Loading parsed DCIR data...\n');
load(parsed_file); % Loads 'pdata', 'data', 'intrim', etc.

data_struct = [];
for i = 1:length(pdata)
    data_struct(i).t = pdata(i).t;
    data_struct(i).V = pdata(i).V;
    data_struct(i).I = pdata(i).I;
    if isfield(pdata(i), 'type')
        data_struct(i).type = char(pdata(i).type);
    end
end

I_1C = 64;              % 1C Current [A]
tolC_main = 0.05;       % 5% tolerance

% Identify ALL discharge 1C (n1C) pulses
n1C_low = -(1+tolC_main) * I_1C;
n1C_high = -(1-tolC_main) * I_1C;

for k = 1:numel(data_struct)
    data_struct(k).avg_I = mean(data_struct(k).I);
    if numel(data_struct(k).t) >= 2
        data_struct(k).duration = max(data_struct(k).t) - min(data_struct(k).t);
    else
        data_struct(k).duration = 0;
    end
end

% A DCIR pulse is typically short (e.g., 60 seconds). 
% Discharge step for SOC change would be much longer (e.g., 6 mins+).
is_n1C_pulse = arrayfun(@(x) ...
    (x.avg_I >= n1C_low && x.avg_I <= n1C_high) && ...
    (x.duration > 10 && x.duration < 100), data_struct);

pulse_indices = find(is_n1C_pulse);
num_pulses = length(pulse_indices);

fprintf('Found %d discharge 1C DCIR pulses.\n', num_pulses);

if num_pulses == 0
    error('No discharge 1C DCIR pulse found.');
end

% Assuming standard RPT, starts at 100% SOC, goes down to 10% SOC or 0% SOC.
% typically 10 pulses for 100, 90... 10% SOC. Or 11 pulses for 100... 0%.
soc_vector = linspace(100, 100 - (100/(num_pulses-1))*(num_pulses-1), num_pulses); 
% Simplification: if 10 pulses, it will be 100, 88.8, 77.7 ... OR we can just use assumed 100:-10...
if num_pulses == 10
    soc_vector = 100:-10:10;
elseif num_pulses == 11
    soc_vector = 100:-10:0;
elseif num_pulses == 9
    soc_vector = 90:-10:10; % Example
else
    % Fallback to linear spacing
    soc_vector = linspace(100, 0, num_pulses);
end

fprintf('Assigned SOC vector: [%s] %%\n', num2str(soc_vector));

%% 2. Electrical Model Parameter Extraction Loop
% Variables: x(1) = R0, x(2) = R1, x(3) = C1
lb_elec = [1e-5, 1e-4, 10]; 
ub_elec = [1e-1, 1e-1, 500000]; % Widened bounds for C1 at low SOC

R0_lut = zeros(num_pulses, 1);
R1_lut = zeros(num_pulses, 1);
C1_lut = zeros(num_pulses, 1);
sse_lut = zeros(num_pulses, 1);

options_ps = optimoptions('particleswarm', 'SwarmSize', 30, 'Display', 'off', 'MaxIterations', 50);

fprintf('\nStarting Particle Swarm Optimization for each SOC...\n');
for i = 1:num_pulses
    idx = pulse_indices(i);
    pulse_data = data_struct(idx);
    
    t_exp = pulse_data.t(:);
    I_exp = pulse_data.I(:);
    V_exp = pulse_data.V(:);
    
    if ~isempty(t_exp)
        t_exp = t_exp - t_exp(1);
    end
    
    Voc_fit = V_exp(1) * ones(size(t_exp)); 
    
    obj_elec = @(x) calc_electrical_error(x, t_exp, I_exp, V_exp, Voc_fit);
    
    [x_elec, fval] = particleswarm(obj_elec, 3, lb_elec, ub_elec, options_ps);
    
    R0_lut(i) = x_elec(1);
    R1_lut(i) = x_elec(2);
    C1_lut(i) = x_elec(3);
    sse_lut(i) = fval;
    
    fprintf('  SOC %.1f%%: R0 = %.4f Ohm, R1 = %.4f Ohm, C1 = %.1f F (SSE: %.4f)\n', ...
        soc_vector(i), R0_lut(i), R1_lut(i), C1_lut(i), sse_lut(i));
end

%% 3. Save Extracted Parameters
save_path = fullfile(pwd, 'Cell_Parameters_Elec_LUT.mat');
save(save_path, 'soc_vector', 'R0_lut', 'R1_lut', 'C1_lut', 'sse_lut');
fprintf('\nSuccessfully saved Electrical LUT to: %s\n', save_path);

%% 4. Visualization (Saved strictly as .fig)
fig_handle = figure('Name', 'Electrical Calibration LUT', 'Position', [100, 100, 1200, 400], 'Color', 'w');

subplot(1,3,1);
plot(soc_vector, R0_lut*1000, '-o', 'LineWidth', 2);
title('R_0 vs SOC'); xlabel('SOC [%]'); ylabel('R_0 [m\Omega]'); grid on;

subplot(1,3,2);
plot(soc_vector, R1_lut*1000, '-o', 'LineWidth', 2);
title('R_1 vs SOC'); xlabel('SOC [%]'); ylabel('R_1 [m\Omega]'); grid on;

subplot(1,3,3);
plot(soc_vector, C1_lut, '-o', 'LineWidth', 2);
title('C_1 vs SOC'); xlabel('SOC [%]'); ylabel('C_1 [F]'); grid on;

fig_save_path = fullfile(pwd, 'Lab_Elec_Calibration_LUT_Results.fig');
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
        
        % V_cell(t) = V_oc + I_cell(t)*R0 + V_rc(t)
        V_sim(k) = Voc(k) + I(k)*R0 + V_rc(k);
    end
    
    V_sim(1) = V_exp(1);
    
    % Sum of Squared Errors
    sse = sum((V_exp - V_sim).^2);
end
