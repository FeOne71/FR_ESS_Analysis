% Lab_Elec_Calibration_LUT_8Cells.m
% Phase 1-A: Extracting R0, R1, C1 Lookup Tables (LUT) across 0-100% SOC
% Uses multiple 1C discharge pulses from ALL 8 cells (Ch09~Ch16) from the RPT test.

clear; clc; close all;

fprintf('--- Starting Phase 1-A: 8-Cell Multi-Pulse Electrical Calibration (LUT) ---\n');

%% 1. Configuration & Data Loading
base_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';
channels = 9:16;
num_cells = length(channels);

I_1C = 64;              % 1C Current [A]
tolC_main = 0.05;       % 5% tolerance
n1C_low = -(1+tolC_main) * I_1C;
n1C_high = -(1-tolC_main) * I_1C;

% We will store all pulses per cell in a structured cell array
% all_pulses{cell_idx, pulse_idx} contains the struct for that pulse
all_pulses = cell(num_cells, 1);

fprintf('Loading parsed DCIR data for Ch09 - Ch16...\n');
for c = 1:num_cells
    ch_num = channels(c);
    parsed_file = fullfile(base_path, sprintf('RPT0_ch%02d_parsed.mat', ch_num));
    
    if ~exist(parsed_file, 'file')
        error('Parsed DCIR file not found: %s', parsed_file);
    end
    
    % Use load function and access fields explicitly to avoid variable pollution
    loaded_data = load(parsed_file, 'pdata');
    pdata_local = loaded_data.pdata;
    
    % Find valid 1C discharge pulses
    valid_pulses_for_cell = {};
    for k = 1:length(pdata_local)
        pd = pdata_local(k);
        avg_I = mean(pd.I);
        
        if length(pd.t) >= 2
            duration = max(pd.t) - min(pd.t);
        else
            duration = 0;
        end
        
        % Check if it's a 1C discharge pulse (typically ~60s duration for DCIR)
        if (avg_I >= n1C_low && avg_I <= n1C_high) && (duration > 10 && duration < 100)
            valid_pulses_for_cell{end+1} = pd; %#ok<SAGROW>
        end
    end
    
    all_pulses{c} = valid_pulses_for_cell;
    fprintf('  Ch%02d: Found %d discharge 1C pulses.\n', ch_num, length(valid_pulses_for_cell));
end

% Check if all cells have the same number of pulses
num_pulses_per_cell = cellfun(@length, all_pulses);
mode_pulses = mode(num_pulses_per_cell);

if any(num_pulses_per_cell ~= mode_pulses)
    warning('Not all cells have the same number of DCIR pulses! Using the minimum common number.');
    min_pulses = min(num_pulses_per_cell);
else
    min_pulses = mode_pulses;
end

if min_pulses == 0
    error('No valid 1C discharge pulses found in the datasets.');
end

num_pulses = min_pulses;
fprintf('\nProceeding with %d DCIR pulses (SOC steps) across all 8 cells.\n', num_pulses);

% Create SOC vector (Assuming RPT test does 100% -> 10% SOC in 10 steps, or 100%->0% in 11 steps)
if num_pulses == 10
    soc_vector = 100:-10:10;
elseif num_pulses == 11
    soc_vector = 100:-10:0;
elseif num_pulses == 9
    soc_vector = 90:-10:10; 
else
    soc_vector = linspace(100, 0, num_pulses);
end

fprintf('Assigned SOC vector: [%s] %%\n', num2str(soc_vector));

%% 2. Electrical Model Parameter Extraction Loop (8 Cells Simultaneous Fit)
% Variables: x(1) = R0, x(2) = R1, x(3) = C1
lb_elec = [1e-5, 1e-4, 10]; 
ub_elec = [1e-1, 1e-1, 500000];

R0_lut = zeros(num_pulses, 1);
R1_lut = zeros(num_pulses, 1);
C1_lut = zeros(num_pulses, 1);
sse_lut = zeros(num_pulses, 1);

options_ps = optimoptions('particleswarm', 'SwarmSize', 30, 'Display', 'off', 'MaxIterations', 50);

fprintf('\nStarting Particle Swarm Optimization for each SOC (Averaging 8 Cells)...\n');
for p = 1:num_pulses
    
    % Collect the data for this SOC from ALL 8 cells
    % We will pass this to the objective function
    soc_pulse_data = cell(num_cells, 1);
    for c = 1:num_cells
        % Extract the p-th pulse for the c-th cell
        pd = all_pulses{c}{p};
        
        % Normalize time
        t_exp = pd.t(:);
        if ~isempty(t_exp)
            t_exp = t_exp - t_exp(1);
        end
        
        soc_pulse_data{c}.t = t_exp;
        soc_pulse_data{c}.I = pd.I(:);
        soc_pulse_data{c}.V = pd.V(:);
        soc_pulse_data{c}.Voc = pd.V(1) * ones(size(t_exp)); % V(1) as Voc
    end
    
    % Objective function minimizes SSE summed across all 8 cells
    obj_elec = @(x) calc_electrical_error_multi_cell(x, soc_pulse_data);
    
    % Run PSO
    [x_elec, fval] = particleswarm(obj_elec, 3, lb_elec, ub_elec, options_ps);
    
    R0_lut(p) = x_elec(1);
    R1_lut(p) = x_elec(2);
    C1_lut(p) = x_elec(3);
    sse_lut(p) = fval;
    
    fprintf('  SOC %3.1f%%: R0 = %6.2f mOhm, R1 = %6.2f mOhm, C1 = %8.1f F (Total SSE 8 Cells: %.4f)\n', ...
        soc_vector(p), R0_lut(p)*1000, R1_lut(p)*1000, C1_lut(p), sse_lut(p));
end

%% 3. Save Extracted Parameters
save_path = fullfile(pwd, 'Cell_Parameters_Elec_LUT_8Cells.mat');
save(save_path, 'soc_vector', 'R0_lut', 'R1_lut', 'C1_lut', 'sse_lut');
fprintf('\nSuccessfully saved Electrical 8-Cell LUT to: %s\n', save_path);

%% 4. Visualization (Saved strictly as .fig)
fig_handle = figure('Name', 'Electrical Calibration 8-Cell LUT', 'Position', [100, 100, 1200, 400], 'Color', 'w');

subplot(1,3,1);
plot(soc_vector, R0_lut*1000, '-bo', 'LineWidth', 2);
title('R_0 vs SOC (8 Cells Avg)'); xlabel('SOC [%]'); ylabel('R_0 [m\Omega]'); grid on;
set(gca, 'XDir', 'reverse'); % Plot high SOC on left, low on right

subplot(1,3,2);
plot(soc_vector, R1_lut*1000, '-ro', 'LineWidth', 2);
title('R_1 vs SOC (8 Cells Avg)'); xlabel('SOC [%]'); ylabel('R_1 [m\Omega]'); grid on;
set(gca, 'XDir', 'reverse');

subplot(1,3,3);
plot(soc_vector, C1_lut, '-go', 'LineWidth', 2);
title('C_1 vs SOC (8 Cells Avg)'); xlabel('SOC [%]'); ylabel('C_1 [F]'); grid on;
set(gca, 'XDir', 'reverse');

fig_save_path = fullfile(pwd, 'Lab_Elec_Calibration_LUT_8Cells.fig');
savefig(fig_handle, fig_save_path);
close(fig_handle);
fprintf('Saved figure strictly as .fig to: %s\n', fig_save_path);

%% ========================================================================
% HELPER FUNCTIONS
function sse_total = calc_electrical_error_multi_cell(x, soc_pulse_data)
    % Calculates total SSE across all cells for a given parameter set x=[R0, R1, C1]
    R0 = x(1);
    R1 = x(2);
    C1 = x(3);
    
    sse_total = 0;
    num_cells = length(soc_pulse_data);
    
    for c = 1:num_cells
        t = soc_pulse_data{c}.t;
        I = soc_pulse_data{c}.I;
        V_exp = soc_pulse_data{c}.V;
        Voc = soc_pulse_data{c}.Voc;
        
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
        
        sse_cell = sum((V_exp - V_sim).^2);
        sse_total = sse_total + sse_cell;
    end
end
