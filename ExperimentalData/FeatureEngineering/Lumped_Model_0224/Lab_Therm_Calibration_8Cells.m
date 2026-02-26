% Lab_Therm_Calibration_8Cells.m
% Phase 1-B: Extracting generalized C_cell and R_th parameters
% Uses continuous 1C charge data from ALL 8 cells (Ch09~Ch16).

clear; clc; close all;

fprintf('--- Starting Phase 1-B: 8-Cell Thermal Calibration ---\n');

%% 1. Configuration & Data Loading
vq_file = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat';
if ~exist(vq_file, 'file')
    error('RPT_VQ_grid.mat not found: %s', vq_file);
end

fprintf('Loading RPT_VQ_grid.mat...\n');
load(vq_file); % Loads 'RPT_VQ_grid'

elec_param_file = fullfile(pwd, 'Cell_Parameters_Elec_LUT_8Cells.mat');
if ~exist(elec_param_file, 'file')
    error('Electrical LUT not found. Run Lab_Elec_Calibration_LUT_8Cells.m first.');
end
fprintf('Loading Electrical LUT...\n');
load(elec_param_file); % Loads soc_vector, R0_lut, R1_lut, C1_lut

channels = 9:16;
num_cells = length(channels);
T_amb = 25; % Assumed chamber temperature

% We will store thermal testing targets in a cell array
all_therm_data = cell(num_cells, 1);

fprintf('Extracting 1C Charge (c1_charge) data for Ch09 - Ch16...\n');
for c = 1:num_cells
    ch_num = channels(c);
    ch_name = sprintf('Ch%02d', ch_num);
    
    if ~isfield(RPT_VQ_grid.cyc0, ch_name)
        error('Channel %s not found in RPT_VQ_grid.cyc0', ch_name);
    end
    
    target_struc = RPT_VQ_grid.cyc0.(ch_name).c1_charge;
    
    t_sec = seconds(target_struc.t_raw);
    t_sec = t_sec - t_sec(1);
    
    I_exp = target_struc.I_raw;
    V_exp = target_struc.V_raw;
    T_exp = target_struc.T1_raw;
    
    % Approximate SOC trajectory mapped to time for this continuous charge
    % Assumes charging from ~0% to ~100% over the duration
    % (Since we just need R0 interpolation)
    soc_approx = linspace(0, 100, length(t_sec))';
    
    % Interpolate R0_lut to get R0 at every time step for Q_gen calculation
    % Note: soc_vector might be descending (eg 100:-10:0), so interp1 works either way
    % if we ensure sample points are distinct.
    % If soc_vector is descending, flip it for interp1
    if soc_vector(1) > soc_vector(end)
        soc_asc = flip(soc_vector);
        R0_asc = flip(R0_lut);
    else
        soc_asc = soc_vector;
        R0_asc = R0_lut;
    end
    
    R0_t = interp1(soc_asc, R0_asc, soc_approx, 'linear', 'extrap');
    
    % Store for objective function
    all_therm_data{c}.t = t_sec;
    all_therm_data{c}.I = I_exp;
    all_therm_data{c}.T1_raw = T_exp;
    all_therm_data{c}.R0_t = R0_t;
end

%% 2. Thermal Model Parameter Extraction Loop (8 Cells Simultaneous Fit)
% Variables: x(1) = C_cell [J/K], x(2) = R_th [K/W]
lb_th = [10,  0.01]; 
ub_th = [5000, 10];  

% Objective Function: Sum of Squared Errors between T_sim and T1_raw across 8 cells
obj_th = @(x) calc_thermal_error_multi_cell(x, all_therm_data, T_amb);

options_ps = optimoptions('particleswarm', 'SwarmSize', 30, 'Display', 'iter', 'MaxIterations', 50);

fprintf('\nStarting Thermal Optimization (C_cell, R_th) over 8 Cells...\n');
[x_th, fval_th] = particleswarm(obj_th, 2, lb_th, ub_th, options_ps);

C_cell_opt = x_th(1);
R_housing_opt = x_th(2); % R_th = R_housing,cell

fprintf('--> 8-Cell Average Thermal Calibration Results:\n');
fprintf('    C_cell = %.2f J/K\n', C_cell_opt);
fprintf('    R_th (R_housing,cell) = %.4f K/W\n', R_housing_opt);
fprintf('    Final Total SSE = %.4f\n', fval_th);

%% 3. Save Extracted Parameters
save_path = fullfile(pwd, 'Cell_Parameters_Therm_8Cells.mat');
save(save_path, 'soc_vector', 'R0_lut', 'R1_lut', 'C1_lut', 'C_cell_opt', 'R_housing_opt');
fprintf('\nSuccessfully saved Final 8-Cell Thermal Parameters to: %s\n', save_path);

%% 4. Visualization (Saved strictly as .fig)
fig_handle = figure('Name', '8-Cell Thermal Calibration', 'Position', [150, 150, 1000, 600], 'Color', 'w');

% We will plot 4 randomly selected channels to show the fit quality
plot_ch_idxs = [1, 3, 5, 8]; % Indexes in channels array (Ch09, Ch11, Ch13, Ch16)
num_subplots = length(plot_ch_idxs);

for i = 1:num_subplots
    c = plot_ch_idxs(i);
    ch_num = channels(c);
    
    t = all_therm_data{c}.t;
    I_exp = all_therm_data{c}.I;
    T_exp = all_therm_data{c}.T1_raw;
    R0_t = all_therm_data{c}.R0_t;
    
    dt = [0; diff(t)];
    N = length(t);
    T_sim = zeros(N, 1);
    T_sim(1) = T_exp(1);
    
    for k = 2:N
        % Q_gen_sim(t) = | I_raw^2 * R0(SOC(t)) |
        Q_gen = abs((I_exp(k-1)^2) * R0_t(k-1)); 
        Q_diss = (T_amb - T_sim(k-1)) / R_housing_opt;
        dT = ((Q_gen + Q_diss) / C_cell_opt) * dt(k);
        T_sim(k) = T_sim(k-1) + dT;
    end
    
    subplot(2, 2, i);
    plot(t/60, T_exp, 'k-', 'LineWidth', 2); hold on;
    plot(t/60, T_sim, 'r--', 'LineWidth', 2);
    title(sprintf('Ch%02d Thermal Fit', ch_num));
    xlabel('Time [Minutes]'); ylabel('Temperature [^\circC]');
    if i == 1, legend('Experimental T1', 'Simulated T_{cell}', 'Location', 'best'); end
    grid on;
end

fig_save_path = fullfile(pwd, 'Lab_Therm_Calibration_8Cells_Results.fig');
savefig(fig_handle, fig_save_path);
close(fig_handle);
fprintf('Saved figure strictly as .fig to: %s\n', fig_save_path);

%% ========================================================================
% HELPER FUNCTIONS
function sse_total = calc_thermal_error_multi_cell(x, all_therm_data, T_amb)
    C_cell = x(1);
    R_th = x(2); 
    
    sse_total = 0;
    num_cells = length(all_therm_data);
    
    for c = 1:num_cells
        t = all_therm_data{c}.t;
        I_exp = all_therm_data{c}.I;
        T_exp = all_therm_data{c}.T1_raw;
        R0_t = all_therm_data{c}.R0_t;
        
        dt = [0; diff(t)];
        N = length(t);
        T_sim = zeros(N, 1);
        T_sim(1) = T_exp(1); 
        
        for k = 2:N
            % Q_gen_sim(t) = | I_raw^2 * R0(SOC(t)) |
            Q_gen = abs((I_exp(k-1)^2) * R0_t(k-1)); 
            Q_diss = (T_amb - T_sim(k-1)) / R_th;
            
            dT = ((Q_gen + Q_diss) / C_cell) * dt(k);
            T_sim(k) = T_sim(k-1) + dT;
        end
        
        sse_cell = sum((T_exp - T_sim).^2);
        sse_total = sse_total + sse_cell;
    end
end
