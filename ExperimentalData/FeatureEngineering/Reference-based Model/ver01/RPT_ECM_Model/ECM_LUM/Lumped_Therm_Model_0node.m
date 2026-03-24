%% Lumped Thermal Model (0-Node) with DCIR 60s Estimation
% Based on Reference SOC Calculation (Global Integration)
% Data Source: Parsed RPT MAT file (e.g., RPT0_ch09_parsed.mat)

clear; clc; close all;

%% 1. Initialization and Parameter Configuration
clc; clear; close all;
warning('off', 'all');

C_nom = 64.0;
I_1C = 64.0;

% Load OCV Grid Global Data
ocvFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat';
if exist(ocvFile, 'file')
    ocv_data = load(ocvFile);
    fprintf('=> Loaded Master OCV Grid Data.\n');
else
    error('OCV File not found: %s', ocvFile);
end

channels = 9:16;
num_channels = length(channels);

% Preallocate arrays for parfor compatibility
res_C_cell = zeros(num_channels, 1);
res_R_th = zeros(num_channels, 1);
res_RMSE = zeros(num_channels, 1);

fprintf('\n=== Starting Multi-Channel Parallel Processing (Ch09~Ch16) ===\n');

parfor iter = 1:num_channels
    chNum = channels(iter);
    ch_str = sprintf('ch%02d', chNum);
    Ch_str_cap = sprintf('Ch%02d', chNum);

    fprintf('>> Processing %s...\n', ch_str);

    % File paths
    dataFile = sprintf('D:\\JCW\\Projects\\KEPCO_ESS_Local\\ExperimentalData\\RPT\\Postprocessing\\Parsed\\RPT0_%s_parsed.mat', ch_str);

    % if ~exist(dataFile, 'file')
    %     error('File not found: %s\nPlease check the path.', dataFile);
    % end

    % Tuning needed according to actual values
    % C_nom = 64;          % Measured reference capacity or 1C current value (Ah base) - Moved to global
    % T_amb = 25;          % Ambient temperature (degC)

    % RPT Analysis Conditions
    % I_1C = 64;           % 1C Current [A] (Assume -64 for discharge, +64 for charge) - Moved to global
    tolC = 0.05;         % 1C identification tolerance (5%)
    V_max = 4.2;         % Upper limit voltage for charge pulse power calculation [V]
    V_min = 3.0;         % Lower limit voltage for discharge pulse power calculation [V]

    %% 2. Data Preparation (Identify and extract only DCIR evaluation steps from parsed pdata array)
    fprintf('== Data Loading: %s ==\n', dataFile);
    loaded_data = load(dataFile); % Explicit assignment required for parfor transparency
    pdata = loaded_data.pdata;
    % 
    % if ~isfield(loaded_data, 'pdata')
    %     error('pdata variable is not found in the loaded file. Check the parsed data format.');
    % end

    %% ---------------------------------------------------------
    % 2.1. Identify Discharge Window
    % ---------------------------------------------------------
    disch_n1C_idx = [];
    disch_c02_idx = [];

    for i = 1:length(pdata)
        if isempty(pdata(i).I), continue; end
        avg_I = mean(pdata(i).I);
        dur = pdata(i).t(end) - pdata(i).t(1);

        if avg_I <= -I_1C*(1-tolC) && avg_I >= -I_1C*(1+tolC) && dur >= 50 && dur <= 70
            disch_n1C_idx(end+1) = i;
        end
        if avg_I <= -12.8*(1-tolC) && avg_I >= -12.8*(1+tolC)
            disch_c02_idx(end+1) = i;
        end
    end

    if isempty(disch_n1C_idx)
        warning('Discharge DCIR 1C pulse steps were not found.');
        dchg_1C_idx = [];
    else
        % Extract starting from the Rest step (-1) immediately before the first 1C discharge pulse to secure initial OCV
        % disch_start = max(1, disch_n1C_idx(1) - 1); % This is now handled by full_valid_idx construction
        dchg_1C_idx = disch_n1C_idx;
        % Find the end of 0.2C discharge, if none exists, use the last 1C pulse as the end
        if ~isempty(disch_c02_idx)
            % disch_end = max(disch_c02_idx); % This is now handled by full_valid_idx construction
            dchg_02C_idx = disch_c02_idx;
        else
            dchg_02C_idx = disch_n1C_idx(end); % If no 0.2C, use last 1C pulse
        end
        % disch_valid_idx = disch_start:disch_end; % This is now handled by full_valid_idx construction
        fprintf('  [Discharge] 1C Pulse idx: '); fprintf('%d ', disch_n1C_idx); fprintf('\n');
        fprintf('  [Discharge] 0.2C idx: '); fprintf('%d ', disch_c02_idx); fprintf('\n');
        % fprintf('  [Discharge Window] Extracting steps %d ~ %d\n', disch_start, disch_end);
    end

    % ---------------------------------------------------------
    % 2.2. Identify Charge Window
    % ---------------------------------------------------------
    charge_n1C_idx = [];
    charge_c02_idx = [];

    for i = 1:length(pdata)
        if isempty(pdata(i).I), continue; end
        avg_I = mean(pdata(i).I);
        dur = pdata(i).t(end) - pdata(i).t(1);

        if avg_I >= I_1C*(1-tolC) && avg_I <= I_1C*(1+tolC) && dur >= 50 && dur <= 70
            charge_n1C_idx(end+1) = i;
        end
        if avg_I >= 12.8*(1-tolC) && avg_I <= 12.8*(1+tolC)
            charge_c02_idx(end+1) = i;
        end
    end

    if isempty(charge_n1C_idx)
        warning('Charge DCIR 1C pulse steps were not found.');
        chg_1C_idx = [];
    else
        % Extract starting from the Rest step (-1) immediately before the first 1C charge pulse to secure initial OCV
        % charge_start = max(1, charge_n1C_idx(1) - 1); % This is now handled by full_valid_idx construction
        chg_1C_idx = charge_n1C_idx;
        if ~isempty(charge_c02_idx)
            % charge_end = max(charge_c02_idx); % This is now handled by full_valid_idx construction
            chg_02C_idx = charge_c02_idx;
        else
            chg_02C_idx = charge_n1C_idx(end); % If no 0.2C, use last 1C pulse
        end
        % charge_valid_idx = charge_start:charge_end; % This is now handled by full_valid_idx construction
        fprintf('  [Charge] 1C Pulse idx: '); fprintf('%d ', charge_n1C_idx); fprintf('\n');
        fprintf('  [Charge] 0.2C idx: '); fprintf('%d ', charge_c02_idx); fprintf('\n');
        % fprintf('  [Charge Window] Extracting steps %d ~ %d\n', charge_start, charge_end);
    end


    %% 3. Integrated Data Extraction (Unify Charge/Discharge windows into a single time series)
    % Assuming charge is performed first and discharge later (Single continuous RPT)
    % start_idx = min(charge_start, disch_start); % Replaced by below
    % end_idx = max(charge_end, disch_end); % Replaced by below
    % fprintf('  [Unified Window] Extracting steps %d ~ %d\n', start_idx, end_idx);

    % 2.3 Extract Windows
    % Ensure indices are valid and sorted
    all_indices = unique([max(1, chg_1C_idx(1)-1):chg_02C_idx(end), max(1, dchg_1C_idx(1)-1):dchg_02C_idx(end)]);
    full_valid_idx = all_indices(all_indices <= length(pdata)); % Ensure indices don't exceed pdata length

    [t_all, I_all, V_all, T_meas_all, SOC_all, dt_all] = extract_window_data(pdata, full_valid_idx, 0.0, C_nom);

    % 2.5 Map Dynamic V_ocv (Overpotential Model)
    ocv_struct = ocv_data.RPT_VQ_grid.cyc0.(Ch_str_cap);

    % Map Dynamic V_ocv based on accumulated Capacity (SOC * C_nom)
    V_ocv_all = zeros(size(t_all));
    for i=1:length(t_all)
        if I_all(i) > 0.1 % 충전(Charge) 구간 (+ 양수 전류)
            q_grid = ocv_struct.OCV_charge.Q;
            v_grid = ocv_struct.OCV_charge.V_grid;
            current_Q = SOC_all(i) * max(q_grid); % 충전 누적량: 채널별 고유 최대 용량으로 SOC 매핑
        else % 방전(Discharge) 및 휴지기 구간 (0A 포함, - 음수 전류)
            q_grid = ocv_struct.OCV_discharge.Q;
            v_grid = ocv_struct.OCV_discharge.V_grid;
            current_Q = (1.0 - SOC_all(i)) * max(q_grid); % 방전 누적량 역전: 채널별 고유 최대 용량으로 SOC 매핑
        end

        % 중복 Q 값 제거 (interp1 조건 충족)
        [q_uniq, idx_uniq] = unique(q_grid);
        v_uniq = v_grid(idx_uniq);

        % Interpolate OCV
        if length(q_uniq) > 1
            V_ocv_all(i) = interp1(q_uniq, v_uniq, current_Q, 'linear', 'extrap');
        else
            V_ocv_all(i) = v_uniq(1); % Handle case with single unique point
        end
    end

    % 3.1 주변 온도(T_amb) 
    T_amb = 25.0;

    % 4. Multi-duration (2s, 10s, 30s, 60s) Global DCIR and Power Calculation
    % (생략: 터미널 출력 최소화, 계산 자체는 유지 필요 시 코드를 유지)
    [R_ch_res, P_ch_res, R_disch_res, P_disch_res, R_series_all] = calc_all_dcir_power(t_all, I_all, V_all, dt_all, SOC_all, I_1C, tolC, V_max, V_min);

    fprintf('  Detected Discharge DCIR points: %d\n', length(R_disch_res.soc));
    fprintf('--------------------------------------------------------------------------------------------------------\n');
    fprintf('  SOC(%%) |  R_2s(mOhm) | R_10s(mOhm) | R_30s(mOhm) | R_60s(mOhm) |  P_2s(W)   |  P_10s(W)  |  P_60s(W)  \n');
    fprintf('--------------------------------------------------------------------------------------------------------\n');
    for i = 1:length(R_disch_res.soc)
        fprintf(' %6.1f | %11.2f | %11.2f | %11.2f | %11.2f | %10.1f | %10.1f | %10.1f\n', ...
            R_disch_res.soc(i)*100, R_disch_res.R2(i)*1000, R_disch_res.R10(i)*1000, R_disch_res.R30(i)*1000, R_disch_res.R60(i)*1000, ...
            P_disch_res.P2(i), P_disch_res.P10(i), P_disch_res.P60(i));
    end

    fprintf('\n  Detected Charge DCIR points: %d\n', length(R_ch_res.soc));
    fprintf('--------------------------------------------------------------------------------------------------------\n');
    fprintf('  SOC(%%) |  R_2s(mOhm) | R_10s(mOhm) | R_30s(mOhm) | R_60s(mOhm) |  P_2s(W)   |  P_10s(W)  |  P_60s(W)  \n');
    fprintf('--------------------------------------------------------------------------------------------------------\n');
    for i = 1:length(R_ch_res.soc)
        fprintf(' %6.1f | %11.2f | %11.2f | %11.2f | %11.2f | %10.1f | %10.1f | %10.1f\n', ...
            R_ch_res.soc(i)*100, R_ch_res.R2(i)*1000, R_ch_res.R10(i)*1000, R_ch_res.R30(i)*1000, R_ch_res.R60(i)*1000, ...
            P_ch_res.P2(i), P_ch_res.P10(i), P_ch_res.P60(i));
    end


    %% 5. Global Optimization (PSO) - Unified Data Based
    fprintf('\n== Lumped Thermal Model (0-Node) Parameter Fitting (Global PSO) ==\n');

    % Output detected SOC
    fprintf('  [Discharge] Detected DCIR SOC (%%): ');
    fprintf('%.1f ', R_disch_res.soc*100);
    fprintf('\n');

    fprintf('  [Charge] Detected DCIR SOC (%%): ');
    fprintf('%.1f ', R_ch_res.soc*100);
    fprintf('\n');

    % Physics-based initial values and bound settings (Based on 64Ah Battery)
    initial_guess = [1000.0, 1.0]; % [C_cell, R_th]

    % C_cell: 100~20000 [J/K], R_th: 0.1~50.0 [K/W]
    lb = [100.0, 0.1];
    ub = [5000.0, 5.0];

    % PSO Options: Outer loop is already parallel, inner must be serial
    % Adding fmincon as a HybridFcn to polish the optimum valley.
    options = optimoptions('particleswarm', ...
        'Display', 'none', ...
        'SwarmSize', 50, ...
        'MaxIterations', 1000, ...
        'InitialSwarmMatrix', initial_guess, ...
        'UseParallel', false, ...
        'HybridFcn', @fmincon);

    fprintf('  [Fitting Unified Time-Series Model (Parallel Computing 2 Parameters)...]\n');
    [opt_all, fval_all] = particleswarm(@(x) obj_therm_model(x, t_all, I_all, V_all, V_ocv_all, T_meas_all, T_amb), 2, lb, ub, options);
    C_cell_opt = opt_all(1); R_th_opt = opt_all(2);
    RMSE_all = fval_all;

    res_C_cell(iter) = C_cell_opt;
    res_R_th(iter) = R_th_opt;
    res_RMSE(iter) = RMSE_all;

    fprintf('  => Optimal Parameters [Unified]: C_cell = %.1f [J/K], R_th = %.3f [K/W]\n  => Unified RMSE = %.3f °C\n', C_cell_opt, R_th_opt, RMSE_all);

    T_pred_all = simulate_thermal_model(t_all, I_all, V_all, V_ocv_all, T_meas_all(1), T_amb, C_cell_opt, R_th_opt);
    fprintf('     (Analysis) Max Measured T (T_meas): %.2f °C, Max Predicted T (T_pred): %.2f °C\n', max(T_meas_all), max(T_pred_all));


    %% 6. Result Visualization
    fprintf('\n== Generating Visualizations ==\n');

    % 6.1 Multi-duration DCIR Resistance Visualization (1 Figure, 2 Subplots)
    fig_R = figure('Name', sprintf('Multi-duration DCIR Profile %s', ch_str), 'Color', [1 1 1], 'Position', [100 100 1200 450], 'Visible', 'off');

    subplot(1,2,1); hold on; grid on;
    if ~isempty(R_ch_res.soc)
        plot(R_ch_res.soc*100, R_ch_res.R2*1000, 'o-', 'LineWidth', 1.5, 'DisplayName', '2s');
        plot(R_ch_res.soc*100, R_ch_res.R10*1000, 's-', 'LineWidth', 1.5, 'DisplayName', '10s');
        plot(R_ch_res.soc*100, R_ch_res.R30*1000, 'd-', 'LineWidth', 1.5, 'DisplayName', '30s');
        plot(R_ch_res.soc*100, R_ch_res.R60*1000, '^-', 'LineWidth', 2.0, 'MarkerSize', 7, 'DisplayName', '60s');
    end
    xlabel('SOC (%)'); ylabel('Resistance (m\Omega)');
    title('Charge DCIR'); legend('Location', 'best');

    subplot(1,2,2); hold on; grid on;
    if ~isempty(R_disch_res.soc)
        plot(R_disch_res.soc*100, R_disch_res.R2*1000, 'o-', 'LineWidth', 1.5, 'DisplayName', '2s');
        plot(R_disch_res.soc*100, R_disch_res.R10*1000, 's-', 'LineWidth', 1.5, 'DisplayName', '10s');
        plot(R_disch_res.soc*100, R_disch_res.R30*1000, 'd-', 'LineWidth', 1.5, 'DisplayName', '30s');
        plot(R_disch_res.soc*100, R_disch_res.R60*1000, 'v-', 'LineWidth', 2.0, 'MarkerSize', 7, 'DisplayName', '60s');
    end
    xlabel('SOC (%)'); ylabel('Resistance (m\Omega)');
    title('Discharge DCIR'); legend('Location', 'best');
    saveas(fig_R, sprintf('Figure_DCIR_%s.fig', ch_str));
    close(fig_R);

    % 6.1.5 Multi-duration Power Capability Visualization (1 Figure, 2 Subplots)
    fig_P = figure('Name', sprintf('Pulse Power Capability Profile %s', ch_str), 'Color', [1 1 1], 'Position', [150 150 1200 450], 'Visible', 'off');

    subplot(1,2,1); hold on; grid on;
    if ~isempty(P_ch_res.soc)
        plot(P_ch_res.soc*100, P_ch_res.P2, 'o-', 'LineWidth', 1.5, 'DisplayName', '2s');
        plot(P_ch_res.soc*100, P_ch_res.P10, 's-', 'LineWidth', 1.5, 'DisplayName', '10s');
        plot(P_ch_res.soc*100, P_ch_res.P30, 'd-', 'LineWidth', 1.5, 'DisplayName', '30s');
        plot(P_ch_res.soc*100, P_ch_res.P60, '^-', 'LineWidth', 2.0, 'MarkerSize', 7, 'DisplayName', '60s');
    end
    xlabel('SOC (%)'); ylabel('Pulse Power Capability (W)');
    title('Charge Power Capability'); legend('Location', 'best');

    subplot(1,2,2); hold on; grid on;
    if ~isempty(P_disch_res.soc)
        plot(P_disch_res.soc*100, P_disch_res.P2, 'o-', 'LineWidth', 1.5, 'DisplayName', '2s');
        plot(P_disch_res.soc*100, P_disch_res.P10, 's-', 'LineWidth', 1.5, 'DisplayName', '10s');
        plot(P_disch_res.soc*100, P_disch_res.P30, 'd-', 'LineWidth', 1.5, 'DisplayName', '30s');
        plot(P_disch_res.soc*100, P_disch_res.P60, 'v-', 'LineWidth', 2.0, 'MarkerSize', 7, 'DisplayName', '60s');
    end
    xlabel('SOC (%)'); ylabel('Pulse Power Capability (W)');
    title('Discharge Power Capability'); legend('Location', 'best');
    saveas(fig_P, sprintf('Figure_Power_%s.fig', ch_str));
    close(fig_P);

    % 6.2 Voltage/Current and SOC Visualization (1 Figure, Unified)
    fig_elec = figure('Name', sprintf('Global Electrical & SOC %s', ch_str), 'Color', [1 1 1], 'Position', [100 200 1000 600], 'Visible', 'off');

    subplot(3,1,1); plot(t_all, I_all, 'b'); grid on; ylabel('Current [A]'); title('Global Electrical States');
    subplot(3,1,2); plot(t_all, V_all, 'k'); grid on; ylabel('Voltage [V]');
    subplot(3,1,3); plot(t_all, SOC_all*100, 'g', 'LineWidth', 1.5); grid on; ylabel('SOC [%]'); xlabel('Time [s]');
    saveas(fig_elec, sprintf('Figure_Electrical_%s.fig', ch_str));
    close(fig_elec);

    % 6.3 Temperature Prediction and Measured Validation Visualization (1 Figure, Unified)
    fig_therm = figure('Name', sprintf('Global Thermal Validation %s', ch_str), 'Color', [1 1 1], 'Position', [150 250 1000 400], 'Visible', 'on');
    plot(t_all/3600, T_meas_all, 'ko-', 'MarkerSize', 2, 'DisplayName', 'Exp.data');
    ylim([20 30]);
    hold on;
    plot(t_all/3600, T_pred_all, 'ro-', 'MarkerSize', 2, 'DisplayName', 'T_{model}');
    ylabel('Temperature [°C]'); grid on;
    xlabel('Time [h]');
    ylim([20 30]);
    % yyaxis right;
    % plot(t_all, R_series_all*1000, 'b-', 'LineWidth', 1, 'DisplayName', 'R_{60}');
    % ylabel('Resistance [m\Omega]');
    title(sprintf('Global Fit %s | C_{heat, cell}: %.1f, R_{thermal}: %.3f | RMSE: %.3f°C', ch_str, C_cell_opt, R_th_opt, RMSE_all));
    legend('Location', 'northwest');

    saveas(fig_therm, sprintf('fig_therm_unified_0node_%s.fig', ch_str));
    close(fig_therm); % 메모리 누수 방지

    % 6.4 3D Objective Function Landscape Visualization 
    % (User explicitly requested to show without saving)
    fprintf('  [3D Landscape] Calculating grid for %s...\n', ch_str);
    c_range = linspace(100, 5000, 40);  
    r_range = linspace(0.1, 100, 40);     
    [C_grid, R_grid] = meshgrid(c_range, r_range);
    RMSE_grid = zeros(size(C_grid));
    
    for rc_i = 1:size(C_grid, 1)
        for rc_j = 1:size(C_grid, 2)
            RMSE_grid(rc_i, rc_j) = obj_therm_model([C_grid(rc_i, rc_j), R_grid(rc_i, rc_j)], t_all, I_all, V_all, V_ocv_all, T_meas_all, T_amb);
        end
    end
    
    % Z축 스케일 폭발 방지 (ODE 수치해석 발산 영역 C_cell=100, dt=Large 방어)
    RMSE_grid_plot = RMSE_grid;
    RMSE_grid_plot(RMSE_grid_plot > 5.0) = 5.0; 
    RMSE_grid_plot(isnan(RMSE_grid_plot)) = 5.0;
    
    fig_3d = figure('Name', sprintf('Thermal Landscape 3D %s', ch_str), 'Color', [1 1 1], 'Position', [200 200 800 600], 'Visible', 'on');
    surf(C_grid, R_grid, RMSE_grid_plot, 'EdgeColor', 'none'); 
    colormap parula; colorbar;
    xlabel('C_{cell} (J/K)', 'FontWeight', 'bold');
    ylabel('R_{th} (K/W)', 'FontWeight', 'bold');
    zlabel('RMSE (°C)', 'FontWeight', 'bold');
    title(sprintf('Obj Func Landscape %s | Opt: C=%.1f, R=%.3f', ch_str, C_cell_opt, R_th_opt));
    view(-45, 30);
    hold on;
    plot3(C_cell_opt, R_th_opt, RMSE_all, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    hold off;

    saveas(fig_3d, sprintf('Figure_grid_%s.fig', ch_str));

    fprintf('>> [%s] Done. C_cell: %.1f, R_th: %.3f\n', ch_str, C_cell_opt, R_th_opt);
end % End of parfor

%% Compile Results into Struct
Thermal_Fit_Results = struct();
for iter = 1:num_channels
    ch_str = sprintf('ch%02d', channels(iter));
    Thermal_Fit_Results.(ch_str).C_cell = res_C_cell(iter);
    Thermal_Fit_Results.(ch_str).R_th = res_R_th(iter);
    Thermal_Fit_Results.(ch_str).RMSE = res_RMSE(iter);
end

fprintf('\n=== Script Execution Completed ===\n');
disp('Thermal Parameter Results:');
disp(Thermal_Fit_Results);
save('Thermal_Fit_Results_Ch09_16.mat', 'Thermal_Fit_Results');


%% [Local Helper Functions]
function [t, I, V, T_meas, SOC, dt] = extract_window_data(pdata, valid_idx, initial_soc, C_nom)
if isempty(valid_idx)
    t=[]; I=[]; V=[]; T_meas=[]; SOC=[]; dt=0.1; return;
end

pdata_filtered = pdata(valid_idx);

I = vertcat(pdata_filtered.I);
V = vertcat(pdata_filtered.V);

if isfield(pdata, 'Batt_Temp') && ~isempty(pdata(1).Batt_Temp)
    T_meas = vertcat(pdata_filtered.Batt_Temp);
end

% Build accumulated time (based on steptime_double)
if isfield(pdata, 'steptime_double') && ~isempty(pdata(1).steptime_double)
    steptime_concat = vertcat(pdata_filtered.steptime_double);
    steptime_concat(1) = 0;

    t = zeros(size(steptime_concat));
    for i=2:length(t)
        dt_val = steptime_concat(i) - steptime_concat(i-1);
        if dt_val < 0
            t(i) = t(i-1) + 0.1;
        else
            t(i) = t(i-1) + dt_val;
        end
    end
else
    t = vertcat(pdata_filtered.t);
    t = t - t(1);
end

dt = median(diff(t));
if isnan(dt) || dt <= 0, dt = 0.1; end

% Eliminate temperature noise (Savitzky-Golay filter, 60-second window, degree 2)
% window_size = round(60 / dt);
% if mod(window_size, 2) == 0, window_size = window_size + 1; end % sgolay 윈도우는 홀수여야 함
% if window_size >= 3
%     T_meas = smoothdata(T_meas, 'sgolay', window_size);
% end

SOC = initial_soc + cumtrapz(t, I) / (3600 * C_nom);
end

function [R_ch_res, P_ch_res, R_disch_res, P_disch_res, R_series] = calc_all_dcir_power(t, I, V, dt, SOC, I_1C, tolC, V_max, V_min)
R60_interp = ones(size(t)) * NaN;

ch_soc=[]; ch_R2=[]; ch_R10=[]; ch_R30=[]; ch_R60=[];
ch_P2=[]; ch_P10=[]; ch_P30=[]; ch_P60=[];

disch_soc=[]; disch_R2=[]; disch_R10=[]; disch_R30=[]; disch_R60=[];
disch_P2=[]; disch_P10=[]; disch_P30=[]; disch_P60=[];

threshold = abs(I_1C) * 0.5;
is_active = abs(I) >= threshold;
edges = find(diff([0; is_active]) == 1);

for k = 1:length(edges)
    idx_active = edges(k);

    % Trace back to find start
    idx_start = idx_active;
    while idx_start > 1 && abs(I(idx_start)) > 1.0
        idx_start = idx_start - 1;
    end
    % Ensure we don't pick a transition point reading by going 1 extra step back if the loop stopped exactly at the boundary
    % Using the reliable rest baseline:
    idx_rest = max(1, idx_start);

    idx_2s = find(t >= t(idx_active) + 2, 1, 'first');
    idx_10s = find(t >= t(idx_active) + 10, 1, 'first');
    idx_30s = find(t >= t(idx_active) + 30, 1, 'first');
    idx_60s = find(t >= t(idx_active) + 59, 1, 'first'); % 59s safe

    if ~isempty(idx_60s) && ~isempty(idx_30s) && ~isempty(idx_10s) && ~isempty(idx_2s)
        % Determine if it's Charge or Discharge based on Current sign
        is_charge = I(idx_60s) > 0;

        if (is_charge && I(idx_60s) >= I_1C*(1-tolC)) || (~is_charge && I(idx_60s) <= -I_1C*(1-tolC))
            V_start = V(idx_rest);
            I_start = I(idx_rest);

            % Reference: dV / dI
            R2 = abs((V(idx_2s) - V_start) / (I(idx_2s) - I_start));
            R10 = abs((V(idx_10s) - V_start) / (I(idx_10s) - I_start));
            R30 = abs((V(idx_30s) - V_start) / (I(idx_30s) - I_start));
            R60 = abs((V(idx_60s) - V_start) / (I(idx_60s) - I_start));

            R60_interp(idx_rest:idx_60s) = R60;

            if ~is_charge
                V_limit = V_min;
                P2 = V_limit * (V_start - V_limit) / R2;
                P10 = V_limit * (V_start - V_limit) / R10;
                P30 = V_limit * (V_start - V_limit) / R30;
                P60 = V_limit * (V_start - V_limit) / R60;

                disch_soc(end+1) = SOC(idx_rest);
                disch_R2(end+1)=R2; disch_R10(end+1)=R10; disch_R30(end+1)=R30; disch_R60(end+1)=R60;
                disch_P2(end+1)=P2; disch_P10(end+1)=P10; disch_P30(end+1)=P30; disch_P60(end+1)=P60;
            else
                V_limit = V_max;
                P2 = V_limit * (V_limit - V_start) / R2;
                P10 = V_limit * (V_limit - V_start) / R10;
                P30 = V_limit * (V_limit - V_start) / R30;
                P60 = V_limit * (V_limit - V_start) / R60;

                ch_soc(end+1) = SOC(idx_rest);
                ch_R2(end+1)=R2; ch_R10(end+1)=R10; ch_R30(end+1)=R30; ch_R60(end+1)=R60;
                ch_P2(end+1)=P2; ch_P10(end+1)=P10; ch_P30(end+1)=P30; ch_P60(end+1)=P60;
            end
        end
    end
end

R_ch_res.soc=ch_soc; R_ch_res.R2=ch_R2; R_ch_res.R10=ch_R10; R_ch_res.R30=ch_R30; R_ch_res.R60=ch_R60;
P_ch_res.soc=ch_soc; P_ch_res.P2=ch_P2; P_ch_res.P10=ch_P10; P_ch_res.P30=ch_P30; P_ch_res.P60=ch_P60;

R_disch_res.soc=disch_soc; R_disch_res.R2=disch_R2; R_disch_res.R10=disch_R10; R_disch_res.R30=disch_R30; R_disch_res.R60=disch_R60;
P_disch_res.soc=disch_soc; P_disch_res.P2=disch_P2; P_disch_res.P10=disch_P10; P_disch_res.P30=disch_P30; P_disch_res.P60=disch_P60;

% Fill R_series array using Zero-order hold method
R_current = 0.0015;
if ~isempty(ch_R60) && ~isempty(disch_R60)
    R_current = mean([ch_R60(1), disch_R60(1)]);
end

R_series = zeros(size(t));
for i = 1:length(t)
    if ~isnan(R60_interp(i))
        R_current = R60_interp(i);
    end
    R_series(i) = R_current;
end
end

function rmse_val = obj_therm_model(x, t, I, V, V_ocv, T_meas, T_amb)
C_cell = x(1);
R_th = x(2);
T_pred = simulate_thermal_model(t, I, V, V_ocv, T_meas(1), T_amb, C_cell, R_th);

% 목적 함수 (Objective Function): RMSE 최소화 (전기적 모델과 동일한 방식 적용)
% N_meas 에 대해 (T_meas - T_pred)^2 의 평균을 구한 뒤 루트(sqrt) 적용
rmse_val = sqrt(mean((T_meas - T_pred).^2));
end

% Lumped Thermal Model Simulation
function T_pred = simulate_thermal_model(t, I, V, V_ocv, T_initial, T_amb, C_cell, R_th)
N = length(t);
T_pred = zeros(N, 1);
T_pred(1) = T_initial;

for i = 1:N-1
    dt = t(i+1) - t(i);
    if dt == 0, dt = 1e-4; end % 강제 타임 루프 오류 방지 (데이터 로깅 오차 대비)

    % Overpotential 기반의 실제 열역학적 발열량 계산 (사용자 요청: Q_gen = Vocv - Vcell 결합 모델)
    Q_gen = abs(I(i) * (V_ocv(i) - V(i)));

    % [수정됨] Forward Euler (dt 발산) 문제를 해결하기 위해 수학적으로 완벽하게 안정적인
    % Exact Analytical Solution (해석적 해/지수 함수)을 적용합니다. 
    % T(t) = T_ss + (T_initial - T_ss) * exp(-dt / tau)
    tau = C_cell * R_th;                  % 열 시정수 (Thermal Time Constant)
    T_ss = T_amb + Q_gen * R_th;          % 해당 발열량에서의 도달 목표 정상상태 온도 (Steady State)
    
    T_pred(i+1) = T_ss + (T_pred(i) - T_ss) * exp(-dt / tau);
end
end
