% DCIR_2RC_Fitting_Comparison.m
% =====================================================================
% 비교: 펄스만 피팅 vs 펄스+휴지기 피팅
% Ch09 cyc0의 첫 5개 펄스로 비교
% =====================================================================
clear; clc; close all;

parsedFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed\RPT0_ch09_parsed.mat';
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4';

d = load(parsedFile);
pdata = d.pdata;
I_1C = 64; tol = 0.05;

% Find discharge 1C pulses
pulse_idx = [];
for k = 1:length(pdata)
    avg_I = mean(pdata(k).I);
    dur = pdata(k).t(end) - pdata(k).t(1);
    if avg_I < -(1-tol)*I_1C && avg_I > -(1+tol)*I_1C && dur > 5 && dur < 70
        pulse_idx(end+1) = k; %#ok<SAGROW>
    end
end

d_dcir = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_Final\DCIR_SOC_data_all_channels_final.mat');
soc_ref = d_dcir.dcir_soc_data.Ch09.cyc0.discharge_table.SOC;

% 2RC model functions
pulse_model = @(p, t) p(1) + p(2)*(1-exp(-t/p(3))) + p(4)*(1-exp(-t/p(5)));
relax_model = @(p, t) p(2)*exp(-t/p(3)) + p(4)*exp(-t/p(5));  % R₀ drops instantly, RC decays

opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 2000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);

n_show = min(6, length(pulse_idx));  % 비교할 펄스 수

% Storage
P_pulse = zeros(n_show, 5);   % pulse-only params
P_both  = zeros(n_show, 5);   % pulse+relax params

fig = figure('Name', '2RC Fitting Comparison', 'Position', [30, 30, 1600, 800], 'Visible', 'off');

for p = 1:n_show
    k = pulse_idx(p);
    
    % Pulse data
    t_p = pdata(k).t(:) - pdata(k).t(1);
    V_p = pdata(k).V(:);
    I_avg = abs(mean(pdata(k).I));
    
    % V_OCV from rest before pulse
    V_OCV = pdata(k-1).V(end);
    
    % R(t) during pulse
    R_pulse = abs(V_p - V_OCV) / I_avg;
    
    % Relaxation data (rest after pulse)
    if k+1 <= length(pdata) && strcmp(char(pdata(k+1).type), 'R')
        t_r_raw = pdata(k+1).t(:) - pdata(k+1).t(1);
        V_r = pdata(k+1).V(:);
        V_end_pulse = V_p(end);  % voltage at end of pulse
        
        % During relaxation: V recovers from V_end_pulse toward V_OCV
        % The RC voltage = (V_OCV - V_r) = remaining overpotential
        % At t_relax=0: RC voltage = V_OCV - V_end_pulse + R₀*I (R₀ already recovered)
        % Actually: dV_relax(t) = V_OCV - V_r(t) = R₁*exp(-t/τ₁) + R₂*exp(-t/τ₂)
        % (R₀ drops instantly at pulse end)
        
        R_relax = abs(V_OCV - V_r) / I_avg;  % normalized by same I for comparison
        % Limit relaxation data to ~120 seconds
        valid_r = t_r_raw <= 120;
        t_r = t_r_raw(valid_r);
        R_relax = R_relax(valid_r);
    else
        t_r = []; R_relax = [];
    end
    
    % --- Method A: Pulse-only fitting ---
    R0_init = R_pulse(1);
    R_inf = R_pulse(end);
    x0 = [R0_init, (R_inf-R0_init)*0.4, 3, (R_inf-R0_init)*0.6, 30];
    lb = [0, 0, 0.1, 0, 1]; ub = [0.01, 0.01, 20, 0.05, 300];
    
    [pA, ~] = lsqcurvefit(pulse_model, x0, t_p, R_pulse, lb, ub, opts);
    R_fitA = pulse_model(pA, t_p);
    P_pulse(p,:) = pA;
    
    % --- Method B: Pulse + Relaxation combined fitting ---
    if ~isempty(t_r) && length(t_r) > 10
        % Combined cost: fit pulse AND relaxation simultaneously
        % Use fmincon to minimize combined error
        cost_func = @(params) combined_cost(params, t_p, R_pulse, t_r, R_relax);
        
        x0B = pA; % start from pulse-only result
        lbB = [0, 0, 0.1, 0, 1]; ubB = [0.01, 0.01, 20, 0.05, 300];
        opts_fmin = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 2000, ...
            'FunctionTolerance', 1e-14, 'StepTolerance', 1e-14);
        [pB, ~] = fmincon(cost_func, x0B, [], [], [], [], lbB, ubB, [], opts_fmin);
    else
        pB = pA;
    end
    P_both(p,:) = pB;
    R_fitB_pulse = pulse_model(pB, t_p);
    
    % --- Plot comparison ---
    subplot(2, n_show, p);
    plot(t_p, R_pulse*1000, 'k.', 'MarkerSize', 1);
    hold on;
    plot(t_p, R_fitA*1000, 'b-', 'LineWidth', 1.5, 'DisplayName', 'A: Pulse Only');
    plot(t_p, R_fitB_pulse*1000, 'r--', 'LineWidth', 1.5, 'DisplayName', 'B: Pulse+Relax');
    title(sprintf('Pulse SOC=%.0f%%', soc_ref(p)), 'FontSize', 10);
    xlabel('t(s)'); ylabel('R(mΩ)');
    legend('FontSize', 6, 'Location', 'best'); grid on;
    
    % Relaxation plot
    subplot(2, n_show, n_show + p);
    if ~isempty(t_r)
        plot(t_r, R_relax*1000, 'k.', 'MarkerSize', 1);
        hold on;
        R_relaxA = relax_model(pA, t_r);
        R_relaxB = relax_model(pB, t_r);
        plot(t_r, R_relaxA*1000, 'b-', 'LineWidth', 1.5, 'DisplayName', 'A: Pulse Only');
        plot(t_r, R_relaxB*1000, 'r--', 'LineWidth', 1.5, 'DisplayName', 'B: Pulse+Relax');
        title(sprintf('Relax SOC=%.0f%%', soc_ref(p)), 'FontSize', 10);
        xlabel('t(s)'); ylabel('R(mΩ)');
        legend('FontSize', 6, 'Location', 'best'); grid on;
    end
end

sgtitle('2RC Fitting: Pulse Only (Blue) vs Pulse+Relaxation (Red)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig, fullfile(saveDir, 'DCIR_2RC_PulseVsRelax_Comparison.fig'));
close(fig);

% --- Summary table ---
fprintf('\n=== Parameter Comparison ===\n');
fprintf('%-6s | %-40s | %-40s\n', 'SOC', 'A: Pulse Only', 'B: Pulse+Relax');
fprintf('%-6s | %8s %8s %6s %8s %6s | %8s %8s %6s %8s %6s\n', ...
    '', 'R0', 'R1', 'τ1', 'R2', 'τ2', 'R0', 'R1', 'τ1', 'R2', 'τ2');
fprintf('%s\n', repmat('-', 1, 95));
for p = 1:n_show
    fprintf('%5.0f%% | %7.3f %7.3f %5.1f %7.3f %5.1f | %7.3f %7.3f %5.1f %7.3f %5.1f\n', ...
        soc_ref(p), ...
        P_pulse(p,1)*1e3, P_pulse(p,2)*1e3, P_pulse(p,3), P_pulse(p,4)*1e3, P_pulse(p,5), ...
        P_both(p,1)*1e3, P_both(p,2)*1e3, P_both(p,3), P_both(p,4)*1e3, P_both(p,5));
end
fprintf('\nDone.\n');

function cost = combined_cost(params, t_pulse, R_pulse, t_relax, R_relax)
    % Pulse model: R0 + R1*(1-exp(-t/tau1)) + R2*(1-exp(-t/tau2))
    R_pred_pulse = params(1) + params(2)*(1-exp(-t_pulse/params(3))) + params(4)*(1-exp(-t_pulse/params(5)));
    err_pulse = R_pulse - R_pred_pulse;
    
    % Relax model: R1*exp(-t/tau1) + R2*exp(-t/tau2)  (R0 drops instantly)
    R_pred_relax = params(2)*exp(-t_relax/params(3)) + params(4)*exp(-t_relax/params(5));
    err_relax = R_relax - R_pred_relax;
    
    % Combined cost (weighted equally)
    cost = sum(err_pulse.^2) + sum(err_relax.^2);
end
