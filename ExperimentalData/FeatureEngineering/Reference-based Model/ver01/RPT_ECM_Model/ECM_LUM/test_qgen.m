% test_qgen.m
clear; clc;

ch_str = 'ch09';
Ch_str_cap = 'Ch09';
dataFile = sprintf('D:\\JCW\\Projects\\KEPCO_ESS_Local\\ExperimentalData\\RPT\\Postprocessing\\Parsed\\RPT0_%s_parsed.mat', ch_str);
loaded_data = load(dataFile);
pdata = loaded_data.pdata;

ocvFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat';
ocv_data = load(ocvFile);
ocv_struct = ocv_data.RPT_VQ_grid.cyc0.(Ch_str_cap);

C_nom = 64.0;
I_1C = 64.0;
tolC = 0.05;

disch_n1C_idx = []; disch_c02_idx = [];
for i = 1:length(pdata)
    if isempty(pdata(i).I), continue; end
    avg_I = mean(pdata(i).I); dur = pdata(i).t(end) - pdata(i).t(1);
    if avg_I <= -I_1C*(1-tolC) && avg_I >= -I_1C*(1+tolC) && dur >= 50 && dur <= 70, disch_n1C_idx(end+1) = i; end
    if avg_I <= -12.8*(1-tolC) && avg_I >= -12.8*(1+tolC), disch_c02_idx(end+1) = i; end
end
dchg_1C_idx = disch_n1C_idx;
if ~isempty(disch_c02_idx), dchg_02C_idx = disch_c02_idx; else dchg_02C_idx = disch_n1C_idx(end); end

charge_n1C_idx = []; charge_c02_idx = [];
for i = 1:length(pdata)
    if isempty(pdata(i).I), continue; end
    avg_I = mean(pdata(i).I); dur = pdata(i).t(end) - pdata(i).t(1);
    if avg_I >= I_1C*(1-tolC) && avg_I <= I_1C*(1+tolC) && dur >= 50 && dur <= 70, charge_n1C_idx(end+1) = i; end
    if avg_I >= 12.8*(1-tolC) && avg_I <= 12.8*(1+tolC), charge_c02_idx(end+1) = i; end
end
chg_1C_idx = charge_n1C_idx;
if ~isempty(charge_c02_idx), chg_02C_idx = charge_c02_idx; else chg_02C_idx = charge_n1C_idx(end); end

all_indices = unique([max(1, chg_1C_idx(1)-1):chg_02C_idx(end), max(1, dchg_1C_idx(1)-1):dchg_02C_idx(end)]);
full_valid_idx = all_indices(all_indices <= length(pdata));

pdata_filtered = pdata(full_valid_idx);
I_all = vertcat(pdata_filtered.I);
V_all = vertcat(pdata_filtered.V);
if isfield(pdata, 'steptime_double') && ~isempty(pdata(1).steptime_double)
    steptime_concat = vertcat(pdata_filtered.steptime_double);
    steptime_concat(1) = 0;
    t_all = zeros(size(steptime_concat));
    for i=2:length(t_all)
        dt_val = steptime_concat(i) - steptime_concat(i-1);
        if dt_val < 0, t_all(i) = t_all(i-1) + 0.1; else t_all(i) = t_all(i-1) + dt_val; end
    end
else
    t_all = vertcat(pdata_filtered.t); t_all = t_all - t_all(1);
end
SOC_all = cumtrapz(t_all, I_all) / (3600 * C_nom);

V_ocv_all = zeros(size(t_all));
for i=1:length(t_all)
    if I_all(i) > 0.1
        q_grid = ocv_struct.OCV_charge.Q; v_grid = ocv_struct.OCV_charge.V_grid;
        current_Q = SOC_all(i) * max(q_grid);
    else
        q_grid = ocv_struct.OCV_discharge.Q; v_grid = ocv_struct.OCV_discharge.V_grid;
        current_Q = (1.0 - SOC_all(i)) * max(q_grid);
    end
    [q_uniq, idx_uniq] = unique(q_grid); v_uniq = v_grid(idx_uniq);
    if length(q_uniq) > 1, V_ocv_all(i) = interp1(q_uniq, v_uniq, current_Q, 'linear', 'extrap'); else V_ocv_all(i) = v_uniq(1); end
end

Q_gen = abs(I_all .* (V_ocv_all - V_all));
Total_Q = trapz(t_all, Q_gen);

fprintf('Max V_ocv = %.3f, Min V_ocv = %.3f\n', max(V_ocv_all), min(V_ocv_all));
fprintf('Max V_all = %.3f, Min V_all = %.3f\n', max(V_all), min(V_all));
fprintf('Max Q_gen (W) = %.3f, Mean Q_gen = %.3f\n', max(Q_gen), mean(Q_gen));
fprintf('Total Heat Energy Generated (J)= %.2f\n', Total_Q);
fprintf('Total Temp Rise (no cooling, C_cell=1000) = %.2f\n', Total_Q / 1000.0);

% Find where max Q_gen occurs
[max_q, max_idx] = max(Q_gen);
fprintf('Max Q occurs at I = %.2f A, V = %.3f V, V_ocv = %.3f V, SOC = %.3f\n', I_all(max_idx), V_all(max_idx), V_ocv_all(max_idx), SOC_all(max_idx));
