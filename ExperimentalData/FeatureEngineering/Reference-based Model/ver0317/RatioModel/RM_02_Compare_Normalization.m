% RM_02_Compare_Normalization.m
% 비교 실험: 동일 전처리에서 3가지 dQ 정규화 방식 비교
% A: Raw dQ (Ah)
% B: dQ / sum(valid dQ)  — 유효 세그먼트 합 대비 비율
% C: dQ / Q_total         — 해당 사이클 총 충전량 대비 비율
%
% 전처리: 모든 방법에 동일 (no smoothing, monotonic + interp1(NaN))
% 모델: GBM만 (가장 안정적 + 비교 명확)

clear; clc; close all;

%% 1. Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
rmDir   = fullfile(verDir, 'RatioModel');
visDir  = fullfile(rmDir, 'Visualization');
rawDir  = fullfile(projDir, 'Rack_raw2mat');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% 2. Load Lab Data
d = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat'));
FM = d.FM;
mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;
num_segs = 12; Q_nom = 64;

%% 3. Prepare Lab features (3 versions)
dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);
dQ_c = FM{:, dQ_c_cols};  dQ_d = FM{:, dQ_d_cols};  % raw
C_eff_c = FM.C_eff_chg;  C_eff_d = FM.C_eff_dch;

% Q_total for each row (sum of ALL dQ segments charge + discharge)
Q_total_c_lab = sum(FM{:, arrayfun(@(i) sprintf('dQ_c_%02d',i), 1:12, 'UniformOutput', false)}, 2, 'omitnan');
Q_total_d_lab = sum(FM{:, arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:12, 'UniformOutput', false)}, 2, 'omitnan');

% A: Raw dQ
X_A = [dQ_c, dQ_d, C_eff_c, C_eff_d];

% B: dQ / sum(valid subset)
dQ_c_B = dQ_c ./ sum(dQ_c, 2, 'omitnan');
dQ_d_B = dQ_d ./ sum(dQ_d, 2, 'omitnan');
X_B = [dQ_c_B, dQ_d_B, C_eff_c, C_eff_d];

% C: dQ / Q_total
dQ_c_C = dQ_c ./ Q_total_c_lab;
dQ_d_C = dQ_d ./ Q_total_d_lab;
X_C = [dQ_c_C, dQ_d_C, C_eff_c, C_eff_d];

y_lab = FM.Static_Capacity / Q_nom * 100;
methods = {'A_RawDQ', 'B_RatioSum', 'C_RatioQtotal'};
X_all = {X_A, X_B, X_C};

%% 4. Train GBM on full Lab data (3 versions)
fprintf('=== Training 3 GBM models ===\n');
Models = struct();
MU = struct(); SIG = struct();

for m = 1:3
    X = X_all{m};
    mu = mean(X, 1, 'omitnan'); sig = std(X, 0, 1, 'omitnan'); sig(sig==0)=1;
    Xs = (X - mu) ./ sig;
    
    mdl = fitrensemble(Xs, y_lab, 'Method', 'LSBoost', 'NumLearningCycles', 100, ...
        'LearnRate', 0.1, 'Learners', templateTree('MaxNumSplits', 10, 'MinLeafSize', 5, 'Surrogate', 'on'));
    
    Models.(methods{m}) = mdl;
    MU.(methods{m}) = mu;
    SIG.(methods{m}) = sig;
    fprintf('  %s trained.\n', methods{m});
end

%% 5. Field Data
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old', datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new', datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new', datetime(2024,9,9);
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new', datetime(2025,7,11);
};

Np = 2; thr_A = Q_nom * 0.05 / Np;
min_chg_sec = [600, 300, 300, 300];
min_dch_sec = [300, 150, 300, 150];
FieldComp = struct();

for k = 1:size(dataFiles,1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3}; base_date=dataFiles{k,4};
    if ~exist(fpath,'file'), continue; end
    fprintf('\n=== %s ===\n', yr);

    S = load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'), t=datetime(D.Time);
        elseif isduration(D.Date_Time), t=base_date+D.Date_Time;
        else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
        if isfield(D,'SOHPct'), raw_soh=D.SOHPct(:); else, raw_soh=nan(size(I_rack)); end
    else
        D=S.Raw;
        if isduration(D.Date_Time), t=base_date+D.Date_Time; else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
        if isfield(D,'SOH_BMS'), raw_soh=D.SOH_BMS(:); else, raw_soh=nan(size(I_rack)); end
    end

    tsec=seconds(t-t(1)); I_cell=I_rack/Np;
    valid_soh=raw_soh; valid_soh(valid_soh<=0)=NaN;
    soh_bms = median(valid_soh, 'omitnan');

    chgSegs=local_find_segments(I_cell>thr_A);
    dchSegs=local_find_segments(I_cell<-thr_A);
    if ~isempty(chgSegs), dur=chgSegs(:,2)-chgSegs(:,1)+1; chgSegs=chgSegs(dur>=min_chg_sec(k),:); end
    if ~isempty(dchSegs), dur=dchSegs(:,2)-dchSegs(:,1)+1; dchSegs=dchSegs(dur>=min_dch_sec(k),:); end
    if isempty(chgSegs), continue; end

    % --- Charge ---
    [~,bi]=max(chgSegs(:,2)-chgSegs(:,1));
    cs=chgSegs(bi,1); ce=chgSegs(bi,2);
    v_c=V_avg(cs:ce); i_c=I_cell(cs:ce); t_c=tsec(cs:ce);
    q_c = cumtrapz(t_c, abs(i_c)) / 3600;
    Q_total_c = q_c(end);  % total charge of this event
    v_mono=v_c(1); q_mono=q_c(1);
    for ii=2:length(v_c)
        if v_c(ii)>v_mono(end), v_mono(end+1)=v_c(ii); q_mono(end+1)=q_c(ii); end
    end
    dQ_chg = nan(1, num_segs);
    if length(v_mono)>1
        QR = interp1(v_mono, q_mono, V_chg, 'linear', NaN);
        dQ_chg = abs(diff(QR));
    end

    % --- Discharge ---
    dQ_dch = nan(1, num_segs); Q_total_d = NaN;
    C_eff_d_val = NaN;
    if ~isempty(dchSegs)
        [~,bi_d]=max(dchSegs(:,2)-dchSegs(:,1));
        ds=dchSegs(bi_d,1); de=dchSegs(bi_d,2);
        v_d=V_avg(ds:de); i_d=I_cell(ds:de); t_d=tsec(ds:de);
        q_d = cumtrapz(t_d, abs(i_d)) / 3600;
        Q_total_d = q_d(end);
        v_mono_d=v_d(1); q_mono_d=q_d(1);
        for ii=2:length(v_d)
            if v_d(ii)<v_mono_d(end), v_mono_d(end+1)=v_d(ii); q_mono_d(end+1)=q_d(ii); end
        end
        v_asc=flip(v_mono_d); q_asc=flip(q_mono_d);
        if length(v_asc)>1
            QR_d = interp1(v_asc, q_asc, V_dch, 'linear', NaN);
            dQ_dch = abs(diff(QR_d));
        end
        C_eff_d_val = mean(abs(i_d)) / Q_nom;
    end
    C_eff_c_val = mean(abs(i_c)) / Q_nom;

    % Extract subsets
    dQ_c_f = dQ_chg(3:12);
    dQ_d_f = dQ_dch(1:11);

    % A: Raw dQ
    X_fA = [dQ_c_f, dQ_d_f, C_eff_c_val, C_eff_d_val];

    % B: dQ / sum(valid)
    dQ_c_fB = dQ_c_f / sum(dQ_c_f, 'omitnan');
    dQ_d_fB = dQ_d_f / sum(dQ_d_f, 'omitnan');
    X_fB = [dQ_c_fB, dQ_d_fB, C_eff_c_val, C_eff_d_val];

    % C: dQ / Q_total
    dQ_c_fC = dQ_c_f / Q_total_c;
    dQ_d_fC = dQ_d_f / Q_total_d;
    X_fC = [dQ_c_fC, dQ_d_fC, C_eff_c_val, C_eff_d_val];

    X_fields = {X_fA, X_fB, X_fC};

    fprintf('  BMS=%.1f%%, Q_total_c=%.1f, Q_total_d=%.1f\n', soh_bms, Q_total_c, Q_total_d);
    fprintf('  %-15s  SOH%%\n', 'Method');
    
    FieldComp.(yr).BMS = soh_bms;
    for m = 1:3
        Xf = X_fields{m};
        Xf_s = (Xf - MU.(methods{m})) ./ SIG.(methods{m});
        Xf_imp = Xf_s; Xf_imp(isnan(Xf_imp)) = 0;
        soh = predict(Models.(methods{m}), Xf_imp);
        FieldComp.(yr).(methods{m}) = soh;
        fprintf('  %-15s  %.1f\n', methods{m}, soh);
    end
end

%% 6. Summary Table
fprintf('\n\n========== NORMALIZATION COMPARISON (GBM, SOH%%) ==========\n');
fprintf('%-6s  %-6s  %-10s %-12s %-12s\n', 'Year', 'BMS', 'A:RawDQ', 'B:dQ/sum', 'C:dQ/Qtotal');
yrs = fieldnames(FieldComp);
for k = 1:length(yrs)
    r = FieldComp.(yrs{k});
    fprintf('%-6s  %-6.1f  %-10.1f %-12.1f %-12.1f\n', ...
        yrs{k}, r.BMS, r.A_RawDQ, r.B_RatioSum, r.C_RatioQtotal);
end

%% 7. Visualization
fig = figure('Position', [50 100 1000 500]);
n = length(yrs); x = 1:n;
bms_arr = arrayfun(@(i) FieldComp.(yrs{i}).BMS, 1:n);
A_arr = arrayfun(@(i) FieldComp.(yrs{i}).A_RawDQ, 1:n);
B_arr = arrayfun(@(i) FieldComp.(yrs{i}).B_RatioSum, 1:n);
C_arr = arrayfun(@(i) FieldComp.(yrs{i}).C_RatioQtotal, 1:n);

hold on; grid on; box on;
plot(x, bms_arr, 'k--s', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'BMS SOH');
plot(x, A_arr, '-o', 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.2 0.2], 'DisplayName', 'A: Raw dQ');
plot(x, B_arr, '-d', 'Color', [0.2 0.6 0.9], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.6 0.9], 'DisplayName', 'B: dQ/sum(valid)');
plot(x, C_arr, '-^', 'Color', [0.1 0.7 0.3], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.1 0.7 0.3], 'DisplayName', 'C: dQ/Q_{total}');

% Annotate
for i = 1:n
    text(x(i)+0.05, bms_arr(i)+0.3, sprintf('%.1f', bms_arr(i)), 'FontSize', 8, 'FontWeight', 'bold');
    text(x(i)+0.05, A_arr(i)-0.8, sprintf('%.1f', A_arr(i)), 'FontSize', 8, 'Color', [0.8 0.2 0.2]);
    text(x(i)+0.05, B_arr(i)+0.3, sprintf('%.1f', B_arr(i)), 'FontSize', 8, 'Color', [0.2 0.6 0.9]);
    text(x(i)+0.05, C_arr(i)-0.8, sprintf('%.1f', C_arr(i)), 'FontSize', 8, 'Color', [0.1 0.7 0.3]);
end

set(gca, 'XTick', x, 'XTickLabel', strrep(yrs, 'Y', ''), 'FontSize', 12);
xlabel('Year', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('SOH (%)', 'FontSize', 13, 'FontWeight', 'bold');
title('GBM Field SOH: Raw dQ vs dQ/sum vs dQ/Q_{total}', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 11);
ylim([80 105]);

saveas(fig, fullfile(visDir, 'RM_Normalization_Comparison.png'));
fprintf('\nSaved: RM_Normalization_Comparison.png\n');
fprintf('=== Comparison Complete! ===\n');

%% Helper
function segs = local_find_segments(mask)
    segs=[]; n=length(mask); i=1;
    while i<=n
        if mask(i), j=i;
            while j<n && mask(j+1), j=j+1; end
            segs=[segs;i,j]; i=j+1;
        else, i=i+1; end
    end
end
