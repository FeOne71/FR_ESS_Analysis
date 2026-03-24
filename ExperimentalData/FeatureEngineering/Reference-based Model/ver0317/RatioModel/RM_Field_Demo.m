% RM_Field_Demo.m
% 필드 시연: 충전/방전 독립 SOH 추정 (마스킹 증강 모델)
% 사전 실행 필요: RM_Charge_Model.m, RM_Discharge_Model.m

clear; clc; close all;

%% 1. Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
rmDir   = fullfile(verDir, 'RatioModel');
visDir  = fullfile(rmDir, 'Visualization');
rawDir  = fullfile(projDir, 'Rack_raw2mat');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;

%% 2. Load trained models
chg_data = load(fullfile(rmDir, 'RM_Charge_Final.mat'));
dch_data = load(fullfile(rmDir, 'RM_Discharge_Final.mat'));

%% 3. Field Data
Q_nom = 64; Np = 2; thr_A = Q_nom * 0.05 / Np;
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old', datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new', datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new', datetime(2024,9,9);
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new', datetime(2025,7,11);
};
min_chg_sec = [600, 300, 300, 300];
min_dch_sec = [300, 150, 300, 150];
FieldResults = struct();

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

    %% --- CHARGE ---
    soh_chg_xgb=NaN; soh_chg_lgb=NaN; C_eff_c_val=NaN;
    if ~isempty(chgSegs)
        [~,bi]=max(chgSegs(:,2)-chgSegs(:,1));
        cs=chgSegs(bi,1); ce=chgSegs(bi,2);
        v_c=V_avg(cs:ce); i_c=I_cell(cs:ce); t_c=tsec(cs:ce);
        q_c=cumtrapz(t_c,abs(i_c))/3600;
        v_mono=v_c(1); q_mono=q_c(1);
        for ii=2:length(v_c), if v_c(ii)>v_mono(end), v_mono(end+1)=v_c(ii); q_mono(end+1)=q_c(ii); end; end
        dQ_chg=nan(1,12);
        if length(v_mono)>1, QR=interp1(v_mono,q_mono,V_chg,'linear',NaN); dQ_chg=abs(diff(QR)); end
        C_eff_c_val = mean(abs(i_c))/Q_nom;

        dQ_c_f = dQ_chg(3:12);  % 10 features
        dQ_c_ratio = dQ_c_f / sum(dQ_c_f, 'omitnan');
        X_chg = [dQ_c_ratio, C_eff_c_val];
        X_chg_s = (X_chg - chg_data.mu_chg) ./ chg_data.sig_chg;

        X_py = py.numpy.array(X_chg_s).reshape(int32(1),int32(-1));
        soh_chg_xgb = double(chg_data.FinalModels.XGB.predict(X_py));
        soh_chg_lgb = double(chg_data.FinalModels.LGB.predict(X_py));

        fprintf('  CHG: Segs=%s, C_eff=%.3f\n', mat2str(find(~isnan(dQ_c_f))), C_eff_c_val);
        fprintf('       XGB=%.1f%%, LGB=%.1f%%\n', soh_chg_xgb, soh_chg_lgb);
    end

    %% --- DISCHARGE ---
    soh_dch_xgb=NaN; soh_dch_lgb=NaN; C_eff_d_val=NaN;
    if ~isempty(dchSegs)
        [~,bi_d]=max(dchSegs(:,2)-dchSegs(:,1));
        ds=dchSegs(bi_d,1); de=dchSegs(bi_d,2);
        v_d=V_avg(ds:de); i_d=I_cell(ds:de); t_d=tsec(ds:de);
        q_d=cumtrapz(t_d,abs(i_d))/3600;
        v_md=v_d(1); q_md=q_d(1);
        for ii=2:length(v_d), if v_d(ii)<v_md(end), v_md(end+1)=v_d(ii); q_md(end+1)=q_d(ii); end; end
        v_asc=flip(v_md); q_asc=flip(q_md);
        dQ_dch=nan(1,12);
        if length(v_asc)>1, QR_d=interp1(v_asc,q_asc,V_dch,'linear',NaN); dQ_dch=abs(diff(QR_d)); end
        C_eff_d_val = mean(abs(i_d))/Q_nom;

        dQ_d_f = dQ_dch(1:11);  % 11 features
        dQ_d_ratio = dQ_d_f / sum(dQ_d_f, 'omitnan');
        X_dch = [dQ_d_ratio, C_eff_d_val];
        X_dch_s = (X_dch - dch_data.mu_dch) ./ dch_data.sig_dch;

        X_py = py.numpy.array(X_dch_s).reshape(int32(1),int32(-1));
        soh_dch_xgb = double(dch_data.FinalModels.XGB.predict(X_py));
        soh_dch_lgb = double(dch_data.FinalModels.LGB.predict(X_py));

        fprintf('  DCH: Segs=%s, C_eff=%.3f\n', mat2str(find(~isnan(dQ_d_f))), C_eff_d_val);
        fprintf('       XGB=%.1f%%, LGB=%.1f%%\n', soh_dch_xgb, soh_dch_lgb);
    end

    FieldResults.(yr).BMS = soh_bms;
    FieldResults.(yr).CHG_XGB = soh_chg_xgb;
    FieldResults.(yr).CHG_LGB = soh_chg_lgb;
    FieldResults.(yr).DCH_XGB = soh_dch_xgb;
    FieldResults.(yr).DCH_LGB = soh_dch_lgb;
    fprintf('  BMS=%.1f%%\n', soh_bms);
end

%% 4. Summary
fprintf('\n\n========== FIELD SOH: CHARGE vs DISCHARGE (Masking-Augmented) ==========\n');
fprintf('%-6s  %-5s  %-9s %-9s  %-9s %-9s\n', 'Year','BMS','CHG_XGB','CHG_LGB','DCH_XGB','DCH_LGB');
yrs = fieldnames(FieldResults);
for k = 1:length(yrs)
    r = FieldResults.(yrs{k});
    fprintf('%-6s  %-5.1f  %-9.1f %-9.1f  %-9.1f %-9.1f\n', ...
        yrs{k}, r.BMS, r.CHG_XGB, r.CHG_LGB, r.DCH_XGB, r.DCH_LGB);
end

%% 5. Visualization
fig = figure('Position', [50 100 1100 500]);
n = length(yrs); x = 1:n;
bms = arrayfun(@(i) FieldResults.(yrs{i}).BMS, 1:n);
cxgb = arrayfun(@(i) FieldResults.(yrs{i}).CHG_XGB, 1:n);
clgb = arrayfun(@(i) FieldResults.(yrs{i}).CHG_LGB, 1:n);
dxgb = arrayfun(@(i) FieldResults.(yrs{i}).DCH_XGB, 1:n);
dlgb = arrayfun(@(i) FieldResults.(yrs{i}).DCH_LGB, 1:n);

hold on; grid on; box on;
plot(x, bms, 'k--s', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'BMS SOH');
plot(x, cxgb, '-o', 'Color', [0.8 0.3 0.1], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.3 0.1], 'DisplayName', 'CHG XGBoost');
plot(x, clgb, '-d', 'Color', [0.9 0.6 0.2], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.9 0.6 0.2], 'DisplayName', 'CHG LightGBM');
plot(x, dxgb, '-^', 'Color', [0.1 0.5 0.8], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.1 0.5 0.8], 'DisplayName', 'DCH XGBoost');
plot(x, dlgb, '-v', 'Color', [0.3 0.7 0.9], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.7 0.9], 'DisplayName', 'DCH LightGBM');

for i=1:n
    text(x(i)+0.05, bms(i)+0.4, sprintf('%.1f', bms(i)), 'FontSize', 8, 'FontWeight', 'bold');
    text(x(i)+0.05, cxgb(i)-0.6, sprintf('%.1f', cxgb(i)), 'FontSize', 7, 'Color', [0.8 0.3 0.1]);
    text(x(i)-0.15, clgb(i)+0.4, sprintf('%.1f', clgb(i)), 'FontSize', 7, 'Color', [0.9 0.6 0.2]);
    text(x(i)+0.05, dxgb(i)-0.6, sprintf('%.1f', dxgb(i)), 'FontSize', 7, 'Color', [0.1 0.5 0.8]);
    text(x(i)-0.15, dlgb(i)+0.4, sprintf('%.1f', dlgb(i)), 'FontSize', 7, 'Color', [0.3 0.7 0.9]);
end

set(gca, 'XTick', x, 'XTickLabel', strrep(yrs, 'Y', ''), 'FontSize', 12);
xlabel('Year', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('SOH (%)', 'FontSize', 13, 'FontWeight', 'bold');
title('Masking-Augmented Model: Charge vs Discharge SOH', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 10);
ylim([80 105]);
saveas(fig, fullfile(visDir, 'RM_Field_ChgDch_Comparison.png'));
fprintf('\nSaved: RM_Field_ChgDch_Comparison.png\n');

save(fullfile(rmDir, 'RM_Field_Results.mat'), 'FieldResults');
fprintf('=== Field Demo Complete! ===\n');

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
