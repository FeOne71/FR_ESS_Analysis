%% Field_Estimation_v6_Charge.m
% =====================================================================
% Phase 6: 현장 데이터 SOH 추정 (Fragmented Charge)
% - 파편화된 현장 데이터에서 v5 (11개) 전압 구간별 dQ 추출
% - 관측되지 않은 구간(NaN)은 Lab 데이터의 평균(Mean)으로 Imputation 처리
%   (GPR/SVM 모델이 NaN 입력 시 예측 불가 방지)
% - ML_Models_v6_Charge 내의 'All', 'Top7', 'Top4' 모델을 모두 적용해보고 결과 비교
% =====================================================================
clear; clc; close all;

%% 1. 모델 및 레퍼런스 데이터 로드
base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                    'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02', 'ML_RandomForest_v4');
mdl_path = fullfile(base_dir, 'ML_Models_v6_Charge.mat');
feat_path = fullfile(base_dir, 'RPT_FeatureLabelExtractor_v5_Charge', 'Feature_Matrix_v5_Charge.mat');

load(mdl_path, 'Final_Models');
% Lab Feature의 Mean 계산 (Imputation 용도)
feat_data = load(feat_path);
X_lab = table2array(feat_data.FeatureTable_v5_Charge(:, 5:15));
mean_X_lab = mean(X_lab, 1, 'omitnan');

% OCV-based Reference SOH
ref_SOH = struct('Y2021', 97.5, 'Y2023', 95.0, 'Y2024', 92.5, 'Y2025', 90.0);

rawDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
files = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old'
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new'
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new'
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new'
};

Np = 2; Q_cell = 64; 
% v5 Equi-Capacity Boundaries
VB = [ 3.40, 3.52, 3.59, 3.64, 3.67, 3.70, 3.75, 3.84, 3.92, 4.01, 4.10, 4.20 ];
N_seg = length(VB)-1;
thr_A = Q_cell * 0.05;

%% 2. 연도별 필드 추정 루프
fprintf('=== Field SOH Estimation (v6 Fragmented Charge) ===\n');
res = struct();
tracks_to_test = {'All', 'Top7'};
mtype = 'GPR';

for k = 1:size(files,1)
    fpath = files{k,1}; yr = files{k,2}; dtype = files{k,3};
    if ~exist(fpath,'file'), fprintf('[Skip] %s not found\n', yr); continue; end
    
    S = load(fpath);
    if strcmp(dtype,'old')
        I_rack = S.Raw.Rack01.DCCurrent_A(:);
        time_arr = S.Raw.Rack01.Time;
        V_cell = S.Raw.Rack01.AverageCV_V(:);
        dt = datetime(time_arr, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        t_s = seconds(dt - dt(1));
    else
        I_rack = S.Raw.DCCurrent(:);
        time_arr = S.Raw.Date_Time;
        V_cell = S.Raw.CVavg(:);
        t_s = seconds(time_arr) - seconds(time_arr(1));
    end
    I_cell = I_rack / Np;
    if isempty(t_s), t_s = (0:length(I_cell)-1)'; end
    
    % ---- 충전 세그먼트 탐색 ----
    chg_mask = I_cell > thr_A;
    d = diff([0; chg_mask; 0]);
    starts = find(d==1); ends = find(d==-1)-1;
    lens = ends - starts + 1;
    if isempty(lens), continue; end
    
    % 제일 긴 충전 구간 선택 (파편화 가정)
    [~, best_idx] = max(lens);
    idx = starts(best_idx):ends(best_idx);
    
    t_seg = t_s(idx); V_seg = V_cell(idx); I_seg = I_cell(idx);
    Q_seg = cumtrapz(t_seg, I_seg)/3600;
    
    % C_eff
    C_eff = mean(abs(I_seg)) / Q_cell;
    
    % dQ 추출 (11개)
    dQ_raw = extract_dQ(V_seg, Q_seg, VB, true);
    
    % 관측 안된 구간(NaN) 개수 확인
    num_nan = sum(isnan(dQ_raw));
    
    % Mean Imputation
    dQ_imp = dQ_raw;
    for s = 1:N_seg
        if isnan(dQ_imp(s))
            dQ_imp(s) = mean_X_lab(s);
        end
    end
    
    fprintf('-- %s --\n', yr);
    fprintf('  V_range: %.3f ~ %.3f V (Missing Segs: %d/%d)\n', min(V_seg), max(V_seg), num_nan, N_seg);
    fprintf('  C_eff:   %.3f C\n', C_eff);
    
    res.(yr).C_eff = C_eff;
    res.(yr).num_nan = num_nan;
    
    % 모델별 평가
    for tidx = 1:length(tracks_to_test)
        trk = tracks_to_test{tidx};
        segs = Final_Models.(trk).(mtype).features;
        feat = [dQ_imp(segs), C_eff];
        mdl = Final_Models.(trk).(mtype).model;
        
        pred_SOH = predict(mdl, feat);
        res.(yr).(trk) = pred_SOH;
        fprintf('  [%s] %s SOH: %.2f%%\n', trk, mtype, pred_SOH);
    end
end

%% 3. 결과 시각화
yrs_valid = fieldnames(res);
n_yrs = length(yrs_valid);
soh_ref = nan(1, n_yrs);
soh_all = nan(1, n_yrs);
soh_top7 = nan(1, n_yrs);

for i = 1:n_yrs
    yr = yrs_valid{i};
    if isfield(ref_SOH, yr), soh_ref(i) = ref_SOH.(yr); end
    if isfield(res.(yr), 'All'), soh_all(i) = res.(yr).All; end
    if isfield(res.(yr), 'Top7'), soh_top7(i) = res.(yr).Top7; end
end

fig = figure('Position',[100,200,800,400],'Name','Field SOH Estimation v6');

% SOH 비교 플롯
hold on; grid on;
b = bar(1:n_yrs, [soh_ref; soh_top7; soh_all]');
b(1).FaceColor = [0.4 0.4 0.4]; b(1).DisplayName = 'Reference (OCV-based)';
b(2).FaceColor = [0.0 0.45 0.74]; b(2).DisplayName = sprintf('Model Top7 (%s)', mtype);
b(3).FaceColor = [0.85 0.32 0.09]; b(3).DisplayName = sprintf('Model All 11 (%s)', mtype);

% Text
for i=1:n_yrs
    yr = yrs_valid{i};
    text(i, soh_all(i)+0.5, sprintf('NaN:%d', res.(yr).num_nan), 'HorizontalAlignment','center','FontSize',8,'Color','k');
end

xticks(1:n_yrs); xticklabels(yrs_valid);
ylabel('SOH (%)'); ylim([80 100]);
title('Field Data SOH Estimation (v6 Fragmented + Imputation)','FontWeight','bold');
legend('Location','southwest');

saveas(fig, fullfile(base_dir, 'Field_Estimation_Results_v6.fig'));
saveas(fig, fullfile(base_dir, 'Field_Estimation_Results_v6.png'));
fprintf('\n>> Saved: Field_Estimation_Results_v6.fig / .png\n');


%% Helper
function dQ = extract_dQ(V, Q, VB, is_charge)
    N = length(VB)-1; dQ = nan(1,N);
    [Vu, uid] = unique(V, 'stable'); Qu = Q(uid);
    mono = true(size(Vu));
    if is_charge
        for ii=2:length(Vu), if Vu(ii)<=Vu(ii-1), mono(ii)=false; end; end
    else
        for ii=2:length(Vu), if Vu(ii)>=Vu(ii-1), mono(ii)=false; end; end
    end
    Vu = Vu(mono); Qu = Qu(mono);
    if isempty(Vu), return; end
    if Vu(1)>Vu(end), Vu=flipud(Vu); Qu=flipud(Qu); end
    if length(Vu)<10, return; end
    
    for s=1:N
        if VB(s)>=min(Vu) && VB(s+1)<=max(Vu)
            q1=interp1(Vu,Qu,VB(s),'linear');
            q2=interp1(Vu,Qu,VB(s+1),'linear');
            if ~isnan(q1)&&~isnan(q2), dQ(s)=abs(q2-q1); end
        end
    end
end
