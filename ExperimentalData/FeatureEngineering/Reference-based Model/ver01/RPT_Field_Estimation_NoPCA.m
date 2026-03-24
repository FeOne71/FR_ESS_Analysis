%% RPT_Field_Estimation_NoPCA.m
% 1. 실험 데이터(Lab Data)로 NoPCA Random Forest 모델 학습
% 2. 필드 데이터(Field Data, 2021~2025) 로드 및 피처 추출 (Master Ruler 적용)
% 3. 학습된 모델을 사용하여 필드 데이터의 SOH, LLI, LAM 추정 및 시각화
% NoPCA 모델은 PCA를 거치지 않은 14개 원본 피처를 직접 사용합니다.

clear; clc; close all;
warning off;

%% ========================================================================
% Section 1: 모델 학습 (실험 데이터 사용)
% ========================================================================
fprintf('=== Section 1: 실험 데이터로 NoPCA_RF 모델 학습 ===\n');

% 1.1 데이터 로드
baseDir = fileparts(mfilename('fullpath')); % 현재 스크립트 위치
labDataPath = fullfile(baseDir, 'Dataset', 'Feature_Matrix_Final.mat');
rulerPath   = fullfile(baseDir, 'Dataset', 'MasterRulers.mat');

if ~exist(labDataPath, 'file')
    error('Lab data file not found: %s', labDataPath);
end
load(labDataPath, 'FeatureTable');
load(rulerPath, 'MasterRulers');

% 1.2 입력/라벨 설정
X_lab = FeatureTable.X_Features; % 정규화 전 원본 피처 (14개)
Y_lab = FeatureTable.Y_Labels;   % [SOH, LLI, LAM]
label_names = {'SOH', 'LLI', 'LAM'};

% NaN 행 제거
nan_rows = any(isnan(X_lab), 2) | any(isnan(Y_lab), 2);
X_lab(nan_rows, :) = [];
Y_lab(nan_rows, :) = [];

% 1.3 정규화 파라미터 계산 (실험 데이터 기준)
mu_X = mean(X_lab, 1, 'omitnan');
sigma_X = std(X_lab, 0, 1, 'omitnan');
sigma_X(sigma_X == 0) = 1;

X_lab_norm = (X_lab - mu_X) ./ sigma_X;

% 1.4 NoPCA RF 모델 학습
% Hyperparameter: 500 trees, 1 min leaf size (RPT_Modeling_NoPCA.m 기준)
RF_Models = cell(1, 3);
rf_trees = 500;
rf_minleaf = 1;

for i = 1:3
    fprintf('  Training NoPCA_RF model for %s...\n', label_names{i});
    t = templateTree('MinLeafSize', rf_minleaf);
    RF_Models{i} = fitrensemble(X_lab_norm, Y_lab(:, i), ...
        'Method', 'Bag', 'NumLearningCycles', rf_trees, 'Learners', t);
end

fprintf('=== 모델 학습 완료 ===\n');


%% ========================================================================
% Section 2: 필드 데이터 로드 및 전처리
% ========================================================================
fprintf('\n=== Section 2: 필드 데이터 로드 (2021~2025) ===\n');

% 경로 설정
rawDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
dataFiles = {
    fullfile(rawDir, 'Old', '2021', '202106', 'Raw_20210603.mat'), 'Y2021';
    fullfile(rawDir, 'New', '2023', '202310', 'Raw_20231016.mat'), 'Y2023';
    fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024';
    fullfile(rawDir, 'New', '2025', '202507', 'Raw_20250711.mat'), 'Y2025'
};

% Master Ruler 설정 (Ch09 기준 - Global Ruler 공유)
ref_ch = 'Ch09'; 
if isfield(MasterRulers, ref_ch)
    V_bounds_chg = MasterRulers.(ref_ch).V_bounds_chg;
    V_bounds_dch = MasterRulers.(ref_ch).V_bounds_dch;
else
    fns = fieldnames(MasterRulers);
    V_bounds_chg = MasterRulers.(fns{1}).V_bounds_chg;
    V_bounds_dch = MasterRulers.(fns{1}).V_bounds_dch;
end

% 필드 데이터 저장용 구조체
FieldData = struct();

for k = 1:size(dataFiles, 1)
    fpath = dataFiles{k, 1};
    year_label = dataFiles{k, 2};
    
    if ~exist(fpath, 'file')
        warning('File not found: %s', fpath);
        continue;
    end
    
    fprintf('  Processing %s...\n', year_label);
    S = load(fpath);
    
    if contains(fpath, 'Old')
        if isfield(S, 'Raw') && isfield(S.Raw, 'Rack01')
            D = S.Raw.Rack01;
            V_avg = D.AverageCV_V(:);
            I_rack = D.DCCurrent_A(:);
        else
            warning('Invalid Old format'); continue;
        end
    else
        D = S.Raw;
        V_avg = D.CVavg(:);
        I_rack = D.DCCurrent(:);
    end
    
    Np = 2;
    I_cell = I_rack / Np;
    full_cap = 64; % Ah
    thr_A = full_cap * 0.05 / Np;
    
    isChg = I_cell > thr_A;
    isDch = I_cell < -thr_A;
    
    chg_segs = find_segments(isChg);
    dch_segs = find_segments(isDch);
    
    % 충전 구간 선택 (가장 긴 것)
    if ~isempty(chg_segs)
        durs = chg_segs(:,2) - chg_segs(:,1);
        [~, max_idx] = max(durs);
        target_chg = chg_segs(max_idx, :);
        idx_s = target_chg(1); idx_e = target_chg(2);
        FieldData.(year_label).Chg.V = V_avg(idx_s:idx_e);
        FieldData.(year_label).Chg.I = I_cell(idx_s:idx_e);
        FieldData.(year_label).Chg.Q = cumtrapz(FieldData.(year_label).Chg.I) / 3600; 
    end
    
    % 방전 구간 선택 (가장 긴 것)
    if ~isempty(dch_segs)
        durs = dch_segs(:,2) - dch_segs(:,1);
        [~, max_idx] = max(durs);
        target_dch = dch_segs(max_idx, :);
        idx_s = target_dch(1); idx_e = target_dch(2);
        FieldData.(year_label).Dch.V = V_avg(idx_s:idx_e);
        FieldData.(year_label).Dch.I = I_cell(idx_s:idx_e);
        FieldData.(year_label).Dch.Q = cumtrapz(abs(FieldData.(year_label).Dch.I)) / 3600; 
    end
end

%% ========================================================================
% Section 3: 피처 추출 및 상태 추정
% ========================================================================
fprintf('\n=== Section 3: 피처 추출 및 상태 추정 ===\n');

years = fieldnames(FieldData);
Results = struct();

% 윈도우 설정 (실험 데이터와 동일하게)
win_chg_min = 3.7; win_chg_max = 3.95;
win_dch_min = 3.75; win_dch_max = 3.88;

for k = 1:length(years)
    yr = years{k};
    if ~isfield(FieldData.(yr), 'Chg') || ~isfield(FieldData.(yr), 'Dch')
        fprintf('  %s: 데이터 부족로 인해 스킵합니다.\n', yr);
        continue;
    end
    
    % 피처 추출
    data_chg.V_grid = FieldData.(yr).Chg.V;
    data_chg.Q = FieldData.(yr).Chg.Q;
    data_dch.V_grid = FieldData.(yr).Dch.V;
    data_dch.Q = FieldData.(yr).Dch.Q;
    
    [fea_chg_dQ, fea_chg_diff] = extract_features_half(data_chg, V_bounds_chg, 'charge', win_chg_min, win_chg_max);
    [fea_dch_dQ, fea_dch_diff] = extract_features_half(data_dch, V_bounds_dch, 'discharge', win_dch_min, win_dch_max);
    
    % 14개 피처 결합 (Lab 모델과 동일한 순서: [Chg_dQ(5), Dch_dQ(5), Chg_Pk(2), Dch_Pk(2)])
    field_features = [fea_chg_dQ, fea_dch_dQ, fea_chg_diff, fea_dch_diff];
    
    if any(isnan(field_features))
        warning('  %s: NaN detected in features.', yr);
        field_features(isnan(field_features)) = 0;
    end
    
    % 정규화 (실험 데이터 mu, sigma 사용)
    field_features_norm = (field_features - mu_X) ./ sigma_X;
    
    % 추정
    pred_rf = zeros(1, 3);
    for i = 1:3
        pred_rf(i) = predict(RF_Models{i}, field_features_norm);
    end
    
    Results.(yr).RF = pred_rf;
    
    fprintf('  %s NoPCA_RF Estimation:\n', yr);
    fprintf('    SOH: %.3f, LLI: %.3f, LAM: %.3f\n', pred_rf);
end

%% ========================================================================
% Section 4: 결과 시각화
% ========================================================================
fprintf('\n=== Section 4: 결과 시각화 ===\n');

yr_list = cellfun(@(x) str2double(strrep(x, 'Y', '')), fieldnames(Results));
num_yrs = length(yr_list);
res_fields = fieldnames(Results);

est_soh = zeros(num_yrs, 1);
est_lli = zeros(num_yrs, 1);
est_lam = zeros(num_yrs, 1);

for k = 1:num_yrs
    yr = res_fields{k};
    est_soh(k) = Results.(yr).RF(1);
    est_lli(k) = Results.(yr).RF(2);
    est_lam(k) = Results.(yr).RF(3);
end

fig = figure('Position', [100, 100, 1200, 400], 'Name', 'Field Estimation (NoPCA_RF)');

% SOH
subplot(1, 3, 1); hold on; grid on;
plot(yr_list, est_soh, '-^', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'r');
title('Field Data SOH (NoPCA\_RF)'); xlabel('Year'); ylabel('SOH (Ah)');
ylim([min(est_soh)*0.95, max(est_soh)*1.05]);

% LLI
subplot(1, 3, 2); hold on; grid on;
plot(yr_list, est_lli, '-^', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'b');
title('Field Data LLI (NoPCA\_RF)'); xlabel('Year'); ylabel('LLI (V shift)');

% LAM
subplot(1, 3, 3); hold on; grid on;
plot(yr_list, est_lam, '-^', 'LineWidth', 2, 'Color', [0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0]);
title('Field Data LAM (NoPCA\_RF)'); xlabel('Year'); ylabel('LAM (Ah/V)');

saveas(fig, fullfile(baseDir, 'Field_Estimation_NoPCA_Result.fig'));
fprintf('Figure saved as Field_Estimation_NoPCA_Result.fig\n');


%% Helper Functions
function segs = find_segments(mask)
    segs = [];
    n = length(mask);
    i = 1;
    while i <= n
        if mask(i)
            j = i;
            while j < n && mask(j+1), j = j + 1; end
            segs = [segs; i, j];
            i = j + 1;
        else
            i = i + 1;
        end
    end
end

function [features_dQ, features_diff] = extract_features_half(data_struct, V_ruler, mode, minV, maxV)
    V_grid = data_struct.V_grid;
    Q_val  = data_struct.Q;
    
    [V_grid, sort_idx] = sort(V_grid);
    Q_val = Q_val(sort_idx);
    
    [V_u, uid] = unique(V_grid); 
    Q_u = Q_val(uid);
    
    if numel(V_u) < 10
         features_dQ = nan(1, length(V_ruler)-1);
         features_diff = [NaN, NaN];
         return;
    end
    
    % 1) dQ Features
    try
        Q_ruler = interp1(V_u, Q_u, V_ruler, 'linear', 'extrap'); 
        features_dQ = abs(diff(Q_ruler));
    catch
        features_dQ = nan(1, length(V_ruler)-1);
    end
    
    % 2) dQ/dV Peak Features (Moving Average smoothing)
    try
        v_min_data = min(V_u);
        v_max_data = max(V_u);
        if v_max_data - v_min_data < 0.01
             features_diff = [0, 0]; return;
        end
        v_interp = v_min_data:0.001:v_max_data;
        v_interp = v_interp(:);
        q_interp = interp1(V_u, Q_u, v_interp, 'linear');
        
        dV = gradient(v_interp);
        dQ = gradient(q_interp);
        dV(dV==0) = NaN;
        dQdV_raw = dQ ./ dV;
        
        window_size = 21;
        dQdV_filt = movmean(dQdV_raw, window_size);
        
        if strcmp(mode, 'charge')
            mask = v_interp >= minV & v_interp <= win_chg_max;
        else
            mask = v_interp <= maxV & v_interp >= minV;
        end
        
        dQdV_win = dQdV_filt(mask);
        V_win = v_interp(mask);
        
        if isempty(dQdV_win)
            features_diff = [NaN, NaN];
        else
            pk_height = max(dQdV_win);
            pk_area = trapz(V_win, dQdV_win);
            features_diff = [pk_height, abs(pk_area)];
        end
    catch
        features_diff = [NaN, NaN];
    end
end
