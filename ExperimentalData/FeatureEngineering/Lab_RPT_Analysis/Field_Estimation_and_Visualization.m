%% Field_Estimation_and_Visualization.m
% 1. 실험 데이터(Lab Data)로 MLR/RF 모델 학습
% 2. 필드 데이터(Field Data, 2021~2025) 로드 및 피처 추출 (Master Ruler 적용)
% 3. 학습된 모델을 사용하여 필드 데이터의 SOH, LLI, LAM 추정 및 시각화

clear; clc; close all;
warning off;

%% ========================================================================
% Section 1: 모델 학습 (실험 데이터 사용)
% ========================================================================
fprintf('=== Section 1: Model Training with Experimental Data ===\n');

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
X_lab = FeatureTable.X_Features; % 정규화 전 원본 피처
Y_lab = FeatureTable.Y_Labels;   % [SOH, LLI, LAM]
label_names = {'SOH', 'LLI', 'LAM'};

% 1.3 정규화 파라미터 계산 (실험 데이터 기준)
mu_X = mean(X_lab, 1, 'omitnan');
sigma_X = std(X_lab, 0, 1, 'omitnan');
sigma_X(sigma_X == 0) = 1;

X_lab_norm = (X_lab - mu_X) ./ sigma_X;

% 1.4 PCA 학습 (95% 설명력)
[coeff_pca, score_lab, ~, ~, explained] = pca(X_lab_norm);
cum_var = cumsum(explained);
n_pc = find(cum_var >= 95, 1, 'first');
fprintf('  PCA: Selected %d PCs (%.2f%% Variance)\n', n_pc, cum_var(n_pc));

X_lab_pc = score_lab(:, 1:n_pc);

% 1.5 MLR 및 RF 모델 학습
MLR_Models = cell(1, 3);
RF_Models  = cell(1, 3);
rf_trees = 100;
rf_minleaf = 3;

for i = 1:3
    fprintf('  Training models for %s...\n', label_names{i});
    % MLR
    MLR_Models{i} = fitlm(X_lab_pc, Y_lab(:, i));
    
    % RF
    t = templateTree('MinLeafSize', rf_minleaf);
    RF_Models{i} = fitrensemble(X_lab_pc, Y_lab(:, i), ...
        'Method', 'Bag', 'NumLearningCycles', rf_trees, 'Learners', t);
end

fprintf('=== Model Training Complete ===\n');


%% ========================================================================
% Section 2: 필드 데이터 로드 및 전처리
% ========================================================================
fprintf('\n=== Section 2: Loading Field Data (2021~2025) ===\n');

% 경로 설정 (사용자 환경에 맞게 수정 필요 시 수정)
rawDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
dataFiles = {
    fullfile(rawDir, 'Old', '2021', '202106', 'Raw_20210603.mat'), 'Y2021';
    fullfile(rawDir, 'New', '2023', '202310', 'Raw_20231016.mat'), 'Y2023';
    fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024';
    fullfile(rawDir, 'New', '2025', '202507', 'Raw_20250711.mat'), 'Y2025'
};

% Feature Extraction 설정
% Master Ruler는 실험 데이터 중 하나(예: Ch09)를 기준으로 사용
% 모든 채널이 동일한 Global Ruler를 공유하므로 Ch09의 Ruler 사용
ref_ch = 'Ch09'; 
if isfield(MasterRulers, ref_ch)
    V_bounds_chg = MasterRulers.(ref_ch).V_bounds_chg;
    V_bounds_dch = MasterRulers.(ref_ch).V_bounds_dch;
else
    % Fallback: 첫 번째 필드 사용
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
    
    fprintf('  Processing %s (%s)...\n', year_label, fpath);
    S = load(fpath);
    
    % 데이터 포맷 파싱 (Reference/Field_Visualize_VIT_byYear.m 참조)
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
    
    % 셀 단위 전류로 변환 (Np=2 가정)
    Np = 2;
    I_cell = I_rack / Np;
    
    % 세그먼트 분리 (0.05C 이상 전류 기준)
    full_cap = 64; % Ah
    thr_A = full_cap * 0.05 / Np;
    
    isChg = I_cell > thr_A;
    isDch = I_cell < -thr_A;
    
    chg_segs = find_segments(isChg);
    dch_segs = find_segments(isDch);
    
    % 유효 세그먼트 필터링 (Capacity 기준, 약 0.5C 충전/방전 구간 찾기)
    % 실험 데이터가 0.1C, 0.5C, 1C 등으로 학습됨.
    % 필드 정기평가는 보통 1C 방전 등을 포함함.
    % 여기서는 가장 긴 구간(용량이 가장 큰 구간)을 메인 평가 구간으로 가정.
    
    % 충전 구간 선택 (가장 긴 것)
    if ~isempty(chg_segs)
        durs = chg_segs(:,2) - chg_segs(:,1);
        [~, max_idx] = max(durs);
        target_chg = chg_segs(max_idx, :);
        
        idx_s = target_chg(1); idx_e = target_chg(2);
        FieldData.(year_label).Chg.V = V_avg(idx_s:idx_e);
        FieldData.(year_label).Chg.I = I_cell(idx_s:idx_e);
        % Ah 계산 (dt=1 assuming 1Hz)
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
        % Ah 계산 (방전은 -I 이므로 절대값 혹은 -적분)
        FieldData.(year_label).Dch.Q = cumtrapz(abs(FieldData.(year_label).Dch.I)) / 3600; 
    end
end

%% ========================================================================
% Section 2.1: 필드 데이터 세그먼트 시각화 (데이터 검증용)
% ========================================================================
fprintf('  Visualizing Selected Segments...\n');
fig_seg = figure('Position', [50, 50, 1200, 600], 'Name', 'Field Data Selected Segments');
years = fieldnames(FieldData);
colors = lines(length(years));

subplot(2,1,1); hold on; grid on; title('Selected Charge Segments (V)');
xlabel('Capacity (Ah)'); ylabel('Voltage (V)');
for k = 1:length(years)
    yr = years{k};
    if isfield(FieldData.(yr), 'Chg')
        plot(FieldData.(yr).Chg.Q, FieldData.(yr).Chg.V, 'LineWidth', 1.5, ...
            'Color', colors(k,:), 'DisplayName', yr);
    end
end
legend('Location', 'best');

subplot(2,1,2); hold on; grid on; title('Selected Discharge Segments (V)');
xlabel('Capacity (Ah)'); ylabel('Voltage (V)');
for k = 1:length(years)
    yr = years{k};
    if isfield(FieldData.(yr), 'Dch')
        plot(FieldData.(yr).Dch.Q, FieldData.(yr).Dch.V, 'LineWidth', 1.5, ...
            'Color', colors(k,:), 'DisplayName', yr);
    end
end
saveas(fig_seg, fullfile(baseDir, 'Field_Selected_Segments.fig'));

%% ========================================================================
% Section 3: 필드 데이터 피처 추출 & 추정
% ========================================================================
fprintf('\n=== Section 3: Feature Extraction & State Estimation ===\n');

years = fieldnames(FieldData);
Results = struct();

% 윈도우 설정 (실험 데이터와 동일하게)
win_chg_min = 3.7; win_chg_max = 3.95;
win_dch_min = 3.75; win_dch_max = 3.88;

for k = 1:length(years)
    yr = years{k};
    if ~isfield(FieldData.(yr), 'Chg') || ~isfield(FieldData.(yr), 'Dch')
        fprintf('  %s: Insufficient charge or discharge data\n', yr);
        continue;
    end
    
    % 데이터 준비
    data_chg.V_grid = FieldData.(yr).Chg.V;
    data_chg.Q = FieldData.(yr).Chg.Q;
    
    data_dch.V_grid = FieldData.(yr).Dch.V;
    data_dch.Q = FieldData.(yr).Dch.Q;
    
    % 피처 추출 (Helper Function 사용)
    [fea_chg_dQ, fea_chg_diff] = extract_features_half(data_chg, V_bounds_chg, 'charge', win_chg_min, win_chg_max);
    [fea_dch_dQ, fea_dch_diff] = extract_features_half(data_dch, V_bounds_dch, 'discharge', win_dch_min, win_dch_max);
    
    % 피처 결합: [Chg_dQ(5), Dch_dQ(5), Chg_PkH, Chg_PkA, Dch_PkH, Dch_PkA] -> 총 14개
    field_features = [fea_chg_dQ, fea_dch_dQ, fea_chg_diff, fea_dch_diff];
    
    if any(isnan(field_features))
        warning('  %s: NaN detected in features. Estimation may be inaccurate.', yr);
        field_features(isnan(field_features)) = 0; % 임시 처리
    end
    
    % 정규화 (실험 데이터 기준)
    field_features_norm = (field_features - mu_X) ./ sigma_X;
    
    % PCA 변환
    field_pc = field_features_norm * coeff_pca(:, 1:n_pc);
    
    % 추정 (MLR & RF)
    pred_mlr = zeros(1, 3);
    pred_rf  = zeros(1, 3);
    
    for i = 1:3
        pred_mlr(i) = predict(MLR_Models{i}, field_pc);
        pred_rf(i)  = predict(RF_Models{i}, field_pc);
    end
    
    Results.(yr).MLR = pred_mlr;
    Results.(yr).RF  = pred_rf;
    
    fprintf('  %s Estimation:\n', yr);
    fprintf('    MLR -> SOH: %.3f, LLI: %.3f, LAM: %.3f\n', pred_mlr);
    fprintf('    RF  -> SOH: %.3f, LLI: %.3f, LAM: %.3f\n', pred_rf);
end

%% ========================================================================
% Section 3.1: 피처 추출 검증 시각화 (dQ/dV + Master Ruler)
% ========================================================================
fprintf('  Visualizing Feature Extraction (dQ/dV)...\n');
fig_feat = figure('Position', [100, 100, 1200, 500], 'Name', 'Feature Extraction Verification');

% Charge dQ/dV
subplot(1,2,1); hold on; grid on; title('Charge dQ/dV & Master Ruler');
xlabel('Voltage (V)'); ylabel('dQ/dV (Ah/V)');
for k = 1:length(years)
    yr = years{k};
    if isfield(FieldData.(yr), 'Chg')
        % Recalculate dQ/dV for plotting
        V = FieldData.(yr).Chg.V; Q = FieldData.(yr).Chg.Q;
        [V_u, uid] = unique(V); Q_u = Q(uid);
        dQdV = gradient(Q_u) ./ gradient(V_u);
        dQdV = movmean(dQdV, 21);
        plot(V_u, dQdV, 'LineWidth', 1.2, 'Color', colors(k,:), 'DisplayName', yr);
    end
end
% Add Master Ruler Lines
ylim_c = ylim;
for i = 1:length(V_bounds_chg)
    xline(V_bounds_chg(i), '--k', 'HandleVisibility', 'off');
end
legend('Location', 'best');

% Discharge dQ/dV
subplot(1,2,2); hold on; grid on; title('Discharge dQ/dV & Master Ruler');
xlabel('Voltage (V)'); ylabel('dQ/dV (Ah/V)');
for k = 1:length(years)
    yr = years{k};
    if isfield(FieldData.(yr), 'Dch')
        V = FieldData.(yr).Dch.V; Q = FieldData.(yr).Dch.Q;
        % Sort descending for discharge plotting stability
        [V_u, uid] = unique(V); Q_u = Q(uid); 
        dQdV = gradient(Q_u) ./ gradient(V_u);
        dQdV = movmean(abs(dQdV), 21); % Absolute for discharge
        plot(V_u, dQdV, 'LineWidth', 1.2, 'Color', colors(k,:), 'DisplayName', yr);
    end
end
% Add Master Ruler Lines
ylim_d = ylim;
for i = 1:length(V_bounds_dch)
    xline(V_bounds_dch(i), '--k', 'HandleVisibility', 'off');
end

saveas(fig_feat, fullfile(baseDir, 'Field_Feature_Extraction_Verify.fig'));

%% ========================================================================
% Section 4: 결과 시각화
% ========================================================================
fprintf('\n=== Section 4: Visualization ===\n');

yr_list = cellfun(@(x) str2double(strrep(x, 'Y', '')), years);
num_yrs = length(yr_list);

est_soh_mlr = zeros(num_yrs, 1); est_soh_rf = zeros(num_yrs, 1);
est_lli_mlr = zeros(num_yrs, 1); est_lli_rf = zeros(num_yrs, 1);
est_lam_mlr = zeros(num_yrs, 1); est_lam_rf = zeros(num_yrs, 1);

for k = 1:num_yrs
    yr = years{k};
    est_soh_mlr(k) = Results.(yr).MLR(1);
    est_soh_rf(k)  = Results.(yr).RF(1);
    
    est_lli_mlr(k) = Results.(yr).MLR(2);
    est_lli_rf(k)  = Results.(yr).RF(2);
    
    est_lam_mlr(k) = Results.(yr).MLR(3);
    est_lam_rf(k)  = Results.(yr).RF(3);
end

fig = figure('Position', [100, 100, 1200, 400], 'Name', 'Field Estimation Result');

% SOH
subplot(1, 3, 1); hold on; grid on;
plot(yr_list, est_soh_mlr, '-o', 'LineWidth', 2, 'DisplayName', 'MLR');
plot(yr_list, est_soh_rf,  '-^', 'LineWidth', 2, 'DisplayName', 'RF');
title('Field Data SOH Estimation'); xlabel('Year'); ylabel('SOH (Ah)');
legend;

% LLI
subplot(1, 3, 2); hold on; grid on;
plot(yr_list, est_lli_mlr, '-o', 'LineWidth', 2, 'DisplayName', 'MLR');
plot(yr_list, est_lli_rf,  '-^', 'LineWidth', 2, 'DisplayName', 'RF');
title('Field Data LLI Estimation'); xlabel('Year'); ylabel('LLI (V shift)');

% LAM
subplot(1, 3, 3); hold on; grid on;
plot(yr_list, est_lam_mlr, '-o', 'LineWidth', 2, 'DisplayName', 'MLR');
plot(yr_list, est_lam_rf,  '-^', 'LineWidth', 2, 'DisplayName', 'RF');
title('Field Data LAM Estimation'); xlabel('Year'); ylabel('LAM (Ah/V)');

saveas(fig, fullfile(baseDir, 'Field_Estimation_Result.fig'));
fprintf('Figure saved.\n');

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
    % RPT_Feature_Extraction_Advanced.m의 로직 복사 (단순화)
    V_grid = data_struct.V_grid;
    Q_val  = data_struct.Q;
    
    % 전압 오름차순 정렬 (필수)
    [V_grid, sort_idx] = sort(V_grid);
    Q_val = Q_val(sort_idx);
    
    % 중복 제거
    [V_u, uid] = unique(V_grid); 
    Q_u = Q_val(uid);
    
    if numel(V_u) < 10
         features_dQ = nan(1, length(V_ruler)-1);
         features_diff = [NaN, NaN];
         return;
    end
    
    % 1) dQ Features
    % Ruler 전압에 해당하는 Q 보간
    try
        Q_ruler = interp1(V_u, Q_u, V_ruler, 'linear', 'extrap'); % Extrap 허용 (필드 데이터 범위 문제 대비)
        features_dQ = abs(diff(Q_ruler));
    catch
        features_dQ = nan(1, length(V_ruler)-1);
    end
    
    % 2) dQ/dV Peak Features
    % 2) dQ/dV Peak Features
    % Preprocessing: Resample to 0.001V grid (mimic RPT_VQ_grid quality)
    try
        v_min = min(V_u);
        v_max = max(V_u);
        if v_max - v_min < 0.01
             features_diff = [0, 0]; return;
        end
        v_interp = v_min:0.001:v_max;
        v_interp = v_interp(:);
        q_interp = interp1(V_u, Q_u, v_interp, 'linear');
        
        % Update V_u, Q_u for gradient calculation
        V_u = v_interp;
        Q_u = q_interp;
    catch
        features_diff = [NaN, NaN]; return;
    end

    dV = gradient(V_u);
    dQ = gradient(Q_u);
    dV(dV==0) = NaN;
    dQdV_raw = dQ ./ dV;
    
    % Smoothing Method: Moving Average (Window=21) - EXACT MATCH to Lab Script
    window_size = 21;
    dQdV_filt = movmean(dQdV_raw, window_size);
    
    if strcmp(mode, 'charge')
        mask = V_u >= minV & V_u <= maxV;
    else
        mask = V_u <= maxV & V_u >= minV;
    end
    
    dQdV_win = dQdV_filt(mask);
    V_win = V_u(mask);
    
    if isempty(dQdV_win)
        features_diff = [NaN, NaN];
    else
        pk_height = max(dQdV_win);
        pk_area = trapz(V_win, dQdV_win);
        features_diff = [pk_height, abs(pk_area)];
    end
end
