%% ML_Feature_Selector.m
% =====================================================================
% 필드 제약사항과 C_eff를 고려하여 1~3개 세그먼트 조합의 RF 성능 평가
% Track A: 넓은 시야 (Seg 5,8,9,10 활용 가능)
% Track B: 좁은 시야 (Seg 7,8,9 필수 활용)
% 모든 조합에는 C_eff(여기서는 data_CrateNum)가 포함됨
% =====================================================================
clear; clc; close all;
rng(42); % 재현성

% 데이터 로드
d_pcc = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\PCC_v5b_data.mat');
X_chg     = d_pcc.dQ_chg_all;  % [N x 11]
y_soh     = d_pcc.data_SOH;    % [N x 1]
c_eff     = d_pcc.data_CrateNum; % [N x 1]
cell_ids  = d_pcc.data_CellID;   % [N x 1] (cell array)

% NaN 제거 (SOH 또는 C_eff가 NaN인 샘플 제외, 피처 NaN 처리는 조합별 수행)
valid_all = ~isnan(y_soh) & ~isnan(c_eff);
X_chg = X_chg(valid_all,:);
y_soh = y_soh(valid_all);
c_eff = c_eff(valid_all);
cell_ids = cell_ids(valid_all);

% Group K-Fold 설정 (K=5)
u_cells = unique(cell_ids);
K = 5;
cv_indices = zeros(size(y_soh));
groups = ceil(linspace(0.01, K, length(u_cells)));
for k=1:K
    test_cells = u_cells(groups==k);
    for ti = 1:length(test_cells)
        cv_indices(strcmp(cell_ids, test_cells{ti})) = k;
    end
end

% 테스트할 세그먼트 풀
pool_A = [10, 9, 5, 8]; % Track A 주요 강력 세그먼트
pool_B = [7, 8, 9];     % Track B 필수 공통 세그먼트

%% ==== COMBINATION GENERATOR ====
combos = {};
combos{end+1} = struct('name','Track A [1] Seg10', 'segs',[10]);
combos{end+1} = struct('name','Track A [2] Seg10+9', 'segs',[10,9]);
combos{end+1} = struct('name','Track A [2] Seg10+5', 'segs',[10,5]);
combos{end+1} = struct('name','Track A [3] Seg10+9+5', 'segs',[10,9,5]);
combos{end+1} = struct('name','Track A [3] Seg10+9+8', 'segs',[10,9,8]);

combos{end+1} = struct('name','Track B [1] Seg9', 'segs',[9]);
combos{end+1} = struct('name','Track B [2] Seg9+8', 'segs',[9,8]);
combos{end+1} = struct('name','Track B [2] Seg9+7', 'segs',[9,7]);
combos{end+1} = struct('name','Track B [3] Seg9+8+7', 'segs',[9,8,7]);

% Baseline
combos{end+1} = struct('name','Baseline (All 11 Segs)', 'segs', 1:11);

%% ==== EVALUATION LOOP ====
fprintf('=== ML Feature Selection: Random Forest (Group 5-Fold) ===\n');
fprintf('모든 모델에 [C_eff] 피처가 기본으로 포함됩니다.\n');
fprintf('%-25s | 피처 수 | R2     | RMSE \n', 'Combination Name');
fprintf(repmat('-', 1, 55)); fprintf('\n');

for c = 1:length(combos)
    cb = combos{c};
    feat_idx = cb.segs;
    
    % 피처 행렬 구성: [선택된 세그먼트 ΔQ, C_eff]
    X_comb = [X_chg(:, feat_idx), c_eff];
    
    % 해당 조합에서 NaN을 가지는 행 제거
    valid_comb = ~any(isnan(X_comb), 2);
    X_comb = X_comb(valid_comb,:);
    y_comb = y_soh(valid_comb);
    cv_idx = cv_indices(valid_comb);
    
    if isempty(X_comb)
        fprintf('%-25s | %d      | N/A    | N/A\n', cb.name, length(feat_idx)+1);
        continue;
    end
    
    y_pred = nan(size(y_comb));
    for k = 1:K
        train_mask = (cv_idx ~= k);
        test_mask  = (cv_idx == k);
        
        X_train = X_comb(train_mask, :); y_train = y_comb(train_mask);
        X_test  = X_comb(test_mask,  :); y_test  = y_comb(test_mask);
        
        if isempty(X_train) || isempty(X_test), continue; end
        
        % Train RF
        model = TreeBagger(50, X_train, y_train, 'Method','regression',...
                           'MinLeafSize', 5);
        y_pred(test_mask) = predict(model, X_test);
    end
    
    % Calculate metrics
    rmse = sqrt(mean((y_comb - y_pred).^2, 'omitnan'));
    y_mean = mean(y_comb, 'omitnan');
    ss_tot = sum((y_comb - y_mean).^2, 'omitnan');
    ss_res = sum((y_comb - y_pred).^2, 'omitnan');
    r2 = 1 - (ss_res / ss_tot);
    
    fprintf('%-25s | %.0f      | %6.3f | %5.3f\n', cb.name, length(feat_idx)+1, r2, rmse);
end
