function [X, Y, cv_group, feature_names, label_names, cellIDs, crateLabels] = RPT_ML_DataLoader(base_dir)
    % RPT_ML_DataLoader: 17개(혹은 18개) 피처와 3개 라벨(SOH, LLI, LAM)을 로딩하고,
    % K-Fold 교차검증 분할 인덱스를 생성하여 반환하는 공통 모듈
    
    %% 1. 데이터 로드
    feature_path = fullfile(base_dir, 'RPT_FeatureLabelExtractor', 'Feature_Matrix_ver02.mat');
    label_path = fullfile(base_dir, 'RPT_FeatureLabelExtractor', 'Label_Matrix_ver02.mat');
    
    fprintf('Loading Feature and Label Matrices from: %s\n', base_dir);
    feat_data = load(feature_path);
    label_data = load(label_path);
    
    FeatureTable = feat_data.FeatureTable_ver02;
    LabelTable = label_data.LabelTable_ver02;
    feature_names = feat_data.feature_names;
    
    % --- Feature Selection: exclude features not used for ML ---
    exclude_feats = {'T_chg_avg', 'T_dch_avg'};
    keep_mask = ~cellfun(@(f) ismember(f, exclude_feats), feature_names);
    feature_names = feature_names(keep_mask);
    fprintf('  -> Excluded features: %s\n', strjoin(exclude_feats, ', '));
    
    % `.mat` 파일 내부에 저장된 label_names를 무시하고 강제로 3개 라벨만 사용하도록 고정
    label_names = {'SOH', 'LLI', 'LAM'};
    
    % 1. Align rows using CellID + Cycle + CrateLabel string match
    [~, ia, ib] = intersect( ...
        strcat(FeatureTable.CellID, num2str(FeatureTable.Cycle), FeatureTable.CrateLabel), ...
        strcat(LabelTable.CellID,   num2str(LabelTable.Cycle),   LabelTable.CrateLabel));
    
    X_all = table2array(FeatureTable(ia, feature_names));
    Y_all = table2array(LabelTable(ib, label_names));
    cellIDs_all = FeatureTable.CellID(ia);
    crateLabels_all = FeatureTable.CrateLabel(ia);  % NEW: C-rate label per row
    
    % 2. Remove NaN rows
    nan_rows = any(isnan([X_all, Y_all]), 2);
    X = X_all(~nan_rows, :);
    Y = Y_all(~nan_rows, :);
    cellIDs = cellIDs_all(~nan_rows);
    crateLabels = crateLabels_all(~nan_rows);  % NEW: C-rate label (NaN rows removed)
    
    fprintf('  -> Extracted %d Valid Samples (after removing NaNs).\n', size(X,1));
    fprintf('  -> Features: %d | Labels: %d\n', length(feature_names), length(label_names));
    
    %% 4. 교차 검증 (Cross-Validation) 인덱스 생성
    cells_unique = unique(cellIDs);
    num_cells = length(cells_unique);
    fprintf('  -> Found %d unique cells for CV Split.\n', num_cells);
    
    % 전략 1: Group K-Fold (셀 단위 분할)
    % N개의 셀을 무작위로 K개의 그룹으로 묶어서 Fold를 형성합니다.
    K_group = min(5, num_cells);
    cv_group = zeros(size(X, 1), 1);
    
    rng(42); % 무작위 분할 결과 고정을 위한 임의의 시드
    shuffled_cells = cells_unique(randperm(num_cells));
    fold_assignment = mod(0:num_cells-1, K_group) + 1;
    
    for i = 1:num_cells
        c_name = shuffled_cells{i};
        f_idx = fold_assignment(i);
        cv_group(strcmp(cellIDs, c_name)) = f_idx;
    end
    
    fprintf('\n[CV Structures Prepared]\n');
    fprintf(' - Group K-Fold: %d Folds (grouped by cell IDs)\n', K_group);
end
