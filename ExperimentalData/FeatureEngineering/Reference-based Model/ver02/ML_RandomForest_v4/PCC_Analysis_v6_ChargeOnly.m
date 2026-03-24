%% PCC_Analysis_v6_ChargeOnly.m
% =====================================================================
% NCM 파편화 충전 용량 (Fragmented Charge Capacity) 분석
% - v5_Charge feature matrix를 로드
% - 각 C-rate 및 통합 (Pooled) 환경에서 11개 세그먼트의 dQ와 SOH간 PCC 계산
% - 등용량 분할의 철학에 따라 상관계수가 고르게 분포되는지 확인
% =====================================================================
clear; clc; close all;

%% 1. 데이터 로드
base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                    'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
feat_path = fullfile(base_dir, 'ML_RandomForest_v4', 'RPT_FeatureLabelExtractor_v5_Charge', 'Feature_Matrix_v5_Charge.mat');
label_path = fullfile(base_dir, 'RPT_FeatureLabelExtractor', 'Label_Matrix_ver02.mat');
saveDir = fullfile(base_dir, 'ML_RandomForest_v4');

fprintf('Loading v5 Charge Features...\n');
feat_data = load(feat_path);
label_data = load(label_path);

FeatureTable = feat_data.FeatureTable_v5_Charge;
LabelTable   = label_data.LabelTable_ver02;

% 행 매칭 (CellID + Cycle + CrateLabel)
[~, ia, ib] = intersect( ...
    strcat(FeatureTable.CellID, num2str(FeatureTable.Cycle), FeatureTable.CrateLabel), ...
    strcat(LabelTable.CellID,   num2str(LabelTable.Cycle),   LabelTable.CrateLabel));

% 필요한 데이터 추출
% FeatureTable 구성: CellID(1), Cycle(2), CrateLabel(3), CrateNum(4), dQ_chg_S1~S11(5~15), Peak_H_chg(16), ...
feature_cols = 5:15; 
X_all = table2array(FeatureTable(ia, feature_cols)); % dQ_chg_S1 ~ dQ_chg_S11
Y_soh = LabelTable.SOH(ib);
C_num = FeatureTable.CrateNum(ia);

% NaN 필터링
valid = ~isnan(Y_soh) & ~any(isnan(X_all), 2);
X_all = X_all(valid, :);
Y_soh = Y_soh(valid);
C_num = C_num(valid);

N_seg = 11;
crates = [0.1, 0.5, 1, 2, 3];
clabels = {'0.1C', '0.5C', '1C', '2C', '3C'};

fprintf('Total valid samples: %d\n', length(Y_soh));

%% 2. PCC 계산
% 2.1 통합 (Pooled)
pcc_pooled = compute_pcc(X_all, Y_soh);

% 2.2 C-rate별
pcc_cr = nan(length(crates), N_seg);
for r = 1:length(crates)
    idx_r = C_num == crates(r);
    pcc_cr(r, :) = compute_pcc(X_all(idx_r, :), Y_soh(idx_r));
end

%% 3. 결과 출력
fprintf('\n=== PCC: ΔQ_chg vs SOH (전 C-rate 풀링) ===\n');
for s = 1:N_seg
    fprintf('  Seg %2d: %+.3f\n', s, pcc_pooled(s));
end

%% 4. 시각화
fig = figure('Name', 'PCC Analysis v6 (Charge-Only NCM)', 'Position', [100, 100, 1200, 400]);

% 4.1 Pooled
subplot(1,2,1);
imagesc(pcc_pooled(:)'); colormap(gca, redblue(256)); clim([-1 1]); colorbar;
xticks(1:N_seg); yticks(1); yticklabels({'SOH'});
xticklabels(arrayfun(@(i) sprintf('S%d',i), 1:N_seg, 'UniformOutput', false));
title('Pooled Data (All C-rates) PCC', 'FontWeight', 'bold');
for s = 1:N_seg
    col = 'w'; if abs(pcc_pooled(s)) < 0.5, col = 'k'; end
    text(s, 1, sprintf('%.2f', pcc_pooled(s)), 'Color', col, 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

% 4.2 C-rate별
subplot(1,2,2);
imagesc(pcc_cr); colormap(gca, redblue(256)); clim([-1 1]); colorbar;
xticks(1:N_seg); yticks(1:5);
xticklabels(arrayfun(@(i) sprintf('S%d',i), 1:N_seg, 'UniformOutput', false));
yticklabels(clabels); xlabel('Segment'); ylabel('C-rate');
title('PCC by C-rate (0.1C ~ 3.0C)', 'FontWeight', 'bold');
for r = 1:5
    for s = 1:N_seg
        if ~isnan(pcc_cr(r,s))
            col = 'w'; if abs(pcc_cr(r,s)) < 0.5, col = 'k'; end
            text(s, r, sprintf('%.2f', pcc_cr(r,s)), 'Color', col, 'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end
end

sgtitle('PCC Analysis: Equi-Capacity Charge Segments (v6 NCM)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig, fullfile(saveDir, 'PCC_Heatmap_v6_ChargeOnly.fig'));
saveas(fig, fullfile(saveDir, 'PCC_Heatmap_v6_ChargeOnly.png'));

fprintf('\n>> Saved Figures: PCC_Heatmap_v6_ChargeOnly.fig / .png\n');


%% Helpers
function pcc = compute_pcc(X, y)
    pcc = nan(1, size(X,2));
    if length(y) < 10, return; end
    for s = 1:size(X,2)
        r = corrcoef(X(:,s), y);
        pcc(s) = r(1,2);
    end
end

function c = redblue(n)
    if nargin < 1, n = 256; end
    half = floor(n/2);
    r = [linspace(0,1,half)'; ones(n-half,1)];
    b = [ones(half,1); linspace(1,0,n-half)'];
    g = [linspace(0,1,half)'; linspace(1,0,n-half)'];
    c = [r g b];
end
