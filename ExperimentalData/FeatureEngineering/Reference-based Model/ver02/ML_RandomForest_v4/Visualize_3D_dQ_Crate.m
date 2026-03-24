%% Visualize_3D_dQ_Crate.m
% =====================================================================
% C-rate별 dQ의 3차원적 분포 시각화 (Segment vs C-rate vs dQ)
% - GPR 모델이 C_eff 내삽(Interpolation)을 어떻게 수행하는지 보여주기 위함
% =====================================================================
clear; clc; close all;

%% 1. 데이터 로드
base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                    'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02', 'ML_RandomForest_v4');
feat_path = fullfile(base_dir, 'RPT_FeatureLabelExtractor_v5_Charge', 'Feature_Matrix_v5_Charge.mat');

fprintf('Loading Feature Matrix...\n');
load(feat_path, 'FeatureTable_v5_Charge');
T = FeatureTable_v5_Charge;

% C-rate 목록
crates = [0.1, 0.5, 1.0, 2.0, 3.0];
N_seg = 11;
N_crates = length(crates);

% 평균 dQ 저장할 행렬 (C-rate x Segments)
mean_dQ = zeros(N_crates, N_seg);

for i = 1:N_crates
    cr = crates(i);
    % 해당 C-rate의 데이터만 추출
    idx = T.CrateNum == cr;
    X_seg = table2array(T(idx, 5:15)); % dQ_chg_S1 ~ dQ_chg_S11
    
    % NaN 무시하고 세그먼트별 평균 계산
    mean_dQ(i, :) = mean(X_seg, 1, 'omitnan');
end

%% 2. 3D 시각화 (Surface)
[X, Y] = meshgrid(1:N_seg, crates);

fig = figure('Name', '3D dQ Surface over C-rates', 'Position', [100, 100, 800, 600]);

% Surface Plot
surf(X, Y, mean_dQ, 'FaceAlpha', 0.8, 'EdgeColor', 'k');
hold on; grid on;

% 각각의 데이터 포인트 (Scatter) 추가하여 명확하게 보이기
scatter3(X(:), Y(:), mean_dQ(:), 50, 'r', 'filled', 'MarkerEdgeColor', 'w');

colormap(gca, turbo);
colorbar;
caxis([0, max(mean_dQ(:))]);

% 축 설정
xlabel('Segment (1 to 11)', 'FontWeight', 'bold');
ylabel('C-rate (C_{eff} from Lab)', 'FontWeight', 'bold');
zlabel('Average \DeltaQ (Ah)', 'FontWeight', 'bold');
title('3D Map of \DeltaQ across Voltage Segments & C-rates', 'FontSize', 14, 'FontWeight', 'bold');

xticks(1:N_seg);
yticks(crates);
set(gca, 'YDir', 'reverse'); % C-rate가 앞으로 나오게
view(-35, 45); % 3D 보기 각도 설정

%% [선택] 현장(Field)의 가상 C_eff 궤도선(Interpolation Line) 추가
% 예: 현장의 C_eff가 1.5C 라고 가정할 때, 이 모델이 상상하는 1.5C 라인을 그림
c_eff_field = 1.5;
% 보간 (Spline)
vq_15 = interp1(crates, mean_dQ, c_eff_field, 'spline');
plot3(1:N_seg, repmat(c_eff_field, 1, N_seg), vq_15, 'g-', 'LineWidth', 4, 'DisplayName', 'Interpolated Field CP (e.g. 1.5C)');

legend('Lab Data Surface', 'Lab Anchor Points', 'Interpolated Field CP (1.5C)', 'Location', 'northeast');

% 저장
saveDir = base_dir;
saveas(fig, fullfile(saveDir, '3D_dQ_Crate_Map_v6.fig'));
saveas(fig, fullfile(saveDir, '3D_dQ_Crate_Map_v6.png'));
fprintf('>> Saved Figures: 3D_dQ_Crate_Map_v6.fig / .png\n');
