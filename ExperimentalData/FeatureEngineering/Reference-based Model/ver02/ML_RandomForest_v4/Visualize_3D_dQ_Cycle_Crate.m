%% Visualize_3D_dQ_Cycle_Crate.m
% =====================================================================
% C-rate 별 & Cycle 별 dQ의 3차원적 분포 시각화 
% - 각 C-rate(0.1C ~ 3.0C)에 대해 X=Segment, Y=Cycle, Z=dQ Surface 출력
% - 배터리가 열화(Cycle 증가)됨에 따라 각 전압 세그먼트의 용량이 
%   어떻게 변하는지 입체적으로 확인
% =====================================================================
clear; clc; close all;

%% 1. 데이터 로드
base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
                    'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02', 'ML_RandomForest_v4');
feat_path = fullfile(base_dir, 'RPT_FeatureLabelExtractor_v5_Charge', 'Feature_Matrix_v5_Charge.mat');

fprintf('Loading Feature Matrix...\n');
load(feat_path, 'FeatureTable_v5_Charge');
T = FeatureTable_v5_Charge;

% C-rate 및 Cycle 목록 확보
crates = [0.1, 0.5, 1.0, 2.0, 3.0];
clabels = {'0.1C', '0.5C', '1.0C', '2.0C', '3.0C'};
unique_cycles = sort(unique(T.Cycle));
N_seg = 11;
N_crates = length(crates);
N_cycles = length(unique_cycles);

%% 2. 3D 시각화 (Subplots for each C-rate)
fig = figure('Name', '3D dQ Surface: Cycle & C-rate', 'Position', [50, 50, 1600, 600]);

% X, Y 그리드 (X = Segment, Y = Cycle)
[X, Y] = meshgrid(1:N_seg, unique_cycles);

for i = 1:N_crates
    cr = crates(i);
    
    % 현재 C-rate에 해당하는 Z(dQ) 매트릭스 생성 (Cycles x Segments)
    Z_dQ = nan(N_cycles, N_seg);
    
    for j = 1:N_cycles
        cyc = unique_cycles(j);
        
        idx = (T.CrateNum == cr) & (T.Cycle == cyc);
        if sum(idx) > 0
            % 여러 셀(Cell09, Cell10 등)이 있을 수 있으므로 평균
            X_seg = table2array(T(idx, 5:15));
            Z_dQ(j, :) = mean(X_seg, 1, 'omitnan');
        end
    end
    
    % Subplot 그리기
    subplot(1, 5, i);
    surf(X, Y, Z_dQ, 'FaceAlpha', 0.85, 'EdgeColor', 'k');
    hold on; grid on;
    
    % Scatter 혼합
    % scatter3(X(:), Y(:), Z_dQ(:), 15, 'r', 'filled', 'MarkerEdgeColor', 'w');
    
    colormap(gca, turbo);
    colorbar('Location', 'southoutside');
    caxis([0, 10]); % 범위를 통일시켜 비교하기 쉽게 (최대 약 9Ah 언저리로 예상)
    
    % 축 설정
    xlabel('Segment (1 to 11)', 'FontWeight', 'bold');
    ylabel('Cycle', 'FontWeight', 'bold');
    zlabel('\DeltaQ (Ah)', 'FontWeight', 'bold');
    title(sprintf('%s (C_{eff} = %.1f)', clabels{i}, cr), 'FontSize', 12, 'FontWeight', 'bold');
    
    xticks(1:2:N_seg);
    yticks(unique_cycles(1:2:end));
    set(gca, 'YDir', 'reverse'); % Cycle 0 (프레쉬) 가 앞으로 오게
    view(-45, 35);
end

sgtitle('3D Map of \DeltaQ Evolution by Cycle for each C-rate (NCM 11 Segments)', 'FontSize', 16, 'FontWeight', 'bold');

% 결과 저장
saveDir = base_dir;
saveas(fig, fullfile(saveDir, '3D_dQ_Cycle_Crate_Map_v6.fig'));
saveas(fig, fullfile(saveDir, '3D_dQ_Cycle_Crate_Map_v6.png'));
fprintf('\n>> Saved Figures: 3D_dQ_Cycle_Crate_Map_v6.fig / .png\n');
