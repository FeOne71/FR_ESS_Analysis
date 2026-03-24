%% Plot_3D_PhysicalMapping.m
% =====================================================================
% Phase 4: 랩 데이터(CC) Physical Mapping (고해상도 3D Surface 시각화)
% - X축: Voltage Segment (1 ~ 11)
% - Y축: C-rate (0.1 ~ 3.0)
% - Z축: dQ (단위 세그먼트 당 용량)
% 사용자 פי드백 반영: 거친 면이나 점 대신 해상도를 높여 부드러운 3D 지형(Topography)으로 렌더링
% 딥러닝 모델이 인식할 비선형 공간(Manifold)을 직관적으로 보여줌
% =====================================================================
clear; clc; close all;

%% 1. 데이터 로드 및 초기화
d_pcc = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\PCC_v5b_data.mat');
X_chg = d_pcc.dQ_chg_all;  
X_dch = d_pcc.dQ_dch_all;  
c_eff = d_pcc.data_CrateNum; 
y_soh = d_pcc.data_SOH;      

N_seg = 11;
% 분석 대상 C-rate (로그 스케일에 가까움)
target_crates = [0.1, 0.5, 1.0, 2.0, 3.0];
N_c = length(target_crates);

% 원본 Grid (11 x 5)
[Seg_Grid, C_Grid] = meshgrid(1:N_seg, target_crates);

Z_Chg = nan(N_c, N_seg); 
Z_Dch = nan(N_c, N_seg); 

%% 2. 데이터 포인트 수집 (C-rate별, 세그먼트별 평균 dQ 추출)
for i = 1:N_c
    cr = target_crates(i);
    idx_c = abs(c_eff - cr) < 0.05;
    
    if sum(idx_c) == 0, continue; end
    
    for s = 1:N_seg
        % Charge
        vec_q_c = X_chg(idx_c, s);
        valid_c = ~isnan(vec_q_c);
        if any(valid_c), Z_Chg(i, s) = mean(vec_q_c(valid_c)); end
        
        % Discharge
        vec_q_d = X_dch(idx_c, s);
        valid_d = ~isnan(vec_q_d);
        if any(valid_d), Z_Dch(i, s) = mean(vec_q_d(valid_d)); end
    end
end

% 빈 공간은 주로 고전압 고-C-rate 구간이므로 0으로 보간하여 지형(Surface)이 끊어지지 않게 유지
Z_Chg(isnan(Z_Chg)) = 0;
Z_Dch(isnan(Z_Dch)) = 0;

%% 3. 고해상도 지형 생성을 위한 보간 (Interpolation)
% 원래 11x5 였던 그리드 공간을 110x50 해상도로 부드럽게 매핑
[Seg_Fine, C_Fine] = meshgrid(linspace(1, N_seg, 110), linspace(target_crates(1), target_crates(end), 50));

Z_Chg_Fine = interp2(Seg_Grid, C_Grid, Z_Chg, Seg_Fine, C_Fine, 'spline');
Z_Dch_Fine = interp2(Seg_Grid, C_Grid, Z_Dch, Seg_Fine, C_Fine, 'spline');

% 스플라인 보간으로 음수가 튀는 것을 방지
Z_Chg_Fine(Z_Chg_Fine < 0) = 0;
Z_Dch_Fine(Z_Dch_Fine < 0) = 0;

%% 4. 고해상도 3D Surface Plotting
fig = figure('Position', [100, 100, 1500, 650], 'Name', 'High-Res 3D Physical Mapping', 'Color', 'w');

% ----------- Subplot 1: Charge -----------
ax1 = subplot(1, 2, 1); hold on; grid on; view(-35, 45);
% Surface 플롯
s1 = surf(Seg_Fine, C_Fine, Z_Chg_Fine);
s1.EdgeColor = 'none';           % 격자선 제거 (매끄러운 표면)
s1.FaceColor = 'interp';         % 색상 부드럽게 보간
s1.FaceLighting = 'gouraud';     % 사실적인 양감 명암 추가
s1.AmbientStrength = 0.4;        
s1.DiffuseStrength = 0.8;

% Contour (등고선) 추가를 통한 가독성 상승
[~, h1] = contour3(Seg_Grid, C_Grid, Z_Chg, 15, 'k', 'LineWidth', 0.5);
h1.ZLocation = 'zmax';           % 등고선을 표면에 붙이거나 위에 표시

colormap(ax1, turbo); 
cb1 = colorbar; 
ylabel(cb1, '\Delta Q Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 11);

% 조명 추가
camlight('headlight'); 
material('dull');

xlabel('Voltage Segments (1 \rightarrow 11)', 'FontWeight', 'bold'); 
ylabel('Current Rate (C-rate)', 'FontWeight', 'bold'); 
zlabel('Capacity \Delta Q (Ah)', 'FontWeight', 'bold');
title('Charge Physical Manifold (\Delta Q_{chg})', 'FontSize', 15, 'FontWeight', 'bold');
xticks(1:11); yticks(target_crates);
zlim([0 max(Z_Chg_Fine(:))*1.1]);


% ----------- Subplot 2: Discharge -----------
ax2 = subplot(1, 2, 2); hold on; grid on; view(-35, 45);
s2 = surf(Seg_Fine, C_Fine, Z_Dch_Fine);
s2.EdgeColor = 'none';           
s2.FaceColor = 'interp';         
s2.FaceLighting = 'gouraud';     
s2.AmbientStrength = 0.4;        
s2.DiffuseStrength = 0.8;

% Contour
[~, h2] = contour3(Seg_Grid, C_Grid, Z_Dch, 15, 'k', 'LineWidth', 0.5);

colormap(ax2, turbo); 
cb2 = colorbar; 
ylabel(cb2, '\Delta Q Capacity (Ah)', 'FontWeight', 'bold', 'FontSize', 11);

% 조명 추가
camlight('headlight'); 
material('dull');

xlabel('Voltage Segments (1 \rightarrow 11)', 'FontWeight', 'bold'); 
ylabel('Current Rate (C-rate)', 'FontWeight', 'bold'); 
zlabel('Capacity \Delta Q (Ah)', 'FontWeight', 'bold');
title('Discharge Physical Manifold (\Delta Q_{dch})', 'FontSize', 15, 'FontWeight', 'bold');
xticks(1:11); yticks(target_crates);
zlim([0 max(Z_Dch_Fine(:))*1.1]);

sgtitle('Deep Learning Physical Mapping Space (High-Res Topography)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);

% Save Figure
saveDir = fileparts(mfilename('fullpath'));
savefig(fig, fullfile(saveDir, 'Physical_Mapping_3D.fig'));
fprintf('>> Saved: Physical_Mapping_3D.fig (High-Res Surface Plot)\n');
