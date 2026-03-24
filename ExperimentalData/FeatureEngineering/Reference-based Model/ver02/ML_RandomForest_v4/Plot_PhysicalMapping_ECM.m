%% Plot_PhysicalMapping_ECM.m
% =====================================================================
% Phase 4: 물리적 매핑 (Physical Mapping) - ECM 보정의 진정한 시각적 증명
% 딥러닝 피처 공간이 아니라, 배터리 물리(Physics) 관점의 매핑 공간입니다.
% 
% - X축: 용량 (Q, Ah)
% - Y축: C-rate (0.1C ~ 3.0C)
% - Z축: 전압 (V)
% 
% (1) ECM 보정 전: Raw 단자 전압 (Vt) -> C-rate가 커질수록 IR Drop으로 인해 찌그러짐
% (2) ECM 보정 후: 보정 개방회로전압 (V_ocv) -> 모든 C-rate가 하나의 평행한 곡면으로 정렬
% =====================================================================
clear; clc; close all;

%% 1. 데이터 로드 준비
% VQ Grid 로드
vq_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat';
fprintf('Loading VQ data...\n');
d_vq = load(vq_path);

% ECM 파라미터 로드
ecm_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\ECM_2RC_AllChannels_AllCycles.mat';
fprintf('Loading ECM parameters...\n');
d_ecm = load(ecm_path);

ch = 'Ch09';
cyc_key = 'cyc0'; % 기준 열역학적 평형 상태 관측용
ch_data = d_vq.RPT_VQ_grid.(cyc_key).(ch);

ecm_chg = d_ecm.All_ECM.(ch).(cyc_key).charge;
ecm_dch = d_ecm.All_ECM.(ch).(cyc_key).discharge;

crates_str = {'C01', 'C05', 'C1', 'C2', 'C3'};
crate_vals = [0.1, 0.5, 1.0, 2.0, 3.0];
N_c = length(crate_vals);

% Capacity Grid (공통 X축으로 보간하기 위한 격자)
Q_grid_chg = linspace(10, 60, 150); % 충전 용량 (Ah)
Q_grid_dch = linspace(10, 60, 150); % 방전 용량 (Ah)

V_raw_chg = nan(N_c, length(Q_grid_chg)); V_ocv_chg = nan(N_c, length(Q_grid_chg));
V_raw_dch = nan(N_c, length(Q_grid_dch)); V_ocv_dch = nan(N_c, length(Q_grid_dch));

% 스무딩 윈도우 크기 (C-rate 1초단위 보정용)
sm_win = 100;

%% 2. 데이터 추출 및 ECM C-rate별 보간
for r = 1:N_c
    % ---------- CHARGE ----------
    f_chg = [crates_str{r} '_charge'];
    if isfield(ch_data, f_chg)
        raw_c = ch_data.(f_chg);
        V_c = raw_c(:,2); I_c = raw_c(:,3); Q_c = raw_c(:,4);
        
        % 중복 Q 제거
        [Q_c_u, idx_u] = unique(Q_c);
        V_c_u = V_c(idx_u); I_c_u = I_c(idx_u);
        
        % V_raw 보간
        V_raw_chg(r, :) = interp1(Q_c_u, V_c_u, Q_grid_chg, 'linear', 'extrap');
        
        % V_ocv (ECM 보정) - Q_0는 임의 공칭값(68) 삽입, Vc만 활용
        Q0 = 68;
        [Vc_c, ~] = ecm_correct(raw_c, ecm_chg, Q0, true, sm_win);
        [Q_c_u2, idx_u2] = unique(raw_c(:,4));
        Vc_c_u = Vc_c(idx_u2);
        V_ocv_chg(r, :) = interp1(Q_c_u2, Vc_c_u, Q_grid_chg, 'linear', 'extrap');
    end
    
    % ---------- DISCHARGE ----------
    f_dch = [crates_str{r} '_discharge'];
    if isfield(ch_data, f_dch)
        raw_d = ch_data.(f_dch);
        V_d = raw_d(:,2); I_d = raw_d(:,3); Q_d = raw_d(:,4);
        
        [Q_d_u, idx_u] = unique(Q_d);
        V_d_u = V_d(idx_u); I_d_u = I_d(idx_u);
        
        V_raw_dch(r, :) = interp1(Q_d_u, V_d_u, Q_grid_dch, 'linear', 'extrap');
        
        [Vc_d, ~] = ecm_correct(raw_d, ecm_dch, Q0, false, sm_win);
        [Q_d_u2, idx_u2] = unique(raw_d(:,4));
        Vc_d_u = Vc_d(idx_u2);
        V_ocv_dch(r, :) = interp1(Q_d_u2, Vc_d_u, Q_grid_dch, 'linear', 'extrap');
    end
end

% 튀는 값이나 스플라인 외삽 오류 자르기
V_raw_chg(V_raw_chg > 4.3 | V_raw_chg < 3.0) = NaN;
V_ocv_chg(V_ocv_chg > 4.3 | V_ocv_chg < 3.0) = NaN;
V_raw_dch(V_raw_dch > 4.3 | V_raw_dch < 3.0) = NaN;
V_ocv_dch(V_ocv_dch > 4.3 | V_ocv_dch < 3.0) = NaN;


%% 3. Surface 3D 시각화 (진정한 Physical Mapping 증명)
[Q_Grid_C, C_Grid_C] = meshgrid(Q_grid_chg, crate_vals);
[Q_Grid_D, C_Grid_D] = meshgrid(Q_grid_dch, crate_vals);

fig = figure('Position', [100, 100, 1600, 700], 'Name', 'ECM Physical Mapping: Raw vs OCV');

% -------- Chg: Raw Vt --------
ax1 = subplot(2,2,1); hold on; grid on; view(-30, 30);
surf(Q_Grid_C, C_Grid_C, V_raw_chg, 'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud');
colormap(ax1, parula); camlight('headlight'); material('dull');
title('Charge: Before ECM (Terminal Voltage V_t)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Capacity (Ah)'); ylabel('C-rate'); zlabel('Voltage (V)');
zlim([3.4 4.25]);

% -------- Chg: V_ocv --------
ax2 = subplot(2,2,3); hold on; grid on; view(-30, 30);
surf(Q_Grid_C, C_Grid_C, V_ocv_chg, 'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud');
colormap(ax2, parula); camlight('headlight'); material('dull');
title('Charge: After ECM (Estimated V_{OCV}) -> Physical Mapped', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Capacity (Ah)'); ylabel('C-rate'); zlabel('Voltage (V)');
zlim([3.4 4.25]);

% -------- Dch: Raw Vt --------
ax3 = subplot(2,2,2); hold on; grid on; view(-30, 30);
surf(Q_Grid_D, C_Grid_D, V_raw_dch, 'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud');
colormap(ax3, parula); camlight('headlight'); material('dull');
title('Discharge: Before ECM (Terminal Voltage V_t)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Capacity (Ah)'); ylabel('C-rate'); zlabel('Voltage (V)');
zlim([3.4 4.25]);

% -------- Dch: V_ocv --------
ax4 = subplot(2,2,4); hold on; grid on; view(-30, 30);
surf(Q_Grid_D, C_Grid_D, V_ocv_dch, 'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud');
colormap(ax4, parula); camlight('headlight'); material('dull');
title('Discharge: After ECM (Estimated V_{OCV}) -> Physical Mapped', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Capacity (Ah)'); ylabel('C-rate'); zlabel('Voltage (V)');
zlim([3.4 4.25]);

sgtitle('The Physical Evidence of ECM: Removing IR Drop to Map Thermodynamic OCV Surface (Cell09, Cyc0)', ...
    'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);

%% 4. 보조 함수 선언
saveDir = fileparts(mfilename('fullpath'));
savefig(fig, fullfile(saveDir, 'Physical_Mapping_ECM.fig'));
fprintf('\n>> Saved: Physical_Mapping_ECM.fig\n');


function [Vc, Q_interp] = ecm_correct(vq_raw, ecm, Qt, isChg, sm_win)
    % ecm_correct: V_raw에서 IR Drop을 빼는 함수 (PCC_Analysis_v5b.m 재사용)
    t = vq_raw(:,1); V = vq_raw(:,2); I = vq_raw(:,3); Q = vq_raw(:,4);
    
    soc = Q / Qt * 100;
    soc(soc < 0) = 0; soc(soc > 100) = 100;
    if ~isChg, soc = 100 - soc; end
    
    r0 = interp1(ecm.soc, ecm.r0, soc, 'linear', 'extrap');
    r1 = interp1(ecm.soc, ecm.r1, soc, 'linear', 'extrap');
    c1 = interp1(ecm.soc, ecm.c1, soc, 'linear', 'extrap');
    
    V1 = zeros(size(t));
    dt = diff(t);
    for k = 1:length(t)-1
        tau = r1(k) * c1(k);
        if tau == 0, V1(k+1) = V1(k); continue; end
        exp_factor = exp(-dt(k) / tau);
        V1(k+1) = V1(k) * exp_factor + r1(k) * I(k) * (1 - exp_factor);
    end
    
    if isChg
        Vc = V - I.*r0 - V1;
    else
        Vc = V + abs(I).*r0 + abs(V1); 
    end
    
    Vc = smoothdata(Vc, 'movmean', sm_win);
    Q_interp = Q;
end
