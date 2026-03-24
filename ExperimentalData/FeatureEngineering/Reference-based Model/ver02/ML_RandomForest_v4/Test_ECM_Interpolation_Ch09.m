%% Test_ECM_Interpolation_Ch09.m
% =====================================================================
% 질문: "C-rate마다 데이터 파이(길이)가 다를 텐데 어떻게 똑같이 자를 수 있나요?"
% 답변 시각화: 내분점 보간법 (Interpolation) 데모 스크립트
% =====================================================================
clear; clc; close all;

%% 1. 데이터 로드 (Cell 09, Cycle 0 데이터 세트)
vq_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat';
ecm_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\ECM_2RC_AllChannels_AllCycles.mat';

load(vq_path, 'RPT_VQ_grid');
load(ecm_path, 'All_ECM');

ch = 'Ch09';
c01_data = RPT_VQ_grid.cyc0.(ch).c01_charge;
c30_data = RPT_VQ_grid.cyc0.(ch).c3_charge;

% ECM 정보
ecm_chg = All_ECM.(ch).cyc0.charge;

%% 2. Raw 전압과 ECM 보정 전압(Vc) 계산
% --- 0.1C ---
V_raw_c01 = double(c01_data.V_raw(:));
I_raw_c01 = double(c01_data.I_raw(:));
t_s_c01   = seconds(c01_data.t_raw(:));
Q_raw_c01 = cumtrapz(t_s_c01, abs(I_raw_c01))/3600;
Q_0 = max(Q_raw_c01);

soc_c01 = Q_raw_c01 / Q_0 * 100;
Ri_c01 = interp1(ecm_chg.SOC, (ecm_chg.R0+ecm_chg.R1+ecm_chg.R2)/1e3, soc_c01, 'linear','extrap');
V_c_c01 = V_raw_c01 - I_raw_c01 .* Ri_c01;

% --- 3.0C ---
V_raw_c30 = double(c30_data.V_raw(:));
I_raw_c30 = double(c30_data.I_raw(:));
t_s_c30   = seconds(c30_data.t_raw(:));
Q_raw_c30 = cumtrapz(t_s_c30, abs(I_raw_c30))/3600;

soc_c30 = Q_raw_c30 / Q_0 * 100;
soc_c30 = max(0, min(100, soc_c30)); % 범위 초과 방지
Ri_c30 = interp1(ecm_chg.SOC, (ecm_chg.R0+ecm_chg.R1+ecm_chg.R2)/1e3, soc_c30, 'linear','extrap');
V_c_c30 = V_raw_c30 - I_raw_c30 .* Ri_c30; 

%% 3. 배열 길이 비교 출력
fprintf('=== 데이터 배열 길이 비교 ===\n');
fprintf('0.1C 데이터 길이 (개수): %d 포인트\n', length(V_c_c01));
fprintf('3.0C 데이터 길이 (개수): %d 포인트\n', length(V_c_c30));
fprintf('* 3.0C는 빨리 충전되므로 센서 포인트 개수가 더 적고 최대 충전 용량도 다릅니다.\n\n');

%% 4. Interpolation 데모 (고정 좌표에서의 추출)
% 기준 전압 경계 (예: 3.65V, 3.75V, 3.85V)
V_target = [3.65, 3.75, 3.85];

% 0.1C 고유 전압에 대한 해당 점의 Q값 찾기 (중복제거 필수)
[V_c_c01_u, uid] = unique(V_c_c01, 'stable'); Q_c01_u = Q_raw_c01(uid);
[V_c_c30_u, uid] = unique(V_c_c30, 'stable'); Q_c30_u = Q_raw_c30(uid);

% 1. 데이터 배열 길이가 달라도, y=f(x) 선상에서 x=V_target의 y좌표(Q)를 내분점으로 추정
Q_target_01 = interp1(V_c_c01_u, Q_c01_u, V_target, 'linear');
Q_target_30 = interp1(V_c_c30_u, Q_c30_u, V_target, 'linear');

fprintf('=== Interpolation 추출 원리 ===\n');
for i=1:3
    fprintf('Target Vc = %.2f V 에 대응하는 용량 (Q):\n', V_target(i));
    fprintf('  0.1C: %.4f Ah\n', Q_target_01(i));
    fprintf('  3.0C: %.4f Ah\n', Q_target_30(i));
end
fprintf('[핵심] 배열(Array) 자체를 자르는게 아니라, 두 연속된 점 사이의 V_target 값을 1차함수로 계산(보간)하여 척출합니다.\n');

%% 5. 시각화
fig = figure('Name', 'Interpolation Concept', 'Position', [100, 100, 1000, 500]);

subplot(1,2,1); hold on; grid on;
plot(Q_raw_c01, V_raw_c01, 'b-', 'LineWidth', 1.5, 'DisplayName', '0.1C Raw');
plot(Q_raw_c30, V_raw_c30, 'r-', 'LineWidth', 1.5, 'DisplayName', '3.0C Raw');
xlabel('Capacity Q (Ah)'); ylabel('Voltage (V)');
title('Raw 전압 (IR Drop 발생 - 오차 큼)');
legend('Location', 'southeast');

subplot(1,2,2); hold on; grid on;
plot(Q_c01_u, V_c_c01_u, 'b-', 'LineWidth', 1.5, 'DisplayName', '0.1C ECM (Vc)');
plot(Q_c30_u, V_c_c30_u, 'r-', 'LineWidth', 1.5, 'DisplayName', '3.0C ECM (Vc)');

% Interpolation Point 표시 (y축이 V이므로 수평선으로 그림)
for i=1:3
    scatter(Q_target_01(i), V_target(i), 80, 'b', 'filled', 'HandleVisibility', 'off');
    scatter(Q_target_30(i), V_target(i), 80, 'r', 'filled', 'HandleVisibility', 'off');
    % 수평선으로 Target 전압 표시
    yline(V_target(i), 'k--', 'HandleVisibility', 'off');
end

xlabel('Capacity Q (Ah)'); ylabel('Compensated Voltage Vc (V)');
title('ECM 보정 전압 (y좌표=V 기준, x좌표=Q 추출)');
legend('Location', 'southeast');

saveDir = fileparts(mfilename('fullpath'));
saveas(fig, fullfile(saveDir, 'Interpolation_Concept_Ch09.fig'));
fprintf('\n>> Saved: Interpolation_Concept_Ch09.fig\n');
