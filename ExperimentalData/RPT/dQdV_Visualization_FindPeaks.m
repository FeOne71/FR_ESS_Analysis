% dQdV_Visualization_FindPeaks.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LLI/LAM dQ/dV Quantification Visualization (Final Integrated Version)
% 🛠️ 0cyc, 200cyc, 400cyc 피크 모두 표시 및 디버깅 테이블 출력
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% Paths and Settings
ocvDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';

% --- 피크 검색 파라미터 (사용자가 디버깅 후 조절 필요) ---
V_PEAK_SEARCH_MIN = 3.3; 
V_PEAK_SEARCH_MAX = 3.8; 
MIN_PEAK_PROMINENCE = 0.001; 

%% Load OCV Data
ocvMatFile = fullfile(ocvDataPath, 'OCV_integrated.mat');
if ~exist(ocvMatFile, 'file'); error('OCV_integrated.mat not found.'); end
load(ocvMatFile, 'OCV_data');
fprintf('OCV_integrated.mat 로드 완료.\n');

% --- Extract OCV and Q Grids ---
Q_grid_0 = OCV_data.q_grid_rpt0; V_ocv_0 = OCV_data.avg_ocv_rpt0; 
Q_grid_200 = OCV_data.q_grid_rpt200; V_ocv_200 = OCV_data.avg_ocv_rpt200; 
Q_grid_400 = OCV_data.q_grid_rpt400; V_ocv_400 = OCV_data.avg_ocv_rpt400; 

% --- Capacity Loss Summary (for display) ---
Cap_0 = OCV_data.mean_capacity_rpt0;
Cap_200 = OCV_data.mean_capacity_rpt200;
Cap_400 = OCV_data.mean_capacity_rpt400;

SOH_Loss_0to200 = (Cap_0 - Cap_200) / Cap_0 * 100;
SOH_Loss_200to400 = (Cap_200 - Cap_400) / Cap_200 * 100;

%% 1. dQ/dV Curve Calculation for all cycles
[dQdV_0, V_mid_0] = calculate_dQdV(Q_grid_0, V_ocv_0);
[dQdV_200, V_mid_200] = calculate_dQdV(Q_grid_200, V_ocv_200);
[dQdV_400, V_mid_400] = calculate_dQdV(Q_grid_400, V_ocv_400);

%% 2. 피크 식별 및 LLI/LAM 정량화

% A. 단일 주기 피크 데이터 추출
[V_peak_0, dQdV_peak_0, Peak_Table_0] = find_single_cycle_peak(V_mid_0, dQdV_0, 0, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);
[V_peak_200, dQdV_peak_200, Peak_Table_200] = find_single_cycle_peak(V_mid_200, dQdV_200, V_peak_0, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);
[V_peak_400, dQdV_peak_400, Peak_Table_400] = find_single_cycle_peak(V_mid_400, dQdV_400, V_peak_200, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);

% B. LLI/LAM 정량화 (0 -> 200cyc)
[LLI_shift_V_200, LAM_loss_rate_200] = quantify_lli_lam(V_peak_0, dQdV_peak_0, V_peak_200, dQdV_peak_200);
% C. LLI/LAM 정량화 (200 -> 400cyc)
[LLI_shift_V_400, LAM_loss_rate_400] = quantify_lli_lam(V_peak_200, dQdV_peak_200, V_peak_400, dQdV_peak_400);


%% 3. 디버깅 출력 (Peak Table)
fprintf('\n=======================================================\n');
fprintf('dQ/dV 피크 식별 결과 (MIN_PROMINENCE = %.4f)\n', MIN_PEAK_PROMINENCE);
fprintf('=======================================================\n');
fprintf('--- 0cyc Primary Peak ---\n');
if ~isempty(Peak_Table_0); disp(Peak_Table_0); else fprintf('피크 없음.\n'); end
fprintf('\n--- 200cyc Primary Peak ---\n');
if ~isempty(Peak_Table_200); disp(Peak_Table_200); else fprintf('피크 없음.\n'); end
fprintf('\n--- 400cyc Primary Peak ---\n');
if ~isempty(Peak_Table_400); disp(Peak_Table_400); else fprintf('피크 없음.\n'); end


%% 4. Visualization
figure('Name', 'dQ/dV Degradation Progression (0cyc, 200cyc, 400cyc)', 'Position', [100 100 1000 600]);
hold on; grid on;

% Plot curves (mAh/V)
plot(V_mid_0, dQdV_0 * 1000, 'b-', 'LineWidth', 2, 'DisplayName', '0cyc (Initial)');
plot(V_mid_200, dQdV_200 * 1000, 'r--', 'LineWidth', 2, 'DisplayName', '200cyc (Aged)');
plot(V_mid_400, dQdV_400 * 1000, 'k:', 'LineWidth', 2, 'DisplayName', '400cyc (Most Aged)');


% --- Highlight Tracked Peaks ---
plot_peak(V_peak_0, dQdV_peak_0, 'bo', '0cyc Peak');
plot_peak(V_peak_200, dQdV_peak_200, 'r^', '200cyc Peak');
plot_peak(V_peak_400, dQdV_peak_400, 'k*', '400cyc Peak');


% --- Annotations and Labels ---
title('dQ/dV Curve Comparison: LLI and LAM Progression (Dynamic Peak Tracking)');
xlabel('Voltage [V]');
ylabel('dQ/dV [mAh/V]'); 
legend('Location', 'best');
xlim([min(V_mid_0) max(V_mid_0)]);

% Summary box for Degradation
dim = [.65 .15 .3 .3];
str = {
    sprintf('\\bf\\color{blue} Degradation Summary (%% of Cap Loss)');
    sprintf('\\color{black}SOH Loss (0->200cyc): \\color{red}%.2f %%', SOH_Loss_0to200);
    sprintf('\\color{black}SOH Loss (200->400cyc): \\color{red}%.2f %%', SOH_Loss_200to400);
    '--- LLI/LAM Analysis (0 \rightarrow 200cyc) ---';
    sprintf('\\color{black}LLI (Shift): \\color{red}%.4f V', LLI_shift_V_200);
    sprintf('\\color{black}LAM (Loss): \\color{red}%.2f %%', LAM_loss_rate_200 * 100);
    '--- LLI/LAM Analysis (200 \rightarrow 400cyc) ---';
    sprintf('\\color{black}LLI (Shift): \\color{red}%.4f V', LLI_shift_V_400);
    sprintf('\\color{black}LAM (Loss): \\color{red}%.2f %%', LAM_loss_rate_400 * 100);
};
annotation('textbox', dim, 'String', str, 'FitBoxToText','on', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'tex');

hold off;

%% ========================================================================
% 🌟 서브 함수 1: dQ/dV 곡선 계산
% ========================================================================

function [dQdV_AhV, V_mid] = calculate_dQdV(Q_grid, V_ocv)
    dQ = diff(Q_grid);
    dV = diff(V_ocv);
    dQdV_AhV = dQ ./ dV;
    dQdV_AhV(abs(dV) < 1e-6) = NaN; 
    V_mid = V_ocv(1:end-1) + dV/2;
end

%% ========================================================================
% 🌟 서브 함수 2: LLI/LAM 정량화 (findpeaks 기반)
% ========================================================================

function [LLI_shift_V, LAM_loss_rate] = quantify_lli_lam(V_peak_start, dQdV_peak_start, V_peak_end, dQdV_peak_end)
    % 피크 값들을 받아서 LLI/LAM 정량화 (LLI는 Ah 단위가 아닌 V shift 단위)
    
    LLI_shift_V = V_peak_start - V_peak_end; 

    if dQdV_peak_start > 1e-6
        LAM_loss_rate = (dQdV_peak_start - dQdV_peak_end) / dQdV_peak_start;
    else
        LAM_loss_rate = NaN;
    end
end

%% ========================================================================
% 🌟 서브 함수 3: 단일 주기 핵심 피크 찾기 (findpeaks)
% ========================================================================

function [V_peak_main, dQdV_peak_main, Peak_Table_full] = find_single_cycle_peak(V_mid, dQdV, V_peak_ref_initial, V_SEARCH_MIN, V_SEARCH_MAX, MIN_PROMINENCE)
    
    % 분석 영역을 전압으로 제한
    mask_start = (V_mid >= V_SEARCH_MIN) & (V_mid <= V_SEARCH_MAX);
    V_search = V_mid(mask_start);
    dQdV_search = dQdV(mask_start);
    
    % findpeaks를 사용하여 피크 찾기
    [pks, locs, ~, prom] = findpeaks(dQdV_search, V_search, 'MinPeakProminence', MIN_PROMINENCE);
    
    if isempty(pks)
        V_peak_main = NaN; dQdV_peak_main = NaN; Peak_Table_full = table();
        return;
    end
    
    % 피크 테이블 생성 및 Prominence 기준으로 정렬
    Peak_Table_full = table(locs', pks', prom', 'VariableNames', {'V_Peak', 'dQdV_Peak', 'Prominence'});
    Peak_Table_full = sortrows(Peak_Table_full, 'Prominence', 'descend');
    
    % --- 핵심 피크 선택 로직 ---
    if V_peak_ref_initial > 0 % 이미 이전 주기의 피크 위치를 알고 있는 경우 (200cyc, 400cyc)
        
        % 이전 주기의 피크 전압과 가장 가까운 피크를 선택하여 추적 (LLI 이동 반영)
        [~, nearest_idx] = min(abs(Peak_Table_full.V_Peak - V_peak_ref_initial));
        
        V_peak_main = Peak_Table_full.V_Peak(nearest_idx);
        dQdV_peak_main = Peak_Table_full.dQdV_Peak(nearest_idx);
        
    else % 0cyc (Initial)의 경우, 가장 두드러진 피크를 선택
        V_peak_main = Peak_Table_full.V_Peak(1);
        dQdV_peak_main = Peak_Table_full.dQdV_Peak(1);
    end
end

%% ========================================================================
% 🌟 서브 함수 4: 피크 시각화 헬퍼
% ========================================================================

function plot_peak(V_peak, dQdV_peak, marker_style, label_text)
    if isfinite(V_peak)
        plot(V_peak, dQdV_peak * 1000, marker_style, 'MarkerSize', 8, 'MarkerFaceColor', marker_style(1));
        % text 함수는 dQdV_peak_start * 1000을 사용해야 mAh/V 단위에 맞습니다.
        text(V_peak * 1.001, dQdV_peak * 1000 * 1.05, label_text, 'FontSize', 9, 'Interpreter', 'none');
    end
end