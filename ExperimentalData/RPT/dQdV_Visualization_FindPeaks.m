% dQdV_Visualization_FindPeaks_v2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dQ/dV 분석 및 LLI/LAM 정량화 스크립트 (RPT_Postprocessing.m 연동 버전)
%
% Description:
% 1. RPT_Postprocessing.m 에서 생성된 OCV_integrated.mat 파일을 로드합니다.
% 2. 파일 내에 존재하는 모든 RPT 주기(0, 200, 400cyc...) 데이터를 동적으로 인식합니다.
% 3. 각 주기에 대해 dQ/dV 곡선을 계산합니다.
% 4. findpeaks 함수를 이용하여 dQ/dV 곡선의 주요 피크를 식별하고, 사이클이 진행됨에 따라 피크를 추적합니다.
% 5. 추적된 피크의 변화(전압 이동, 높이 감소)를 통해 LLI와 LAM을 정량화합니다.
% 6. 모든 dQ/dV 곡선과 정량화 결과를 하나의 그래프에 시각화하여 열화 경향을 명확하게 보여줍니다.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
warning off;

%% ========================================================================
%  1. 설정 및 데이터 로딩
% =========================================================================

% --- 경로 설정 ---
% RPT_Postprocessing.m 에서 OCV_integrated.mat 파일이 저장된 폴더 경로를 지정하세요.
ocvDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';

% --- 피크 검색 파라미터 (분석할 배터리 종류에 따라 사용자가 조절) ---
V_PEAK_SEARCH_MIN = 3.3;  % 피크를 탐색할 최소 전압 [V]
V_PEAK_SEARCH_MAX = 3.8;  % 피크를 탐색할 최대 전압 [V]
MIN_PEAK_PROMINENCE = 0.001; % 피크로 인정할 최소 돌출값 (노이즈 제거용)

% --- OCV 데이터 로딩 ---
ocvMatFile = fullfile(ocvDataPath, 'OCV_integrated.mat');
if ~exist(ocvMatFile, 'file')
    error('지정된 경로에 OCV_integrated.mat 파일이 없습니다. RPT_Postprocessing.m을 먼저 실행하세요.');
end
load(ocvMatFile, 'OCV_data');
fprintf('OCV_integrated.mat 로드 완료.\n');


%% ========================================================================
%  2. 데이터 동적 추출 및 구성
%  (OCV_data 구조체에서 모든 사이클 데이터를 자동으로 읽어옵니다)
% =========================================================================

% OCV_data 구조체의 필드 이름을 분석하여 사용 가능한 사이클 목록을 찾습니다.
all_fields = fieldnames(OCV_data);
q_grid_fields = all_fields(startsWith(all_fields, 'q_grid_rpt'));

% 'q_grid_rpt' 다음의 숫자(사이클)를 추출하여 cell array로 만듭니다.
cycle_keys_str = cellfun(@(s) s(11:end), q_grid_fields, 'UniformOutput', false);

% 문자열인 사이클 번호를 숫자로 변환하여 오름차순으로 정렬합니다. (e.g., '0', '1000', '200' -> 0, 200, 1000)
[~, sort_idx] = sort(cellfun(@str2double, cycle_keys_str));
cycle_keys = cycle_keys_str(sort_idx); % 정렬된 사이클 키 ('0', '200', '400', '600')

fprintf('분석 대상 사이클: %s\n', strjoin(cycle_keys, ', '));

% 각 사이클별 데이터를 저장할 구조체 초기화
Q_data = struct();
V_data = struct();
Cap_data = struct();

% 정렬된 사이클 순서대로 OCV_data에서 Q, V, Capacity 값을 추출하여 새로운 구조체에 저장
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    Q_data.(['c' key]) = OCV_data.(['q_grid_rpt' key]);
    V_data.(['c' key]) = OCV_data.(['avg_ocv_rpt' key]);
    Cap_data.(['c' key]) = OCV_data.(['mean_capacity_rpt' key]);
end

%% ========================================================================
%  3. dQ/dV 계산 및 열화 정량화
% =========================================================================

% --- 1. 모든 사이클에 대해 dQ/dV 곡선 계산 ---
dQdV_results = struct();
V_mid_points = struct();
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    [dQdV_results.(['c' key]), V_mid_points.(['c' key])] = calculate_dQdV(Q_data.(['c' key]), V_data.(['c' key]));
end
fprintf('모든 사이클의 dQ/dV 계산 완료.\n');

% --- 2. 피크 식별 및 LLI/LAM 정량화 ---
V_peaks = struct();
dQdV_peaks = struct();
Peak_Tables = struct();
LLI_LAM_results = struct();

% 첫 번째 사이클 (기준)의 피크 찾기
base_key = cycle_keys{1};
[V_peaks.(['c' base_key]), dQdV_peaks.(['c' base_key]), Peak_Tables.(['c' base_key])] = ...
    find_main_peak(V_mid_points.(['c' base_key]), dQdV_results.(['c' base_key]), 0, V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);

% 이후 사이클들의 피크를 '추적'하며 LLI/LAM 계산
for i = 2:length(cycle_keys)
    prev_key = cycle_keys{i-1};
    curr_key = cycle_keys{i};
    
    % 이전 사이클의 피크 위치(V_peaks)를 참조하여 현재 사이클의 피크를 찾음 (피크 추적)
    [V_peaks.(['c' curr_key]), dQdV_peaks.(['c' curr_key]), Peak_Tables.(['c' curr_key])] = ...
        find_main_peak(V_mid_points.(['c' curr_key]), dQdV_results.(['c' curr_key]), V_peaks.(['c' prev_key]), V_PEAK_SEARCH_MIN, V_PEAK_SEARCH_MAX, MIN_PEAK_PROMINENCE);

    % LLI/LAM 정량화
    [lli_V, lam_rate] = quantify_lli_lam(V_peaks.(['c' prev_key]), dQdV_peaks.(['c' prev_key]), V_peaks.(['c' curr_key]), dQdV_peaks.(['c' curr_key]));
    LLI_LAM_results.(sprintf('LLI_V_%s_to_%s', prev_key, curr_key)) = lli_V;
    LLI_LAM_results.(sprintf('LAM_rate_%s_to_%s', prev_key, curr_key)) = lam_rate;
end
fprintf('피크 추적 및 LLI/LAM 정량화 완료.\n');

%% ========================================================================
%  4. 결과 출력 및 시각화
% =========================================================================

% --- 1. 디버깅용 피크 테이블 출력 ---
fprintf('\n=======================================================\n');
fprintf('dQ/dV 피크 식별 결과 (MinPeakProminence = %.4f)\n', MIN_PEAK_PROMINENCE);
fprintf('=======================================================\n');
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    fprintf('\n--- %s cyc Primary Peak ---\n', key);
    if ~isempty(Peak_Tables.(['c' key])); disp(Peak_Tables.(['c' key])); else; fprintf('피크를 찾지 못했습니다.\n'); end
end

% --- 2. 종합 그래프 시각화 ---
figure('Name', 'dQ/dV Degradation Analysis', 'Position', [100 100 1200 700]);
hold on; grid on;

% dQ/dV 곡선 플롯
colors = lines(length(cycle_keys)); % 사이클 개수에 맞춰 자동으로 색상 생성
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    plot(V_mid_points.(['c' key]), dQdV_results.(['c' key]) * 1000, ...
        'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', sprintf('%s cyc', key));
end

% 추적된 피크 위치 표시
for i = 1:length(cycle_keys)
    key = cycle_keys{i};
    plot_peak(V_peaks.(['c' key]), dQdV_peaks.(['c' key]), 'o', colors(i,:), sprintf('%s cyc Peak', key));
end

% 그래프 제목 및 레이블
title('dQ/dV Curve Comparison: Degradation Progression');
xlabel('Voltage [V]');
ylabel('dQ/dV [mAh/V]');
legend('Location', 'northeast');
xlim([V_PEAK_SEARCH_MIN-0.1 V_PEAK_SEARCH_MAX+0.1]);

% --- 3. 정량화 결과 요약 텍스트 박스 추가 ---
summary_str = {'\bfDegradation Analysis Summary'};
summary_str{end+1} = ' ';

for i = 2:length(cycle_keys)
    prev_key = cycle_keys{i-1};
    curr_key = cycle_keys{i};
    
    % SOH 감소율
    soh_loss = (Cap_data.(['c' prev_key]) - Cap_data.(['c' curr_key])) / Cap_data.(['c' base_key]) * 100;
    summary_str{end+1} = sprintf('\\bf\\color[rgb]{0,0.4,0.8}--- [%s → %s cyc] ---', prev_key, curr_key);
    summary_str{end+1} = sprintf('\\color{red}Total Capacity Loss: %.2f %%', soh_loss);
    
    % LLI / LAM
    lli_V = LLI_LAM_results.(sprintf('LLI_V_%s_to_%s', prev_key, curr_key));
    lam_rate = LLI_LAM_results.(sprintf('LAM_rate_%s_to_%s', prev_key, curr_key));
    summary_str{end+1} = sprintf('LLI (Voltage Shift): %.4f V', lli_V);
    summary_str{end+1} = sprintf('LAM (Peak-Height Loss): %.2f %%', lam_rate * 100);
end

annotation('textbox', [0.15, 0.15, 0.3, 0.3], 'String', summary_str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'tex');

hold off;

%% ========================================================================
%                         🌟 서브 함수 (Helper Functions) 🌟
% =========================================================================

% 서브 함수 1: dQ/dV 곡선 계산
function [dQdV_AhV, V_mid] = calculate_dQdV(Q_grid_Ah, V_ocv)
    % Q와 V 데이터로부터 dQ/dV를 계산합니다.
    dQ = diff(Q_grid_Ah);
    dV = diff(V_ocv);
    
    % dV가 0에 가까우면 dQ/dV가 무한대로 튀는 것을 방지
    valid_indices = abs(dV) > 1e-6; 
    
    dQdV_AhV = NaN(size(dV)); % 결과를 NaN으로 초기화
    dQdV_AhV(valid_indices) = dQ(valid_indices) ./ dV(valid_indices);
    
    % dQ/dV 값에 해당하는 중간 전압 계산
    V_mid = V_ocv(1:end-1) + dV/2;
end

% 서브 함수 2: 주요 피크 찾기
function [V_peak_main, dQdV_peak_main, Peak_Table] = find_main_peak(V_mid, dQdV, V_peak_ref, V_SEARCH_MIN, V_SEARCH_MAX, MIN_PROMINENCE)
    % dQ/dV 곡선에서 가장 의미있는 피크를 찾습니다.
    
    % 유효한 dQ/dV 값과 전압 범위 내의 데이터만 필터링
    valid_mask = ~isnan(dQdV) & (V_mid >= V_SEARCH_MIN) & (V_mid <= V_SEARCH_MAX);
    V_search = V_mid(valid_mask);
    dQdV_search = dQdV(valid_mask);
    
    if isempty(V_search)
        V_peak_main = NaN; dQdV_peak_main = NaN; Peak_Table = table();
        return;
    end
    
    % findpeaks 함수로 피크 검색
    [pks, locs, ~, prom] = findpeaks(dQdV_search, V_search, 'MinPeakProminence', MIN_PROMINENCE);
    
    if isempty(pks)
        V_peak_main = NaN; dQdV_peak_main = NaN; Peak_Table = table();
        return;
    end
    
    Peak_Table = table(locs', pks', prom', 'VariableNames', {'V_Peak', 'dQdV_Peak', 'Prominence'});
    Peak_Table = sortrows(Peak_Table, 'Prominence', 'descend'); % 중요도 순으로 정렬
    
    % --- 핵심 피크 선택 로직 ---
    if V_peak_ref > 0 % 이전 사이클의 피크 위치(V_peak_ref) 정보가 있는 경우
        % 이전 피크와 전압이 가장 가까운 피크를 현재 피크로 선택 (피크 추적)
        [~, nearest_idx] = min(abs(Peak_Table.V_Peak - V_peak_ref));
        V_peak_main = Peak_Table.V_Peak(nearest_idx);
        dQdV_peak_main = Peak_Table.dQdV_Peak(nearest_idx);
    else % 첫 사이클 (기준)인 경우
        % 가장 높은 dQ/dV 값을 가진 피크를 메인 피크로 선택 (prominence가 아닌 절대값 기준)
        [~, max_idx] = max(Peak_Table.dQdV_Peak);
        V_peak_main = Peak_Table.V_Peak(max_idx);
        dQdV_peak_main = Peak_Table.dQdV_Peak(max_idx);
    end
end

% 서브 함수 3: LLI/LAM 정량화
function [LLI_shift_V, LAM_loss_rate] = quantify_lli_lam(V_peak_start, dQdV_peak_start, V_peak_end, dQdV_peak_end)
    % 두 시점의 피크 정보를 바탕으로 LLI와 LAM을 계산합니다.
    
    % LLI: 피크의 전압 이동량 [V]
    LLI_shift_V = V_peak_start - V_peak_end; 

    % LAM: 피크 높이의 감소율 [%]
    if dQdV_peak_start > 1e-6
        LAM_loss_rate = (dQdV_peak_start - dQdV_peak_end) / dQdV_peak_start;
    else
        LAM_loss_rate = NaN; % 기준 피크 높이가 0에 가까우면 계산 불가
    end
end

% 서브 함수 4: 피크 시각화
function plot_peak(V_peak, dQdV_peak, marker, color, label_text)
    % 그래프에 피크 위치와 레이블을 표시합니다.
    if isfinite(V_peak) && isfinite(dQdV_peak)
        plot(V_peak, dQdV_peak * 1000, marker, 'MarkerSize', 10, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k');
        % 텍스트를 피크 바로 오른쪽에 위치
        text(V_peak + 0.01, dQdV_peak * 1000, label_text, 'FontSize', 9, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end