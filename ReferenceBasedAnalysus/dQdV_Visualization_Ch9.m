%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOC_Profile_Visualization_Ch9.m
% Ch9의 0cyc, 200cyc, 400cyc 주행부하 SOC 앵커링 결과 시각화
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% A. 환경 및 경로 설정 (Master Script와 동일)
inputFolder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
ocvDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\RPT\Postprocessing\OCV_integrated'; 

% --- 설정 변수 ---
Cnom_Ah = 64; 
targetChannel = 'Ch9';
cycleTypes = {'0cyc', '200cyc', '400cyc'};
SOCLevels = {'SOC90', 'SOC70', 'SOC50'};
profileName = 'DC1'; % 대표 프로파일 DC1만 시각화

%% B. OCV 데이터 로드
ocvMatFile = fullfile(ocvDataPath, 'OCV_integrated.mat');
if ~exist(ocvMatFile, 'file'); error('OCV_integrated.mat 파일을 찾을 수 없습니다.'); end
load(ocvMatFile, 'OCV_data');
fprintf('OCV_integrated.mat 로드 완료.\n');

%% C. 시각화 루프 (Subplot 사용)

% 3x3 Subplot 구조: Row=SOC Level (90, 70, 50), Col=Cycle Type (0, 200, 400)
figure('Name', sprintf('%s SOC Profile Visualization (Profile %s)', targetChannel, profileName), 'Position', [100 100 1400 800]);

for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    inputFileName = fullfile(inputFolder, sprintf('parsedDriveCycle_%s_filtered.mat', cycleType));
    
    if ~exist(inputFileName, 'file'); continue; end
    load(inputFileName, ['parsedDriveCycle_', cycleType]);
    eval(['currentData = parsedDriveCycle_', cycleType, ';']);
    
    % 해당 사이클의 통합 OCV 데이터 선택
    ocv_struct_name = ['OCV_integrated_', cycleType];
    eval(['OCV_data_current = OCV_data.', ocv_struct_name, ';']);
    
    % 채널 이름 정리 (ch9_Drive_0cyc -> Ch9)
    channelFieldName = sprintf('ch%s_Drive_%s', targetChannel(3:end), cycleType);
    
    if ~isfield(currentData, channelFieldName); continue; end
    channelData = currentData.(channelFieldName);
    
    for socIdx = 1:length(SOCLevels)
        socLevel = SOCLevels{socIdx};
        
        if ~isfield(channelData, socLevel); continue; end
        if ~isfield(channelData.(socLevel), profileName); continue; end
        
        rawProfile = channelData.(socLevel).(profileName);
        
        % 데이터 추출
        V = rawProfile.V;
        I = rawProfile.I;
        t_s = rawProfile.t; % Time [s] (Duration type from parser)

        % 🛠️ 시간 스텝 계산 및 duration to double 변환
        dt_sec = seconds(mean(diff(t_s))); 
        t_double = seconds(t_s) - seconds(t_s(1)); % Time from start [s]
        
        % 🛠️ OCV 앵커링 기반 SOC 프로파일 계산
        try
            [SOC_full, V_ocv_initial, V_ocv_final] = calculate_anchored_soc(V, I, t_s, dt_sec, OCV_data_current);
        catch ME
            fprintf('[오류] %s - %s 계산 실패: %s\n', cycleType, socLevel, ME.message);
            SOC_full = NaN;
        end
        
        % Subplot 인덱스 계산 (Row: SOC level, Col: Cycle Type)
        subplot_idx = (socIdx - 1) * 3 + cycleIdx;
        subplot(3, 3, subplot_idx);
        
        if any(~isnan(SOC_full))
            % SOC Profile 플롯
            plot(t_double, SOC_full, 'LineWidth', 2, 'Color', [0.1 0.4 0.7]);
            
            % 앵커링 지점 표시
            initial_idx = find(abs(diff(SOC_full)) > 1e-6, 1, 'first');
            final_idx = find(abs(diff(SOC_full)) > 1e-6, 1, 'last') + 1;

            if isempty(initial_idx); initial_idx = 1; end
            if isempty(final_idx); final_idx = length(SOC_full); end

            plot(t_double(initial_idx), SOC_full(initial_idx), 'go', 'MarkerSize', 8, 'DisplayName', 'SOC Start');
            plot(t_double(final_idx), SOC_full(final_idx), 'rs', 'MarkerSize', 8, 'DisplayName', 'SOC End');

            % 레이블 및 제목 설정
            ylim([min(soc_grid) - 5, max(soc_grid) + 5]); % OCV_data의 SOC grid 사용
            title(sprintf('%s (%s): %.1f -> %.1f%%', cycleType, socLevel, SOC_full(initial_idx), SOC_full(final_idx)), 'FontSize', 10);
            ylabel('SOC [%]');
            xlabel('Time [s]');
            grid on;
        else
            title(sprintf('%s (%s): 데이터 없음', cycleType, socLevel), 'FontSize', 10);
        end
    end
end
sgtitle(sprintf('Ch9 SOC Profile Analysis (Anchored SOC Method) for Profile %s', profileName), 'FontSize', 14, 'FontWeight', 'bold');


%% D. Sub-Functions (Master Script에서 사용된 로직 복사)

% Sub-Function 1: dQ/dV Analyzer에서 사용된 dQ/dV 파싱의 기초
function [dQdV_AhV, V_mid] = calculate_dQdV_raw(Q_grid, V_ocv)
    dQ = diff(Q_grid); dV = diff(V_ocv);
    dQdV_AhV = dQ ./ dV; dQdV_AhV(abs(dV) < 1e-6) = NaN; 
    V_mid = V_ocv(1:end-1) + dV/2;
end

% Sub-Function 2: OCV 앵커링 기반 SOC 계산 (Master Script에서 복사)
function [SOC_full, V_ocv_initial, V_ocv_final] = calculate_anchored_soc(V_measured, I_measured, t_measured_duration, dt_sec, OCV_data)

    % t_measured_duration은 duration 타입일 수 있으므로 seconds()로 변환
    t_measured = seconds(t_measured_duration);
    
    % --- OCV 데이터 추출 ---
    soc_grid = OCV_data.SOC_grid;
    ocv_values = OCV_data.V_avg_SOC;
    [ocv_values_sorted, uniqueIdx] = unique(ocv_values, 'stable');
    soc_grid_sorted = soc_grid(uniqueIdx);
    inverse_OCV_func = @(voltage) interp1(ocv_values_sorted, soc_grid_sorted, voltage, 'linear');
    
    % --- 휴지기 찾기 ---
    rest_mask = abs(I_measured) <= 2.0;
    min_rest_time_points = floor(5990 * 0.95); % 약 10분 휴지기의 95% (안정성 확보)
    
    % 연속된 휴지기 기간 찾기 (min_rest_time_points는 dt_sec에 따라 달라지지 않도록 포인트 수로 유지)
    rest_periods = []; in_rest = false; rest_start = 0;
    for i = 1:length(rest_mask)
        if rest_mask(i) && ~in_rest; rest_start = i; in_rest = true;
        elseif ~rest_mask(i) && in_rest
            rest_duration_points = (i-1) - rest_start + 1;
            if rest_duration_points >= min_rest_time_points; rest_periods = [rest_periods; rest_start, i-1]; end
            in_rest = false;
        end
    end
    if in_rest; rest_duration_points = length(rest_mask) - rest_start + 1;
        if rest_duration_points >= min_rest_time_points; rest_periods = [rest_periods; rest_start, length(rest_mask)]; end
    end
    
    if size(rest_periods, 1) < 2
        error('시작 및 끝 휴지기 구간을 모두 찾을 수 없습니다.');
    end
    
    % --- SOC 앵커링 ---
    initial_rest_end = rest_periods(1, 2);
    final_rest_end = rest_periods(end, 2);
    
    V_ocv_initial = V_measured(initial_rest_end);
    SOC_initial = inverse_OCV_func(V_ocv_initial);
    
    V_ocv_final = V_measured(final_rest_end);
    SOC_final = inverse_OCV_func(V_ocv_final);
    
    % --- 전류 적산 및 SOC 프로파일 계산 ---
    N = length(I_measured); SOC_full = zeros(N, 1);
    
    cumulative_current_integral = zeros(N, 1);
    for i = (initial_rest_end + 1):final_rest_end
        cumulative_current_integral(i) = cumulative_current_integral(i-1) + I_measured(i) * dt_sec;
    end
    total_current_integration = cumulative_current_integral(final_rest_end);
    
    if abs(total_current_integration) < 1e-6
        SOC_full(:) = SOC_initial; return;
    end
    
    % SOC(t) = SOC1 + (SOC2-SOC1) * (∫₀ᵗ I*dt) / (∫₀ᵉⁿᵈ I*dt) 적용
    for i = 1:N
        if i <= initial_rest_end; SOC_full(i) = SOC_initial;
        elseif i >= final_rest_end; SOC_full(i) = SOC_final;
        else
            SOC_full(i) = SOC_initial + (SOC_final - SOC_initial) * (cumulative_current_integral(i) / total_current_integration);
        end
    end
    
end