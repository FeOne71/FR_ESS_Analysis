%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOC_Profile_Visualization_Ch9.m
% Ch9의 0cyc, 200cyc, 400cyc 주행부하 SOC 앵커링 결과 시각화
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% A. 환경 및 경로 설정 (Master Script와 동일)
inputFolder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
ocvDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated'; 

% --- 설정 변수 ---
targetChannels = {'Ch10'};
cycleTypes = {'0cyc', '200cyc', '400cyc'};
SOCLevels = {'SOC90', 'SOC70', 'SOC50'};
profileName = 'DC1'; % 대표 프로파일 DC1만 시각화

%% B. OCV 데이터 로드
ocvMatFile = fullfile(ocvDataPath, 'OCV_integrated.mat');
if ~exist(ocvMatFile, 'file'); error('OCV_integrated.mat 파일을 찾을 수 없습니다.'); end
load(ocvMatFile, 'OCV_data');
fprintf('OCV_integrated.mat 로드 완료.\n');

if ~isfield(OCV_data, 'OCV_integrated_0cyc')
    error('OCV_data 구조체에 OCV_integrated_0cyc 필드가 없어 SOC_grid를 찾을 수 없습니다.');
end
soc_grid = OCV_data.OCV_integrated_0cyc.SOC_grid; 

%% C. 시각화 루프 (CycleType별 Figure 생성)

for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    inputFileName = fullfile(inputFolder, sprintf('parsedDriveCycle_%s_filtered.mat', cycleType));
    
    if ~exist(inputFileName, 'file'); continue; end
    load(inputFileName, ['parsedDriveCycle_', cycleType]);
    eval(['currentData = parsedDriveCycle_', cycleType, ';']);
    
    ocv_struct_name = ['OCV_integrated_', cycleType];
    if ~isfield(OCV_data, ocv_struct_name); continue; end
    eval(['OCV_data_current = OCV_data.', ocv_struct_name, ';']);
    
    % CycleType별 새로운 Figure 생성 (3x1 Subplot)
    figure('Name', sprintf('%s Cycle - SOC Profile Analysis', cycleType), 'Position', [100 100 800 1000]); % Figure 크기 조정
    sgtitle(sprintf('%s Cycle SOC/Current Profile Analysis (Profile %s)', cycleType, profileName), 'FontSize', 12, 'FontWeight', 'bold');
    
    for socIdx = 1:length(SOCLevels)
        socLevel = SOCLevels{socIdx};
        subplot(3, 1, socIdx);
        
        valid_channels_found = 0;
        
        % --- 플롯 초기화 ---
        cla; % Subplot 초기화
        yyaxis left; hold on; grid on;
        
        % 🛠️ SOC Y축 범위 동적 설정
        SOC_Nominal = str2double(regexp(socLevel, '\d+', 'match', 'once')); % 'SOC90' -> 90
        ylim_min = SOC_Nominal - 10;
        ylim_max = SOC_Nominal + 10;
        ylim([ylim_min, ylim_max]);
        ylabel('SOC [%] (OCV Anchored)', 'Color', 'b');
        
        yyaxis right; hold on;
        ylabel('Current [A]', 'Color', 'r');
        
        % --- 채널 루프: 모든 채널 순회 및 오버레이 플롯 ---
        for chIdx = 1:length(targetChannels)
            channel = targetChannels{chIdx};
            ch_num = channel(3:end);
            channelFieldName = sprintf('ch%s_Drive_%s', ch_num, cycleType); 
            
            if isfield(currentData, channelFieldName)
                channelData = currentData.(channelFieldName);
                if isfield(channelData, socLevel) && isfield(channelData.(socLevel), profileName)
                    rawProfile = channelData.(socLevel).(profileName);
                    
                    V = rawProfile.V; I = rawProfile.I; t_s = rawProfile.t; 
                    dt_sec = seconds(mean(diff(t_s))); 
                    t_double = seconds(t_s) - seconds(t_s(1)); 
                    
                    try
                        [SOC_full, ~, ~, initial_idx, final_idx] = calculate_anchored_soc_v2(V, I, t_s, dt_sec, OCV_data_current);
                    catch ME
                        continue;
                    end
                    
                    if any(~isnan(SOC_full))
                        % 🛠️ FIX 1: SOC Profile (yyaxis left)
                        yyaxis left;
                        plot(t_double, SOC_full, 'LineWidth', 1, 'DisplayName', channel);
                        
                        % 🛠️ FIX 2: Current Profile (yyaxis right)
                        yyaxis right;
                        plot(t_double, I, 'r-', 'LineWidth', 0.5, 'HandleVisibility','off'); % 전류는 얇게 표시
                        
                        % 앵커링 지점 표시 (Left Axis)
                        yyaxis left;
                        plot(t_double(initial_idx), SOC_full(initial_idx), 'go', 'MarkerSize', 5, 'HandleVisibility','off');
                        plot(t_double(final_idx), SOC_full(final_idx), 'ms', 'MarkerSize', 5, 'HandleVisibility','off');
                        
                        valid_channels_found = valid_channels_found + 1;
                    end
                end
            end
        end % End Channel Loop
        
        % --- 최종 제목 및 레이블 설정 ---
        yyaxis right; % Right axis max/min 설정 (전류 범위 확보)
        current_max = max(abs(I));
        if current_max < 10; current_max = 10; end % 최소 전류 범위 확보
        ylim([-current_max, current_max]);
        
        yyaxis left;
        title(sprintf('%s Level (Start: %.1f%%)', socLevel, SOC_Nominal), 'FontSize', 10);
        xlabel('Time [s]');
        legend('Location', 'bestoutside', 'FontSize', 8);
        hold off;
    end % End SOC Level Loop
end % End Cycle Type Loop

%% D. Sub-Functions (로직 복사 및 수정)

% Sub-Function 1: Raw dQ/dV Calculation (used by Y-Label, not this visualization)
function [dQdV_AhV, V_mid] = calculate_dQdV_raw(Q_grid, V_ocv)
    dQ = diff(Q_grid); dV = diff(V_ocv);
    dQdV_AhV = dQ ./ dV; dQdV_AhV(abs(dV) < 1e-6) = NaN; 
    V_mid = V_ocv(1:end-1) + dV/2;
end

% Sub-Function 2: OCV 앵커링 기반 SOC 계산 (V2)
function [SOC_full, V_ocv_initial, V_ocv_final, initial_rest_end, final_rest_end] = calculate_anchored_soc_v2(V_measured, I_measured, t_measured_duration, dt_sec, OCV_data)

    t_measured = seconds(t_measured_duration);
    
    % --- OCV 데이터 추출 ---
    soc_grid = OCV_data.SOC_grid;
    ocv_values = OCV_data.V_avg_SOC;
    [ocv_values_sorted, uniqueIdx] = unique(ocv_values, 'stable');
    soc_grid_sorted = soc_grid(uniqueIdx);
    inverse_OCV_func = @(voltage) interp1(ocv_values_sorted, soc_grid_sorted, voltage, 'linear');
    
    % --- 휴지기 찾기 ---
    rest_mask = abs(I_measured) <= 2.0;
    min_rest_time_points = 5990; 
    
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