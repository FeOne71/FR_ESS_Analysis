%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract_Lab_Features_Anchored_SOC.m
% 🛠️ 최종 수정: duration 타입 오류 해결
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% --- 경로 설정 (사용자 제공 정보 기반) ---
inputFolder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
ocvDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated'; 
outputFolder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\features_anchored_soc';

if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end

% --- 설정 변수 ---
Cnom_Ah = 64; % 공칭 용량 사용
cycleTypes = {'0cyc', '200cyc', '400cyc'};
SOCLevels = {'SOC90', 'SOC70', 'SOC50'};

% --- OCV 데이터 로드 ---
ocvMatFile = fullfile(ocvDataPath, 'OCV_integrated.mat');
if ~exist(ocvMatFile, 'file')
    error('OCV_integrated.mat 파일을 찾을 수 없습니다: %s\n', ocvMatFile);
end
load(ocvMatFile, 'OCV_data');
fprintf('OCV_integrated.mat 로드 완료.\n');

% 최종 결과 테이블 초기화
featureTable_Lab = table();

fprintf('랩 데이터 주행부하 특징 추출 시작...\n');

%% 2. 데이터 로드 및 특징 추출 루프
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    inputFileName = fullfile(inputFolder, sprintf('parsedDriveCycle_%s_filtered.mat', cycleType));
    
    if ~exist(inputFileName, 'file'); continue; end
    
    fprintf('\n=== %s 데이터 로드 및 처리 시작 ===\n', cycleType);
    load(inputFileName, ['parsedDriveCycle_', cycleType]);
    
    ocv_struct_name = ['OCV_integrated_', cycleType];
    if ~isfield(OCV_data, ocv_struct_name)
         fprintf('경고: %s 에 대한 통합 OCV 데이터가 없습니다. 스킵합니다.\n', cycleType);
         continue;
    end
    eval(['OCV_data_current = OCV_data.', ocv_struct_name, ';']);
    
    eval(['currentData = parsedDriveCycle_', cycleType, ';']);
    channelFields = fieldnames(currentData);
    
    for chIdx = 1:length(channelFields)
        channelFieldName = channelFields{chIdx};
        channelData = currentData.(channelFieldName);
        
        for socIdx = 1:length(SOCLevels)
            socLevel = SOCLevels{socIdx};
            
            if ~isfield(channelData, socLevel); continue; end
            
            profileFields = fieldnames(channelData.(socLevel));
            
            for pIdx = 1:length(profileFields)
                profileName = profileFields{pIdx};
                rawProfile = channelData.(socLevel).(profileName);
                
                % 데이터 추출 및 Power 계산
                V = rawProfile.V;
                I = rawProfile.I;
                t = rawProfile.t; % Time is duration type here, based on error
                
                % Power 계산 (P = V * I / 1000) [kW]
                P_kW = (V .* I) / 1000;
                
                % 🛠️ 오류 수정: duration 타입을 double(초)로 변환
                dt_sec = seconds(mean(diff(t))); 
                dt_hr = dt_sec / 3600;
                
                % fs_hz 계산 (이제 dt_sec은 double 타입)
                fs_hz = 1 / dt_sec;
                
                % 🛠️ SOC 프로파일 재계산 (OCV 앵커링 로직 사용)
                try
                    [SOC_full, V_ocv_initial, V_ocv_final] = calculate_anchored_soc(V, I, t, dt_sec, OCV_data_current);
                catch ME
                    fprintf('  [오류] %s - %s - %s: SOC 계산 실패. 스킵: %s\n', cycleType, socLevel, profileName, ME.message);
                    continue;
                end
                
                % 12가지 특징 추출 함수 호출
                features = extract_duty_cycle_features(P_kW, SOC_full, dt_hr, fs_hz);
                
                % 결과 테이블에 추가
                newRow = struct2table(features);
                
                % 식별 변수 추가
                newRow.Channel = {channelFieldName};
                newRow.CycleType = {cycleType};
                newRow.SOCLevel = {socLevel};
                newRow.ProfileName = {profileName};
                newRow.V_OCV_Start = V_ocv_initial;
                newRow.V_OCV_End = V_ocv_final;
                
                featureTable_Lab = [featureTable_Lab; newRow];
            end
        end
    end
end

%% 3. 결과 저장
saveFileName = fullfile(outputFolder, 'Lab_DutyCycle_Features_X_Lab.mat');
save(saveFileName, 'featureTable_Lab');
fprintf('\n=== 랩 데이터 특징 추출 완료 ===\n');
fprintf('총 %d개의 주행부하 데이터 포인트가 추출되었습니다.\n', size(featureTable_Lab, 1));
fprintf('결과 파일: %s\n', saveFileName);

%% ========================================================================
% 🌟 서브 함수 1: 12가지 특징 추출 (extract_duty_cycle_features)
% ========================================================================

function [features] = extract_duty_cycle_features(P, SOE, dt_hr, fs_hz)

% I. 시간 영역 지표 (6개)
P_dis = P; P_dis(P_dis > -1e-6) = 0; % P < 0 방전
P_ch  = P; P_ch(P_ch < 1e-6) = 0;   % P > 0 충전

features.E_dis = abs(sum(P_dis)) * dt_hr;
features.E_ch = sum(P_ch) * dt_hr;   

features.P_dis_max = abs(min(P_dis)); 
features.P_ch_max = max(P_ch); 

features.SOE_mean = mean(SOE);
features.SOE_std = std(SOE);

% II. 주파수 영역 지표 (6개) - Hz 단위
[f_max_dis, f_10_dis, f_90_dis] = calculate_psd_metrics(abs(P_dis), fs_hz); 
features.f_max_dis = f_max_dis;
features.f_10_dis = f_10_dis;
features.f_90_dis = f_90_dis;

[f_max_ch, f_10_ch, f_90_ch] = calculate_psd_metrics(P_ch, fs_hz); 
features.f_max_ch = f_max_ch;
features.f_10_ch = f_10_ch;
features.f_90_ch = f_90_ch;

end

%% ========================================================================
% 🌟 서브 함수 2: PSD 특징 추출 (calculate_psd_metrics)
% ========================================================================

function [f_max, f_10, f_90] = calculate_psd_metrics(P_event, fs_hz)

events = {}; 
in_rest = true; 
current_event = [];
P_event(isnan(P_event)) = 0;

for i = 1:length(P_event)
    if P_event(i) > 1e-6 
        current_event = [current_event, P_event(i)];
        in_rest = false;
    elseif ~in_rest && P_event(i) <= 1e-6 
        if length(current_event) > 1 
            centered_profile = [current_event, -current_event];
            events{end+1} = centered_profile;
        end
        current_event = [];
        in_rest = true;
    end
end
if ~isempty(current_event) && length(current_event) > 1
    centered_profile = [current_event, -current_event];
    events{end+1} = centered_profile;
end


if isempty(events)
    f_max = 0; f_10 = 0; f_90 = 0; 
    return;
end

combined_profile = cell2mat(events);

nfft = 2^nextpow2(length(combined_profile));
window = hamming(floor(length(combined_profile)/4)); 
noverlap = floor(length(window) * 0.5);

[Pxx, F] = pwelch(combined_profile, window, noverlap, nfft, fs_hz); 

Pxx = Pxx(2:end);
F = F(2:end);

if isempty(F)
     f_max = 0; f_10 = 0; f_90 = 0;
     return;
end

total_power = sum(Pxx);
cumulative_power = cumsum(Pxx); 

[~, max_idx] = max(Pxx);
f_max = F(max_idx);

idx_10 = find(cumulative_power >= 0.1 * total_power, 1, 'first');
f_10 = F(idx_10);

idx_90 = find(cumulative_power >= 0.9 * total_power, 1, 'first');
f_90 = F(idx_90);

end

%% ========================================================================
% 🌟 서브 함수 3: OCV 앵커링 기반 SOC 계산 (calculate_anchored_soc)
% ========================================================================

function [SOC_full, V_ocv_initial, V_ocv_final] = calculate_anchored_soc(V_measured, I_measured, t_measured, dt_sec, OCV_data)

    % --- OCV 데이터 추출 ---
    soc_grid = OCV_data.SOC_grid;
    ocv_values = OCV_data.V_avg_SOC;
    
    % Inverse OCV function (Voltage → SOC) 생성
    [ocv_values_sorted, uniqueIdx] = unique(ocv_values, 'stable');
    soc_grid_sorted = soc_grid(uniqueIdx);
    inverse_OCV_func = @(voltage) interp1(ocv_values_sorted, soc_grid_sorted, voltage, 'linear');
    
    % --- 휴지기 찾기 ---
    rest_mask = abs(I_measured) <= 3.2; % |I| <= 2A를 휴지기로 정의
    min_rest_time_points = dt_sec * 600; 
    
    rest_periods = [];
    in_rest = false;
    rest_start = 0;
    
    for i = 1:length(rest_mask)
        if rest_mask(i) && ~in_rest
            rest_start = i;
            in_rest = true;
        elseif ~rest_mask(i) && in_rest
            rest_duration_points = (i-1) - rest_start + 1;
            if rest_duration_points >= min_rest_time_points
                rest_periods = [rest_periods; rest_start, i-1];
            end
            in_rest = false;
        end
    end
    
    if in_rest 
        rest_duration_points = length(rest_mask) - rest_start + 1;
        if rest_duration_points >= min_rest_time_points
            rest_periods = [rest_periods; rest_start, length(rest_mask)];
        end
    end
    
    if size(rest_periods, 1) < 2
        error('시작 및 끝 휴지기 구간을 모두 찾을 수 없습니다. SOC 계산 불가.');
    end
    
    % --- SOC 앵커링 ---
    initial_rest_end = rest_periods(1, 2);
    final_rest_end = rest_periods(end, 2);
    
    % SOC1 (초기)
    V_ocv_initial = V_measured(initial_rest_end);
    SOC_initial = inverse_OCV_func(V_ocv_initial);
    
    % SOC2 (최종)
    V_ocv_final = V_measured(final_rest_end);
    SOC_final = inverse_OCV_func(V_ocv_final);
    
    % --- 전류 적산 및 SOC 프로파일 계산 ---
    N = length(I_measured);
    SOC_full = zeros(N, 1);
    
    % dt_sec을 사용하여 누적 전류 적분 (A·s)
    cumulative_current_integral = zeros(N, 1);
    
    for i = (initial_rest_end + 1):final_rest_end
        cumulative_current_integral(i) = cumulative_current_integral(i-1) + I_measured(i) * dt_sec;
    end
    
    % 전체 전류 적산 (A·s)
    total_current_integration = cumulative_current_integral(final_rest_end);
    
    if abs(total_current_integration) < 1e-6
        SOC_full(:) = SOC_initial;
        return;
    end
    
    % SOC(t) = SOC1 + (SOC2-SOC1) * (∫₀ᵗ I*dt) / (∫₀ᵉⁿᵈ I*dt) 적용
    for i = 1:N
        if i <= initial_rest_end
            SOC_full(i) = SOC_initial;
        elseif i >= final_rest_end
            SOC_full(i) = SOC_final;
        else
            SOC_full(i) = SOC_initial + (SOC_final - SOC_initial) * (cumulative_current_integral(i) / total_current_integration);
        end
    end
    
end