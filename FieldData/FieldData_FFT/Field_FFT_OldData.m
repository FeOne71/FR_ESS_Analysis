%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field_FFT_OldData.m
% ESS 운전 주기 특징 추출 및 시각화 (P < 0: 방전, P > 0: 충전 규칙 적용)
% 최종 수정: PSD 플롯을 일자별 Subplot (충전/방전) 형식으로 변경 (단위: Hz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% File D
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_FFT','FFTMetric_Results_OldData_Hz');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters (사용자 정의 변수)
Cnom = 128;
Cnom_cell = 64;                   % Rack nominal Capacity (Ah)
idle_thr = Cnom_cell*0.05;         % Idle threshold [charge, discharge] (A)
Ns = 17*14;    % 238s
Np = 2;        % 2p

%% Load Old Data
yearList = {'2021'};
rackNames_all = {'Rack01'};

%% Process each year
for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    fprintf('Processing year: %s\n', year);
    yearPath = fullfile(dataDir, year);
    
    monthDirs = dir(fullfile(yearPath, '20*'));
    
    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        fprintf('Processing month: %s\n', monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, '*.mat'));
        
        [~, idx] = sort({matFiles.name});
        matFiles = matFiles(idx);
        
        % 월별 테이블 초기화
        monthly_table = table();
        month_name = monthDirs(m).name;
        
        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            fprintf('Loading file: %s\n', matFiles(f).name);
            
            try
                 load(matFilePath);
            catch ME
                 fprintf('Error loading file %s: %s\n', matFiles(f).name, ME.message);
                 continue;
            end
            
            for rack_idx = 1:length(rackNames_all)
                rackName = rackNames_all{rack_idx};
                fprintf('  %s: ', rackName);
                
                if ~isfield(Raw, rackName)
                    fprintf('0 chgPeaks 0 dchPeaks (no data)\n');
                    continue;
                end
                
                rackData = Raw.(rackName);
                t_dt = rackData.Time; 
                soc = rackData.SOCPct;
                P = rackData.DCPower_kW; % P < 0 Discharge, P > 0 Charge

                t = datetime(t_dt);
                
                % =========================================================
                % 💡 특징 추출 및 시각화 로직 💡
                % =========================================================
                
                % 1. 유효 데이터 필터링 및 시간 간격 계산
                if length(t) < 2
                    fprintf('Data too short for analysis.\n');
                    continue;
                end
                
                P_valid = P(~isnan(P));
                t_valid = t(~isnan(P));
                soc_valid = soc(~isnan(P));
                
                if length(t_valid) < 2
                    fprintf('Insufficient valid data points.\n');
                    continue;
                end
                
                dt_sec = mean(seconds(diff(t_valid))); 
                dt_hr = dt_sec / 3600;  
                fs_hz = 1 / dt_sec;         

                % 2. 데이터 분리
                P_dis = abs(P_valid); P_dis(P_valid > -1e-6) = 0; % 방전 (P < 0) - 양수화
                P_ch  = P_valid; P_ch(P_valid < 1e-6) = 0;   % 충전 (P > 0)

                % 3. 12가지 특징 추출 함수 호출 (시간 영역 + 주파수 영역 3가지 지표)
                try
                    % 주파수 영역 지표만 필요하므로 calculate_psd_metrics를 별도로 호출하지 않고,
                    % extract_duty_cycle_features를 사용하여 6가지 지표를 테이블에 저장합니다.
                    features = extract_duty_cycle_features(P_valid, soc_valid, dt_hr, fs_hz);
                catch ME
                    fprintf('Error during feature extraction for %s: %s\n', matFiles(f).name, ME.message);
                    continue;
                end

                % 4. PSD 시각화를 위한 전체 데이터 (Pxx, F, f_max, f_10, f_90) 추출
                [f_max_dis, f_10_dis, f_90_dis, Pxx_dis, F_dis] = calculate_psd_metrics(P_dis, fs_hz);
                [f_max_ch, f_10_ch, f_90_ch, Pxx_ch, F_ch] = calculate_psd_metrics(P_ch, fs_hz);
                
                [~, name, ~] = fileparts(matFiles(f).name);
                DateStr = name(5:12); 
                
                % 5. Subplot으로 PSD 시각화
                plot_combined_psd_analysis(DateStr, Pxx_dis, F_dis, f_max_dis, f_10_dis, f_90_dis, ...
                                           Pxx_ch, F_ch, f_max_ch, f_10_ch, f_90_ch);


                % 6. 결과를 테이블로 변환하고 월별 테이블에 추가
                fNames = fieldnames(features);
                T_row = table('Size', [1, length(fNames)+2], ...
                              'VariableTypes', [{'datetime'}, {'string'}, repmat({'double'}, 1, length(fNames))], ...
                              'VariableNames', [{'Date'}, {'RackName'}, fNames']);
                
                T_row.Date(1) = datetime(DateStr, 'InputFormat', 'yyyyMMdd'); 
                T_row.RackName(1) = rackName;

                for k = 1:length(fNames)
                    T_row.(fNames{k})(1) = features.(fNames{k});
                end
                
                % 🛠️ calculate_psd_metrics에서 Pxx, F를 반환하도록 수정했기 때문에
                % extract_duty_cycle_features에서 반환하는 f_max_dis, f_10_dis, f_90_dis 등의 값을
                % 계산 시 사용한 값으로 덮어쓰기 합니다. (정확성 유지)
                T_row.f_max_dis = f_max_dis; T_row.f_10_dis = f_10_dis; T_row.f_90_dis = f_90_dis;
                T_row.f_max_ch = f_max_ch; T_row.f_10_ch = f_10_ch; T_row.f_90_ch = f_90_ch;


                monthly_table = [monthly_table; T_row];

                fprintf('Metrics calculated successfully.\n');

            end % End of rack loop
        end % End of file loop
        
        % 7. 월별 결과 저장 및 박스 플롯 시각화
        if ~isempty(monthly_table)
            saveFileName = fullfile(saveDir, ['FFTMetrics_', year, '_', month_name, '_Hz.mat']); % 저장 파일명에 Hz 명시
            save(saveFileName, 'monthly_table');
            fprintf('Saved monthly metrics to: %s\n', saveFileName);

            % 박스 플롯 시각화 (변경 없음)
            
            % 시간/SOE 지표 (1~6)
            figure('Name', ['Time Domain Metrics: ', month_name]);
            time_metrics = {'E_dis', 'E_ch', 'P_dis_max', 'P_ch_max', 'SOE_mean', 'SOE_std'};
            boxplot(monthly_table{:, time_metrics}, 'Labels', time_metrics); 
            title(['Time Domain Metrics Distribution for ', month_name]);
            ylabel('Value');
            grid on;

            % 주파수 지표 (7~12)
            figure('Name', ['Frequency Domain Metrics: ', month_name, ' (Hz)']);
            freq_metrics = {'f_max_dis', 'f_10_dis', 'f_90_dis', 'f_max_ch', 'f_10_ch', 'f_90_ch'};
            boxplot(monthly_table{:, freq_metrics}, 'Labels', freq_metrics); 
            title(['Frequency Domain Metrics Distribution for ', month_name, ' (Hz)']);
            ylabel('Frequency (Hz)'); 
            grid on;
        end
    end % End of month loop
end % End of year loop
disp('Processing complete.');

%% ========================================================================
% 🌟 서브 함수 1: 12가지 특징 추출 (extract_duty_cycle_features)
% ========================================================================

function [features] = extract_duty_cycle_features(P, SOE, dt_hr, fs_hz)
% ESS 운전 주기의 12가지 특징을 추출합니다. 

% I. 시간 영역 지표 (6개)
P_dis = P; P_dis(P_dis > -1e-6) = 0; % 방전 (P < 0)
P_ch  = P; P_ch(P_ch < 1e-6) = 0;   % 충전 (P > 0)

features.E_dis = abs(sum(P_dis)) * dt_hr;
features.E_ch = sum(P_ch) * dt_hr; 
features.P_dis_max = abs(min(P_dis));
features.P_ch_max = max(P_ch);
features.SOE_mean = mean(SOE);
features.SOE_std = std(SOE);

% II. 주파수 영역 지표 (6개)
% PSD 커브는 필요 없으므로 3가지 지표만 받습니다.
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
% 🌟 서브 함수 2: PSD 특징 추출 및 커브 반환 (calculate_psd_metrics)
% 🛠️ Pxx와 F를 추가로 반환하도록 수정
% ========================================================================

function [f_max, f_10, f_90, Pxx, F] = calculate_psd_metrics(P_event, fs_hz)
% PSD 특징과 플롯 데이터를 추출합니다.

% 1. 개별 이벤트 식별 및 평균 중심화 
events = {}; 
in_event = false;
current_event = [];
P_event(isnan(P_event)) = 0;

for i = 1:length(P_event)
    if P_event(i) > 1e-6 
        current_event = [current_event, P_event(i)];
        in_event = true;
    elseif in_event 
        if ~isempty(current_event)
            centered_profile = [current_event, -current_event];
            events{end+1} = centered_profile;
        end
        current_event = [];
        in_event = false;
    end
end
if in_event && ~isempty(current_event)
    centered_profile = [current_event, -current_event];
    events{end+1} = centered_profile;
end

if isempty(events)
    f_max = 0; f_10 = 0; f_90 = 0; Pxx = 0; F = 0; % 데이터가 없으면 0 반환
    return;
end

% 2. 모든 이벤트 프로파일을 하나로 연결
combined_profile = cell2mat(events);

% 3. Welch Method를 사용하여 PSD 추정 (pwelch)
nfft = 2^nextpow2(length(combined_profile));
window = hamming(floor(length(combined_profile)/4)); 
noverlap = floor(length(window) * 0.5);

[Pxx, F] = pwelch(combined_profile, window, noverlap, nfft, fs_hz); 

% DC 성분 (F=0)을 제외하고 분석 시작
Pxx = Pxx(2:end);
F = F(2:end);

if isempty(F)
     f_max = 0; f_10 = 0; f_90 = 0; Pxx = 0; F = 0;
     return;
end

% 4. PSD 지표 추출
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
% 🌟 서브 함수 3: Subplot에 PSD 그래프 그리기 (plot_subplot_psd)
% ========================================================================

function plot_subplot_psd(F, Pxx, f_max, f_10, f_90, title_str)
% 하나의 subplot에 PSD 플롯을 그립니다.

    if isscalar(Pxx) && Pxx == 0
        text(0.5, 0.5, 'No significant power events', 'HorizontalAlignment', 'center', 'Color', 'r');
        title(title_str);
        return;
    end

    loglog(F, Pxx, 'b', 'LineWidth', 1.5);
    grid on;
    title(title_str);
    xlabel('Frequency (Hz)', 'Interpreter', 'none'); 
    ylabel('Power Spectral Density (Pxx)', 'Interpreter', 'none'); 
    
    % f_max, f_10, f_90 선 표시
    hold on;
    
    y_min = min(Pxx) * 0.1;
    y_max = max(Pxx) * 10;
    
    % f_max 표시
    loglog([f_max, f_max], [y_min, y_max], 'r--', 'LineWidth', 1);
    text(f_max * 1.1, max(Pxx) * 0.5, ['fmax: ', num2str(f_max, '%.4f')], 'Color', 'r', 'FontSize', 10, 'Interpreter', 'none');
    
    % f_10 표시
    loglog([f_10, f_10], [y_min, y_max], 'g:', 'LineWidth', 1);
    text(f_10 * 1.1, max(Pxx) * 0.3, ['f10%: ', num2str(f_10, '%.4f')], 'Color', 'g', 'FontSize', 10, 'Interpreter', 'none');
    
    % f_90 표시
    loglog([f_90, f_90], [y_min, y_max], 'm:', 'LineWidth', 1);
    text(f_90 * 1.1, max(Pxx) * 0.7, ['f90%: ', num2str(f_90, '%.4f')], 'Color', 'm', 'FontSize', 10, 'Interpreter', 'none');
    
    hold off;
    
end

%% ========================================================================
% 🌟 서브 함수 4: 충전/방전 PSD 통합 플롯 함수 (plot_combined_psd_analysis)
% ========================================================================

function plot_combined_psd_analysis(DateStr, Pxx_dis, F_dis, f_max_dis, f_10_dis, f_90_dis, ...
                                    Pxx_ch, F_ch, f_max_ch, f_10_ch, f_90_ch)
    
    figure('Name', ['PSD Analysis - ', DateStr], 'Position', [100 100 800 600]);
    sgtitle(['Power Spectral Density Analysis (', DateStr, ') - Unit: Hz']); % 전체 제목
    
    % Subplot 1: Discharge
    subplot(2, 1, 1);
    plot_subplot_psd(F_dis, Pxx_dis, f_max_dis, f_10_dis, f_90_dis, 'Discharge PSD');
    
    % Subplot 2: Charge
    subplot(2, 1, 2);
    plot_subplot_psd(F_ch, Pxx_ch, f_max_ch, f_10_ch, f_90_ch, 'Charge PSD');

end