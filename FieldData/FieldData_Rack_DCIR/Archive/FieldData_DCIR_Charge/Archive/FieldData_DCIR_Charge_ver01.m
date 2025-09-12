%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filed Data DCIR Charge ver01.m
% ESS 충전 이벤트 및 저항 분석 (BSC_Charge 기반)
% Step 1: Folder traversal
% Step 2: BSC_Charge 기반 이벤트 검출
% Step 3: DCIR 계산
% Step 4: Visualization
% 
% Ver02 Updates:
% - Added Power-based filtering (max_P, min_P, max_P_std)
% - Extended Condition #3: Current + Power magnitude check
% - Extended Condition #4: Current + Power stability check
% - Power data stored in P_seq and P_seg_ridIdle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average SOC(%)	Highest SOC(%)	Lowest SOC(%)	Highest SOC Pos.(R)	Lowest SOC Pos.(R)	Average SOH(%)	Average C.V. Sum(V)	DC Current(A)	Highest DC Current(A)	Lowest DC Current(A)	Highest DC Current Pos.(R)	Lowest DC Current Pos.(R)	DC Chg. Current Limit(A)	DC Dchg. Current Limit(A)	DC Power(kW)

clc; clear; close all;

%% Directory
% dataDir  = 'D:\JCW\KENTECH\Projects\KEPCO\ESS_Data_Preprocessing';
dataDir  = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\raw2mat_ver04';
yearList = {'2023'}; % ,'2023','2023'};
saveDir  = fullfile('G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\FieldData_DCIR_Charge');

if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

%% Variables
C_nom        = 1024;          % Ah
min_charge_duration = 30;     % [s] - Charging duration
max_I        = 0.20 * C_nom;  % Max charge current               [256.0 (A)]
min_I        = 0.10 * C_nom;  % min charge current               [102.4 (A)]
std_I        = 0.01 * C_nom;  % Allowed std for charging current [51.2  (A)]
max_P        = 250;           % Max charge power [kW]
min_P        = 150;           % Min charge power [kW]
max_P_std    = 50;            % Max power standard deviation [kW]
dt_sec       = 1;             % Sampling time [sec]

%% Step 1: Raw Files Traversal
fprintf('Step 1: Detecting charging events using BSC_Charge...\n');
eventStruct = struct();
eventCount = 0;

total_events    = 0;
filtered_events = 0;

for y = 1:length(yearList)
    year = yearList{y};
    yearPath = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));

    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, 'Raw_*.mat'));

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            [~, name, ~] = fileparts(matFiles(f).name);
            dateStr = extractAfter(name, 'Raw_');
            fieldName = ['date_' dateStr];
            
            load(matFilePath); 
            t          = Raw.Time;
            I          = Raw.DCCurrent_A;
            V          = Raw.Total_AverageCVSum_V;
            soc        = Raw.Total_AverageSOC;
            T_batt     = Raw.Total_AverageMT_degC;
            bsc_charge = Raw.Charge;  % string array
            status     = Raw.Status;
            dc_power   = Raw.DCPower_kW;
            
%% Step 2: Charging Event

            % 1) Define "Idle" and "Charging" status
            is_idle = strcmp(bsc_charge, 'Idle');
            is_charging = strcmp(bsc_charge, 'Charging');
            
            % Detect "Idle" -> "Charging"
            idle_to_charge = find(is_idle(1:end-1) & is_charging(2:end));
            total_events = total_events + length(idle_to_charge);
            fprintf('File: %s - Found %d potential events\n', matFiles(f).name, length(idle_to_charge));
            
            fprintf('  Total Idle points: %d\n', sum(is_idle));
            fprintf('  Total Charging points: %d\n', sum(is_charging));
            fprintf('  Idle→Charging transitions: %d\n', length(idle_to_charge));

            step2_count = 0;
            step3_count = 0;
            step4_count = 0;
            step5_count = 0;
            step6_count = 0;
            
            for i = 1:length(idle_to_charge)
                idx1 = idle_to_charge(i);  % "Idle"의 마지막 시점
                start_charge_idx = idx1 + 1;  % "Charging" 시작점
                
                % 2) "Charging" 구간 끝점 찾기
                charge_end_idx = start_charge_idx;
                while charge_end_idx <= length(bsc_charge) && strcmp(bsc_charge(charge_end_idx), 'Charging')
                    charge_end_idx = charge_end_idx + 1;
                end
                charge_end_idx = charge_end_idx - 1;  % "Charging"의 마지막 시점
                
                % 이벤트 전체 구간 (Idle 시작점부터 Charging 끝점까지)
                start_idx = idx1;
                end_idx = charge_end_idx;
                
                step2_count = step2_count + 1;
                
                % Condition #1 "Charging" 구간 지속 시간 확인 (5초 이상)
                charge_duration = charge_end_idx - start_charge_idx + 1;
                if charge_duration < min_charge_duration % 5sec
                    fprintf('  Event %d filtered: Charge duration too short (%ds)\n', i, charge_duration);
                    continue;
                end
                
                step3_count = step3_count + 1;
                
                t_seg   = t(start_idx:end_idx);
                I_seg   = I(start_idx:end_idx);
                V_seg   = V(start_idx:end_idx);
                soc_seg = soc(start_idx:end_idx);
                T_seg   = T_batt(start_idx:end_idx);
                bsc_seg = bsc_charge(start_idx:end_idx);
                P_seg   = dc_power(start_idx:end_idx);
                
                duration_sec = seconds(t_seg(end) - t_seg(1));
                deltaI = max(I_seg) - min(I_seg);
                I_max = max(abs(I_seg));
                
                fprintf('Event %d: Duration=%.1fs, Charge_Duration=%ds, deltaI=%.2fA, I_max=%.2fA\n', i, duration_sec, charge_duration, deltaI, I_max);
                
                % idx2는 "Charging"의 마지막 시점
                idx2 = charge_end_idx - start_idx + 1;  % 세그먼트 내에서의 상대적 위치
                
                % Charging 구간의 시작점 (세그먼트 내 상대위치)
                charge_start_in_segment = start_charge_idx - start_idx + 1;

                % 5) 전류 및 출력 크기 (확장 조건)
                max_current = max(I_seg(3:idx2));
                min_current = min(I_seg(3:idx2));
                max_power   = max(P_seg(3:idx2));
                min_power   = min(P_seg(3:idx2));

                fprintf('  Event %d: max_current=%.1fA (limit=%.1fA), max_power=%.1fkW (limit=%.1fkW)\n', ...
                    i, max_current, max_I, max_power, max_P);
                fprintf('  Event %d: min_current=%.1fA (limit=%.1fA), min_power=%.1fkW (limit=%.1fkW)\n', ...
                    i, min(I_seg(5:idx2)), min_I, min_power, min_P);

                % 전류 크기 (Charging 구간만 체크)
                current_ok = (max_current <= max_I) && (min(I_seg(5:idx2)) >= min_I);
                
                % Condition #3b: 출력 크기 (Charging 구간만 체크)
                power_ok = (max_power <= max_P) && (min_power >= min_P);
                
                % Condition #2 충전 전류 min Max 
                % if max_current > max_I || min_current < min_I
                if any(I_seg(3:idx2) > max_I) || any(I_seg(3:idx2) < min_I)
                    fprintf('    → FILTERED: Max current too high (%.1fA > %.1fA)\n', max(I_seg(3:idx2)), max_I);
                    fprintf('    → FILTERED: min current too low  (%.1fA < %.1fA)\n', min(I_seg(3:idx2)), min_I);
                     continue;
                 end
                 
                % if max_power > max_P || min_power < min_P
                %     fprintf('    → FILTERED: Max power too high (%.1fA > %.1fA)\n', max_power, max_P);
                %     fprintf('    → FILTERED: min power too low  (%.1fA < %.1fA)\n', min_power, min_P);
                %     continue;
                % end

                step4_count = step4_count + 1;                               

                % Condition #3: 전류 및 출력 안정성 (Charging 구간만 체크)
                current_std = std(I_seg(3:idx2));
                power_std   = std(P_seg(3:idx2));
                
                current_stable = current_std <= std_I;
                power_stable   = power_std <= max_P_std;
                
                fprintf('  Event %d: current_std=%.1fA (limit=%.1fA), power_std=%.1fkW (limit=%.1fkW)\n', ...
                    i, current_std, C_nom*0.1, power_std, max_P_std);

                if ~current_stable %|| ~power_stable
                    fprintf('    → FILTERED: Current instability (std=%.1fA > %.1fA)\n', current_std, C_nom*0.1);
                continue;    
                end
                step5_count = step5_count + 1;

                % if ~power_stable
                %     fprintf('    → FILTERED: Current instability (std=%.1fA > %.1fA)\n', current_std, C_nom*0.1);
                %     continue;
                % end
                % 
                % step6_count = step6_count + 1;


%% Step 3: DCIR Calculation

                % Save Events as Struct
                eventCount = eventCount + 1;
                evtName = sprintf('event%d', eventCount);
                eventStruct.(evtName).start_idx = start_idx;
                eventStruct.(evtName).end_idx = end_idx;
                eventStruct.(evtName).idx1 = idx1;  % Idle의 마지막 시점
                eventStruct.(evtName).idx2 = idx2;  % Charging의 마지막 시점 (세그먼트 내 상대위치)
                eventStruct.(evtName).charge_duration = charge_duration;

                eventStruct.(evtName).t_seq   = t_seg;
                eventStruct.(evtName).I_seq   = I_seg;
                eventStruct.(evtName).V_seq   = V_seg;
                eventStruct.(evtName).T_seq   = T_seg;
                eventStruct.(evtName).soc_seq = soc_seg;
                eventStruct.(evtName).bsc_seq = bsc_seg;
                eventStruct.(evtName).P_seq   = P_seg;

                eventStruct.(evtName).t_seg_ridIdle = t_seg(1:idx2);
                eventStruct.(evtName).I_seg_ridIdle = I_seg(1:idx2);
                eventStruct.(evtName).V_seg_ridIdle = V_seg(1:idx2);
                eventStruct.(evtName).T_seg_ridIdle = T_seg(1:idx2);
                eventStruct.(evtName).soc_seq_ridIdle = soc_seg(1:idx2);
                eventStruct.(evtName).P_seg_ridIdle = P_seg(1:idx2);

                % DCIR 계산: 시간별
                dt_list = [1, 3, 5, 10, 30, 50];
                for d = 1:length(dt_list)
                    dt_sec = dt_list(d);
                    idx_dt = 1 + dt_sec;

                    if idx_dt <= length(I_seg)
                        V1 = V_seg(1);
                        V2 = V_seg(idx_dt);
                        I1 = I_seg(1);
                        I2 = I_seg(idx_dt);
                        dV = V2 - V1;
                        dI = I2 - I1;

                        if (dI) > 0 && dV > 0
                            dcir_val = (dV / dI) * 1000; % [mOhm]
                        else
                            dcir_val = NaN;
                        end
                    else
                        V1 = NaN; V2 = NaN; I1 = NaN; I2 = NaN;
                        dcir_val = NaN;
                    end

                    fieldName = sprintf('DCIR_%ds', dt_sec);

                    eventStruct.(evtName).(fieldName).val = dcir_val;
                    eventStruct.(evtName).(fieldName).V1 = V1;
                    eventStruct.(evtName).(fieldName).V2 = V2;
                    eventStruct.(evtName).(fieldName).I1 = I1;
                    eventStruct.(evtName).(fieldName).I2 = I2;
                    eventStruct.(evtName).(fieldName).dV = dV;
                    eventStruct.(evtName).(fieldName).dI = dI;
                end              

                % Charging 끝점 기준 DCIR 계산
                if idx2 <= length(I_seg)
                    V1_r = V_seg(1);
                    V2_r = V_seg(idx2);
                    I1_r = I_seg(1);
                    I2_r = I_seg(idx2);
                    dV_ramp = V2_r - V1_r;
                    dI_ramp = I2_r - I1_r;

                    if (dV_ramp) > 0 && (dI_ramp) > 0
                        dcir_ramp = (dV_ramp / dI_ramp) * 1000; % [mOhm]
                    else
                        dcir_ramp = NaN;
                    end
                else
                    V1_r = NaN; V2_r = NaN; I1_r = NaN; I2_r = NaN;
                    dcir_ramp = NaN;
                end

                eventStruct.(evtName).DCIR_ramp.val = dcir_ramp;
                eventStruct.(evtName).DCIR_ramp.V1  = V1_r;
                eventStruct.(evtName).DCIR_ramp.V2  = V2_r;
                eventStruct.(evtName).DCIR_ramp.I1  = I1_r;
                eventStruct.(evtName).DCIR_ramp.I2  = I2_r;
                eventStruct.(evtName).DCIR_ramp.dV_ramp  = dV_ramp;
                eventStruct.(evtName).DCIR_ramp.dI_ramp  = dI_ramp;

                filtered_events = filtered_events + 1;
            end
            
            % 파일별 필터링 결과 요약
            fprintf('  Step summary:\n');
            fprintf('    After charging end detection: %d\n', step2_count);
            fprintf('    After duration filter (≥%ds): %d\n', min_charge_duration, step3_count);
            fprintf('    After current magnitude check: %d\n', step4_count);
            fprintf('    After current stability check: %d\n', step5_count);
            % fprintf('    After power stability check: %d\n', step6_count);
            fprintf('    Final valid events from this file: %d\n', step5_count);
        end
    end
end

save(fullfile(saveDir, 'all_chg_events_2023.mat'), 'eventStruct');
fprintf('Step 1 complete. Events saved.\n');
fprintf('Total potential events detected: %d\n', total_events);
fprintf('Events that passed all filters: %d\n', filtered_events);
fprintf('Final stored events: %d\n', eventCount);
fprintf('\n');

%% Step 4: Visualization

%% fig1: Raw Voltage vs. Time
eventNames = fieldnames(eventStruct); 

fprintf('Creating Fig1 with %d events:\n', length(eventNames));

figure; hold on;
for i = 1:length(eventNames)
    evt = eventStruct.(eventNames{i});
    
    if ~isfield(evt, 'V_seq') || isempty(evt.V_seq)
        fprintf('Event %s: Skipped (empty V_seq)\n', eventNames{i});
        continue;
    end

    t_evt = evt.t_seq;
    V_evt = evt.V_seq;
    t_rel = seconds(t_evt - t_evt(1));
    
    fprintf('Event %s: %d data points, duration=%.1fs\n', ...
        eventNames{i}, length(V_evt), t_rel(end));

    plot(t_rel, V_evt, '-', 'LineWidth', 1.2);
end
xlabel('Time [s]');
ylabel('Voltage [V]');
ylim([870 920]);
title('Fig 1b. Voltage vs. Time');
grid on;

saveas(gcf, fullfile(saveDir, 'fig1b_Voltage_vs_Time_2023.fig'));

% Current vs. Time
figure; hold on;
for i = 1:length(eventNames)
    evt = eventStruct.(eventNames{i});
    
    if ~isfield(evt, 'V_seq') || isempty(evt.V_seq)
        fprintf('Event %s: Skipped (empty V_seq)\n', eventNames{i});
        continue;
    end

    t_evt = evt.t_seq;
    V_evt = evt.V_seq;
    I_evt = evt.I_seq;
    t_rel = seconds(t_evt - t_evt(1));
    
    fprintf('Event %s: %d data points, duration=%.1fs\n', ...
        eventNames{i}, length(V_evt), t_rel(end));

    plot(t_rel, I_evt, '-', 'LineWidth', 1.2);
end
xlabel('Time [s]');
ylabel('Current [A]');
title(sprintf('Fig 1a. Current vs. Time (%s)', yearList{1}));
grid on;

saveas(gcf, fullfile(saveDir, 'fig1a_Current_vs_Time_2023.fig'));

%% DCIR 플롯
dcir_fields = {'DCIR_1s', 'DCIR_3s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s', 'DCIR_50s', 'DCIR_ramp'};
dcir_labels = {'1s', '3s', '5s', '10s', '30s', '50s', 'End'};

eventNames = fieldnames(eventStruct);
nEvents = length(eventNames);

figure('Name', 'DCIR per Event', 'Position', [100, 100, 1200, 900]);

for k = 1:length(dcir_fields)
    dcir_vals = NaN(nEvents,1);
    
    for i = 1:nEvents
        evt = eventStruct.(eventNames{i});
        if isfield(evt, dcir_fields{k}) && isfield(evt.(dcir_fields{k}), 'val')
            dcir_vals(i) = evt.(dcir_fields{k}).val;
        end
    end
    
    subplot(4,2,k); hold on;
    valid_idx = ~isnan(dcir_vals);
    x = find(valid_idx);
    y = dcir_vals(valid_idx);

    % o- plot
    plot(x, y, 'o-', 'Color', '#0072BD', 'LineWidth', 1.5, 'MarkerSize', 5);

    % Average
    dcir_mean = mean(y);
    yline(dcir_mean, '--', 'Color', '#EE4C97', 'LineWidth', 1.5);
    annotation_str = sprintf('Mean = %.2f mΩ', dcir_mean);
    text(0.97, 0.95, annotation_str, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'Color', '#EE4C97', 'FontWeight', 'bold');

    title(sprintf('DCIR at %s', dcir_labels{k}));
    xlabel('Event Index');
    ylabel('[mΩ]');
    grid on;
end

sgtitle('Charging Event DCIR (2023.06)');
saveas(gcf, fullfile(saveDir, 'fig_2023_ver04_DCIR_subplots.fig'));

%% 표준편차 기반 이상치 제거 플롯 추가
fprintf('\n=== 표준편차 기반 이상치 제거 분석 ===\n');

% 이상치 제거 기준 설정
sigma_threshold = 1;  % 2-sigma 기준
dcir_fields = {'DCIR_1s', 'DCIR_3s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s', 'DCIR_50s', 'DCIR_ramp'};
dcir_labels = {'1s', '3s', '5s', '10s', '30s', '50s', 'End'};

% 전체 통계 저장용
outlier_stats = struct();

figure('Name', sprintf('DCIR Outlier Detection (%dσ)', sigma_threshold), 'Position', [50, 50, 1600, 1000]);

for k = 1:length(dcir_fields)
    dcir_vals = NaN(nEvents,1);
    
    % 데이터 수집
    for i = 1:nEvents
        evt = eventStruct.(eventNames{i});
        if isfield(evt, dcir_fields{k}) && isfield(evt.(dcir_fields{k}), 'val')
            dcir_vals(i) = evt.(dcir_fields{k}).val;
        end
    end
    
    % 유효한 데이터와 NaN 데이터 분리
    valid_idx = ~isnan(dcir_vals);
    valid_data = dcir_vals(valid_idx);
    event_indices = find(valid_idx);
    
    % NaN 값이 있는 위치
    nan_idx = find(isnan(dcir_vals));
    
    if length(valid_data) < 3
        fprintf('Warning: %s has insufficient data (%d points)\n', dcir_fields{k}, length(valid_data));
        continue;
    end
    
    % 통계 계산
    data_mean = mean(valid_data);
    data_std = std(valid_data);
    
    % 이상치 검출 (σ 기준)
    outlier_mask = abs(valid_data - data_mean) > (sigma_threshold * data_std);
    outliers = valid_data(outlier_mask);
    normal_data = valid_data(~outlier_mask);
    
    % 이상치 제거 후 통계
    clean_mean = mean(normal_data);
    clean_std = std(normal_data);
    
    % 통계 저장
    outlier_stats.(dcir_fields{k}).original_count = length(valid_data);
    outlier_stats.(dcir_fields{k}).outlier_count = length(outliers);
    outlier_stats.(dcir_fields{k}).outlier_percentage = (length(outliers) / length(valid_data)) * 100;
    outlier_stats.(dcir_fields{k}).original_mean = data_mean;
    outlier_stats.(dcir_fields{k}).original_std = data_std;
    outlier_stats.(dcir_fields{k}).clean_mean = clean_mean;
    outlier_stats.(dcir_fields{k}).clean_std = clean_std;
    
    % 플롯 생성 - 더 큰 subplot 크기로 조정
    subplot(4,2,k); hold on;
    
    % % subplot 간격 조정 (더 큰 subplot 크기)
    % pos = get(gca, 'Position');
    % pos(3) = pos(3) * 1.1;  % 너비 10% 증가
    % pos(4) = pos(4) * 1.1;  % 높이 10% 증가
    % set(gca, 'Position', pos);
    % 
    
    % 원본 데이터 (회색 점) - 실제 이벤트 인덱스 사용
    plot(event_indices, valid_data, 'o', 'Color', [0.7 0.7 0.7], 'MarkerSize', 4, 'DisplayName', 'Original Data');
    
    % 정상 데이터 (파란색 점)
    normal_event_indices = event_indices(~outlier_mask);
    plot(normal_event_indices, normal_data, 'o-', 'Color', '#0072BD', 'MarkerSize', 5, 'DisplayName', 'Normal Data');
    
    % 이상치 (빨간색 X)
    if ~isempty(outliers)
        outlier_event_ind%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field Data DCIR Charge Current Clustering
% ESS 충전 이벤트 및 저항 분석 (전류별 군집화)
% Step 1: Folder traversal
% Step 2: BSC_Charge 기반 이벤트 검출
% Step 3: DCIR 계산 (시간별 + 차이값)
% Step 4: 전류별 군집화 및 시각화
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir  = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\raw2mat_ver04';
yearList = {'2021'}; 
saveDir  = fullfile('G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\FieldData_DCIR_Charge\Results');

if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

%% Variables
C_nom        = 1024;          % Ah
min_charge_duration = 10;      % [s] - Charging duration (5초 이상)
max_P_std    = 10;            % Max power standard deviation [kW]
max_I_std    = C_nom *0.01;
dt_sec       = 1;             % Sampling time [sec]

% 충전 전류 구간별 분류 설정
current_bins = [50, 100, 150, 200, 250, 300];  % [A] - 충전 전류 구간 경계
current_labels = {'range_50_100A', 'range_100_150A', 'range_150_200A', 'range_200_250A', 'range_250_300A'};  % 구간 라벨

%% Step 1: Raw Files Traversal
fprintf('=====Step 1: Loading Field Data=====\n');
eventStruct = struct();
eventCount = 0;

total_events    = 0;
filtered_events = 0;

for y = 1:length(yearList)
    year = yearList{y};
    yearPath = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));

    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, 'Raw_*.mat'));

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            [~, name, ~] = fileparts(matFiles(f).name);
            dateStr = extractAfter(name, 'Raw_');
            fieldName = ['date_' dateStr];
            
            load(matFilePath); 
            t          = Raw.Time;
            I          = Raw.DCCurrent_A;
            V          = Raw.Total_AverageCVSum_V;
            soc        = Raw.Total_AverageSOC;
            T_batt     = Raw.Total_AverageMT_degC;
            bsc_charge = Raw.Charge;  % string array
            status     = Raw.Status;
            dc_power   = Raw.DCPower_kW;
            
%% Step 2: Searching Idle to Charging Transition
            fprintf('=====Step 2: Searching Idle to Charging Transition =====\n');
            % 1) Define "Idle" and "Charging" status
            is_idle = strcmp(bsc_charge, 'Idle');
            is_charging = strcmp(bsc_charge, 'Charging');
            
            % Detect "Idle" -> "Charging"
            idle_to_charge = find(is_idle(1:end-1) & is_charging(2:end));
            total_events = total_events + length(idle_to_charge);
            fprintf('File: %s - Found %d potential events\n', matFiles(f).name, length(idle_to_charge));
            
            fprintf('  Total Idle points: %d\n', sum(is_idle));
            fprintf('  Total Charging points: %d\n', sum(is_charging));
            fprintf('  Idle→Charging transitions: %d\n', length(idle_to_charge));

            step2_count = 0;
            step3_count = 0;
            step4_count = 0;

%% Step 3: Events Filtering              
            fprintf('=====Step 3: Filtering Events =====\n');
            for i = 1:length(idle_to_charge)
                idx1 = idle_to_charge(i);  % "Idle"의 마지막 시점
                start_charge_idx = idx1 + 1;  % "Charging" 시작점
                
                % 2) "Charging" 구간 끝점 찾기
                charge_end_idx = start_charge_idx;
                while charge_end_idx <= length(bsc_charge) && strcmp(bsc_charge(charge_end_idx), 'Charging')
                    charge_end_idx = charge_end_idx + 1;
                end
                charge_end_idx = charge_end_idx - 1;  % "Charging"의 마지막 시점
                
                % 이벤트 전체 구간 (Idle 시작점부터 Charging 끝점까지)
                start_idx = idx1;
                end_idx = charge_end_idx;
                
                step2_count = step2_count + 1;
                
                % Condition #1 "Charging" 구간 지속 시간 확인 (5초 이상)
                charge_duration = charge_end_idx - start_charge_idx + 1;
                if charge_duration < min_charge_duration
                    fprintf('  Event %d filtered: Charge duration too short (%ds)\n', i, charge_duration);
                    continue;
                end
                
                step3_count = step3_count + 1;
                
                t_seg   = t(start_idx:end_idx);
                I_seg   = I(start_idx:end_idx);
                V_seg   = V(start_idx:end_idx);
                soc_seg = soc(start_idx:end_idx);
                T_seg   = T_batt(start_idx:end_idx);
                bsc_seg = bsc_charge(start_idx:end_idx);
                P_seg   = dc_power(start_idx:end_idx);
                
                duration_sec = seconds(t_seg(end) - t_seg(1));
                deltaI = max(I_seg) - min(I_seg);
                I_max = max(abs(I_seg));
                
                fprintf('Event %d: Duration=%.1fs, Charge_Duration=%ds, deltaI=%.2fA, I_max=%.2fA\n', i, duration_sec, charge_duration, deltaI, I_max);
                
                % idx2는 "Charging"의 마지막 시점
                idx2 = charge_end_idx - start_idx + 1;  % 세그먼트 내에서의 상대적 위치
                
                % Charging 구간의 시작점 (세그먼트 내 상대위치)
                charge_start_in_segment = start_charge_idx - start_idx + 1;

                % Condition #3: 출력 안정성 (Charging 구간만 체크)
                power_std = std(P_seg(3:idx2));
                current_std = std(I_seg(3:idx2));
                
                fprintf('  Event %d: power_std=%.1fkW (limit=%.1fkW)\n', i, power_std, max_P_std);

                if power_std >= max_P_std || current_std >= max_I_std
                    fprintf('    → FILTERED: Power & Current instability (std=%.1fkW > %.1fkW) | (std=%.1fA > %.1fA)\n', power_std, max_P_std, current_std, max_I_std);
                    continue;    
                end
                step4_count = step4_count + 1;

%% Step 4: DCIR Calculation
                fprintf('=====Step 4: DCIR Calculation =====\n');
                % Save Events as Struct
                eventCount = eventCount + 1;
                evtName = sprintf('event%d', eventCount);
                eventStruct.(evtName).start_idx = start_idx;
                eventStruct.(evtName).end_idx = end_idx;
                eventStruct.(evtName).idx1 = idx1;  % Idle의 마지막 시점
                eventStruct.(evtName).idx2 = idx2;  % Charging의 마지막 시점 (세그먼트 내 상대위치)
                eventStruct.(evtName).charge_duration = charge_duration;

                eventStruct.(evtName).t_seq   = t_seg;
                eventStruct.(evtName).I_seq   = I_seg;
                eventStruct.(evtName).V_seq   = V_seg;
                eventStruct.(evtName).T_seq   = T_seg;
                eventStruct.(evtName).soc_seq = soc_seg;
                eventStruct.(evtName).bsc_seq = bsc_seg;
                eventStruct.(evtName).P_seq   = P_seg;

                eventStruct.(evtName).t_seg_ridIdle = t_seg(1:idx2);
                eventStruct.(evtName).I_seg_ridIdle = I_seg(1:idx2);
                eventStruct.(evtName).V_seg_ridIdle = V_seg(1:idx2);
                eventStruct.(evtName).T_seg_ridIdle = T_seg(1:idx2);
                eventStruct.(evtName).soc_seq_ridIdle = soc_seg(1:idx2);
                eventStruct.(evtName).P_seg_ridIdle = P_seg(1:idx2);

                % DCIR 계산: 시간별
                dt_list = [1, 3, 5, 10, 30, 50];
                for d = 1:length(dt_list)
                    dt_sec = dt_list(d);
                    idx_dt = 1 + dt_sec;

                    if idx_dt <= length(I_seg)
                        V1 = V_seg(1);
                        V2 = V_seg(idx_dt);
                        I1 = I_seg(1);
                        I2 = I_seg(idx_dt);
                        dV = V2 - V1;
                        dI = I2 - I1;

                        if (dI) > 0 && dV > 0
                            dcir_val = (dV / dI) * 1000; % [mOhm]
                        else
                            dcir_val = NaN;
                        end
                    else
                        V1 = NaN; V2 = NaN; I1 = NaN; I2 = NaN;
                        dcir_val = NaN;
                    end

                    fieldName = sprintf('DCIR_%ds', dt_sec);

                    eventStruct.(evtName).(fieldName).val = dcir_val;
                    eventStruct.(evtName).(fieldName).V1 = V1;
                    eventStruct.(evtName).(fieldName).V2 = V2;
                    eventStruct.(evtName).(fieldName).I1 = I1;
                    eventStruct.(evtName).(fieldName).I2 = I2;
                    eventStruct.(evtName).(fieldName).dV = dV;
                    eventStruct.(evtName).(fieldName).dI = dI;
                end              

                % Charging 끝점 기준 DCIR 계산
                if idx2 <= length(I_seg)
                    V1_r = V_seg(1);
                    V2_r = V_seg(idx2);
                    I1_r = I_seg(1);
                    I2_r = I_seg(idx2);
                    dV_ramp = V2_r - V1_r;
                    dI_ramp = I2_r - I1_r;

                    if (dV_ramp) > 0 && (dI_ramp) > 0
                        dcir_ramp = (dV_ramp / dI_ramp) * 1000; % [mOhm]
                    else
                        dcir_ramp = NaN;
                    end
                else
                    V1_r = NaN; V2_r = NaN; I1_r = NaN; I2_r = NaN;
                    dcir_ramp = NaN;
                end

                eventStruct.(evtName).DCIR_ramp.val = dcir_ramp;
                eventStruct.(evtName).DCIR_ramp.V1  = V1_r;
                eventStruct.(evtName).DCIR_ramp.V2  = V2_r;
                eventStruct.(evtName).DCIR_ramp.I1  = I1_r;
                eventStruct.(evtName).DCIR_ramp.I2  = I2_r;
                eventStruct.(evtName).DCIR_ramp.dV_ramp  = dV_ramp;
                eventStruct.(evtName).DCIR_ramp.dI_ramp  = dI_ramp;

%% Step 5: DCIR Difference Calculation
                fprintf('=====Step 5: DCIR Difference Calculation Ts-1s =====\n');
                % 추가 DCIR 차이값 계산
                if isfield(eventStruct.(evtName), 'DCIR_5s') && isfield(eventStruct.(evtName), 'DCIR_1s')
                    if ~isnan(eventStruct.(evtName).DCIR_5s.val) && ~isnan(eventStruct.(evtName).DCIR_1s.val)
                        eventStruct.(evtName).DCIR_diff_5s_1s = eventStruct.(evtName).DCIR_5s.val - eventStruct.(evtName).DCIR_1s.val;
                    else
                        eventStruct.(evtName).DCIR_diff_5s_1s = NaN;
                    end
                end
                
                if isfield(eventStruct.(evtName), 'DCIR_10s') && isfield(eventStruct.(evtName), 'DCIR_1s')
                    if ~isnan(eventStruct.(evtName).DCIR_10s.val) && ~isnan(eventStruct.(evtName).DCIR_1s.val)
                        eventStruct.(evtName).DCIR_diff_10s_1s = eventStruct.(evtName).DCIR_10s.val - eventStruct.(evtName).DCIR_1s.val;
                    else
                        eventStruct.(evtName).DCIR_diff_10s_1s = NaN;
                    end
                end
                
                if isfield(eventStruct.(evtName), 'DCIR_30s') && isfield(eventStruct.(evtName), 'DCIR_1s')
                    if ~isnan(eventStruct.(evtName).DCIR_30s.val) && ~isnan(eventStruct.(evtName).DCIR_1s.val)
                        eventStruct.(evtName).DCIR_diff_30s_1s = eventStruct.(evtName).DCIR_30s.val - eventStruct.(evtName).DCIR_1s.val;
                    else
                        eventStruct.(evtName).DCIR_diff_30s_1s = NaN;
                    end
                end
                
                if isfield(eventStruct.(evtName), 'DCIR_50s') && isfield(eventStruct.(evtName), 'DCIR_1s')
                    if ~isnan(eventStruct.(evtName).DCIR_50s.val) && ~isnan(eventStruct.(evtName).DCIR_1s.val)
                        eventStruct.(evtName).DCIR_diff_50s_1s = eventStruct.(evtName).DCIR_50s.val - eventStruct.(evtName).DCIR_1s.val;
                    else
                        eventStruct.(evtName).DCIR_diff_50s_1s = NaN;
                    end
                end

                filtered_events = filtered_events + 1;
            end
            
            % 파일별 필터링 결과 요약
            fprintf('  Step summary:\n');
            fprintf('    After charging end detection: %d\n', step2_count);
            fprintf('    After duration filter (≥%ds): %d\n', min_charge_duration, step3_count);
            fprintf('    After power stability check: %d\n', step4_count);
            fprintf('    Final valid events from this file: %d\n', step4_count);
        end
    end
end

save(fullfile(saveDir, 'all_chg_events_current_clustering_2023.mat'), 'eventStruct');
fprintf('Step 1 complete. Events saved.\n');
fprintf('Total potential events detected: %d\n', total_events);
fprintf('Events that passed all filters: %d\n', filtered_events);
fprintf('Final stored events: %d\n', eventCount);
fprintf('\n');

%% Step 6: 전류별 군집화 및 시각화
fprintf('=====Step 6: Current-based Clustering =====\n');
% 전류 구간별 이벤트 분류
current_clusters = struct();
for i = 1:length(current_labels)
    current_clusters.(current_labels{i}) = [];
end

eventNames = fieldnames(eventStruct);
fprintf('Step 4: Current-based clustering analysis...\n');

for i = 1:length(eventNames)
    evt = eventStruct.(eventNames{i});
    
    % 평균 충전 전류 계산
    if isfield(evt, 'I_seg_ridIdle') && ~isempty(evt.I_seg_ridIdle)
        avg_current = mean(evt.I_seg_ridIdle(3:end));  % 충전 구간 평균 전류
        
        % 전류 구간 결정
        current_bin_idx = find(avg_current >= current_bins, 1, 'last');
        if ~isempty(current_bin_idx) && current_bin_idx < length(current_bins)
            current_range = current_labels{current_bin_idx};
            current_clusters.(current_range) = [current_clusters.(current_range), i];
            
            % 이벤트에 전류 정보 저장
            eventStruct.(eventNames{i}).avg_current = avg_current;
            eventStruct.(eventNames{i}).current_range = current_range;
        end
    end
end

% 구간별 통계 출력
fprintf('\nCurrent-based clustering results:\n');
for i = 1:length(current_labels)
    cluster_size = length(current_clusters.(current_labels{i}));
    fprintf('  %s: %d events\n', current_labels{i}, cluster_size);
end

%% Step 7: Current-based Visualization
fprintf('=====Step 7: Current-based Visualization =====\n');

% Display labels for plots (with proper formatting)
display_labels = {'50-100A', '100-150A', '150-200A', '200-250A', '250-300A'};

% 1. Current label별 시간-전류, 시간-전압 그래프
figure('Name', 'Current vs Time and Voltage vs Time by Current Range', 'Position', [50, 50, 1200, 800]);

for i = 1:length(current_labels)
    cluster_events = current_clusters.(current_labels{i});
    if isempty(cluster_events)
        continue;
    end
    
    % 각 current label별로 최대 5개 이벤트만 표시
    num_events_to_show = min(5, length(cluster_events));
    
    subplot(2, 3, i); hold on;
    
    for j = 1:num_events_to_show
        evt_idx = cluster_events(j);
        evt = eventStruct.(eventNames{evt_idx});
        
        if isfield(evt, 't_seg_ridIdle') && isfield(evt, 'I_seg_ridIdle')
            t_data = evt.t_seg_ridIdle;
            I_data = evt.I_seg_ridIdle;
            
            % 시간을 초 단위로 변환
            t_seconds = seconds(t_data - t_data(1));
            
            plot(t_seconds, I_data, 'LineWidth', 1.5);
        end
    end
    
    title(sprintf('Current vs Time: %s (%d events)', display_labels{i}, length(cluster_events)));
    xlabel('Time [s]');
    ylabel('Current [A]');
    grid on;
    legend('Location', 'best');
end

sgtitle('Current vs Time by Current Range (2023)');
saveas(gcf, fullfile(saveDir, 'fig_2023_Current_vs_Time_by_range.fig'));

% Voltage vs Time 그래프
figure('Name', 'Voltage vs Time by Current Range', 'Position', [100, 100, 1200, 800]);

for i = 1:length(current_labels)
    cluster_events = current_clusters.(current_labels{i});
    if isempty(cluster_events)
        continue;
    end
    
    % 각 current label별로 최대 5개 이벤트만 표시
    num_events_to_show = min(5, length(cluster_events));
    
    subplot(2, 3, i); hold on;
    
    for j = 1:num_events_to_show
        evt_idx = cluster_events(j);
        evt = eventStruct.(eventNames{evt_idx});
        
        if isfield(evt, 't_seg_ridIdle') && isfield(evt, 'V_seg_ridIdle')
            t_data = evt.t_seg_ridIdle;
            V_data = evt.V_seg_ridIdle;
            
            % 시간을 초 단위로 변환
            t_seconds = seconds(t_data - t_data(1));
            
            plot(t_seconds, V_data, 'LineWidth', 1.5);
        end
    end
    
    title(sprintf('Voltage vs Time: %s (%d events)', display_labels{i}, length(cluster_events)));
    xlabel('Time [s]');
    ylabel('Voltage [V]');
    grid on;
    legend('Location', 'best');
end

sgtitle('Voltage vs Time by Current Range (2023)');
saveas(gcf, fullfile(saveDir, 'fig_2023_Voltage_vs_Time_by_range.fig'));

% DCIR 필드 및 라벨
% 기존 DCIR 분석 시각화 부분을 모두 대체

dcir_fields = {'DCIR_1s', 'DCIR_3s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s', 'DCIR_50s'};
dcir_labels = {'1s', '3s', '5s', '10s', '30s', '50s'};
dcir_diff_fields = {'DCIR_diff_3s_1s', 'DCIR_diff_5s_1s', 'DCIR_diff_10s_1s', 'DCIR_diff_30s_1s', 'DCIR_diff_50s_1s'};
dcir_diff_labels = {'3s-1s', '5s-1s', '10s-1s', '30s-1s', '50s-1s'};

for i = 1:length(current_labels)
    cluster_events = current_clusters.(current_labels{i});
    if isempty(cluster_events)
        continue;
    end

    figure('Name', sprintf('DCIR & Diff - %s', display_labels{i}), 'Position', [200, 200, 1400, 800]);
    % 6+5=11 subplot, 3행 4열, 마지막 한 칸 비움
    for k = 1:length(dcir_fields)
        subplot(3, 4, k); hold on;
        vals = nan(1, length(cluster_events));
        for j = 1:length(cluster_events)
            evt_idx = cluster_events(j);
            evt = eventStruct.(eventNames{evt_idx});
            if isfield(evt, dcir_fields{k}) && isfield(evt.(dcir_fields{k}), 'val')
                vals(j) = evt.(dcir_fields{k}).val;
            end
        end
        plot(1:length(vals), vals, 'o', 'Color', '#0072BD', 'MarkerSize', 1, 'DisplayName', dcir_labels{k});
        m = nanmean(vals);
        s = nanstd(vals);
        % 1시그마 밖 값 x로 표시
        outlier_idx = find(abs(vals - m) > s);
        outliers = vals(outlier_idx);
        plot(outlier_idx, outliers, 'x', 'Color', '#D95319', 'MarkerSize', 6);
        yline(m, '--g', 'LineWidth', 2, 'DisplayName', 'Mean');
        yline(m+s, ':r', 'LineWidth', 2, 'DisplayName', '+1σ');
        yline(m-s, ':r', 'LineWidth', 2, 'DisplayName', '-1σ');
        title(['DCIR ' dcir_labels{k}]);
        xlabel('Event Index');
        ylabel('[mΩ]');
        grid on;
    end

    for k = 1:length(dcir_diff_fields)
        subplot(3, 4, 6+k); hold on;
        vals = nan(1, length(cluster_events));
        for j = 1:length(cluster_events)
            evt_idx = cluster_events(j);
            evt = eventStruct.(eventNames{evt_idx});
            if isfield(evt, dcir_diff_fields{k})
                vals(j) = evt.(dcir_diff_fields{k});
            end
        end
        plot(1:length(vals), vals, 'o', 'Color', '#0072BD', 'MarkerSize', 1, 'DisplayName', dcir_diff_labels{k});
        m = nanmean(vals);
        s = nanstd(vals);
        % 1시그마 밖 값 x로 표시
        outlier_idx = find(abs(vals - m) > s);
        outliers = vals(outlier_idx);
        plot(outlier_idx, outliers, 'x', 'Color', '#D95319', 'MarkerSize', 6);
        yline(m, '--g', 'LineWidth', 2, 'DisplayName', 'Mean');
        yline(m+s, ':r', 'LineWidth', 2, 'DisplayName', '+1σ');
        yline(m-s, ':r', 'LineWidth', 2, 'DisplayName', '-1σ');
        title(['DCIR Diff ' dcir_diff_labels{k}]);
        xlabel('Event Index');
        ylabel('[mΩ]');
        grid on;
    end
    % 마지막 subplot(3,4,12)은 비워둠

    sgtitle(sprintf('DCIR & Diff by Event (%s, %d events)', display_labels{i}, length(cluster_events)));
    saveas(gcf, fullfile(saveDir, sprintf('fig_2023_DCIR_Diff_vs_Event_%s.fig', current_labels{i})));
end

fprintf('\n=== Event Preprocessing Completed ===\n'); ices = event_indices(outlier_mask);
        plot(outlier_event_indices, outliers, 'x', 'Color', '#D95319', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Outliers');
        
        % 이상치 좌표에 값과 이벤트 번호 표기
        for j = 1:length(outlier_event_indices)
            text(outlier_event_indices(j), outliers(j), sprintf('E%d: %.1f', outlier_event_indices(j), outliers(j)), ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
                'FontSize', 8, 'Color', '#D95319', 'FontWeight', 'bold');
        end
    end
    
    % 통계선 그리기
    xlim_vals = xlim;
    
    % 원본 평균선 (점선)
    yline(data_mean, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', sprintf('Original Mean (%.2f mΩ)', data_mean));
    
    % 정상 데이터 평균선 (실선)
    yline(clean_mean, '-', 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', sprintf('Clean Mean (%.2f mΩ)', clean_mean));
    
    % σ 경계선
    yline(data_mean + sigma_threshold * data_std, ':', 'Color', '#D95319', 'LineWidth', 1, 'DisplayName', sprintf('+%dσ (%.2f mΩ)', sigma_threshold, data_std));
    yline(data_mean - sigma_threshold * data_std, ':', 'Color', '#D95319', 'LineWidth', 1, 'DisplayName', sprintf('-%dσ (%.2f mΩ)', sigma_threshold, data_std));
    
    % 제목 및 레이블
    title(sprintf('DCIR %s (Outliers: %d/%d, %.1f%% removed)', dcir_labels{k}, length(outliers), length(valid_data), (length(outliers)/length(valid_data))*100));
    xlabel('Event Index');
    ylabel('[mΩ]');
    
    % x축 설정: 전체 이벤트 범위, 간격 5
    xticks(1:5:nEvents);  % 1, 6, 11, 16, 21, ...
    xlim([0.5, nEvents + 0.5]);
    
    % NaN 값 표시 (Legend에만)
    if ~isempty(nan_idx)
        % NaN 이벤트 리스트 생성 (줄당 최대 6개)
        nan_events_str = '';
        for j = 1:length(nan_idx)
            if j == 1
                nan_events_str = sprintf('E%d', nan_idx(j));
            else
                % 6개마다 줄바꿈
                if mod(j-1, 6) == 0
                    nan_events_str = sprintf('%s,\nE%d', nan_events_str, nan_idx(j));
                else
                    nan_events_str = sprintf('%s, E%d', nan_events_str, nan_idx(j));
                end
            end
        end
        
        % Legend에 NaN 정보 표시 (빨간색 X 마커)
        plot(NaN, NaN, 'x', 'Color', 'red', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'NaN');
        
        % NaN 이벤트 리스트를 별도 legend 항목으로 표시
        plot(NaN, NaN, 'w.', 'MarkerSize', 0.1, 'DisplayName', nan_events_str);
    end
    
    grid on;
    
    legend('Location', 'eastoutside', 'FontSize', 8);
    
    % 통계 정보 출력
    fprintf('%s: %d outliers (%.1f%%) | Original: μ=%.2f, σ=%.2f | Clean: μ=%.2f, σ=%.2f\n', ...
        dcir_fields{k}, length(outliers), (length(outliers)/length(valid_data))*100, ...
        data_mean, data_std, clean_mean, clean_std);
end

sgtitle(sprintf('DCIR Outlier Removal using %dσ (2023.06)', sigma_threshold));
saveas(gcf, fullfile(saveDir, 'fig_2023_DCIR_outlier_removal.fig'));

% 이상치 통계 저장
save(fullfile(saveDir, 'outlier_stats_2023.mat'), 'outlier_stats');
fprintf('\nOutlier statistics saved to: outlier_stats_2023.mat\n');

%% 이상치 제거 전후 비교 히스토그램
figure('Name', 'DCIR Distribution Comparison', 'Position', [200, 200, 1400, 800]);

for k = 1:length(dcir_fields)
    dcir_vals = NaN(nEvents,1);
    
    % 데이터 수집
    for i = 1:nEvents
        evt = eventStruct.(eventNames{i});
        if isfield(evt, dcir_fields{k}) && isfield(evt.(dcir_fields{k}), 'val')
            dcir_vals(i) = evt.(dcir_fields{k}).val;
        end
    end
    
    % 유효한 데이터만 추출
    valid_idx = ~isnan(dcir_vals);
    valid_data = dcir_vals(valid_idx);
    
    if length(valid_data) < 3
        continue;
    end
    
    % 통계 계산
    data_mean = mean(valid_data);
    data_std = std(valid_data);
    
    % 이상치 검출
    outlier_mask = abs(valid_data - data_mean) > (sigma_threshold * data_std);
    normal_data = valid_data(~outlier_mask);
    
    % 히스토그램 플롯
    subplot(4,2,k); hold on;
    
    % 원본 데이터 히스토그램 (투명한 회색)
    histogram(valid_data, 10, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'black', 'DisplayName', 'Original');
    
    % 정상 데이터 히스토그램 (파란색)
    histogram(normal_data, 10, 'FaceColor', '#0072BD', 'FaceAlpha', 0.7, 'EdgeColor', 'blue', 'DisplayName', 'Clean');
    
    % 통계선
    xline(data_mean, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'Original Mean');
    xline(mean(normal_data), '-', 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', 'Clean Mean');
    
    title(sprintf('DCIR %s Distribution', dcir_labels{k}));
    xlabel('[mΩ]');
    ylabel('Frequency');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
end

sgtitle('DCIR Distribution: Original vs Clean Data (2023.06)');
saveas(gcf, fullfile(saveDir, 'fig_2023_DCIR_distribution_comparison.fig'));

fprintf('\n=== 이상치 제거 완료 ===\n');
