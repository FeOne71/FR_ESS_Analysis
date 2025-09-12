%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field Data DCIR Discharge ver01.m
% ESS 방전 이벤트 및 저항 분석 (BSC_Discharge 기반)
% Step 1: Folder traversal
% Step 2: BSC_Discharge 기반 이벤트 검출
% Step 3: DCIR 계산
% Step 4: Visualization
% 
% Ver01 Updates:
% - Added Power-based filtering (max_P, min_P, max_P_std)
% - Extended Condition #3: Current + Power magnitude check
% - Extended Condition #4: Current + Power stability check
% - Power data stored in P_seq and P_seg_ridIdle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average SOC(%)	Highest SOC(%)	Lowest SOC(%)	Highest SOC Pos.(R)	Lowest SOC Pos.(R)	Average SOH(%)	Average C.V. Sum(V)	DC Current(A)	Highest DC Current(A)	Lowest DC Current(A)	Highest DC Current Pos.(R)	Lowest DC Current Pos.(R)	DC Chg. Current Limit(A)	DC Dchg. Current Limit(A)	DC Power(kW)

clc; clear; close all;

%% Directory
dataDir  = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\raw2mat_ver04';
yearList = {'2021', '2022', '2023'}; % ,'2023','2023'};
saveDir  = fullfile('G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData_DCIR_Charge');

if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

%% Variables
C_nom        = 1024;          % Ah
min_discharge_duration = 30;  % [s] - Discharging duration
max_I        = 0.20 * C_nom;  % Max discharge current               [256.0 (A)]
min_I        = 0.10 * C_nom;  % min discharge current               [102.4 (A)]
std_I        = 0.01 * C_nom;  % Allowed std for discharging current [51.2  (A)]
max_P        = 250;           % Max discharge power [kW]
min_P        = 150;           % Min discharge power [kW]
max_P_std    = 50;            % Max power standard deviation [kW]
dt_sec       = 1;             % Sampling time [sec]

%% Step 1: Raw Files Traversal
fprintf('Step 1: Detecting discharging events using BSC_Discharge...\n');
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
            bsc_discharge = Raw.Discharge;  % string array
            status     = Raw.Status;
            dc_power   = Raw.DCPower_kW;
            
%% Step 2: Discharging Event Detection
            fprintf('=====Step 2: Searching Idle to Discharging Transition =====\n');
            % 1) Define "Idle" and "Discharging" status
            is_idle = strcmp(bsc_discharge, 'Idle');
            is_discharging = strcmp(bsc_discharge, 'Discharging');
            
            % Detect "Idle" -> "Discharging"
            idle_to_discharge = find(is_idle(1:end-1) & is_discharging(2:end));
            total_events = total_events + length(idle_to_discharge);
            fprintf('File: %s - Found %d potential events\n', matFiles(f).name, length(idle_to_discharge));
            
            fprintf('  Total Idle points: %d\n', sum(is_idle));
            fprintf('  Total Discharging points: %d\n', sum(is_discharging));
            fprintf('  Idle→Discharging transitions: %d\n', length(idle_to_discharge));

            step2_count = 0;
            step3_count = 0;
            step4_count = 0;
            step5_count = 0;
            step6_count = 0;

%% Step 3: Events Filtering              
            fprintf('=====Step 3: Filtering Events =====\n');
            for i = 1:length(idle_to_discharge)
                idx1 = idle_to_discharge(i);  % "Idle"의 마지막 시점
                start_discharge_idx = idx1 + 1;  % "Discharging" 시작점
                
                % 2) "Discharging" 구간 끝점 찾기
                discharge_end_idx = start_discharge_idx;
                while discharge_end_idx <= length(bsc_discharge) && strcmp(bsc_discharge(discharge_end_idx), 'Discharging')
                    discharge_end_idx = discharge_end_idx + 1;
                end
                discharge_end_idx = discharge_end_idx - 1;  % "Discharging"의 마지막 시점
                
                % 이벤트 전체 구간 (Idle 시작점부터 Discharging 끝점까지)
                start_idx = idx1;
                end_idx = discharge_end_idx;
                
                step2_count = step2_count + 1;
                
                % Condition #1 "Discharging" 구간 지속 시간 확인 (30초 이상)
                discharge_duration = discharge_end_idx - start_discharge_idx + 1;
                if discharge_duration < min_discharge_duration % 30sec
                    fprintf('  Event %d filtered: Discharge duration too short (%ds)\n', i, discharge_duration);
                    continue;
                end
                
                step3_count = step3_count + 1;
                
                t_seg   = t(start_idx:end_idx);
                I_seg   = I(start_idx:end_idx);
                V_seg   = V(start_idx:end_idx);
                soc_seg = soc(start_idx:end_idx);
                T_seg   = T_batt(start_idx:end_idx);
                bsc_seg = bsc_discharge(start_idx:end_idx);
                P_seg   = dc_power(start_idx:end_idx);
                
                duration_sec = seconds(t_seg(end) - t_seg(1));
                deltaI = max(I_seg) - min(I_seg);
                I_max = max(abs(I_seg));
                
                fprintf('Event %d: Duration=%.1fs, Discharge_Duration=%ds, deltaI=%.2fA, I_max=%.2fA\n', i, duration_sec, discharge_duration, deltaI, I_max);
                
                % Condition #2: 전류 크기 확인 (Discharging 구간만 체크)
                idx2 = discharge_end_idx - start_idx + 1;  % Discharging 구간 내 상대 위치
                
                max_current = max(abs(I_seg(1:idx2)));
                min_current = min(abs(I_seg(1:idx2)));
                
                current_ok = (max_current <= max_I) && (min_current >= min_I);
                
                fprintf('  Event %d: max_current=%.1fA (limit=%.1fA), min_current=%.1fA (limit=%.1fA)\n', ...
                    i, max_current, max_I, min_current, min_I);
                
                if ~current_ok
                    if max_current > max_I
                        fprintf('    → FILTERED: max current too high (%.1fA > %.1fA)\n', max_current, max_I);
                    elseif min_current < min_I
                        fprintf('    → FILTERED: min current too low (%.1fA < %.1fA)\n', min_current, min_I);
                    end
                    continue;
                end
                
                % Condition #2-1: 출력 크기 확인 (Discharging 구간만 체크)
                max_power = max(abs(P_seg(1:idx2)));
                min_power = min(abs(P_seg(1:idx2)));
                
                power_ok = (max_power <= max_P) && (min_power >= min_P);
                
                fprintf('  Event %d: max_power=%.1fkW (limit=%.1fkW), min_power=%.1fkW (limit=%.1fkW)\n', ...
                    i, max_power, max_P, min_power, min_P);
                
                if ~power_ok
                    if max_power > max_P
                        fprintf('    → FILTERED: max power too high (%.1fkW > %.1fkW)\n', max_power, max_P);
                    elseif min_power < min_P
                        fprintf('    → FILTERED: min power too low (%.1fkW < %.1fkW)\n', min_power, min_P);
                    end
                    continue;
                end

                step4_count = step4_count + 1;                               

                % Condition #3: 전류 및 출력 안정성 (Discharging 구간만 체크)
                current_std = std(I_seg(3:idx2));
                power_std   = std(P_seg(3:idx2));
                
                current_stable = current_std <= std_I;
                power_stable   = power_std <= max_P_std;
                
                fprintf('  Event %d: current_std=%.1fA (limit=%.1fA), power_std=%.1fkW (limit=%.1fkW)\n', ...
                    i, current_std, std_I, power_std, max_P_std);

                if ~current_stable %|| ~power_stable
                    fprintf('    → FILTERED: Current instability (std=%.1fA > %.1fA)\n', current_std, std_I);
                    continue;    
                end
                step5_count = step5_count + 1;

                % if ~power_stable
                %     fprintf('    → FILTERED: Power instability (std=%.1fkW > %.1fkW)\n', power_std, max_P_std);
                %     continue;
                % end
                % 
                % step6_count = step6_count + 1;

%% Step 4: DCIR Calculation
                fprintf('=====Step 4: DCIR Calculation =====\n');
                % Save Events as Struct
                eventCount = eventCount + 1;
                evtName = sprintf('event%d', eventCount);
                eventStruct.(evtName).start_idx = start_idx;
                eventStruct.(evtName).end_idx = end_idx;
                eventStruct.(evtName).idx1 = idx1;  % Idle의 마지막 시점
                eventStruct.(evtName).idx2 = idx2;  % Discharging의 마지막 시점 (세그먼트 내 상대위치)
                eventStruct.(evtName).discharge_duration = discharge_duration;

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

                        if (dI) < 0 && dV < 0  % Discharge: current decreases, voltage decreases
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

                % Discharging 끝점 기준 DCIR 계산
                if idx2 <= length(I_seg)
                    V1_r = V_seg(1);
                    V2_r = V_seg(idx2);
                    I1_r = I_seg(1);
                    I2_r = I_seg(idx2);
                    dV_ramp = V2_r - V1_r;
                    dI_ramp = I2_r - I1_r;

                    if (dV_ramp) < 0 && (dI_ramp) < 0  % Discharge: voltage and current decrease
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

                fprintf('  Event %d: DCIR_1s=%.2fmΩ, DCIR_5s=%.2fmΩ, DCIR_10s=%.2fmΩ\n', ...
                    i, eventStruct.(evtName).DCIR_1s.val, eventStruct.(evtName).DCIR_5s.val, eventStruct.(evtName).DCIR_10s.val);
                
                filtered_events = filtered_events + 1;
            end
            
            fprintf('  File %s: Step2=%d, Step3=%d, Step4=%d, Step5=%d, Step6=%d\n', ...
                matFiles(f).name, step2_count, step3_count, step4_count, step5_count, step6_count);
        end
    end
end

%% Step 6: Save Results
fprintf('=====Step 6: Save Results =====\n');
fprintf('Total events detected: %d\n', total_events);
fprintf('Filtered events: %d\n', filtered_events);
fprintf('Event count: %d\n', eventCount);

% Save event structure
save(fullfile(saveDir, 'all_discharge_events.mat'), 'eventStruct');
fprintf('Event structure saved to: %s\n', fullfile(saveDir, 'all_discharge_events.mat'));

%% Step 7: Basic Visualization
fprintf('=====Step 7: Basic Visualization =====\n');

% DCIR values extraction
dcir_1s = [];
dcir_3s = [];
dcir_5s = [];
dcir_10s = [];
dcir_30s = [];
dcir_50s = [];

event_names = fieldnames(eventStruct);
for i = 1:length(event_names)
    evt = eventStruct.(event_names{i});
    
    if isfield(evt, 'DCIR_1s') && ~isnan(evt.DCIR_1s.val)
        dcir_1s = [dcir_1s, evt.DCIR_1s.val];
    end
    if isfield(evt, 'DCIR_3s') && ~isnan(evt.DCIR_3s.val)
        dcir_3s = [dcir_3s, evt.DCIR_3s.val];
    end
    if isfield(evt, 'DCIR_5s') && ~isnan(evt.DCIR_5s.val)
        dcir_5s = [dcir_5s, evt.DCIR_5s.val];
    end
    if isfield(evt, 'DCIR_10s') && ~isnan(evt.DCIR_10s.val)
        dcir_10s = [dcir_10s, evt.DCIR_10s.val];
    end
    if isfield(evt, 'DCIR_30s') && ~isnan(evt.DCIR_30s.val)
        dcir_30s = [dcir_30s, evt.DCIR_30s.val];
    end
    if isfield(evt, 'DCIR_50s') && ~isnan(evt.DCIR_50s.val)
        dcir_50s = [dcir_50s, evt.DCIR_50s.val];
    end
end

% Create boxplot
figure('Name', 'DCIR Distribution - Discharge Events', 'Position', [100, 100, 1200, 800]);
subplot(2, 3, 1); boxplot(dcir_1s); title('DCIR 1s'); ylabel('DCIR [mΩ]');
subplot(2, 3, 2); boxplot(dcir_3s); title('DCIR 3s'); ylabel('DCIR [mΩ]');
subplot(2, 3, 3); boxplot(dcir_5s); title('DCIR 5s'); ylabel('DCIR [mΩ]');
subplot(2, 3, 4); boxplot(dcir_10s); title('DCIR 10s'); ylabel('DCIR [mΩ]');
subplot(2, 3, 5); boxplot(dcir_30s); title('DCIR 30s'); ylabel('DCIR [mΩ]');
subplot(2, 3, 6); boxplot(dcir_50s); title('DCIR 50s'); ylabel('DCIR [mΩ]');

sgtitle('DCIR Distribution for Discharge Events', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(saveDir, 'fig_DCIR_distribution_discharge.fig'));

% Statistics
fprintf('DCIR Statistics:\n');
fprintf('1s:  Mean=%.2f, Std=%.2f, Count=%d\n', mean(dcir_1s), std(dcir_1s), length(dcir_1s));
fprintf('3s:  Mean=%.2f, Std=%.2f, Count=%d\n', mean(dcir_3s), std(dcir_3s), length(dcir_3s));
fprintf('5s:  Mean=%.2f, Std=%.2f, Count=%d\n', mean(dcir_5s), std(dcir_5s), length(dcir_5s));
fprintf('10s: Mean=%.2f, Std=%.2f, Count=%d\n', mean(dcir_10s), std(dcir_10s), length(dcir_10s));
fprintf('30s: Mean=%.2f, Std=%.2f, Count=%d\n', mean(dcir_30s), std(dcir_30s), length(dcir_30s));
fprintf('50s: Mean=%.2f, Std=%.2f, Count=%d\n', mean(dcir_50s), std(dcir_50s), length(dcir_50s));

fprintf('Discharge event detection and DCIR calculation completed successfully.\n'); 