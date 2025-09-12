%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Postprocessing - Generate and Save Figures by Category
% NEW VERSION: SOC calculation using trapz/64Ah
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% File Directory
folderPath = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data\RPT';
savePath   = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\RPT_Postprocessing';

% Change save path to current directory to resolve path issues
currentDir = pwd;
savePath = fullfile(currentDir, 'RPT_Postprocessing_Results_fixedCap');

channels = {'Ch9', 'Ch10', 'Ch11', 'Ch12', 'Ch13', 'Ch14', 'Ch15', 'Ch16'};
rpt_cycles = {'0cyc', '200cyc'};
rpt_cycles_fields = {'cyc0', 'cyc200'}; % For field names

Time      = [];
totalTime = [];

% DCIR-SOC 데이터 저장을 위한 구조체
dcir_soc_data = struct();

%% 3. DCIR vs SOC (NEW VERSION)
fprintf('\n=== SOC-DCIR Processing (NEW VERSION) ===\n');
dcir_socFolder = fullfile(savePath, 'DCIR');
if ~exist(dcir_socFolder, 'dir')
    mkdir(dcir_socFolder);
end

% DCIR step conditions
dcir_charge_step = 14;   % 충전 DCIR 펄스
dcir_discharge_step = 21; % 방전 DCIR 펄스

% SOC calculation steps
soc_charge_step = 16;     % 충전 SOC 적산 구간
soc_discharge_step = 23;  % 방전 SOC 적산 구간

% Fixed capacity for SOC calculation
fixed_capacity = 64; % Ah

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};

    figure('Name', sprintf('DCIR vs SOC - %s (NEW)', channel), 'Position', [100 100 1200 800]);
    set(gcf, 'Visible', 'off');

    % Subplot 1: Charge DCIR
    subplot(2,1,1);
    hold on;

    % Subplot 2: Discharge DCIR
    subplot(2,1,2);
    hold on;

    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
        filepath = fullfile(folderPath, filename);

        T = readtable(filepath);

        fprintf('Channel: %s, RPT: %s\n', channel, rpt_cycle);

        % === 충전 DCIR vs SOC 처리 ===
        % StepIdx=14에서 DCIR 값들을 순서대로 추출
        charge_dcir_mask = T.StepIndex == dcir_charge_step;
        charge_dcir_data = T(charge_dcir_mask, :);

        if height(charge_dcir_data) > 0
            % DCIR cycle들을 순서대로 정렬
            charge_dcir_cycles = unique(charge_dcir_data.CycleIndex);
            charge_dcir_cycles = sort(charge_dcir_cycles);
            fprintf('Charge DCIR cycles (ordered): ');
            fprintf('%d ', charge_dcir_cycles);
            fprintf('\n');

            % StepIdx=16에서 SOC 계산을 위한 데이터
            charge_soc_mask = T.StepIndex == soc_charge_step;
            charge_soc_data = T(charge_soc_mask, :);

            if height(charge_soc_data) > 0
                % 각 DCIR cycle별로 순서대로 처리
                for i = 1:length(charge_dcir_cycles)
                    cycle_idx = charge_dcir_cycles(i);

                    % 해당 cycle의 DCIR 데이터
                    cycle_dcir_mask = charge_dcir_data.CycleIndex == cycle_idx;
                    cycle_dcir = charge_dcir_data(cycle_dcir_mask, :);

                    if height(cycle_dcir) > 0
                        % DCIR 값 계산 - 정확한 시점에서 계산
                        % 1분 펄스에서 0.1s, 5s, 10s, 30s, 60s 시점의 DCIR 추출
                        dcir_time_points = [0.1, 5, 10, 30, 60]; % seconds
                        dcir_values_01s = NaN;
                        dcir_values_5s = NaN;
                        dcir_values_10s = NaN;
                        dcir_values_30s = NaN;
                        dcir_values_60s = NaN;

                        % 시간 간격 저장 및 상대시간 계산
                        time_intervals = diff(seconds(cycle_dcir.TotalTime));
                        relative_times = [0; cumsum(time_intervals)]; % 0초부터 시작하는 상대시간

                        for time_idx = 1:length(dcir_time_points)
                            target_time = dcir_time_points(time_idx);

                            % 정확한 시점의 데이터 찾기
                            exact_idx = find(abs(relative_times - target_time) < 0.01, 1);

                            if exact_idx <= height(cycle_dcir)
                                % DCIR 펄스 시작값 (V1, I1)
                                v1 = cycle_dcir.Voltage_V_(1);  % 펄스 시작 전압
                                i1 = cycle_dcir.Current_A_(1);  % 펄스 시작 전류
                                
                                % 해당 시점값 (V2, I2)
                                v2 = cycle_dcir.Voltage_V_(exact_idx);  % 해당 시점 전압
                                i2 = cycle_dcir.Current_A_(exact_idx);  % 해당 시점 전류
                                
                                % DCIR 계산 - 단순화된 공식
                                % DCIR = (V2 - V1) / I2
                                dcir_value = (v2 - v1) / i2;
                                
                                % 시간별로 정확히 저장
                                if abs(target_time - 0.1) < 0.01
                                    dcir_values_01s = dcir_value;
                                elseif abs(target_time - 5) < 0.01
                                    dcir_values_5s = dcir_value;
                                elseif abs(target_time - 10) < 0.01
                                    dcir_values_10s = dcir_value;
                                elseif abs(target_time - 30) < 0.01
                                    dcir_values_30s = dcir_value;
                                elseif abs(target_time - 60) < 0.01
                                    dcir_values_60s = dcir_value;
                                end
                            end
                        end

                        % SOC 계산 - 순서에 따라
                        if i == 1
                            % 첫 번째 DCIR → SOC 0%
                            soc_percent = 0;
                        else
                            % 두 번째부터는 실제 적분 계산
                            % 해당 cycle의 SOC 구간 데이터
                            cycle_soc_mask = charge_soc_data.CycleIndex == cycle_idx;
                            cycle_soc = charge_soc_data(cycle_soc_mask, :);

                            if height(cycle_soc) > 0
                                % Cycle 2부터 (현재 cycle-1)까지의 누적 적분 계산
                                total_integral = 0;
                                for prev_cycle = 2:cycle_idx-1
                                    prev_cycle_mask = charge_soc_data.CycleIndex == prev_cycle;
                                    prev_cycle_data = charge_soc_data(prev_cycle_mask, :);
                                    
                                    if height(prev_cycle_data) > 0
                                        % 시간을 hours로 변환하여 적분
                                        time_hours = seconds(prev_cycle_data.TotalTime - prev_cycle_data.TotalTime(1)) / 3600;
                                        cycle_integral = trapz(time_hours, prev_cycle_data.Current_A_);
                                        total_integral = total_integral + cycle_integral;
                                    end
                                end
                                
                                % 디버깅 출력
                                fprintf('  RPT%s - Cycle %d: Total=%.6f, Fixed capacity=%.1f Ah, Ratio=%.3f%%\n', ...
                                    rpt_cycle, charge_dcir_cycles(i), total_integral, fixed_capacity, (total_integral/fixed_capacity)*100);
                                
                                % 충전 SOC = 0 + trapz(total_current) / 64Ah
                                soc_percent = (total_integral / fixed_capacity) * 100;
                            else
                                soc_percent = 0;
                            end
                        end

                        fprintf('Charge DCIR #%d (Cycle %d): SOC = %.1f%%, DCIR values calculated\n', ...
                            i, cycle_idx, soc_percent);

                        % RPT cycle별 필드명 설정
                        if strcmp(rpt_cycle, '0cyc')
                            rpt_field = 'cyc0';
                        else
                            rpt_field = 'cyc200';
                        end

                        % DCIR-SOC 데이터 구조체 초기화
                        if ~isfield(dcir_soc_data, channel)
                            dcir_soc_data.(channel) = struct();
                        end
                        
                        if ~isfield(dcir_soc_data.(channel), rpt_field)
                            dcir_soc_data.(channel).(rpt_field) = struct();
                        end
                        
                        % 충전 DCIR-SOC 데이터 필드 초기화
                        if ~isfield(dcir_soc_data.(channel).(rpt_field), 'charge_soc')
                            dcir_soc_data.(channel).(rpt_field).charge_soc = [];
                            dcir_soc_data.(channel).(rpt_field).charge_dcir_01s = [];
                            dcir_soc_data.(channel).(rpt_field).charge_dcir_5s = [];
                            dcir_soc_data.(channel).(rpt_field).charge_dcir_10s = [];
                            dcir_soc_data.(channel).(rpt_field).charge_dcir_30s = [];
                            dcir_soc_data.(channel).(rpt_field).charge_dcir_60s = [];
                            dcir_soc_data.(channel).(rpt_field).charge_cycles = [];
                        end

                        % DCIR-SOC 데이터 저장
                        fprintf('  DCIR values check: 01s=%.6f, 5s=%.6f, 10s=%.6f, 30s=%.6f, 60s=%.6f\n', ...
                            dcir_values_01s, dcir_values_5s, dcir_values_10s, dcir_values_30s, dcir_values_60s);
                        
                        % UPDATED와 동일하게 조건 없이 저장 (NaN 값도 포함)
                        dcir_soc_data.(channel).(rpt_field).charge_soc = [dcir_soc_data.(channel).(rpt_field).charge_soc; soc_percent];
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_01s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_01s; dcir_values_01s];
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_5s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_5s; dcir_values_5s];
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_10s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_10s; dcir_values_10s];
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_30s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_30s; dcir_values_30s];
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_60s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_60s; dcir_values_60s];
                        dcir_soc_data.(channel).(rpt_field).charge_cycles = [dcir_soc_data.(channel).(rpt_field).charge_cycles; cycle_idx];
                        fprintf('  Data saved for cycle %d (including NaN values)\n', cycle_idx);
                    end
                end
            end
            
            % === 방전 DCIR vs SOC 처리 ===
            % StepIdx=21에서 DCIR 값들을 순서대로 추출
            discharge_dcir_mask = T.StepIndex == dcir_discharge_step;
            discharge_dcir_data = T(discharge_dcir_mask, :);

            if height(discharge_dcir_data) > 0
                % DCIR cycle들을 순서대로 정렬 (16 이상만)
                discharge_dcir_cycles = unique(discharge_dcir_data.CycleIndex);
                discharge_dcir_cycles = discharge_dcir_cycles(discharge_dcir_cycles >= 16);
                discharge_dcir_cycles = sort(discharge_dcir_cycles);
                fprintf('Discharge DCIR cycles (ordered): ');
                fprintf('%d ', discharge_dcir_cycles);
                fprintf('\n');

                % StepIdx=23에서 SOC 계산을 위한 데이터
                discharge_soc_mask = T.StepIndex == soc_discharge_step;
                discharge_soc_data = T(discharge_soc_mask, :);

                if height(discharge_soc_data) > 0
                    % 각 DCIR cycle별로 순서대로 처리
                    for i = 1:length(discharge_dcir_cycles)
                        cycle_idx = discharge_dcir_cycles(i);

                        % 해당 cycle의 DCIR 데이터
                        cycle_dcir_mask = discharge_dcir_data.CycleIndex == cycle_idx;
                        cycle_dcir = discharge_dcir_data(cycle_dcir_mask, :);

                        if height(cycle_dcir) > 0
                            % DCIR 값 계산 - 정확한 시점에서 계산
                            % 1분 펄스에서 0.1s, 5s, 10s, 30s, 60s 시점의 DCIR 추출
                            dcir_time_points = [0.1, 5, 10, 30, 60]; % seconds
                            dcir_values_01s = NaN;
                            dcir_values_5s = NaN;
                            dcir_values_10s = NaN;
                            dcir_values_30s = NaN;
                            dcir_values_60s = NaN;

                            % 시간 간격 저장 및 상대시간 계산
                            time_intervals = diff(seconds(cycle_dcir.TotalTime));
                            relative_times = [0; cumsum(time_intervals)]; % 0초부터 시작하는 상대시간

                            for time_idx = 1:length(dcir_time_points)
                                target_time = dcir_time_points(time_idx);

                                % 정확한 시점의 데이터 찾기
                                exact_idx = find(abs(relative_times - target_time) < 0.01, 1);

                                if exact_idx <= height(cycle_dcir)
                                    % DCIR 펄스 시작값 (V1, I1)
                                    v1 = cycle_dcir.Voltage_V_(1);  % 펄스 시작 전압
                                    i1 = cycle_dcir.Current_A_(1);  % 펄스 시작 전류
                                    
                                    % 해당 시점값 (V2, I2)
                                    v2 = cycle_dcir.Voltage_V_(exact_idx);  % 해당 시점 전압
                                    i2 = cycle_dcir.Current_A_(exact_idx);  % 해당 시점 전류
                                    
                                    % DCIR 계산 - 단순화된 공식
                                    % DCIR = (V2 - V1) / I2
                                    dcir_value = (v2 - v1) / i2;
                                    
                                    % 시간별로 정확히 저장
                                    if abs(target_time - 0.1) < 0.01
                                        dcir_values_01s = dcir_value;
                                    elseif abs(target_time - 5) < 0.01
                                        dcir_values_5s = dcir_value;
                                    elseif abs(target_time - 10) < 0.01
                                        dcir_values_10s = dcir_value;
                                    elseif abs(target_time - 30) < 0.01
                                        dcir_values_30s = dcir_value;
                                    elseif abs(target_time - 60) < 0.01
                                        dcir_values_60s = dcir_value;
                                    end
                                end
                            end

                            % SOC 계산 - 순서에 따라
                            if i == 1
                                % 첫 번째 DCIR → SOC 100%
                                soc_percent = 100;
                            else
                                % 두 번째부터는 실제 적분 계산
                                % 해당 cycle의 SOC 구간 데이터
                                cycle_soc_mask = discharge_soc_data.CycleIndex == cycle_idx;
                                cycle_soc = discharge_soc_data(cycle_soc_mask, :);

                                if height(cycle_soc) > 0
                                    % Cycle 16부터 (현재 cycle-1)까지의 누적 적분 계산
                                    total_integral = 0;
                                    for prev_cycle = 16:cycle_idx-1
                                        prev_cycle_mask = discharge_soc_data.CycleIndex == prev_cycle;
                                        prev_cycle_data = discharge_soc_data(prev_cycle_mask, :);
                                        
                                        if height(prev_cycle_data) > 0
                                            % 시간을 hours로 변환하여 적분 (방전은 절댓값)
                                            time_hours = seconds(prev_cycle_data.TotalTime - prev_cycle_data.TotalTime(1)) / 3600;
                                            cycle_integral = abs(trapz(time_hours, prev_cycle_data.Current_A_));
                                            total_integral = total_integral + cycle_integral;
                                        end
                                    end
                                    
                                    % 디버깅 출력
                                    fprintf('  RPT%s - Cycle %d: Total=%.6f, Fixed capacity=%.1f Ah, Ratio=%.3f%%, SOC=%.1f%%\n', ...
                                        rpt_cycle, discharge_dcir_cycles(i), total_integral, fixed_capacity, (total_integral/fixed_capacity)*100, 100-(total_integral/fixed_capacity)*100);
                                    
                                    % 방전 SOC = 1 - trapz(total_current) / 64Ah
                                    soc_percent = 100 - (total_integral / fixed_capacity) * 100;
                                else
                                    soc_percent = 100;
                                end
                            end

                            fprintf('Discharge DCIR #%d (Cycle %d): SOC = %.1f%%, DCIR values calculated\n', ...
                                i, cycle_idx, soc_percent);

                            % RPT cycle별 필드명 설정
                            if strcmp(rpt_cycle, '0cyc')
                                rpt_field = 'cyc0';
                            else
                                rpt_field = 'cyc200';
                            end

                            % 방전 DCIR-SOC 데이터 필드 초기화 (필드가 없을 때만)
                            if ~isfield(dcir_soc_data.(channel).(rpt_field), 'discharge_soc')
                                dcir_soc_data.(channel).(rpt_field).discharge_soc = [];
                                dcir_soc_data.(channel).(rpt_field).discharge_dcir_01s = [];
                                dcir_soc_data.(channel).(rpt_field).discharge_dcir_5s = [];
                                dcir_soc_data.(channel).(rpt_field).discharge_dcir_10s = [];
                                dcir_soc_data.(channel).(rpt_field).discharge_dcir_30s = [];
                                dcir_soc_data.(channel).(rpt_field).discharge_dcir_60s = [];
                                dcir_soc_data.(channel).(rpt_field).discharge_cycles = [];
                            end

                            % DCIR-SOC 데이터 저장
                            fprintf('  Discharge DCIR values check: 01s=%.6f, 5s=%.6f, 10s=%.6f, 30s=%.6f, 60s=%.6f\n', ...
                                dcir_values_01s, dcir_values_5s, dcir_values_10s, dcir_values_30s, dcir_values_60s);
                            
                            % UPDATED와 동일하게 조건 없이 저장 (NaN 값도 포함)
                            dcir_soc_data.(channel).(rpt_field).discharge_soc = [dcir_soc_data.(channel).(rpt_field).discharge_soc; soc_percent];
                            dcir_soc_data.(channel).(rpt_field).discharge_dcir_01s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_01s; dcir_values_01s];
                            dcir_soc_data.(channel).(rpt_field).discharge_dcir_5s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_5s; dcir_values_5s];
                            dcir_soc_data.(channel).(rpt_field).discharge_dcir_10s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_10s; dcir_values_10s];
                            dcir_soc_data.(channel).(rpt_field).discharge_dcir_30s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_30s; dcir_values_30s];
                            dcir_soc_data.(channel).(rpt_field).discharge_dcir_60s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_60s; dcir_values_60s];
                            dcir_soc_data.(channel).(rpt_field).discharge_cycles = [dcir_soc_data.(channel).(rpt_field).discharge_cycles; cycle_idx];
                            fprintf('  Discharge data saved for cycle %d (including NaN values)\n', cycle_idx);
                        end
                    end
                end
            end
        end
    end

    % === 저장된 데이터를 사용한 시각화 ===
    fprintf('\n=== Plotting from saved data ===\n');
    
    % 데이터 구조체 디버깅
    fprintf('=== Data Structure Debug ===\n');
    if isfield(dcir_soc_data, channel)
        fprintf('Channel %s exists in dcir_soc_data\n', channel);
        if isfield(dcir_soc_data.(channel), 'cyc0')
            fprintf('cyc0 exists for channel %s\n', channel);
            if isfield(dcir_soc_data.(channel).cyc0, 'charge_soc')
                fprintf('cyc0.charge_soc length: %d\n', length(dcir_soc_data.(channel).cyc0.charge_soc));
                if ~isempty(dcir_soc_data.(channel).cyc0.charge_soc)
                    fprintf('cyc0.charge_soc values: ');
                    fprintf('%.1f ', dcir_soc_data.(channel).cyc0.charge_soc);
                    fprintf('\n');
                    fprintf('cyc0.charge_dcir_60s values: ');
                    fprintf('%.6f ', dcir_soc_data.(channel).cyc0.charge_dcir_60s);
                    fprintf('\n');
                end
            end
        end
        if isfield(dcir_soc_data.(channel), 'cyc200')
            fprintf('cyc200 exists for channel %s\n', channel);
            if isfield(dcir_soc_data.(channel).cyc200, 'charge_soc')
                fprintf('cyc200.charge_soc length: %d\n', length(dcir_soc_data.(channel).cyc200.charge_soc));
                if ~isempty(dcir_soc_data.(channel).cyc200.charge_soc)
                    fprintf('cyc200.charge_soc values: ');
                    fprintf('%.1f ', dcir_soc_data.(channel).cyc200.charge_soc);
                    fprintf('\n');
                    fprintf('cyc200.charge_dcir_60s values: ');
                    fprintf('%.6f ', dcir_soc_data.(channel).cyc200.charge_dcir_60s);
                    fprintf('\n');
                end
            end
        end
    else
        fprintf('Channel %s does NOT exist in dcir_soc_data\n', channel);
    end
    
    % 충전 DCIR vs SOC 플롯
    subplot(2,1,1);
    hold on;
    
    % RPT0 데이터 플롯
    if isfield(dcir_soc_data.(channel), 'cyc0') && isfield(dcir_soc_data.(channel).cyc0, 'charge_soc')
        if ~isempty(dcir_soc_data.(channel).cyc0.charge_soc)
            % NaN 값 제거
            valid_idx = ~isnan(dcir_soc_data.(channel).cyc0.charge_dcir_60s);
            if any(valid_idx)
                plot(dcir_soc_data.(channel).cyc0.charge_soc(valid_idx), dcir_soc_data.(channel).cyc0.charge_dcir_60s(valid_idx), ...
                    'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'DisplayName', 'RPT0 60s');
                fprintf('RPT0 Charge: %d valid points plotted\n', sum(valid_idx));
            end
        end
    end
    
    % RPT200 데이터 플롯
    if isfield(dcir_soc_data.(channel), 'cyc200') && isfield(dcir_soc_data.(channel).cyc200, 'charge_soc')
        if ~isempty(dcir_soc_data.(channel).cyc200.charge_soc)
            % NaN 값 제거
            valid_idx = ~isnan(dcir_soc_data.(channel).cyc200.charge_dcir_60s);
            if any(valid_idx)
                plot(dcir_soc_data.(channel).cyc200.charge_soc(valid_idx), dcir_soc_data.(channel).cyc200.charge_dcir_60s(valid_idx), ...
                    'o', 'Color', 'r', 'MarkerFaceColor', 'none', 'MarkerSize', 6, 'DisplayName', 'RPT200 60s');
                fprintf('RPT200 Charge: %d valid points plotted\n', sum(valid_idx));
            end
        end
    end
    
    % 방전 DCIR vs SOC 플롯
    subplot(2,1,2);
    hold on;
    
    % RPT0 데이터 플롯
    if isfield(dcir_soc_data.(channel), 'cyc0') && isfield(dcir_soc_data.(channel).cyc0, 'discharge_soc')
        if ~isempty(dcir_soc_data.(channel).cyc0.discharge_soc)
            % NaN 값 제거
            valid_idx = ~isnan(dcir_soc_data.(channel).cyc0.discharge_dcir_60s);
            if any(valid_idx)
                plot(dcir_soc_data.(channel).cyc0.discharge_soc(valid_idx), dcir_soc_data.(channel).cyc0.discharge_dcir_60s(valid_idx), ...
                    'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'DisplayName', 'RPT0 60s');
                fprintf('RPT0 Discharge: %d valid points plotted\n', sum(valid_idx));
            end
        end
    end
    
    % RPT200 데이터 플롯
    if isfield(dcir_soc_data.(channel), 'cyc200') && isfield(dcir_soc_data.(channel).cyc200, 'discharge_soc')
        if ~isempty(dcir_soc_data.(channel).cyc200.discharge_soc)
            % NaN 값 제거
            valid_idx = ~isnan(dcir_soc_data.(channel).cyc200.discharge_dcir_60s);
            if any(valid_idx)
                plot(dcir_soc_data.(channel).cyc200.discharge_soc(valid_idx), dcir_soc_data.(channel).cyc200.discharge_dcir_60s(valid_idx), ...
                    'o', 'Color', 'r', 'MarkerFaceColor', 'none', 'MarkerSize', 6, 'DisplayName', 'RPT200 60s');
                fprintf('RPT200 Discharge: %d valid points plotted\n', sum(valid_idx));
            end
        end
    end

    % Plot formatting
    subplot(2,1,1);
    xlabel('SOC [%]');
    ylabel('DCIR [Ω]');
    title(sprintf('%s - Charge DCIR vs SOC (60s DCIR)', channel));
    legend('Location', 'best');
    grid on;
    xlim([0 100]);

    subplot(2,1,2);
    xlabel('SOC [%]');
    ylabel('DCIR [Ω]');
    title(sprintf('%s - Discharge DCIR vs SOC (60s DCIR)', channel));
    legend('Location', 'best');
    grid on;
    xlim([0 100]);

    % Save figure - 보이는 상태로 저장
    figName = fullfile(dcir_socFolder, sprintf('%s_DCIR_SOC_fixedCap.fig', channel));
    set(gcf, 'Visible', 'on');  % figure를 보이게 설정
    savefig(gcf, figName);
    fprintf('Saved: %s\n', figName);
end

% Save DCIR-SOC data to mat file
fprintf('\n=== Saving DCIR-SOC data to MAT file ===\n');
matFileName = fullfile(savePath, 'DCIR_SOC_data_NEW.mat');
save(matFileName, 'dcir_soc_data');
fprintf('Saved: %s\n', matFileName);

fprintf('\n=== RPT Postprocessing (NEW VERSION) Completed ===\n'); 