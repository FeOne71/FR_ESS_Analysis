%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Postprocessing - Generate and Save Figures by Category
% UPDATED VERSION: New DCIR calculation logic applied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% File Directory
folderPath = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data\RPT';
savePath   = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\RPT_Postprocessing';

% Change save path to current directory to resolve path issues
currentDir = pwd;
savePath = fullfile(currentDir, 'RPT_Postprocessing_Results_calcCap');

channels = {'Ch9'}; %, 'Ch10', 'Ch11', 'Ch12', 'Ch13', 'Ch14', 'Ch15', 'Ch16'};
rpt_cycles = {'0cyc', '200cyc'};
rpt_cycles_fields = {'cyc0', 'cyc200'}; % For field names

Time      = [];
totalTime = [];

%% 3. DCIR vs SOC
fprintf('\n=== SOC-DCIR Processing (UPDATED) ===\n');
dcir_socFolder = fullfile(savePath, 'DCIR');
if ~exist(dcir_socFolder, 'dir')
    mkdir(dcir_socFolder);
end

% DCIR step conditions5
dcir_charge_step = 14;   % 충전 DCIR 펄스
dcir_discharge_step = 21; % 방전 DCIR 펄스

% SOC calculation steps
soc_charge_step = 16;     % 충전 SOC 적산 구간
soc_discharge_step = 23;  % 방전 SOC 적산 구간

% DCIR-SOC 데이터 저장을 위한 구조체
dcir_soc_data = struct();

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};

    figure('Name', sprintf('DCIR vs SOC - %s (UPDATED)', channel), 'Position', [100 100 1200 800]);
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
                % 1. StepIdx=14 (DCIR)에서 cycle 2부터 있는 인덱스까지 적분값 계산
                dcir_cycles_mask = (T.StepIndex == dcir_charge_step) & (T.CycleIndex >= 2);
                dcir_cycles_data = T(dcir_cycles_mask, :);
                dcir_cycles = unique(dcir_cycles_data.CycleIndex);
                dcir_cycles = sort(dcir_cycles);
                
                % StepIdx=14의 각 cycle별 적분값 계산
                dcir_cycle_integrals = [];
                for cycle_idx = 1:length(dcir_cycles)
                    current_cycle = dcir_cycles(cycle_idx);
                    cycle_mask = (T.StepIndex == dcir_charge_step) & (T.CycleIndex == current_cycle);
                    cycle_data = T(cycle_mask, :);
                    
                    if height(cycle_data) > 0
                        % 각 cycle의 개별 적분
                        cycle_time = seconds(cycle_data.TotalTime - cycle_data.TotalTime(1));
                        cycle_current = cycle_data.Current_A_;
                        cycle_integral = trapz(cycle_time, cycle_current);
                        dcir_cycle_integrals = [dcir_cycle_integrals; cycle_integral];
                    end
                end
                
                % 2. StepIdx=16에서 모든 cycle의 적분값 계산
                step16_mask = (T.StepIndex == soc_charge_step);
                step16_data = T(step16_mask, :);
                step16_cycles = unique(step16_data.CycleIndex);
                step16_cycles = sort(step16_cycles);
                
                step16_cycle_integrals = [];
                for cycle_idx = 1:length(step16_cycles)
                    current_cycle = step16_cycles(cycle_idx);
                    cycle_mask = (T.StepIndex == soc_charge_step) & (T.CycleIndex == current_cycle);
                    cycle_data = T(cycle_mask, :);
                    
                    if height(cycle_data) > 0
                        % 각 cycle의 개별 적분
                        cycle_time = seconds(cycle_data.TotalTime - cycle_data.TotalTime(1));
                        cycle_current = cycle_data.Current_A_;
                        cycle_integral = trapz(cycle_time, cycle_current);
                        step16_cycle_integrals = [step16_cycle_integrals; cycle_integral];
                    end
                end
                
                % 3. StepIdx=19에서 제일 마지막 시점의 capacity
                step19_mask = (T.StepIndex == 19);
                step19_data = T(step19_mask, :);
                if height(step19_data) > 0
                    final_capacity = step19_data.Capacity_Ah_(end);  % 제일 마지막 capacity
                else
                    final_capacity = 0;
                end
                
                % 전체 capacity 값 = StepIdx=14 capacity 합 + StepIdx=16 capacity 합 + StepIdx=19 마지막 capacity
                step14_capacity_total = 0;
                step16_capacity_total = 0;
                
                % StepIdx=14 capacity 합 계산
                for cyc = 2:16
                    cyc_mask = (T.StepIndex == dcir_charge_step) & (T.CycleIndex == cyc);
                    cyc_data = T(cyc_mask, :);
                    if height(cyc_data) > 0
                        step14_capacity_total = step14_capacity_total + cyc_data.Capacity_Ah_(end);
                    end
                end
                
                % StepIdx=16 capacity 합 계산
                for cyc = 2:16
                    cyc_mask = (T.StepIndex == soc_charge_step) & (T.CycleIndex == cyc);
                    cyc_data = T(cyc_mask, :);
                    if height(cyc_data) > 0
                        step16_capacity_total = step16_capacity_total + cyc_data.Capacity_Ah_(end);
                    end
                end
                
                total_integral = step14_capacity_total + step16_capacity_total + final_capacity;
                
                % 디버깅: 각 StepIdx와 cycle별 capacity 값들
                fprintf('  Debug - StepIdx=14 cycle별 capacity (Ah):\n');
                for cyc = 2:16
                    cyc_mask = (T.StepIndex == dcir_charge_step) & (T.CycleIndex == cyc);
                    cyc_data = T(cyc_mask, :);
                    if height(cyc_data) > 0
                        fprintf('    Cycle %d: %.6f Ah\n', cyc, cyc_data.Capacity_Ah_(end));
                    end
                end
                
                fprintf('  Debug - StepIdx=16 cycle별 capacity (Ah):\n');
                for cyc = 2:16
                    cyc_mask = (T.StepIndex == soc_charge_step) & (T.CycleIndex == cyc);
                    cyc_data = T(cyc_mask, :);
                    if height(cyc_data) > 0
                        fprintf('    Cycle %d: %.6f Ah\n', cyc, cyc_data.Capacity_Ah_(end));
                    end
                end
                
                fprintf('  Debug - StepIdx=19 final capacity: %.6f Ah\n', final_capacity);
                fprintf('  Debug - StepIdx=14 capacity total: %.6f Ah\n', step14_capacity_total);
                fprintf('  Debug - StepIdx=16 capacity total: %.6f Ah\n', step16_capacity_total);
                fprintf('  Debug - Total capacity: %.6f Ah\n', total_integral);

                % 각 DCIR cycle별로 순서대로 처리
                for i = 1:length(charge_dcir_cycles)
                    cycle_idx = charge_dcir_cycles(i);

                    % 해당 cycle의 DCIR 데이터
                    cycle_dcir_mask = charge_dcir_data.CycleIndex == cycle_idx;
                    cycle_dcir = charge_dcir_data(cycle_dcir_mask, :);

                    if height(cycle_dcir) > 0
                        % DCIR 값 계산 - 시간별로 정확히 저장
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
                            exact_idx = find(relative_times == target_time, 1);

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
                                if target_time == 0.1
                                    dcir_values_01s = dcir_value;
                                elseif target_time == 5
                                    dcir_values_5s = dcir_value;
                                elseif target_time == 10
                                    dcir_values_10s = dcir_value;
                                elseif target_time == 30
                                    dcir_values_30s = dcir_value;
                                elseif target_time == 60
                                    dcir_values_60s = dcir_value;
                                end
                            end
                        end

                        % SOC 계산 - 수정된 방식
                        if i == 1
                            % 첫 번째 DCIR → SOC 0%
                            soc_percent = 0;
                        else
                            % SOC 계산 - 누적 방식
                            if i == 1
                                soc_percent = 0;  % SOC1 = 0
                            else
                                % 현재 cycle의 StepIdx=14와 StepIdx=16 capacity 값
                                current_cycle = charge_dcir_cycles(i);
                                
                                % StepIdx=14, 현재 cycle의 capacity
                                step14_current_mask = (T.StepIndex == dcir_charge_step) & (T.CycleIndex == current_cycle);
                                step14_current_data = T(step14_current_mask, :);
                                if height(step14_current_data) > 0
                                    step14_current_capacity = step14_current_data.Capacity_Ah_(end);
                                else
                                    step14_current_capacity = 0;
                                end
                                
                                % StepIdx=16, 현재 cycle의 capacity
                                step16_current_mask = (T.StepIndex == soc_charge_step) & (T.CycleIndex == current_cycle);
                                step16_current_data = T(step16_current_mask, :);
                                if height(step16_current_data) > 0
                                    step16_current_capacity = step16_current_data.Capacity_Ah_(end);
                                else
                                    step16_current_capacity = 0;
                                end
                                
                                % 분자 = 현재 cycle의 StepIdx=14 + StepIdx=16 capacity
                                numerator = step14_current_capacity + step16_current_capacity;
                                
                                % 분모 = StepIdx=14 capacity 합 + StepIdx=16 capacity 합 + StepIdx=19
                                denominator = step14_capacity_total + step16_capacity_total + final_capacity;
                                
                                % 현재 cycle의 SOC 증가량
                                current_soc_increment = (numerator / denominator) * 100;
                                
                                % 누적 SOC 계산
                                if i == 2
                                    soc_percent = current_soc_increment;  % SOC2 = 첫 번째 증가량
                                else
                                    % 이전 SOC + 현재 증가량
                                    soc_percent = soc_percent + current_soc_increment;
                                end
                            end
                            
                            % 디버깅 출력
                            fprintf('  RPT%s - Cycle %d: StepIdx=14=%.6f, StepIdx=16=%.6f, Total=%.6f, SOC=%.1f%%\n', ...
                                rpt_cycle, cycle_idx, step14_current_capacity, step16_current_capacity, denominator, soc_percent);
                            
                            % SOC는 이미 계산됨
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
                        dcir_soc_data.(channel).(rpt_field).charge_soc = [dcir_soc_data.(channel).(rpt_field).charge_soc; soc_percent];
                        
                        % 시간별로 정확히 저장된 DCIR 값들 저장
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_01s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_01s; dcir_values_01s];
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_5s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_5s; dcir_values_5s];
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_10s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_10s; dcir_values_10s];
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_30s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_30s; dcir_values_30s];
                        dcir_soc_data.(channel).(rpt_field).charge_dcir_60s = [dcir_soc_data.(channel).(rpt_field).charge_dcir_60s; dcir_values_60s];
                        
                        dcir_soc_data.(channel).(rpt_field).charge_cycles = [dcir_soc_data.(channel).(rpt_field).charge_cycles; cycle_idx];
                    end
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
                % 원본 데이터 직접 사용 (시간 간격이 동일하므로 보간 불필요)
                all_cycles_mask = (T.StepIndex == soc_discharge_step);
                all_cycles_data = T(all_cycles_mask, :);
                all_cycles = unique(all_cycles_data.CycleIndex);
                all_cycles = sort(all_cycles);

                % 각 cycle의 개별 적분값 계산
                all_cycle_integrals = [];
                for cycle_idx = 1:length(all_cycles)
                    current_cycle = all_cycles(cycle_idx);
                    cycle_mask = all_cycles_data.CycleIndex == current_cycle;
                    cycle_data = all_cycles_data(cycle_mask, :);
                    
                    if height(cycle_data) > 0
                        % 각 cycle의 개별 적분
                        cycle_time = seconds(cycle_data.TotalTime - cycle_data.TotalTime(1));
                        cycle_current = cycle_data.Current_A_;
                        cycle_integral = trapz(cycle_time, cycle_current);
                        all_cycle_integrals = [all_cycle_integrals; cycle_integral];
                    end
                end
                
                % 전체 적분값
                total_integral = sum(all_cycle_integrals);

                % 각 DCIR cycle별로 순서대로 처리
                for i = 1:length(discharge_dcir_cycles)
                    cycle_idx = discharge_dcir_cycles(i);

                    % 해당 cycle의 DCIR 데이터
                    cycle_dcir_mask = discharge_dcir_data.CycleIndex == cycle_idx;
                    cycle_dcir = discharge_dcir_data(cycle_dcir_mask, :);

                    if height(cycle_dcir) > 0
                        % DCIR 값 계산 - 시간별로 정확히 저장
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
                            exact_idx = find(relative_times == target_time, 1);

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
                                if target_time == 0.1
                                    dcir_values_01s = dcir_value;
                                elseif target_time == 5
                                    dcir_values_5s = dcir_value;
                                elseif target_time == 10
                                    dcir_values_10s = dcir_value;
                                elseif target_time == 30
                                    dcir_values_30s = dcir_value;
                                elseif target_time == 60
                                    dcir_values_60s = dcir_value;
                                end
                            end
                        end

                        % SOC 계산 - 각 cycle 개별 적분 후 합산
                        if i == 1
                            % 첫 번째 DCIR → SOC 100%
                            soc_percent = 100;
                        else
                            % DCIR cycle 순서에 맞춰 누적 적분 계산
                            cumulative_integral = sum(all_cycle_integrals(1:i-1));
                            
                            % 디버깅 출력
                            fprintf('  RPT%s - Cycle %d: Cumulative=%.6f, Total=%.6f, SOC=%.1f%%\n', ...
                                rpt_cycle, cycle_idx, cumulative_integral, total_integral, 100-(cumulative_integral / total_integral) * 100);
                            
                            % 방전 SOC = 100% - (누적 적분 / 전체 적분)
                            soc_percent = 100 - (cumulative_integral / total_integral) * 100;
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
                        dcir_soc_data.(channel).(rpt_field).discharge_soc = [dcir_soc_data.(channel).(rpt_field).discharge_soc; soc_percent];
                        
                        % 시간별로 정확히 저장된 DCIR 값들 저장
                        dcir_soc_data.(channel).(rpt_field).discharge_dcir_01s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_01s; dcir_values_01s];
                        dcir_soc_data.(channel).(rpt_field).discharge_dcir_5s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_5s; dcir_values_5s];
                        dcir_soc_data.(channel).(rpt_field).discharge_dcir_10s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_10s; dcir_values_10s];
                        dcir_soc_data.(channel).(rpt_field).discharge_dcir_30s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_30s; dcir_values_30s];
                        dcir_soc_data.(channel).(rpt_field).discharge_dcir_60s = [dcir_soc_data.(channel).(rpt_field).discharge_dcir_60s; dcir_values_60s];
                        
                        dcir_soc_data.(channel).(rpt_field).discharge_cycles = [dcir_soc_data.(channel).(rpt_field).discharge_cycles; cycle_idx];
                    end
                end
            end
        end
    end

    % === 저장된 데이터를 사용한 시각화 ===
    fprintf('\n=== Plotting from saved data ===\n');
    
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
    figName = fullfile(dcir_socFolder, sprintf('%s_DCIR_SOC_calcCap.fig', channel));
    set(gcf, 'Visible', 'on');  % figure를 보이게 설정
    savefig(gcf, figName);
    fprintf('Saved: %s\n', figName);
end

% Save DCIR-SOC data to mat file
fprintf('\n=== Saving DCIR-SOC data to MAT file ===\n');
matFileName = fullfile(savePath, 'DCIR_SOC_data_UPDATED.mat');
save(matFileName, 'dcir_soc_data');
fprintf('Saved: %s\n', matFileName);

fprintf('\n=== RPT Postprocessing (UPDATED) Completed ===\n'); 
fprintf('\n=== RPT Postprocessing (UPDATED) Completed ===\n'); 