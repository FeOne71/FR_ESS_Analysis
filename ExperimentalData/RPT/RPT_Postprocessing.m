%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Postprocessing - Generate and Save Figures by Category
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% File Directory
folderPath = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data\RPT';
savePath   = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\RPT_Postprocessing';

% Change save path to current directory to resolve path issues
currentDir = pwd;
savePath = fullfile(currentDir, 'RPT_Postprocessing_Results_0723');

channels = {'Ch9'}; % , 'Ch10', 'Ch11', 'Ch12', 'Ch13', 'Ch14', 'Ch15', 'Ch16'};
rpt_cycles = {'0cyc', '200cyc'};
rpt_cycles_fields = {'cyc0', 'cyc200'}; % For field names

Time      = [];
totalTime = [];

%% 1. Static Capacity
% fprintf('=== Static Capacity Processing ===\n');
% staticCapacityFolder = fullfile(savePath, 'StaticCapacity');
% if ~exist(staticCapacityFolder, 'dir')
%     mkdir(staticCapacityFolder);
% end
%
% for ch_idx = 1:length(channels)
%     channel = channels{ch_idx};
%     figure('Name', sprintf('Static Capacity - %s', channel), 'Position', [100 100 1200 800]);
%     set(gcf, 'Visible', 'off');
%
%     % Static Capacity conditions: (stepidx, cycleidx) = (1,1), (1,3) - charge, (3,1), (3,3) - discharge
%     step_cycle_conditions = {[1,1], [1,3], [3,1], [3,3]}; % charge, charge, discharge, discharge
%     condition_names = {'Charge RPT0', 'Charge RPT200', 'Discharge RPT0', 'Discharge RPT200'};
%     colors = {'b', 'b--', 'r', 'r--'};
%
%     hold on;
%     legend_entries = {};
%
%     for rpt_idx = 1:length(rpt_cycles)
%         rpt_cycle = rpt_cycles{rpt_idx};
%         filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
%         filepath = fullfile(folderPath, filename);
%
%         T = readtable(filepath);
%
%         % Charge (step 1)
%         for cycle_idx = [1, 3]
%             static_idx = (T{:,2} == 1) & (T{:,4} == cycle_idx);
%             if sum(static_idx) > 0
%                 capacity_data = T{static_idx, 9};
%                 voltage_data = T{static_idx, 8};
%                 final_capacity = capacity_data(end);
%
%                 if cycle_idx == 1
%                     if strcmp(rpt_cycle, '0cyc')
%                         plot(capacity_data, voltage_data, 'b-', 'LineWidth', 2);
%                         legend_entries{end+1} = sprintf('Charge RPT0 (Q=%.2f Ah)', final_capacity);
%                     else
%                         plot(capacity_data, voltage_data, 'b--', 'LineWidth', 2);
%                         legend_entries{end+1} = sprintf('Charge RPT200 (Q=%.2f Ah)', final_capacity);
%                     end
%                 else
%                     if strcmp(rpt_cycle, '0cyc')
%                         plot(capacity_data, voltage_data, 'r-', 'LineWidth', 2);
%                         legend_entries{end+1} = sprintf('Discharge RPT0 (Q=%.2f Ah)', final_capacity);
%                     else
%                         plot(capacity_data, voltage_data, 'r--', 'LineWidth', 2);
%                         legend_entries{end+1} = sprintf('Discharge RPT200 (Q=%.2f Ah)', final_capacity);
%                     end
%                 end
%             end
%         end
%
%         % Discharge (step 3)
%         for cycle_idx = [1, 3]
%             static_idx = (T{:,2} == 3) & (T{:,4} == cycle_idx);
%             if sum(static_idx) > 0
%                 capacity_data = T{static_idx, 9};
%                 voltage_data = T{static_idx, 8};
%                 final_capacity = capacity_data(end);
%
%                 if cycle_idx == 1
%                     if strcmp(rpt_cycle, '0cyc')
%                         plot(capacity_data, voltage_data, 'g-', 'LineWidth', 2);
%                         legend_entries{end+1} = sprintf('Discharge RPT0 (Q=%.2f Ah)', final_capacity);
%                     else
%                         plot(capacity_data, voltage_data, 'g--', 'LineWidth', 2);
%                         legend_entries{end+1} = sprintf('Discharge RPT200 (Q=%.2f Ah)', final_capacity);
%                     end
%                 else
%                     if strcmp(rpt_cycle, '0cyc')
%                         plot(capacity_data, voltage_data, 'm-', 'LineWidth', 2);
%                         legend_entries{end+1} = sprintf('Discharge RPT0 (Q=%.2f Ah)', final_capacity);
%                     else
%                         plot(capacity_data, voltage_data, 'm--', 'LineWidth', 2);
%                         legend_entries{end+1} = sprintf('Discharge RPT200 (Q=%.2f Ah)', final_capacity);
%                     end
%                 end
%             end
%         end
%     end
%
%     xlabel('Capacity [Ah]');
%     ylabel('Voltage [V]');
%     title(sprintf('%s - Static Capacity', channel));
%     legend(legend_entries, 'Location', 'best');
%     grid on;
%     xlim([0 70]);
%     ylim([2.5 4.5]);
%
%     % Save figure
%     figName = fullfile(staticCapacityFolder, sprintf('%s_StaticCapacity.fig', channel));
%     savefig(gcf, figName);
%     fprintf('Saved: %s\n', figName);
% end

%% 2. OCV
% fprintf('\n=== OCV Processing ===\n');
% ocvFolder = fullfile(savePath, 'OCV');
% if ~exist(ocvFolder, 'dir')
%     mkdir(ocvFolder);
% end
%
% % OCV conditions: charge (8,2), discharge (10,2)
% ocv_conditions = {'charge', 'discharge'};
% ocv_steps = [8, 10];
%
% % Structure to store OCV data for each channel
% channel_ocv_data = struct();
%
% for ch_idx = 1:length(channels)
%     channel = channels{ch_idx};
%
%     % Figure 1: Individual channel OCV
%     figure('Name', sprintf('OCV - %s', channel), 'Position', [100 100 1200 800]);
%     set(gcf, 'Visible', 'off');
%
%     for rpt_idx = 1:length(rpt_cycles)
%         rpt_cycle = rpt_cycles{rpt_idx};
%         filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
%         filepath = fullfile(folderPath, filename);
%
%         T = readtable(filepath);
%
%         subplot(2,1,rpt_idx);
%         hold on;
%
%         charge_ocv = [];
%         discharge_ocv = [];
%
%         for ocv_idx = 1:length(ocv_conditions)
%             step_idx = ocv_steps(ocv_idx);
%             ocv_data = (T{:,2} == step_idx) & (T{:,4} == 2);
%
%             capacity_data = T{ocv_data, 9};
%             voltage_data = T{ocv_data, 8};
%
%             % Divide capacity into 0~100 (101 points)
%             soc_data = linspace(0, 100, 101);
%
%             % Reverse discharge OCV so voltage increases with SOC
%             if strcmp(ocv_conditions{ocv_idx}, 'discharge')
%                 voltage_data = flipud(voltage_data);
%             end
%
%             % Interpolate voltage_data to match soc_data length
%             voltage_interp = interp1(linspace(0, 100, length(voltage_data)), voltage_data, soc_data, 'linear');
%
%             if strcmp(ocv_conditions{ocv_idx}, 'charge')
%                 charge_ocv = voltage_interp;
%                 plot(soc_data, voltage_interp, 'b-', 'LineWidth', 2, 'DisplayName', 'Charge');
%             else
%                 discharge_ocv = voltage_interp;
%                 plot(soc_data, voltage_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'Discharge');
%             end
%         end
%
%         % Calculate and display average charge/discharge OCV
%         if ~isempty(charge_ocv) && ~isempty(discharge_ocv)
%             avg_ocv = (charge_ocv + discharge_ocv) / 2;
%             plot(soc_data, avg_ocv, 'g-', 'LineWidth', 2, 'DisplayName', 'Average');
%
%             % Save data for dQdV calculation
%             if ~isfield(channel_ocv_data, channel)
%                 channel_ocv_data.(channel) = struct();
%             end
%
%             % Field name mapping
%             if strcmp(rpt_cycle, '0cyc')
%                 field_name = 'cyc0';
%             else
%                 field_name = 'cyc200';
%             end
%
%             channel_ocv_data.(channel).(field_name) = struct();
%             channel_ocv_data.(channel).(field_name).soc = soc_data;
%             channel_ocv_data.(channel).(field_name).avg_ocv = avg_ocv;
%
%             % Save final capacity value of OCV step
%             charge_mask = ocv_data & (T{:,2} == ocv_steps(1));
%             discharge_mask = ocv_data & (T{:,2} == ocv_steps(2));
%             charge_capacity = T{charge_mask, 9};
%             discharge_capacity = T{discharge_mask, 9};
%             if ~isempty(charge_capacity) && ~isempty(discharge_capacity)
%                 avg_capacity = (charge_capacity(end) + discharge_capacity(end)) / 2;
%                 channel_ocv_data.(channel).(field_name).capacity = avg_capacity;
%             end
%         end
%
%         xlabel('SOC [%]');
%         ylabel('Voltage [V]');
%         title(sprintf('%s - %s OCV', channel, rpt_cycle));
%         legend('Location', 'best');
%         grid on;
%         xlim([0 100]);
%     end
%
%     % Save figure
%     figName = fullfile(ocvFolder, sprintf('%s_OCV.fig', channel));
%     savefig(gcf, figName);
%     fprintf('Saved: %s\n', figName);
% end
%
% % Figure 2: Average OCV for each cycle separately
% figure('Name', 'Average OCV by Cycle', 'Position', [100 100 1200 800]);
% hold on;
%
% % RPT0 평균 OCV 계산
% all_ocv_data_rpt0 = [];
% for ch_idx = 1:length(channels)
%     channel = channels{ch_idx};
%     rpt_cycle = '0cyc';
%         filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
%         filepath = fullfile(folderPath, filename);
%
%         T = readtable(filepath);
%
%     % Combine charge and discharge OCV for RPT0
%         for ocv_idx = 1:length(ocv_conditions)
%             step_idx = ocv_steps(ocv_idx);
%             ocv_data = (T{:,2} == step_idx) & (T{:,4} == 2);
%
%             capacity_data = T{ocv_data, 9};
%             voltage_data = T{ocv_data, 8};
%
%         % Divide capacity into 0~100
%             soc_data = linspace(0, 100, length(capacity_data));
%
%         % Reverse discharge OCV so voltage increases with SOC
%             if strcmp(ocv_conditions{ocv_idx}, 'discharge')
%                 voltage_data = flipud(voltage_data);
%             end
%
%             % Interpolate to common SOC grid
%             soc_grid = 0:1:100;
%             voltage_interp = interp1(soc_data, voltage_data, soc_grid, 'linear');
%         all_ocv_data_rpt0 = [all_ocv_data_rpt0; voltage_interp];
%     end
% end
%
% % RPT200 평균 OCV 계산
% all_ocv_data_rpt200 = [];
% for ch_idx = 1:length(channels)
%     channel = channels{ch_idx};
%     rpt_cycle = '200cyc';
%     filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
%     filepath = fullfile(folderPath, filename);
%
%     T = readtable(filepath);
%
%     % Combine charge and discharge OCV for RPT200
%     for ocv_idx = 1:length(ocv_conditions)
%         step_idx = ocv_steps(ocv_idx);
%         ocv_data = (T{:,2} == step_idx) & (T{:,4} == 2);
%
%         capacity_data = T{ocv_data, 9};
%         voltage_data = T{ocv_data, 8};
%
%         % Divide capacity into 0~100
%         soc_data = linspace(0, 100, length(capacity_data));
%
%         % Reverse discharge OCV so voltage increases with SOC
%         if strcmp(ocv_conditions{ocv_idx}, 'discharge')
%             voltage_data = flipud(voltage_data);
%         end
%
%         % Interpolate to common SOC grid
%         soc_grid = 0:1:100;
%         voltage_interp = interp1(soc_data, voltage_data, soc_grid, 'linear');
%         all_ocv_data_rpt200 = [all_ocv_data_rpt200; voltage_interp];
%     end
% end
%
% % Calculate average OCV for each cycle
% avg_ocv_rpt0 = mean(all_ocv_data_rpt0, 1, 'omitnan');
% avg_ocv_rpt200 = mean(all_ocv_data_rpt200, 1, 'omitnan');
% soc_grid = 0:1:100;
%
% % Grid와 전압을 모두 오름차순으로 정렬
% [soc_grid_sorted, sort_idx] = sort(soc_grid);
% avg_ocv_rpt0_sorted = avg_ocv_rpt0(sort_idx);
% avg_ocv_rpt200_sorted = avg_ocv_rpt200(sort_idx);
%
% % Plot both RPT0 and RPT200 average OCV
% plot(soc_grid_sorted, avg_ocv_rpt0_sorted, 'b-', 'LineWidth', 3, 'DisplayName', 'RPT0 Average');
% plot(soc_grid_sorted, avg_ocv_rpt200_sorted, 'r-', 'LineWidth', 3, 'DisplayName', 'RPT200 Average');
% xlabel('SOC [%]');
% ylabel('Voltage [V]');
% title('Average OCV by Cycle (All Channels)');
% legend('Location', 'best');
% grid on;
% xlim([0 100]);
%
% % Save average OCV
% figName = fullfile(ocvFolder, 'Average_OCV_by_Cycle.fig');
% savefig(gcf, figName);
% set(gcf, 'Visible', 'off');
% fprintf('Saved: %s\n', figName);
%
% % Create OCV functions with sorted data
% OCV_func_rpt0 = @(soc_query) interp1(soc_grid_sorted, avg_ocv_rpt0_sorted, soc_query, 'linear');
% OCV_func_rpt200 = @(soc_query) interp1(soc_grid_sorted, avg_ocv_rpt200_sorted, soc_query, 'linear');

%% 3. DCIR vs SOC
fprintf('\n=== SOC-DCIR Processing ===\n');
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

% DCIR-SOC 데이터 저장을 위한 구조체
dcir_soc_data = struct();

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};

    figure('Name', sprintf('DCIR vs SOC - %s', channel), 'Position', [100 100 1200 800]);
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
                % SOC 계산을 위한 연속 시간축 생성 (한 번만)
                all_cycles_mask = (T.StepIndex == soc_charge_step);
                all_cycles_data = T(all_cycles_mask, :);
                all_cycles = unique(all_cycles_data.CycleIndex);
                all_cycles = sort(all_cycles);

                % 연속 시간축 생성
                continuous_time = [];
                continuous_current = [];
                time_offset = 0;

                for cycle_soc_idx = 1:length(all_cycles)
                    cycle_soc_idx_val = all_cycles(cycle_soc_idx);
                    cycle_soc_mask = all_cycles_data.CycleIndex == cycle_soc_idx_val;
                    cycle_soc_data = all_cycles_data(cycle_soc_mask, :);

                    if height(cycle_soc_data) > 0
                        % 현재 cycle의 시간 범위
                        cycle_start_time = cycle_soc_data.TotalTime(1);
                        cycle_end_time = cycle_soc_data.TotalTime(end);
                        cycle_duration = seconds(cycle_end_time - cycle_start_time);

                        % 연속 시간축에 추가
                        cycle_continuous_time = (0:0.1:cycle_duration) + time_offset;
                        continuous_time = [continuous_time, cycle_continuous_time];

                        % 전류 데이터 보간하여 추가
                        cycle_current_interp = interp1(seconds(cycle_soc_data.TotalTime - cycle_start_time), ...
                            cycle_soc_data.Current_A_, cycle_continuous_time - time_offset, 'linear');
                        continuous_current = [continuous_current, cycle_current_interp];

                        % 다음 cycle을 위한 시간 오프셋
                        time_offset = time_offset + cycle_duration;
                    end
                end

                % 누적 적분 계산 (cumtrapz 사용)
                cumulative_current = cumtrapz(continuous_time, continuous_current);
                total_current = cumulative_current(end);

                % 각 DCIR cycle별로 순서대로 처리
                for i = 1:length(charge_dcir_cycles)
                    cycle_idx = charge_dcir_cycles(i);

                    % 해당 cycle의 DCIR 데이터
                    cycle_dcir_mask = charge_dcir_data.CycleIndex == cycle_idx;
                    cycle_dcir = charge_dcir_data(cycle_dcir_mask, :);

                    if height(cycle_dcir) > 0
                        % DCIR 값 계산 - 특정 시점에서만 계산
                        % 1분 펄스에서 1s, 5s, 10s, 30s 시점의 DCIR 추출
                        dcir_time_points = [1, 5, 10, 30]; % seconds
                        dcir_values = [];

                        for time_idx = 1:length(dcir_time_points)
                            target_time = dcir_time_points(time_idx);

                            % 해당 시점에 가장 가까운 데이터 찾기
                            pulse_times = seconds(cycle_dcir.TotalTime - cycle_dcir.TotalTime(1));
                            [~, closest_idx] = min(abs(pulse_times - target_time));

                            if closest_idx <= height(cycle_dcir)
                                voltage = cycle_dcir.Voltage_V_(closest_idx);
                                current = cycle_dcir.Current_A_(closest_idx);
                                dcir_value = abs(voltage / current);
                                dcir_values = [dcir_values; dcir_value];
                            end
                        end

                        % SOC 계산 - 순서에 따라
                        if i == 1
                            % 첫 번째 DCIR → SOC 0%
                            soc_percent = 0;
                        else
                            % 두 번째부터는 실제 적분 계산
                            % 해당 cycle의 SOC 구간 시작 시간 찾기
                            cycle_soc_mask = charge_soc_data.CycleIndex == cycle_idx;
                            cycle_soc = charge_soc_data(cycle_soc_mask, :);

                            if height(cycle_soc) > 0
                                % 해당 cycle의 SOC 구간 시작 시간
                                cycle_soc_start_time = cycle_soc.TotalTime(1);

                                % 연속 시간축에서 해당 cycle 시작 위치 찾기
                                dcir_continuous_time = 0;
                                for cycle_soc_idx = 1:length(all_cycles)
                                    cycle_soc_idx_val = all_cycles(cycle_soc_idx);
                                    cycle_soc_mask = all_cycles_data.CycleIndex == cycle_soc_idx_val;
                                    cycle_soc_data = all_cycles_data(cycle_soc_mask, :);

                                    if height(cycle_soc_data) > 0
                                        cycle_start_time = cycle_soc_data.TotalTime(1);
                                        cycle_end_time = cycle_soc_data.TotalTime(end);
                                        cycle_duration = seconds(cycle_end_time - cycle_start_time);

                                        if cycle_soc_idx_val == cycle_idx
                                            % 해당 cycle의 시작 시간을 연속 시간축에서 찾기
                                            break;
                                        else
                                            dcir_continuous_time = dcir_continuous_time + cycle_duration;
                                        end
                                    end
                                end

                                % 해당 cycle 시작 시점까지의 적분 (trapz 사용)
                                [~, cycle_start_idx] = min(abs(continuous_time - dcir_continuous_time));
                                cycle_start_idx = min(cycle_start_idx, length(cumulative_current));

                                                                % 부분 적분 (trapz 사용)
                                partial_current = trapz(continuous_time(1:cycle_start_idx), continuous_current(1:cycle_start_idx));
                                
                                % 디버깅 출력
                                fprintf('  RPT%s - Cycle %d: Partial=%.6f, Total=%.6f, Ratio=%.3f%%\n', ...
                                    rpt_cycle, cycle_idx, partial_current, total_current, (partial_current/total_current)*100);
                                
                                % 충전 SOC = ∫(0→t) I dt / ∫(0→end) I dt
                                soc_percent = (partial_current / total_current) * 100;
                            else
                                soc_percent = 0;
                            end
                        end

                        fprintf('Charge DCIR #%d (Cycle %d): SOC = %.1f%%, DCIR values = %d\n', ...
                            i, cycle_idx, soc_percent, length(dcir_values));

                        % DCIR 값들을 같은 SOC로 플롯
                        soc_array = repmat(soc_percent, size(dcir_values));

                        subplot(2,1,1);
                        if strcmp(rpt_cycle, '0cyc')
                            plot(soc_array, dcir_values, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'DisplayName', sprintf('RPT0 DCIR #%d', i));
                        else
                            plot(soc_array, dcir_values, 'bo', 'MarkerFaceColor', 'none', 'MarkerSize', 5, 'DisplayName', sprintf('RPT200 DCIR #%d', i));
                        end
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
                % SOC 계산을 위한 연속 시간축 생성 (한 번만)
                all_cycles_mask = (T.StepIndex == soc_discharge_step);
                all_cycles_data = T(all_cycles_mask, :);
                all_cycles = unique(all_cycles_data.CycleIndex);
                all_cycles = sort(all_cycles);

                % 연속 시간축 생성
                continuous_time = [];
                continuous_current = [];
                time_offset = 0;

                for cycle_soc_idx = 1:length(all_cycles)
                    cycle_soc_idx_val = all_cycles(cycle_soc_idx);
                    cycle_soc_mask = all_cycles_data.CycleIndex == cycle_soc_idx_val;
                    cycle_soc_data = all_cycles_data(cycle_soc_mask, :);

                    if height(cycle_soc_data) > 0
                        % 현재 cycle의 시간 범위
                        cycle_start_time = cycle_soc_data.TotalTime(1);
                        cycle_end_time = cycle_soc_data.TotalTime(end);
                        cycle_duration = seconds(cycle_end_time - cycle_start_time);

                        % 연속 시간축에 추가
                        cycle_continuous_time = (0:0.1:cycle_duration) + time_offset;
                        continuous_time = [continuous_time, cycle_continuous_time];

                        % 전류 데이터 보간하여 추가
                        cycle_current_interp = interp1(seconds(cycle_soc_data.TotalTime - cycle_start_time), ...
                            cycle_soc_data.Current_A_, cycle_continuous_time - time_offset, 'linear');
                        continuous_current = [continuous_current, cycle_current_interp];

                        % 다음 cycle을 위한 시간 오프셋
                        time_offset = time_offset + cycle_duration;
                    end
                end

                % 누적 적분 계산 (cumtrapz 사용)
                cumulative_current = cumtrapz(continuous_time, continuous_current);
                total_current = cumulative_current(end);

                % 각 DCIR cycle별로 순서대로 처리
                for i = 1:length(discharge_dcir_cycles)
                    cycle_idx = discharge_dcir_cycles(i);

                    % 해당 cycle의 DCIR 데이터
                    cycle_dcir_mask = discharge_dcir_data.CycleIndex == cycle_idx;
                    cycle_dcir = discharge_dcir_data(cycle_dcir_mask, :);

                    if height(cycle_dcir) > 0
                        % DCIR 값 계산 - 특정 시점에서만 계산
                        % 1분 펄스에서 1s, 5s, 10s, 30s 시점의 DCIR 추출
                        dcir_time_points = [1, 5, 10, 30]; % seconds
                        dcir_values = [];

                        for time_idx = 1:length(dcir_time_points)
                            target_time = dcir_time_points(time_idx);

                            % 해당 시점에 가장 가까운 데이터 찾기
                            pulse_times = seconds(cycle_dcir.TotalTime - cycle_dcir.TotalTime(1));
                            [~, closest_idx] = min(abs(pulse_times - target_time));

                            if closest_idx <= height(cycle_dcir)
                                voltage = cycle_dcir.Voltage_V_(closest_idx);
                                current = cycle_dcir.Current_A_(closest_idx);
                                dcir_value = abs(voltage / current);
                                dcir_values = [dcir_values; dcir_value];
                            end
                        end

                        % SOC 계산 - 순서에 따라
                        if i == 1
                            % 첫 번째 DCIR → SOC 100%
                            soc_percent = 100;
                        else
                            % 두 번째부터는 실제 적분 계산
                            % 해당 cycle의 SOC 구간 시작 시간 찾기
                            cycle_soc_mask = discharge_soc_data.CycleIndex == cycle_idx;
                            cycle_soc = discharge_soc_data(cycle_soc_mask, :);

                            if height(cycle_soc) > 0
                                % 해당 cycle의 SOC 구간 시작 시간
                                cycle_soc_start_time = cycle_soc.TotalTime(1);

                                % 연속 시간축에서 해당 cycle 시작 위치 찾기
                                dcir_continuous_time = 0;
                                for cycle_soc_idx = 1:length(all_cycles)
                                    cycle_soc_idx_val = all_cycles(cycle_soc_idx);
                                    cycle_soc_mask = all_cycles_data.CycleIndex == cycle_soc_idx_val;
                                    cycle_soc_data = all_cycles_data(cycle_soc_mask, :);

                                    if height(cycle_soc_data) > 0
                                        cycle_start_time = cycle_soc_data.TotalTime(1);
                                        cycle_end_time = cycle_soc_data.TotalTime(end);
                                        cycle_duration = seconds(cycle_end_time - cycle_start_time);

                                        if cycle_soc_idx_val == cycle_idx
                                            % 해당 cycle의 시작 시간을 연속 시간축에서 찾기
                                            break;
                                        else
                                            dcir_continuous_time = dcir_continuous_time + cycle_duration;
                                        end
                                    end
                                end

                                % 해당 cycle 시작 시점까지의 적분 (trapz 사용)
                                [~, cycle_start_idx] = min(abs(continuous_time - dcir_continuous_time));
                                cycle_start_idx = min(cycle_start_idx, length(cumulative_current));

                                                                % 부분 적분 (trapz 사용)
                                partial_current = trapz(continuous_time(1:cycle_start_idx), continuous_current(1:cycle_start_idx));
                                
                                % 디버깅 출력
                                fprintf('  RPT%s - Cycle %d: Partial=%.6f, Total=%.6f, Ratio=%.3f%%, SOC=%.1f%%\n', ...
                                    rpt_cycle, cycle_idx, partial_current, total_current, (partial_current/total_current)*100, 100-(partial_current/total_current)*100);
                                
                                % 방전 SOC = 1 - ∫(0→t) I dt / ∫(0→end) I dt
                                soc_percent = 100 - (partial_current / total_current) * 100;
                            else
                                soc_percent = 100;
                            end
                        end

                        fprintf('Discharge DCIR #%d (Cycle %d): SOC = %.1f%%, DCIR values = %d\n', ...
                            i, cycle_idx, soc_percent, length(dcir_values));

                        % DCIR 값들을 같은 SOC로 플롯
                        soc_array = repmat(soc_percent, size(dcir_values));

                        subplot(2,1,2);
                        if strcmp(rpt_cycle, '0cyc')
                            plot(soc_array, dcir_values, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', sprintf('RPT0 DCIR #%d', i));
                        else
                            plot(soc_array, dcir_values, 'ro', 'MarkerFaceColor', 'none', 'MarkerSize', 5, 'DisplayName', sprintf('RPT200 DCIR #%d', i));
                        end
                    end
                end
            end
        end
    end

    % Plot formatting
    subplot(2,1,1);
    xlabel('SOC [%]');
    ylabel('DCIR [Ω]');
    title(sprintf('%s - Charge DCIR vs SOC', channel));
    % legend('Location', 'best');
    grid on;
    xlim([0 100]);

    subplot(2,1,2);
    xlabel('SOC [%]');
    ylabel('DCIR [Ω]');
    title(sprintf('%s - Discharge DCIR vs SOC', channel));
    % legend('Location', 'best');
    grid on;
    xlim([0 100]);

    % Save figure
    figName = fullfile(dcir_socFolder, sprintf('%s_DCIR_SOC.fig', channel));
    savefig(gcf, figName);
    set(gcf, 'Visible', 'on');
    fprintf('Saved: %s\n', figName);
end

%% 4. C-rate
% fprintf('\n=== C-rate Processing ===\n');
% crateFolder = fullfile(savePath, 'C_rate');
% if ~exist(crateFolder, 'dir')
%     mkdir(crateFolder);
% end
%
% % C-rate conditions - 각 채널별 정확한 cycle idx
% crate_charge_0cyc = {[28,30], [32,30], [36,30], [40,30], [44,30]}; % 0cyc 충전
% crate_discharge_0cyc = {[48,32], [52,32], [56,32], [60,32], [64,32]}; % 0cyc 방전
%
% % 200cyc는 대부분 30, 일부 채널만 다름
% crate_charge_200cyc = {[28,30], [32,30], [36,30], [40,30], [44,30]}; % 200cyc 충전 (기본)
% crate_discharge_200cyc = {[48,30], [52,30], [56,30], [60,30], [64,30]}; % 200cyc 방전 (기본)
%
% % 특별한 채널들 처리
% ch9_200cyc_charge = {[28,31], [32,31], [36,31], [40,31], [44,31]};
% ch9_200cyc_discharge = {[48,31], [52,31], [56,31], [60,31], [64,31]};
%
% ch15_200cyc_charge = {[28,31], [32,31], [36,31], [40,31], [44,31]};
% ch15_200cyc_discharge = {[48,31], [52,31], [56,31], [60,31], [64,31]};
%
% ch16_0cyc_charge   = {[28,32], [32,32], [36,32], [40,32], [44,32]};
% ch16_0cyc_discharge = {[48,32], [52,32], [56,32], [60,32], [64,32]};
%
%
% ch16_200cyc_discharge = {[48,30], [52,30], [56,30], [60,30], [64,30]};
%
% crate_names = {'0.1C', '0.5C', '1C', '2C', '3C'};
%
% for ch_idx = 1:length(channels)
%     channel = channels{ch_idx};
%
%     figure('Name', sprintf('C-rate - %s', channel), 'Position', [100 100 1200 800]);
%     set(gcf, 'Visible', 'off');
%
%     for rpt_idx = 1:length(rpt_cycles)
%         rpt_cycle = rpt_cycles{rpt_idx};
%         filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
%         filepath = fullfile(folderPath, filename);
%
%         T = readtable(filepath);
%
%         subplot(2,1,rpt_idx);
%         hold on;
%
%         % 0cyc와 200cyc에 따라 다른 cycle idx 사용
%         if strcmp(rpt_cycle, '0cyc')
%             charge_conditions = crate_charge_0cyc;
%             discharge_conditions = crate_discharge_0cyc;
%         else % 200cyc
%             if strcmp(channel, 'Ch9')
%                 charge_conditions = ch9_200cyc_charge;
%                 discharge_conditions = ch9_200cyc_discharge;
%             elseif strcmp(channel, 'Ch15')
%                 charge_conditions = ch15_200cyc_charge;
%                 discharge_conditions = ch15_200cyc_discharge;
%             else
%                 charge_conditions = crate_charge_200cyc;
%                 discharge_conditions = crate_discharge_200cyc;
%             end
%         end
%
%         % Charge C-rate
%         for crate_idx = 1:length(charge_conditions)
%             step_cycle = charge_conditions{crate_idx};
%             crate_idx_data = (T{:,2} == step_cycle(1)) & (T{:,4} == step_cycle(2));
%
%             if sum(crate_idx_data) > 0
%                 capacity_data = T{crate_idx_data, 9};
%                 voltage_data = T{crate_idx_data, 8};
%                 plot(capacity_data, voltage_data, 'LineWidth', 2, 'DisplayName', sprintf('Charge %s', crate_names{crate_idx}));
%             end
%         end
%
%         % Discharge C-rate
%         for crate_idx = 1:length(discharge_conditions)
%             step_cycle = discharge_conditions{crate_idx};
%             crate_idx_data = (T{:,2} == step_cycle(1)) & (T{:,4} == step_cycle(2));
%
%             if sum(crate_idx_data) > 0
%                 capacity_data = T{crate_idx_data, 9};
%                 voltage_data = T{crate_idx_data, 8};
%                 plot(capacity_data, voltage_data, '--', 'LineWidth', 2, 'DisplayName', sprintf('    Discharge %s', crate_names{crate_idx}));
%             end
%         end
%
%         xlabel('Capacity [Ah]');
%         ylabel('Voltage [V]');
%         title(sprintf('%s - %s C-rate', channel, rpt_cycle));
%         legend('Location', 'best');
%         grid on;
%     end
%
%     % Save figure
%     figName = fullfile(crateFolder, sprintf('%s_Crate.fig', channel));
%     savefig(gcf, figName);
%     set(gcf, 'Visible', 'off');
%     fprintf('Saved: %s\n', figName);
% % end

%% 5. dQdV Processing
% fprintf('\n=== dQdV Processing ===\n');
% dqdvFolder = fullfile(savePath, 'dQdV');
% if ~exist(dqdvFolder, 'dir')
%     mkdir(dqdvFolder);
% end
%
% for ch_idx = 1:length(channels)
%     channel = channels{ch_idx};
%
%     figure('Name', sprintf('dQdV - %s', channel), 'Position', [100 100 800 600]);
%     hold on;
%
%     for rpt_idx = 1:length(rpt_cycles)
%         rpt_cycle = rpt_cycles{rpt_idx};
%
%         % OCV 섹션에서 계산한 평균 OCV 데이터 사용
%         % 필드명 매핑
%         if strcmp(rpt_cycle, '0cyc')
%             field_name = 'cyc0';
%         else
%             field_name = 'cyc200';
%         end
%
%         if isfield(channel_ocv_data, channel) && isfield(channel_ocv_data.(channel), field_name)
%             soc_data = channel_ocv_data.(channel).(field_name).soc;
%             avg_ocv = channel_ocv_data.(channel).(field_name).avg_ocv;
%             capacity = channel_ocv_data.(channel).(field_name).capacity;
%
%             % SOC를 실제 용량으로 변환
%             capacity_data = soc_data * capacity / 100; % SOC%를 실제 용량[Ah]로 변환
%
%             % dQdV = dQ/dV 계산
%             dqdv = diff(capacity_data) ./ diff(avg_ocv) * 1000; % mAh/V 단위
%             voltage_center = (avg_ocv(1:end-1) + avg_ocv(2:end)) / 2; % 중간값
%
%             if strcmp(rpt_cycle, '0cyc')
%                 plot(voltage_center, dqdv, 'b-', 'LineWidth', 2, 'DisplayName', '0cyc');
%             else
%                 plot(voltage_center, dqdv, 'r-', 'LineWidth', 2, 'DisplayName', '200cyc');
%             end
%         end
%     end
%
%     xlabel('Voltage [V]');
%     ylabel('dQdV [mAh/V]');
%     title(sprintf('%s - dQdV Comparison', channel));
%     legend('Location', 'best');
%     grid on;
%
%     % Save figure
%     figName = fullfile(dqdvFolder, sprintf('%s_dQdV.fig', channel));
%     savefig(gcf, figName);
%     set(gcf, 'Visible', 'on');
%     fprintf('Saved: %s\n', figName);
% end

%% Save all data to mat file
% fprintf('\n=== Saving to MAT file ===\n');
% matFileName = fullfile(savePath, 'RPT_processed_data.mat');
%
% % Create structure to save all processed data
% RPT_data = struct();
% RPT_data.channels = channels;
% RPT_data.rpt_cycles = rpt_cycles;
% RPT_data.ocv_function = OCV_func;
% RPT_data.processing_date = datestr(now);
%
% save(matFileName, 'RPT_data');
% fprintf('Saved: %s\n', matFileName);
%
% fprintf('\n=== RPT Postprocessing Completed ===\n');