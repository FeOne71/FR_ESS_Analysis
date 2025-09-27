%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Postprocessing - Generate and Save Figures by Category
% Chg/Dch OCV > Avg OCV > Iterate 8Ch 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% File Directory
ExperimentalDataPath = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data\RPT';
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';
if ~exist(saveDir,'dir'); mkdir(saveDir); end

channels = {'Ch09', 'Ch10', 'Ch11', 'Ch12', 'Ch13', 'Ch14', 'Ch15', 'Ch16'};
rpt_cycles = {'0cyc', '200cyc', '400cyc'};

%% OCV Processing
fprintf('\n=== OCV Processing ===\n');

% OCV conditions: charge (8,2), discharge (10,2)
ocv_conditions = {'charge', 'discharge'};
ocv_steps = [8, 10];

% Structure to store OCV data for each channel
channel_ocv_data = struct();
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    channel_ocv_data.(channel) = struct();
end

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};

    % Figure 1: Individual channel OCV
    figure('Name', sprintf('OCV - %s', channel), 'Position', [100 100 1200 800]);
    set(gcf, 'Visible', 'off');

    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
        filepath = fullfile(ExperimentalDataPath, filename);

        T = readtable(filepath);

        subplot(3,1,rpt_idx);
        hold on;

        charge_ocv = [];
        discharge_ocv = [];

        for ocv_idx = 1:length(ocv_conditions)
            step_idx = ocv_steps(ocv_idx);
            ocv_data = (T{:,2} == step_idx) & (T{:,4} == 2);

            capacity_data = T{ocv_data, 9};
            voltage_data = T{ocv_data, 8};

            % Divide capacity into 0~100 (101 points)
            soc_data = linspace(0, 100, 101);

            % Reverse discharge OCV so voltage increases with SOC
            if strcmp(ocv_conditions{ocv_idx}, 'discharge')
                voltage_data = flipud(voltage_data);
            end

            % Interpolate voltage_data to match soc_data length
            voltage_interp = interp1(linspace(0, 100, length(voltage_data)), voltage_data, soc_data, 'linear');

            if strcmp(ocv_conditions{ocv_idx}, 'charge')
                charge_ocv = voltage_interp;
                plot(soc_data, voltage_interp, 'b-', 'LineWidth', 2, 'DisplayName', 'Charge');
            else
                discharge_ocv = voltage_interp;
                plot(soc_data, voltage_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'Discharge');
            end
        end

        % Calculate and display average charge/discharge OCV
        avg_ocv = (charge_ocv + discharge_ocv) / 2;
        plot(soc_data, avg_ocv, 'g-', 'LineWidth', 2, 'DisplayName', 'Average');

        % Field name mapping
        if strcmp(rpt_cycle, '0cyc')
            field_name = 'cyc0';
        elseif strcmp(rpt_cycle, '200cyc')
            field_name = 'cyc200';
        else
            field_name = 'cyc400';
        end

        % Save soc, avg_ocv, capacity (use last capacity of each step)
        charge_mask = (T{:,2} == ocv_steps(1)) & (T{:,4} == 2);
        discharge_mask = (T{:,2} == ocv_steps(2)) & (T{:,4} == 2);
        channel_ocv_data.(channel).(field_name).soc = soc_data;
        channel_ocv_data.(channel).(field_name).avg_ocv = avg_ocv;
        channel_ocv_data.(channel).(field_name).capacity = (T{find(charge_mask,1,'last'), 9} + T{find(discharge_mask,1,'last'), 9})/2;

        xlabel('SOC [%]');
        ylabel('Voltage [V]');
        title(sprintf('%s - %s OCV', channel, rpt_cycle));
        legend('Location', 'best');
        grid on;
        xlim([0 100]);
    end

    % Save figure
    figName = fullfile(saveDir, sprintf('%s_OCV.fig', channel));
    savefig(gcf, figName);
    close(gcf);  % 그래프 창 닫기
    fprintf('Saved: %s\n', figName);
end

% Figure 2: Average OCV for each cycle separately (Approach A)
figure('Name', 'Average OCV by Cycle', 'Position', [100 100 1200 800]);
hold on;

% Collect per-channel averaged OCV (charge/discharge averaged per channel)
per_channel_avg_rpt0 = [];
per_channel_avg_rpt200 = [];
per_channel_avg_rpt400 = [];
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    per_channel_avg_rpt0 = [per_channel_avg_rpt0; channel_ocv_data.(channel).cyc0.avg_ocv];
    per_channel_avg_rpt200 = [per_channel_avg_rpt200; channel_ocv_data.(channel).cyc200.avg_ocv];
    per_channel_avg_rpt400 = [per_channel_avg_rpt400; channel_ocv_data.(channel).cyc400.avg_ocv];
end

% Average across channels
avg_ocv_rpt0 = mean(per_channel_avg_rpt0, 1, 'omitnan');
avg_ocv_rpt200 = mean(per_channel_avg_rpt200, 1, 'omitnan');
avg_ocv_rpt400 = mean(per_channel_avg_rpt400, 1, 'omitnan');
soc_grid = 0:1:100;

% Sorted grid (already ascending, but keep consistent)
[soc_grid_sorted, sort_idx] = sort(soc_grid);
avg_ocv_rpt0_sorted = avg_ocv_rpt0(sort_idx);
avg_ocv_rpt200_sorted = avg_ocv_rpt200(sort_idx);
avg_ocv_rpt400_sorted = avg_ocv_rpt400(sort_idx);

% Plot all three RPT average OCV
plot(soc_grid_sorted, avg_ocv_rpt0_sorted, 'b-', 'LineWidth', 3, 'DisplayName', 'RPT0 Average');
plot(soc_grid_sorted, avg_ocv_rpt200_sorted, 'r-', 'LineWidth', 3, 'DisplayName', 'RPT200 Average');
plot(soc_grid_sorted, avg_ocv_rpt400_sorted, 'g-', 'LineWidth', 3, 'DisplayName', 'RPT400 Average');
xlabel('SOC [%]');
ylabel('Voltage [V]');
title('Average OCV by Cycle (Approach A: per-channel avg → channel avg)');
legend('Location', 'best');
grid on;
xlim([0 100]);

% Save average OCV
figName = fullfile(saveDir, 'Average_OCV_by_Cycle.fig');
savefig(gcf, figName);
close(gcf);  % 그래프 창 닫기
fprintf('Saved: %s\n', figName);

% Create OCV functions with sorted data
OCV_func_rpt0 = @(soc_query) interp1(soc_grid_sorted, avg_ocv_rpt0_sorted, soc_query, 'linear');
OCV_func_rpt200 = @(soc_query) interp1(soc_grid_sorted, avg_ocv_rpt200_sorted, soc_query, 'linear');
OCV_func_rpt400 = @(soc_query) interp1(soc_grid_sorted, avg_ocv_rpt400_sorted, soc_query, 'linear');

% Compute mean capacity across channels for RPT0, RPT200, and RPT400 (for OCV-Q)
cap_list_rpt0 = [];
cap_list_rpt200 = [];
cap_list_rpt400 = [];
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    cap_list_rpt0 = [cap_list_rpt0, channel_ocv_data.(channel).cyc0.capacity];
    cap_list_rpt200 = [cap_list_rpt200, channel_ocv_data.(channel).cyc200.capacity];
    cap_list_rpt400 = [cap_list_rpt400, channel_ocv_data.(channel).cyc400.capacity];
end
mean_capacity_rpt0 = mean(cap_list_rpt0, 'omitnan');
mean_capacity_rpt200 = mean(cap_list_rpt200, 'omitnan');
mean_capacity_rpt400 = mean(cap_list_rpt400, 'omitnan');

% OCV as a function of charge amount Q (Ah)
q_grid_rpt0 = mean_capacity_rpt0 .* (soc_grid_sorted ./ 100);
q_grid_rpt200 = mean_capacity_rpt200 .* (soc_grid_sorted ./ 100);
q_grid_rpt400 = mean_capacity_rpt400 .* (soc_grid_sorted ./ 100);
OCV_Q_func_rpt0 = @(q_query) interp1(q_grid_rpt0, avg_ocv_rpt0_sorted, q_query, 'linear');
OCV_Q_func_rpt200 = @(q_query) interp1(q_grid_rpt200, avg_ocv_rpt200_sorted, q_query, 'linear');
OCV_Q_func_rpt400 = @(q_query) interp1(q_grid_rpt400, avg_ocv_rpt400_sorted, q_query, 'linear');

% Save OCV arrays and metadata to MAT for later use
OCV_data = struct();
OCV_data.soc_grid = soc_grid_sorted;
OCV_data.avg_ocv_rpt0 = avg_ocv_rpt0_sorted;
OCV_data.avg_ocv_rpt200 = avg_ocv_rpt200_sorted;
OCV_data.avg_ocv_rpt400 = avg_ocv_rpt400_sorted;
OCV_data.mean_capacity_rpt0 = mean_capacity_rpt0;
OCV_data.mean_capacity_rpt200 = mean_capacity_rpt200;
OCV_data.mean_capacity_rpt400 = mean_capacity_rpt400;
OCV_data.q_grid_rpt0 = q_grid_rpt0;
OCV_data.q_grid_rpt200 = q_grid_rpt200;
OCV_data.q_grid_rpt400 = q_grid_rpt400;
OCV_data.OCV_SOC_func_rpt0 = OCV_func_rpt0;
OCV_data.OCV_SOC_func_rpt200 = OCV_func_rpt200;
OCV_data.OCV_SOC_func_rpt400 = OCV_func_rpt400;
OCV_data.OCV_Q_func_rpt0 = OCV_Q_func_rpt0;
OCV_data.OCV_Q_func_rpt200 = OCV_Q_func_rpt200;
OCV_data.OCV_Q_func_rpt400 = OCV_Q_func_rpt400;

%% 8개 채널 통합 OCV 함수 생성 (OCV_integrated.m 로직)
fprintf('\n=== Creating Integrated OCV Functions ===\n');

% 0cyc, 200cyc, 400cyc: 8개 채널의 OCV를 8x101 행렬로 쌓기
all_V_OCV_0cyc = [];
all_V_OCV_200cyc = [];
all_V_OCV_400cyc = [];
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    all_V_OCV_0cyc = [all_V_OCV_0cyc; channel_ocv_data.(channel).cyc0.avg_ocv];
    all_V_OCV_200cyc = [all_V_OCV_200cyc; channel_ocv_data.(channel).cyc200.avg_ocv];
    all_V_OCV_400cyc = [all_V_OCV_400cyc; channel_ocv_data.(channel).cyc400.avg_ocv];
end

% 8채널 평균 OCV 계산
V_avg_SOC_0cyc = mean(all_V_OCV_0cyc, 1, 'omitnan'); % 1x101
V_avg_SOC_200cyc = mean(all_V_OCV_200cyc, 1, 'omitnan'); % 1x101
V_avg_SOC_400cyc = mean(all_V_OCV_400cyc, 1, 'omitnan'); % 1x101

% 통합 OCV 함수 생성
OCV_SOC_func_0cyc = @(soc_query) interp1(soc_grid_sorted, V_avg_SOC_0cyc, soc_query, 'linear');
OCV_SOC_func_200cyc = @(soc_query) interp1(soc_grid_sorted, V_avg_SOC_200cyc, soc_query, 'linear');
OCV_SOC_func_400cyc = @(soc_query) interp1(soc_grid_sorted, V_avg_SOC_400cyc, soc_query, 'linear');

% 통합 구조체 생성
OCV_integrated_struct_0cyc.SOC_grid = soc_grid_sorted;
OCV_integrated_struct_0cyc.V_avg_SOC = V_avg_SOC_0cyc;
OCV_integrated_struct_0cyc.OCV_SOC_func = OCV_SOC_func_0cyc;
OCV_integrated_struct_0cyc.mean_capacity = mean_capacity_rpt0;

OCV_integrated_struct_200cyc.SOC_grid = soc_grid_sorted;
OCV_integrated_struct_200cyc.V_avg_SOC = V_avg_SOC_200cyc;
OCV_integrated_struct_200cyc.OCV_SOC_func = OCV_SOC_func_200cyc;
OCV_integrated_struct_200cyc.mean_capacity = mean_capacity_rpt200;

OCV_integrated_struct_400cyc.SOC_grid = soc_grid_sorted;
OCV_integrated_struct_400cyc.V_avg_SOC = V_avg_SOC_400cyc;
OCV_integrated_struct_400cyc.OCV_SOC_func = OCV_SOC_func_400cyc;
OCV_integrated_struct_400cyc.mean_capacity = mean_capacity_rpt400;

% 통합 OCV 시각화
figure('Name', 'Integrated SOC-OCV Curve', 'Position', [100 100 1200 800]);
plot(soc_grid_sorted, V_avg_SOC_0cyc, 'bo-', 'LineWidth', 2, 'DisplayName', 'RPT0 (0cyc)');
hold on;
plot(soc_grid_sorted, V_avg_SOC_200cyc, 'ro-', 'LineWidth', 2, 'DisplayName', 'RPT200 (200cyc)');
plot(soc_grid_sorted, V_avg_SOC_400cyc, 'go-', 'LineWidth', 2, 'DisplayName', 'RPT400 (400cyc)');
xlabel('SOC [%]');
ylabel('OCV [V]');
title('Integrated SOC-OCV Curve (8-cell Average)');
grid on;
xlim([0 100]);
legend('Location', 'best');

% 통합 그래프 저장
figName = fullfile(saveDir, 'Integrated_OCV_8cell_Average.fig');
savefig(gcf, figName);
close(gcf);  % 그래프 창 닫기
fprintf('Saved: %s\n', figName);

% 통합 데이터를 OCV_data에 추가
OCV_data.OCV_integrated_0cyc = OCV_integrated_struct_0cyc;
OCV_data.OCV_integrated_200cyc = OCV_integrated_struct_200cyc;
OCV_data.OCV_integrated_400cyc = OCV_integrated_struct_400cyc;

matSaveFile = fullfile(saveDir, 'OCV_integrated.mat');
save(matSaveFile, 'OCV_data');
fprintf('Saved OCV data MAT: %s\n', matSaveFile);

fprintf('\n=== RPT Postprocessing Completed ===\n');
