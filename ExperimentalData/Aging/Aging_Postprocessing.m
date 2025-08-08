%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aging Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% 경로 설정
folderPath = 'C:\Users\Chulwon Jung\내 드라이브\JCW_GDrive\Experiment Data\Aging';
saveFolder_fig = fullfile(folderPath, 'figures');
saveFolder_mat = fullfile(folderPath, 'mat');

if ~exist(saveFolder_fig, 'dir'); mkdir(saveFolder_fig); end
if ~exist(saveFolder_mat, 'dir'); mkdir(saveFolder_mat); end

% 파일 목록 수동 지정
fileNames = {
    'Ch09_1C1C_4pcd_V70_0to200cyc.csv';
    'Ch10_1C1C_4cpd_Vmax_0to200cyc.csv';
    'Ch11_1C1C_8cpd_V70_0to200cyc.csv';
    'Ch12_1C1C_8cpd_Vmax_0to200cyc.csv';
    'Ch13_3C1C_4pcd_V70_0to200cyc.csv';
    'Ch14_3C1C_4cpd_Vmax_0to200cyc.csv';
    'Ch15_3C1C_8cpd_V70_0to200cyc.csv';
    'Ch16_3C1C_8cpd_Vmax_0to200cyc.csv'
};

% 색상 지정
color_voltage = '#0073C2';  % 파랑
color_current = '#E18727';  % 주황
color_capacity = '#20854E'; % 초록

for i = 1:length(fileNames)
    filename = fileNames{i};
    filepath = fullfile(folderPath, filename);

    if exist(filepath, 'file') ~= 2
        fprintf('❌ 파일 없음: %s\n', filepath);
        continue;
    end

    T = readtable(filepath);

    stepIndex   = T{:,2};
    stepType    = T{:,3};
    cycleIndex  = T{:,4};
    totalTime   = T{:,6};  % [s]
    current     = T{:,7};  % [A]
    voltage     = T{:,8};  % [V]
    capacity    = T{:,9};  % [Ah]
    time_rel    = seconds(totalTime - totalTime(1));  % 0초 기준

    % 채널 이름
    [~, baseName, ~] = fileparts(filename);
    channelLabel = extractBetween(baseName, 'Ch', '_');
    channelName = sprintf('Channel %s', channelLabel{1});

    parsed.stepIndex = stepIndex;
    parsed.stepType = stepType;
    parsed.cycleIndex = cycleIndex;
    parsed.totalTime = totalTime;
    parsed.current = current;
    parsed.voltage = voltage;
    parsed.capacity = capacity;

    matname = fullfile(saveFolder_mat, sprintf('%s_degradation.mat', baseName));
    save(matname, 'parsed');

    % === Figure ===
    figure('Name', channelName, 'NumberTitle', 'off');

    subplot(2,1,1);
    yyaxis left
    plot(time_rel, voltage, 'Color', color_voltage, 'LineWidth', 1.5);
    ylabel('Voltage [V]');
    ylim([min(voltage)*0.98, max(voltage)*1.02]);

    yyaxis right
    plot(time_rel, current, 'Color', color_current, 'LineWidth', 1.5);
    ylabel('Current [A]');
    ylim([min(current)*1.2, max(current)*1.2]);

    xlabel('Time [s]');
    title(sprintf('%s - Voltage & Current vs Time', channelName));
    grid on;


    % 2. 시간–용량
    subplot(2,1,2); hold on;

    plot(time_rel, capacity, 'Color', color_capacity, 'LineWidth', 1.5);

    stepType = T{:,3};  
    maxCycle = max(cycleIndex);  
    mask_Qend = strcmp(stepType, 'CC DChg') & (cycleIndex == maxCycle);

    Qend = capacity(find(mask_Qend, 1, 'last'));
    legendStr = sprintf('Cycle %d - Final Q = %.2f Ah', maxCycle, Qend);

    legend(legendStr, 'Location', 'best');
    xlabel('Time [s]');
    ylabel('Capacity [Ah]');
    title(sprintf('%s - Capacity vs Time', channelName));
    grid on;

    %% Figure 2: 시간 vs 용량 
    figure('Name', [channelName ' - Capacity vs Time'], 'Color', 'w');

    plot(time_rel, capacity, 'Color', color_capacity, 'LineWidth', 1.5); hold on;

    % Qend
    stepType = T{:,3}; 
    maxCycle = max(cycleIndex);
    mask_Qend = strcmp(stepType, 'CC DChg') & (cycleIndex == maxCycle);
    Qend = capacity(find(mask_Qend, 1, 'last'));

    legendStr = sprintf('Cycle %d - Final Q = %.2f Ah', maxCycle, Qend);
    legend(legendStr, 'Location', 'best');

    xlabel('Time [s]');
    ylabel('Capacity [Ah]');
    ylim([0 70]);
    title(sprintf('%s - Capacity vs Time', channelName));
    grid on;

end