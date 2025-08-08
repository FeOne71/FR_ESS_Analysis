%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drive Cycle Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% 경로 설정
folderPath = 'C:\Users\Chulwon Jung\내 드라이브\JCW_GDrive\Experiment Data\Drive Cycle';
saveFolder = 'C:\Users\Chulwon Jung\내 드라이브\JCW_GDrive\Experiment Data\Drive Cycle\figures';
if ~exist(saveFolder, 'dir'); mkdir(saveFolder); end

% 파일 목록 수동 입력
fileNames = {
    'Ch9_Drive_0cyc.csv';
    'Ch10_Drive_0cyc.csv';
    'Ch11_Drive_0cyc.csv';
    'Ch12_Drive_0cyc.csv';
    'Ch13_Drive_0cyc.csv';
    'Ch14_Drive_0cyc.csv';
    'Ch15_Drive_0cyc.csv';
    'Ch16_Drive_0cyc.csv'
};

% 색상
color_voltage = '#0073C2';
color_current = '#E18727';

% 저장 폴더 생성
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

for i = 1:length(fileNames)
    filename = fileNames{i};
    filepath = fullfile(folderPath, filename);

    if exist(filepath, 'file') ~= 2
        fprintf('파일 없음: %s\n', filepath);
        continue;
    end

    T = readtable(filepath);

    stepIndex   = T{:,2};
    cycleIndex  = T{:,4};
    time        = T{:,5};  % [s]
    totalTime   = T{:,6};  % [s]
    current     = T{:,7};  % [A]
    voltage     = T{:,8};  % [V]
    capacity    = T{:,9};  % [Ah]

    % 시간 기준 정렬 (0 기준)
    time_rel = seconds(totalTime - totalTime(1));

    % 채널 이름 추출
    [~, baseName, ~] = fileparts(filename);
    channelLabel = extractBetween(baseName, 'Ch', '_');
    channelName = sprintf('Channel %s', channelLabel{1});

    % === 데이터 구조체 저장 ===
    parsed.stepIndex   = stepIndex;
    parsed.cycleIndex  = cycleIndex;
    parsed.time        = time;
    parsed.totalTime   = totalTime;
    parsed.current     = current;
    parsed.voltage     = voltage;
    parsed.capacity    = capacity;

    matname = fullfile(saveFolder, sprintf('%s_real_load.mat', baseName));
    save(matname, 'parsed');

    % === Figure ===
    figure('Name', channelName, 'NumberTitle', 'off');

    yyaxis left
    plot(parsed.totalTime, parsed.voltage, 'Color', color_voltage, 'LineWidth', 1.5);
    ylabel('Voltage [V]');
    ylim([min(parsed.voltage)*0.98, max(parsed.voltage)*1.02]);

    yyaxis right
    plot(parsed.totalTime, parsed.current, 'Color', color_current, 'LineWidth', 1.5);
    ylabel('Current [A]');
    ylim([min(parsed.current)*1.1, max(parsed.current)*1.1]);

    xlabel('Time [s]');
    % title(sprintf('%s - Time vs Voltage/Current (Static Capacity)', baseName));
    title(sprintf('%s - Time vs Voltage/Current (Drive Cycle)', channelName));
    grid on;

    % === Figure 저장 ===
    fig_fig = fullfile(saveFolder, sprintf('%s_real_load.fig', baseName));
    savefig(gcf, fig_fig);
end