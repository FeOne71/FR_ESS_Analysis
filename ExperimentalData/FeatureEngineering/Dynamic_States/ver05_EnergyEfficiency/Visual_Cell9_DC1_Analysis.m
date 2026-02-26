%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual_Cell9_DC1_Analysis.m
% - 9번 셀(Ch09) 주행부하 DC1 시각화 분석
% - V-Q 그래프 (전압 vs 누적 용량)
% - E-t 그래프 (누적 에너지 vs 시간)
% - 결과는 스크립트명 폴더에 .fig 및 .pdf로 저장
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% Configuration
% File paths
driveCycleDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
scriptName = 'Visual_Cell9_DC1_Analysis';
outputDir = fullfile(pwd, scriptName);
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% RPT cycles
rptCycles = {'0cyc', '200cyc', '400cyc', '600cyc', '800cyc', '1000cyc'};
targetSOC = 'SOC90';
targetDC = 'DC1';

% Color map for cycles
colors = lines(length(rptCycles));

%% Initialize Data Storage
allData = struct();

%% Process Each Cycle
fprintf('\n=== Loading and Processing Data ===\n');

for i = 1:length(rptCycles)
    cycle = rptCycles{i};
    fileName = sprintf('parsedDriveCycle_%s_filtered.mat', cycle);
    filePath = fullfile(driveCycleDir, fileName);
    
    if ~exist(filePath, 'file')
        fprintf('  Warning: File not found: %s\n', fileName);
        continue;
    end
    
    fprintf('  Processing %s...\n', cycle);
    data = load(filePath);
    varName = sprintf('parsedDriveCycle_%s', cycle);
    ch9FieldName = sprintf('ch9_Drive_%s', cycle);
    
    if ~isfield(data, varName) || ~isfield(data.(varName), ch9FieldName)
        fprintf('  Warning: Ch09 data not found in %s\n', fileName);
        continue;
    end
    
    ch9Data = data.(varName).(ch9FieldName);
    
    if ~isfield(ch9Data, targetSOC) || ~isfield(ch9Data.(targetSOC), targetDC)
        fprintf('  Warning: %s %s not found in Ch09 data\n', targetSOC, targetDC);
        continue;
    end
    
    dcData = ch9Data.(targetSOC).(targetDC);
    
    % Extract V, I, t
    V = dcData.V;
    I = dcData.I;
    t = dcData.t;
    
    if isa(t, 'duration')
        t_sec = seconds(t);
    else
        t_sec = t;
    end
    
    % Remove NaN
    validIdx = ~isnan(V) & ~isnan(I) & ~isnan(t_sec);
    V = V(validIdx);
    I = I(validIdx);
    t_sec = t_sec(validIdx);
    
    % Normalize time
    t_sec = t_sec - t_sec(1);
    dt = [0; diff(t_sec)];
    
    % Ah Calculation (Capacity Q)
    % Q = integral |I| dt / 3600
    Q = cumtrapz(t_sec, abs(I)) / 3600;
    
    % Wh Calculation (Energy E)
    % E = integral |V*I| dt / 3600
    E = cumtrapz(t_sec, abs(V .* I)) / 3600;
    
    % Store (Use a valid field name by prefixing with 'cyc_')
    validFieldName = ['cyc_', cycle];
    allData.(validFieldName).t = t_sec;
    allData.(validFieldName).V = V;
    allData.(validFieldName).Q = Q;
    allData.(validFieldName).E = E;
end

%% Visualization 1: V-Q Graph
fprintf('\n=== Generating V-Q Graph ===\n');
fig_vq = figure('Name', 'Voltage vs Cumulative Capacity', 'Position', [100, 100, 800, 600]);
hold on;

for i = 1:length(rptCycles)
    cycle = rptCycles{i};
    validFieldName = ['cyc_', cycle];
    if isfield(allData, validFieldName)
        plot(allData.(validFieldName).Q, allData.(validFieldName).V, 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', cycle);
    end
end

xlabel('Cumulative Capacity [Ah]', 'FontSize', 12);
ylabel('Voltage [V]', 'FontSize', 12);
title(sprintf('V-Q Curve for %s %s %s', 'Ch09', targetSOC, targetDC), 'FontSize', 14);
legend('Location', 'best');
grid on;
grid minor;
hold off;

% Save V-Q Graph
saveas(fig_vq, fullfile(outputDir, 'VQ_Graph.fig'));
saveas(fig_vq, fullfile(outputDir, 'VQ_Graph.pdf'));

%% Visualization 2: E-t Graph
fprintf('\n=== Generating E-t Graph ===\n');
fig_et = figure('Name', 'Cumulative Energy vs Time', 'Position', [100, 100, 800, 600]);
hold on;

for i = 1:length(rptCycles)
    cycle = rptCycles{i};
    validFieldName = ['cyc_', cycle];
    if isfield(allData, validFieldName)
        plot(allData.(validFieldName).t / 3600, allData.(validFieldName).E, 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', cycle);
    end
end

xlabel('Time [hr]', 'FontSize', 12);
ylabel('Cumulative Energy [Wh]', 'FontSize', 12);
title(sprintf('Cumulative Energy vs Time for %s %s %s', 'Ch09', targetSOC, targetDC), 'FontSize', 14);
legend('Location', 'best');
grid on;
grid minor;
hold off;

% Save E-t Graph
saveas(fig_et, fullfile(outputDir, 'Et_Graph.fig'));
saveas(fig_et, fullfile(outputDir, 'Et_Graph.pdf'));

fprintf('\n=== Visualization Completed ===\n');
fprintf('Results saved in: %s\n', outputDir);

%% Report Key Statistics
fprintf('\n=== Summary Statistics for Cell 9 DC1 ===\n');
fprintf('%-10s | %-12s | %-12s | %-12s\n', 'Cycle', 'Final Q [Ah]', 'Final E [Wh]', 'Avg V [V]');
fprintf('------------------------------------------------------------\n');

for i = 1:length(rptCycles)
    cycle = rptCycles{i};
    validFieldName = ['cyc_', cycle];
    if isfield(allData, validFieldName)
        finalQ = allData.(validFieldName).Q(end);
        finalE = allData.(validFieldName).E(end);
        avgV = mean(allData.(validFieldName).V);
        fprintf('%-10s | %-12.4f | %-12.4f | %-12.4f\n', cycle, finalQ, finalE, avgV);
    end
end
