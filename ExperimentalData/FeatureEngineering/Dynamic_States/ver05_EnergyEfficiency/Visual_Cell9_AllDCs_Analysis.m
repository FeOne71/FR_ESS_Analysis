%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual_Cell9_AllDCs_Analysis.m
% - 9번 셀(Ch09) 주행부하 8개(DC1~DC8) 통합 시각화 분석
% - V-Q 그래프 (전압 vs 누적 용량) - 2x4 Subplots
% - E-t 그래프 (누적 에너지 vs 시간) - 2x4 Subplots
% - 결과는 스크립트명 폴더에 .fig 및 .pdf로 저장
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% Configuration
% File paths
driveCycleDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
scriptName = 'Visual_Cell9_AllDCs_Analysis';
outputDir = fullfile(pwd, scriptName);
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% RPT cycles
rptCycles = {'0cyc', '200cyc', '400cyc', '600cyc', '800cyc', '1000cyc'};
targetSOC = 'SOC90';
dcTypes = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

% Color map for cycles
colors = lines(length(rptCycles));

%% Initialize Data Storage
allData = struct();

%% Process Data
fprintf('\n=== Loading and Processing Data for All DCs ===\n');

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
    validFieldName = ['cyc_', cycle];
    
    for d = 1:length(dcTypes)
        dcType = dcTypes{d};
        
        if ~isfield(ch9Data, targetSOC) || ~isfield(ch9Data.(targetSOC), dcType)
            continue;
        end
        
        dcData = ch9Data.(targetSOC).(dcType);
        
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
        
        if isempty(V)
            continue;
        end
        
        % Normalize time
        t_sec = t_sec - t_sec(1);
        
        % Ah/Wh Calculation
        Q = cumtrapz(t_sec, abs(I)) / 3600;
        E = cumtrapz(t_sec, abs(V .* I)) / 3600;
        
        % Store
        allData.(validFieldName).(dcType).t = t_sec;
        allData.(validFieldName).(dcType).V = V;
        allData.(validFieldName).(dcType).Q = Q;
        allData.(validFieldName).(dcType).E = E;
    end
end

%% Visualization 1: V-Q Subplots
fprintf('\n=== Generating V-Q Subplots ===\n');
fig_vq = figure('Name', 'Voltage vs Cumulative Capacity (DC1~DC8)', 'Position', [100, 50, 1600, 900]);

for d = 1:length(dcTypes)
    dcType = dcTypes{d};
    subplot(2, 4, d);
    hold on;
    
    for i = 1:length(rptCycles)
        cycle = rptCycles{i};
        validFieldName = ['cyc_', cycle];
        if isfield(allData, validFieldName) && isfield(allData.(validFieldName), dcType)
            plot(allData.(validFieldName).(dcType).Q, allData.(validFieldName).(dcType).V, ...
                'Color', colors(i,:), 'LineWidth', 1, 'DisplayName', cycle);
        end
    end
    
    xlabel('Cum. Q [Ah]');
    ylabel('Voltage [V]');
    title(dcType, 'FontWeight', 'bold');
    grid on;
    if d == 1; legend('Location', 'best', 'FontSize', 8); end
end
sgtitle(sprintf('V-Q Curve Comparison (Cell 9, %s, DC1-DC8)', targetSOC));

% Save
saveas(fig_vq, fullfile(outputDir, 'VQ_AllDCs_Subplots.fig'));

%% Visualization 2: E-t Subplots
fprintf('\n=== Generating E-t Subplots ===\n');
fig_et = figure('Name', 'Cumulative Energy vs Time (DC1~DC8)', 'Position', [100, 50, 1600, 900]);

for d = 1:length(dcTypes)
    dcType = dcTypes{d};
    subplot(2, 4, d);
    hold on;
    
    for i = 1:length(rptCycles)
        cycle = rptCycles{i};
        validFieldName = ['cyc_', cycle];
        if isfield(allData, validFieldName) && isfield(allData.(validFieldName), dcType)
            plot(allData.(validFieldName).(dcType).t / 3600, allData.(validFieldName).(dcType).E, ...
                'Color', colors(i,:), 'LineWidth', 1, 'DisplayName', cycle);
        end
    end
    
    xlabel('Time [hr]');
    ylabel('Cum. Energy [Wh]');
    title(dcType, 'FontWeight', 'bold');
    grid on;
    if d == 1; legend('Location', 'best', 'FontSize', 8); end
end
sgtitle(sprintf('E-t Curve Comparison (Cell 9, %s, DC1-DC8)', targetSOC));

% Save
saveas(fig_et, fullfile(outputDir, 'Et_AllDCs_Subplots.fig'));

fprintf('\n=== Visualization Completed ===\n');
fprintf('Results saved in: %s\n', outputDir);
