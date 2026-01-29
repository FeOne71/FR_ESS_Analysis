%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DriveCycle_Visualization.m
% - 시각화 주행부하 8개 (DC1~DC8)
% - Figure 1: 서브플롯 형태, x축 시간, y축 전압
% - Figure 2: 서브플롯 형태, x축 시간, y축 전류
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% Configuration
% File paths
driveCycleDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\ver05_EnergyEfficiency\Results';
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% Target channel
targetChannel = 'Ch09';

% RPT cycles
rptCycles = {'0cyc', '200cyc', '400cyc', '600cyc', '800cyc', '1000cyc'};
rptCycleNums = [0, 200, 400, 600, 800, 1000];

% SOC levels
socLevels = {'SOC90', 'SOC70', 'SOC50'};

% Drive cycle types
dcTypes = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

% Visualization settings
targetCycle = '0cyc';  % 시각화할 사이클 선택 (변경 가능)
targetSOC = 'SOC90';   % 시각화할 SOC 레벨 선택 (변경 가능)

%% Find drive cycle data files
fprintf('\n=== Finding Drive Cycle Data Files ===\n');
driveCycleFiles = dir(fullfile(driveCycleDir, 'parsedDriveCycle_*cyc_filtered.mat'));

if isempty(driveCycleFiles)
    error('No parsedDriveCycle files found in %s. Please check the path.', driveCycleDir);
end

fprintf('Found %d drive cycle files in: %s\n', length(driveCycleFiles), driveCycleDir);

%% Load target cycle data
fprintf('\n=== Loading Target Cycle Data ===\n');
targetFileName = sprintf('parsedDriveCycle_%s_filtered.mat', targetCycle);
targetFilePath = fullfile(driveCycleDir, targetFileName);

if ~exist(targetFilePath, 'file')
    error('Target file not found: %s', targetFilePath);
end

fprintf('Loading: %s\n', targetFileName);
driveData = load(targetFilePath);

% Find the variable name
varName = sprintf('parsedDriveCycle_%s', targetCycle);
if ~isfield(driveData, varName)
    error('Variable %s not found in file', varName);
end

parsedData = driveData.(varName);

% Find Ch09 data
ch9FieldName = sprintf('ch9_Drive_%s', targetCycle);
if ~isfield(parsedData, ch9FieldName)
    error('Ch09 data not found: %s', ch9FieldName);
end

ch9Data = parsedData.(ch9FieldName);

% Check if target SOC level exists
if ~isfield(ch9Data, targetSOC)
    error('SOC level %s not found in Ch09 data', targetSOC);
end

socData = ch9Data.(targetSOC);

fprintf('Successfully loaded data for %s %s\n', targetCycle, targetSOC);

%% Prepare data for visualization
fprintf('\n=== Preparing Data for Visualization ===\n');
dcData = struct();

for dcIdx = 1:length(dcTypes)
    dcType = dcTypes{dcIdx};
    
    if ~isfield(socData, dcType)
        fprintf('  Warning: %s not found in data. Skipping.\n', dcType);
        continue;
    end
    
    dcInfo = socData.(dcType);
    
    % Check required fields
    if ~isfield(dcInfo, 'V') || ~isfield(dcInfo, 'I') || ~isfield(dcInfo, 't')
        fprintf('  Warning: %s missing required fields (V, I, t). Skipping.\n', dcType);
        continue;
    end
    
    V = dcInfo.V;
    I = dcInfo.I;
    t = dcInfo.t;
    
    % Convert time to seconds if needed
    if isa(t, 'duration')
        t_sec = seconds(t);
    else
        t_sec = t;
    end
    
    % Remove NaN and ensure column vectors
    validIdx = ~isnan(V) & ~isnan(I) & ~isnan(t_sec);
    if sum(validIdx) < 2
        fprintf('  Warning: %s has insufficient valid data. Skipping.\n', dcType);
        continue;
    end
    
    V = V(validIdx);
    I = I(validIdx);
    t_sec = t_sec(validIdx);
    
    % Normalize time to start from 0
    t_sec = t_sec - t_sec(1);
    
    % Convert to minutes for better readability
    t_min = t_sec / 60;
    
    % Store data
    dcData.(dcType).t_min = t_min;
    dcData.(dcType).V = V;
    dcData.(dcType).I = I;
    
    fprintf('  %s: %d data points, duration: %.2f min\n', dcType, length(t_min), max(t_min));
end

%% Create Figure 1: Voltage vs Time (8 subplots)
fprintf('\n=== Creating Figure 1: Voltage vs Time ===\n');
fig1 = figure('Name', sprintf('Voltage vs Time - %s %s', targetCycle, targetSOC), ...
    'Position', [100, 100, 1800, 1000]);
set(gcf, 'Visible', 'off');

for dcIdx = 1:length(dcTypes)
    dcType = dcTypes{dcIdx};
    
    subplot(2, 4, dcIdx);
    
    if ~isfield(dcData, dcType)
        text(0.5, 0.5, sprintf('%s\nNo Data', dcType), ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        axis off;
        continue;
    end
    
    t_vis = dcData.(dcType).t_min;
    V_vis = dcData.(dcType).V;
    
    plot(t_vis, V_vis, 'b-', 'LineWidth', 1.5);
    xlabel('Time [min]', 'FontSize', 10);
    ylabel('Voltage [V]', 'FontSize', 10);
    title(sprintf('%s', dcType), 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    grid minor;
end

sgtitle(sprintf('Voltage vs Time - %s %s', targetCycle, targetSOC), ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save figure
figDir = fullfile(outputDir, 'Figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

fig1Name = sprintf('DriveCycle_Voltage_%s_%s', targetCycle, targetSOC);
savefig(fig1, fullfile(figDir, [fig1Name, '.fig']));
saveas(fig1, fullfile(figDir, [fig1Name, '.png']));
close(fig1);
fprintf('Saved figure: %s\n', fig1Name);

%% Create Figure 2: Current vs Time (8 subplots)
fprintf('\n=== Creating Figure 2: Current vs Time ===\n');
fig2 = figure('Name', sprintf('Current vs Time - %s %s', targetCycle, targetSOC), ...
    'Position', [100, 100, 1800, 1000]);
set(gcf, 'Visible', 'off');

for dcIdx = 1:length(dcTypes)
    dcType = dcTypes{dcIdx};
    
    subplot(2, 4, dcIdx);
    
    if ~isfield(dcData, dcType)
        text(0.5, 0.5, sprintf('%s\nNo Data', dcType), ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        axis off;
        continue;
    end
    
    t_vis = dcData.(dcType).t_min;
    I_vis = dcData.(dcType).I;
    
    plot(t_vis, I_vis, 'r-', 'LineWidth', 1.5);
    xlabel('Time [min]', 'FontSize', 10);
    ylabel('Current [A]', 'FontSize', 10);
    title(sprintf('%s', dcType), 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    grid minor;
end

sgtitle(sprintf('Current vs Time - %s %s', targetCycle, targetSOC), ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save figure
fig2Name = sprintf('DriveCycle_Current_%s_%s', targetCycle, targetSOC);
savefig(fig2, fullfile(figDir, [fig2Name, '.fig']));
saveas(fig2, fullfile(figDir, [fig2Name, '.png']));
close(fig2);
fprintf('Saved figure: %s\n', fig2Name);

fprintf('\n=== Drive Cycle Visualization Completed ===\n');
