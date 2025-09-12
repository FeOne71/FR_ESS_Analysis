%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldData_RPT_Comparison.m
% Compare RPT resistance across multiple years (2021, 2023, 2025)
% - Load results from FieldData_RPT_Old and FieldData_RPT_New
% - Visualize R10s resistance trends across years
% - Create comparison table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

%% Paths
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\FieldData_RPT_results';

% Data files
dataFiles = {
    'FieldData_RPT_Old_20210607.mat',     % 2021
    'FieldData_RPT_New_20231016.mat',     % 2023  
    'FieldData_RPT_New_20250711.mat'      % 2025
};

% Corresponding dates
dates = {
    '2021/06/07',
    '2023/10/16', 
    '2025/07/11'
};

% Colors for each year
yearColors = [
    0.0000 0.4470 0.7410;  % 2021 - blue
    0.8500 0.3250 0.0980;  % 2023 - orange
    0.4940 0.1840 0.5560   % 2025 - purple
];

yearMarkers = {'o', 's', '^'};
yearNames = {'2021', '2023', '2025'};

%% Load data
Results = cell(1, numel(dataFiles));
for i = 1:numel(dataFiles)
    filePath = fullfile(saveDir, dataFiles{i});
    if exist(filePath, 'file')
        S = load(filePath);
        Results{i} = S.Result;
        fprintf('Loaded: %s\n', dataFiles{i});
    else
        fprintf('Warning: File not found: %s\n', filePath);
        Results{i} = [];
    end
end

%% Extract R1s, R5s, R10s data for comparison
% Use Rack01 for all years (assuming single rack in new data)
rackName = 'Rack01';
R1s_data = cell(1, numel(Results));
R5s_data = cell(1, numel(Results));
R10s_data = cell(1, numel(Results));
validYears = [];

for i = 1:numel(Results)
    if ~isempty(Results{i}) && isfield(Results{i}, rackName)
        if isfield(Results{i}.(rackName), 'R_1s_mOhm') && ...
           isfield(Results{i}.(rackName), 'R_5s_mOhm') && ...
           isfield(Results{i}.(rackName), 'R_10s_mOhm')
            R1s_data{i} = Results{i}.(rackName).R_1s_mOhm;
            R5s_data{i} = Results{i}.(rackName).R_5s_mOhm;
            R10s_data{i} = Results{i}.(rackName).R_10s_mOhm;
            validYears = [validYears, i];
        end
    end
end

if isempty(validYears)
    error('No valid R10s data found in any of the loaded files');
end

%% Create comparison figure
fig = figure('Name','RPT Resistance Comparison Across Years','NumberTitle','off');
set(gcf,'Position',[100 100 1200 600]);

% Use tiledlayout for better control
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Left tile: resistance plot
ax1 = nexttile(tl,1);
hold(ax1,'on'); grid(ax1,'on');

% Plot R10s for each valid year
for i = validYears
    R_data = R10s_data{i};
    if ~isempty(R_data) && any(~isnan(R_data))
        valid_pulses = find(~isnan(R_data));
        if ~isempty(valid_pulses)
            plot(ax1, valid_pulses, R_data(valid_pulses), '-', ...
                'Color', yearColors(i,:), 'LineWidth', 2, ...
                'Marker', yearMarkers{i}, 'MarkerSize', 8, ...
                'DisplayName', yearNames{i});
        end
    end
end

xlabel(ax1,'Pulse Index'); ylabel(ax1,'R_{10s} [m\Omega]');
title(ax1,'R10s Resistance Comparison Across Years');
xticks(ax1,1:10); xlim(ax1,[0.5 10.5]);
legend(ax1,'Location','best');
set(ax1,'FontSize',12);

% Right tile: comparison table - R10s only
ax2 = nexttile(tl,2);
posTbl = ax2.Position; delete(ax2);

% Prepare table data - R10s only
maxPulses = 10;
varNames = cell(1, maxPulses + 1);  % Pulse rows + Average
tblData = cell(maxPulses + 1, numel(validYears));
colNames = cell(1, numel(validYears));

% Headers - one column per year
for j = 1:numel(validYears)
    colNames{j} = dates{validYears(j)};
end

% Data rows for each pulse - R10s only
for p = 1:maxPulses
    varNames{p} = sprintf('Pulse %d', p);
    for j = 1:numel(validYears)
        i = validYears(j);
        
        % R10s data only
        if ~isempty(R10s_data{i}) && p <= numel(R10s_data{i})
            if ~isnan(R10s_data{i}(p))
                tblData{p,j} = sprintf('%.4f', abs(R10s_data{i}(p)));
            else
                tblData{p,j} = 'N/A';
            end
        else
            tblData{p,j} = 'N/A';
        end
    end
end

% Add average row - R10s only
varNames{maxPulses+1} = 'Average';
for j = 1:numel(validYears)
    i = validYears(j);
    
    % Calculate R10s average only
    R10s_avg = mean(R10s_data{i}(~isnan(R10s_data{i})));
    
    if ~isnan(R10s_avg)
        tblData{maxPulses+1,j} = sprintf('%.4f', abs(R10s_avg));
    else
        tblData{maxPulses+1,j} = 'N/A';
    end
end

% Create uitable
uit = uitable(fig, 'Data', tblData, 'ColumnName', colNames, 'RowName', varNames, ...
             'Units','normalized', 'Position', posTbl);
uit.ColumnWidth = 'auto';

set(findall(fig,'-property','FontSize'),'FontSize',12);

%% Print summary statistics
fprintf('\n=== R10s Resistance Summary ===\n');
fprintf('%-12s %-10s %-10s %-10s %-10s\n', 'Year', 'Mean', 'Std', 'Min', 'Max');
fprintf('%s\n', repmat('-', 1, 60));

for j = 1:numel(validYears)
    i = validYears(j);
    R_data = R10s_data{i};
    valid_data = R_data(~isnan(R_data));
    
    if ~isempty(valid_data)
        fprintf('%-12s %-10.2f %-10.2f %-10.2f %-10.2f\n', ...
            yearNames{i}, mean(valid_data), std(valid_data), ...
            min(valid_data), max(valid_data));
    else
        fprintf('%-12s %-10s %-10s %-10s %-10s\n', ...
            yearNames{i}, 'N/A', 'N/A', 'N/A', 'N/A');
    end
end

fprintf('\nComparison completed successfully!\n');
