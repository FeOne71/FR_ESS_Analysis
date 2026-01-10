%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 03_DriveCycle_DataSplit.m
% 4-Fold Cross Validation 데이터 분할
% 
% 목적: 
% - Summary Table을 4-Fold CV에 맞게 분할
% - 셀 기반 층화 분할 (Stratified Unit-based Split)
% - Group A (Mild): Ch09-12, Group B (Harsh): Ch13-16
% - 각 Fold: Test [Group A 1개 + Group B 1개], Train [나머지 6개]
%
% 입력:
% - DriveCycle_Summary_Table.mat
%
% 출력:
% - DriveCycle_Fold1_Train.mat, DriveCycle_Fold1_Test.mat
% - DriveCycle_Fold2_Train.mat, DriveCycle_Fold2_Test.mat
% - DriveCycle_Fold3_Train.mat, DriveCycle_Fold3_Test.mat
% - DriveCycle_Fold4_Train.mat, DriveCycle_Fold4_Test.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

fprintf('=== Drive Cycle Data Splitting (4-Fold CV) ===\n');

%% Configuration - User Settings
% =========================================================================
inputFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01\DriveCycle_Summary_Table.mat';
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Data groups
groupA_channels = [9, 10, 11, 12];  % Mild Aging (1C 충전)
groupB_channels = [13, 14, 15, 16]; % Harsh Aging (3C 충전)

% 4-Fold definition
fold_definitions = {
    struct('test_channels', [9, 13], 'fold_num', 1),  % Fold 1: Test [Ch09, Ch13]
    struct('test_channels', [10, 14], 'fold_num', 2),  % Fold 2: Test [Ch10, Ch14]
    struct('test_channels', [11, 15], 'fold_num', 3),  % Fold 3: Test [Ch11, Ch15]
    struct('test_channels', [12, 16], 'fold_num', 4)   % Fold 4: Test [Ch12, Ch16]
};
% =========================================================================

%% Load Summary Table
fprintf('\n=== Loading Summary Table ===\n');
if ~exist(inputFile, 'file')
    error('Summary table not found: %s\nPlease run DriveCycle_DataAggregation_02.m first', inputFile);
end

load(inputFile, 'summaryTable');
fprintf('Loaded summary table: %d rows, %d columns\n', height(summaryTable), width(summaryTable));

%% Display data statistics
fprintf('\n=== Data Statistics ===\n');
fprintf('Cycles: %s\n', mat2str(unique(summaryTable.Cycle)'));
fprintf('Channels: %s\n', mat2str(unique(summaryTable.Channel)'));
fprintf('DC Profiles: %s\n', strjoin(unique(summaryTable.DC_Profile), ', '));
fprintf('Event Types: %s\n', strjoin(unique(summaryTable.EventType), ', '));

fprintf('\nGroup A (Mild Aging): Channels %s\n', mat2str(groupA_channels));
fprintf('Group B (Harsh Aging): Channels %s\n', mat2str(groupB_channels));

%% Create 4-Fold splits
fprintf('\n=== Creating 4-Fold Splits ===\n');

fold_results = struct();

for fold_idx = 1:length(fold_definitions)
    fold_info = fold_definitions{fold_idx};
    fold_num = fold_info.fold_num;
    test_channels = fold_info.test_channels;
    
    fprintf('\n--- Fold %d ---\n', fold_num);
    fprintf('Test channels: %s\n', mat2str(test_channels));
    
    % Split data
    testIdx = ismember(summaryTable.Channel, test_channels);
    trainIdx = ~testIdx;
    
    trainData = summaryTable(trainIdx, :);
    testData = summaryTable(testIdx, :);
    
    % Verify train channels
    train_channels = unique(trainData.Channel);
    fprintf('Train channels: %s (%d channels)\n', mat2str(train_channels), length(train_channels));
    
    % Statistics
    fprintf('Train set: %d rows\n', height(trainData));
    fprintf('Test set: %d rows\n', height(testData));
    
    % Check data distribution
    fprintf('Train cycles: %s\n', mat2str(unique(trainData.Cycle)'));
    fprintf('Test cycles: %s\n', mat2str(unique(testData.Cycle)'));
    
    % Store results
    fold_results.(sprintf('Fold%d', fold_num)) = struct();
    fold_results.(sprintf('Fold%d', fold_num)).trainData = trainData;
    fold_results.(sprintf('Fold%d', fold_num)).testData = testData;
    fold_results.(sprintf('Fold%d', fold_num)).test_channels = test_channels;
    fold_results.(sprintf('Fold%d', fold_num)).train_channels = train_channels;
    
    % Save individual fold files
    savePath_train = fullfile(outputDir, sprintf('DriveCycle_Fold%d_Train.mat', fold_num));
    savePath_test = fullfile(outputDir, sprintf('DriveCycle_Fold%d_Test.mat', fold_num));
    
    save(savePath_train, 'trainData', 'test_channels', 'train_channels', 'fold_num');
    save(savePath_test, 'testData', 'test_channels', 'train_channels', 'fold_num');
    
    fprintf('Saved: %s\n', savePath_train);
    fprintf('Saved: %s\n', savePath_test);
end

%% Save fold definitions
savePath_folds = fullfile(outputDir, 'DriveCycle_4Fold_Definitions.mat');
save(savePath_folds, 'fold_definitions', 'groupA_channels', 'groupB_channels', 'fold_results');
fprintf('\nSaved fold definitions: %s\n', savePath_folds);

%% Summary
fprintf('\n=== Split Summary ===\n');
fprintf('Total folds: %d\n', length(fold_definitions));
fprintf('Test set size per fold: 2 channels\n');
fprintf('Train set size per fold: 6 channels\n');

%% Generate Data Split Visualization
fprintf('\n=== Generating Data Split Visualization ===\n');
figuresDir = fullfile(outputDir, 'figures', 'DataSplit');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

% Each fold distribution visualization
for f = 1:length(fold_definitions)
    fold_name = sprintf('Fold%d', f);
    if isfield(fold_results, fold_name)
        tr = fold_results.(fold_name).trainData;
        ts = fold_results.(fold_name).testData;
        
        fig = figure('Name', sprintf('Data Split Distribution - %s', fold_name), ...
            'Position', [100, 100, 1200, 800], 'Visible', 'on');
        
        % 1. Cycle별 분포
        subplot(2, 3, 1);
        histogram(tr.Cycle, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Train');
        hold on;
        histogram(ts.Cycle, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'Test');
        xlabel('Cycle', 'FontSize', 12);
        ylabel('Count', 'FontSize', 12);
        title('Data Distribution by Cycle', 'FontSize', 14);
        legend('Location', 'best');
        grid on;
        
        % 2. Channel별 분포
        subplot(2, 3, 2);
        train_channels = unique(tr.Channel);
        test_channels = unique(ts.Channel);
        all_channels = unique([train_channels; test_channels]);
        train_counts = histcounts(tr.Channel, [all_channels; max(all_channels)+1]);
        test_counts = histcounts(ts.Channel, [all_channels; max(all_channels)+1]);
        bar_data = [train_counts(:), test_counts(:)];
        b = bar(all_channels, bar_data);
        b(1).FaceColor = 'b';
        b(1).FaceAlpha = 0.5;
        b(2).FaceColor = 'r';
        b(2).FaceAlpha = 0.5;
        xlabel('Channel', 'FontSize', 12);
        ylabel('Count', 'FontSize', 12);
        title('Data Distribution by Channel', 'FontSize', 14);
        legend('Train', 'Test', 'Location', 'best');
        grid on;
        
        % 3. Capacity 분포 (Label Distribution)
        subplot(2, 3, 3);
        if ismember('Capacity_C3', tr.Properties.VariableNames)
            histogram(tr.Capacity_C3, 20, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Train');
            hold on;
            histogram(ts.Capacity_C3, 20, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'Test');
            xlabel('Capacity (Ah)', 'FontSize', 12);
            ylabel('Count', 'FontSize', 12);
            title('Label (Capacity) Distribution', 'FontSize', 14);
            legend('Location', 'best');
            grid on;
        end
        
        % 4. SOC별 분포
        subplot(2, 3, 4);
        if ismember('SOC', tr.Properties.VariableNames)
            train_soc = tr.SOC(~isnan(tr.SOC));
            test_soc = ts.SOC(~isnan(ts.SOC));
            if ~isempty(train_soc) && ~isempty(test_soc)
                all_soc = unique([train_soc; test_soc]);
                train_soc_counts = histcounts(train_soc, [all_soc; max(all_soc)+1]);
                test_soc_counts = histcounts(test_soc, [all_soc; max(all_soc)+1]);
                bar_data = [train_soc_counts(:), test_soc_counts(:)];
                b = bar(all_soc, bar_data);
                b(1).FaceColor = 'b';
                b(1).FaceAlpha = 0.5;
                b(2).FaceColor = 'r';
                b(2).FaceAlpha = 0.5;
                xlabel('SOC (%)', 'FontSize', 12);
                ylabel('Count', 'FontSize', 12);
                title('Data Distribution by SOC', 'FontSize', 14);
                legend('Train', 'Test', 'Location', 'best');
                grid on;
            end
        end
        
        % 5. EventType별 분포
        subplot(2, 3, 5);
        if ismember('EventType', tr.Properties.VariableNames)
            train_et = categorical(tr.EventType);
            test_et = categorical(ts.EventType);
            train_et_counts = countcats(train_et);
            test_et_counts = countcats(test_et);
            categories_et = categories(train_et);
            if isempty(categories_et)
                categories_et = categories(test_et);
            end
            bar_data = [train_et_counts(:), test_et_counts(:)];
            b = bar(categorical(categories_et), bar_data);
            b(1).FaceColor = 'b';
            b(1).FaceAlpha = 0.5;
            b(2).FaceColor = 'r';
            b(2).FaceAlpha = 0.5;
            xlabel('Event Type', 'FontSize', 12);
            ylabel('Count', 'FontSize', 12);
            title('Data Distribution by EventType', 'FontSize', 14);
            legend('Train', 'Test', 'Location', 'best');
            grid on;
        end
        
        % 6. Data size summary
        subplot(2, 3, 6);
        axis off;
        summary_text = sprintf('Fold %d Summary\n\n', f);
        summary_text = [summary_text, sprintf('Train: %d rows\n', height(tr))];
        summary_text = [summary_text, sprintf('Test: %d rows\n', height(ts))];
        summary_text = [summary_text, sprintf('\nTest Channels: %s\n', mat2str(fold_results.(fold_name).test_channels))];
        summary_text = [summary_text, sprintf('Train Channels: %s\n', mat2str(fold_results.(fold_name).train_channels))];
        text(0.1, 0.5, summary_text, 'FontSize', 12, 'VerticalAlignment', 'middle');
        
        sgtitle(sprintf('Data Split Summary - %s (Test Channels: %s)', fold_name, ...
            mat2str(fold_results.(fold_name).test_channels)), 'FontSize', 16);
        
        savefig(fullfile(figuresDir, sprintf('Split_Dist_%s.fig', fold_name)));
        close(fig);
        fprintf('  Saved: Split_Dist_%s.fig\n', fold_name);
    end
end
fprintf('Data split visualizations saved to: %s\n', figuresDir);

%% Generate Data Split Diagram (8 Cells Icon)
fprintf('\n=== Generating Data Split Diagram ===\n');

% Create a schematic diagram showing 8 cells split into Train (6) and Test (2)
fig_diagram = figure('Name', 'Data Split Diagram', 'Position', [200, 200, 1200, 600], 'Visible', 'on');
axis off;

% Define cell positions
cell_width = 0.08;
cell_height = 0.15;
cell_spacing = 0.02;
start_x = 0.1;
start_y = 0.5;

% Draw Train box
train_box_x = 0.1;
train_box_y = 0.2;
train_box_width = 0.35;
train_box_height = 0.6;
rectangle('Position', [train_box_x, train_box_y, train_box_width, train_box_height], ...
          'EdgeColor', 'b', 'LineWidth', 3, 'FaceColor', [0.9, 0.9, 1]);
text(train_box_x + train_box_width/2, train_box_y + train_box_height + 0.05, ...
     'Train Set (6 Cells)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'b');

% Draw Test box
test_box_x = 0.55;
test_box_y = 0.2;
test_box_width = 0.35;
test_box_height = 0.6;
rectangle('Position', [test_box_x, test_box_y, test_box_width, test_box_height], ...
          'EdgeColor', 'r', 'LineWidth', 3, 'FaceColor', [1, 0.9, 0.9]);
text(test_box_x + test_box_width/2, test_box_y + test_box_height + 0.05, ...
     'Test Set (2 Cells)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'r');

% Draw cells in Train box (6 cells)
train_channels = [];
for f = 1:length(fold_definitions)
    fold_name = sprintf('Fold%d', f);
    if isfield(fold_results, fold_name)
        train_channels = fold_results.(fold_name).train_channels;
        break; % Use first fold as example
    end
end

if isempty(train_channels)
    train_channels = [9, 10, 11, 12, 13, 14]; % Default
end

train_cell_x = train_box_x + 0.05;
train_cell_y = train_box_y + train_box_height - 0.1;
for i = 1:min(6, length(train_channels))
    ch = train_channels(i);
    cell_x = train_cell_x + mod(i-1, 3) * (cell_width + cell_spacing);
    cell_y = train_cell_y - floor((i-1)/3) * (cell_height + cell_spacing);
    
    % Draw cell rectangle
    rectangle('Position', [cell_x, cell_y, cell_width, cell_height], ...
              'EdgeColor', 'b', 'LineWidth', 2, 'FaceColor', [0.8, 0.8, 1]);
    text(cell_x + cell_width/2, cell_y + cell_height/2, ...
         sprintf('Ch%d', ch), 'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% Draw cells in Test box (2 cells)
test_channels = [];
for f = 1:length(fold_definitions)
    fold_name = sprintf('Fold%d', f);
    if isfield(fold_results, fold_name)
        test_channels = fold_results.(fold_name).test_channels;
        break; % Use first fold as example
    end
end

if isempty(test_channels)
    test_channels = [15, 16]; % Default
end

test_cell_x = test_box_x + 0.1;
test_cell_y = train_box_y + train_box_height - 0.1;
for i = 1:min(2, length(test_channels))
    ch = test_channels(i);
    cell_x = test_cell_x + (i-1) * (cell_width + cell_spacing);
    cell_y = test_cell_y;
    
    % Draw cell rectangle
    rectangle('Position', [cell_x, cell_y, cell_width, cell_height], ...
              'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', [1, 0.8, 0.8]);
    text(cell_x + cell_width/2, cell_y + cell_height/2, ...
         sprintf('Ch%d', ch), 'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% Add title
title('4-Fold Cross Validation Data Split Strategy', 'FontSize', 18, 'FontWeight', 'bold');

% Add legend/description
desc_text = sprintf('Strategy: Each fold uses 2 cells for testing (1 from Group A + 1 from Group B)\n');
desc_text = [desc_text, sprintf('Group A (Mild): Ch09-12 (1C charge)\n')];
desc_text = [desc_text, sprintf('Group B (Harsh): Ch13-16 (3C charge)\n')];
desc_text = [desc_text, sprintf('Example shown: Fold 1 (Test: Ch09, Ch13)')];
text(0.5, 0.1, desc_text, 'FontSize', 12, 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'black');

savefig(fullfile(figuresDir, 'Data_Split_Diagram.fig'));
close(fig_diagram);
fprintf('  Saved: Data_Split_Diagram.fig\n');

fprintf('\n=== Data Splitting Complete ===\n');
