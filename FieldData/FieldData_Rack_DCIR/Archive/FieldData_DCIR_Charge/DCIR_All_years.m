%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCIR All Years Visualization
% Multi-year DCIR analysis and visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Load data
dataPath = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_DCIR_Charge\AutoResults_Charge\all_chg_events_current_clustering_all_years.mat';

load(dataPath);

%% Create Figure directories
saveDir = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_DCIR_Charge\AutoResults_Charge\Figure';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Create year-specific directories for individual event plots
for year = {'2021', '2022', '2023'}
    yearDir = fullfile(saveDir, year{1});
    if ~exist(yearDir, 'dir')
        mkdir(yearDir);
    end
end

%% Variables
C_nom = 1024; % Ah
yearList = {'2021', '2022', '2023'};
yearColors = [0 0.451 0.761; 0.937 0.753 0; 0.804 0.325 0.298]; % 2021=파랑(#0073C2), 2022=노랑(#EFC000), 2023=빨강(#CD534C) (RGB)



%% Get all cluster labels
all_labels = fieldnames(global_eventStruct);

%% Generate display labels
display_labels = {};
for i = 1:length(all_labels)
    label = all_labels{i};
    if contains(label, 'cluster_')
        current_val = regexp(label, 'cluster_(\d+)A', 'tokens');
        if ~isempty(current_val)
            current_num = str2double(current_val{1}{1});
            c_rate = current_num / C_nom;
            display_labels{i} = sprintf('%.2fC', c_rate);
        else
            display_labels{i} = label;
        end
    else
        display_labels{i} = label;
    end
end

%% Combined Current and Voltage vs Time plots (Year-specific)
for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    year_str = sprintf('year_%s', year);
    
    % Create figure for this year
    figure('Name', sprintf('Current & Voltage vs Time - %s', year), 'Position', [100, 100, 1000, 600]);
    
    subplot_idx = 1;
    for i = 1:length(all_labels)
        label = all_labels{i};
        cluster_data = global_eventStruct.(label);
        
        if isempty(cluster_data) || ~isfield(cluster_data, year_str)
            continue;
        end
        
        year_data = cluster_data.(year_str);
        event_names = fieldnames(year_data);
        num_events = length(event_names);
        
        if num_events == 0
            continue;
        end
        
        % Current vs Time subplot
        subplot(2, 3, subplot_idx); hold on;
        for event_name = event_names'
            evt = year_data.(event_name{1});
            if isfield(evt, 't_seq') && isfield(evt, 'I_seq')
                t_data = evt.t_seq;
                I_data = evt.I_seq;
                t_seconds = seconds(t_data - t_data(1));
                plot(t_seconds, I_data, 'LineWidth', 1.5);
            end
        end
        title(sprintf('Current vs Time: %s (%d events)', display_labels{i}, num_events));
        xlabel('Time [s]');
        ylabel('Current [A]');
        grid on;
        
        % Voltage vs Time subplot
        subplot(2, 3, subplot_idx + 3); hold on;
        for event_name = event_names'
            evt = year_data.(event_name{1});
            if isfield(evt, 't_seq') && isfield(evt, 'V_seq')
                t_data = evt.t_seq;
                V_data = evt.V_seq;
                t_seconds = seconds(t_data - t_data(1));
                plot(t_seconds, V_data, 'LineWidth', 1.5);
            end
        end
        title(sprintf('Voltage vs Time: %s (%d events)', display_labels{i}, num_events));
        xlabel('Time [s]');
        ylabel('Voltage [V]');
        grid on;
        
        subplot_idx = subplot_idx + 1;
        
        % Only show up to 3 clusters per figure
        if subplot_idx > 3
            break;
        end
    end
    
    sgtitle(sprintf('Current & Voltage vs Time - %s', year));
    saveas(gcf, fullfile(saveDir, sprintf('fig_%s_Current_Voltage_vs_Time.fig', year)));
end

%% DCIR Individual Event Plots and Box Plots
dcir_fields = {'DCIR_1s', 'DCIR_3s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s', 'DCIR_50s'};
dcir_labels = {'1s', '3s', '5s', '10s', '30s', '50s'};

% Group similar clusters (within 10A difference) - avoid overlapping groups
cluster_groups = {};
used_clusters = [];

% First, extract all cluster values and sort them
cluster_values = [];
cluster_names = {};
for i = 1:length(all_labels)
    val = regexp(all_labels{i}, 'cluster_(\d+)A', 'tokens');
    if ~isempty(val)
        cluster_values = [cluster_values, str2double(val{1}{1})];
        cluster_names{end+1} = all_labels{i};
    end
end

% Sort by cluster values
[sorted_values, sort_idx] = sort(cluster_values);
sorted_names = cluster_names(sort_idx);

% Group consecutive clusters within 10A difference
i = 1;
while i <= length(sorted_names)
    if any(strcmp(used_clusters, sorted_names{i}))
        i = i + 1;
        continue;
    end
    
    current_group = {sorted_names{i}};
    used_clusters = [used_clusters, sorted_names{i}];
    current_val = sorted_values(i);
    
    % Find consecutive clusters within 10A
    j = i + 1;
    while j <= length(sorted_names) && abs(sorted_values(j) - current_val) <= 10
        if ~any(strcmp(used_clusters, sorted_names{j}))
            current_group{end+1} = sorted_names{j};
            used_clusters = [used_clusters, sorted_names{j}];
        end
        j = j + 1;
    end
    
    if length(current_group) >= 2
        cluster_groups{end+1} = current_group;
    end
    
    i = j;
end

% Create individual event plots for each cluster (saved in year-specific folders)
for i = 1:length(all_labels)
    label = all_labels{i};
    cluster_data = global_eventStruct.(label);
    
    if isempty(cluster_data)
        continue;
    end
    
    % Create individual event plots for each year
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_str = sprintf('year_%s', year);
        
        if isfield(cluster_data, year_str)
            year_data = cluster_data.(year_str);
            event_names = fieldnames(year_data);
            
            if ~isempty(event_names)
                % Create individual event plot for this year and cluster
                figure('Name', sprintf('DCIR Individual Events - %s %s', year, display_labels{i}), 'Position', [100, 100, 1000, 600]);
                
                for k = 1:length(dcir_fields)
                    subplot(2, 3, k); hold on;
                    
                    event_numbers = [];
                    dcir_values = [];
                    
                    for event_idx = 1:length(event_names)
                        event_name = event_names{event_idx};
                        evt = year_data.(event_name);
                        if isfield(evt, dcir_fields{k}) && isfield(evt.(dcir_fields{k}), 'val')
                            val = evt.(dcir_fields{k}).val;
                            % Include all events, even NaN values
                            event_numbers = [event_numbers, event_idx];
                            if ~isnan(val)
                                dcir_values = [dcir_values, val];
                            else
                                dcir_values = [dcir_values, NaN];
                            end
                        end
                    end
                    
                    if ~isempty(event_numbers)
                        % Plot valid data points
                        valid_idx = ~isnan(dcir_values);
                        if any(valid_idx)
                            plot(event_numbers(valid_idx), dcir_values(valid_idx), 'o', 'MarkerSize', 5, 'MarkerFaceColor', yearColors(year_idx, :));
                        end
                        
                        % Plot NaN values as red X markers
                        nan_idx = isnan(dcir_values);
                        if any(nan_idx)
                            plot(event_numbers(nan_idx), zeros(1, sum(nan_idx)), 'rx', 'MarkerSize', 6);
                        end
                        
                        title(sprintf('DCIR %s - %s %s (%d events)', dcir_labels{k}, year, display_labels{i}, length(event_numbers)), 'FontSize', 14, 'FontWeight', 'bold');
                        xlabel('Event Number', 'FontSize', 12, 'FontWeight', 'bold');
                        ylabel('DCIR [mΩ]', 'FontSize', 12, 'FontWeight', 'bold');
                        set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
                        grid on;
                    end
                end
                
                sgtitle(sprintf('DCIR Individual Events - %s %s', year, display_labels{i}), 'FontSize', 16, 'FontWeight', 'bold');
                saveas(gcf, fullfile(saveDir, year, sprintf('fig_%s_DCIR_Individual_%s.fig', year, label)));
            end
        end
    end
end

% Create box plots for similar cluster groups
for group_idx = 1:length(cluster_groups)
    group = cluster_groups{group_idx};
    
    % Find the representative cluster name (middle value)
    group_values = [];
    for i = 1:length(group)
        val = regexp(group{i}, 'cluster_(\d+)A', 'tokens');
        if ~isempty(val)
            group_values = [group_values, str2double(val{1}{1})];
        end
    end
    [~, mid_idx] = min(abs(group_values - mean(group_values)));
    representative_cluster = group{mid_idx};
    
    % Create DCIR box plot figure for similar clusters
    figure('Name', sprintf('DCIR Box Plot - Similar Clusters (%.0fA)', mean(group_values)), 'Position', [100, 100, 1000, 600]);
    
    for k = 1:length(dcir_fields)
        subplot(2, 3, k); hold on;
        
        % Collect data for each year and cluster in the group
        all_data = [];
        group_labels = {};
        
        for year_idx = 1:length(yearList)
            year = yearList{year_idx};
            year_str = sprintf('year_%s', year);
            
            for cluster_idx = 1:length(group)
                cluster_label = group{cluster_idx};
                if isfield(global_eventStruct, cluster_label) && isfield(global_eventStruct.(cluster_label), year_str)
                    cluster_data = global_eventStruct.(cluster_label).(year_str);
                    year_values = [];
                    
                    for event_name = fieldnames(cluster_data)'
                        evt = cluster_data.(event_name{1});
                        if isfield(evt, dcir_fields{k}) && isfield(evt.(dcir_fields{k}), 'val')
                            val = evt.(dcir_fields{k}).val;
                            if ~isnan(val)
                                year_values = [year_values, val];
                            end
                        end
                    end
                    
                    if ~isempty(year_values)
                        all_data = [all_data, year_values];
                        group_labels = [group_labels, repmat({year}, 1, length(year_values))];
                    end
                end
            end
        end
        
        % Create enhanced box plot
        if ~isempty(all_data)
            % Create custom box plot with better styling
            unique_years = unique(group_labels);
            box_data = cell(1, length(unique_years));
            for y = 1:length(unique_years)
                year_idx = strcmp(group_labels, unique_years{y});
                box_data{y} = all_data(year_idx);
            end
            
            % Create box plot with custom properties
            bp = boxplot(all_data, group_labels, 'Colors', yearColors(1:length(unique_years),:), ...
                'BoxStyle', 'filled', 'MedianStyle', 'line', 'OutlierSize', 10, 'Width', 0.8);
            
            % Customize box plot appearance
            set(findobj(bp, 'type', 'line'), 'LineWidth', 2.5);
            set(findobj(bp, 'type', 'line', 'Tag', 'Median'), 'LineWidth', 4);
            set(findobj(bp, 'type', 'line', 'Tag', 'Outliers'), 'MarkerSize', 12);
            
            % Make boxes more visible
            boxes = findobj(bp, 'type', 'patch');
            for b = 1:length(boxes)
                set(boxes(b), 'FaceAlpha', 0.8, 'EdgeColor', 'black', 'LineWidth', 2);
            end
            
            % Add statistics text
            for y = 1:length(unique_years)
                if ~isempty(box_data{y})
                    mean_val = mean(box_data{y});
                    std_val = std(box_data{y});
                    text(y, max(box_data{y}) + max(all_data) * 0.05, sprintf('μ=%.1f\nσ=%.1f', mean_val, std_val), ...
                        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
                end
            end
            
            title(sprintf('DCIR %s - Similar Clusters (%.0fA)', dcir_labels{k}, mean(group_values)), 'FontSize', 14, 'FontWeight', 'bold');
            ylabel('DCIR [mΩ]', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
            set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
            grid on;
            ylim([0, max(all_data) * 1.15]);
        end
    end
    
    sgtitle(sprintf('DCIR Box Plots - Similar Clusters (%.0fA)', mean(group_values)), 'FontSize', 16, 'FontWeight', 'bold');
    saveas(gcf, fullfile(saveDir, sprintf('fig_DCIR_BoxPlot_Similar_%.0fA.fig', mean(group_values))));
end

%% DCIR Difference Individual Event Plots and Box Plots
dcir_diff_fields = {'DCIR_diff_5s_1s', 'DCIR_diff_10s_1s', 'DCIR_diff_30s_1s', 'DCIR_diff_50s_1s'};
dcir_diff_labels = {'5s-1s', '10s-1s', '30s-1s', '50s-1s'};

% Create individual event plots for each cluster (saved in year-specific folders)
for i = 1:length(all_labels)
    label = all_labels{i};
    cluster_data = global_eventStruct.(label);
    
    if isempty(cluster_data)
        continue;
    end
    
    % Create individual event plots for each year
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_str = sprintf('year_%s', year);
        
        if isfield(cluster_data, year_str)
            year_data = cluster_data.(year_str);
            event_names = fieldnames(year_data);
            
            if ~isempty(event_names)
                % Create individual event plot for this year and cluster
                figure('Name', sprintf('DCIR Diff Individual Events - %s %s', year, display_labels{i}), 'Position', [100, 100, 1000, 600]);
                
                for k = 1:length(dcir_diff_fields)
                    subplot(2, 2, k); hold on;
                    
                    event_numbers = [];
                    dcir_diff_values = [];
                    
                    for event_idx = 1:length(event_names)
                        event_name = event_names{event_idx};
                        evt = year_data.(event_name);
                        if isfield(evt, dcir_diff_fields{k})
                            val = evt.(dcir_diff_fields{k});
                            % Include all events, even NaN values
                            event_numbers = [event_numbers, event_idx];
                            if ~isnan(val)
                                dcir_diff_values = [dcir_diff_values, val];
                            else
                                dcir_diff_values = [dcir_diff_values, NaN];
                            end
                        end
                    end
                    
                    if ~isempty(event_numbers)
                        % Plot valid data points
                        valid_idx = ~isnan(dcir_diff_values);
                        if any(valid_idx)
                            plot(event_numbers(valid_idx), dcir_diff_values(valid_idx), 'o', 'MarkerSize', 5, 'MarkerFaceColor', yearColors(year_idx, :));
                        end
                        
                        % Plot NaN values as red X markers
                        nan_idx = isnan(dcir_diff_values);
                        if any(nan_idx)
                            plot(event_numbers(nan_idx), zeros(1, sum(nan_idx)), 'rx', 'MarkerSize', 6);
                        end
                        
                        title(sprintf('DCIR Diff %s - %s %s (%d events)', dcir_diff_labels{k}, year, display_labels{i}, length(event_numbers)), 'FontSize', 14, 'FontWeight', 'bold');
                        xlabel('Event Number', 'FontSize', 12, 'FontWeight', 'bold');
                        ylabel('DCIR Difference [mΩ]', 'FontSize', 12, 'FontWeight', 'bold');
                        set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
                        grid on;
                    end
                end
                
                sgtitle(sprintf('DCIR Diff Individual Events - %s %s', year, display_labels{i}), 'FontSize', 16, 'FontWeight', 'bold');
                saveas(gcf, fullfile(saveDir, year, sprintf('fig_%s_DCIR_Diff_Individual_%s.fig', year, label)));
            end
        end
    end
end

% Create box plots for similar cluster groups (DCIR Difference)
for group_idx = 1:length(cluster_groups)
    group = cluster_groups{group_idx};
    
    % Find the representative cluster name (middle value)
    group_values = [];
    for i = 1:length(group)
        val = regexp(group{i}, 'cluster_(\d+)A', 'tokens');
        if ~isempty(val)
            group_values = [group_values, str2double(val{1}{1})];
        end
    end
    [~, mid_idx] = min(abs(group_values - mean(group_values)));
    representative_cluster = group{mid_idx};
    
    % Create DCIR Difference box plot figure for similar clusters
    figure('Name', sprintf('DCIR Diff Box Plot - Similar Clusters (%.0fA)', mean(group_values)), 'Position', [100, 100, 1000, 600]);
    
    for k = 1:length(dcir_diff_fields)
        subplot(2, 2, k); hold on;
        
        % Collect data for each year and cluster in the group
        all_data = [];
        group_labels = {};
        
        for year_idx = 1:length(yearList)
            year = yearList{year_idx};
            year_str = sprintf('year_%s', year);
            
            for cluster_idx = 1:length(group)
                cluster_label = group{cluster_idx};
                if isfield(global_eventStruct, cluster_label) && isfield(global_eventStruct.(cluster_label), year_str)
                    cluster_data = global_eventStruct.(cluster_label).(year_str);
                    year_values = [];
                    
                    for event_name = fieldnames(cluster_data)'
                        evt = cluster_data.(event_name{1});
                        if isfield(evt, dcir_diff_fields{k})
                            val = evt.(dcir_diff_fields{k});
                            if ~isnan(val)
                                year_values = [year_values, val];
                            end
                        end
                    end
                    
                    if ~isempty(year_values)
                        all_data = [all_data, year_values];
                        group_labels = [group_labels, repmat({year}, 1, length(year_values))];
                    end
                end
            end
        end
        
        % Create enhanced box plot
        if ~isempty(all_data)
            % Create custom box plot with better styling
            unique_years = unique(group_labels);
            box_data = cell(1, length(unique_years));
            for y = 1:length(unique_years)
                year_idx = strcmp(group_labels, unique_years{y});
                box_data{y} = all_data(year_idx);
            end
            
            % Create box plot with custom properties
            bp = boxplot(all_data, group_labels, 'Colors', yearColors(1:length(unique_years),:), ...
                'BoxStyle', 'filled', 'MedianStyle', 'line', 'OutlierSize', 10, 'Width', 0.8);
            
            % Customize box plot appearance
            set(findobj(bp, 'type', 'line'), 'LineWidth', 2.5);
            set(findobj(bp, 'type', 'line', 'Tag', 'Median'), 'LineWidth', 4);
            set(findobj(bp, 'type', 'line', 'Tag', 'Outliers'), 'MarkerSize', 12);
            
            % Make boxes more visible
            boxes = findobj(bp, 'type', 'patch');
            for b = 1:length(boxes)
                set(boxes(b), 'FaceAlpha', 0.8, 'EdgeColor', 'black', 'LineWidth', 2);
            end
            
            % Add statistics text
            for y = 1:length(unique_years)
                if ~isempty(box_data{y})
                    mean_val = mean(box_data{y});
                    std_val = std(box_data{y});
                    text(y, max(box_data{y}) + max(all_data) * 0.05, sprintf('μ=%.1f\nσ=%.1f', mean_val, std_val), ...
                        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
                end
            end
            
            title(sprintf('DCIR Diff %s - Similar Clusters (%.0fA)', dcir_diff_labels{k}, mean(group_values)), 'FontSize', 14, 'FontWeight', 'bold');
            ylabel('DCIR Difference [mΩ]', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
            set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
            grid on;
            ylim([min(all_data) * 0.9, max(all_data) * 1.15]);
        end
    end
    
    sgtitle(sprintf('DCIR Diff Box Plots - Similar Clusters (%.0fA)', mean(group_values)), 'FontSize', 16, 'FontWeight', 'bold');
    saveas(gcf, fullfile(saveDir, sprintf('fig_DCIR_Diff_BoxPlot_Similar_%.0fA.fig', mean(group_values))));
end

fprintf('All years visualization completed successfully.\n'); 

%% Year-specific Current and Voltage vs Time plots (All events combined)
fprintf('Creating year-specific Current and Voltage vs Time plots...\n');

for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    year_str = sprintf('year_%s', year);
    
    % Collect all events for this year
    all_t_current = {};
    all_I_current = {};
    all_t_voltage = {};
    all_V_voltage = {};
    
    % Collect data from all clusters for this year
    for i = 1:length(all_labels)
        label = all_labels{i};
        cluster_data = global_eventStruct.(label);
        
        if isempty(cluster_data) || ~isfield(cluster_data, year_str)
            continue;
        end
        
        year_data = cluster_data.(year_str);
        event_names = fieldnames(year_data);
        
        for event_name = event_names'
            evt = year_data.(event_name{1});
            if isfield(evt, 't_seq') && isfield(evt, 'I_seq') && isfield(evt, 'V_seq')
                t_data = evt.t_seq;
                I_data = evt.I_seq;
                V_data = evt.V_seq;
                
                % Normalize time to start from 0
                t_normalized = seconds(t_data - t_data(1));
                
                all_t_current{end+1} = t_normalized;
                all_I_current{end+1} = I_data;
                all_t_voltage{end+1} = t_normalized;
                all_V_voltage{end+1} = V_data;
            end
        end
    end
    
    % Create figure for this year
    figure('Name', sprintf('Current and Voltage vs Time - %s', year), 'Position', [100, 100, 1200, 600]);
    
    % Current vs Time subplot (1,2,1)
    subplot(1, 2, 1); hold on;
    for i = 1:length(all_t_current)
        plot(all_t_current{i}, all_I_current{i}, 'LineWidth', 1.5, 'LineStyle', '-');
    end
    title(sprintf('Current vs Time - %s (All Events)', year), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Current [A]', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
    grid on;
    
    % Voltage vs Time subplot (1,2,2)
    subplot(1, 2, 2); hold on;
    for i = 1:length(all_t_voltage)
        plot(all_t_voltage{i}, all_V_voltage{i}, 'LineWidth', 1.5, 'LineStyle', '-');
    end
    title(sprintf('Voltage vs Time - %s (All Events)', year), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Voltage [V]', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
    grid on;
    
    sgtitle(sprintf('Current and Voltage vs Time - %s (All Events)', year), 'FontSize', 16, 'FontWeight', 'bold');
    saveas(gcf, fullfile(saveDir, sprintf('fig_%s_Current_Voltage_vs_Time_AllEvents.fig', year)));
    
    fprintf('Created figure for year %s with %d events\n', year, length(all_t_current));
end

fprintf('Year-specific Current and Voltage vs Time plots completed successfully.\n'); 