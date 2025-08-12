%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCIR Outlier Removal and Conditional Plots
% Remove outliers using 1σ criterion and create conditional plots
% Only create cluster-specific plots when same cluster exists in multiple years
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Load original data
dataPath = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_DCIR_Charge\AutoResults_Charge\all_chg_events_current_clustering_all_years.mat';

load(dataPath);

%% Create directories
saveDir = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_DCIR_Charge\AutoResults_Charge\Figure';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

figureDir = fullfile(saveDir, 'Figure');
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

%% Variables and Settings
yearList = {'2021', '2022', '2023'};
yearColors = [0 0.451 0.761; 0.937 0.753 0; 0.804 0.325 0.298]; % 2021=파랑, 2022=노랑, 2023=빨강
dcir_fields = {'DCIR_1s', 'DCIR_3s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s', 'DCIR_50s'};
dcir_labels = {'1s', '3s', '5s', '10s', '30s', '50s'};
dcir_diff_fields = {'DCIR_diff_5s_1s', 'DCIR_diff_10s_1s', 'DCIR_diff_30s_1s', 'DCIR_diff_50s_1s'};
dcir_diff_labels = {'5s-1s', '10s-1s', '30s-1s', '50s-1s'};

%% Get all cluster labels and generate display labels
all_labels = fieldnames(global_eventStruct);
display_labels = generate_display_labels(all_labels);

%% Helper Functions
function display_labels = generate_display_labels(all_labels)
    display_labels = {};
    for i = 1:length(all_labels)
        label = all_labels{i};
        if contains(label, 'cluster_')
            if contains(label, '0_')
                c_rate_str = regexp(label, 'cluster_0_(\d+)C', 'tokens');
                if ~isempty(c_rate_str)
                    c_rate_val = str2double(c_rate_str{1}{1}) / 100;
                    display_labels{i} = sprintf('%.2fC', c_rate_val);
                else
                    display_labels{i} = label;
                end
            else
                c_rate_str = regexp(label, 'cluster_(\d+)C', 'tokens');
                if ~isempty(c_rate_str)
                    c_rate_val = str2double(c_rate_str{1}{1});
                    display_labels{i} = sprintf('%.0fC', c_rate_val);
                else
                    display_labels{i} = label;
                end
            end
        else
            display_labels{i} = label;
        end
    end
end

function [year_data_all, year_labels, outlier_stats] = collect_year_data(cluster_data, field_name, yearList)
    year_data_all = {};
    year_labels = {};
    outlier_stats = struct('removed', 0, 'total', 0);
    
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_str = sprintf('year_%s', year);
        
        if isfield(cluster_data, year_str)
            year_data = cluster_data.(year_str);
            year_values = [];
            
            for event_name = fieldnames(year_data)'
                evt = year_data.(event_name{1});
                if isfield(evt, field_name)
                    if isstruct(evt.(field_name)) && isfield(evt.(field_name), 'val')
                        val = evt.(field_name).val;
                    else
                        val = evt.(field_name);
                    end
                    % Only include non-NaN values (outliers have been set to NaN)
                    if ~isnan(val)
                        year_values = [year_values, val];
                    end
                end
            end
            
            if ~isempty(year_values)
                year_data_all{end+1} = year_values;
                year_labels{end+1} = year;
            end
        end
    end
    
    % Calculate outlier statistics
    total_count = 0;
    removed_count = 0;
    
    % Count total and removed data points
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_str = sprintf('year_%s', year);
        
        if isfield(cluster_data, year_str)
            year_data = cluster_data.(year_str);
            
            for event_name = fieldnames(year_data)'
                evt = year_data.(event_name{1});
                if isfield(evt, field_name)
                    if isstruct(evt.(field_name)) && isfield(evt.(field_name), 'val')
                        val = evt.(field_name).val;
                    else
                        val = evt.(field_name);
                    end
                    if ~isnan(val)
                        total_count = total_count + 1;
                    else
                        removed_count = removed_count + 1;
                    end
                end
            end
        end
    end
    
    % Count filtered data points (after outlier removal)
    filtered_total = sum(cellfun(@length, year_data_all));
    
    outlier_stats.total = total_count + removed_count;
    outlier_stats.removed = removed_count;
end

function create_enhanced_boxplot(data, labels, title_str, ylabel_str, yearColors, save_path, outlier_stats)
    if isempty(data) || all(isnan(data))
        return;
    end
    
    valid_idx = ~isnan(data);
    filtered_data = data(valid_idx);
    filtered_labels = labels(valid_idx);
    
    if isempty(filtered_data) || ~all(isfinite(filtered_data))
        return;
    end
    
    unique_years = unique(filtered_labels);
    box_data = cell(1, length(unique_years));
    for y = 1:length(unique_years)
        year_idx = strcmp(filtered_labels, unique_years{y});
        box_data{y} = filtered_data(year_idx);
    end
    
    if ~all(cellfun(@(x) ~isempty(x) && all(isfinite(x)), box_data))
        return;
    end
    
    try
        bp = boxplot(filtered_data, filtered_labels, 'BoxStyle', 'outline', 'MedianStyle', 'line');
    catch ME
        fprintf('Error in boxplot: %s\n', ME.message);
        return;
    end
    
    % Customize appearance
    set(findobj(bp, 'type', 'line'), 'LineWidth', 2, 'Color', 'black');
    set(findobj(bp, 'type', 'line', 'Tag', 'Median'), 'LineWidth', 3, 'Color', 'red');
    set(gca, 'XTick', 1:length(unique_years), 'XTickLabel', unique_years);
    
    boxes = findobj(bp, 'type', 'patch');
    for b = 1:length(boxes)
        set(boxes(b), 'FaceColor', 'none', 'EdgeColor', 'black', 'LineWidth', 2);
    end
    
    % Add individual data points
    hold on;
    for y = 1:length(unique_years)
        year_idx = strcmp(labels, unique_years{y});
        year_data = data(year_idx);
        x_pos = y * ones(size(year_data));
        scatter(x_pos, year_data, 50, yearColors(y, :), 'o', 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    % Calculate ANOVA p-value if we have multiple years
    anova_p_value = NaN;
    if length(unique_years) >= 2
        try
            % Prepare data for ANOVA
            groups = [];
            values = [];
            for y = 1:length(unique_years)
                if ~isempty(box_data{y})
                    groups = [groups; repmat(y, length(box_data{y}), 1)];
                    values = [values; box_data{y}(:)];
                end
            end
            
                            if length(unique(groups)) >= 2 && length(values) >= 3
                    % Debug information
                    fprintf('\n=== ANOVA Debug - %s ===\n', title_str);
                    
                    % Print detailed group information
                    unique_groups = unique(groups);
                    fprintf('Groups found: ');
                    for g = 1:length(unique_groups)
                        group_idx = unique_groups(g);
                        group_values = values(groups == group_idx);
                        fprintf('Group %d (n=%d, μ=%.4f, σ=%.4f) ', group_idx, length(group_values), mean(group_values), std(group_values));
                    end
                    fprintf('\n');
                    
                    % Print all values for debugging
                    fprintf('All values: [');
                    for i = 1:length(values)
                        fprintf('%.4f', values(i));
                        if i < length(values)
                            fprintf(', ');
                        end
                    end
                    fprintf(']\n');
                    
                    fprintf('Group labels: [');
                    for i = 1:length(groups)
                        fprintf('%d', groups(i));
                        if i < length(groups)
                            fprintf(', ');
                        end
                    end
                    fprintf(']\n');
                    
                    % Perform one-way ANOVA
                    [~, tbl, ~] = anova1(values, groups, 'off');
                    anova_p_value = tbl{2, 6}; % p-value is in row 2, column 6
                    
                    % Print complete ANOVA table
                    fprintf('Complete ANOVA Table:\n');
                    fprintf('Source\t\tSS\t\t\tDF\t\tMS\t\t\tF\t\t\tProb>F\n');
                    fprintf('Groups\t\t%.6f\t%d\t\t%.6f\t%.6f\t%.6f\n', tbl{2,2}, tbl{2,3}, tbl{2,4}, tbl{2,5}, tbl{2,6});
                    fprintf('Error\t\t%.6f\t%d\t\t%.6f\n', tbl{3,2}, tbl{3,3}, tbl{3,4});
                    fprintf('Total\t\t%.6f\t%d\n', tbl{4,2}, tbl{4,3});
                    
                    fprintf('ANOVA result: p-value = %.6f\n', anova_p_value);
                    
                    % Manual calculation verification
                    fprintf('\n--- Manual Calculation Verification ---\n');
                    grand_mean = mean(values);
                    fprintf('Grand mean = %.4f\n', grand_mean);
                    
                    % Between-group sum of squares
                    ss_between = 0;
                    for g = 1:length(unique_groups)
                        group_idx = unique_groups(g);
                        group_values = values(groups == group_idx);
                        group_mean = mean(group_values);
                        n_group = length(group_values);
                        ss_between = ss_between + n_group * (group_mean - grand_mean)^2;
                    end
                    fprintf('SS_between = %.6f\n', ss_between);
                    
                    % Within-group sum of squares
                    ss_within = 0;
                    for g = 1:length(unique_groups)
                        group_idx = unique_groups(g);
                        group_values = values(groups == group_idx);
                        group_mean = mean(group_values);
                        for i = 1:length(group_values)
                            ss_within = ss_within + (group_values(i) - group_mean)^2;
                        end
                    end
                    fprintf('SS_within = %.6f\n', ss_within);
                    
                    % Degrees of freedom
                    df_between = length(unique_groups) - 1;
                    df_within = length(values) - length(unique_groups);
                    fprintf('DF_between = %d, DF_within = %d\n', df_between, df_within);
                    
                    % Mean squares
                    ms_between = ss_between / df_between;
                    ms_within = ss_within / df_within;
                    fprintf('MS_between = %.6f, MS_within = %.6f\n', ms_between, ms_within);
                    
                    % F-statistic
                    f_stat = ms_between / ms_within;
                    fprintf('F-statistic = %.6f\n', f_stat);
                    fprintf('=== End ANOVA Debug ===\n\n');
            end
        catch ME
            fprintf('ANOVA calculation error: %s\n', ME.message);
            anova_p_value = NaN;
        end
    end
    
    % Add statistics text with ANOVA p-value
    for y = 1:length(unique_years)
        if ~isempty(box_data{y})
            mean_val = mean(box_data{y});
            std_val = std(box_data{y});
            text(y, max(box_data{y}) + max(data) * 0.05, sprintf('μ=%.1f\nσ=%.1f', mean_val, std_val), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    
    % Add Kruskal-Wallis p-value to title if available
    if ~isnan(anova_p_value)
        if anova_p_value < 0.001
            p_value_str = 'p < 0.001';
        elseif anova_p_value < 0.01
            p_value_str = sprintf('p = %.3f', anova_p_value);
        else
            p_value_str = sprintf('p = %.3f', anova_p_value);
        end
        
                    % Add p-value to title
            title_with_p = sprintf('%s (ANOVA: %s)', title_str, p_value_str);
            title(title_with_p, 'FontSize', 14, 'FontWeight', 'bold');
    else
        title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % Add outlier statistics subtitle
    if nargin >= 7 && ~isempty(outlier_stats)
        subtitle_str = sprintf('Outliers removed: %d/%d', outlier_stats.removed, outlier_stats.total);
        subtitle(subtitle_str, 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
    end
    
    ylabel(ylabel_str, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
    grid on;
    ylim([0, max(data) * 1.15]);
end

% function create_dcir_histogram(year_data_all, year_labels, title_str, yearColors, save_path)
%     if isempty(year_data_all)
%         return;
%     end
% 
%     all_values = [year_data_all{:}];
%     if isempty(all_values) || ~all(isfinite(all_values))
%         return;
%     end
% 
%     min_val = min(all_values);
%     max_val = max(all_values);
%     bin_edges = linspace(min_val, max_val, 20);
% 
%     for y = 1:length(year_data_all)
%         [counts, edges] = histcounts(year_data_all{y}, bin_edges);
%         bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
% 
%         bar(bin_centers, counts, 1, 'FaceColor', yearColors(y, :), 'FaceAlpha', 0.7, 'EdgeColor', 'none', ...
%             'DisplayName', sprintf('%s (Mean: %.3f, Std: %.3f)', year_labels{y}, mean(year_data_all{y}), std(year_data_all{y})));
%         hold on;
% 
%         % Add normal distribution curve
%         mu = mean(year_data_all{y});
%         sigma = std(year_data_all{y});
%         x_norm = linspace(min_val, max_val, 100);
%         y_norm = normpdf(x_norm, mu, sigma) * length(year_data_all{y}) * (max_val - min_val) / 20;
%         plot(x_norm, y_norm, 'Color', yearColors(y, :), 'LineWidth', 2, 'DisplayName', sprintf('Normal dist. %s', year_labels{y}));
%     end
% 
%     title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
%     ylabel('Frequency', 'FontSize', 12, 'FontWeight', 'bold');
%     xlabel('DCIR [mΩ]', 'FontSize', 12, 'FontWeight', 'bold');
%     set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
%     grid on;
%     legend('Location', 'best');
% end

function create_individual_events_plot(year_data, field_name, field_label, year, cluster_label, yearColors, year_idx, figureDir)
    event_names = fieldnames(year_data);
    if isempty(event_names)
        return;
    end
    
    figure('Name', sprintf('DCIR Individual Events - %s %s', year, cluster_label), 'Position', [100, 100, 1200, 800]);
    
    for k = 1:length(field_name)
        subplot(2, 3, k); hold on;
        
        event_numbers = [];
        values = [];
        outlier_flags = [];
        
        % Collect original data for statistics
        original_values = [];
        for event_idx = 1:length(event_names)
            event_name = event_names{event_idx};
            evt = year_data.(event_name);
            if isfield(evt, field_name{k})
                if isstruct(evt.(field_name{k})) && isfield(evt.(field_name{k}), 'val')
                    val = evt.(field_name{k}).val;
                else
                    val = evt.(field_name{k});
                end
                if ~isnan(val)
                    original_values = [original_values, val];
                end
            end
        end
        
        % Calculate statistics and plot lines
        if ~isempty(original_values)
            mean_val = mean(original_values);
            std_val = std(original_values);
            lower_bound = mean_val - std_val;
            upper_bound = mean_val + std_val;
            
            x_range = [1, length(event_names)];
            plot(x_range, [mean_val, mean_val], '--', 'Color', 'red', 'LineWidth', 2, 'DisplayName', sprintf('Mean: %.2f', mean_val));
            plot(x_range, [lower_bound, lower_bound], '--', 'Color', 'black', 'LineWidth', 1, 'DisplayName', sprintf('-1σ: %.2f', lower_bound));
            plot(x_range, [upper_bound, upper_bound], '--', 'Color', 'black', 'LineWidth', 1, 'DisplayName', sprintf('+1σ: %.2f', upper_bound));
        end
        
        % Plot individual events
        for event_idx = 1:length(event_names)
            event_name = event_names{event_idx};
            evt = year_data.(event_name);
            if isfield(evt, field_name{k})
                if isstruct(evt.(field_name{k})) && isfield(evt.(field_name{k}), 'val')
                    val = evt.(field_name{k}).val;
                else
                    val = evt.(field_name{k});
                end
                event_numbers = [event_numbers, event_idx];
                
                if ~isnan(val)
                    values = [values, val];
                    if val < lower_bound || val > upper_bound
                        outlier_flags = [outlier_flags, 1];
                    else
                        outlier_flags = [outlier_flags, 0];
                    end
                else
                    values = [values, NaN];
                    outlier_flags = [outlier_flags, 2];
                end
            end
        end
        
        if ~isempty(event_numbers)
            % Plot normal data points
            normal_idx = outlier_flags == 0;
            if any(normal_idx)
                plot(event_numbers(normal_idx), values(normal_idx), 'o', 'MarkerSize', 8, 'MarkerFaceColor', yearColors(year_idx, :), 'MarkerEdgeColor', yearColors(year_idx, :), 'DisplayName', year);
            end
            
            % Plot outliers as red X
            outlier_idx = outlier_flags == 1;
            if any(outlier_idx)
                plot(event_numbers(outlier_idx), values(outlier_idx), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Outliers');
            end
            
            % Plot NaN values as black X
            nan_idx = outlier_flags == 2;
            if any(nan_idx)
                plot(event_numbers(nan_idx), zeros(1, sum(nan_idx)), 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'NaN');
            end
            
            title(sprintf('DCIR %s - %s %s (%d events)', field_label{k}, year, cluster_label, length(event_numbers)), 'FontSize', 14, 'FontWeight', 'bold');
            xlabel('Event Number', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('DCIR [mΩ]', 'FontSize', 12, 'FontWeight', 'bold');
            set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
            grid on;
            legend('Location', 'best');
        end
    end
    
    sgtitle(sprintf('DCIR Individual Events - %s %s (Outlier Removed)', year, cluster_label), 'FontSize', 16, 'FontWeight', 'bold');
    saveas(gcf, fullfile(figureDir, sprintf('fig_%s_DCIR_Individual_%s.fig', year, cluster_label)));
end

%% Create outlier-removed structure and remove outliers
outlier_removed_struct = global_eventStruct;

fprintf('Removing outliers using 1σ criterion with year-specific baselines...\n');
for i = 1:length(all_labels)
    label = all_labels{i};
    cluster_data = global_eventStruct.(label);
    
    if isempty(cluster_data)
        continue;
    end
    
    % Remove outliers year by year (independent for each year)
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_str = sprintf('year_%s', year);
        
        if ~isfield(cluster_data, year_str)
            continue;
        end
        
        year_data = cluster_data.(year_str);
        event_names = fieldnames(year_data);
        
        if isempty(event_names)
            continue;
        end
        
        % Collect all DCIR values for this year only
        all_dcir_values = struct();
        for k = 1:length(dcir_fields)
            all_dcir_values.(dcir_fields{k}) = [];
        end
        
        for event_idx = 1:length(event_names)
            event_name = event_names{event_idx};
            evt = year_data.(event_name);
            
            for k = 1:length(dcir_fields)
                if isfield(evt, dcir_fields{k}) && isfield(evt.(dcir_fields{k}), 'val')
                    val = evt.(dcir_fields{k}).val;
                    if ~isnan(val)
                        all_dcir_values.(dcir_fields{k}) = [all_dcir_values.(dcir_fields{k}), val];
                    end
                end
            end
        end
        
        % Calculate year-specific statistics and remove outliers
        for k = 1:length(dcir_fields)
            field_name = dcir_fields{k};
            values = all_dcir_values.(field_name);
            
            if length(values) >= 3
                year_mean = mean(values);
                year_std = std(values);
                year_lower_bound = year_mean - year_std;
                year_upper_bound = year_mean + year_std;
                
                fprintf('Cluster %s, Year %s, Field %s: μ=%.2f, σ=%.2f, bounds=[%.2f, %.2f]\n', ...
                    label, year, field_name, year_mean, year_std, year_lower_bound, year_upper_bound);
                
                for event_idx = 1:length(event_names)
                    event_name = event_names{event_idx};
                    evt = year_data.(event_name);
                    
                    if isfield(evt, field_name) && isfield(evt.(field_name), 'val')
                        val = evt.(field_name).val;
                        if ~isnan(val) && (val < year_lower_bound || val > year_upper_bound)
                            outlier_removed_struct.(label).(year_str).(event_name).(field_name).val = NaN;
                            fprintf('Outlier removed: %s %s %s %s = %.2f (year-specific)\n', label, year, event_name, field_name, val);
                        end
                    end
                end
            end
        end
    end
end

%% Save outlier-removed data
save(fullfile(saveDir, 'all_chg_events_outlier_removed.mat'), 'outlier_removed_struct');
fprintf('Outlier-removed data saved.\n');

%% Find clusters that exist in multiple years
multi_year_clusters = {};
for i = 1:length(all_labels)
    label = all_labels{i};
    cluster_data = outlier_removed_struct.(label);
    
    if isempty(cluster_data)
        continue;
    end
    
    year_count = 0;
    for year_idx = 1:length(yearList)
        year_str = sprintf('year_%s', yearList{year_idx});
        if isfield(cluster_data, year_str)
            year_data = cluster_data.(year_str);
            if ~isempty(fieldnames(year_data))
                year_count = year_count + 1;
            end
        end
    end
    
    if year_count >= 2
        multi_year_clusters{end+1} = label;
    end
end

fprintf('Found %d clusters that exist in multiple years.\n', length(multi_year_clusters));

%% Create plots for multi-year clusters only
for cluster_idx = 1:length(multi_year_clusters)
    label = multi_year_clusters{cluster_idx};
    cluster_data = outlier_removed_struct.(label);
    display_label = display_labels{strcmp(all_labels, label)};
    
    % Create DCIR box plots
    figure('Name', sprintf('DCIR Box Plot - %s (Outlier Removed)', display_label), 'Position', [100, 100, 1200, 800]);
    for k = 1:length(dcir_fields)
        subplot(2, 3, k);
        [year_data_all, year_labels, outlier_stats] = collect_year_data(cluster_data, dcir_fields{k}, yearList);
        
        if ~isempty(year_data_all)
            all_data = [year_data_all{:}];
            group_labels = [];
            for y = 1:length(year_labels)
                group_labels = [group_labels, repmat({year_labels{y}}, 1, length(year_data_all{y}))];
            end
            
            create_enhanced_boxplot(all_data, group_labels, sprintf('DCIR %s - %s', dcir_labels{k}, display_label), 'DCIR [mΩ]', yearColors, '', outlier_stats);
        end
    end
    sgtitle(sprintf('DCIR Box Plots - %s (Outlier Removed)', display_label), 'FontSize', 16, 'FontWeight', 'bold');
    saveas(gcf, fullfile(figureDir, sprintf('fig_DCIR_BoxPlot_%s_OutlierRemoved.fig', label)));
    
    % Create DCIR Difference box plots
    figure('Name', sprintf('DCIR Diff Box Plot - %s (Outlier Removed)', display_label), 'Position', [100, 100, 1000, 600]);
    for k = 1:length(dcir_diff_fields)
        subplot(2, 2, k);
        [year_data_all, year_labels, outlier_stats] = collect_year_data(cluster_data, dcir_diff_fields{k}, yearList);
        
        if ~isempty(year_data_all)
            all_data = [year_data_all{:}];
            group_labels = [];
            for y = 1:length(year_labels)
                group_labels = [group_labels, repmat({year_labels{y}}, 1, length(year_data_all{y}))];
            end
            
            create_enhanced_boxplot(all_data, group_labels, sprintf('DCIR Diff %s - %s', dcir_diff_labels{k}, display_label), 'DCIR Difference [mΩ]', yearColors, '', outlier_stats);
        end
    end
    sgtitle(sprintf('DCIR Diff Box Plots - %s (Outlier Removed)', display_label), 'FontSize', 16, 'FontWeight', 'bold');
    saveas(gcf, fullfile(figureDir, sprintf('fig_DCIR_Diff_BoxPlot_%s_OutlierRemoved.fig', label)));
    
    % % Create DCIR distribution histogram
    % figure('Name', sprintf('DCIR Distribution Histogram - %s', display_label), 'Position', [100, 100, 1200, 800]);
    % for k = 1:length(dcir_fields)
    %     subplot(2, 3, k);
    %     [year_data_all, year_labels, outlier_stats] = collect_year_data(cluster_data, dcir_fields{k}, yearList);
    %     create_dcir_histogram(year_data_all, year_labels, sprintf('DCIR Distribution %s - %s', dcir_labels{k}, display_label), yearColors, '');
    % end
    % sgtitle(sprintf('DCIR Distribution Histograms - %s (Outlier Removed)', display_label), 'FontSize', 16, 'FontWeight', 'bold');
    % saveas(gcf, fullfile(figureDir, sprintf('fig_DCIR_Distribution_Histogram_%s_OutlierRemoved.fig', label)));
end

%% Create individual event plots for all clusters
fprintf('Creating individual event plots...\n');
for i = 1:length(all_labels)
    label = all_labels{i};
    cluster_data = outlier_removed_struct.(label);
    display_label = display_labels{i};
    
    if isempty(cluster_data)
        continue;
    end
    
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_str = sprintf('year_%s', year);
        
        if isfield(cluster_data, year_str)
            year_data = cluster_data.(year_str);
            create_individual_events_plot(year_data, dcir_fields, dcir_labels, year, display_label, yearColors, year_idx, figureDir);
        end
    end
end

fprintf('Outlier removal and conditional plotting completed successfully.\n');
fprintf('Results saved in: %s\n', saveDir); 