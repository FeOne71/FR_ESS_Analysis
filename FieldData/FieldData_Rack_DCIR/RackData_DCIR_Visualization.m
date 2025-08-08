%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rack Data DCIR Visualization Script
% Basic MATLAB functions only - no fancy graphics
% Visualizes DCIR trends, current distributions, and rack comparisons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Load data
dataPath = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_Rack_DCIR\DCIR_Charge_0728\all_chg_events_current_clustering_all_years.mat';
load(dataPath);

%% Parameters
rackNames = {'Rack01', 'Rack02', 'Rack03', 'Rack04', 'Rack05', 'Rack06', 'Rack07', 'Rack08'};
years = {'2021', '2022', '2023'};
dcir_times = {'1s', '5s', '10s', '30s', '50s'};


%% Figure 1: Current Distribution by Rack
figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);

for rack_idx = 1:length(rackNames)
    rackName = rackNames{rack_idx};
    
    subplot(2, 4, rack_idx);
    
    % Collect all current data for this rack
    all_currents = [];
    cluster_labels = {};
    
    clusters = fieldnames(global_eventStruct.(rackName));
    
    for cluster_idx = 1:length(clusters)
        cluster_name = clusters{cluster_idx};
        
        for year_idx = 1:length(years)
            year = years{year_idx};
            year_str = sprintf('year_%s', year);
            
            if ~isfield(global_eventStruct.(rackName).(cluster_name), year_str)
                continue;
            end
            
            events = fieldnames(global_eventStruct.(rackName).(cluster_name).(year_str));
            
            for event_idx = 1:length(events)
                event_name = events{event_idx};
                evt = global_eventStruct.(rackName).(cluster_name).(year_str).(event_name);
                
                if ~isnan(evt.avg_current)
                    all_currents = [all_currents, evt.avg_current];
                    cluster_labels{end+1} = cluster_name;
                end
            end
        end
    end
    
    if ~isempty(all_currents)
        % Create histogram
        histogram(all_currents, 20, 'FaceAlpha', 0.7);
        title(sprintf('%s Current Distribution', rackName), 'FontSize', 12);
        xlabel('Average Current (A)', 'FontSize', 10);
        ylabel('Frequency', 'FontSize', 10);
        grid on;
    end
end

sgtitle('Current Distribution by Rack', 'FontSize', 14);

%% Figure 1-2. Power Distribution by Rack
figure(2)
set(gcf, 'Position', [100, 100, 800, 600]);

for rack_idx = 1:length(rackNames)
    rackName = rackNames{rack_idx};
    
    subplot(2, 4, rack_idx);
    
    % Collect all power data for this rack
    all_powers = [];
    cluster_labels = {};
    
    clusters = fieldnames(global_eventStruct.(rackName));
    
    for cluster_idx = 1:length(clusters)
        cluster_name = clusters{cluster_idx};
        
        for year_idx = 1:length(years)
            year = years{year_idx};
            year_str = sprintf('year_%s', year);
            
            if ~isfield(global_eventStruct.(rackName).(cluster_name), year_str)
                continue;
            end
            
            events = fieldnames(global_eventStruct.(rackName).(cluster_name).(year_str));
            
            for event_idx = 1:length(events)
                event_name = events{event_idx};
                evt = global_eventStruct.(rackName).(cluster_name).(year_str).(event_name);
                
                % Calculate average power from P_seq data
                    avg_power = mean(evt.P_seq);
                    all_powers = [all_powers, avg_power];
                    cluster_labels{end+1} = cluster_name;
            end
        end
    end
    
    if ~isempty(all_powers)
        % Create histogram
        histogram(all_powers, 20, 'FaceAlpha', 0.7);
        title(sprintf('%s Power Distribution', rackName), 'FontSize', 12);
        xlabel('Average Power (kW)', 'FontSize', 10);
        ylabel('Frequency', 'FontSize', 10);
        grid on;
    end
end

sgtitle('Power Distribution by Rack', 'FontSize', 14);

%% Figure 3: Time-Current and Time-Voltage Profiles by Rack
for rack_idx = 1:length(rackNames)
    rackName = rackNames{rack_idx};
    
    % Find clusters for this rack that have data in ALL years
    clusters = fieldnames(global_eventStruct.(rackName));
    valid_clusters = {};
    cluster_scores = [];
    
    for cluster_idx = 1:length(clusters)
        cluster_name = clusters{cluster_idx};
        
        % Check if this cluster has data in ALL years
        has_all_years = true;
        total_events = 0;
        
        for year_idx = 1:length(years)
            year = years{year_idx};
            year_str = sprintf('year_%s', year);
            
            if ~isfield(global_eventStruct.(rackName).(cluster_name), year_str)
                has_all_years = false;
                break;
            end
            
            events = fieldnames(global_eventStruct.(rackName).(cluster_name).(year_str));
            if isempty(events)
                has_all_years = false;
                break;
            end
            total_events = total_events + length(events);
        end
        
        % Only include clusters that have data in ALL years
        if has_all_years && total_events > 0
            valid_clusters{end+1} = cluster_name;
            cluster_scores(end+1) = total_events;
        end
    end
    
    if isempty(valid_clusters)
        fprintf('Rack %s: No clusters with data in all years\n', rackName);
        continue;
    end
    
    % Get the cluster with most events for this rack
    [~, max_idx] = max(cluster_scores);
    best_cluster = valid_clusters{max_idx};
    
    fprintf('Rack %s: Best cluster %s (%d events)\n', rackName, best_cluster, cluster_scores(max_idx));
    
    % Create figure for this rack
    figure(3 + rack_idx);
    set(gcf, 'Position', [100, 100, 800, 600]);
    
    % First row: Current plots by year (subplot 1, 2, 3)
    for year_idx = 1:length(years)
        year = years{year_idx};
        year_str = sprintf('year_%s', year);
        
        subplot(2, 3, year_idx);
        hold on;
        
        if ~isfield(global_eventStruct.(rackName).(best_cluster), year_str)
            title(sprintf('%s %s Current (%s) - No Data', rackName, best_cluster, year), 'FontSize', 10);
            continue;
        end
        
        events = fieldnames(global_eventStruct.(rackName).(best_cluster).(year_str));
        
        if isempty(events)
            title(sprintf('%s %s Current (%s) - No Data', rackName, best_cluster, year), 'FontSize', 10);
            continue;
        end
        
        plot_count = 0;
        for event_idx = 1:length(events)
            event_name = events{event_idx};
            evt = global_eventStruct.(rackName).(best_cluster).(year_str).(event_name);
            
            % Convert time data to seconds
            t_data = evt.t_seq;
            if isstring(t_data) || iscell(t_data)
                % Convert string datetime to datetime object
                t_datetime = datetime(t_data);
                t_seconds = seconds(t_datetime - t_datetime(1));
            else
                % Handle non-datetime time data
                t_seconds = double(t_data) - double(t_data(1));
            end
            
            % Plot current
            plot(t_seconds, evt.I_seq, 'LineWidth', 1);
            plot_count = plot_count + 1;
        end
        
        title(sprintf('%s %s Current (%s) - %d events', rackName, best_cluster, year, length(events)), 'FontSize', 10);
        xlabel('Time (s)', 'FontSize', 8);
        ylabel('Current (A)', 'FontSize', 8);
        grid on;
        
        if plot_count == 0
            text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        end
    end
    
    % Second row: Voltage plots by year (subplot 4, 5, 6)
    for year_idx = 1:length(years)
        year = years{year_idx};
        year_str = sprintf('year_%s', year);
        
        subplot(2, 3, year_idx + 3);
        hold on;
        
        if ~isfield(global_eventStruct.(rackName).(best_cluster), year_str)
            title(sprintf('%s %s Voltage (%s) - No Data', rackName, best_cluster, year), 'FontSize', 10);
            continue;
        end
        
        events = fieldnames(global_eventStruct.(rackName).(best_cluster).(year_str));
        
        if isempty(events)
            title(sprintf('%s %s Voltage (%s) - No Data', rackName, best_cluster, year), 'FontSize', 10);
            continue;
        end
        
        plot_count = 0;
        for event_idx = 1:length(events)
            event_name = events{event_idx};
            evt = global_eventStruct.(rackName).(best_cluster).(year_str).(event_name);
            
            % Convert time data to seconds
            t_data = evt.t_seq;
            if isstring(t_data) || iscell(t_data)
                % Convert string datetime to datetime object
                t_datetime = datetime(t_data);
                t_seconds = seconds(t_datetime - t_datetime(1));
            else
                % Handle non-datetime time data
                t_seconds = double(t_data) - double(t_data(1));
            end
            
            % Plot voltage
            plot(t_seconds, evt.V_seq, 'LineWidth', 1);
            plot_count = plot_count + 1;
        end
        
        title(sprintf('%s %s Voltage (%s) - %d events', rackName, best_cluster, year, length(events)), 'FontSize', 10);
        xlabel('Time (s)', 'FontSize', 8);
        ylabel('Voltage (V)', 'FontSize', 8);
        grid on;
        
        if plot_count == 0
            text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        end
    end
    
    sgtitle(sprintf('Time-Current and Time-Voltage Profiles for %s', rackName), 'FontSize', 14);
end

%% Figure 4: DCIR Boxplots by Year for each Rack
for rack_idx = 1:length(rackNames)
    rackName = rackNames{rack_idx};
    
    % Find clusters for this rack that have data in ALL years
    clusters = fieldnames(global_eventStruct.(rackName));
    valid_clusters = {};
    cluster_scores = [];
    
    for cluster_idx = 1:length(clusters)
        cluster_name = clusters{cluster_idx};
        
        % Check if this cluster has data in ALL years
        has_all_years = true;
        total_events = 0;
        
        for year_idx = 1:length(years)
            year = years{year_idx};
            year_str = sprintf('year_%s', year);
            
            if ~isfield(global_eventStruct.(rackName).(cluster_name), year_str)
                has_all_years = false;
                break;
            end
            
            events = fieldnames(global_eventStruct.(rackName).(cluster_name).(year_str));
            if isempty(events)
                has_all_years = false;
                break;
            end
            total_events = total_events + length(events);
        end
        
        % Only include clusters that have data in ALL years
        if has_all_years && total_events > 0
            valid_clusters{end+1} = cluster_name;
            cluster_scores(end+1) = total_events;
        end
    end
    
    if isempty(valid_clusters)
        fprintf('Rack %s: No clusters with data in all years\n', rackName);
        continue;
    end
    
    % Get the cluster with most events for this rack
    [~, max_idx] = max(cluster_scores);
    best_cluster = valid_clusters{max_idx};
    
    % Create figure for DCIR boxplots
    figure(12 + rack_idx);
    set(gcf, 'Position', [100, 100, 1200, 800]);
    
    % DCIR time points
    dcir_times = [1, 3, 5, 10, 30];
    
    % Create subplots for each year
    for year_idx = 1:length(years)
        year = years{year_idx};
        year_str = sprintf('year_%s', year);
        
        subplot(1, 3, year_idx);
        
        if ~isfield(global_eventStruct.(rackName).(best_cluster), year_str)
            title(sprintf('%s %s DCIR (%s) - No Data', rackName, best_cluster, year), 'FontSize', 12);
            continue;
        end
        
        events = fieldnames(global_eventStruct.(rackName).(best_cluster).(year_str));
        
        if isempty(events)
            title(sprintf('%s %s DCIR (%s) - No Data', rackName, best_cluster, year), 'FontSize', 12);
            continue;
        end
        
        % Collect DCIR data for each time point
        dcir_data = cell(length(dcir_times), 1);
        
        for event_idx = 1:length(events)
            event_name = events{event_idx};
            evt = global_eventStruct.(rackName).(best_cluster).(year_str).(event_name);
            
            % Extract DCIR values for each time point
            for time_idx = 1:length(dcir_times)
                time_str = sprintf('DCIR_%ds', dcir_times(time_idx));
                if isfield(evt, time_str) && ~isnan(evt.(time_str).val)
                    dcir_data{time_idx} = [dcir_data{time_idx}, evt.(time_str).val];
                end
            end
        end
        
        % Create boxplot with original data (no outlier removal)
        if ~all(cellfun(@isempty, dcir_data))
            % Prepare data for boxplot
            all_dcir_values = [];
            all_time_labels = [];
            
            for time_idx = 1:length(dcir_times)
                if ~isempty(dcir_data{time_idx})
                    all_dcir_values = [all_dcir_values; dcir_data{time_idx}'];
                    all_time_labels = [all_time_labels; repmat(time_idx, length(dcir_data{time_idx}), 1)];
                end
            end
            
            if ~isempty(all_dcir_values)
                % Create labels only for time points that have data
                available_labels = {};
                time_indices_with_data = [];
                for time_idx = 1:length(dcir_times)
                    if ~isempty(dcir_data{time_idx})
                        available_labels{end+1} = sprintf('%ds', dcir_times(time_idx));
                        time_indices_with_data = [time_indices_with_data, time_idx];
                    end
                end
                
                boxplot(all_dcir_values, all_time_labels, 'Labels', available_labels);
                
                % Add mean values as red dots
                hold on;
                for i = 1:length(time_indices_with_data)
                    time_idx = time_indices_with_data(i);
                    if ~isempty(dcir_data{time_idx})
                        mean_val = mean(dcir_data{time_idx});
                        plot(i, mean_val, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                    end
                end
                hold off;
                
                title(sprintf('%s %s DCIR (%s) - %d events', rackName, best_cluster, year, length(events)), 'FontSize', 12);
                xlabel('Time', 'FontSize', 10);
                ylabel('DCIR (mΩ)', 'FontSize', 10);
                ylim([0 1.8]);
                grid on;
            else
                title(sprintf('%s %s DCIR (%s) - No Valid Data', rackName, best_cluster, year), 'FontSize', 12);
            end
        % else
        %     title(sprintf('%s %s DCIR (%s) - No Valid Data', rackName, best_cluster, year), 'FontSize', 12);
        % end
        end
    end
    
    % sgtitle(sprintf('DCIR Boxplots for %s', rackName), 'FontSize', 14);
end

%% Save figures
saveDir = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_Rack_DCIR\DCIR_Charge_0728\Figures';

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

saveas(figure(1), fullfile(saveDir, 'DCIR_Trends_by_Rack.fig'));

saveas(figure(2), fullfile(saveDir, 'Power_Distribution_by_Rack.fig'));

% Save DCIR boxplot figures for each rack
for rack_idx = 1:length(rackNames)
    rackName = rackNames{rack_idx};
    saveas(figure(12 + rack_idx), fullfile(saveDir, sprintf('DCIR_Boxplots_%s.fig', rackName)));
end

fprintf('Visualization completed. Figures saved to %s\n', saveDir);