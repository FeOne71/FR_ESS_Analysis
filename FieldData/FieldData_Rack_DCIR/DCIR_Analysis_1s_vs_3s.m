%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCIR Analysis: 1s vs 3s Comparison
% Analyze why 1s DCIR values are larger than 3s DCIR values
% Focus on Rack01 data analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Load data
dataPath = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_Rack_DCIR\DCIR_Charge_0728\all_chg_events_current_clustering_all_years.mat';
load(dataPath);

%% Analysis parameters
rackName = 'Rack01';
min_events_per_cluster = 1;  % Include all events (minimum 1 event)
dcir_1s_field = 'DCIR_1s';
dcir_2s_field = 'DCIR_3s';  % 1s vs 3s comparison

%% Extract Rack01 data
if ~isfield(global_eventStruct, rackName)
    error('Rack01 data not found in the loaded structure');
end

rack01_data = global_eventStruct.(rackName);
cluster_names = fieldnames(rack01_data);

fprintf('=== DCIR 1s vs 3s Analysis for %s ===\n', rackName);
fprintf('Total clusters found: %d\n', length(cluster_names));

%% Initialize analysis results
analysis_results = struct();
year_analysis = struct();  % Year-specific analysis
total_events = 0;
anomalous_events = 0;

%% Analyze each cluster
for c = 1:length(cluster_names)
    cluster_name = cluster_names{c};
    cluster_data = rack01_data.(cluster_name);
    year_names = fieldnames(cluster_data);
    
    fprintf('\n--- Cluster: %s ---\n', cluster_name);
    
    cluster_events = [];
    cluster_analysis = struct();
    
    % Collect all events from all years for this cluster
    for y = 1:length(year_names)
        year_name = year_names{y};
        year_data = cluster_data.(year_name);
        event_names = fieldnames(year_data);
        
        % Initialize year-specific data
        if ~isfield(year_analysis, year_name)
            year_analysis.(year_name) = struct();
        end
        if ~isfield(year_analysis.(year_name), cluster_name)
            year_analysis.(year_name).(cluster_name) = [];
        end
        
        for e = 1:length(event_names)
            event_name = event_names{e};
            event_data = year_data.(event_name);
            
            % Extract DCIR values
            dcir_1s = event_data.(dcir_1s_field).val;
            dcir_3s = event_data.(dcir_2s_field).val;
            
            % Check if both values are valid
            if ~isnan(dcir_1s) && ~isnan(dcir_3s)
                total_events = total_events + 1;
                
                % Calculate difference
                dcir_diff = dcir_1s - dcir_3s;
                
                % Store event info
                event_info = struct();
                event_info.event_name = event_name;
                event_info.year = year_name;
                event_info.dcir_1s = dcir_1s;
                event_info.dcir_3s = dcir_3s;
                event_info.dcir_diff = dcir_diff;
                event_info.is_anomalous = dcir_diff > 0;
                
                % Extract voltage and current data for detailed analysis
                event_info.V1_1s = event_data.(dcir_1s_field).V1;
                event_info.V2_1s = event_data.(dcir_1s_field).V2;
                event_info.I1_1s = event_data.(dcir_1s_field).I1;
                event_info.I2_1s = event_data.(dcir_1s_field).I2;
                event_info.dV_1s = event_data.(dcir_1s_field).dV;
                event_info.dI_1s = event_data.(dcir_1s_field).dI;
                
                event_info.V1_3s = event_data.(dcir_2s_field).V1;
                event_info.V2_3s = event_data.(dcir_2s_field).V2;
                event_info.I1_3s = event_data.(dcir_2s_field).I1;
                event_info.I2_3s = event_data.(dcir_2s_field).I2;
                event_info.dV_3s = event_data.(dcir_2s_field).dV;
                event_info.dI_3s = event_data.(dcir_2s_field).dI;
                
                % Extract time series data for trend analysis
                event_info.t_seq = event_data.t_seq;
                event_info.I_seq = event_data.I_seq;
                event_info.V_seq = event_data.V_seq;
                event_info.soc_seq = event_data.soc_seq;
                event_info.T_seq = event_data.T_seq;
                
                cluster_events = [cluster_events, event_info];
                year_analysis.(year_name).(cluster_name) = [year_analysis.(year_name).(cluster_name), event_info];
                
                if event_info.is_anomalous
                    anomalous_events = anomalous_events + 1;
                end
            end
        end
    end
    
    % Analyze this cluster (include all events)
    if length(cluster_events) >= min_events_per_cluster
        % Calculate statistics
        dcir_1s_vals = [cluster_events.dcir_1s];
        dcir_3s_vals = [cluster_events.dcir_3s];
        dcir_diffs = [cluster_events.dcir_diff];
        
        cluster_analysis.total_events = length(cluster_events);
        cluster_analysis.anomalous_events = sum([cluster_events.is_anomalous]);
        cluster_analysis.anomalous_ratio = cluster_analysis.anomalous_events / cluster_analysis.total_events;
        cluster_analysis.mean_dcir_1s = mean(dcir_1s_vals);
        cluster_analysis.mean_dcir_3s = mean(dcir_3s_vals);
        cluster_analysis.mean_diff = mean(dcir_diffs);
        cluster_analysis.std_diff = std(dcir_diffs);
        
        fprintf('Total events: %d\n', cluster_analysis.total_events);
        fprintf('Anomalous events: %d (%.1f%%)\n', cluster_analysis.anomalous_events, cluster_analysis.anomalous_ratio * 100);
        fprintf('Mean DCIR 1s: %.2f mΩ\n', cluster_analysis.mean_dcir_1s);
        fprintf('Mean DCIR 3s: %.2f mΩ\n', cluster_analysis.mean_dcir_3s);
        fprintf('Mean difference: %.2f mΩ\n', cluster_analysis.mean_diff);
        
        % Store cluster analysis
        analysis_results.(cluster_name) = cluster_analysis;
        analysis_results.(cluster_name).events = cluster_events;
    else
        fprintf('No events for analysis (%d events)\n', length(cluster_events));
    end
end

%% Overall statistics
fprintf('\n=== Overall Statistics ===\n');
fprintf('Total valid events: %d\n', total_events);
fprintf('Total anomalous events: %d (%.1f%%)\n', anomalous_events, anomalous_events/total_events * 100);

%% Detailed analysis of anomalous events
fprintf('\n=== Detailed Analysis of Anomalous Events ===\n');

% Collect all anomalous events
all_anomalous_events = [];
cluster_names_analysis = fieldnames(analysis_results);

for c = 1:length(cluster_names_analysis)
    cluster_name = cluster_names_analysis{c};
    if isfield(analysis_results.(cluster_name), 'events')
        cluster_events = analysis_results.(cluster_name).events;
        anomalous_in_cluster = cluster_events([cluster_events.is_anomalous]);
        all_anomalous_events = [all_anomalous_events, anomalous_in_cluster];
    end
end

fprintf('Found %d anomalous events for detailed analysis\n', length(all_anomalous_events));

%% Analyze patterns in anomalous events
if length(all_anomalous_events) > 0
    % Pattern 1: Check if anomalous events have specific voltage/current characteristics
    fprintf('\n--- Pattern Analysis ---\n');
    
    % Analyze voltage and current changes
    dV_1s_anomalous = [all_anomalous_events.dV_1s];
    dV_3s_anomalous = [all_anomalous_events.dV_3s];
    dI_1s_anomalous = [all_anomalous_events.dI_1s];
    dI_3s_anomalous = [all_anomalous_events.dI_3s];
    
    fprintf('Anomalous events - Mean dV (1s): %.4f V\n', mean(dV_1s_anomalous));
    fprintf('Anomalous events - Mean dV (3s): %.4f V\n', mean(dV_3s_anomalous));
    fprintf('Anomalous events - Mean dI (1s): %.4f A\n', mean(dI_1s_anomalous));
    fprintf('Anomalous events - Mean dI (3s): %.4f A\n', mean(dI_3s_anomalous));
    
    % Check if dV/dI ratios are consistent
    dcir_calc_1s = abs(dV_1s_anomalous ./ dI_1s_anomalous) * 1000;
    dcir_calc_3s = abs(dV_3s_anomalous ./ dI_3s_anomalous) * 1000;
    
    fprintf('Calculated DCIR 1s (from dV/dI): %.2f mΩ\n', mean(dcir_calc_1s));
    fprintf('Calculated DCIR 3s (from dV/dI): %.2f mΩ\n', mean(dcir_calc_3s));
    
    % Pattern 2: Check timing of measurements
    fprintf('\n--- Timing Analysis ---\n');
    
    % Analyze if the issue is related to measurement timing
    for i = 1:min(5, length(all_anomalous_events))  % Show first 5 examples
        event = all_anomalous_events(i);
        fprintf('Event %s (%s):\n', event.event_name, event.year);
        fprintf('  1s: V1=%.3fV, V2=%.3fV, I1=%.3fA, I2=%.3fA, dV=%.4fV, dI=%.4fA\n', ...
            event.V1_1s, event.V2_1s, event.I1_1s, event.I2_1s, event.dV_1s, event.dI_1s);
        fprintf('  3s: V1=%.3fV, V2=%.3fV, I1=%.3fA, I2=%.3fA, dV=%.4fV, dI=%.4fA\n', ...
            event.V1_3s, event.V2_3s, event.I1_3s, event.I2_3s, event.dV_3s, event.dI_3s);
        fprintf('  DCIR 1s: %.2f mΩ, DCIR 3s: %.2f mΩ, Diff: %.2f mΩ\n', ...
            event.dcir_1s, event.dcir_3s, event.dcir_diff);
        fprintf('\n');
    end
end

%% Year-specific Analysis and Visualization
fprintf('\n=== Year-specific Analysis ===\n');

% Get all years
all_years = fieldnames(year_analysis);
fprintf('Years available: %s\n', strjoin(all_years, ', '));

% Create separate figures for each year
for year_idx = 1:length(all_years)
    year_name = all_years{year_idx};
    fprintf('\n--- Analyzing %s ---\n', year_name);
    
    % Get year-specific data
    year_data = year_analysis.(year_name);
    year_clusters = fieldnames(year_data);
    
    % Debug: 실제 클러스터 이름 확인
    fprintf('  Available clusters in %s: %s\n', year_name, strjoin(year_clusters, ', '));
    
    % Calculate year-specific statistics
    year_total_events = 0;
    year_anomalous_events = 0;
    year_analysis_results = struct();
    
    for c = 1:length(year_clusters)
        cluster_name = year_clusters{c};
        cluster_events = year_data.(cluster_name);
        
        if length(cluster_events) >= min_events_per_cluster
            dcir_1s_vals = [cluster_events.dcir_1s];
            dcir_3s_vals = [cluster_events.dcir_3s];
            dcir_diffs = [cluster_events.dcir_diff];
            
            year_total_events = year_total_events + length(cluster_events);
            year_anomalous_events = year_anomalous_events + sum([cluster_events.is_anomalous]);
            
            year_analysis_results.(cluster_name).total_events = length(cluster_events);
            year_analysis_results.(cluster_name).anomalous_events = sum([cluster_events.is_anomalous]);
            year_analysis_results.(cluster_name).anomalous_ratio = sum([cluster_events.is_anomalous]) / length(cluster_events);
            year_analysis_results.(cluster_name).mean_dcir_1s = mean(dcir_1s_vals);
            year_analysis_results.(cluster_name).mean_dcir_3s = mean(dcir_3s_vals);
            year_analysis_results.(cluster_name).mean_diff = mean(dcir_diffs);
            year_analysis_results.(cluster_name).events = cluster_events;
            
            fprintf('  %s: %d events, %d anomalous (%.1f%%)\n', ...
                cluster_name, length(cluster_events), sum([cluster_events.is_anomalous]), ...
                sum([cluster_events.is_anomalous])/length(cluster_events)*100);
        else
            fprintf('  %s: No events\n', cluster_name);
        end
    end
    
    fprintf('  %s Total: %d events, %d anomalous (%.1f%%)\n', ...
        year_name, year_total_events, year_anomalous_events, ...
        year_anomalous_events/year_total_events*100);
    
    % Create visualization for this year
    if year_total_events > 0
        figure('Name', sprintf('DCIR 1s vs 3s Analysis - %s', year_name), 'Position', [100, 100, 1200, 800]);
        
        % Subplot 1: Scatter plot of DCIR 1s vs 3s
        subplot(2, 3, 1);
        all_dcir_1s = [];
        all_dcir_3s = [];
        all_colors = [];
        
        year_cluster_names = fieldnames(year_analysis_results);
        for c = 1:length(year_cluster_names)
            cluster_name = year_cluster_names{c};
            if isfield(year_analysis_results.(cluster_name), 'events')
                cluster_events = year_analysis_results.(cluster_name).events;
                dcir_1s = [cluster_events.dcir_1s];
                dcir_3s = [cluster_events.dcir_3s];
                
                if ~isempty(dcir_1s)
                    all_dcir_1s = [all_dcir_1s, dcir_1s];
                    all_dcir_3s = [all_dcir_3s, dcir_3s];
                    all_colors = [all_colors, repmat(c, 1, length(dcir_1s))];
                end
            end
        end
        
        if ~isempty(all_dcir_1s)
            scatter(all_dcir_1s, all_dcir_3s, 50, all_colors, 'filled', 'MarkerFaceAlpha', 0.6);
            hold on;
            plot([0, max(all_dcir_1s)], [0, max(all_dcir_1s)], 'r--', 'LineWidth', 2);
            xlabel('DCIR 1s (mΩ)');
            ylabel('DCIR 3s (mΩ)');
            title(sprintf('DCIR 1s vs 3s Comparison - %s', year_name));
            grid on;
            colorbar;
        else
            text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
            title(sprintf('DCIR 1s vs 3s Comparison - %s', year_name));
        end
        
        % Subplot 2: Histogram of differences
        subplot(2, 3, 2);
        all_diffs = [];
        for c = 1:length(year_cluster_names)
            cluster_name = year_cluster_names{c};
            if isfield(year_analysis_results.(cluster_name), 'events')
                cluster_events = year_analysis_results.(cluster_name).events;
                diffs = [cluster_events.dcir_diff];
                if ~isempty(diffs)
                    all_diffs = [all_diffs, diffs];
                end
            end
        end
        
        if ~isempty(all_diffs)
            histogram(all_diffs, 30, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
            hold on;
            xline(0, 'r--', 'LineWidth', 2);
            xlabel('DCIR 1s - DCIR 3s (mΩ)');
            ylabel('Frequency');
            title(sprintf('Distribution of DCIR Differences - %s', year_name));
            grid on;
        else
            text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
            title(sprintf('Distribution of DCIR Differences - %s', year_name));
        end
        
        % Subplot 3: Box plot by cluster
        subplot(2, 3, 3);
        cluster_diffs = [];
        cluster_labels = {};
        group_labels = [];
        
        for c = 1:length(year_cluster_names)
            cluster_name = year_cluster_names{c};
            if isfield(year_analysis_results.(cluster_name), 'events')
                cluster_events = year_analysis_results.(cluster_name).events;
                diffs = [cluster_events.dcir_diff];
                if ~isempty(diffs)
                    cluster_diffs = [cluster_diffs, diffs];
                    cluster_labels{end+1} = cluster_name;
                    group_labels = [group_labels, repmat(c, 1, length(diffs))];
                end
            end
        end
        
        if ~isempty(cluster_diffs)
            boxplot(cluster_diffs, group_labels, 'Labels', cluster_labels);
            ylabel('DCIR 1s - DCIR 3s (mΩ)');
            title(sprintf('DCIR Differences by Cluster - %s', year_name));
            grid on;
            xtickangle(45);
        else
            text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
            title(sprintf('DCIR Differences by Cluster - %s', year_name));
        end
        
        % Subplot 4: Anomalous ratio by cluster
        subplot(2, 3, 4);
        anomalous_ratios = [];
        cluster_names_plot = {};
        
        for c = 1:length(year_cluster_names)
            cluster_name = year_cluster_names{c};
            if isfield(year_analysis_results.(cluster_name), 'anomalous_ratio')
                anomalous_ratios = [anomalous_ratios, year_analysis_results.(cluster_name).anomalous_ratio];
                cluster_names_plot{end+1} = cluster_name;
            end
        end
        
        if ~isempty(anomalous_ratios)
            bar(anomalous_ratios);
            set(gca, 'XTickLabel', cluster_names_plot);
            ylabel('Anomalous Ratio');
            title(sprintf('Anomalous Events Ratio by Cluster - %s', year_name));
            grid on;
            xtickangle(45);
        else
            text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
            title(sprintf('Anomalous Events Ratio by Cluster - %s', year_name));
        end
        
        % Subplot 5: Mean DCIR values by cluster (Clustering_Auto와 동일한 클러스터 사용)
        subplot(2, 3, 5);
        
        % Clustering_Auto와 동일하게 가장 많은 이벤트를 가진 클러스터 찾기
        cluster_names_available = fieldnames(year_analysis_results);
        best_cluster_name = '';
        max_events = 0;
        
        for c = 1:length(cluster_names_available)
            cluster_name = cluster_names_available{c};
            if isfield(year_analysis_results.(cluster_name), 'total_events')
                if year_analysis_results.(cluster_name).total_events > max_events
                    max_events = year_analysis_results.(cluster_name).total_events;
                    best_cluster_name = cluster_name;
                end
            end
        end
        
        if ~isempty(best_cluster_name) && isfield(year_analysis_results, best_cluster_name)
            % Clustering_Auto와 동일한 방식으로 원본 데이터에서 DCIR 값 추출
            target_cluster_events = year_analysis_results.(best_cluster_name).events;
            
            % Clustering_Auto와 동일한 조건: ~isnan(evt.(time_str).val)만 체크
            dcir_1s_target = [];
            dcir_3s_target = [];
            
            for evt_idx = 1:length(target_cluster_events)
                evt = target_cluster_events(evt_idx);
                % Clustering_Auto와 동일한 조건 사용
                if ~isnan(evt.dcir_1s)  % 1s만 체크 (Clustering_Auto와 동일)
                    dcir_1s_target = [dcir_1s_target, evt.dcir_1s];
                end
                if ~isnan(evt.dcir_3s)  % 3s만 체크 (Clustering_Auto와 동일)
                    dcir_3s_target = [dcir_3s_target, evt.dcir_3s];
                end
            end
            
            % Clustering_Auto와 동일한 평균 계산
            target_mean_1s = mean(dcir_1s_target);
            target_mean_3s = mean(dcir_3s_target);
            
            % Bar chart로 표시
            x_pos = [1, 2];
            bar(x_pos, [target_mean_1s, target_mean_3s]);
            set(gca, 'XTickLabel', {'1s', '3s'});
            ylabel('Mean DCIR (mΩ)');
            title(sprintf('Target Cluster Mean DCIR Values - %s (%s)', year_name, best_cluster_name));
            grid on;
            
            % 평균값 텍스트로 표시
            text(1, target_mean_1s + 0.05, sprintf('%.2f', target_mean_1s), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            text(2, target_mean_3s + 0.05, sprintf('%.2f', target_mean_3s), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        else
            text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
            title(sprintf('Target Cluster Mean DCIR Values - %s', year_name));
        end
        
        % Subplot 6: dV vs dI analysis for anomalous events
        subplot(2, 3, 6);
        all_anomalous_events_year = [];
        for c = 1:length(year_cluster_names)
            cluster_name = year_cluster_names{c};
            if isfield(year_analysis_results.(cluster_name), 'events')
                cluster_events = year_analysis_results.(cluster_name).events;
                anomalous_in_cluster = cluster_events([cluster_events.is_anomalous]);
                all_anomalous_events_year = [all_anomalous_events_year, anomalous_in_cluster];
            end
        end
        
        if length(all_anomalous_events_year) > 0
            dV_1s_anom = [all_anomalous_events_year.dV_1s];
            dI_1s_anom = [all_anomalous_events_year.dI_1s];
            dV_3s_anom = [all_anomalous_events_year.dV_3s];
            dI_3s_anom = [all_anomalous_events_year.dI_3s];
            
            scatter(dV_1s_anom, dI_1s_anom, 50, 'b', 'filled', 'DisplayName', '1s');
            hold on;
            scatter(dV_3s_anom, dI_3s_anom, 50, 'r', 'filled', 'DisplayName', '3s');
            xlabel('ΔV (V)');
            ylabel('ΔI (A)');
            title(sprintf('ΔV vs ΔI for Anomalous Events - %s', year_name));
            legend;
            grid on;
        else
            text(0.5, 0.5, 'No anomalous events', 'HorizontalAlignment', 'center', 'FontSize', 12);
            title(sprintf('ΔV vs ΔI for Anomalous Events - %s', year_name));
        end
        
        sgtitle(sprintf('DCIR Analysis: 1s vs 3s Comparison - %s (%s)', rackName, year_name));
    end
end

%% Overall Visualization (Original)
fprintf('\n=== Creating Overall Visualizations ===\n');

% Create figure for DCIR comparison
figure('Name', 'DCIR 1s vs 3s Analysis - Overall', 'Position', [100, 100, 1200, 800]);

% Subplot 1: Scatter plot of DCIR 1s vs 3s
subplot(2, 3, 1);
all_dcir_1s = [];
all_dcir_3s = [];
all_colors = [];

for c = 1:length(cluster_names_analysis)
    cluster_name = cluster_names_analysis{c};
    if isfield(analysis_results.(cluster_name), 'events')
        cluster_events = analysis_results.(cluster_name).events;
        dcir_1s = [cluster_events.dcir_1s];
        dcir_3s = [cluster_events.dcir_3s];
        
        if ~isempty(dcir_1s)
            all_dcir_1s = [all_dcir_1s, dcir_1s];
            all_dcir_3s = [all_dcir_3s, dcir_3s];
            all_colors = [all_colors, repmat(c, 1, length(dcir_1s))];
        end
    end
end

if ~isempty(all_dcir_1s)
    scatter(all_dcir_1s, all_dcir_3s, 50, all_colors, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    plot([0, max(all_dcir_1s)], [0, max(all_dcir_1s)], 'r--', 'LineWidth', 2);
    xlabel('DCIR 1s (mΩ)');
    ylabel('DCIR 3s (mΩ)');
    title('DCIR 1s vs 3s Comparison');
    grid on;
    colorbar;
else
    text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('DCIR 1s vs 3s Comparison');
end

% Subplot 2: Histogram of differences
subplot(2, 3, 2);
all_diffs = [];
for c = 1:length(cluster_names_analysis)
    cluster_name = cluster_names_analysis{c};
    if isfield(analysis_results.(cluster_name), 'events')
        cluster_events = analysis_results.(cluster_name).events;
        diffs = [cluster_events.dcir_diff];
        if ~isempty(diffs)
            all_diffs = [all_diffs, diffs];
        end
    end
end

if ~isempty(all_diffs)
    histogram(all_diffs, 30, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
    hold on;
    xline(0, 'r--', 'LineWidth', 2);
    xlabel('DCIR 1s - DCIR 3s (mΩ)');
    ylabel('Frequency');
    title('Distribution of DCIR Differences');
    grid on;
else
    text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('Distribution of DCIR Differences');
end

% Subplot 3: Box plot by cluster
subplot(2, 3, 3);
cluster_diffs = [];
cluster_labels = {};
group_labels = [];

for c = 1:length(cluster_names_analysis)
    cluster_name = cluster_names_analysis{c};
    if isfield(analysis_results.(cluster_name), 'events')
        cluster_events = analysis_results.(cluster_name).events;
        diffs = [cluster_events.dcir_diff];
        if ~isempty(diffs)
            cluster_diffs = [cluster_diffs, diffs];
            cluster_labels{end+1} = cluster_name;
            group_labels = [group_labels, repmat(c, 1, length(diffs))];
        end
    end
end

if ~isempty(cluster_diffs)
    boxplot(cluster_diffs, group_labels, 'Labels', cluster_labels);
    ylabel('DCIR 1s - DCIR 3s (mΩ)');
    title('DCIR Differences by Cluster');
    grid on;
    xtickangle(45);
else
    text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('DCIR Differences by Cluster');
end

% Subplot 4: Anomalous ratio by cluster
subplot(2, 3, 4);
anomalous_ratios = [];
cluster_names_plot = {};

for c = 1:length(cluster_names_analysis)
    cluster_name = cluster_names_analysis{c};
    if isfield(analysis_results.(cluster_name), 'anomalous_ratio')
        anomalous_ratios = [anomalous_ratios, analysis_results.(cluster_name).anomalous_ratio];
        cluster_names_plot{end+1} = cluster_name;
    end
end

if ~isempty(anomalous_ratios)
    bar(anomalous_ratios);
    set(gca, 'XTickLabel', cluster_names_plot);
    ylabel('Anomalous Ratio');
    title('Anomalous Events Ratio by Cluster');
    grid on;
    xtickangle(45);
else
    text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('Anomalous Events Ratio by Cluster');
end

% Subplot 5: Mean DCIR values by cluster (Clustering_Auto와 동일한 클러스터 사용)
subplot(2, 3, 5);

% Clustering_Auto와 동일하게 가장 많은 이벤트를 가진 클러스터 찾기
cluster_names_available = fieldnames(analysis_results);
best_cluster_name = '';
max_events = 0;

for c = 1:length(cluster_names_available)
    cluster_name = cluster_names_available{c};
    if isfield(analysis_results.(cluster_name), 'total_events')
        if analysis_results.(cluster_name).total_events > max_events
            max_events = analysis_results.(cluster_name).total_events;
            best_cluster_name = cluster_name;
        end
    end
end

if ~isempty(best_cluster_name) && isfield(analysis_results, best_cluster_name)
    % Clustering_Auto와 동일한 방식으로 원본 데이터에서 DCIR 값 추출
    target_cluster_events = analysis_results.(best_cluster_name).events;
    
    % Clustering_Auto와 동일한 조건: ~isnan(evt.(time_str).val)만 체크
    dcir_1s_target = [];
    dcir_3s_target = [];
    
    for evt_idx = 1:length(target_cluster_events)
        evt = target_cluster_events(evt_idx);
        % Clustering_Auto와 동일한 조건 사용
        if ~isnan(evt.dcir_1s)  % 1s만 체크 (Clustering_Auto와 동일)
            dcir_1s_target = [dcir_1s_target, evt.dcir_1s];
        end
        if ~isnan(evt.dcir_3s)  % 3s만 체크 (Clustering_Auto와 동일)
            dcir_3s_target = [dcir_3s_target, evt.dcir_3s];
        end
    end
    
    % Clustering_Auto와 동일한 평균 계산
    target_mean_1s = mean(dcir_1s_target);
    target_mean_3s = mean(dcir_3s_target);
    
    % Bar chart로 표시
    x_pos = [1, 2];
    bar(x_pos, [target_mean_1s, target_mean_3s]);
    set(gca, 'XTickLabel', {'1s', '3s'});
    ylabel('Mean DCIR (mΩ)');
    title(sprintf('Target Cluster Mean DCIR Values (%s)', best_cluster_name));
    grid on;
    
    % 평균값 텍스트로 표시
    text(1, target_mean_1s + 0.05, sprintf('%.2f', target_mean_1s), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(2, target_mean_3s + 0.05, sprintf('%.2f', target_mean_3s), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
else
    text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('Target Cluster Mean DCIR Values');
end

% Subplot 6: dV vs dI analysis for anomalous events
subplot(2, 3, 6);
if length(all_anomalous_events) > 0
    dV_1s_anom = [all_anomalous_events.dV_1s];
    dI_1s_anom = [all_anomalous_events.dI_1s];
    dV_3s_anom = [all_anomalous_events.dV_3s];
    dI_3s_anom = [all_anomalous_events.dI_3s];
    
    scatter(dV_1s_anom, dI_1s_anom, 50, 'b', 'filled', 'DisplayName', '1s');
    hold on;
    scatter(dV_3s_anom, dI_3s_anom, 50, 'r', 'filled', 'DisplayName', '3s');
    xlabel('ΔV (V)');
    ylabel('ΔI (A)');
    title('ΔV vs ΔI for Anomalous Events');
    legend;
    grid on;
end

sgtitle(sprintf('DCIR Analysis: 1s vs 3s Comparison - %s', rackName));

%% Save analysis results
save(fullfile('DCIR_Analysis_Results.mat'), 'analysis_results', 'all_anomalous_events');

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results saved to: DCIR_Analysis_Results.mat\n');
fprintf('Visualization completed.\n'); 