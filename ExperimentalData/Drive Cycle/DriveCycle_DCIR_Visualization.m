%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab Data - Drive Cycle DCIR Visualization
% Loads DriveCycle_Summary_Table.mat and creates integrated plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Drive Cycle DCIR Visualization ===\n');

%% Load Summary Table
if ~exist('DriveCycle_Summary_Table.mat', 'file')
    fprintf('ERROR: DriveCycle_Summary_Table.mat not found!\n');
    fprintf('Please run DriveCycle_DCIR_Analysis.m first to generate the summary table.\n');
    return;
end

load('DriveCycle_Summary_Table.mat', 'summaryTable');
fprintf('Loaded summary table: %d rows\n', height(summaryTable));

%% Create figures directory
if ~exist('figures', 'dir')
    mkdir('figures');
end
if ~exist('figures/Integrated_Capacity_Rchg', 'dir')
    mkdir('figures/Integrated_Capacity_Rchg');
        end

%% Define time intervals to plot
timeIntervals = {'R_1s', 'R_10s', 'R_30s'};
timeIntervalLabels = {'Rchg 1s', 'Rchg 10s', 'Rchg 30s'};

%% Get unique DC Profiles and Channels for color/marker mapping
uniqueDCProfiles = unique(summaryTable.DC_Profile);
uniqueChannels = unique(summaryTable.Channel);
    
% Sort DC profiles (DC1, DC2, ...)
dc_nums = zeros(length(uniqueDCProfiles), 1);
for d = 1:length(uniqueDCProfiles)
    dc_str = uniqueDCProfiles{d};
    num_match = regexp(dc_str, 'DC(\d+)', 'tokens');
    if ~isempty(num_match)
        dc_nums(d) = str2double(num_match{1}{1});
    else
        dc_nums(d) = 999;
    end
end
[~, dc_sort_idx] = sort(dc_nums);
sortedDCProfiles = uniqueDCProfiles(dc_sort_idx);

% Sort channels
sortedChannels = sort(uniqueChannels);

%% Define color map for DC Profiles (8 colors)
dcColors = lines(length(sortedDCProfiles));

%% Define marker styles for Channels (8 different markers)
channelMarkers = {'o', 's', '^', 'd', 'v', '>', '<', 'p'};
if length(sortedChannels) > length(channelMarkers)
    % If more channels than markers, cycle through markers
    channelMarkers = repmat(channelMarkers, 1, ceil(length(sortedChannels) / length(channelMarkers)));
end

%% Create 3 figures (one for each time interval)
for fig_idx = 1:length(timeIntervals)
    timeInterval = timeIntervals{fig_idx};
    timeLabel = timeIntervalLabels{fig_idx};
    
    % Check if column exists
    if ~ismember(timeInterval, summaryTable.Properties.VariableNames)
        fprintf('WARNING: Column %s not found in summary table. Skipping.\n', timeInterval);
        continue;
end

    % Create figure
    fig = figure('Name', sprintf('Capacity vs %s', timeLabel), ...
                 'Position', [100 + (fig_idx-1)*50, 100 + (fig_idx-1)*50, 1200, 800]);
    
    hold on;
    
    % Plot data for each DC Profile and Channel combination
    legendEntries = {};
    legendHandles = [];
    
    for dc_idx = 1:length(sortedDCProfiles)
        dcProfile = sortedDCProfiles{dc_idx};
        dcColor = dcColors(dc_idx, :);
        
        for ch_idx = 1:length(sortedChannels)
            channel = sortedChannels(ch_idx);
            marker = channelMarkers{ch_idx};
            
            % Filter data for this DC Profile and Channel
            mask = strcmp(summaryTable.DC_Profile, dcProfile) & ...
                   summaryTable.Channel == channel & ...
                   ~isnan(summaryTable.(timeInterval)) & ...
                   ~isnan(summaryTable.Capacity);
            
            if sum(mask) > 0
                capData = summaryTable.Capacity(mask);
                rchgData = summaryTable.(timeInterval)(mask);
        
                % Plot scatter points
                h = scatter(capData, rchgData, 100, dcColor, marker, ...
                           'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1.5);
                
                % Add to legend (only once per DC Profile)
                if ch_idx == 1
                    legendEntries{end+1} = dcProfile;
                    legendHandles(end+1) = h;
                end
            end
        end
    end
    
    % Add legend
    if ~isempty(legendHandles)
        legend(legendHandles, legendEntries, 'Location', 'best', 'FontSize', 10);
    end
    
    % Labels and title
    xlabel('Capacity (Ah)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(sprintf('%s (mΩ)', timeLabel), 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Capacity vs %s (All Channels & DC Profiles)', timeLabel), ...
          'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Add text annotation for channel markers (if space allows)
    if length(sortedChannels) <= 8
        markerLegend = cell(length(sortedChannels), 1);
        for ch_idx = 1:length(sortedChannels)
            markerLegend{ch_idx} = sprintf('Ch%d: %s', sortedChannels(ch_idx), channelMarkers{ch_idx});
        end
        text(0.02, 0.98, strjoin(markerLegend, '\n'), ...
             'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
    
    % Save figure
    saveas(fig, sprintf('figures/Integrated_Capacity_Rchg/Capacity_vs_%s.fig', timeInterval));
    fprintf('Saved: figures/Integrated_Capacity_Rchg/Capacity_vs_%s.fig\n', timeInterval);
end

fprintf('\n=== Visualization Complete ===\n');
fprintf('All figures saved to figures/Integrated_Capacity_Rchg/\n');
fprintf('Total figures created: %d\n', length(timeIntervals));
