%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 01_DriveCycle_EventAnalysis.m
% 이벤트 검출 및 Rchg 계산 (분석 및 검증)
% 
% 목적: 
% - Drive Cycle 데이터에서 Idle -> Load 전환 구간 검출
% - 필터링 조건 적용 (min_duration, max_P_std, max_I_std)
% - 각 이벤트에서 Rchg 계산 (1s, 3s, 5s, 10s, 30s, 60s)
% - debug_mode = true 시 필터링된 이벤트의 전류 그래프 저장
%
% 출력:
% - Lab_DC_DCIR_*cyc_Events.mat (각 사이클별 결과)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Configuration - User Settings
% =========================================================================
% 데이터 경로 설정
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
% 현재 스크립트가 있는 폴더의 Results 폴더에 저장
scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir, 'Results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

selectedChannels = [];  % Empty array: analyze all channels (Ch9~Ch16)
selectedCycles = {'0cyc','200cyc','400cyc','600cyc','800cyc'};  
selectedSOC = {'SOC90','SOC70','SOC50'};  

% 디버그 모드
debug_mode = false;  % true: 필터링된 이벤트의 전류 그래프를 저장
% =========================================================================

fprintf('=== Drive Cycle Event Analysis ===\n');
fprintf('Data directory: %s\n', dataDir);

% Display user selections
fprintf('\n=== User Selections ===\n');
if isempty(selectedChannels)
    fprintf('Channels: ALL (9-16)\n');
else
    fprintf('Channels: %s\n', mat2str(selectedChannels));
end
if isempty(selectedCycles)
    fprintf('Cycles: ALL available cycles\n');
else
    fprintf('Cycles: %s\n', strjoin(selectedCycles, ', '));
end
if isempty(selectedSOC)
    fprintf('SOC: ALL (SOC90, SOC70, SOC50)\n');
else
    fprintf('SOC: %s\n', strjoin(selectedSOC, ', '));
end

%% Auto-detect available cycles
fprintf('\n=== Auto-detecting available cycles ===\n');
matFiles = dir(fullfile(dataDir, 'parsedDriveCycle_*cyc_filtered.mat'));
availableCycles = {};

for i = 1:length(matFiles)
    filename = matFiles(i).name;
    match = regexp(filename, 'parsedDriveCycle_(\d+cyc)_filtered\.mat', 'tokens');
    if ~isempty(match)
        cycleType = match{1}{1};
        if ~ismember(cycleType, availableCycles)
            availableCycles{end+1} = cycleType;
        end
    end
end

% Sort cycles numerically
cycleNumbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), availableCycles);
[sortedNumbers, sortIdx] = sort(cycleNumbers);
availableCycles = availableCycles(sortIdx);

fprintf('Available cycles: %s\n', strjoin(availableCycles, ', '));

% Filter cycles based on user selection
if isempty(selectedCycles)
    cycleTypes = availableCycles;
    fprintf('Processing all available cycles\n');
else
    cycleTypes = intersect(selectedCycles, availableCycles);
    if isempty(cycleTypes)
        error('No valid cycles found in selectedCycles! Available: %s', strjoin(availableCycles, ', '));
    end
    fprintf('Processing selected cycles: %s\n', strjoin(cycleTypes, ', '));
end

% Debug mode: Process only first cycle (all channels/SOC use same drive cycle)
if debug_mode
    if ~isempty(cycleTypes)
        cycleTypes = cycleTypes(1);  % Keep only first cycle
    end
end

%% Settings (64Ah cell)
dt_list = [1, 3, 5, 10, 30, 60];
Cnom = 64;    % Battery capacity [Ah]
Pnom = 0.235; % Nominal Power    [kW] 3.68V * 64Ah = 235.52W
current_threshold = Cnom * 0.01;  % Idle current threshold [A]
target_smoothing_time = 3;        % Target smoothing duration [s] for moving average (time-based, not sample-based)

% Charging Event Detection Parameters
min_duration_chg  = 10;            % Minimum driving duration for charging events [s] (relaxed from 10s based on debug analysis)
max_P_std_chg     = Pnom * 0.1;   % Maximum power standard deviation for charging events [kW] (relaxed from 0.05 to match discharging)
max_I_std_chg     = 2;          % Maximum current standard deviation for charging events [A] (relaxed from 10A based on debug analysis)
max_CV_chg        = 2;         % Maximum coefficient of variation for charging events (relaxed from 0.70 to match discharging)

% Discharging Event Detection Parameters
min_duration_dchg = 30;           % Minimum driving duration for discharging events [s]
max_P_std_dchg    = Pnom * 0.1;   % Maximum power standard deviation for discharging events [kW]
max_I_std_dchg    = 0.1;          % Maximum current standard deviation for discharging events [A]
max_CV_dchg       = 0.1;         % Maximum coefficient of variation for discharging events

fprintf('\n=== Analysis Parameters ===\n');
fprintf('Current threshold: %.2f A\n', current_threshold);
fprintf('Moving average smoothing time: %.1f s (window size calculated based on dt)\n', target_smoothing_time);
fprintf('\n--- Charging Event Parameters ---\n');
fprintf('Min duration: %d s\n', min_duration_chg);
fprintf('Max P std: %.3f kW\n', max_P_std_chg);
fprintf('Max I std: %.2f A (%.1f%%)\n', max_I_std_chg, (max_I_std_chg/Cnom)*100);
fprintf('Max CV: %.1f%%\n', max_CV_chg*100);
fprintf('\n--- Discharging Event Parameters ---\n');
fprintf('Min duration: %d s\n', min_duration_dchg);
fprintf('Max P std: %.3f kW\n', max_P_std_dchg);
fprintf('Max I std: %.2f A (%.1f%%)\n', max_I_std_dchg, (max_I_std_dchg/Cnom)*100);
fprintf('Max CV: %.1f%%\n', max_CV_dchg*100);

% Create debug folder if debug mode is enabled
if debug_mode
    if ~exist('figures', 'dir')
        mkdir('figures');
    end
    if ~exist('figures/debug', 'dir')
        mkdir('figures/debug');
    end
end

%% Process each cycle
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    fprintf('\n\n');
    fprintf('========================================\n');
    fprintf('=== Processing %s ===\n', cycleType);
    fprintf('========================================\n');
    
    % Load data for this cycle
    dataFile = fullfile(dataDir, sprintf('parsedDriveCycle_%s_filtered.mat', cycleType));
    
    if ~exist(dataFile, 'file')
        fprintf('WARNING: File not found: %s\n', dataFile);
        fprintf('Skipping %s cycle...\n', cycleType);
        continue;
    end
    
    fprintf('Loading: %s\n', dataFile);
    try
        load(dataFile);
    catch ME
        fprintf('ERROR: Failed to load file: %s\n', dataFile);
        fprintf('Error message: %s\n', ME.message);
        continue;
    end
    
    % Get data variable name
    dataVarName = sprintf('parsedDriveCycle_%s', cycleType);
    if ~exist(dataVarName, 'var')
        fprintf('WARNING: Variable %s not found in file\n', dataVarName);
        fprintf('Skipping %s cycle...\n', cycleType);
        continue;
    end
    
    eval(sprintf('data_var = %s;', dataVarName));
    
    % Check data structure
    if ~isstruct(data_var)
        fprintf('WARNING: Data is not a structure\n');
        fprintf('Skipping %s cycle...\n', cycleType);
        continue;
    end
    
    % Initialize result structure for this cycle
    resultStructName = sprintf('Lab_DC_DCIR_%s', cycleType);
    eval(sprintf('%s = struct();', resultStructName));
    
    % Initialize filtering statistics structure for this cycle
    % Structure: filterStats.(channelName).(socLevel).(dcProfile).(rejection_type) = count
    filterStats = struct();
    
    % Get available channels
    allChannels = fieldnames(data_var);
    fprintf('Available channels: %s\n', strjoin(allChannels, ', '));
    
    % Filter channels based on user selection
    if isempty(selectedChannels)
        channels = allChannels;
        fprintf('Processing all available channels\n');
    else
        % Extract channel numbers from selectedChannels
        if isnumeric(selectedChannels)
            selectedChannelNums = selectedChannels;
        else
            selectedChannelNums = [];
            for i = 1:length(selectedChannels)
                match = regexp(selectedChannels{i}, 'ch(\d+)_', 'tokens');
                if ~isempty(match)
                    selectedChannelNums = [selectedChannelNums, str2double(match{1}{1})];
                end
            end
        end
        
        % Match channels by number
        channels = {};
        for i = 1:length(allChannels)
            chName = allChannels{i};
            match = regexp(chName, 'ch(\d+)_', 'tokens');
            if ~isempty(match)
                chNum = str2double(match{1}{1});
                if ismember(chNum, selectedChannelNums)
                    channels{end+1} = chName;
                end
            end
        end
        
        if isempty(channels)
            fprintf('WARNING: No valid channels found in selectedChannels!\n');
            fprintf('Available channels: %s\n', strjoin(allChannels, ', '));
            fprintf('Skipping %s cycle...\n', cycleType);
            continue;
        end
        fprintf('Processing selected channels: %s\n', strjoin(channels, ', '));
    end

    for ch_idx = 1:length(channels)
        channelName = channels{ch_idx};
        
        % Check if channel exists in data
        if ~isfield(data_var, channelName)
            fprintf('\nWARNING: Channel %s not found in data, skipping...\n', channelName);
            continue;
        end
        
        channel_data = data_var.(channelName);
        
        % Check if channel_data is valid
        if ~isstruct(channel_data)
            fprintf('\nWARNING: Channel %s data is not a structure, skipping...\n', channelName);
            continue;
        end
        
        all_soc_levels = fieldnames(channel_data);
        
        if isempty(all_soc_levels)
            fprintf('\nWARNING: Channel %s has no SOC levels, skipping...\n', channelName);
            continue;
        end
        
        % Filter SOC levels based on user selection
        if isempty(selectedSOC)
            soc_levels = all_soc_levels;
        else
            soc_levels = intersect(selectedSOC, all_soc_levels);
            if isempty(soc_levels)
                fprintf('\nWARNING: Channel %s has no matching SOC levels!\n', channelName);
                fprintf('Available SOC: %s\n', strjoin(all_soc_levels, ', '));
                fprintf('Selected SOC: %s\n', strjoin(selectedSOC, ', '));
                continue;
            end
        end
    
        fprintf('\n=== Processing %s ===\n', channelName);
        fprintf('Available SOC levels: %s\n', strjoin(all_soc_levels, ', '));
        if ~isempty(selectedSOC)
            fprintf('Selected SOC levels: %s\n', strjoin(soc_levels, ', '));
        end
    
        chg_struct_name = [channelName '_ChgEvent'];
        dchg_struct_name = [channelName '_DchgEvent'];
        eval(sprintf('if ~isfield(%s, chg_struct_name); %s.(chg_struct_name) = struct(); end', resultStructName, resultStructName));
        eval(sprintf('if ~isfield(%s, dchg_struct_name); %s.(dchg_struct_name) = struct(); end', resultStructName, resultStructName));
    
        for soc_idx = 1:length(soc_levels)
            soc_level = soc_levels{soc_idx};
            soc_data = channel_data.(soc_level);
            profiles = fieldnames(soc_data);
            
            fprintf('  %s: ', soc_level);
            for i = 1:length(profiles)
                fprintf('%s ', profiles{i});
            end
            fprintf('\n');
            
            eval(sprintf('if ~isfield(%s.(chg_struct_name), soc_level); %s.(chg_struct_name).(soc_level) = struct(); end', resultStructName, resultStructName));
            eval(sprintf('if ~isfield(%s.(dchg_struct_name), soc_level); %s.(dchg_struct_name).(soc_level) = struct(); end', resultStructName, resultStructName));
            
            for prof_idx = 1:length(profiles)
                profile_name = profiles{prof_idx};
                profile_data = soc_data.(profile_name);
                
                eval(sprintf('if ~isfield(%s.(chg_struct_name).(soc_level), ''%s''); %s.(chg_struct_name).(soc_level).(''%s'') = struct(); end', resultStructName, profile_name, resultStructName, profile_name));
                eval(sprintf('if ~isfield(%s.(dchg_struct_name).(soc_level), ''%s''); %s.(dchg_struct_name).(soc_level).(''%s'') = struct(); end', resultStructName, profile_name, resultStructName, profile_name));
                
                V = profile_data.V;
                I = profile_data.I;
                t = profile_data.t;
                totalTime = profile_data.totalTime;
                stepIndex = profile_data.stepIndex;
                
                if isa(t, 'duration')
                    t_seconds = seconds(t);
                else
                    t_seconds = t;
                end
                
                % Calculate sampling interval (dt) for this DC profile
                % This is important because ch15, ch16 have different sampling rates:
                % - Before 600 cycles: 1 second interval
                % - After 600 cycles: 0.1 second interval
                if length(t_seconds) > 1
                    dt = median(diff(t_seconds));  % Use median to handle any outliers
                else
                    dt = 1.0;  % Default to 1 second if only one point
                end
            
                
                %% Apply moving average filter to remove noise
                % Calculate window size based on target smoothing time (time-based, not sample-based)
                % This ensures same smoothing duration regardless of dt
                % dt=1s: 5 points = 5 seconds, dt=0.1s: 50 points = 5 seconds
                moving_avg_window = max(1, round(target_smoothing_time / dt));
                
                if length(I) >= moving_avg_window
                    I = movmean(I, moving_avg_window);
                    V = movmean(V, moving_avg_window);
                else
                    % If data is shorter than window, use all data
                    I = movmean(I, length(I));
                    V = movmean(V, length(V));
                end
                
                %% Step 1: Idle -> Load transition detection
                is_idle = abs(I) < current_threshold;
                is_driving = abs(I) >= current_threshold;
                
                % Find all idle -> driving transitions
                idle_to_driving = find(is_idle(1:end-1) & is_driving(2:end));
                
                % Initialize statistics for this DC profile
                if ~isfield(filterStats, channelName)
                    filterStats.(channelName) = struct();
                end
                if ~isfield(filterStats.(channelName), soc_level)
                    filterStats.(channelName).(soc_level) = struct();
                end
                if ~isfield(filterStats.(channelName).(soc_level), profile_name)
                    filterStats.(channelName).(soc_level).(profile_name) = struct();
                    filterStats.(channelName).(soc_level).(profile_name).total_transitions = 0;
                    filterStats.(channelName).(soc_level).(profile_name).rejected_short_duration = 0;
                    filterStats.(channelName).(soc_level).(profile_name).rejected_direction = 0;
                    filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window = 0;
                    filterStats.(channelName).(soc_level).(profile_name).rejected_I_std = 0;
                    filterStats.(channelName).(soc_level).(profile_name).rejected_CV = 0;
                    filterStats.(channelName).(soc_level).(profile_name).rejected_P_std = 0;
                    filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_after_cut = 0;
                end
                
                if isempty(idle_to_driving)
                    fprintf('    %s: No transitions detected → SKIP\n', profile_name);
                    continue;  % Skip if no transition
                end
                
                %% Step 2: Event Detection (charging/discharging split)
                chg_event_count = 0;
                dchg_event_count = 0;
                rejected_short_duration = 0;
                rejected_direction = 0;
                rejected_too_short_window = 0;
                rejected_filtering = 0;
                rejected_too_short_after_cut = 0;
                
                filterStats.(channelName).(soc_level).(profile_name).total_transitions = length(idle_to_driving);
                
                % Debug: Count charging vs discharging transitions
                num_charging_transitions = 0;
                num_discharging_transitions = 0;
                
                for i = 1:length(idle_to_driving)
                    idx1 = idle_to_driving(i);
                    start_driving_idx = idx1 + 1;
                    
                    % Determine event type from current direction at start
                    event_type = sign(I(start_driving_idx));  % 1: charge, -1: discharge
                    if event_type == 0
                        continue;  % Skip if current is exactly zero
                    end
                    
                    % Count transition types for debugging
                    if event_type > 0
                        num_charging_transitions = num_charging_transitions + 1;
                    else
                        num_discharging_transitions = num_discharging_transitions + 1;
                    end
                    
                    % Find end index of driving with same current direction
                    driving_end_idx = start_driving_idx;
                    while driving_end_idx <= length(I)
                        % Check if still driving (abs >= threshold) AND same direction
                        if abs(I(driving_end_idx)) >= current_threshold && sign(I(driving_end_idx)) == event_type
                        driving_end_idx = driving_end_idx + 1;
                        else
                            break;
                        end
                    end
                    driving_end_idx = driving_end_idx - 1;
                    
                    start_idx = idx1;
                    end_idx = driving_end_idx;
                    
                    % Check driving duration (actual time)
                    driving_time = t_seconds(driving_end_idx) - t_seconds(start_driving_idx);
                    
                    % Select appropriate parameters based on event type
                    if event_type > 0
                        % Charging event
                        min_duration = min_duration_chg;
                        max_P_std = max_P_std_chg;
                        max_I_std = max_I_std_chg;
                        max_CV = max_CV_chg;
                    else
                        % Discharging event
                        min_duration = min_duration_dchg;
                        max_P_std = max_P_std_dchg;
                        max_I_std = max_I_std_dchg;
                        max_CV = max_CV_dchg;
                    end
                    
                    if driving_time < min_duration
                        rejected_short_duration = rejected_short_duration + 1;
                        filterStats.(channelName).(soc_level).(profile_name).rejected_short_duration = ...
                            filterStats.(channelName).(soc_level).(profile_name).rejected_short_duration + 1;
                        % Debug plotting: Rejected due to short duration
                        if debug_mode
                            t_seg_temp = t_seconds(start_idx:end_idx);
                            I_seg_temp = I(start_idx:end_idx);
                            figure('Visible', 'off');
                            plot(t_seg_temp - t_seg_temp(1), I_seg_temp, 'b-', 'LineWidth', 1.5);
                            xlabel('Time (s)');
                            ylabel('Current (A)');
                            title(sprintf('Rejected: Short Duration (Val=%.2fs, Min=%.2fs)', driving_time, min_duration));
                            grid on;
                            % Ensure debug folder exists
                            debugDir = fullfile(outputDir, 'figures', 'debug');
                            if ~exist(debugDir, 'dir')
                                mkdir(debugDir);
                            end
                            filename = fullfile(debugDir, sprintf('%s_%s_%s_event%d_ShortDuration.fig', ...
                                channelName, soc_level, profile_name, i));
                            saveas(gcf, filename);
                            close(gcf);
                        end
                        continue;
                    end
                    
                    % Extract segment
                    t_seg = t_seconds(start_idx:end_idx);
                    I_seg = I(start_idx:end_idx);
                    V_seg = V(start_idx:end_idx);
                    
                    % Verify current direction consistency in the driving segment only
                    % Exclude first index (idle state): check from start_driving_idx to end_idx
                    I_driving = I(start_driving_idx:end_idx);  % Driving segment only (exclude idle state)
                    
                    % Charging event: majority of currents should be positive
                    % Discharging event: majority of currents should be negative
                    % Note: After moving average filtering, some points near zero may occur
                    % So we check majority direction rather than all points
                    if event_type > 0
                        % Charging event: majority of currents should be > 0
                        % Allow some tolerance for filtered data near zero
                        positive_ratio = sum(I_driving > 0) / length(I_driving);
                        if positive_ratio < 0.5  % At least 50% of points should be positive (relaxed from 70%)
                            rejected_direction = rejected_direction + 1;
                            filterStats.(channelName).(soc_level).(profile_name).rejected_direction = ...
                                filterStats.(channelName).(soc_level).(profile_name).rejected_direction + 1;
                            % Current is not predominantly positive in charging event, skip
                            if debug_mode
                                figure('Visible', 'off');
                                plot(t_seg - t_seg(1), I_seg, 'b-', 'LineWidth', 1.5);
                                hold on;
                                plot(t_seg(2:end) - t_seg(1), I_driving, 'r-', 'LineWidth', 2);
                                xlabel('Time (s)');
                                ylabel('Current (A)');
                                title(sprintf('Rejected: Charging Event Has Low Positive Ratio (%.1f%%, Min=%.2fA)', ...
                                    positive_ratio*100, min(I_driving)));
                                grid on;
                                legend('Full Segment', 'Driving Segment', 'Location', 'best');
                                debugDir = fullfile(outputDir, 'figures', 'debug');
                                if ~exist(debugDir, 'dir')
                                    mkdir(debugDir);
                                end
                                filename = fullfile(debugDir, sprintf('%s_%s_%s_event%d_ChargingNonPositive.fig', ...
                                    channelName, soc_level, profile_name, i));
                                saveas(gcf, filename);
                                close(gcf);
                            end
                            continue;
                        end
                    elseif event_type < 0
                        % Discharging event: all currents must be < 0
                        if any(I_driving >= 0)
                            rejected_direction = rejected_direction + 1;
                            filterStats.(channelName).(soc_level).(profile_name).rejected_direction = ...
                                filterStats.(channelName).(soc_level).(profile_name).rejected_direction + 1;
                            % Current is not negative in discharging event, skip
                            if debug_mode
                                figure('Visible', 'off');
                                plot(t_seg - t_seg(1), I_seg, 'b-', 'LineWidth', 1.5);
                                hold on;
                                plot(t_seg(2:end) - t_seg(1), I_driving, 'r-', 'LineWidth', 2);
                                xlabel('Time (s)');
                                ylabel('Current (A)');
                                title(sprintf('Rejected: Discharging Event Has Non-Negative Current (Max=%.2fA)', max(I_driving)));
                                grid on;
                                legend('Full Segment', 'Driving Segment', 'Location', 'best');
                                debugDir = fullfile(outputDir, 'figures', 'debug');
                                if ~exist(debugDir, 'dir')
                                    mkdir(debugDir);
                                end
                                filename = fullfile(debugDir, sprintf('%s_%s_%s_event%d_DischargingNonNegative.fig', ...
                                    channelName, soc_level, profile_name, i));
                                saveas(gcf, filename);
                                close(gcf);
                            end
                            continue;
                        end
                    end
                    
                    P_seg = V_seg .* I_seg / 1000;  % [kW]
                    
                    % dt is already calculated at DC profile level (above)
                    % idx2: end of driving segment
                    idx2 = driving_end_idx - start_idx + 1;
                    
                    %% Step 3: Stability check
                    % Use time-based window: same time range regardless of sampling interval
                    % This ensures same time range for 1s and 0.1s interval data
                    % Note: t_seg(1) is idle state, t_seg(2) is driving start
                    % Charging: 3-10 seconds from driving start (t_seg(2))
                    % Discharging: 2-12 seconds from driving start (t_seg(2))
                    if event_type > 0
                        % Charging event: 3-10 seconds from driving start (fixed window)
                        stable_start_time = 3;   % Start time [s] from driving start
                        stable_end_time = 10;   % End time [s] from driving start
                        
                        % Calculate indices based on actual time values (not dt-based round)
                        % t_seg(1) is idle state, t_seg(2) is driving start (index 2)
                        % Use actual time to find indices: ensures exact time range regardless of dt variations
                        if length(t_seg) >= 2
                            driving_start_idx = 2;  % Driving starts at index 2
                            driving_start_time = t_seg(2);  % Actual time of driving start
                            
                            % Find indices by actual time values (more accurate than round(time/dt))
                            target_start_time = driving_start_time + stable_start_time;
                            target_end_time = driving_start_time + stable_end_time;
                            
                            % Find first index >= target time (ensures we use time >= target, not just closest)
                            idx_start_candidates = find(t_seg >= target_start_time, 1);
                            idx_end_candidates = find(t_seg >= target_end_time, 1);
                            
                            % If no index found (target time beyond data), use last valid index
                            if isempty(idx_start_candidates)
                                stable_start_idx = length(t_seg);
                            else
                                stable_start_idx = idx_start_candidates;
                            end
                            
                            if isempty(idx_end_candidates)
                                stable_end_idx = length(t_seg);
                            else
                                stable_end_idx = idx_end_candidates;
                            end
                            
                            % Ensure indices are at least driving_start_idx and valid
                            stable_start_idx = max(driving_start_idx, min(stable_start_idx, length(t_seg)));
                            stable_end_idx = max(stable_start_idx + 1, min(stable_end_idx, length(t_seg)));
                        else
                            % Not enough points, reject this event
                            stable_start_idx = [];
                            stable_end_idx = [];
                        end
                        
                        % Check if event is long enough for stable window (time-based, not point-based)
                        % Need at least stable_end_time seconds from driving start
                        % driving_time is already calculated above: t_seconds(driving_end_idx) - t_seconds(start_driving_idx)
                        % This check must come FIRST, before point-based checks, to ensure same filtering for different dt
                        if driving_time < stable_end_time
                            rejected_too_short_window = rejected_too_short_window + 1;
                            filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window = ...
                                filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window + 1;
                            % Event is shorter than required time, skip
                            if debug_mode
                                t_seg_temp = t_seconds(start_idx:end_idx);
                                I_seg_temp = I(start_idx:end_idx);
                                figure('Visible', 'off');
                                plot(t_seg_temp - t_seg_temp(1), I_seg_temp, 'b-', 'LineWidth', 1.5);
                                xlabel('Time (s)');
                                ylabel('Current (A)');
                                title(sprintf('Rejected: Too Short for Stable Window (Time=%.2fs, Need>=%.2fs)', ...
                                    driving_time, stable_end_time));
                                grid on;
                                debugDir = fullfile(outputDir, 'figures', 'debug');
                                if ~exist(debugDir, 'dir')
                                    mkdir(debugDir);
                                end
                                filename = fullfile(debugDir, sprintf('%s_%s_%s_event%d_TooShortForStableWindow.fig', ...
                                    channelName, soc_level, profile_name, i));
                                saveas(gcf, filename);
                                close(gcf);
                            end
                            continue;
                        end
                        
                        % Adjust indices to valid range if needed (time check already passed)
                        % If calculated indices exceed array, adjust to available range
                        % This ensures we use available data even if event is exactly at minimum time
                        stable_start_idx = max(2, min(stable_start_idx, length(t_seg)));
                        stable_end_idx = max(stable_start_idx + 1, min(stable_end_idx, length(t_seg)));
                        
                        % Ensure stable_end_idx is at least stable_start_idx + minimum points for std calculation
                        % Need at least 3 points for meaningful std calculation
                        min_points_for_std = 3;
                        if stable_end_idx < stable_start_idx + min_points_for_std - 1
                            % If adjusted window is too small, try to extend stable_end_idx if possible
                            stable_end_idx = min(stable_start_idx + min_points_for_std - 1, length(t_seg));
                        end
                        
                        % Final check: ensure stable_start_idx < stable_end_idx (at least 1 point difference)
                        if stable_start_idx >= stable_end_idx
                            % This should not happen if time check passed, but handle gracefully
                            rejected_too_short_window = rejected_too_short_window + 1;
                            filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window = ...
                                filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window + 1;
                            continue;
                        end
                        
                        % Find indices corresponding to stable window (time-based, adjusted to available range)
                        stable_indices = stable_start_idx:stable_end_idx;
                        
                        if length(stable_indices) >= 3  % Need at least 3 points for std
                            % Calculate statistics for stable window
                            I_stable = I_seg(stable_indices);
                            power_std   = std(P_seg(stable_indices));
                            current_std = std(I_stable);   % Current variation (std) [A]
                            current_mean = mean(I_stable); % Average current [A]
                            current_range = max(I_stable) - min(I_stable); % Current range (max - min) [A]
                            
                            % Coefficient of Variation (CV = std/mean)
                            if abs(current_mean) > 1e-6  % Avoid division by zero
                                current_CV = current_std / abs(current_mean);  % CV (dimensionless)
                            else
                                current_CV = Inf;  % If mean is too small, CV is infinite
                            end
                        else
                            % This should not happen if event is >= 12 seconds
                            % But handle it gracefully
                            if debug_mode
                                fprintf('WARNING: Stable window has < 3 points for %s %s %s event %d\n', ...
                                    channelName, soc_level, profile_name, i);
                            end
                            power_std   = NaN;
                            current_std = NaN;
                            current_mean = NaN;
                            current_range = NaN;
                            current_CV = NaN;
                        end
                    else
                        % Discharging event: Find stable regions from 3 seconds to event end
                        % Use sliding window to find all stable regions, then use the last one
                        if length(t_seg) >= 2
                            driving_start_idx = 2;  % Driving starts at index 2
                            driving_start_time = t_seg(2);  % Actual time of driving start
                            
                            % Start checking from 3 seconds after driving start
                            stable_start_time = 3;  % Start time [s] from driving start
                            target_start_time = driving_start_time + stable_start_time;
                            idx_start_candidates = find(t_seg >= target_start_time, 1);
                            
                            if isempty(idx_start_candidates)
                                % Event is too short (less than 3 seconds), reject
                                rejected_too_short_window = rejected_too_short_window + 1;
                                filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window = ...
                                    filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window + 1;
                                continue;
                            end
                            
                            check_start_idx = idx_start_candidates;
                            check_start_idx = max(driving_start_idx, min(check_start_idx, length(t_seg)));
                            
                            % Use sliding window to find stable regions
                            % Window duration: 3 seconds (time-based)
                            window_duration = 3;  % Window duration [s]
                            if length(t_seg) > 1
                                window_num_points = round(window_duration / dt);  % Convert to points based on dt
                            else
                                window_num_points = 3;  % Default
                            end
                            
                            % Minimum points required for stability check
                            min_check_duration = 2;  % Minimum check duration [s] for stability
                            min_window_points = max(3, round(min_check_duration / dt));  % Time-based, minimum 3 points
                            
                            % Find all stable regions from check_start_idx to end
                            stable_regions = [];  % Store [start_idx, end_idx] pairs
                            current_stable_start = [];
                            
                            for check_idx = check_start_idx:length(t_seg)
                                % Define window ending at check_idx
                                window_start_idx = max(2, check_idx - window_num_points + 1);
                                
                                if window_start_idx < check_start_idx
                                    continue;  % Window starts before check_start_idx
                                end
                                
                                window_indices = window_start_idx:check_idx;
                                
                                if length(window_indices) < min_window_points
                                    continue;  % Not enough points
                                end
                                
                                % Check stability criteria in this window
                                I_window = I_seg(window_indices);
                                P_window = P_seg(window_indices);
                                
                                window_power_std = std(P_window);
                                window_current_std = std(I_window);
                                window_current_mean = mean(I_window);
                                
                                if abs(window_current_mean) > 1e-6
                                    window_CV = window_current_std / abs(window_current_mean);
                                else
                                    window_CV = Inf;
                                end
                                
                                % Check if this window is stable
                                is_stable = true;
                                if window_power_std >= max_P_std || ...
                                   window_current_std >= max_I_std || ...
                                   window_CV >= max_CV
                                    is_stable = false;
                                end
                                
                                % Also check current direction and threshold
                                I_current = I_seg(check_idx);
                                if sign(I_current) ~= event_type || abs(I_current) < current_threshold
                                    is_stable = false;
                                end
                                
                                if is_stable
                                    % Start or continue stable region
                                    if isempty(current_stable_start)
                                        current_stable_start = window_start_idx;
                                    end
                                else
                                    % End of stable region
                                    if ~isempty(current_stable_start)
                                        stable_regions(end+1, :) = [current_stable_start, check_idx - 1];
                                        current_stable_start = [];
                                    end
                                end
                            end
                            
                            % If still in a stable region at the end, close it
                            if ~isempty(current_stable_start)
                                stable_regions(end+1, :) = [current_stable_start, length(t_seg)];
                            end
                            
                            % Use the last stable region for statistics
                            if ~isempty(stable_regions)
                                last_stable_region = stable_regions(end, :);
                                stable_start_idx = last_stable_region(1);
                                stable_end_idx = last_stable_region(2);
                                
                                stable_indices = stable_start_idx:stable_end_idx;
                                
                                if length(stable_indices) >= 3
                                    % Calculate statistics for last stable window
                                    I_stable = I_seg(stable_indices);
                                    power_std   = std(P_seg(stable_indices));
                                    current_std = std(I_stable);   % Current variation (std) [A]
                                    current_mean = mean(I_stable); % Average current [A]
                                    current_range = max(I_stable) - min(I_stable); % Current range (max - min) [A]
                                    
                                    % Coefficient of Variation (CV = std/mean)
                                    if abs(current_mean) > 1e-6  % Avoid division by zero
                                        current_CV = current_std / abs(current_mean);  % CV (dimensionless)
                                    else
                                        current_CV = Inf;  % If mean is too small, CV is infinite
                                    end
                                else
                                    % No valid stable region found
                                    rejected_too_short_window = rejected_too_short_window + 1;
                                    filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window = ...
                                        filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window + 1;
                                    continue;
                                end
                            else
                                % No stable region found
                                rejected_too_short_window = rejected_too_short_window + 1;
                                filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window = ...
                                    filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window + 1;
                                continue;
                            end
                        else
                            % Not enough points, reject this event
                            rejected_too_short_window = rejected_too_short_window + 1;
                            filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window = ...
                                filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_window + 1;
                            continue;
                        end
                    end
                    
                    % Check multiple criteria: power std, current std, and CV
                    reject_reasons = {};
                    rejected_I_std = false;
                    rejected_CV = false;
                    rejected_P_std = false;
                    
                    if power_std >= max_P_std
                        reject_reasons{end+1} = sprintf('High P Std (%.3f >= %.3f)', power_std, max_P_std);
                        rejected_P_std = true;
                        filterStats.(channelName).(soc_level).(profile_name).rejected_P_std = ...
                            filterStats.(channelName).(soc_level).(profile_name).rejected_P_std + 1;
                    end
                    if current_std >= max_I_std
                        reject_reasons{end+1} = sprintf('High I Std (%.3f >= %.3f)', current_std, max_I_std);
                        rejected_I_std = true;
                        filterStats.(channelName).(soc_level).(profile_name).rejected_I_std = ...
                            filterStats.(channelName).(soc_level).(profile_name).rejected_I_std + 1;
                    end
                    if current_CV >= max_CV
                        reject_reasons{end+1} = sprintf('High CV (%.3f >= %.3f)', current_CV, max_CV);
                        rejected_CV = true;
                        filterStats.(channelName).(soc_level).(profile_name).rejected_CV = ...
                            filterStats.(channelName).(soc_level).(profile_name).rejected_CV + 1;
                    end
                    
                    if ~isempty(reject_reasons)
                        rejected_filtering = rejected_filtering + 1;
                        % Debug plotting: Rejected due to filtering criteria
                        if debug_mode
                            figure('Visible', 'off');
                            plot(t_seg - t_seg(1), I_seg, 'b-', 'LineWidth', 1.5);
                            hold on;
                            % Mark stable window
                            stable_t = t_seg(stable_indices) - t_seg(1);
                            stable_I = I_seg(stable_indices);
                            plot(stable_t, stable_I, 'r-', 'LineWidth', 2);
                            xlabel('Time (s)');
                            ylabel('Current (A)');
                            title_str = sprintf('Rejected: %s', strjoin(reject_reasons, ', '));
                            title(title_str, 'FontSize', 10);
                            if event_type > 0
                                legend('Full Event', sprintf('Stable Window (points %d-%d)', stable_start_idx, stable_end_idx), 'Location', 'best');
                            else
                                legend('Full Event', sprintf('Stable Window (points %d-%d)', stable_start_idx, stable_end_idx), 'Location', 'best');
                            end
                            grid on;
                            % Add statistics text
                            stats_text = sprintf('I_std=%.3f (limit=%.3f)\nI_range=%.3f\nCV=%.3f (limit=%.3f)', ...
                                current_std, max_I_std, current_range, current_CV, max_CV);
                            text(0.02, 0.98, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                                'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
                            hold off;
                            % Ensure debug folder exists
                            debugDir = fullfile(outputDir, 'figures', 'debug');
                            if ~exist(debugDir, 'dir')
                                mkdir(debugDir);
                            end
                            filename = fullfile(debugDir, sprintf('%s_%s_%s_event%d_Rejected.fig', ...
                                channelName, soc_level, profile_name, i));
                            saveas(gcf, filename);
                            close(gcf);
                        end
                        continue;
                    end                   
                    
                    %% Step 4: Find stable region end point using sliding window
                    if event_type > 0
                        % Charging event: Find unstable region and cut before it
                        % Use stable window average current as reference
                        I_ref = current_mean;  % Reference current from stable window
                        I_ref_tolerance = 0.10;  % ±10% tolerance
                        I_ref_lower = I_ref * (1 - I_ref_tolerance);
                        I_ref_upper = I_ref * (1 + I_ref_tolerance);
                        
                        % Find unstable region start point (forward search)
                        unstable_start_idx = length(t_seg) + 1;  % Default: no unstable region found
                        
                        % Use time-based window: same time range regardless of sampling interval
                        % Window duration: 3 seconds
                        window_duration = 3;  % Window duration [s]
                        % Calculate window_num_points based on actual dt
                        if length(t_seg) > 1
                            window_num_points = round(window_duration / dt);  % Convert to points based on dt
                        else
                            window_num_points = 3;  % Default
                        end
                        
                        % Start checking from 3 seconds from driving start
                        check_start_time = 3;  % Start time [s] from driving start
                        
                        % Calculate index based on actual time values
                        if length(t_seg) >= 2
                            driving_start_idx = 2;  % Driving starts at index 2
                            driving_start_time = t_seg(2);  % Actual time of driving start
                            
                            % Check if event has enough time first
                            if driving_time < check_start_time
                                continue;  % Event too short in time, skip
                            end
                            
                            % Find index by actual time value
                            target_check_time = driving_start_time + check_start_time;
                            idx_check_candidates = find(t_seg >= target_check_time, 1);
                            if isempty(idx_check_candidates)
                                check_start_idx = length(t_seg);
                            else
                                check_start_idx = idx_check_candidates;
                            end
                            check_start_idx = max(driving_start_idx, min(check_start_idx, length(t_seg)));
                        else
                            % Not enough points, skip this event
                            continue;
                        end
                        
                        % Calculate min_start_idx once before loop (time-based, using actual time)
                        min_start_time = 3;  % Minimum start time [s] from driving start
                        target_min_time = driving_start_time + min_start_time;
                        idx_min_candidates = find(t_seg >= target_min_time, 1);
                        if isempty(idx_min_candidates)
                            min_start_idx = length(t_seg);
                        else
                            min_start_idx = idx_min_candidates;
                        end
                        min_start_idx = max(driving_start_idx, min(min_start_idx, length(t_seg)));
                        
                        if check_start_idx < length(t_seg)
                            % Forward search: Find where instability starts
                            for check_idx = check_start_idx:length(t_seg)
                                % Define window ending at check_idx
                                window_start_idx = max(2, check_idx - window_num_points + 1);
                                
                                % Window should not start before the minimum start point
                                if window_start_idx < min_start_idx
                                    continue;  % Window starts before minimum time, skip
                                end
                                
                                % Find indices in this window
                                window_indices = window_start_idx:check_idx;
                                
                                % Minimum points required
                                min_check_duration = 2;  % Minimum check duration [s] for stability
                                min_window_points = max(3, round(min_check_duration / dt));  % Time-based, minimum 3 points
                                
                                if length(window_indices) < min_window_points
                                    continue;  % Not enough points
                                end
                                
                                % Check termination conditions
                                I_current = I_seg(check_idx);
                                
                                % Condition 1: Current direction changed
                                if sign(I_current) ~= event_type
                                    unstable_start_idx = check_idx;
                                    break;
                                end
                                
                                % Condition 2: Current below threshold
                                if abs(I_current) < current_threshold
                                    unstable_start_idx = check_idx;
                                    break;
                                end
                                
                                % Check stability criteria in this window
                                I_window = I_seg(window_indices);
                                P_window = P_seg(window_indices);
                                
                                window_power_std = std(P_window);
                                window_current_std = std(I_window);
                                window_current_mean = mean(I_window);
                                
                                if abs(window_current_mean) > 1e-6
                                    window_CV = window_current_std / abs(window_current_mean);
                                else
                                    window_CV = Inf;
                                end
                                
                                % Check if this window is unstable
                                is_unstable = false;
                                if window_power_std >= max_P_std || ...
                                   window_current_std >= max_I_std || ...
                                   window_CV >= max_CV
                                    is_unstable = true;
                                end
                                
                                % Also check if current is outside ±10% of reference
                                if I_current < I_ref_lower || I_current > I_ref_upper
                                    is_unstable = true;
                                end
                                
                                % Check current change rate (slope)
                                if ~is_unstable
                                    % Charging: current should not be too much higher than reference
                                    if I_current > I_ref + 0.2
                                        is_unstable = true;
                                    end
                                    
                                    % Also check recent trend
                                    if check_idx > 1
                                        lookback_duration = 5;  % Look back 5 seconds [s]
                                        num_lookback_points = max(3, round(lookback_duration / dt));
                                        prev_start_idx = max(1, check_idx - num_lookback_points);
                                        prev_indices = prev_start_idx:(check_idx - 1);
                                        if length(prev_indices) >= 3
                                            I_prev_avg = mean(I_seg(prev_indices));
                                            I_trend_diff = abs(I_current - I_prev_avg);
                                            if I_trend_diff > 0.2
                                                is_unstable = true;
                                            end
                                        end
                                    end
                                end
                                
                                % If unstable, immediately mark as unstable region start
                                if is_unstable
                                    unstable_start_idx = check_idx;
                                    break;
                                end
                            end
                        end
                        
                        % Backward search: Find last stable point before unstable region
                        if unstable_start_idx <= length(t_seg)
                            stable_end_idx = [];  % Initialize as empty
                            
                            search_back_start = unstable_start_idx - 1;
                            min_search_back_time = 3;  % [s] from driving start
                            if length(t_seg) >= 2
                                driving_start_time = t_seg(2);
                                target_min_search_time = driving_start_time + min_search_back_time;
                                min_search_back_idx = find(t_seg >= target_min_search_time, 1);
                                if isempty(min_search_back_idx)
                                    min_search_back_idx = 2;
                                end
                            else
                                min_search_back_idx = 2;
                            end
                            
                            % Search backward from unstable region to find last stable point
                            for check_idx = search_back_start:-1:min_search_back_idx
                                window_start_idx = max(2, check_idx - window_num_points + 1);
                                
                                if window_start_idx < 2
                                    break;
                                end
                                
                                window_indices = window_start_idx:check_idx;
                                
                                min_check_duration = 2;
                                min_window_points = max(3, round(min_check_duration / dt));
                                
                                if length(window_indices) < min_window_points
                                    continue;
                                end
                                
                                % Check stability criteria in this window
                                I_window = I_seg(window_indices);
                                P_window = P_seg(window_indices);
                                
                                window_power_std = std(P_window);
                                window_current_std = std(I_window);
                                window_current_mean = mean(I_window);
                                
                                if abs(window_current_mean) > 1e-6
                                    window_CV = window_current_std / abs(window_current_mean);
                                else
                                    window_CV = Inf;
                                end
                                
                                % Check if this window is stable
                                is_stable = true;
                                if window_power_std >= max_P_std || ...
                                   window_current_std >= max_I_std || ...
                                   window_CV >= max_CV
                                    is_stable = false;
                                end
                                
                                I_current = I_seg(check_idx);
                                if I_current < I_ref_lower || I_current > I_ref_upper
                                    is_stable = false;
                                end
                                
                                if sign(I_current) ~= event_type || abs(I_current) < current_threshold
                                    is_stable = false;
                                end
                                
                                if is_stable
                                    stable_end_idx = check_idx;
                                    break;
                                end
                            end
                            
                            % If no stable point found in backward search
                            if isempty(stable_end_idx)
                                stable_end_idx = unstable_start_idx - 1;
                                if stable_end_idx < check_start_idx
                                    stable_end_idx = check_start_idx;
                                end
                            end
                        else
                            % No unstable region found, use full event
                            stable_end_idx = length(t_seg);
                        end
                        
                        % Final validation: ensure stable_end_idx is within valid range
                        if stable_end_idx < 2
                            stable_end_idx = 2;
                        end
                        if stable_end_idx > length(t_seg)
                            stable_end_idx = length(t_seg);
                        end
                    else
                        % Discharging event: Use the last stable region end point found in Step 3
                        % stable_end_idx is already set in Step 3 as the end of the last stable region
                        % No need to search for unstable region - just use the last stable region end
                        % stable_end_idx is already correctly set from Step 3
                    end
                    
                    % Store original full event for visualization
                    t_seg_full = t_seg;  % Original full event
                    I_seg_full = I_seg;  % Original full event
                    
                    % Cut event to stable region
                    t_seg_cut = t_seg(1:stable_end_idx);
                    I_seg_cut = I_seg(1:stable_end_idx);
                    V_seg_cut = V_seg(1:stable_end_idx);
                    P_seg_cut = P_seg(1:stable_end_idx);
                    
                    % Store stable region end point time for visualization (relative to original event start)
                    % This is the actual time value at stable_end_idx in t_seg_full, relative to t_seg_full(1)
                    stable_region_end_time = t_seg_full(stable_end_idx) - t_seg_full(1);  % Time [s] relative to event start
                    % Also store the index for more accurate visualization
                    stable_region_end_idx = stable_end_idx;  % Index in t_seg_full (same as t_full)
                    
                    % Check minimum duration after cutting (use appropriate min_duration based on event_type)
                    % Note: min_duration is already set based on event_type above
                    % For charging events, use a more lenient threshold after cutting since they tend to be shorter
                    % and cutting unstable regions may reduce duration significantly
                    driving_time_cut = t_seg_cut(end) - t_seg_cut(1);
                    if event_type > 0
                        % Charging event: require at least 5 seconds after cutting (relaxed from 8s based on debug analysis)
                        min_duration_after_cut = 5;
                    else
                        % Discharging event: require at least 5 seconds after cutting (same as charging)
                        min_duration_after_cut = 5;
                    end
                    
                    if driving_time_cut < min_duration_after_cut
                        rejected_too_short_after_cut = rejected_too_short_after_cut + 1;
                        filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_after_cut = ...
                            filterStats.(channelName).(soc_level).(profile_name).rejected_too_short_after_cut + 1;
                        % Event too short after cutting, skip
                        if debug_mode
                            figure('Visible', 'off');
                            plot(t_seg_full - t_seg_full(1), I_seg_full, 'b-', 'LineWidth', 1.5);
                            hold on;
                            plot(t_seg_cut - t_seg_full(1), I_seg_cut, 'r-', 'LineWidth', 2);
                            xlabel('Time (s)');
                            ylabel('Current (A)');
                            title(sprintf('Rejected: Too Short After Cutting (Duration=%.2fs, Min=%.2fs)', ...
                                driving_time_cut, min_duration));
                            legend('Full Event', 'Cut Event', 'Location', 'best');
                            grid on;
                            debugDir = fullfile(outputDir, 'figures', 'debug');
                            if ~exist(debugDir, 'dir')
                                mkdir(debugDir);
                            end
                            filename = fullfile(debugDir, sprintf('%s_%s_%s_event%d_TooShortAfterCut.fig', ...
                                channelName, soc_level, profile_name, i));
                            saveas(gcf, filename);
                            close(gcf);
                        end
                        continue;
                    end
                    
                    % Update segment variables to use cut version (for analysis)
                    t_seg = t_seg_cut;
                    I_seg = I_seg_cut;
                    V_seg = V_seg_cut;
                    P_seg = P_seg_cut;
                    driving_time = driving_time_cut;
                    end_idx = start_idx + stable_end_idx - 1;  % Update end_idx
                    
                    % Charging/discharging split based on event_type (determined at start)
                    % event_type was already determined at the beginning of the loop
                    if event_type > 0
                        % Charging event
                        chg_event_count = chg_event_count + 1;
                        evtName = sprintf('event%d', chg_event_count);
                        eval(sprintf('target_struct = %s.(chg_struct_name);', resultStructName));
                    elseif event_type < 0
                        % Discharging event
                        dchg_event_count = dchg_event_count + 1;
                        evtName = sprintf('event%d', dchg_event_count);
                        eval(sprintf('target_struct = %s.(dchg_struct_name);', resultStructName));
                    else
                        continue;
                    end
                    
                    % Store event in the appropriate structure
                    target_struct.(soc_level).(profile_name).(evtName).channel = channelName;
                    target_struct.(soc_level).(profile_name).(evtName).soc_level = soc_level;
                    target_struct.(soc_level).(profile_name).(evtName).profile_name = profile_name;
                    target_struct.(soc_level).(profile_name).(evtName).stepIndex = stepIndex;
                    target_struct.(soc_level).(profile_name).(evtName).event_number = str2double(evtName(6:end));
                    target_struct.(soc_level).(profile_name).(evtName).transition_idx = idx1;
                    target_struct.(soc_level).(profile_name).(evtName).driving_duration = driving_time;
                    target_struct.(soc_level).(profile_name).(evtName).t = t_seg;  % Cut event (for analysis)
                    target_struct.(soc_level).(profile_name).(evtName).I = I_seg;  % Cut event (for analysis)
                    target_struct.(soc_level).(profile_name).(evtName).V = V_seg;
                    target_struct.(soc_level).(profile_name).(evtName).P = P_seg;
                    target_struct.(soc_level).(profile_name).(evtName).t_full = t_seg_full;  % Original full event (for visualization)
                    target_struct.(soc_level).(profile_name).(evtName).I_full = I_seg_full;  % Original full event (for visualization)
                    target_struct.(soc_level).(profile_name).(evtName).stable_region_end_time = stable_region_end_time;  % Stable region end point time [s] relative to event start
                    target_struct.(soc_level).(profile_name).(evtName).stable_region_end_idx = stable_region_end_idx;  % Stable region end point index in t_full
                    target_struct.(soc_level).(profile_name).(evtName).I_std = current_std;      % Current std in stable window [A] (Chg: points 3-12, Dchg: points 2-12)
                    target_struct.(soc_level).(profile_name).(evtName).I_mean = current_mean;     % Average current in stable window [A] (Chg: points 3-12, Dchg: points 2-12)
                    target_struct.(soc_level).(profile_name).(evtName).I_range = current_range;  % Current range (max-min) in stable window [A] (Chg: points 3-12, Dchg: points 2-12)
                    target_struct.(soc_level).(profile_name).(evtName).I_CV = current_CV;        % Coefficient of variation (std/mean) in stable window (Chg: points 3-12, Dchg: points 2-12)
                    target_struct.(soc_level).(profile_name).(evtName).P_std = power_std;
                    
                    %% === [NEW] Paper-based Feature Extraction ===
                    % 논문에서 제시한 16개 + a 특징 추출
                    
                    % 1. 기본 데이터 준비
                    % 시간 및 용량 계산
                    t_rel = t_seg - t_seg(1); % 상대 시간 (0부터 시작)
                    Q_seg = cumtrapz(t_rel, I_seg) / 3600; % [Ah] 누적 용량 변화
                    
                    % 노이즈 제거 (dQ/dV 계산을 위해 필수)
                    % Window size는 데이터 길이에 따라 유동적으로 (최소 5포인트)
                    smooth_window = max(5, round(length(V_seg)/20)); 
                    V_smooth = smoothdata(V_seg, 'gaussian', smooth_window);
                    Q_smooth = smoothdata(Q_seg, 'gaussian', smooth_window);
                    
                    % 2. dQ/dV 계산 (미분)
                    dQ = diff(Q_smooth);
                    dV = diff(V_smooth);
                    
                    % dV가 0에 가까우면 Inf 발생하므로 처리
                    valid_diff = abs(dV) > 1e-5; 
                    if sum(valid_diff) > 5
                        dQdV = dQ(valid_diff) ./ dV(valid_diff);
                    else
                        dQdV = zeros(size(dQ)); % 예외 처리
                    end
                    
                    % --- Feature Group A: Geometric (기하학적 특징) ---
                    feat_Time = t_rel(end); % 충전/방전 지속 시간
                    feat_EOCV = V_seg(end); % 종료 전압
                    feat_CapThroughput = abs(Q_seg(end)); % 총 처리 용량
                    
                    % Slope (기울기 계산 - Linear Fit)
                    if length(t_rel) > 1
                        p_V = polyfit(t_rel, V_seg, 1);
                        feat_Slope_V = p_V(1); % 전압 기울기
                        
                        p_I = polyfit(t_rel, I_seg, 1);
                        feat_Slope_I = p_I(1); % 전류 기울기 (CC가 아닌 경우 중요)
                    else
                        feat_Slope_V = 0; feat_Slope_I = 0;
                    end
                    
                    % --- Feature Group B: Statistics (통계적 특징) ---
                    % 전압(V) 분포 특성
                    feat_Mean_V = mean(V_seg);
                    feat_Var_V = var(V_seg);
                    feat_Skew_V = skewness(V_seg);
                    feat_Kurt_V = kurtosis(V_seg);
                    
                    % 전류(I) 분포 특성 (FR 패턴의 격렬함 반영)
                    feat_Mean_I = mean(I_seg);
                    feat_Var_I = var(I_seg);
                    feat_Skew_I = skewness(I_seg);
                    feat_Kurt_I = kurtosis(I_seg);
                    
                    % 용량 변화(Delta Q) 분포 특성
                    feat_Skew_Q = skewness(Q_seg);
                    feat_Kurt_Q = kurtosis(Q_seg);
                    
                    % --- Feature Group C: Electrochemical (dQ/dV 통계) ---
                    % 논문의 핵심: Peak를 찾는게 아니라 분포를 본다
                    if ~isempty(dQdV) && sum(~isinf(dQdV)) > 0
                        valid_dQdV = dQdV(~isinf(dQdV));
                        feat_Mean_dQdV = mean(valid_dQdV);
                        feat_Var_dQdV = var(valid_dQdV); % **핵심**: 노화될수록 뭉개짐(Variance 변화)
                        feat_Skew_dQdV = skewness(valid_dQdV);
                        feat_Kurt_dQdV = kurtosis(valid_dQdV);
                    else
                        feat_Mean_dQdV=NaN; feat_Var_dQdV=NaN; 
                        feat_Skew_dQdV=NaN; feat_Kurt_dQdV=NaN;
                    end

                    % --- [Struct 저장] ---
                    % target_struct에 필드 추가
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Time = feat_Time;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_EOCV = feat_EOCV;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Cap = feat_CapThroughput;
                    
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Slope_V = feat_Slope_V;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Slope_I = feat_Slope_I;
                    
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Mean_V = feat_Mean_V;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Var_V = feat_Var_V;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Skew_V = feat_Skew_V;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Kurt_V = feat_Kurt_V;
                    
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Mean_I = feat_Mean_I;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Var_I = feat_Var_I;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Skew_I = feat_Skew_I;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Kurt_I = feat_Kurt_I;
                    
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Skew_Q = feat_Skew_Q;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Kurt_Q = feat_Kurt_Q;
                    
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Mean_dQdV = feat_Mean_dQdV;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Var_dQdV = feat_Var_dQdV;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Skew_dQdV = feat_Skew_dQdV;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Kurt_dQdV = feat_Kurt_dQdV;
                    
                    %% [추가] dV/dQ Features (방전 데이터 맞춤)
                    % 논문의 DV Analysis (Differential Voltage) 적용
                    % 방전 평탄 구간(Plateau)에서 dQ/dV보다 훨씬 안정적입니다.
                    
                    % dV/dQ 계산 (분모가 용량 변화량이므로 0이 될 확률이 매우 낮음)
                    % dQ가 0인 경우(충전 멈춤)만 예외 처리
                    valid_dQ = abs(dQ) > 1e-6;
                    
                    if sum(valid_dQ) > 5
                        dVdQ = dV(valid_dQ) ./ dQ(valid_dQ);
                        
                        % 이상치 제거 (너무 큰 값은 잘라냄 - Clipping)
                        % dV/dQ는 저항과 유사하므로 수천 단위로 튀면 안 됩니다.
                        limit_val = 50; 
                        dVdQ(dVdQ > limit_val) = limit_val;
                        dVdQ(dVdQ < -limit_val) = -limit_val;
                        
                        % 통계량 추출
                        feat_Mean_dVdQ = mean(dVdQ);
                        feat_Var_dVdQ = var(dVdQ);       % 저항 변동성 (핵심 지표)
                        feat_Skew_dVdQ = skewness(dVdQ); % 저항 증가 패턴의 치우침
                        feat_Kurt_dVdQ = kurtosis(dVdQ);
                    else
                        feat_Mean_dVdQ = NaN; feat_Var_dVdQ = NaN;
                        feat_Skew_dVdQ = NaN; feat_Kurt_dVdQ = NaN;
                    end
                    
                    % --- [Struct 저장] ---
                    % dV/dQ (DV Analysis)
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Mean_dVdQ = feat_Mean_dVdQ;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Var_dVdQ = feat_Var_dVdQ;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Skew_dVdQ = feat_Skew_dVdQ;
                    target_struct.(soc_level).(profile_name).(evtName).Feat_Kurt_dVdQ = feat_Kurt_dVdQ;
                    
                    for dt_idx = 1:length(dt_list)
                        dt_sec = dt_list(dt_idx);
                        if length(t_seg) > 1
                            dt = t_seg(2) - t_seg(1); % Sample interval
                            idx_dt = round(dt_sec / dt + 1); % Integer index
                        else
                            idx_dt = 1;
                        end
                        if idx_dt <= length(I_seg)
                            V1 = V_seg(1);
                            V2 = V_seg(idx_dt);
                            I1 = I_seg(1);
                            I2 = I_seg(idx_dt);
                            dV = V2 - V1;
                            dI = I2 - I1;
                            if abs(dI) > 1e-6  % Avoid division by zero
                                dcir_val = (dV / dI) * 1000;
                            else
                                dcir_val = NaN;
                            end
                        else
                            V1 = NaN; V2 = NaN; I1 = NaN; I2 = NaN;
                            dV = NaN; dI = NaN;
                            dcir_val = NaN;
                        end
                        fieldName = sprintf('DCIR_%ds', dt_sec);
                        target_struct.(soc_level).(profile_name).(evtName).(fieldName).val = dcir_val;
                        target_struct.(soc_level).(profile_name).(evtName).(fieldName).V1 = V1;
                        target_struct.(soc_level).(profile_name).(evtName).(fieldName).V2 = V2;
                        target_struct.(soc_level).(profile_name).(evtName).(fieldName).I1 = I1;
                        target_struct.(soc_level).(profile_name).(evtName).(fieldName).I2 = I2;
                        target_struct.(soc_level).(profile_name).(evtName).(fieldName).dV = dV;
                        target_struct.(soc_level).(profile_name).(evtName).(fieldName).dI = dI;
                    end
                    % Save back to main structure
                    if mean(I_seg) > 0
                        eval(sprintf('%s.(chg_struct_name) = target_struct;', resultStructName));
                    else
                        eval(sprintf('%s.(dchg_struct_name) = target_struct;', resultStructName));
                    end
                end
                
                % Print DC profile summary
                if chg_event_count > 0 || dchg_event_count > 0
                    fprintf('    %s: %d charging, %d discharging events\n', profile_name, chg_event_count, dchg_event_count);
                else
                    fprintf('    %s: No valid events detected\n', profile_name);
                end
            end
        end
    end
    
    % Save results for this cycle
    savePath = fullfile(outputDir, sprintf('Lab_DC_DCIR_%s_Events.mat', cycleType));
    eval(sprintf('save(''%s'', ''%s'');', savePath, resultStructName));
    fprintf('\n=== %s cycle processing complete ===\n', cycleType);
    fprintf('Saved: %s\n', savePath);
    
    % Detailed debug summary sections removed
end

fprintf('\n\n========================================\n');
fprintf('=== Event Analysis Complete ===\n');
fprintf('========================================\n');
fprintf('All results saved to: %s\n', outputDir);
fprintf('Event files: Lab_DC_DCIR_*cyc_Events.mat\n');

%% Event Count Summary Table
% =========================================================================
fprintf('\n\n========================================\n');
fprintf('=== Event Count Summary (by Cycle, SOC, DC) ===\n');
fprintf('========================================\n');

% SOC 레벨 및 DC 프로파일 리스트
socLevels = {'SOC90', 'SOC70', 'SOC50'};
dcProfiles = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

% 각 사이클별로 통계 수집
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    eventFile = fullfile(outputDir, sprintf('Lab_DC_DCIR_%s_Events.mat', cycleType));
    
    if ~exist(eventFile, 'file')
        fprintf('\n[%s] Event file not found: %s\n', cycleType, eventFile);
        continue;
    end
    
    % 이벤트 파일 로드
    load(eventFile);
    resultStructName = sprintf('Lab_DC_DCIR_%s', cycleType);
    eval(sprintf('eventData = %s;', resultStructName));
    
    % 채널 목록 가져오기 (이미 ChgEvent/DchgEvent 포함된 전체 필드 이름)
    channelNames = fieldnames(eventData);
    
    fprintf('\n');
    fprintf('========================================\n');
    fprintf('Cycle: %s\n', cycleType);
    fprintf('========================================\n');
    fprintf('Available channel fields: %s\n', strjoin(channelNames, ', '));
    
    % SOC별, DC별로 이벤트 개수 카운트
    for socIdx = 1:length(socLevels)
        socLevel = socLevels{socIdx};
        
        fprintf('\n--- %s ---\n', socLevel);
        fprintf('%-6s | %-10s | %-10s | %-10s\n', 'DC', 'Charge', 'Discharge', 'Total');
        fprintf('------|------------|------------|------------\n');
        
        for dcIdx = 1:length(dcProfiles)
            dcProfile = dcProfiles{dcIdx};
            
            chgCount = 0;
            dchgCount = 0;
            
            % 모든 채널에서 이벤트 개수 카운트
            for chIdx = 1:length(channelNames)
                channelFieldName = channelNames{chIdx};
                
                % 채널 필드 이름이 ChgEvent 또는 DchgEvent로 끝나는지 확인
                % Note: Support both _DchEvent (old) and _DchgEvent (new) for backward compatibility
                if contains(channelFieldName, '_ChgEvent', 'IgnoreCase', true)
                    % 충전 이벤트 카운트
                    if isfield(eventData, channelFieldName)
                        chgData = eventData.(channelFieldName);
                        if isfield(chgData, socLevel) && isfield(chgData.(socLevel), dcProfile)
                            allFields = fieldnames(chgData.(socLevel).(dcProfile));
                            eventFields = allFields(startsWith(allFields, 'event'));
                            chgCount = chgCount + length(eventFields);
                        end
                    end
                elseif contains(channelFieldName, '_DchgEvent', 'IgnoreCase', true) || contains(channelFieldName, '_DchEvent', 'IgnoreCase', true)
                    % 방전 이벤트 카운트
                    if isfield(eventData, channelFieldName)
                        dchgData = eventData.(channelFieldName);
                        if isfield(dchgData, socLevel) && isfield(dchgData.(socLevel), dcProfile)
                            allFields = fieldnames(dchgData.(socLevel).(dcProfile));
                            eventFields = allFields(startsWith(allFields, 'event'));
                            dchgCount = dchgCount + length(eventFields);
                        end
                    end
                end
            end
            
            totalCount = chgCount + dchgCount;
            fprintf('%-6s | %10d | %10d | %10d\n', dcProfile, chgCount, dchgCount, totalCount);
        end
        
        % SOC별 총합
        socChgTotal = 0;
        socDchgTotal = 0;
        for dcIdx = 1:length(dcProfiles)
            dcProfile = dcProfiles{dcIdx};
            for chIdx = 1:length(channelNames)
                channelFieldName = channelNames{chIdx};
                if contains(channelFieldName, '_ChgEvent', 'IgnoreCase', true)
                    if isfield(eventData, channelFieldName)
                        chgData = eventData.(channelFieldName);
                        if isfield(chgData, socLevel) && isfield(chgData.(socLevel), dcProfile)
                            allFields = fieldnames(chgData.(socLevel).(dcProfile));
                            eventFields = allFields(startsWith(allFields, 'event'));
                            socChgTotal = socChgTotal + length(eventFields);
                        end
                    end
                elseif contains(channelFieldName, '_DchgEvent', 'IgnoreCase', true) || contains(channelFieldName, '_DchEvent', 'IgnoreCase', true)
                    if isfield(eventData, channelFieldName)
                        dchgData = eventData.(channelFieldName);
                        if isfield(dchgData, socLevel) && isfield(dchgData.(socLevel), dcProfile)
                            allFields = fieldnames(dchgData.(socLevel).(dcProfile));
                            eventFields = allFields(startsWith(allFields, 'event'));
                            socDchgTotal = socDchgTotal + length(eventFields);
                        end
                    end
                end
            end
        end
        fprintf('------|------------|------------|------------\n');
        fprintf('%-6s | %10d | %10d | %10d\n', 'TOTAL', socChgTotal, socDchgTotal, socChgTotal + socDchgTotal);
    end
    
    % 사이클별 총합
    cycleChgTotal = 0;
    cycleDchgTotal = 0;
    for socIdx = 1:length(socLevels)
        socLevel = socLevels{socIdx};
        for dcIdx = 1:length(dcProfiles)
            dcProfile = dcProfiles{dcIdx};
            for chIdx = 1:length(channelNames)
                channelFieldName = channelNames{chIdx};
                if contains(channelFieldName, '_ChgEvent', 'IgnoreCase', true)
                    if isfield(eventData, channelFieldName)
                        chgData = eventData.(channelFieldName);
                        if isfield(chgData, socLevel) && isfield(chgData.(socLevel), dcProfile)
                            allFields = fieldnames(chgData.(socLevel).(dcProfile));
                            eventFields = allFields(startsWith(allFields, 'event'));
                            cycleChgTotal = cycleChgTotal + length(eventFields);
                        end
                    end
                elseif contains(channelFieldName, '_DchgEvent', 'IgnoreCase', true) || contains(channelFieldName, '_DchEvent', 'IgnoreCase', true)
                    if isfield(eventData, channelFieldName)
                        dchgData = eventData.(channelFieldName);
                        if isfield(dchgData, socLevel) && isfield(dchgData.(socLevel), dcProfile)
                            allFields = fieldnames(dchgData.(socLevel).(dcProfile));
                            eventFields = allFields(startsWith(allFields, 'event'));
                            cycleDchgTotal = cycleDchgTotal + length(eventFields);
                        end
                    end
                end
            end
        end
    end
    fprintf('\n');
    fprintf('========================================\n');
    fprintf('Cycle %s TOTAL: Charge=%d, Discharge=%d, Total=%d\n', ...
        cycleType, cycleChgTotal, cycleDchgTotal, cycleChgTotal + cycleDchgTotal);
    fprintf('========================================\n');
end

fprintf('\n=== Event Count Summary Complete ===\n');

%% Comprehensive Filtering Condition Table (All Channels)
% =========================================================================
fprintf('\n\n========================================\n');
fprintf('=== 필터링 조건표 (전체 채널) ===\n');
fprintf('========================================\n');

% 각 사이클별로 필터링 통계 집계
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    
    fprintf('\n');
    fprintf('========================================\n');
    fprintf('Cycle: %s\n', cycleType);
    fprintf('========================================\n');
    
    % filterStats에서 전체 채널 통계 집계
    if ~exist('filterStats', 'var') || isempty(filterStats)
        fprintf('  필터링 통계 데이터가 없습니다.\n');
        continue;
    end
    
    % 전체 전환 수 및 필터링 통계 집계 (전체 채널)
    totalTransitions_all = 0;
    totalRejectedShortDuration = 0;
    totalRejectedDirection = 0;
    totalRejectedTooShortWindow = 0;
    totalRejectedIStd = 0;
    totalRejectedCV = 0;
    totalRejectedPStd = 0;
    totalRejectedTooShortAfterCut = 0;
    
    % filterStats는 채널별, SOC별, DC별로 구조화되어 있음
    % 모든 채널에 대해 집계
    allChannelNames = fieldnames(filterStats);
    for chIdx = 1:length(allChannelNames)
        channelName = allChannelNames{chIdx};
        
        if ~isfield(filterStats, channelName)
            continue;
        end
        
        for socIdx = 1:length(socLevels)
            socLevel = socLevels{socIdx};
            
            if ~isfield(filterStats.(channelName), socLevel)
                continue;
            end
            
            for dcIdx = 1:length(dcProfiles)
                dcProfile = dcProfiles{dcIdx};
                
                if ~isfield(filterStats.(channelName).(socLevel), dcProfile)
                    continue;
                end
                
                stats = filterStats.(channelName).(socLevel).(dcProfile);
                
                % 전체 전환 수 집계
                totalTransitions_all = totalTransitions_all + stats.total_transitions;
                totalRejectedShortDuration = totalRejectedShortDuration + stats.rejected_short_duration;
                totalRejectedDirection = totalRejectedDirection + stats.rejected_direction;
                totalRejectedTooShortWindow = totalRejectedTooShortWindow + stats.rejected_too_short_window;
                totalRejectedIStd = totalRejectedIStd + stats.rejected_I_std;
                totalRejectedCV = totalRejectedCV + stats.rejected_CV;
                totalRejectedPStd = totalRejectedPStd + stats.rejected_P_std;
                totalRejectedTooShortAfterCut = totalRejectedTooShortAfterCut + stats.rejected_too_short_after_cut;
            end
        end
    end
    
    % 실제 최종 이벤트 수는 Event Count Summary에서 계산된 값 사용
    eventFile = fullfile(outputDir, sprintf('Lab_DC_DCIR_%s_Events.mat', cycleType));
    if exist(eventFile, 'file')
        load(eventFile);
        resultStructName = sprintf('Lab_DC_DCIR_%s', cycleType);
        eval(sprintf('eventData = %s;', resultStructName));
        channelNames = fieldnames(eventData);
        
        % 최종 검출된 충전/방전 이벤트 수 계산
        finalChgCount = 0;
        finalDchgCount = 0;
        for socIdx = 1:length(socLevels)
            socLevel = socLevels{socIdx};
            for dcIdx = 1:length(dcProfiles)
                dcProfile = dcProfiles{dcIdx};
                for chIdx = 1:length(channelNames)
                    channelFieldName = channelNames{chIdx};
                    if contains(channelFieldName, '_ChgEvent', 'IgnoreCase', true)
                        if isfield(eventData, channelFieldName)
                            chgData = eventData.(channelFieldName);
                            if isfield(chgData, socLevel) && isfield(chgData.(socLevel), dcProfile)
                                allFields = fieldnames(chgData.(socLevel).(dcProfile));
                                eventFields = allFields(startsWith(allFields, 'event'));
                                finalChgCount = finalChgCount + length(eventFields);
                            end
                        end
                    elseif contains(channelFieldName, '_DchgEvent', 'IgnoreCase', true) || contains(channelFieldName, '_DchEvent', 'IgnoreCase', true)
                        if isfield(eventData, channelFieldName)
                            dchgData = eventData.(channelFieldName);
                            if isfield(dchgData, socLevel) && isfield(dchgData.(socLevel), dcProfile)
                                allFields = fieldnames(dchgData.(socLevel).(dcProfile));
                                eventFields = allFields(startsWith(allFields, 'event'));
                                finalDchgCount = finalDchgCount + length(eventFields);
                            end
                        end
                    end
                end
            end
        end
    else
        finalChgCount = 0;
        finalDchgCount = 0;
    end
    
    if totalTransitions_all > 0
        % 충전 이벤트 필터링 조건표
        fprintf('\n--- 충전 이벤트 필터링 조건표 ---\n');
        fprintf('전체 전환 수: %d\n', totalTransitions_all);
        fprintf('최종 검출된 충전 이벤트 수: %d\n', finalChgCount);
        fprintf('\n');
        fprintf('%-5s | %-40s | %-50s | %12s | %12s | %8s\n', ...
            '단계', '필터링 조건', '설명', '제거 개수', '검출된 이벤트 수', '제거율');
        fprintf('-----|------------------------------------------|--------------------------------------------------|--------------|--------------|--------\n');
        
        % 각 단계에서 제거된 개수 계산 (누적)
        cumulative_rejected = 0;
        
        % Step 1: Short Duration
        step1_rejected = totalRejectedShortDuration;
        cumulative_rejected = cumulative_rejected + step1_rejected;
        step1_remaining = totalTransitions_all - cumulative_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            1, sprintf('최소 지속시간 (< %.1f초)', min_duration_chg), ...
            '이벤트 지속시간이 최소 기준보다 짧은 경우 제거', ...
            step1_rejected, '-', (step1_rejected/totalTransitions_all)*100);
        
        % Step 2: Direction
        step2_rejected = totalRejectedDirection;
        cumulative_rejected = cumulative_rejected + step2_rejected;
        step2_remaining = totalTransitions_all - cumulative_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            2, '충전 방향 확인 (< 50%% 양의 전류)', ...
            '양의 전류 비율이 50%% 미만인 경우 제거 (방전 이벤트)', ...
            step2_rejected, '-', (step2_rejected/totalTransitions_all)*100);
        
        % Step 3: Too Short Window
        step3_rejected = totalRejectedTooShortWindow;
        cumulative_rejected = cumulative_rejected + step3_rejected;
        step3_remaining = totalTransitions_all - cumulative_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            3, '안정 구간 길이 확인 (< 10.0초)', ...
            '안정 구간(3~10초) 길이가 부족한 경우 제거', ...
            step3_rejected, '-', (step3_rejected/totalTransitions_all)*100);
        
        % Step 4-6: I Std, CV, P Std (동시 확인)
        step4_rejected = totalRejectedIStd;
        step5_rejected = totalRejectedCV;
        step6_rejected = totalRejectedPStd;
        % 동시 확인되므로 중복 제거 개수는 최대값 사용
        max_stability_rejected = max([step4_rejected, step5_rejected, step6_rejected]);
        cumulative_rejected = cumulative_rejected + max_stability_rejected;
        step6_remaining = totalTransitions_all - cumulative_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            4, sprintf('전류 표준편차 (>= %.2f A)', max_I_std_chg), ...
            sprintf('안정 구간 전류 표준편차가 %.2f A 이상인 경우 제거', max_I_std_chg), ...
            step4_rejected, '-', (step4_rejected/totalTransitions_all)*100);
        
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            5, sprintf('변동계수 (>= %.1f%%)', max_CV_chg*100), ...
            sprintf('CV = std(I)/|mean(I)| >= %.1f%%인 경우 제거', max_CV_chg*100), ...
            step5_rejected, '-', (step5_rejected/totalTransitions_all)*100);
        
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            6, sprintf('전력 표준편차 (>= %.3f kW)', max_P_std_chg), ...
            sprintf('안정 구간 전력 표준편차가 %.3f kW 이상인 경우 제거', max_P_std_chg), ...
            step6_rejected, '-', (step6_rejected/totalTransitions_all)*100);
        
        fprintf('-----|------------------------------------------|--------------------------------------------------|--------------|--------------|--------\n');
        fprintf('  참고: 단계 4-6 (전류 표준편차, 변동계수, 전력 표준편차)는 동시에 확인됩니다.\n');
        fprintf('        각 제거 개수는 해당 조건을 만족하지 못한 이벤트 수를 나타냅니다.\n');
        fprintf('        하나의 이벤트가 여러 조건을 만족하지 못할 수 있어 제거 개수 합계가 중복될 수 있습니다.\n');
        fprintf('        중간 단계의 "검출된 이벤트 수"는 충전/방전 구분이 불가능하여 표시하지 않습니다.\n');
        fprintf('-----|------------------------------------------|--------------------------------------------------|--------------|--------------|--------\n');
        
        % Step 7: Too Short After Cut
        step7_rejected = totalRejectedTooShortAfterCut;
        cumulative_rejected = cumulative_rejected + step7_rejected;
        step7_remaining = totalTransitions_all - cumulative_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            7, '컷팅 후 최소 길이 확인 (< 5.0초)', ...
            '이벤트 시작 부분을 자른 후 길이가 5초 미만인 경우 제거', ...
            step7_rejected, '-', (step7_rejected/totalTransitions_all)*100);
        
        fprintf('-----|------------------------------------------|--------------------------------------------------|--------------|--------------|--------\n');
        fprintf('%-5s | %-40s | %-50s | %12d | %12d | %7.1f%%\n', ...
            '최종', '최종 검출된 충전 이벤트 수', ...
            '모든 필터링 조건을 통과한 최종 충전 이벤트 수', ...
            totalTransitions_all - finalChgCount, finalChgCount, ...
            ((totalTransitions_all - finalChgCount)/totalTransitions_all)*100);
        fprintf('========================================\n');
        
        % 방전 이벤트 필터링 조건표
        fprintf('\n--- 방전 이벤트 필터링 조건표 ---\n');
        fprintf('전체 전환 수: %d\n', totalTransitions_all);
        fprintf('최종 검출된 방전 이벤트 수: %d\n', finalDchgCount);
        fprintf('\n');
        fprintf('%-5s | %-40s | %-50s | %12s | %12s | %8s\n', ...
            '단계', '필터링 조건', '설명', '제거 개수', '검출된 이벤트 수', '제거율');
        fprintf('-----|------------------------------------------|--------------------------------------------------|--------------|--------------|--------\n');
        
        % 각 단계에서 제거된 개수 계산 (누적)
        cumulative_rejected = 0;
        
        % Step 1: Short Duration
        step1_rejected = totalRejectedShortDuration;
        cumulative_rejected = cumulative_rejected + step1_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            1, sprintf('최소 지속시간 (< %.1f초)', min_duration_dchg), ...
            '이벤트 지속시간이 최소 기준보다 짧은 경우 제거', ...
            step1_rejected, '-', (step1_rejected/totalTransitions_all)*100);
        
        % Step 2: Direction
        step2_rejected = totalRejectedDirection;
        cumulative_rejected = cumulative_rejected + step2_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            2, '방전 방향 확인 (모두 음의 전류)', ...
            '모든 전류가 음수가 아닌 경우 제거 (충전 이벤트)', ...
            step2_rejected, '-', (step2_rejected/totalTransitions_all)*100);
        
        % Step 3: Too Short Window
        step3_rejected = totalRejectedTooShortWindow;
        cumulative_rejected = cumulative_rejected + step3_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            3, '안정 구간 길이 확인 (< 3초)', ...
            '안정 구간(이벤트 끝까지) 길이가 3초 미만인 경우 제거', ...
            step3_rejected, '-', (step3_rejected/totalTransitions_all)*100);
        
        % Step 4-6: I Std, CV, P Std (동시 확인)
        step4_rejected = totalRejectedIStd;
        step5_rejected = totalRejectedCV;
        step6_rejected = totalRejectedPStd;
        % 동시 확인되므로 중복 제거 개수는 최대값 사용
        max_stability_rejected = max([step4_rejected, step5_rejected, step6_rejected]);
        cumulative_rejected = cumulative_rejected + max_stability_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            4, sprintf('전류 표준편차 (>= %.2f A)', max_I_std_dchg), ...
            sprintf('안정 구간 전류 표준편차가 %.2f A 이상인 경우 제거', max_I_std_dchg), ...
            step4_rejected, '-', (step4_rejected/totalTransitions_all)*100);
        
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            5, sprintf('변동계수 (>= %.1f%%)', max_CV_dchg*100), ...
            sprintf('CV = std(I)/|mean(I)| >= %.1f%%인 경우 제거', max_CV_dchg*100), ...
            step5_rejected, '-', (step5_rejected/totalTransitions_all)*100);
        
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            6, sprintf('전력 표준편차 (>= %.3f kW)', max_P_std_dchg), ...
            sprintf('안정 구간 전력 표준편차가 %.3f kW 이상인 경우 제거', max_P_std_dchg), ...
            step6_rejected, '-', (step6_rejected/totalTransitions_all)*100);
        
        fprintf('-----|------------------------------------------|--------------------------------------------------|--------------|--------------|--------\n');
        fprintf('  참고: 단계 4-6 (전류 표준편차, 변동계수, 전력 표준편차)는 동시에 확인됩니다.\n');
        fprintf('        각 제거 개수는 해당 조건을 만족하지 못한 이벤트 수를 나타냅니다.\n');
        fprintf('        하나의 이벤트가 여러 조건을 만족하지 못할 수 있어 제거 개수 합계가 중복될 수 있습니다.\n');
        fprintf('        중간 단계의 "검출된 이벤트 수"는 충전/방전 구분이 불가능하여 표시하지 않습니다.\n');
        fprintf('-----|------------------------------------------|--------------------------------------------------|--------------|--------------|--------\n');
        
        % Step 7: Too Short After Cut
        step7_rejected = totalRejectedTooShortAfterCut;
        cumulative_rejected = cumulative_rejected + step7_rejected;
        fprintf('%-5d | %-40s | %-50s | %12d | %12s | %7.1f%%\n', ...
            7, '컷팅 후 최소 길이 확인 (< 5.0초)', ...
            '이벤트 시작 부분을 자른 후 길이가 5초 미만인 경우 제거', ...
            step7_rejected, '-', (step7_rejected/totalTransitions_all)*100);
        
        fprintf('-----|------------------------------------------|--------------------------------------------------|--------------|--------------|--------\n');
        fprintf('%-5s | %-40s | %-50s | %12d | %12d | %7.1f%%\n', ...
            '최종', '최종 검출된 방전 이벤트 수', ...
            '모든 필터링 조건을 통과한 최종 방전 이벤트 수', ...
            totalTransitions_all - finalDchgCount, finalDchgCount, ...
            ((totalTransitions_all - finalDchgCount)/totalTransitions_all)*100);
        fprintf('========================================\n');
    end
end

%% Visualization: SOC별 충전/방전 이벤트 전류 프로파일 (SOC90, 70, 50)
% =========================================================================
fprintf('\n\n========================================\n');
fprintf('=== Creating SOC별 Charge/Discharge Current Profiles ===\n');
fprintf('========================================\n');

% figures 폴더 미리 생성
figuresBaseDir = fullfile(outputDir, 'figures');
if ~exist(figuresBaseDir, 'dir')
    mkdir(figuresBaseDir);
    fprintf('Created base figures directory: %s\n', figuresBaseDir);
end

% SOC 레벨 리스트
targetSOCLevels = {'SOC90', 'SOC70', 'SOC50'};
eventTypes = {'Charge', 'Discharge'};

% 0cyc 데이터 로드
if ismember('0cyc', cycleTypes)
    eventFile = fullfile(outputDir, 'Lab_DC_DCIR_0cyc_Events.mat');
    fprintf('Looking for event file: %s\n', eventFile);
    if exist(eventFile, 'file')
        fprintf('Loading event file...\n');
        load(eventFile);
        resultStructName = 'Lab_DC_DCIR_0cyc';
        eval(sprintf('eventData = %s;', resultStructName));
        fprintf('Event data loaded successfully.\n');
        
        % 모든 채널 찾기
        channelNames = fieldnames(eventData);
        fprintf('Available channels: %s\n', strjoin(channelNames, ', '));
        
        % 각 SOC 레벨과 이벤트 타입별로 figure 생성
        for socIdx = 1:length(targetSOCLevels)
            socLevel = targetSOCLevels{socIdx};
            
            for eventTypeIdx = 1:length(eventTypes)
                eventType = eventTypes{eventTypeIdx};
                
                % 이벤트 타입에 맞는 구조체 이름 찾기
                if strcmpi(eventType, 'Charge')
                    eventStructSuffix = 'ChgEvent';
                else
                    eventStructSuffix = 'DchgEvent';
                end
                
                % Figure 생성
                fig = figure('Name', sprintf('0cyc %s %s Events - Current vs Time', socLevel, eventType), ...
                    'Position', [100, 100, 1400, 800], 'Visible', 'on');
                
                % DC 프로파일 리스트 (DC1~DC8)
                dcProfiles = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};
                numProfiles = length(dcProfiles);
                cols = 4;
                rows = ceil(numProfiles / cols);
                
                % 모든 채널에서 데이터 수집
                allEventsFound = false;
                
                for p = 1:numProfiles
                    profileName = dcProfiles{p};
                    subplot(rows, cols, p);
                    
                    eventCount = 0;
                    hold on;
                    
                    % 모든 채널에서 이벤트 찾기
                    for chIdx = 1:length(channelNames)
                        channelFieldName = channelNames{chIdx};
                        
                        % 채널 필드 이름이 이벤트 타입과 일치하는지 확인
                        % channelFieldName은 이미 ch9_Drive_0cyc_ChgEvent 형식
                        % Note: Support both _DchEvent (old) and _DchgEvent (new) for backward compatibility
                        if strcmpi(eventType, 'Charge')
                            if ~contains(channelFieldName, '_ChgEvent', 'IgnoreCase', true)
                                continue;  % 충전 이벤트가 아니면 스킵
                            end
                        else  % Discharge
                            if ~contains(channelFieldName, '_DchgEvent', 'IgnoreCase', true) && ~contains(channelFieldName, '_DchEvent', 'IgnoreCase', true)
                                continue;  % 방전 이벤트가 아니면 스킵
                            end
                        end
                        
                        % 이벤트 데이터 확인
                        if isfield(eventData, channelFieldName)
                            chgData = eventData.(channelFieldName);
                            
                            if isfield(chgData, socLevel) && isfield(chgData.(socLevel), profileName)
                                % 이벤트 필드 찾기
                                allFields = fieldnames(chgData.(socLevel).(profileName));
                                eventFields = allFields(startsWith(allFields, 'event'));
                                
                                % 여러 이벤트를 위한 색상 팔레트
                                colorPalette = lines(length(eventFields));
                                if isempty(colorPalette)
                                    colorPalette = [0 0 1];  % 기본 파란색
                                end
                                
                                for e = 1:length(eventFields)
                                    evtName = eventFields{e};
                                    
                                    % Use full event data for visualization if available
                                    if isfield(chgData.(socLevel).(profileName).(evtName), 't_full') && ...
                                       isfield(chgData.(socLevel).(profileName).(evtName), 'I_full')
                                        t_evt = chgData.(socLevel).(profileName).(evtName).t_full;
                                        I_evt = chgData.(socLevel).(profileName).(evtName).I_full;
                                    elseif isfield(chgData.(socLevel).(profileName).(evtName), 't') && ...
                                           isfield(chgData.(socLevel).(profileName).(evtName), 'I')
                                        t_evt = chgData.(socLevel).(profileName).(evtName).t;
                                        I_evt = chgData.(socLevel).(profileName).(evtName).I;
                                    else
                                        continue;
                                    end
                                    
                                    if ~isempty(t_evt) && ~isempty(I_evt) && length(t_evt) == length(I_evt)
                                        % Check if data is valid (not all NaN or Inf)
                                        if any(~isnan(t_evt)) && any(~isnan(I_evt)) && any(~isinf(t_evt)) && any(~isinf(I_evt))
                                            % 시간을 0부터 시작하도록 정규화
                                            t_norm = t_evt - t_evt(1);
                                            
                                            % 각 이벤트마다 다른 색상 사용
                                            if size(colorPalette, 1) >= e
                                                plot(t_norm, I_evt, '-', 'LineWidth', 1.5, 'Color', colorPalette(e, :));
                                            else
                                                plot(t_norm, I_evt, '-', 'LineWidth', 1.5);
                                            end
                                            
                                            % Mark stable region end point with a dot
                                            % First try to use index if available (more accurate)
                                            if isfield(chgData.(socLevel).(profileName).(evtName), 'stable_region_end_idx')
                                                stable_end_idx = chgData.(socLevel).(profileName).(evtName).stable_region_end_idx;
                                                if stable_end_idx > 0 && stable_end_idx <= length(I_evt) && stable_end_idx <= length(t_norm)
                                                    plot(t_norm(stable_end_idx), I_evt(stable_end_idx), 'ro', ...
                                                        'MarkerSize', 3, 'MarkerFaceColor', 'r', 'LineWidth', 2);
                                                end
                                            % Fallback to time-based search if index not available
                                            elseif isfield(chgData.(socLevel).(profileName).(evtName), 'stable_region_end_time')
                                                stable_end_t = chgData.(socLevel).(profileName).(evtName).stable_region_end_time;
                                                if stable_end_t > 0 && stable_end_t <= t_norm(end)
                                                    % Find index closest to stable_end_t in normalized time
                                                    [~, stable_end_idx] = min(abs(t_norm - stable_end_t));
                                                    if stable_end_idx > 0 && stable_end_idx <= length(I_evt)
                                                        plot(t_norm(stable_end_idx), I_evt(stable_end_idx), 'ro', ...
                                                            'MarkerSize', 3, 'MarkerFaceColor', 'r', 'LineWidth', 2);
                                                    end
                                                end
                                            end
                                            
                                            eventCount = eventCount + 1;
                                            allEventsFound = true;
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    % 이벤트가 없으면 "no data" 표시
                    if eventCount == 0
                        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'middle', 'FontSize', 14, ...
                            'Color', 'red', 'Units', 'normalized');
                    end
                    
                    xlabel('Time (s)');
                    ylabel('Current (A)');
                    title(sprintf('%s (%s, %d events)', profileName, eventType, eventCount));
                    grid on;
                    hold off;
                end
                
                sgtitle(sprintf('0cyc %s - %s Events (Current vs Time, All Channels)', socLevel, eventType), ...
                    'FontSize', 14, 'FontWeight', 'bold');
                
                % 이벤트가 하나라도 있으면 저장, 없으면 스킵
                if allEventsFound
                    % 저장
                    figDir = fullfile(outputDir, 'figures');
                    savePath = fullfile(figDir, sprintf('0cyc_%s_%s_CurrentTime.fig', socLevel, eventType));
                    fprintf('Saving %s %s figure to: %s\n', socLevel, eventType, savePath);
                    try
                        saveas(fig, savePath);
                        fprintf('✓ Successfully saved: %s\n', savePath);
                        set(fig, 'Visible', 'on');
                    catch ME
                        fprintf('✗ Error saving figure: %s\n', ME.message);
                        fprintf('  Error details: %s\n', getReport(ME));
                    end
                else
                    % 이벤트가 없으면 figure 닫기
                    fprintf('No %s events found for %s. Skipping figure creation.\n', eventType, socLevel);
                    close(fig);
                end
            end
        end
    else
        fprintf('0cyc event file not found: %s\n', eventFile);
    end
else
    fprintf('0cyc not in selected cycles. Skipping SOC visualization.\n');
end
% 
% DriveCycle_DataAggregation_02
% DriveCycle_Visualization_03
% DriveCycle_CorrelationAnalysis_04


