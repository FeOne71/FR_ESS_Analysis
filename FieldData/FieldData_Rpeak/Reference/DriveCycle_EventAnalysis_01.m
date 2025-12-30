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
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

selectedChannels = [];  % Empty array: analyze all channels (Ch9~Ch16)
selectedCycles = {'0cyc','200cyc','400cyc','600cyc','800cyc'};  
selectedSOC = {'SOC50'};  

% 디버그 모드
debug_mode = false;  % true: 필터링된 이벤트의 전류 그래프를 저장 (느림, 필요시만 true로 설정)
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

%% Settings (64Ah cell)
dt_list = [1, 3, 5, 10, 30, 60];
Cnom = 64;    % Battery capacity [Ah]
Pnom = 0.235; % Nominal Power    [kW] 3.68V * 64Ah = 235.52W
current_threshold = Cnom * 0.02;  % Idle current threshold [A]
min_duration      = 60;           % Minimum driving duration [s]
max_P_std         = Pnom * 0.05;  % Maximum power standard deviation [kW] (cell level)
max_I_std         = Cnom * 0.05;  % Maximum current standard deviation [A] (2%)
max_I_range       = 7.0;          % Maximum current range (max - min) [A] in 2~60s window
max_CV            = 0.50;          % Maximum coefficient of variation (std/mean) in 2~60s window (10%)

fprintf('\n=== Analysis Parameters ===\n');
fprintf('Current threshold: %.2f A\n', current_threshold);
fprintf('Min duration: %d s\n', min_duration);
fprintf('Max P std: %.3f kW\n', max_P_std);
fprintf('Max I std: %.2f A (%.1f%%)\n', max_I_std, (max_I_std/Cnom)*100);
fprintf('Max I range: %.2f A\n', max_I_range);
fprintf('Max CV: %.1f%%\n', max_CV*100);
fprintf('Debug mode: %s\n', mat2str(debug_mode));

% Create debug folder if debug mode is enabled
if debug_mode
    if ~exist('figures', 'dir')
        mkdir('figures');
    end
    if ~exist('figures/debug', 'dir')
        mkdir('figures/debug');
    end
    fprintf('Debug plots will be saved to: figures/debug/\n');
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
        dchg_struct_name = [channelName '_DchEvent'];
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
                
                %% Step 1: Idle -> Load transition detection
                is_idle = abs(I) < current_threshold;
                is_driving = abs(I) >= current_threshold;
                
                idle_to_driving = find(is_idle(1:end-1) & is_driving(2:end));
                
                if isempty(idle_to_driving)
                    fprintf('    %s: No transitions detected → SKIP\n', profile_name);
                    continue;  % Skip if no transition
                end
                
                %% Step 2: Event Detection (charging/discharging split)
                chg_event_count = 0;
                dchg_event_count = 0;
                for i = 1:length(idle_to_driving)
                    idx1 = idle_to_driving(i);
                    start_driving_idx = idx1 + 1;
                    
                    % Determine event type from current direction at start
                    event_type = sign(I(start_driving_idx));  % 1: charge, -1: discharge
                    if event_type == 0
                        continue;  % Skip if current is exactly zero
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
                    
                    if driving_time < min_duration
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
                    
                    % Charging event: all currents must be positive
                    % Discharging event: all currents must be negative
                    if event_type > 0
                        % Charging event: all currents must be > 0
                        if any(I_driving <= 0)
                            % Current is not positive in charging event, skip
                            if debug_mode
                                figure('Visible', 'off');
                                plot(t_seg - t_seg(1), I_seg, 'b-', 'LineWidth', 1.5);
                                hold on;
                                plot(t_seg(2:end) - t_seg(1), I_driving, 'r-', 'LineWidth', 2);
                                xlabel('Time (s)');
                                ylabel('Current (A)');
                                title(sprintf('Rejected: Charging Event Has Non-Positive Current (Min=%.2fA)', min(I_driving)));
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
                    
                    % idx2: end of driving segment
                    idx2 = driving_end_idx - start_idx + 1;
                    
                    %% Step 3: Stability check
                    % Use time-based indexing: fixed 2~60 second window
                    % This ensures consistent time window regardless of event duration
                    t_start_stable = t_seg(1) + 2;   % 2 seconds after event start
                    t_end_stable = t_seg(1) + 60;    % 60 seconds after event start (fixed)
                    
                    % Check if event is long enough for 2~60s window
                    if t_seg(end) < t_end_stable
                        % Event is shorter than 60 seconds, skip
                        if debug_mode
                            t_seg_temp = t_seconds(start_idx:end_idx);
                            I_seg_temp = I(start_idx:end_idx);
                            figure('Visible', 'off');
                            plot(t_seg_temp - t_seg_temp(1), I_seg_temp, 'b-', 'LineWidth', 1.5);
                            xlabel('Time (s)');
                            ylabel('Current (A)');
                            title(sprintf('Rejected: Too Short for Stable Window (Duration=%.2fs, Need>=60s)', driving_time));
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
                    
                    % Find indices corresponding to stable time window (2~60 seconds)
                    stable_mask = (t_seg >= t_start_stable) & (t_seg <= t_end_stable);
                    stable_indices = find(stable_mask);
                    
                    if length(stable_indices) >= 3  % Need at least 3 points for std
                        % Calculate statistics for 2~60s window
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
                        % This should not happen if event is >= 60 seconds
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
                    
                    % Check multiple criteria: power std, current std, current range, and CV
                    reject_reasons = {};
                    if power_std >= max_P_std
                        reject_reasons{end+1} = sprintf('High P Std (%.3f >= %.3f)', power_std, max_P_std);
                    end
                    if current_std >= max_I_std
                        reject_reasons{end+1} = sprintf('High I Std (%.3f >= %.3f)', current_std, max_I_std);
                    end
                    if current_range >= max_I_range
                        reject_reasons{end+1} = sprintf('High I Range (%.3f >= %.3f)', current_range, max_I_range);
                    end
                    if current_CV >= max_CV
                        reject_reasons{end+1} = sprintf('High CV (%.3f >= %.3f)', current_CV, max_CV);
                    end
                    
                    if ~isempty(reject_reasons)
                        % Debug plotting: Rejected due to filtering criteria
                        if debug_mode
                            figure('Visible', 'off');
                            plot(t_seg - t_seg(1), I_seg, 'b-', 'LineWidth', 1.5);
                            hold on;
                            % Mark stable window (2~60s)
                            stable_t = t_seg(stable_indices) - t_seg(1);
                            stable_I = I_seg(stable_indices);
                            plot(stable_t, stable_I, 'r-', 'LineWidth', 2);
                            xlabel('Time (s)');
                            ylabel('Current (A)');
                            title_str = sprintf('Rejected: %s', strjoin(reject_reasons, ', '));
                            title(title_str, 'FontSize', 10);
                            legend('Full Event', 'Stable Window (2-60s)', 'Location', 'best');
                            grid on;
                            % Add statistics text
                            stats_text = sprintf('I_std=%.3f (limit=%.3f)\nI_range=%.3f (limit=%.3f)\nCV=%.3f (limit=%.3f)', ...
                                current_std, max_I_std, current_range, max_I_range, current_CV, max_CV);
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
                    target_struct.(soc_level).(profile_name).(evtName).t = t_seg;
                    target_struct.(soc_level).(profile_name).(evtName).I = I_seg;
                    target_struct.(soc_level).(profile_name).(evtName).V = V_seg;
                    target_struct.(soc_level).(profile_name).(evtName).P = P_seg;
                    target_struct.(soc_level).(profile_name).(evtName).I_std = current_std;      % Current std in 2~60s window [A]
                    target_struct.(soc_level).(profile_name).(evtName).I_mean = current_mean;     % Average current in 2~60s window [A]
                    target_struct.(soc_level).(profile_name).(evtName).I_range = current_range;  % Current range (max-min) in 2~60s window [A]
                    target_struct.(soc_level).(profile_name).(evtName).I_CV = current_CV;        % Coefficient of variation (std/mean) in 2~60s window
                    target_struct.(soc_level).(profile_name).(evtName).P_std = power_std;
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
end

fprintf('\n\n========================================\n');
fprintf('=== Event Analysis Complete ===\n');
fprintf('========================================\n');
fprintf('All results saved to: %s\n', outputDir);
fprintf('Event files: Lab_DC_DCIR_*cyc_Events.mat\n');
if debug_mode
    fprintf('Debug plots saved to: %s\n', fullfile(outputDir, 'figures', 'debug'));
end

%% Visualization: 0cyc DC별 전류-시간 그래프
% =========================================================================
fprintf('\n\n========================================\n');
fprintf('=== Creating 0cyc DC별 Current-Time Plots ===\n');
fprintf('========================================\n');

% figures 폴더 미리 생성
figuresBaseDir = fullfile(outputDir, 'figures');
if ~exist(figuresBaseDir, 'dir')
    mkdir(figuresBaseDir);
    fprintf('Created base figures directory: %s\n', figuresBaseDir);
end

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
        
        % 9번 채널의 충전 이벤트 찾기
        channelNames = fieldnames(eventData);
        fprintf('Available channels: %s\n', strjoin(channelNames, ', '));
        
        % 9번 채널의 충전 이벤트 필드 찾기 (ch9_Drive_0cyc_ChgEvent 형식)
        ch9ChgEventName = '';
        for i = 1:length(channelNames)
            if contains(channelNames{i}, 'ch9', 'IgnoreCase', true) && contains(channelNames{i}, 'ChgEvent')
                ch9ChgEventName = channelNames{i};
                break;
            end
        end
        
        if isempty(ch9ChgEventName)
            fprintf('ERROR: ch9 Charge Event not found. Available fields: %s\n', strjoin(channelNames, ', '));
        else
            fprintf('Using channel: %s\n', ch9ChgEventName);
            chgData = eventData.(ch9ChgEventName);
            
            % SOC 레벨 확인 (첫 번째 SOC 사용)
            socLevels = fieldnames(chgData);
            if ~isempty(socLevels)
                socLevel = socLevels{1};
                fprintf('Using SOC level: %s\n', socLevel);
                
                % figures 폴더 확인 (이미 생성되어 있음)
                figDir = fullfile(outputDir, 'figures');
                fprintf('Using figures directory: %s\n', figDir);
                
                % Charge 이벤트 그래프 생성 (무조건 생성)
                fprintf('\n--- Processing Charge Events (Ch9 only) ---\n');
                fig_charge = figure('Name', sprintf('0cyc Ch9 Charge Events - Current vs Time'), ...
                    'Position', [100, 100, 1400, 800], 'Visible', 'on');
                
                if isfield(chgData, socLevel)
                    fprintf('Charge data found for %s\n', socLevel);
                    chgProfiles = fieldnames(chgData.(socLevel));
                    fprintf('Found %d charge profiles: %s\n', length(chgProfiles), strjoin(chgProfiles, ', '));
                else
                    fprintf('Charge data not found for %s - creating empty plot\n', socLevel);
                    chgProfiles = {};
                end
                
                % DC 프로파일이 없으면 기본 프로파일 리스트 사용 (DC1~DC8)
                if isempty(chgProfiles)
                    chgProfiles = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};
                end
                
                numProfiles = length(chgProfiles);
                cols = 4;
                rows = ceil(numProfiles / cols);
                
                for p = 1:numProfiles
                    profileName = chgProfiles{p};
                    subplot(rows, cols, p);
                    
                    eventCount = 0;
                    if isfield(chgData, socLevel) && isfield(chgData.(socLevel), profileName)
                        % 이벤트 필드 찾기
                        allFields = fieldnames(chgData.(socLevel).(profileName));
                        eventFields = allFields(startsWith(allFields, 'event'));
                        fprintf('  Profile %s: Found %d event fields\n', profileName, length(eventFields));
                        
                        hold on;
                        % 여러 이벤트를 위한 색상 팔레트 (파란색 계열)
                        colorPalette = lines(length(eventFields));  % 각 이벤트마다 다른 색상
                        if isempty(colorPalette)
                            colorPalette = [0 0 1];  % 기본 파란색
                        end
                        
                        for e = 1:length(eventFields)
                            evtName = eventFields{e};
                            if isfield(chgData.(socLevel).(profileName).(evtName), 't') && ...
                               isfield(chgData.(socLevel).(profileName).(evtName), 'I')
                                t_evt = chgData.(socLevel).(profileName).(evtName).t;
                                I_evt = chgData.(socLevel).(profileName).(evtName).I;
                                if ~isempty(t_evt) && ~isempty(I_evt)
                                    % 시간을 0부터 시작하도록 정규화
                                    t_norm = t_evt - t_evt(1);
                                    % 각 이벤트마다 다른 색상 사용
                                    if size(colorPalette, 1) >= e
                                        plot(t_norm, I_evt, '-', 'LineWidth', 1.5, 'Color', colorPalette(e, :));
                                    else
                                        plot(t_norm, I_evt, 'b-', 'LineWidth', 1.5);
                                    end
                                    eventCount = eventCount + 1;
                                end
                            end
                        end
                    else
                        fprintf('  Profile %s: Not found in data\n', profileName);
                    end
                    
                    % 이벤트가 없으면 "no data" 표시
                    if eventCount == 0
                        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'middle', 'FontSize', 14, ...
                            'Color', 'red', 'Units', 'normalized');
                    end
                    
                    fprintf('  Profile %s: Plotted %d events\n', profileName, eventCount);
                    xlabel('Time (s)');
                    ylabel('Current (A)');
                    title(sprintf('%s (Charge, %d events)', profileName, eventCount));
                    grid on;
                    hold off;
                end
                
                sgtitle(sprintf('0cyc Ch9 %s - Charge Events (Current vs Time)', socLevel), ...
                    'FontSize', 14, 'FontWeight', 'bold');
                
                % 저장 (figures 폴더에 직접 저장) - 무조건 저장
                savePath_charge = fullfile(figDir, sprintf('0cyc_Ch9_%s_Charge_CurrentTime.fig', socLevel));
                fprintf('Saving charge figure to: %s\n', savePath_charge);
                try
                    saveas(fig_charge, savePath_charge);
                    fprintf('✓ Successfully saved: %s\n', savePath_charge);
                    set(fig_charge, 'Visible', 'on');
                catch ME
                    fprintf('✗ Error saving charge figure: %s\n', ME.message);
                    fprintf('  Error details: %s\n', getReport(ME));
                end
            else
                fprintf('No SOC levels found in charge data.\n');
            end
        end
    else
        fprintf('0cyc event file not found: %s\n', eventFile);
    end
else
    fprintf('0cyc not in selected cycles. Skipping visualization.\n');
end

fprintf('\n=== 0cyc Visualization Complete ===\n');

DriveCycle_DataAggregation_02
DriveCycle_Visualization_03
DriveCycle_CorrelationAnalysis_04


