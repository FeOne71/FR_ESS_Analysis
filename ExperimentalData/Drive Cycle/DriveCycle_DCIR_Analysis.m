%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab Data - Drive Cycle DCIR Analysis (RPT0)
% Event detection and time-based DCIR calculation from real load profile data
% Logic adapted from FieldData_DCIR_Charge_CurrentClustering_Auto.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Configuration - User Settings
% =========================================================================
% 데이터 경로 설정
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
capacityDataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Capacity_Trend_Figures\Capacity_Data_Static.mat';

selectedChannels = [];  % Empty array: analyze all channels (Ch9~Ch16)
selectedCycles = {'0cyc','200cyc','400cyc','600cyc','800cyc'};  
selectedSOC = {'SOC50'};  

% 상관분석 옵션
correlationByDC = true;  % true: DC별로 상관분석, false: 모든 DC 평균하여 상관분석

% 디버그 모드
debug_mode = false;  % true: 필터링된 이벤트의 전류 그래프를 저장 (느림, 필요시만 true로 설정)
% =========================================================================

fprintf('=== Lab Drive Cycle DCIR Analysis ===\n');
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

%% Load Capacity Data
fprintf('\n=== Loading Capacity Data ===\n');
if exist(capacityDataFile, 'file')
    fprintf('Loading: %s\n', capacityDataFile);
    load(capacityDataFile);
    if exist('allChannelsCapacity', 'var')
        fprintf('Capacity data loaded successfully\n');
        capacityDataAvailable = true;
    else
        fprintf('WARNING: allChannelsCapacity variable not found in file\n');
        capacityDataAvailable = false;
    end
else
    fprintf('WARNING: Capacity data file not found: %s\n', capacityDataFile);
    fprintf('Correlation analysis will be skipped\n');
    capacityDataAvailable = false;
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
min_duration      = 10;           % Minimum driving duration [s]
max_P_std         = Pnom * 0.03;  % Maximum power standard deviation [kW] (cell level)
max_I_std         = Cnom * 0.025; % Maximum current standard deviation [A]

fprintf('\n=== Analysis Parameters ===\n');
fprintf('Current threshold: %.2f A\n', current_threshold);
fprintf('Min duration: %d s\n', min_duration);
fprintf('Max P std: %.3f kW\n', max_P_std);
fprintf('Max I std: %.2f A\n', max_I_std);

%% Process each cycle
allCycleResults = struct();

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
        else
            fprintf('Processing all SOC levels: %s\n', strjoin(soc_levels, ', '));
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
                
                % Find end index of driving
                driving_end_idx = start_driving_idx;
                while driving_end_idx <= length(is_driving) && is_driving(driving_end_idx)
                    driving_end_idx = driving_end_idx + 1;
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
                        if ~exist('figures/debug', 'dir')
                            mkdir('figures/debug');
                        end
                        filename = sprintf('figures/debug/%s_%s_%s_event%d_ShortDuration.fig', ...
                            channelName, soc_level, profile_name, i);
                        saveas(gcf, filename);
                        close(gcf);
                    end
                    continue;
                end
                
                % Extract segment
                t_seg = t_seconds(start_idx:end_idx);
                I_seg = I(start_idx:end_idx);
                V_seg = V(start_idx:end_idx);
                
                P_seg = V_seg .* I_seg / 1000;  % [kW]
                
                % idx2: end of driving segment
                idx2 = driving_end_idx - start_idx + 1;
                
                %% Step 3: Stability check
                power_std   = std(P_seg(3:idx2-2));
                current_std = std(I_seg(3:idx2-2));
                
                if power_std >= max_P_std || current_std >= max_I_std
                    % Debug plotting: Rejected due to high std
                    if debug_mode
                        figure('Visible', 'off');
                        plot(t_seg - t_seg(1), I_seg, 'b-', 'LineWidth', 1.5);
                        xlabel('Time (s)');
                        ylabel('Current (A)');
                        if power_std >= max_P_std && current_std >= max_I_std
                            title(sprintf('Rejected: High Std (P_std=%.3f, I_std=%.3f)', power_std, current_std));
                        elseif power_std >= max_P_std
                            title(sprintf('Rejected: High P Std (Val=%.3f, Limit=%.3f)', power_std, max_P_std));
                        else
                            title(sprintf('Rejected: High I Std (Val=%.3f, Limit=%.3f)', current_std, max_I_std));
                        end
                        grid on;
                        % Ensure debug folder exists
                        if ~exist('figures/debug', 'dir')
                            mkdir('figures/debug');
                        end
                        filename = sprintf('figures/debug/%s_%s_%s_event%d_HighStd.fig', ...
                            channelName, soc_level, profile_name, i);
                        saveas(gcf, filename);
                        close(gcf);
                    end
                    continue;
                end                   
                
                % Charging/discharging split
                mean_current = mean(I_seg);
                if mean_current > 0
                    chg_event_count = chg_event_count + 1;
                    evtName = sprintf('event%d', chg_event_count);
                    eval(sprintf('target_struct = %s.(chg_struct_name);', resultStructName));
                elseif mean_current < 0
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
                target_struct.(soc_level).(profile_name).(evtName).I_std = current_std;
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
    eval(sprintf('save(''Lab_DC_DCIR_%s_Events.mat'', ''%s'');', cycleType, resultStructName));
    eval(sprintf('cycleResult = %s;', resultStructName));
    % Convert cycle type to valid field name (add prefix if starts with number)
    if ~isempty(regexp(cycleType, '^\d', 'once'))
        validFieldName = ['cycle_' cycleType];
    else
        validFieldName = cycleType;
    end
    allCycleResults.(validFieldName) = cycleResult;
    fprintf('\n=== %s cycle processing complete ===\n', cycleType);
end

%% Summary
fprintf('\n\n========================================\n');
fprintf('=== Analysis Summary ===\n');
fprintf('========================================\n');

for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    resultStructName = sprintf('Lab_DC_DCIR_%s', cycleType);
    
    if ~exist(resultStructName, 'var')
        continue;
    end
    
    fprintf('\n--- %s Cycle ---\n', cycleType);
    eval(sprintf('resultData = %s;', resultStructName));
    
    allFields = fieldnames(resultData);
    fprintf('Available fields in result: %s\n', strjoin(allFields, ', '));
total_events = 0;
    
    % Extract unique channels from field names
    % Fields are already named as ch16_Drive_0cyc_ChgEvent and ch16_Drive_0cyc_DchEvent
    channels = {};
    chg_structs = {};
    dchg_structs = {};
    
    for i = 1:length(allFields)
        fieldName = allFields{i};
        % Extract channel name (remove _ChgEvent or _DchEvent suffix)
        if contains(fieldName, '_ChgEvent')
            channelName = strrep(fieldName, '_ChgEvent', '');
            if ~ismember(channelName, channels)
                channels{end+1} = channelName;
                chg_structs{end+1} = fieldName;  % Store the actual field name
            end
        elseif contains(fieldName, '_DchEvent')
            channelName = strrep(fieldName, '_DchEvent', '');
            if ~ismember(channelName, channels)
                % Check if channel already exists
                chIdx = find(strcmp(channels, channelName), 1);
                if isempty(chIdx)
                    channels{end+1} = channelName;
                    chg_structs{end+1} = '';
                    dchg_structs{end+1} = fieldName;
                else
                    dchg_structs{chIdx} = fieldName;
                end
            end
        end
    end
    
for ch_idx = 1:length(channels)
    channelName = channels{ch_idx};
        
        % Get actual struct names
        if ch_idx <= length(chg_structs)
            chg_struct_name = chg_structs{ch_idx};
        else
    chg_struct_name = [channelName '_ChgEvent'];
        end
        
        if ch_idx <= length(dchg_structs) && ~isempty(dchg_structs{ch_idx})
            dchg_struct_name = dchg_structs{ch_idx};
        else
    dchg_struct_name = [channelName '_DchEvent'];
        end
    
    chg_events = 0;
    dchg_events = 0;
        chg_dc_counts = containers.Map();  % DC별 충전 이벤트 개수
        dchg_dc_counts = containers.Map(); % DC별 방전 이벤트 개수
        
        if ~isempty(chg_struct_name) && isfield(resultData, chg_struct_name)
            chg_data = resultData.(chg_struct_name);
            all_soc_levels = fieldnames(chg_data);
            
            % Filter SOC levels based on user selection
            if isempty(selectedSOC)
                soc_levels = all_soc_levels;
            else
                soc_levels = intersect(selectedSOC, all_soc_levels);
            end
            
            for soc_idx = 1:length(soc_levels)
                soc_level = soc_levels{soc_idx};
                if isfield(chg_data, soc_level)
                    profiles = fieldnames(chg_data.(soc_level));
                    for prof_idx = 1:length(profiles)
                        profile_name = profiles{prof_idx};
                        if isfield(chg_data.(soc_level), profile_name)
                            events = fieldnames(chg_data.(soc_level).(profile_name));
                            event_count = length(events);
                            chg_events = chg_events + event_count;
                            if isKey(chg_dc_counts, profile_name)
                                chg_dc_counts(profile_name) = chg_dc_counts(profile_name) + event_count;
                            else
                                chg_dc_counts(profile_name) = event_count;
                            end
                        end
                    end
            end
        end
    end
    
        if ~isempty(dchg_struct_name) && isfield(resultData, dchg_struct_name)
            dchg_data = resultData.(dchg_struct_name);
            all_soc_levels = fieldnames(dchg_data);
            
            % Filter SOC levels based on user selection
            if isempty(selectedSOC)
                soc_levels = all_soc_levels;
            else
                soc_levels = intersect(selectedSOC, all_soc_levels);
            end
            
            for soc_idx = 1:length(soc_levels)
                soc_level = soc_levels{soc_idx};
                if isfield(dchg_data, soc_level)
                    profiles = fieldnames(dchg_data.(soc_level));
                    for prof_idx = 1:length(profiles)
                        profile_name = profiles{prof_idx};
                        if isfield(dchg_data.(soc_level), profile_name)
                            events = fieldnames(dchg_data.(soc_level).(profile_name));
                            event_count = length(events);
                            dchg_events = dchg_events + event_count;
                            if isKey(dchg_dc_counts, profile_name)
                                dchg_dc_counts(profile_name) = dchg_dc_counts(profile_name) + event_count;
                            else
                                dchg_dc_counts(profile_name) = event_count;
                            end
                        end
                    end
            end
        end
    end
    
        if chg_events > 0 || dchg_events > 0
            fprintf('  %s:\n', channelName);
            if chg_events > 0
                fprintf('    Charging: Total %d events', chg_events);
                if ~isempty(chg_dc_counts)
                    dc_keys = keys(chg_dc_counts);
                    [sorted_dc, sort_idx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match')), dc_keys));
                    sorted_dc_keys = dc_keys(sort_idx);
                    fprintf(' [');
                    for dc_idx = 1:length(sorted_dc_keys)
                        dc_name = sorted_dc_keys{dc_idx};
                        if dc_idx > 1, fprintf(', '); end
                        fprintf('%s:%d', dc_name, chg_dc_counts(dc_name));
                    end
                    fprintf(']\n');
                else
                    fprintf('\n');
                end
            end
            if dchg_events > 0
                fprintf('    Discharging: Total %d events', dchg_events);
                if ~isempty(dchg_dc_counts)
                    dc_keys = keys(dchg_dc_counts);
                    [sorted_dc, sort_idx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match')), dc_keys));
                    sorted_dc_keys = dc_keys(sort_idx);
                    fprintf(' [');
                    for dc_idx = 1:length(sorted_dc_keys)
                        dc_name = sorted_dc_keys{dc_idx};
                        if dc_idx > 1, fprintf(', '); end
                        fprintf('%s:%d', dc_name, dchg_dc_counts(dc_name));
                    end
                    fprintf(']\n');
                else
                    fprintf('\n');
                end
            end
    total_events = total_events + chg_events + dchg_events;
        end
    end
    
    fprintf('Total events for %s: %d\n', cycleType, total_events);
end

%% Cycle-by-Cycle DCIR Trend Analysis
fprintf('\n\n========================================\n');
fprintf('=== Cycle-by-Cycle DCIR Trend Analysis ===\n');
fprintf('========================================\n');

% Initialize cycle trend structure
cycleTrendResults = struct();

% Extract cycle numbers
cycleNumbers = [];
for i = 1:length(cycleTypes)
    cycleType = cycleTypes{i};
    cycleNum = str2double(regexp(cycleType, '\d+', 'match'));
    cycleNumbers = [cycleNumbers, cycleNum];
end
[sortedCycleNums, sortIdx] = sort(cycleNumbers);
sortedCycleTypes = cycleTypes(sortIdx);

fprintf('Analyzing cycles: %s\n', mat2str(sortedCycleNums));

% Process each channel
for cycleIdx = 1:length(sortedCycleTypes)
    cycleType = sortedCycleTypes{cycleIdx};
    
    % Convert cycle type to valid field name (add prefix if starts with number)
    if ~isempty(regexp(cycleType, '^\d', 'once'))
        validFieldName = ['cycle_' cycleType];
    else
        validFieldName = cycleType;
    end
    
    if ~isfield(allCycleResults, validFieldName)
        continue;
    end
    
    resultData = allCycleResults.(validFieldName);
    allChannels = fieldnames(resultData);
    
    for ch_idx = 1:length(allChannels)
        channelName = allChannels{ch_idx};
        
        % Extract channel number
        match = regexp(channelName, 'ch(\d+)_', 'tokens');
        if isempty(match)
            continue;
        end
        chNum = str2double(match{1}{1});
        
        % Use unified channel name (channel number only)
        unifiedChannelName = sprintf('ch%d', chNum);
        
        % Initialize channel structure if needed
        if ~isfield(cycleTrendResults, unifiedChannelName)
            cycleTrendResults.(unifiedChannelName) = struct();
            cycleTrendResults.(unifiedChannelName).DC = struct();
        end
        
        % channelName is already the full name like 'ch16_Drive_0cyc_ChgEvent'
        % Extract base name for structure access
        if contains(channelName, '_ChgEvent')
            baseName = strrep(channelName, '_ChgEvent', '');
            chg_struct_name = [baseName '_ChgEvent'];
        elseif contains(channelName, '_DchEvent')
            baseName = strrep(channelName, '_DchEvent', '');
            chg_struct_name = [baseName '_ChgEvent'];
        else
            chg_struct_name = [channelName '_ChgEvent'];
        end
        
        % Collect DCIR values for each DC profile (DC1~DC8) and time interval
        % Use only charging events for Rchg
        if isfield(resultData, chg_struct_name)
            chgData = resultData.(chg_struct_name);
            socLevels = fieldnames(chgData);
            for soc_idx = 1:length(socLevels)
                socLevel = socLevels{soc_idx};
                profiles = fieldnames(chgData.(socLevel));
                
                for prof_idx = 1:length(profiles)
                    profile = profiles{prof_idx};  % DC1, DC2, ..., DC8
                    
                    % Initialize DC structure if needed
                    if ~isfield(cycleTrendResults.(unifiedChannelName).DC, profile)
                        cycleTrendResults.(unifiedChannelName).DC.(profile) = struct();
                    end
                    
                    % Collect DCIR values for each time interval
                    for dt_idx = 1:length(dt_list)
                        dt_sec = dt_list(dt_idx);
                        fieldName = sprintf('DCIR_%ds', dt_sec);
                        
                        % Initialize field if needed
                        if ~isfield(cycleTrendResults.(unifiedChannelName).DC.(profile), fieldName)
                            cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).cycles = [];
                            cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).mean = [];
                            cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).std = [];
                            cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).event_counts = [];
                        end
                        
                        % Collect DCIR values from all events in this profile (charging events only)
                        dcirVals = [];
                        events = fieldnames(chgData.(socLevel).(profile));
                        eventCount = 0;
                        for evt_idx = 1:length(events)
                            event = events{evt_idx};
                            if isfield(chgData.(socLevel).(profile).(event), fieldName)
                                if isfield(chgData.(socLevel).(profile).(event).(fieldName), 'val')
                                    val = chgData.(socLevel).(profile).(event).(fieldName).val;
                                    if ~isnan(val)
                                        dcirVals = [dcirVals, val];
                                        eventCount = eventCount + 1;
                                    end
                                end
                            end
                        end
                        
                        % Store statistics
                        if ~isempty(dcirVals)
                            cycleNum = str2double(regexp(cycleType, '\d+', 'match'));
                            cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).cycles = ...
                                [cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).cycles, cycleNum];
                            cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).mean = ...
                                [cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).mean, mean(dcirVals)];
                            cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).std = ...
                                [cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).std, std(dcirVals)];
                            cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).event_counts = ...
                                [cycleTrendResults.(unifiedChannelName).DC.(profile).(fieldName).event_counts, eventCount];
                        end
                    end
                end
            end
        end
    end
end

% Save results (visualization is done in DriveCycle_DCIR_Visualization.m)
save('DriveCycle_DCIR_CycleTrend.mat', 'cycleTrendResults');
fprintf('\nCycle trend results saved to: DriveCycle_DCIR_CycleTrend.mat\n');

% Summary table (aggregate across all DC profiles)
fprintf('\n=== Cycle Trend Summary (Slope: mΩ/cycle) ===\n');
fprintf('%-12s', 'Channel');
for dt_idx = 1:length(dt_list)
    fprintf('  %-10s', sprintf('Rchg_%ds', dt_list(dt_idx)));
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 12 + 12*length(dt_list)));

% Get all channels from cycleTrendResults
trendChannels = fieldnames(cycleTrendResults);

for ch_idx = 1:length(trendChannels)
    unifiedChannelName = trendChannels{ch_idx};
    match = regexp(unifiedChannelName, 'ch(\d+)', 'tokens');
    if ~isempty(match)
        chNum = str2double(match{1}{1});
        chDisplayName = sprintf('Ch%d', chNum);
    else
        chDisplayName = unifiedChannelName;
    end
    
    fprintf('%-12s', chDisplayName);
    for dt_idx = 1:length(dt_list)
        dt_sec = dt_list(dt_idx);
        fieldName = sprintf('DCIR_%ds', dt_sec);
        
        % Aggregate data across all DC profiles
        allCycles = [];
        allMeans = [];
        
        if isfield(cycleTrendResults.(unifiedChannelName), 'DC')
            dcProfiles = fieldnames(cycleTrendResults.(unifiedChannelName).DC);
            for dc_idx = 1:length(dcProfiles)
                dcProfile = dcProfiles{dc_idx};
                if isfield(cycleTrendResults.(unifiedChannelName).DC.(dcProfile), fieldName)
                    cycles = cycleTrendResults.(unifiedChannelName).DC.(dcProfile).(fieldName).cycles;
                    means = cycleTrendResults.(unifiedChannelName).DC.(dcProfile).(fieldName).mean;
                    if ~isempty(cycles)
                        allCycles = [allCycles, cycles];
                        allMeans = [allMeans, means];
                    end
                end
            end
        end
        
        if ~isempty(allCycles) && length(unique(allCycles)) >= 2
            % Group by cycle and calculate mean across all DC profiles
            uniqueCycles = unique(allCycles);
            aggregatedMeans = [];
            for cyc_idx = 1:length(uniqueCycles)
                cyc = uniqueCycles(cyc_idx);
                idx = find(allCycles == cyc);
                aggregatedMeans = [aggregatedMeans, mean(allMeans(idx))];
            end
            % Calculate slope (mΩ per cycle)
            p = polyfit(uniqueCycles, aggregatedMeans, 1);
            fprintf('  %+8.4f', p(1));  % Slope in mΩ/cycle
        else
            fprintf('  %10s', 'N/A');
        end
    end
    fprintf('\n');
end

%% Capacity-DCIR Correlation Analysis
if capacityDataAvailable && exist('allChannelsCapacity', 'var')
    fprintf('\n\n========================================\n');
    fprintf('=== Capacity-DCIR Correlation Analysis ===\n');
    fprintf('========================================\n');
    
    % Initialize correlation results structure
    correlationResults = struct();
    
    % Extract cycle numbers from capacity data
    if isfield(allChannelsCapacity, 'Ch09')
        capacityCycles = allChannelsCapacity.Ch09.cycles;
        fprintf('Capacity cycles available: %s\n', mat2str(capacityCycles));
    else
        fprintf('WARNING: Cannot find Ch09 in capacity data\n');
        capacityCycles = [];
    end
    
    % Map cycle types to cycle numbers
    cycleTypeToNumber = containers.Map();
    for i = 1:length(cycleTypes)
        cycleType = cycleTypes{i};
        cycleNum = str2double(regexp(cycleType, '\d+', 'match'));
        cycleTypeToNumber(cycleType) = cycleNum;
    end
    
    % Process each channel (filter based on user selection)
    allChannelNames = {'Ch09', 'Ch10', 'Ch11', 'Ch12', 'Ch13', 'Ch14', 'Ch15', 'Ch16'};
    
    % Filter channels based on user selection
    if isempty(selectedChannels)
        channelNames = allChannelNames;
    else
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
        channelNames = {};
        for i = 1:length(allChannelNames)
            chName = allChannelNames{i};
            chNum = str2double(chName(3:end));
            if ismember(chNum, selectedChannelNums)
                channelNames{end+1} = chName;
            end
        end
    end
    
    fprintf('Processing channels for correlation: %s\n', strjoin(channelNames, ', '));
    
    for ch_idx = 1:length(channelNames)
        chName = channelNames{ch_idx};
        chNum = str2double(chName(3:end));
        channelFieldName = sprintf('ch%d_Drive_', chNum);
        
        % Check if capacity data exists for this channel
        if ~isfield(allChannelsCapacity, chName)
            fprintf('\nSkipping %s: No capacity data\n', chName);
            continue;
        end
        
        channelCapacity = allChannelsCapacity.(chName);
        capacityValues = channelCapacity.capacity;
        capacityCycleNumbers = channelCapacity.cycles;
        
        fprintf('\n--- %s Analysis ---\n', chName);
        fprintf('Capacity values: %s [Ah]\n', mat2str(capacityValues, 4));
        fprintf('Cycle numbers: %s\n', mat2str(capacityCycleNumbers));
        
        % Initialize channel correlation structure
        correlationResults.(chName) = struct();
        correlationResults.(chName).cycles = capacityCycleNumbers;
        correlationResults.(chName).capacity = capacityValues;
        correlationResults.(chName).DCIR = struct();
        
        % Process each cycle
        for cycleIdx = 1:length(cycleTypes)
            cycleType = cycleTypes{cycleIdx};
            cycleNum = cycleTypeToNumber(cycleType);
            
            % Find matching capacity index
            capIdx = find(capacityCycleNumbers == cycleNum, 1);
            if isempty(capIdx)
                fprintf('  %s: No capacity data for cycle %d\n', cycleType, cycleNum);
                continue;
            end
            
            % Get DCIR results for this cycle
            % Convert cycle type to valid field name (add prefix if starts with number)
            if ~isempty(regexp(cycleType, '^\d', 'once'))
                validFieldName = ['cycle_' cycleType];
            else
                validFieldName = cycleType;
            end
            
            if ~isfield(allCycleResults, validFieldName)
                fprintf('  %s: No DCIR results found\n', cycleType);
                continue;
            end
            
            dcirData = allCycleResults.(validFieldName);
            
            % Find channel in DCIR data
            channelPattern = sprintf('%s*', channelFieldName);
            dcirChannels = fieldnames(dcirData);
            matchedChannel = [];
            for i = 1:length(dcirChannels)
                if contains(dcirChannels{i}, channelFieldName)
                    matchedChannel = dcirChannels{i};
                    break;
                end
            end
            
            if isempty(matchedChannel)
                fprintf('  %s: Channel not found in DCIR data\n', cycleType);
                continue;
            end
            
            % Extract DCIR values for different time intervals
            % matchedChannel is already the full name like 'ch16_Drive_0cyc_ChgEvent'
            % Extract base name for structure access
            if contains(matchedChannel, '_ChgEvent')
                baseName = strrep(matchedChannel, '_ChgEvent', '');
                chgStructName = [baseName '_ChgEvent'];
                dchgStructName = [baseName '_DchEvent'];
            elseif contains(matchedChannel, '_DchEvent')
                baseName = strrep(matchedChannel, '_DchEvent', '');
                chgStructName = [baseName '_ChgEvent'];
                dchgStructName = [baseName '_DchEvent'];
            else
                chgStructName = [matchedChannel '_ChgEvent'];
                dchgStructName = [matchedChannel '_DchEvent'];
            end
            
            % Store DCIR values per SOC level for this cycle
            % Convert cycle type to valid field name (add prefix if starts with number)
            if ~isempty(regexp(cycleType, '^\d', 'once'))
                validCycleField = ['cycle_' cycleType];
            else
                validCycleField = cycleType;
            end
            
            % Get available SOC levels from charging events
            availableSocLevels = {};
            if isfield(dcirData, chgStructName)
                availableSocLevels = fieldnames(dcirData.(chgStructName));
            elseif isfield(dcirData, dchgStructName)
                availableSocLevels = fieldnames(dcirData.(dchgStructName));
            end
            
            % Filter SOC levels based on user selection
            if isempty(selectedSOC)
                socLevels = availableSocLevels;
            else
                socLevels = intersect(selectedSOC, availableSocLevels);
            end
            
            % Initialize SOC structure if needed
            if ~isfield(correlationResults.(chName), 'SOC')
                correlationResults.(chName).SOC = struct();
            end
            
            % Process each SOC level separately
            for soc_idx = 1:length(socLevels)
                socLevel = socLevels{soc_idx};
                
                if ~isfield(correlationResults.(chName).SOC, socLevel)
                    correlationResults.(chName).SOC.(socLevel) = struct();
                    correlationResults.(chName).SOC.(socLevel).DCIR = struct();
                end
                
                % Store DCIR values per DC profile for this cycle and SOC level
                % Get available DC profiles
                availableProfiles = {};
                if isfield(dcirData, chgStructName) && isfield(dcirData.(chgStructName), socLevel)
                    availableProfiles = fieldnames(dcirData.(chgStructName).(socLevel));
                elseif isfield(dcirData, dchgStructName) && isfield(dcirData.(dchgStructName), socLevel)
                    availableProfiles = fieldnames(dcirData.(dchgStructName).(socLevel));
                end
                
                % Initialize DC structure if needed
                if ~isfield(correlationResults.(chName).SOC.(socLevel), 'DC')
                    correlationResults.(chName).SOC.(socLevel).DC = struct();
                end
                
                % Process each DC profile separately if correlationByDC is true, otherwise aggregate
                if correlationByDC
                    % DC별로 개별 저장
                    for prof_idx = 1:length(availableProfiles)
                        profile = availableProfiles{prof_idx};
                        
                        if ~isfield(correlationResults.(chName).SOC.(socLevel).DC, profile)
                            correlationResults.(chName).SOC.(socLevel).DC.(profile) = struct();
                            correlationResults.(chName).SOC.(socLevel).DC.(profile).DCIR = struct();
                        end
                        
                        dcirValues = struct();
                        for dt_idx = 1:length(dt_list)
                            dt_sec = dt_list(dt_idx);
                            fieldName = sprintf('DCIR_%ds', dt_sec);
                            
                            % Collect DCIR values from this DC profile only
                            dcirVals = [];
                            
                            % Charging events
                            if isfield(dcirData, chgStructName) && isfield(dcirData.(chgStructName), socLevel) && isfield(dcirData.(chgStructName).(socLevel), profile)
                                events = fieldnames(dcirData.(chgStructName).(socLevel).(profile));
                                for evt_idx = 1:length(events)
                                    event = events{evt_idx};
                                    if isfield(dcirData.(chgStructName).(socLevel).(profile).(event), fieldName)
                                        if isfield(dcirData.(chgStructName).(socLevel).(profile).(event).(fieldName), 'val')
                                            val = dcirData.(chgStructName).(socLevel).(profile).(event).(fieldName).val;
                                            if ~isnan(val)
                                                dcirVals = [dcirVals, val];
                                            end
                                        end
                                    end
                                end
                            end
                            
                            % Discharging events
                            if isfield(dcirData, dchgStructName) && isfield(dcirData.(dchgStructName), socLevel) && isfield(dcirData.(dchgStructName).(socLevel), profile)
                                events = fieldnames(dcirData.(dchgStructName).(socLevel).(profile));
                                for evt_idx = 1:length(events)
                                    event = events{evt_idx};
                                    if isfield(dcirData.(dchgStructName).(socLevel).(profile).(event), fieldName)
                                        if isfield(dcirData.(dchgStructName).(socLevel).(profile).(event).(fieldName), 'val')
                                            val = dcirData.(dchgStructName).(socLevel).(profile).(event).(fieldName).val;
                                            if ~isnan(val)
                                                dcirVals = [dcirVals, val];
                                            end
                                        end
                                    end
                                end
                            end
                            
                            if ~isempty(dcirVals)
                                dcirValues.(fieldName).mean = mean(dcirVals);
                                dcirValues.(fieldName).std = std(dcirVals);
                                dcirValues.(fieldName).median = median(dcirVals);
                                dcirValues.(fieldName).all = dcirVals;
                            else
                                dcirValues.(fieldName).mean = NaN;
                                dcirValues.(fieldName).std = NaN;
                                dcirValues.(fieldName).median = NaN;
                                dcirValues.(fieldName).all = [];
                            end
                        end
                        
                        % Store DCIR values for this cycle, SOC level, and DC profile
                        correlationResults.(chName).SOC.(socLevel).DC.(profile).DCIR.(validCycleField) = dcirValues;
                    end
                else
                    % 모든 DC 프로파일 평균 (기존 방식)
                    dcirValues = struct();
                    for dt_idx = 1:length(dt_list)
                        dt_sec = dt_list(dt_idx);
                        fieldName = sprintf('DCIR_%ds', dt_sec);
                        
                        % Collect DCIR values from all DC profiles for this SOC level
                        dcirVals = [];
                        
                        % Charging events
                        if isfield(dcirData, chgStructName) && isfield(dcirData.(chgStructName), socLevel)
                            chgData = dcirData.(chgStructName).(socLevel);
                            profiles = fieldnames(chgData);
                            for prof_idx = 1:length(profiles)
                                profile = profiles{prof_idx};
                                events = fieldnames(chgData.(profile));
                                for evt_idx = 1:length(events)
                                    event = events{evt_idx};
                                    if isfield(chgData.(profile).(event), fieldName)
                                        if isfield(chgData.(profile).(event).(fieldName), 'val')
                                            val = chgData.(profile).(event).(fieldName).val;
                                            if ~isnan(val)
                                                dcirVals = [dcirVals, val];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        
                        % Discharging events
                        if isfield(dcirData, dchgStructName) && isfield(dcirData.(dchgStructName), socLevel)
                            dchgData = dcirData.(dchgStructName).(socLevel);
                            profiles = fieldnames(dchgData);
                            for prof_idx = 1:length(profiles)
                                profile = profiles{prof_idx};
                                events = fieldnames(dchgData.(profile));
                                for evt_idx = 1:length(events)
                                    event = events{evt_idx};
                                    if isfield(dchgData.(profile).(event), fieldName)
                                        if isfield(dchgData.(profile).(event).(fieldName), 'val')
                                            val = dchgData.(profile).(event).(fieldName).val;
                                            if ~isnan(val)
                                                dcirVals = [dcirVals, val];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        
                        if ~isempty(dcirVals)
                            dcirValues.(fieldName).mean = mean(dcirVals);
                            dcirValues.(fieldName).std = std(dcirVals);
                            dcirValues.(fieldName).median = median(dcirVals);
                            dcirValues.(fieldName).all = dcirVals;
                        else
                            dcirValues.(fieldName).mean = NaN;
                            dcirValues.(fieldName).std = NaN;
                            dcirValues.(fieldName).median = NaN;
                            dcirValues.(fieldName).all = [];
                        end
                    end
                    
                    % Store DCIR values for this cycle and SOC level (all DCs averaged)
                    correlationResults.(chName).SOC.(socLevel).DCIR.(validCycleField) = dcirValues;
                end
            end
        end
        
        % Calculate correlation for each SOC level and time interval
        for soc_idx = 1:length(socLevels)
            socLevel = socLevels{soc_idx};
            
            if ~isfield(correlationResults.(chName).SOC, socLevel)
                continue;
            end
            
            if correlationByDC
                % DC별 상관분석
                if ~isfield(correlationResults.(chName).SOC.(socLevel), 'DC')
                    continue;
                end
                
                % Get available DC profiles
                availableProfiles = fieldnames(correlationResults.(chName).SOC.(socLevel).DC);
                if isempty(availableProfiles)
                    continue;
                end
                
                % Sort DC profiles (DC1, DC2, ...)
                dc_nums = zeros(length(availableProfiles), 1);
                for d = 1:length(availableProfiles)
                    dc_str = availableProfiles{d};
                    num_match = regexp(dc_str, 'DC(\d+)', 'tokens');
                    if ~isempty(num_match)
                        dc_nums(d) = str2double(num_match{1}{1});
                    else
                        dc_nums(d) = 999;
                    end
                end
                [~, dc_sort_idx] = sort(dc_nums);
                sortedProfiles = availableProfiles(dc_sort_idx);
                
                % Process each DC profile
                for dc_idx = 1:length(sortedProfiles)
                    profile = sortedProfiles{dc_idx};
                    
                    fprintf('\n  === Capacity-Rchg Correlation Analysis (%s - %s - %s) ===\n', chName, socLevel, profile);
                    fprintf('  Time Interval |    R     |   R²    |  Slope (mΩ/Ah) |  Intercept  |  n  | Cycles\n');
                    fprintf('  %s\n', repmat('-', 1, 85));
                    
                    for dt_idx = 1:length(dt_list)
                        dt_sec = dt_list(dt_idx);
                        fieldName = sprintf('DCIR_%ds', dt_sec);
                        displayName = sprintf('Rchg_%ds', dt_sec);
                        
                        % Collect capacity and DCIR data across cycles for this SOC level and DC profile
                        capData = [];
                        dcirData = [];
                        cycleNumbers = [];
                        
                        for cycleIdx = 1:length(cycleTypes)
                            cycleType = cycleTypes{cycleIdx};
                            cycleNum = cycleTypeToNumber(cycleType);
                            capIdx = find(capacityCycleNumbers == cycleNum, 1);
                            
                            % Convert cycle type to valid field name
                            if ~isempty(regexp(cycleType, '^\d', 'once'))
                                validCycleField = ['cycle_' cycleType];
                            else
                                validCycleField = cycleType;
                            end
                            
                            if ~isempty(capIdx) && isfield(correlationResults.(chName).SOC.(socLevel).DC.(profile).DCIR, validCycleField)
                                if isfield(correlationResults.(chName).SOC.(socLevel).DC.(profile).DCIR.(validCycleField), fieldName)
                                    dcirMean = correlationResults.(chName).SOC.(socLevel).DC.(profile).DCIR.(validCycleField).(fieldName).mean;
                                    if ~isnan(dcirMean)
                                        capData = [capData, capacityValues(capIdx)];
                                        dcirData = [dcirData, dcirMean];
                                        cycleNumbers = [cycleNumbers, cycleNum];
                                    end
                                end
                            end
                        end
                        
                        if length(capData) >= 3
                            % Calculate correlation coefficient
                            R = corrcoef(capData, dcirData);
                            corrCoeff = R(1, 2);
                            
                            % Linear regression
                            p = polyfit(capData, dcirData, 1);
                            fitLine = polyval(p, capData);
                            
                            % Initialize Correlation structure for this SOC level and DC profile if needed
                            if ~isfield(correlationResults.(chName).SOC.(socLevel).DC.(profile), 'Correlation')
                                correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation = struct();
                            end
                            
                            % Store correlation results
                            correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation.(fieldName).capacity = capData;
                            correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation.(fieldName).DCIR = dcirData;
                            correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation.(fieldName).cycles = cycleNumbers;
                            correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation.(fieldName).R = corrCoeff;
                            correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation.(fieldName).R_squared = corrCoeff^2;
                            correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation.(fieldName).slope = p(1);
                            correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation.(fieldName).intercept = p(2);
                            
                            % Format cycles string
                            cycles_str = sprintf('%d', cycleNumbers(1));
                            for cyc_idx = 2:length(cycleNumbers)
                                cycles_str = sprintf('%s,%d', cycles_str, cycleNumbers(cyc_idx));
                            end
                            
                            fprintf('  %-13s | %7.3f | %7.3f | %14.2f | %11.2f | %3d | %s\n', ...
                                displayName, corrCoeff, corrCoeff^2, p(1), p(2), length(capData), cycles_str);
                        else
                            fprintf('  %-13s | %7s | %7s | %14s | %11s | %3d | %s\n', ...
                                displayName, 'N/A', 'N/A', 'N/A', 'N/A', length(capData), 'N/A');
                        end
                    end
                end
            else
                % 모든 DC 평균 상관분석 (기존 방식)
                if ~isfield(correlationResults.(chName).SOC.(socLevel), 'DCIR')
                    continue;
                end
                
                fprintf('\n  === Capacity-Rchg Correlation Analysis (%s - %s) ===\n', chName, socLevel);
                fprintf('  Time Interval |    R     |   R²    |  Slope (mΩ/Ah) |  Intercept  |  n  | Cycles\n');
                fprintf('  %s\n', repmat('-', 1, 85));
                
                for dt_idx = 1:length(dt_list)
                    dt_sec = dt_list(dt_idx);
                    fieldName = sprintf('DCIR_%ds', dt_sec);
                    displayName = sprintf('Rchg_%ds', dt_sec);
                    
                    % Collect capacity and DCIR data across cycles for this SOC level
                    capData = [];
                    dcirData = [];
                    cycleNumbers = [];
                    
                    for cycleIdx = 1:length(cycleTypes)
                        cycleType = cycleTypes{cycleIdx};
                        cycleNum = cycleTypeToNumber(cycleType);
                        capIdx = find(capacityCycleNumbers == cycleNum, 1);
                        
                        % Convert cycle type to valid field name
                        if ~isempty(regexp(cycleType, '^\d', 'once'))
                            validCycleField = ['cycle_' cycleType];
                        else
                            validCycleField = cycleType;
                        end
                        
                        if ~isempty(capIdx) && isfield(correlationResults.(chName).SOC.(socLevel).DCIR, validCycleField)
                            if isfield(correlationResults.(chName).SOC.(socLevel).DCIR.(validCycleField), fieldName)
                                dcirMean = correlationResults.(chName).SOC.(socLevel).DCIR.(validCycleField).(fieldName).mean;
                                if ~isnan(dcirMean)
                                    capData = [capData, capacityValues(capIdx)];
                                    dcirData = [dcirData, dcirMean];
                                    cycleNumbers = [cycleNumbers, cycleNum];
                                end
                            end
                        end
                    end
                    
                    if length(capData) >= 3
                        % Calculate correlation coefficient
                        R = corrcoef(capData, dcirData);
                        corrCoeff = R(1, 2);
                        
                        % Linear regression
                        p = polyfit(capData, dcirData, 1);
                        fitLine = polyval(p, capData);
                        
                        % Initialize Correlation structure for this SOC level if needed
                        if ~isfield(correlationResults.(chName).SOC.(socLevel), 'Correlation')
                            correlationResults.(chName).SOC.(socLevel).Correlation = struct();
                        end
                        
                        % Store correlation results
                        correlationResults.(chName).SOC.(socLevel).Correlation.(fieldName).capacity = capData;
                        correlationResults.(chName).SOC.(socLevel).Correlation.(fieldName).DCIR = dcirData;
                        correlationResults.(chName).SOC.(socLevel).Correlation.(fieldName).cycles = cycleNumbers;
                        correlationResults.(chName).SOC.(socLevel).Correlation.(fieldName).R = corrCoeff;
                        correlationResults.(chName).SOC.(socLevel).Correlation.(fieldName).R_squared = corrCoeff^2;
                        correlationResults.(chName).SOC.(socLevel).Correlation.(fieldName).slope = p(1);
                        correlationResults.(chName).SOC.(socLevel).Correlation.(fieldName).intercept = p(2);
                        
                        % Format cycles string
                        cycles_str = sprintf('%d', cycleNumbers(1));
                        for cyc_idx = 2:length(cycleNumbers)
                            cycles_str = sprintf('%s,%d', cycles_str, cycleNumbers(cyc_idx));
                        end
                        
                        fprintf('  %-13s | %7.3f | %7.3f | %14.2f | %11.2f | %3d | %s\n', ...
                            displayName, corrCoeff, corrCoeff^2, p(1), p(2), length(capData), cycles_str);
                    else
                        fprintf('  %-13s | %7s | %7s | %14s | %11s | %3d | %s\n', ...
                            displayName, 'N/A', 'N/A', 'N/A', 'N/A', length(capData), 'N/A');
                    end
                end
            end
        end
    end
    
    % Save correlation results
    save('DriveCycle_Capacity_DCIR_Correlation.mat', 'correlationResults');
    fprintf('\nCorrelation results saved to: DriveCycle_Capacity_DCIR_Correlation.mat\n');
    
    % Summary table
    if correlationByDC
        % DC별 Summary table
        fprintf('\n=== Correlation Summary (R values) - DC별 ===\n');
        for ch_idx = 1:length(channelNames)
            chName = channelNames{ch_idx};
            if ~isfield(correlationResults, chName) || ~isfield(correlationResults.(chName), 'SOC')
                continue;
            end
            
            socLevels = fieldnames(correlationResults.(chName).SOC);
            for soc_idx = 1:length(socLevels)
                socLevel = socLevels{soc_idx};
                if ~isfield(correlationResults.(chName).SOC.(socLevel), 'DC')
                    continue;
                end
                
                availableProfiles = fieldnames(correlationResults.(chName).SOC.(socLevel).DC);
                if isempty(availableProfiles)
                    continue;
                end
                
                % Sort DC profiles
                dc_nums = zeros(length(availableProfiles), 1);
                for d = 1:length(availableProfiles)
                    dc_str = availableProfiles{d};
                    num_match = regexp(dc_str, 'DC(\d+)', 'tokens');
                    if ~isempty(num_match)
                        dc_nums(d) = str2double(num_match{1}{1});
                    else
                        dc_nums(d) = 999;
                    end
                end
                [~, dc_sort_idx] = sort(dc_nums);
                sortedProfiles = availableProfiles(dc_sort_idx);
                
                fprintf('\n%s - %s:\n', chName, socLevel);
                fprintf('%-8s', 'DC');
                for dt_idx = 1:length(dt_list)
                    fprintf('  %-8s', sprintf('Rchg_%ds', dt_list(dt_idx)));
                end
                fprintf('\n');
                fprintf('%s\n', repmat('-', 1, 8 + 9*length(dt_list)));
                
                for dc_idx = 1:length(sortedProfiles)
                    profile = sortedProfiles{dc_idx};
                    fprintf('%-8s', profile);
                    for dt_idx = 1:length(dt_list)
                        dt_sec = dt_list(dt_idx);
                        fieldName = sprintf('DCIR_%ds', dt_sec);
                        if isfield(correlationResults.(chName).SOC.(socLevel).DC.(profile), 'Correlation') && ...
                           isfield(correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation, fieldName)
                            R = correlationResults.(chName).SOC.(socLevel).DC.(profile).Correlation.(fieldName).R;
                            fprintf('  %8.3f', R);
                        else
                            fprintf('  %8s', 'N/A');
                        end
                    end
                    fprintf('\n');
                end
            end
        end
    else
        % 기존 Summary table (모든 DC 평균)
        fprintf('\n=== Correlation Summary (R values) ===\n');
        fprintf('%-10s', 'Channel');
        for dt_idx = 1:length(dt_list)
            fprintf('  %-8s', sprintf('Rchg_%ds', dt_list(dt_idx)));
        end
        fprintf('\n');
        fprintf('%s\n', repmat('-', 1, 10 + 9*length(dt_list)));
        
        for ch_idx = 1:length(channelNames)
            chName = channelNames{ch_idx};
            if isfield(correlationResults, chName) && isfield(correlationResults.(chName), 'SOC')
                socLevels = fieldnames(correlationResults.(chName).SOC);
                for soc_idx = 1:length(socLevels)
                    socLevel = socLevels{soc_idx};
                    if isfield(correlationResults.(chName).SOC.(socLevel), 'Correlation')
                        fprintf('%-10s', sprintf('%s-%s', chName, socLevel));
                        for dt_idx = 1:length(dt_list)
                            dt_sec = dt_list(dt_idx);
                            fieldName = sprintf('DCIR_%ds', dt_sec);
                            if isfield(correlationResults.(chName).SOC.(socLevel).Correlation, fieldName)
                                R = correlationResults.(chName).SOC.(socLevel).Correlation.(fieldName).R;
                                fprintf('  %8.3f', R);
                            else
                                fprintf('  %8s', 'N/A');
                            end
                        end
                        fprintf('\n');
                    end
                end
            elseif isfield(correlationResults, chName) && isfield(correlationResults.(chName), 'Correlation')
                fprintf('%-10s', chName);
                for dt_idx = 1:length(dt_list)
                    dt_sec = dt_list(dt_idx);
                    fieldName = sprintf('DCIR_%ds', dt_sec);
                    if isfield(correlationResults.(chName).Correlation, fieldName)
                        R = correlationResults.(chName).Correlation.(fieldName).R;
                        fprintf('  %8.3f', R);
                    else
                        fprintf('  %8s', 'N/A');
                    end
                end
                fprintf('\n');
            end
        end
    end
    
    % Correlation visualization is done in DriveCycle_DCIR_Visualization.m
    
else
    fprintf('\nSkipping correlation analysis: Capacity data not available\n');
end

%% Create Summary Table
fprintf('\n=== Creating Summary Table ===\n');

% Initialize table data
tableData = [];

% Iterate through all cycles
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    
    % Extract cycle number from cycle type
    match = regexp(cycleType, '(\d+)cyc', 'tokens');
    if ~isempty(match)
        cycleNum = str2double(match{1}{1});
    else
        cycleNum = NaN;
    end
    
    % Convert cycle type to valid field name
    if ~isempty(regexp(cycleType, '^\d', 'once'))
        validCycleField = ['cycle_' cycleType];
    else
        validCycleField = cycleType;
    end
    
    % Load cycle data
    dataFile = fullfile(dataDir, sprintf('parsedDriveCycle_%s_filtered.mat', cycleType));
    if ~exist(dataFile, 'file')
        continue;
    end
    
    dataVarName = sprintf('parsedDriveCycle_%s', cycleType);
    load(dataFile, dataVarName);
    cycleData = eval(dataVarName);
    
    % Get capacity for this cycle
    capIdx = find(capacityCycleNumbers == cycleNum, 1);
    if isempty(capIdx)
        capacity = NaN;
    else
        capacity = capacityValues(capIdx);
    end
    
        % Iterate through all channels
        for ch_idx = 1:length(channelNames)
            chName = channelNames{ch_idx};
            
            % Extract channel number from channel name
            match = regexp(chName, 'ch(\d+)', 'tokens', 'ignorecase');
            if ~isempty(match)
                chNum = str2double(match{1}{1});
            else
                chNum = NaN;
            end
        
        % Check if this channel exists in cycle data
        channelFieldName = sprintf('ch%d', chNum);
        chgStructName = sprintf('%s_Drive_%s_ChgEvent', channelFieldName, cycleType);
        dchgStructName = sprintf('%s_Drive_%s_DchEvent', channelFieldName, cycleType);
        
        if ~isfield(cycleData, chgStructName) && ~isfield(cycleData, dchgStructName)
            continue;
        end
        
        % Get available SOC levels
        availableSocLevels = {};
        if isfield(cycleData, chgStructName)
            availableSocLevels = fieldnames(cycleData.(chgStructName));
        elseif isfield(cycleData, dchgStructName)
            availableSocLevels = fieldnames(cycleData.(dchgStructName));
        end
        
        % Filter SOC levels based on user selection
        if isempty(selectedSOC)
            socLevels = availableSocLevels;
        else
            socLevels = intersect(selectedSOC, availableSocLevels);
        end
        
        % Iterate through SOC levels
        for soc_idx = 1:length(socLevels)
            socLevel = socLevels{soc_idx};
            
            % Get available DC profiles
            availableProfiles = {};
            if isfield(cycleData, chgStructName) && isfield(cycleData.(chgStructName), socLevel)
                availableProfiles = fieldnames(cycleData.(chgStructName).(socLevel));
            elseif isfield(cycleData, dchgStructName) && isfield(cycleData.(dchgStructName), socLevel)
                availableProfiles = fieldnames(cycleData.(dchgStructName).(socLevel));
            end
            
            % Iterate through DC profiles
            for prof_idx = 1:length(availableProfiles)
                profile = availableProfiles{prof_idx};
                
                % Collect Rchg values for all time intervals
                rchg_vals = struct();
                for dt_idx = 1:length(dt_list)
                    dt_sec = dt_list(dt_idx);
                    fieldName = sprintf('DCIR_%ds', dt_sec);
                    rchg_vals.(sprintf('R_%ds', dt_sec)) = [];
                end
                
                % Collect from charging events
                if isfield(cycleData, chgStructName) && isfield(cycleData.(chgStructName), socLevel) && ...
                   isfield(cycleData.(chgStructName).(socLevel), profile)
                    events = fieldnames(cycleData.(chgStructName).(socLevel).(profile));
                    for evt_idx = 1:length(events)
                        event = events{evt_idx};
                        for dt_idx = 1:length(dt_list)
                            dt_sec = dt_list(dt_idx);
                            fieldName = sprintf('DCIR_%ds', dt_sec);
                            if isfield(cycleData.(chgStructName).(socLevel).(profile).(event), fieldName)
                                if isfield(cycleData.(chgStructName).(socLevel).(profile).(event).(fieldName), 'val')
                                    val = cycleData.(chgStructName).(socLevel).(profile).(event).(fieldName).val;
                                    if ~isnan(val)
                                        rchg_vals.(sprintf('R_%ds', dt_sec)) = [rchg_vals.(sprintf('R_%ds', dt_sec)), val];
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Collect from discharging events
                if isfield(cycleData, dchgStructName) && isfield(cycleData.(dchgStructName), socLevel) && ...
                   isfield(cycleData.(dchgStructName).(socLevel), profile)
                    events = fieldnames(cycleData.(dchgStructName).(socLevel).(profile));
                    for evt_idx = 1:length(events)
                        event = events{evt_idx};
                        for dt_idx = 1:length(dt_list)
                            dt_sec = dt_list(dt_idx);
                            fieldName = sprintf('DCIR_%ds', dt_sec);
                            if isfield(cycleData.(dchgStructName).(socLevel).(profile).(event), fieldName)
                                if isfield(cycleData.(dchgStructName).(socLevel).(profile).(event).(fieldName), 'val')
                                    val = cycleData.(dchgStructName).(socLevel).(profile).(event).(fieldName).val;
                                    if ~isnan(val)
                                        rchg_vals.(sprintf('R_%ds', dt_sec)) = [rchg_vals.(sprintf('R_%ds', dt_sec)), val];
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Calculate mean Rchg values for this DC profile
                rowData = struct();
                rowData.Cycle = cycleNum;
                rowData.Channel = chNum;
                rowData.DC_Profile = profile;
                rowData.Capacity = capacity;
                
                for dt_idx = 1:length(dt_list)
                    dt_sec = dt_list(dt_idx);
                    fieldName = sprintf('R_%ds', dt_sec);
                    if ~isempty(rchg_vals.(fieldName))
                        rowData.(fieldName) = mean(rchg_vals.(fieldName));
                    else
                        rowData.(fieldName) = NaN;
                    end
                end
                
                % Add to table data
                tableData = [tableData; rowData];
            end
        end
    end
end

% Convert to table
if ~isempty(tableData)
    summaryTable = struct2table(tableData);
    
    % Reorder columns: Cycle, Channel, DC_Profile, Capacity, R_1s, R_3s, R_5s, R_10s, R_30s, R_60s
    colOrder = {'Cycle', 'Channel', 'DC_Profile', 'Capacity'};
    for dt_idx = 1:length(dt_list)
        dt_sec = dt_list(dt_idx);
        colOrder{end+1} = sprintf('R_%ds', dt_sec);
    end
    summaryTable = summaryTable(:, colOrder);
    
    % Save table
    save('DriveCycle_Summary_Table.mat', 'summaryTable');
    fprintf('Summary table created: %d rows\n', height(summaryTable));
    fprintf('Saved to: DriveCycle_Summary_Table.mat\n');
    
    % Display first few rows
    fprintf('\nFirst 5 rows of summary table:\n');
    disp(summaryTable(1:min(5, height(summaryTable)), :));
else
    fprintf('WARNING: No data found for summary table\n');
end

%% Start Visualization
fprintf('\nStarting visualization script...\n');
% DriveCycle_DCIR_Visualization; 
