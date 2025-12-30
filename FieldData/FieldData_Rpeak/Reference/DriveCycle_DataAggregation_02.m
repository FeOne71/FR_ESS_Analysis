%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 02_DriveCycle_DataAggregation.m
% 데이터 취합 및 Summary Table 생성
% 
% 목적: 
% - 01번 스크립트가 생성한 Lab_DC_DCIR_*cyc_Events.mat 파일들을 로드
% - Capacity_Data_Static.mat에서 용량 데이터 로드
% - 모든 데이터를 하나의 Table로 취합
% - 컬럼: Cycle, Channel, DC_Profile, Capacity, R_1s, R_3s, R_5s, R_10s, R_30s, R_60s
%
% 입력:
% - Lab_DC_DCIR_*cyc_Events.mat (01번 스크립트 출력)
% - Capacity_Data_Static.mat (용량 데이터)
%
% 출력:
% - DriveCycle_Summary_Table.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

fprintf('=== Drive Cycle Data Aggregation ===\n');

%% Configuration - User Settings
% =========================================================================
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
capacityDataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Capacity_Trend_Figures\Capacity_Data_Static.mat';
inputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\Events Analysis v01';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

selectedCycles = {'0cyc','200cyc','400cyc','600cyc','800cyc'};  
selectedSOC = {'SOC50'};  
% =========================================================================

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
    fprintf('Capacity column will be NaN\n');
    capacityDataAvailable = false;
end

%% Auto-detect available cycle result files
fprintf('\n=== Auto-detecting available cycle result files ===\n');
matFiles = dir(fullfile(inputDir, 'Lab_DC_DCIR_*cyc_Events.mat'));
availableCycles = {};

for i = 1:length(matFiles)
    filename = matFiles(i).name;
    match = regexp(filename, 'Lab_DC_DCIR_(\d+cyc)_Events\.mat', 'tokens');
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

%% Settings
dt_list = [1, 3, 5, 10, 30, 60];

%% Initialize table data
fprintf('\n=== Creating Summary Table ===\n');
tableData = [];

%% Extract channel names from first cycle (assuming all cycles have same channels)
if ~isempty(cycleTypes)
    firstCycleFile = fullfile(inputDir, sprintf('Lab_DC_DCIR_%s_Events.mat', cycleTypes{1}));
    if exist(firstCycleFile, 'file')
        load(firstCycleFile);
        resultVarName = sprintf('Lab_DC_DCIR_%s', cycleTypes{1});
        eval(sprintf('firstCycleData = %s;', resultVarName));
        allFields = fieldnames(firstCycleData);
        
        % Extract unique channel names
        % Use 'Ch09', 'Ch10', ... format to match allChannelsCapacity structure
        channelNames = {};
        for i = 1:length(allFields)
            fieldName = allFields{i};
            if contains(fieldName, '_ChgEvent')
                channelName = strrep(fieldName, '_ChgEvent', '');
                % Extract channel number (e.g., 'ch13_Drive_0cyc' -> 13)
                match = regexp(channelName, 'ch(\d+)', 'tokens', 'ignorecase');
                if ~isempty(match)
                    chNum = str2double(match{1}{1});
                    chName = sprintf('Ch%02d', chNum);  % Ch09, Ch10, ... format
                    if ~ismember(chName, channelNames)
                        channelNames{end+1} = chName;
                    end
                end
            end
        end
        % Sort by channel number
        chNums = cellfun(@(x) str2double(x(3:end)), channelNames);
        [~, sortIdx] = sort(chNums);
        channelNames = channelNames(sortIdx);
        fprintf('Detected channels: %s\n', strjoin(channelNames, ', '));
    end
end

%% Iterate through all cycles
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    
    % Extract cycle number from cycle type
    match = regexp(cycleType, '(\d+)cyc', 'tokens');
    if ~isempty(match)
        cycleNum = str2double(match{1}{1});
    else
        cycleNum = NaN;
    end
    
    fprintf('\n--- Processing %s (Cycle %d) ---\n', cycleType, cycleNum);
    
    % Load cycle result file
    resultFile = fullfile(inputDir, sprintf('Lab_DC_DCIR_%s_Events.mat', cycleType));
    if ~exist(resultFile, 'file')
        fprintf('WARNING: Result file not found: %s\n', resultFile);
        continue;
    end
    
    load(resultFile);
    resultVarName = sprintf('Lab_DC_DCIR_%s', cycleType);
    eval(sprintf('cycleData = %s;', resultVarName));
    
    % Get capacity for this cycle (if available)
    if capacityDataAvailable
        % Get capacity for all channels (will be matched per channel later)
        % For now, we'll get it per channel in the loop below
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
        
        % Get capacity for this channel and cycle
        % chName is in 'Ch09' format, which matches allChannelsCapacity structure
        if capacityDataAvailable && isfield(allChannelsCapacity, chName)
            channelCapacity = allChannelsCapacity.(chName);
            capacityValues = channelCapacity.capacity;
            capacityCycleNumbers = channelCapacity.cycles;
            capIdx = find(capacityCycleNumbers == cycleNum, 1);
            if ~isempty(capIdx)
                capacity = capacityValues(capIdx);
            else
                capacity = NaN;
                fprintf('  WARNING: No capacity data for %s cycle %d\n', chName, cycleNum);
            end
        else
            capacity = NaN;
            if capacityDataAvailable
                fprintf('  WARNING: Channel %s not found in capacity data\n', chName);
            end
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
                
                % Process charging and discharging events separately
                for eventType = {'Charge', 'Discharge'}
                    eventTypeStr = eventType{1};
                    
                    % Collect Rchg values for all time intervals
                    rchg_vals = struct();
                    for dt_idx = 1:length(dt_list)
                        dt_sec = dt_list(dt_idx);
                        rchg_vals.(sprintf('R_%ds', dt_sec)) = [];
                    end
                    
                    % Collect from charging events
                    if strcmp(eventTypeStr, 'Charge')
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
                    end
                    
                    % Collect from discharging events
                    if strcmp(eventTypeStr, 'Discharge')
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
                    end
                    
                    % Calculate mean Rchg values for this DC profile and event type
                    % Only add row if there is at least one valid Rchg value
                    hasData = false;
                    for dt_idx = 1:length(dt_list)
                        dt_sec = dt_list(dt_idx);
                        fieldName = sprintf('R_%ds', dt_sec);
                        if ~isempty(rchg_vals.(fieldName))
                            hasData = true;
                            break;
                        end
                    end
                    
                    if hasData
                        rowData = struct();
                        rowData.Cycle = cycleNum;
                        rowData.Channel = chNum;
                        rowData.DC_Profile = profile;
                        rowData.EventType = eventTypeStr;
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
    end
end

%% Convert to table
if ~isempty(tableData)
    summaryTable = struct2table(tableData);
    
    % Reorder columns: Cycle, Channel, DC_Profile, EventType, Capacity, R_1s, R_3s, R_5s, R_10s, R_30s, R_60s
    colOrder = {'Cycle', 'Channel', 'DC_Profile', 'EventType', 'Capacity'};
    for dt_idx = 1:length(dt_list)
        dt_sec = dt_list(dt_idx);
        colOrder{end+1} = sprintf('R_%ds', dt_sec);
    end
    summaryTable = summaryTable(:, colOrder);
    
    % Save table
    savePath = fullfile(outputDir, 'DriveCycle_Summary_Table.mat');
    save(savePath, 'summaryTable');
    fprintf('\n=== Summary Table Created ===\n');
    fprintf('Total rows: %d\n', height(summaryTable));
    fprintf('Saved to: %s\n', savePath);
    
    % Display first few rows
    fprintf('\nFirst 5 rows of summary table:\n');
    disp(summaryTable(1:min(5, height(summaryTable)), :));
    
    % Display summary statistics
    fprintf('\nSummary Statistics:\n');
    fprintf('  Cycles: %s\n', mat2str(unique(summaryTable.Cycle)'));
    fprintf('  Channels: %s\n', mat2str(unique(summaryTable.Channel)'));
    fprintf('  DC Profiles: %s\n', strjoin(unique(summaryTable.DC_Profile), ', '));
    if ismember('EventType', summaryTable.Properties.VariableNames)
        fprintf('  Event Types: %s\n', strjoin(unique(summaryTable.EventType), ', '));
    end
else
    fprintf('WARNING: No data found for summary table\n');
end

fprintf('\n=== Data Aggregation Complete ===\n');

