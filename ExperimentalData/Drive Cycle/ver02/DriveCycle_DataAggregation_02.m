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
rptDcirDataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_Final\DCIR_SOC_data_all_channels_final.mat';
% 현재 스크립트가 있는 폴더의 Results 폴더 사용
scriptDir = fileparts(mfilename('fullpath'));
inputDir = fullfile(scriptDir, 'Results');
outputDir = fullfile(scriptDir, 'Results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

selectedCycles = {'0cyc','200cyc','400cyc','600cyc','800cyc'};  
selectedSOC = {'SOC50', 'SOC70', 'SOC90'};  % SOC50, SOC70, SOC90 모두 포함  
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
    fprintf('Capacity columns will be NaN\n');
    capacityDataAvailable = false;
end

%% Load RPT-DCIR Data
fprintf('\n=== Loading RPT-DCIR Data ===\n');
rptDataAvailable = false;
if ~isempty(rptDcirDataFile) && exist(rptDcirDataFile, 'file')
    fprintf('Loading: %s\n', rptDcirDataFile);
    load(rptDcirDataFile);
    if exist('dcir_soc_data', 'var')
        fprintf('RPT-DCIR data loaded successfully\n');
        rptDataAvailable = true;
    else
        fprintf('WARNING: dcir_soc_data variable not found in file\n');
        fprintf('Please check the variable name in the RPT-DCIR data file\n');
        rptDataAvailable = false;
    end
else
    if isempty(rptDcirDataFile)
        fprintf('INFO: RPT-DCIR data file path not specified\n');
        fprintf('RPT-DCIR columns will be NaN\n');
    else
        fprintf('WARNING: RPT-DCIR data file not found: %s\n', rptDcirDataFile);
        fprintf('RPT-DCIR columns will be NaN\n');
    end
    rptDataAvailable = false;
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

%% RPT SOC bins: 10% intervals (0, 10, 20, ..., 100)
rptSocBins = 0:10:100;  % [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

%% Initialize table data
fprintf('\n=== Creating Summary Table ===\n');
% Use Cell Array for better performance (O(N) instead of O(N^2))
dataList = {};     % Cell Array to store structs
rowCounter = 0;    % Counter for rows

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
        capacity_C3 = NaN;
        capacity_OCV = NaN;
        if capacityDataAvailable && isfield(allChannelsCapacity, chName)
            channelCapacity = allChannelsCapacity.(chName);
            capacityValues_C3 = channelCapacity.capacity;  % C/3 용량
            capacityValues_OCV = channelCapacity.capacity_ocv;  % OCV 용량
            capacityCycleNumbers = channelCapacity.cycles;
            capIdx = find(capacityCycleNumbers == cycleNum, 1);
            if ~isempty(capIdx)
                capacity_C3 = capacityValues_C3(capIdx);
                capacity_OCV = capacityValues_OCV(capIdx);
            else
                fprintf('  WARNING: No capacity data for %s cycle %d\n', chName, cycleNum);
            end
        else
            if capacityDataAvailable
                fprintf('  WARNING: Channel %s not found in capacity data\n', chName);
            end
        end
        
        % Get RPT-DCIR data for this channel and cycle
        rptData = struct();
        if rptDataAvailable && isfield(dcir_soc_data, chName)
            % Cycle name mapping: 0 -> cyc0, 200 -> cyc200, etc.
            cycleFieldName = sprintf('cyc%d', cycleNum);
            if isfield(dcir_soc_data.(chName), cycleFieldName)
                cycleRptData = dcir_soc_data.(chName).(cycleFieldName);
                
                % Extract charge table data
                if isfield(cycleRptData, 'charge_table')
                    rptData.charge_table = cycleRptData.charge_table;
                end
                
                % Extract discharge table data
                if isfield(cycleRptData, 'discharge_table')
                    rptData.discharge_table = cycleRptData.discharge_table;
                end
            end
        end
        
        % Check if this channel exists in cycle data
        channelFieldName = sprintf('ch%d', chNum);
        chgStructName = sprintf('%s_Drive_%s_ChgEvent', channelFieldName, cycleType);
        dchgStructName = sprintf('%s_Drive_%s_DchgEvent', channelFieldName, cycleType);  % Fixed: DchEvent -> DchgEvent (match EventAnalysis)
        
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
                    
                    % Process each event individually (not averaging)
                    % Get list of events for this DC profile and event type
                    events = {};
                    if strcmp(eventTypeStr, 'Charge')
                        if isfield(cycleData, chgStructName) && isfield(cycleData.(chgStructName), socLevel) && ...
                           isfield(cycleData.(chgStructName).(socLevel), profile)
                            events = fieldnames(cycleData.(chgStructName).(socLevel).(profile));
                        end
                    else % Discharge
                        if isfield(cycleData, dchgStructName) && isfield(cycleData.(dchgStructName), socLevel) && ...
                           isfield(cycleData.(dchgStructName).(socLevel), profile)
                            events = fieldnames(cycleData.(dchgStructName).(socLevel).(profile));
                        end
                    end
                    
                    % Process each event individually
                    for evt_idx = 1:length(events)
                        event = events{evt_idx};
                        
                        % Get the appropriate structure
                        if strcmp(eventTypeStr, 'Charge')
                            eventStruct = cycleData.(chgStructName).(socLevel).(profile).(event);
                        else
                            eventStruct = cycleData.(dchgStructName).(socLevel).(profile).(event);
                        end
                        
                        % Check if this event has at least one valid resistance value
                        hasData = false;
                        for dt_idx = 1:length(dt_list)
                            dt_sec = dt_list(dt_idx);
                            fieldName = sprintf('DCIR_%ds', dt_sec);
                            if isfield(eventStruct, fieldName) && isfield(eventStruct.(fieldName), 'val')
                                val = eventStruct.(fieldName).val;
                                if ~isnan(val)
                                    hasData = true;
                                    break;
                                end
                            end
                        end
                        
                        % Only create row if event has at least one valid resistance value
                        if hasData
                            rowData = struct();
                            rowData.Cycle = cycleNum;
                            rowData.Channel = chNum;
                            rowData.DC_Profile = profile;
                            rowData.EventType = eventTypeStr;
                            
                            % Extract SOC value from socLevel (e.g., 'SOC50' -> 50)
                            socValue = NaN;
                            if ~isempty(socLevel)
                                socMatch = regexp(socLevel, 'SOC(\d+)', 'tokens');
                                if ~isempty(socMatch)
                                    socValue = str2double(socMatch{1}{1});
                                end
                            end
                            rowData.SOC = socValue;  % SOC value as numeric (50, 70, 90)
                            
                            rowData.Capacity_C3 = capacity_C3;
                            rowData.Capacity_OCV = capacity_OCV;
                            
                            % Add Drive Cycle DCIR values (individual event values, not mean)
                            for dt_idx = 1:length(dt_list)
                                dt_sec = dt_list(dt_idx);
                                fieldName = sprintf('DCIR_%ds', dt_sec);
                                rFieldName = sprintf('R_%ds', dt_sec);
                                if isfield(eventStruct, fieldName) && isfield(eventStruct.(fieldName), 'val')
                                    val = eventStruct.(fieldName).val;
                                    if ~isnan(val)
                                        rowData.(rFieldName) = val;
                                    else
                                        rowData.(rFieldName) = NaN;
                                    end
                                else
                                    rowData.(rFieldName) = NaN;
                                end
                            end
                            
                            % [NEW] Automatically add all Feature fields
                            evtFields = fieldnames(eventStruct);
                            for f_idx = 1:length(evtFields)
                                fname = evtFields{f_idx};
                                if startsWith(fname, 'Feat_')
                                    rowData.(fname) = eventStruct.(fname);
                                end
                            end
                        
                            % Add RPT-DCIR data
                            % Initialize all RPT columns (10% SOC bins) with NaN first
                            for socBin = rptSocBins
                                socStr = sprintf('SOC%d', socBin);
                                for dt_idx = 1:length(dt_list)
                                    dt_sec = dt_list(dt_idx);
                                    rowData.(sprintf('RPT_Chg_%s_R%d', socStr, dt_sec)) = NaN;
                                    rowData.(sprintf('RPT_Dch_%s_R%d', socStr, dt_sec)) = NaN;
                                end
                            end
                            
                            % Extract SOC levels and DCIR values from RPT data
                            % Group by 10% SOC bins
                            if ~isempty(fieldnames(rptData))
                                % Process charge table
                                if isfield(rptData, 'charge_table') && ~isempty(rptData.charge_table)
                                    chargeTable = rptData.charge_table;
                                    if istable(chargeTable) && height(chargeTable) > 0 && ismember('SOC', chargeTable.Properties.VariableNames)
                                    % For each SOC bin (0, 10, 20, ..., 100)
                                    for binIdx = 1:length(rptSocBins)
                                        socBin = rptSocBins(binIdx);
                                        socStr = sprintf('SOC%d', socBin);
                                        
                                        % Determine SOC range for this bin
                                        if socBin == 0
                                            socRange = [0, 9.99];  % 0-9.99%
                                        elseif socBin == 100
                                            socRange = [100, 100];  % Exactly 100%
                                        else
                                            socRange = [socBin, socBin + 9.99];  % e.g., 10-19.99%
                                        end
                                        
                                        % Find rows within this SOC range
                                        socRows = chargeTable.SOC >= socRange(1) & chargeTable.SOC <= socRange(2);
                                        
                                        if any(socRows)  % If there are data points in this range
                                            % Add RPT charge DCIR values for each time interval
                                            for dt_idx = 1:length(dt_list)
                                                dt_sec = dt_list(dt_idx);
                                                rptFieldName = sprintf('RPT_Chg_%s_R%d', socStr, dt_sec);
                                                rColName = sprintf('R%d_mOhm', dt_sec);
                                                
                                                if ismember(rColName, chargeTable.Properties.VariableNames)
                                                    rValues = chargeTable.(rColName)(socRows);
                                                    rValues = rValues(~isnan(rValues));
                                                    if ~isempty(rValues)
                                                        rowData.(rptFieldName) = mean(rValues);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            
                                % Process discharge table
                                if isfield(rptData, 'discharge_table') && ~isempty(rptData.discharge_table)
                                    dischargeTable = rptData.discharge_table;
                                    if istable(dischargeTable) && height(dischargeTable) > 0 && ismember('SOC', dischargeTable.Properties.VariableNames)
                                    % For each SOC bin (0, 10, 20, ..., 100)
                                    for binIdx = 1:length(rptSocBins)
                                        socBin = rptSocBins(binIdx);
                                        socStr = sprintf('SOC%d', socBin);
                                        
                                        % Determine SOC range for this bin
                                        if socBin == 0
                                            socRange = [0, 9.99];  % 0-9.99%
                                        elseif socBin == 100
                                            socRange = [100, 100];  % Exactly 100%
                                        else
                                            socRange = [socBin, socBin + 9.99];  % e.g., 10-19.99%
                                        end
                                        
                                        % Find rows within this SOC range
                                        socRows = dischargeTable.SOC >= socRange(1) & dischargeTable.SOC <= socRange(2);
                                        
                                        if any(socRows)  % If there are data points in this range
                                            % Add RPT discharge DCIR values for each time interval
                                            for dt_idx = 1:length(dt_list)
                                                dt_sec = dt_list(dt_idx);
                                                rptFieldName = sprintf('RPT_Dch_%s_R%d', socStr, dt_sec);
                                                rColName = sprintf('R%d_mOhm', dt_sec);
                                                
                                                if ismember(rColName, dischargeTable.Properties.VariableNames)
                                                    rValues = dischargeTable.(rColName)(socRows);
                                                    rValues = rValues(~isnan(rValues));
                                                    if ~isempty(rValues)
                                                        rowData.(rptFieldName) = mean(rValues);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            end
                            
                            % Add to data list (Cell Array - O(1) operation, very fast)
                            rowCounter = rowCounter + 1;
                            dataList{rowCounter, 1} = rowData;
                        end
                    end
                end
            end
        end
    end
end

%% Convert to table (루프 종료 후 한 번에 변환)
fprintf('\n=== Finalizing Data Structure (Merging %d rows)... ===\n', rowCounter);

if ~isempty(dataList)
    % 1. Collect all field names from all rows
    allFields = {};
    for i = 1:length(dataList)
        allFields = union(allFields, fieldnames(dataList{i}));
    end
    fprintf('  Total unique fields: %d\n', length(allFields));
    
    % 2. Ensure all structs have the same fields (fill missing fields with NaN)
    fprintf('  Unifying structure fields...\n');
    for i = 1:length(dataList)
        currentFields = fieldnames(dataList{i});
        missingFields = setdiff(allFields, currentFields);
        for f = 1:length(missingFields)
            dataList{i}.(missingFields{f}) = NaN;
        end
    end
    
    % 3. Convert Cell Array to Struct Array (very fast - single operation)
    fprintf('  Converting to structure array...\n');
    dataStructArray = [dataList{:}];
    
    % 4. Convert to table
    fprintf('  Converting to table...\n');
    summaryTable = struct2table(dataStructArray);
    
    % Reorder columns: Cycle, Channel, DC_Profile, EventType, SOC, Capacity_C3, Capacity_OCV, 
    % Drive Cycle DCIR (R_1s, R_3s, ...), RPT-DCIR columns
    colOrder = {'Cycle', 'Channel', 'DC_Profile', 'EventType', 'SOC', 'Capacity_C3', 'Capacity_OCV'};
    
    % Add Drive Cycle DCIR columns
    for dt_idx = 1:length(dt_list)
        dt_sec = dt_list(dt_idx);
        colOrder{end+1} = sprintf('R_%ds', dt_sec);
    end
    
    % Add RPT-DCIR columns (if they exist)
    % Sort by SOC bin (0, 10, 20, ..., 100), then by time interval
    allCols = summaryTable.Properties.VariableNames;
    rptCols = allCols(startsWith(allCols, 'RPT_'));
    if ~isempty(rptCols)
        % Separate charge and discharge
        chgCols = rptCols(contains(rptCols, 'RPT_Chg_'));
        dchCols = rptCols(contains(rptCols, 'RPT_Dch_'));
        
        % Sort charge columns: by SOC bin (0, 10, 20, ..., 100), then by time (1, 3, 5, 10, 30, 60)
        chgColsSorted = {};
        for socBin = rptSocBins
            socStr = sprintf('SOC%d', socBin);
            for dt_sec = dt_list
                colName = sprintf('RPT_Chg_%s_R%d', socStr, dt_sec);
                if ismember(colName, chgCols)
                    chgColsSorted{end+1} = colName;
                end
            end
        end
        
        % Sort discharge columns: by SOC bin (0, 10, 20, ..., 100), then by time (1, 3, 5, 10, 30, 60)
        dchColsSorted = {};
        for socBin = rptSocBins
            socStr = sprintf('SOC%d', socBin);
            for dt_sec = dt_list
                colName = sprintf('RPT_Dch_%s_R%d', socStr, dt_sec);
                if ismember(colName, dchCols)
                    dchColsSorted{end+1} = colName;
                end
            end
        end
        
        colOrder = [colOrder, chgColsSorted, dchColsSorted];
    end
    
    % [FIX] Feature 컬럼들을 colOrder에 추가 (이 부분이 빠져서 삭제되었던 것임)
    allTableCols = summaryTable.Properties.VariableNames;
    featCols = allTableCols(startsWith(allTableCols, 'Feat_')); % Feat_로 시작하는 모든 컬럼 찾기
    colOrder = [colOrder, featCols]; % 기존 순서 뒤에 Feature 추가
    
    % Only include columns that exist in the table
    colOrder = colOrder(ismember(colOrder, summaryTable.Properties.VariableNames));
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
        unique_event_types = unique(summaryTable.EventType);
        fprintf('  Event Types: %s\n', strjoin(unique_event_types, ', '));
        
        % Count by EventType
        fprintf('\n  Data count by EventType:\n');
        for i = 1:length(unique_event_types)
            et = unique_event_types(i);
            count = sum(strcmpi(summaryTable.EventType, et));
            fprintf('    %s: %d rows\n', string(et), count);
        end
        
        % Check if both Charge and Discharge exist
        has_charge = any(strcmpi(unique_event_types, 'Charge'));
        has_discharge = any(strcmpi(unique_event_types, 'Discharge'));
        
        if ~has_charge
            fprintf('\n  WARNING: No Charge events found in summary table!\n');
        end
        if ~has_discharge
            fprintf('\n  WARNING: No Discharge events found in summary table!\n');
            fprintf('  This may be because:\n');
            fprintf('    1. Discharge events were not detected in DriveCycle_EventAnalysis_01.m\n');
            fprintf('    2. Discharge events were filtered out by strict criteria\n');
            fprintf('    3. Discharge events exist but were not aggregated properly\n');
        end
    end
else
    fprintf('WARNING: No data found for summary table\n');
end



fprintf('\n=== Data Aggregation Complete ===\n');

