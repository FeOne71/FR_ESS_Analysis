%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 02_FieldData_DataAggregation.m
% 데이터 취합 및 Summary Table 생성 (필드 데이터용)
% 
% 목적: 
% - 01번 스크립트가 생성한 FieldData_*_Events.mat 파일들을 로드
% - SOH 매핑 추가 (Year → SOH: 2021→95%, 2022→90%, 2023→87%, 2024→85%, 2025→78%)
% - 모든 데이터를 하나의 Table로 취합
% - 컬럼: Year, Rack, EventType, Group, SOH, R_chg_1s, R_chg_3s, R_chg_5s, R_chg_10s, R_chg_30s, R_chg_60s
%
% 입력:
% - FieldData_*_Events.mat (01번 스크립트 출력)
%
% 출력:
% - FieldData_Summary_Table.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

fprintf('=== Field Data Aggregation ===\n');

%% Configuration - User Settings
% =========================================================================
% 스크립트 위치를 기준으로 상대 경로 설정
scriptDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(scriptDir, 'EventsResults');
inputDir = fullfile(resultsDir, 'FieldData_Analysis');
outputDir = fullfile(resultsDir, 'FieldData_Analysis');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% SOH 매핑 (Year → SOH %)
soh_map = containers.Map({'Y2021', 'Y2022', 'Y2023', 'Y2024', 'Y2025'}, ...
                         {95, 90, 87, 85, 78});

% 필터링: 모든 연도에 데이터가 있는 그룹만 사용
% NOTE: 02 스크립트에서는 필터링하지 않고 모든 그룹을 포함시킵니다.
% 필터링은 03, 04 스크립트에서 수행됩니다.
require_all_years = false;  % false: 모든 그룹 포함 (필터링은 03, 04에서 수행)
% =========================================================================

%% Auto-detect available event files
fprintf('\n=== Auto-detecting available event files ===\n');
matFiles = dir(fullfile(inputDir, 'FieldData_*_Events.mat'));
availableFiles = {};

for i = 1:length(matFiles)
    filename = matFiles(i).name;
    availableFiles{end+1} = filename;
end

if isempty(availableFiles)
    error('No FieldData_*_Events.mat files found in %s', inputDir);
end

fprintf('Found %d event file(s):\n', length(availableFiles));
for i = 1:length(availableFiles)
    fprintf('  %s\n', availableFiles{i});
end

%% Initialize table data
fprintf('\n=== Creating Summary Table ===\n');
tableData = [];

%% Process each event file
for fileIdx = 1:length(availableFiles)
    filename = availableFiles{fileIdx};
    filePath = fullfile(inputDir, filename);
    
    fprintf('\n--- Processing %s ---\n', filename);
    
    % 파일명에서 Rack과 EventType 추출
    % 예: FieldData_Rack01_charge_Events.mat
    match = regexp(filename, 'FieldData_(Rack\d+)_(\w+)_Events\.mat', 'tokens');
    if isempty(match)
        fprintf('WARNING: Cannot parse filename format: %s\n', filename);
        continue;
    end
    
    rackName = match{1}{1};
    eventType = match{1}{2};
    
    fprintf('  Rack: %s, EventType: %s\n', rackName, eventType);
    
    % 데이터 로드
    try
        load(filePath);
    catch ME
        fprintf('ERROR: Failed to load file: %s\n', filePath);
        fprintf('Error message: %s\n', ME.message);
        continue;
    end
    
    % 변수명 추출 (FieldData_Rack01_charge 형식)
    varName = sprintf('FieldData_%s_%s', rackName, eventType);
    if ~exist(varName, 'var')
        fprintf('WARNING: Variable %s not found in file\n', varName);
        continue;
    end
    
    eval(sprintf('eventData = %s;', varName));
    
    % 그룹별로 데이터 추출
    groupNames = fieldnames(eventData);
    fprintf('  Found %d groups\n', length(groupNames));
    
    % 사용 가능한 모든 연도 확인 (SOH 맵에서)
    all_required_years = keys(soh_map);
    all_required_years = sort(all_required_years);
    
    for g_idx = 1:length(groupNames)
        groupName_safe = groupNames{g_idx};
        % Convert safe field name back to original group name (for display)
        % Reverse of 01 script: _N -> -, and _ between numbers -> .
        % Example: V3_77_I12 -> V3.77 I12 (SOC excluded from grouping)
        groupName = strrep(groupName_safe, '_N', '-');  % _N을 -로 변환
        % Replace _ between digits with . (e.g., V3_77 -> V3.77)
        % Use regex to find pattern: digit_ digit and replace with digit. digit
        groupName = regexprep(groupName, '(\d)_(\d)', '$1.$2');
        groupData = eventData.(groupName_safe);
        
        % 연도별로 데이터 추출
        yearNames = fieldnames(groupData);
        yearNames = yearNames(startsWith(yearNames, 'Y'));
        
        % 모든 연도에 데이터가 있는지 확인
        if require_all_years
            has_all_years = true;
            missing_years = {};
            for y_idx = 1:length(all_required_years)
                req_year = all_required_years{y_idx};
                if ~any(strcmp(yearNames, req_year))
                    has_all_years = false;
                    missing_years{end+1} = req_year;
                end
            end
            
            if ~has_all_years
                fprintf('    Skipping group %s: missing years %s\n', groupName, strjoin(missing_years, ', '));
                continue;
            end
        end
        
        for y_idx = 1:length(yearNames)
            yearName = yearNames{y_idx};
            yearData = groupData.(yearName);
            
            % SOH 가져오기
            soh_value = NaN;
            if isKey(soh_map, yearName)
                soh_value = soh_map(yearName);
            end
            
            % 이벤트 필드 찾기
            eventFields = fieldnames(yearData);
            eventIndices = startsWith(eventFields, 'event');
            eventFields = eventFields(eventIndices);
            
            % 각 이벤트별로 행 추가
            for e_idx = 1:length(eventFields)
                eventName = eventFields{e_idx};
                eventStruct = yearData.(eventName);
                
                % Rchg 값 추출
                R_chg_1s = NaN;
                R_chg_3s = NaN;
                R_chg_5s = NaN;
                R_chg_10s = NaN;
                R_chg_30s = NaN;
                R_chg_60s = NaN;
                
                if isfield(eventStruct, 'R_chg_1s')
                    R_chg_1s = eventStruct.R_chg_1s;
                end
                if isfield(eventStruct, 'R_chg_3s')
                    R_chg_3s = eventStruct.R_chg_3s;
                end
                if isfield(eventStruct, 'R_chg_5s')
                    R_chg_5s = eventStruct.R_chg_5s;
                end
                if isfield(eventStruct, 'R_chg_10s')
                    R_chg_10s = eventStruct.R_chg_10s;
                end
                if isfield(eventStruct, 'R_chg_30s')
                    R_chg_30s = eventStruct.R_chg_30s;
                end
                if isfield(eventStruct, 'R_chg_60s')
                    R_chg_60s = eventStruct.R_chg_60s;
                end
                
                % 연도 숫자 추출 (Y2021 → 2021)
                yearNum = str2double(yearName(2:end));
                
                % 행 데이터 생성
                rowData = struct();
                rowData.Year = yearNum;
                rowData.Rack = rackName;
                rowData.EventType = eventType;
                rowData.Group = groupName;
                rowData.SOH = soh_value;
                rowData.R_chg_1s = R_chg_1s;
                rowData.R_chg_3s = R_chg_3s;
                rowData.R_chg_5s = R_chg_5s;
                rowData.R_chg_10s = R_chg_10s;
                rowData.R_chg_30s = R_chg_30s;
                rowData.R_chg_60s = R_chg_60s;
                
                % 타임스탬프 저장 (선택적)
                if isfield(eventStruct, 'timestamp')
                    rowData.Timestamp = eventStruct.timestamp;
                else
                    rowData.Timestamp = NaT;
                end
                
                % 테이블에 추가
                tableData = [tableData; rowData];
            end
        end
    end
end

%% Convert to table
if ~isempty(tableData)
    summaryTable = struct2table(tableData);
    
    % 컬럼 순서 정렬
    colOrder = {'Year', 'Rack', 'EventType', 'Group', 'SOH', 'R_chg_1s', 'R_chg_3s', 'R_chg_5s', 'R_chg_10s', 'R_chg_30s', 'R_chg_60s', 'Timestamp'};
    existingCols = ismember(colOrder, summaryTable.Properties.VariableNames);
    colOrder = colOrder(existingCols);
    summaryTable = summaryTable(:, colOrder);
    
    % Save table
    savePath = fullfile(outputDir, 'FieldData_Summary_Table.mat');
    save(savePath, 'summaryTable');
    fprintf('\n=== Summary Table Created ===\n');
    fprintf('Total rows: %d\n', height(summaryTable));
    fprintf('Saved to: %s\n', savePath);
    
    % Display first few rows
    fprintf('\nFirst 5 rows of summary table:\n');
    disp(summaryTable(1:min(5, height(summaryTable)), :));
    
    % Display summary statistics
    fprintf('\nSummary Statistics:\n');
    fprintf('  Years: %s\n', mat2str(unique(summaryTable.Year)'));
    fprintf('  Racks: %s\n', strjoin(unique(summaryTable.Rack), ', '));
    fprintf('  Event Types: %s\n', strjoin(unique(summaryTable.EventType), ', '));
    fprintf('  Groups: %d unique groups\n', length(unique(summaryTable.Group)));
    
    % 그룹별 이벤트 수
    fprintf('\nTop 10 groups by event count:\n');
    try
        groupCounts = groupsummary(summaryTable, 'Group', @(x) length(x));
    catch
        % Fallback: 직접 계산
        uniqueGroups = unique(summaryTable.Group);
        groupCounts = table(uniqueGroups, zeros(length(uniqueGroups), 1), ...
                           'VariableNames', {'Group', 'Count'});
        for i = 1:length(uniqueGroups)
            groupCounts.Count(i) = sum(strcmp(summaryTable.Group, uniqueGroups{i}));
        end
    end
    if istable(groupCounts)
        if ismember('GroupCount_Group', groupCounts.Properties.VariableNames)
            groupCounts = sortrows(groupCounts, 'GroupCount_Group', 'descend');
            for i = 1:min(10, height(groupCounts))
                fprintf('  %s: %d events\n', groupCounts.Group{i}, groupCounts.GroupCount_Group(i));
            end
        elseif ismember('Count', groupCounts.Properties.VariableNames)
            groupCounts = sortrows(groupCounts, 'Count', 'descend');
            for i = 1:min(10, height(groupCounts))
                fprintf('  %s: %d events\n', groupCounts.Group{i}, groupCounts.Count(i));
            end
        end
    end
    
    % SOH 분포 확인
    fprintf('\nSOH distribution:\n');
    sohStats = groupsummary(summaryTable, 'SOH', {'mean', 'std'}, {'R_chg_1s', 'R_chg_3s', 'R_chg_5s', 'R_chg_10s', 'R_chg_30s', 'R_chg_60s'});
    disp(sohStats);
    
else
    fprintf('WARNING: No data found for summary table\n');
end

fprintf('\n=== Data Aggregation Complete ===\n');

