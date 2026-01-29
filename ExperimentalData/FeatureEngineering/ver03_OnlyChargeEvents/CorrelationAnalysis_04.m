%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 04_CorrelationAnalysis.m
% 목적: 조합별 시간별 저항값(Rchg/Rdchg)의 사이클 간 p-value 계산
% 
% 입력:
%   - Results/<paramLabel>/Lab_DC_Events_Features_*cyc.mat (02번 스크립트 출력)
%
% 출력:
%   - PValue_Summary_<type>.txt (각 조합별 p-value 요약)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== P-Value Analysis (Time-Resistance by Cycle) ===\n');

%% 설정
baseResultsDir = fullfile(pwd, 'Results');

% 이벤트 타입 선택: 'Charge', 'Discharge', 'Both'
event_type_selection = 'Both';  % 'Charge': 충전만, 'Discharge': 방전만, 'Both': 모두

% 이벤트 타입 선택 검증
if ~ismember(event_type_selection, {'Charge', 'Discharge', 'Both'})
    error('event_type_selection must be ''Charge'', ''Discharge'', or ''Both''');
end

fprintf('Base results directory: %s\n', baseResultsDir);
fprintf('Event type selection: %s\n', event_type_selection);

% 파라미터 조합 (Simple_EventDetection_01/FeatureExtraction_02와 동일)
min_duration_list = [10, 30, 60];
max_I_std_list = [2.0, 1.0, 0.5];
max_I_range_list = 0;

paramSets = struct('min_duration', {}, 'max_I_std', {}, 'max_I_range', {});
setIdx = 0;
for md = min_duration_list
    for ms = max_I_std_list
        for mr = max_I_range_list
            setIdx = setIdx + 1;
            paramSets(setIdx).min_duration = md;
            paramSets(setIdx).max_I_std = ms;
            paramSets(setIdx).max_I_range = mr;
        end
    end
end

%% 조합별 처리
for p = 1:length(paramSets)
    min_duration = paramSets(p).min_duration;
    max_I_std = paramSets(p).max_I_std;
    max_I_range = paramSets(p).max_I_range;
    paramLabel = sprintf('min%d_std%s_rng%s_%s', ...
        min_duration, num2str(max_I_std), num2str(max_I_range), event_type_selection);
    paramLabel = regexprep(paramLabel, '\.', 'p');
    paramLabel = regexprep(paramLabel, '[^a-zA-Z0-9_]', '_');
    
    inputDir = fullfile(baseResultsDir, paramLabel);
    outputDir = inputDir;
    
    fprintf('\n=== ParamSet: %s ===\n', paramLabel);
    fprintf('Input directory: %s\n', inputDir);
    fprintf('Output directory: %s\n', outputDir);
    
    %% 파일 찾기
    matFiles = dir(fullfile(inputDir, 'Lab_DC_Events_Features_*cyc.mat'));
    if isempty(matFiles)
        fprintf('  WARNING: No feature files found in %s\n', inputDir);
        continue;
    end
    
    fprintf('Found %d feature files\n', length(matFiles));
    
    %% 데이터 집계
    fprintf('\n=== Aggregating Features ===\n');
    
    % 집계된 데이터를 저장할 셀 배열
    aggregatedData = {};
    
    for i = 1:length(matFiles)
        fileName = matFiles(i).name;
        filePath = fullfile(inputDir, fileName);
    
    % 사이클 번호 추출
    token = regexp(fileName, 'Lab_DC_Events_Features_(\d+cyc)\.mat', 'tokens');
    if isempty(token)
        fprintf('  WARNING: Cannot parse cycle from %s. Skipping.\n', fileName);
        continue;
    end
    cycleStr = token{1}{1};
    cycleNum = str2double(regexp(cycleStr, '\d+', 'match', 'once'));
    
    fprintf('  Processing %s (Cycle %d)...\n', fileName, cycleNum);
    
    % 파일 로드
    load(filePath);
    
    % 변수명 찾기
    vars = who('-file', filePath);
    dataVarName = '';
    for v = 1:length(vars)
        if contains(vars{v}, 'Lab_DC_DCIR_') || contains(vars{v}, 'Lab_DC_Events_Features_')
            dataVarName = vars{v};
            break;
        end
    end
    
    if isempty(dataVarName)
        fprintf('    WARNING: Cannot find data variable. Skipping.\n');
        continue;
    end
    
    data = eval(dataVarName);
    
    % 채널 필드 필터링 (선택된 타입만)
    allChFields = fieldnames(data);
    if strcmp(event_type_selection, 'Charge')
        chFields = allChFields(contains(allChFields, '_Charge'));
    elseif strcmp(event_type_selection, 'Discharge')
        chFields = allChFields(contains(allChFields, '_Discharge'));
    else % 'Both'
        chFields = allChFields;
    end
    
    if isempty(chFields)
        fprintf('    WARNING: No channels found for type ''%s''. Skipping.\n', event_type_selection);
        continue;
    end
    
    % 각 채널 처리
    for chIdx = 1:length(chFields)
        chField = chFields{chIdx};
        
        % 채널 번호 추출
        chMatch = regexp(chField, '^(ch\d+)_', 'tokens', 'once');
        if isempty(chMatch)
            continue;
        end
        chNumStr = chMatch{1};
        chNum = str2double(regexp(chNumStr, '\d+', 'match', 'once'));
        
        % 이벤트 타입 추출
        if contains(chField, '_Charge')
            eventType = 'Charge';
        elseif contains(chField, '_Discharge')
            eventType = 'Discharge';
        else
            eventType = 'Unknown';
        end
        
        if ~isfield(data, chField) || ~isstruct(data.(chField))
            continue;
        end
        
        socs = fieldnames(data.(chField));
        
        for s = 1:length(socs)
            socName = socs{s};
            if ~isfield(data.(chField), socName) || ~isstruct(data.(chField).(socName))
                continue;
            end
            
            profs = fieldnames(data.(chField).(socName));
            
            for p = 1:length(profs)
                profName = profs{p};
                if ~isfield(data.(chField).(socName), profName) || ~isstruct(data.(chField).(socName).(profName))
                    continue;
                end
                
                events = fieldnames(data.(chField).(socName).(profName));
                
                for e = 1:length(events)
                    evtName = events{e};
                    if ~isfield(data.(chField).(socName).(profName), evtName)
                        continue;
                    end
                    
                    evtData = data.(chField).(socName).(profName).(evtName);
                    
                    % 이벤트 번호 추출
                    evtNumMatch = regexp(evtName, 'event(\d+)', 'tokens', 'once');
                    if isempty(evtNumMatch)
                        evtNum = NaN;
                    else
                        evtNum = str2double(evtNumMatch{1});
                    end
                    
                    % 피쳐 추출 (시간별 저항만)
                    rowData = struct();
                    rowData.Cycle = cycleNum;
                    rowData.Channel = chNum;
                    rowData.SOC = socName;
                    rowData.Profile = profName;
                    rowData.EventType = eventType;
                    rowData.EventNumber = evtNum;
                    
                    % Rchg/Rdchg 필드 추출 (시간별 저항) - 모든 필드를 항상 추가 (필드 일치를 위해)
                    timePoints = [1, 3, 5, 10, 30, 60];
                    for t = 1:length(timePoints)
                        rchgField = sprintf('Rchg_%ds', timePoints(t));
                        rdchgField = sprintf('Rdchg_%ds', timePoints(t));
                        
                        % Charge 이벤트는 Rchg_*만, Discharge는 Rdchg_*만, 나머지는 NaN
                        if strcmp(eventType, 'Charge')
                            if isfield(evtData, rchgField)
                                rowData.(rchgField) = evtData.(rchgField);
                            else
                                rowData.(rchgField) = NaN;
                            end
                            rowData.(rdchgField) = NaN; % Discharge 필드는 항상 NaN
                        elseif strcmp(eventType, 'Discharge')
                            rowData.(rchgField) = NaN; % Charge 필드는 항상 NaN
                            if isfield(evtData, rdchgField)
                                rowData.(rdchgField) = evtData.(rdchgField);
                            else
                                rowData.(rdchgField) = NaN;
                            end
                        else % 'Both' 또는 기타
                            if isfield(evtData, rchgField)
                                rowData.(rchgField) = evtData.(rchgField);
                            else
                                rowData.(rchgField) = NaN;
                            end
                            if isfield(evtData, rdchgField)
                                rowData.(rdchgField) = evtData.(rdchgField);
                            else
                                rowData.(rdchgField) = NaN;
                            end
                        end
                    end
                    
                    % 셀 배열에 추가
                    aggregatedData{end+1} = rowData;
                end
            end
        end
    end
end

fprintf('  Total events aggregated: %d\n', length(aggregatedData));

%% 구조체 배열을 테이블로 변환
fprintf('\n=== Converting to Table ===\n');

if isempty(aggregatedData)
    error('No data aggregated. Please check if feature files contain the selected event type.');
end

% 구조체 배열 생성
dataStruct = [aggregatedData{:}];

% 테이블로 변환
summaryTable = struct2table(dataStruct);

fprintf('  Table created: %d rows, %d columns\n', height(summaryTable), width(summaryTable));
fprintf('  Columns: %s\n', strjoin(summaryTable.Properties.VariableNames, ', '));

% 기본 통계
fprintf('\n=== Data Summary ===\n');
fprintf('  Cycles: %s\n', mat2str(unique(summaryTable.Cycle)'));
fprintf('  Channels: %s\n', mat2str(unique(summaryTable.Channel)'));
fprintf('  Event Types: %s\n', strjoin(unique(summaryTable.EventType), ', '));
fprintf('  Profiles: %s\n', strjoin(unique(summaryTable.Profile), ', '));

%% AllChannels p-value (Charge/Discharge, SOC별)
fprintf('\n=== P-Value Summary (AllChannels) ===\n');
resistanceFields = {'Rchg_1s','Rchg_3s','Rchg_5s','Rchg_10s','Rchg_30s','Rchg_60s', ...
                    'Rdchg_1s','Rdchg_3s','Rdchg_5s','Rdchg_10s','Rdchg_30s','Rdchg_60s'};
availableResFields = resistanceFields(ismember(resistanceFields, summaryTable.Properties.VariableNames));
pvalueResults = {};
if isempty(availableResFields)
    fprintf('  WARNING: No Rchg/Rdchg fields found in summary table.\n');
else
    socList = unique(summaryTable.SOC);
    eventTypes = unique(summaryTable.EventType);
    for s = 1:length(socList)
        socName = socList{s};
        for e = 1:length(eventTypes)
            evtType = eventTypes{e};
            if strcmp(evtType, 'Charge')
                fieldList = availableResFields(startsWith(availableResFields, 'Rchg_'));
            elseif strcmp(evtType, 'Discharge')
                fieldList = availableResFields(startsWith(availableResFields, 'Rdchg_'));
            else
                continue;
            end
            
            if isempty(fieldList)
                continue;
            end
            
            for f = 1:length(fieldList)
                featName = fieldList{f};
                mask = strcmp(summaryTable.SOC, socName) & strcmp(summaryTable.EventType, evtType);
                mask = mask & ~isnan(summaryTable.Cycle) & ~isnan(summaryTable.(featName));
                vals = summaryTable.(featName)(mask);
                groups = summaryTable.Cycle(mask);
                
                pVal = NaN;
                if numel(vals) >= 3 && numel(unique(groups)) >= 2
                    pVal = kruskalwallis(vals, groups, 'off');
                end
                
                fprintf('[PVAL][%s][SOC=%s][AllChannels] %s: %s\n', ...
                    evtType, socName, featName, pvalue_str(pVal));
                
                % p-value 결과 저장
                pvalueResults{end+1} = struct('EventType', evtType, 'SOC', socName, ...
                    'ResistanceField', featName, 'PValue', pVal);
            end
        end
    end
    
    % p-value 결과를 테이블로 변환하여 저장
    if ~isempty(pvalueResults)
        pvalueTable = struct2table([pvalueResults{:}]);
        pvaluePath = fullfile(outputDir, sprintf('PValue_Summary_%s.txt', event_type_selection));
        writetable(pvalueTable, pvaluePath, 'Delimiter', '\t');
        fprintf('\n  Saved p-value summary to: %s\n', pvaluePath);
    end
end

    fprintf('\n=== P-Value Analysis Complete ===\n');
    fprintf('P-value summary saved to: %s\n', outputDir);
end

%% 로컬 함수
function outStr = pvalue_str(pVal)
    if isnan(pVal)
        outStr = 'NaN';
    else
        outStr = sprintf('%.4g', pVal);
    end
end
