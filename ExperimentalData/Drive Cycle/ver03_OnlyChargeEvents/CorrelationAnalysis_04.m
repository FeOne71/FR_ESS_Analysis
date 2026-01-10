%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 04_CorrelationAnalysis.m
% 목적: 모든 채널/사이클의 피쳐 값을 한 테이블로 집계하고 상관관계 분석
% 
% 입력:
%   - Lab_DC_Events_Features_*cyc.mat (02번 스크립트 출력)
%
% 출력:
%   - 피쳐 집계 테이블 (Summary Table)
%   - 산점도: Cycle vs Feat_Var_dVdQ, Feat_Slope_V, DCIR 등
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Correlation Analysis (Feature Aggregation) ===\n');

%% 설정
inputDir = fullfile(pwd, 'Results');
outputDir = fullfile(pwd, 'Results');

% 이벤트 타입 선택: 'Charge', 'Discharge', 'Both'
event_type_selection = 'Discharge';  % 'Charge': 충전만, 'Discharge': 방전만, 'Both': 모두

% 이벤트 타입 선택 검증
if ~ismember(event_type_selection, {'Charge', 'Discharge', 'Both'})
    error('event_type_selection must be ''Charge'', ''Discharge'', or ''Both''');
end

fprintf('Input directory: %s\n', inputDir);
fprintf('Output directory: %s\n', outputDir);
fprintf('Event type selection: %s\n', event_type_selection);

%% 파일 찾기
matFiles = dir(fullfile(inputDir, 'Lab_DC_Events_Features_*cyc.mat'));
if isempty(matFiles)
    error('No feature files found in %s\nPlease run 02_FeatureExtraction.m first', inputDir);
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
                    
                    % 피쳐 추출
                    rowData = struct();
                    rowData.Cycle = cycleNum;
                    rowData.Channel = chNum;
                    rowData.SOC = socName;
                    rowData.Profile = profName;
                    rowData.EventType = eventType;
                    rowData.EventNumber = evtNum;
                    
                    % 주요 피쳐 추출
                    if isfield(evtData, 'Feat_Var_dVdQ')
                        rowData.Feat_Var_dVdQ = evtData.Feat_Var_dVdQ;
                    else
                        rowData.Feat_Var_dVdQ = NaN;
                    end
                    
                    if isfield(evtData, 'Feat_Slope_V')
                        rowData.Feat_Slope_V = evtData.Feat_Slope_V;
                    else
                        rowData.Feat_Slope_V = NaN;
                    end
                    
                    % DCIR 값들 추출
                    dcir_fields = {'DCIR_1s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s'};
                    for d = 1:length(dcir_fields)
                        dcir_field = dcir_fields{d};
                        if isfield(evtData, dcir_field) && isstruct(evtData.(dcir_field))
                            if isfield(evtData.(dcir_field), 'val')
                                rowData.(dcir_field) = evtData.(dcir_field).val;
                            else
                                rowData.(dcir_field) = NaN;
                            end
                        else
                            rowData.(dcir_field) = NaN;
                        end
                    end
                    
                    % 추가 피쳐들 (선택사항)
                    if isfield(evtData, 'Feat_dQ')
                        rowData.Feat_dQ = evtData.Feat_dQ;
                    else
                        rowData.Feat_dQ = NaN;
                    end
                    
                    if isfield(evtData, 'Feat_Mean_dVdQ')
                        rowData.Feat_Mean_dVdQ = evtData.Feat_Mean_dVdQ;
                    else
                        rowData.Feat_Mean_dVdQ = NaN;
                    end
                    
                    if isfield(evtData, 'duration')
                        rowData.Duration = evtData.duration;
                    else
                        rowData.Duration = NaN;
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

%% 테이블 저장
savePath = fullfile(outputDir, sprintf('Feature_Summary_Table_%s.mat', event_type_selection));
save(savePath, 'summaryTable', 'event_type_selection');
fprintf('\n  Saved summary table to: %s\n', savePath);

%% 산점도 생성
fprintf('\n=== Creating Scatter Plots ===\n');

figuresDir = fullfile(outputDir, 'figures', 'CorrelationAnalysis');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

% 주요 피쳐 목록 (DCIR과 일반 피쳐 구분)
dcirFeatures = {'DCIR_1s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s'};
dcirLabels = {'DCIR 1s (mΩ)', 'DCIR 5s (mΩ)', 'DCIR 10s (mΩ)', 'DCIR 30s (mΩ)'};
generalFeatures = {'Feat_Var_dVdQ', 'Feat_Slope_V'};
generalLabels = {'Var dV/dQ (V/Ah)²', 'Voltage Slope (V/s)'};

% 사용 가능한 피쳐 확인
availableDcirFeatures = {};
availableDcirLabels = {};
for f = 1:length(dcirFeatures)
    if ismember(dcirFeatures{f}, summaryTable.Properties.VariableNames)
        availableDcirFeatures{end+1} = dcirFeatures{f};
        availableDcirLabels{end+1} = dcirLabels{f};
    end
end

availableGeneralFeatures = {};
availableGeneralLabels = {};
for f = 1:length(generalFeatures)
    if ismember(generalFeatures{f}, summaryTable.Properties.VariableNames)
        availableGeneralFeatures{end+1} = generalFeatures{f};
        availableGeneralLabels{end+1} = generalLabels{f};
    end
end

fprintf('  Available DCIR features: %s\n', strjoin(availableDcirLabels, ', '));
fprintf('  Available general features: %s\n', strjoin(availableGeneralLabels, ', '));

% 채널 목록
uniqueChannels = unique(summaryTable.Channel);
uniqueChannels = sort(uniqueChannels);
fprintf('  Channels: %s\n', mat2str(uniqueChannels'));

%% 일반 피쳐: Cycle vs 피쳐 (통합 그래프)
for f = 1:length(availableGeneralFeatures)
    featName = availableGeneralFeatures{f};
    featLabel = availableGeneralLabels{f};
    
    % 유효한 데이터만 선택
    validMask = ~isnan(summaryTable.Cycle) & ~isnan(summaryTable.(featName));
    cycleData = summaryTable.Cycle(validMask);
    featData = summaryTable.(featName)(validMask);
    
    if length(cycleData) < 3
        fprintf('    Skipping %s: insufficient data (n=%d)\n', featLabel, length(cycleData));
        continue;
    end
    
    % 산점도 생성
    fig = figure('Name', sprintf('Cycle vs %s', featLabel), ...
                 'Position', [100, 100, 1200, 800], 'Visible', 'on');
    
    % 채널별로 색상 구분
    channelData = summaryTable.Channel(validMask);
    channelColors = lines(length(uniqueChannels));
    
    hold on;
    for chIdx = 1:length(uniqueChannels)
        chNum = uniqueChannels(chIdx);
        chMask = channelData == chNum;
        if any(chMask)
            scatter(cycleData(chMask), featData(chMask), 80, ...
                    channelColors(chIdx, :), 'filled', 'MarkerFaceAlpha', 0.7, ...
                    'DisplayName', sprintf('Ch%d', chNum));
        end
    end
    
    % 선형 회귀선 (전체 데이터)
    if length(cycleData) >= 2
        p_fit = polyfit(cycleData, featData, 1);
        cycleRange = [min(cycleData), max(cycleData)];
        fitLine = polyval(p_fit, cycleRange);
        plot(cycleRange, fitLine, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Linear Fit (All)');
        
        % 상관계수 계산
        [R, P] = corrcoef(cycleData, featData);
        r_val = R(1, 2);
        p_val = P(1, 2);
        r_squared = r_val^2;
    else
        r_val = NaN;
        p_val = NaN;
        r_squared = NaN;
    end
    
    xlabel('Cycle', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel(featLabel, 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('Cycle vs %s (%s Events)\nn=%d', featLabel, event_type_selection, length(cycleData)), ...
          'FontSize', 15, 'FontWeight', 'bold');
    
    % 통계 정보 텍스트
    if ~isnan(r_val)
        statsText = sprintf('R = %.4f\nR² = %.4f\np = %.4e\nn = %d', ...
                          r_val, r_squared, p_val, length(cycleData));
        if p_val < 0.001
            sigText = '***';
        elseif p_val < 0.01
            sigText = '**';
        elseif p_val < 0.05
            sigText = '*';
        else
            sigText = 'ns';
        end
        statsText = sprintf('%s\nSig: %s', statsText, sigText);
        
        text(0.98, 0.98, statsText, ...
             'Units', 'normalized', ...
             'HorizontalAlignment', 'right', ...
             'VerticalAlignment', 'top', ...
             'FontSize', 11, ...
             'BackgroundColor', 'white', ...
             'EdgeColor', 'black', ...
             'LineWidth', 1.5);
    end
    
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    hold off;
    
    % 저장 (fig만)
    saveFileName = sprintf('Scatter_Cycle_vs_%s_%s', ...
                           strrep(featName, 'Feat_', ''), event_type_selection);
    saveFileName = regexprep(saveFileName, '[^a-zA-Z0-9_]', '_');
    savePath_fig = fullfile(figuresDir, [saveFileName, '.fig']);
    
    saveas(fig, savePath_fig);
    
    fprintf('    Saved: %s (n=%d, R=%.4f, R²=%.4f, p=%.4e)\n', ...
            featLabel, length(cycleData), r_val, r_squared, p_val);
    
    close(fig);
end

%% DCIR 피쳐: 채널별로 별도 그래프 생성
for f = 1:length(availableDcirFeatures)
    featName = availableDcirFeatures{f};
    featLabel = availableDcirLabels{f};
    
    % 각 채널별로 그래프 생성
    for chIdx = 1:length(uniqueChannels)
        chNum = uniqueChannels(chIdx);
        
        % 해당 채널의 데이터만 선택
        chMask = summaryTable.Channel == chNum;
        validMask = chMask & ~isnan(summaryTable.Cycle) & ~isnan(summaryTable.(featName));
        cycleData = summaryTable.Cycle(validMask);
        featData = summaryTable.(featName)(validMask);
        
        if length(cycleData) < 2
            continue;  % 데이터가 부족한 채널은 스킵
        end
        
        % 산점도 생성
        fig = figure('Name', sprintf('Cycle vs %s (Ch%d)', featLabel, chNum), ...
                     'Position', [100 + chIdx*50, 100 + f*50, 1000, 700], 'Visible', 'on');
        
        hold on;
        
        % 스캐터 플롯
        scatter(cycleData, featData, 100, 'filled', 'MarkerFaceAlpha', 0.7, ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1.5);
        
        % 선형 회귀선
        if length(cycleData) >= 2
            p_fit = polyfit(cycleData, featData, 1);
            cycleRange = [min(cycleData), max(cycleData)];
            fitLine = polyval(p_fit, cycleRange);
            plot(cycleRange, fitLine, 'r--', 'LineWidth', 2.5, 'DisplayName', 'Linear Fit');
            
            % 상관계수 계산
            [R, P] = corrcoef(cycleData, featData);
            r_val = R(1, 2);
            p_val = P(1, 2);
            r_squared = r_val^2;
        else
            r_val = NaN;
            p_val = NaN;
            r_squared = NaN;
        end
        
        xlabel('Cycle', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel(featLabel, 'FontSize', 13, 'FontWeight', 'bold');
        title(sprintf('Cycle vs %s - Channel %d (%s Events)\nn=%d', ...
              featLabel, chNum, event_type_selection, length(cycleData)), ...
              'FontSize', 15, 'FontWeight', 'bold');
        
        % 통계 정보 텍스트
        if ~isnan(r_val)
            statsText = sprintf('R = %.4f\nR² = %.4f\np = %.4e\nn = %d', ...
                              r_val, r_squared, p_val, length(cycleData));
            if p_val < 0.001
                sigText = '***';
            elseif p_val < 0.01
                sigText = '**';
            elseif p_val < 0.05
                sigText = '*';
            else
                sigText = 'ns';
            end
            statsText = sprintf('%s\nSig: %s', statsText, sigText);
            
            text(0.98, 0.98, statsText, ...
                 'Units', 'normalized', ...
                 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'top', ...
                 'FontSize', 11, ...
                 'BackgroundColor', 'white', ...
                 'EdgeColor', 'black', ...
                 'LineWidth', 1.5);
        end
        
        grid on;
        hold off;
        
        % 저장 (fig만)
        saveFileName = sprintf('Scatter_Cycle_vs_%s_Ch%d_%s', ...
                               featName, chNum, event_type_selection);
        saveFileName = regexprep(saveFileName, '[^a-zA-Z0-9_]', '_');
        savePath_fig = fullfile(figuresDir, [saveFileName, '.fig']);
        
        saveas(fig, savePath_fig);
        
        fprintf('    Saved: %s - Ch%d (n=%d, R=%.4f, R²=%.4f, p=%.4e)\n', ...
                featLabel, chNum, length(cycleData), r_val, r_squared, p_val);
        
        close(fig);
    end
end

fprintf('\n=== Correlation Analysis Complete ===\n');
fprintf('Summary table and scatter plots saved to: %s\n', outputDir);
