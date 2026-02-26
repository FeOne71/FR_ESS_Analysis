%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FancyBoxplot_Resistance.m
% 목적: min10std2 조합의 충전/방전 시간별 저항을 fancy한 박스플롯으로 시각화
%
% 입력:
%   - Results/min10_std2_rng0_Both/Lab_DC_Events_Features_*cyc.mat
%
% 출력:
%   - Results/figures/FancyBoxplots/*.fig, *.png
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Fancy Boxplot Visualization for Resistance ===\n');

%% 설정
baseResultsDir = fullfile(pwd, 'Results');
paramLabel = 'min10_std2_rng0_Both';  % 또는 'min10std2' 등 실제 폴더명에 맞게 수정
inputDir = fullfile(baseResultsDir, paramLabel);

% daboxplot 경로 추가
daboxplotPath = 'C:\Users\Chulwon Jung\Downloads\frank-pk-DataViz-3.2.3.0\daboxplot';
if exist(daboxplotPath, 'dir')
    addpath(daboxplotPath);
    fprintf('daboxplot path added: %s\n', daboxplotPath);
else
    warning('daboxplot path not found. Using default boxplot.');
end

% 출력 폴더
figuresBaseDir = fullfile(baseResultsDir, 'figures', 'FancyBoxplots');
if ~exist(figuresBaseDir, 'dir')
    mkdir(figuresBaseDir);
end

% 시간 포인트 및 저항 필드
timePoints = [1, 3, 5, 10, 30, 60];
rFieldsCharge = arrayfun(@(t) sprintf('Rchg_%ds', t), timePoints, 'UniformOutput', false);
rFieldsDischarge = arrayfun(@(t) sprintf('Rdchg_%ds', t), timePoints, 'UniformOutput', false);

% 사이클 목록
cyclesWanted = {'0cyc','200cyc','400cyc','600cyc','800cyc','1000cyc'};

% 시각화 옵션
figVisible = 'on';  % 'on' 또는 'off'

fprintf('Input directory: %s\n', inputDir);

%% MAT 파일 찾기
matFiles = dir(fullfile(inputDir, 'Lab_DC_Events_Features_*cyc.mat'));
if isempty(matFiles)
    fprintf('ERROR: No feature files found in %s\n', inputDir);
    return;
end

fprintf('Found %d feature files\n', length(matFiles));

%% 데이터 집계
fprintf('\n=== Aggregating Data ===\n');

% 집계 구조체: aggCharge.(soc).(cycle).(rField), aggDischarge.(soc).(cycle).(rField)
aggCharge = struct();
aggDischarge = struct();

% 저항 차이 집계 구조체: aggChargeDiff.(soc).(cycle).(diffField)
% diffField: 'R3_R1', 'R5_R1', 'R10_R1', 'R30_R1', 'R60_R1'
aggChargeDiff = struct();
aggDischargeDiff = struct();

for i = 1:length(matFiles)
    fileName = matFiles(i).name;
    filePath = fullfile(inputDir, fileName);
    
    % 사이클 번호 추출
    token = regexp(fileName, 'Lab_DC_Events_Features_(\d+cyc)\.mat', 'tokens', 'once');
    if isempty(token)
        continue;
    end
    cycleName = token{1};
    cycleField = make_valid_cycle_field(cycleName);
    
    fprintf('  Loading %s...\n', fileName);
    
    % 데이터 로드
    dataStruct = load(filePath);
    vars = fieldnames(dataStruct);
    dataVarName = '';
    for v = 1:length(vars)
        if contains(vars{v}, 'Lab_DC_DCIR_') || contains(vars{v}, 'Lab_DC_Events_Features_')
            dataVarName = vars{v};
            break;
        end
    end
    if isempty(dataVarName)
        dataVarName = vars{1};
    end
    data = dataStruct.(dataVarName);
    
    % 채널 분류
    chFields = fieldnames(data);
    chargeCh = chFields(contains(chFields, '_Charge'));
    dischargeCh = chFields(contains(chFields, '_Discharge'));
    
    % 충전 이벤트 처리
    for c = 1:length(chargeCh)
        chName = chargeCh{c};
        if ~isfield(data, chName) || ~isstruct(data.(chName))
            continue;
        end
        socs = fieldnames(data.(chName));
        for s = 1:length(socs)
            socName = socs{s};
            if ~isfield(data.(chName), socName) || ~isstruct(data.(chName).(socName))
                continue;
            end
            profs = fieldnames(data.(chName).(socName));
            for pIdx = 1:length(profs)
                profName = profs{pIdx};
                if ~isfield(data.(chName).(socName), profName) || ~isstruct(data.(chName).(socName).(profName))
                    continue;
                end
                evts = fieldnames(data.(chName).(socName).(profName));
                for e = 1:length(evts)
                    evt = data.(chName).(socName).(profName).(evts{e});
                    if ~isstruct(evt)
                        continue;
                    end
                    
                    % 각 시간별 저항값 수집
                    rVals = struct();
                    hasR1 = false;
                    r1Val = [];
                    
                    for rIdx = 1:length(rFieldsCharge)
                        rField = rFieldsCharge{rIdx};
                        if ~isfield(evt, rField)
                            continue;
                        end
                        rVal = evt.(rField);
                        if ~isfinite(rVal) || rVal == 0
                            continue;
                        end
                        rVals.(rField) = rVal;
                        
                        % R1 값 저장
                        if strcmp(rField, 'Rchg_1s')
                            hasR1 = true;
                            r1Val = rVal;
                        end
                    end
                    
                    % 각 시간별 저항값 집계
                    for rIdx = 1:length(rFieldsCharge)
                        rField = rFieldsCharge{rIdx};
                        if isfield(rVals, rField)
                            if ~isfield(aggCharge, socName)
                                aggCharge.(socName) = struct();
                            end
                            if ~isfield(aggCharge.(socName), cycleField)
                                aggCharge.(socName).(cycleField) = struct();
                            end
                            if ~isfield(aggCharge.(socName).(cycleField), rField)
                                aggCharge.(socName).(cycleField).(rField) = [];
                            end
                            aggCharge.(socName).(cycleField).(rField) = ...
                                [aggCharge.(socName).(cycleField).(rField); rVals.(rField)];
                        end
                    end
                    
                    % 저항 차이 계산 (R3-R1, R5-R1, R10-R1, R30-R1, R60-R1)
                    if hasR1 && ~isempty(r1Val)
                        diffFields = {'Rchg_5s', 'Rchg_10s', 'Rchg_30s', 'Rchg_60s'};  % R3-R1 제외
                        diffFieldNames = {'R5_R1', 'R10_R1', 'R30_R1', 'R60_R1'};
                        
                        for dIdx = 1:length(diffFields)
                            rField = diffFields{dIdx};
                            diffFieldName = diffFieldNames{dIdx};
                            
                            if isfield(rVals, rField)
                                diffVal = rVals.(rField) - r1Val;
                                
                                if ~isfield(aggChargeDiff, socName)
                                    aggChargeDiff.(socName) = struct();
                                end
                                if ~isfield(aggChargeDiff.(socName), cycleField)
                                    aggChargeDiff.(socName).(cycleField) = struct();
                                end
                                if ~isfield(aggChargeDiff.(socName).(cycleField), diffFieldName)
                                    aggChargeDiff.(socName).(cycleField).(diffFieldName) = [];
                                end
                                aggChargeDiff.(socName).(cycleField).(diffFieldName) = ...
                                    [aggChargeDiff.(socName).(cycleField).(diffFieldName); diffVal];
                            end
                        end
                    end
                end
            end
        end
    end
    
    % 방전 이벤트 처리
    for c = 1:length(dischargeCh)
        chName = dischargeCh{c};
        if ~isfield(data, chName) || ~isstruct(data.(chName))
            continue;
        end
        socs = fieldnames(data.(chName));
        for s = 1:length(socs)
            socName = socs{s};
            if ~isfield(data.(chName), socName) || ~isstruct(data.(chName).(socName))
                continue;
            end
            profs = fieldnames(data.(chName).(socName));
            for pIdx = 1:length(profs)
                profName = profs{pIdx};
                if ~isfield(data.(chName).(socName), profName) || ~isstruct(data.(chName).(socName).(profName))
                    continue;
                end
                evts = fieldnames(data.(chName).(socName).(profName));
                for e = 1:length(evts)
                    evt = data.(chName).(socName).(profName).(evts{e});
                    if ~isstruct(evt)
                        continue;
                    end
                    
                    % 각 시간별 저항값 수집
                    rVals = struct();
                    hasR1 = false;
                    r1Val = [];
                    
                    for rIdx = 1:length(rFieldsDischarge)
                        rField = rFieldsDischarge{rIdx};
                        if ~isfield(evt, rField)
                            continue;
                        end
                        rVal = evt.(rField);
                        if ~isfinite(rVal) || rVal == 0
                            continue;
                        end
                        rVals.(rField) = rVal;
                        
                        % R1 값 저장
                        if strcmp(rField, 'Rdchg_1s')
                            hasR1 = true;
                            r1Val = rVal;
                        end
                    end
                    
                    % 각 시간별 저항값 집계
                    for rIdx = 1:length(rFieldsDischarge)
                        rField = rFieldsDischarge{rIdx};
                        if isfield(rVals, rField)
                            if ~isfield(aggDischarge, socName)
                                aggDischarge.(socName) = struct();
                            end
                            if ~isfield(aggDischarge.(socName), cycleField)
                                aggDischarge.(socName).(cycleField) = struct();
                            end
                            if ~isfield(aggDischarge.(socName).(cycleField), rField)
                                aggDischarge.(socName).(cycleField).(rField) = [];
                            end
                            aggDischarge.(socName).(cycleField).(rField) = ...
                                [aggDischarge.(socName).(cycleField).(rField); rVals.(rField)];
                        end
                    end
                    
                    % 저항 차이 계산 (R3-R1, R5-R1, R10-R1, R30-R1, R60-R1)
                    if hasR1 && ~isempty(r1Val)
                        diffFields = {'Rdchg_5s', 'Rdchg_10s', 'Rdchg_30s', 'Rdchg_60s'};  % R3-R1 제외
                        diffFieldNames = {'R5_R1', 'R10_R1', 'R30_R1', 'R60_R1'};
                        
                        for dIdx = 1:length(diffFields)
                            rField = diffFields{dIdx};
                            diffFieldName = diffFieldNames{dIdx};
                            
                            if isfield(rVals, rField)
                                diffVal = rVals.(rField) - r1Val;
                                
                                if ~isfield(aggDischargeDiff, socName)
                                    aggDischargeDiff.(socName) = struct();
                                end
                                if ~isfield(aggDischargeDiff.(socName), cycleField)
                                    aggDischargeDiff.(socName).(cycleField) = struct();
                                end
                                if ~isfield(aggDischargeDiff.(socName).(cycleField), diffFieldName)
                                    aggDischargeDiff.(socName).(cycleField).(diffFieldName) = [];
                                end
                                aggDischargeDiff.(socName).(cycleField).(diffFieldName) = ...
                                    [aggDischargeDiff.(socName).(cycleField).(diffFieldName); diffVal];
                            end
                        end
                    end
                end
            end
        end
    end
end

fprintf('Data aggregation complete.\n');

%% 저항값 0 체크
fprintf('\n=== Checking for Zero Resistance Values ===\n');
zeroCountCharge = 0;
zeroCountDischarge = 0;

% 충전 저항값 체크
socNamesCharge = fieldnames(aggCharge);
for s = 1:length(socNamesCharge)
    socName = socNamesCharge{s};
    cycles = fieldnames(aggCharge.(socName));
    for c = 1:length(cycles)
        cycleField = cycles{c};
        rFields = fieldnames(aggCharge.(socName).(cycleField));
        for r = 1:length(rFields)
            rField = rFields{r};
            rVals = aggCharge.(socName).(cycleField).(rField);
            zeroVals = rVals == 0;
            if any(zeroVals)
                zeroCountCharge = zeroCountCharge + sum(zeroVals);
                fprintf('  [Charge] %s | %s | %s: %d zero values found\n', ...
                    socName, cycleField, rField, sum(zeroVals));
            end
        end
    end
end

% 방전 저항값 체크
socNamesDischarge = fieldnames(aggDischarge);
for s = 1:length(socNamesDischarge)
    socName = socNamesDischarge{s};
    cycles = fieldnames(aggDischarge.(socName));
    for c = 1:length(cycles)
        cycleField = cycles{c};
        rFields = fieldnames(aggDischarge.(socName).(cycleField));
        for r = 1:length(rFields)
            rField = rFields{r};
            rVals = aggDischarge.(socName).(cycleField).(rField);
            zeroVals = rVals == 0;
            if any(zeroVals)
                zeroCountDischarge = zeroCountDischarge + sum(zeroVals);
                fprintf('  [Discharge] %s | %s | %s: %d zero values found\n', ...
                    socName, cycleField, rField, sum(zeroVals));
            end
        end
    end
end

if zeroCountCharge == 0 && zeroCountDischarge == 0
    fprintf('  No zero resistance values found.\n');
else
    fprintf('  Total zero values: Charge=%d, Discharge=%d\n', ...
        zeroCountCharge, zeroCountDischarge);
end

%% Fancy 박스플롯 생성
fprintf('\n=== Creating Fancy Boxplots ===\n');

% 충전 저항 시각화
socNamesCharge = fieldnames(aggCharge);
for s = 1:length(socNamesCharge)
    socName = socNamesCharge{s};
    
    % SOC 70만 visible on, 나머지는 off
    if contains(socName, '70') || strcmp(socName, 'SOC70')
        currentVisible = 'on';
    else
        currentVisible = 'off';
    end
    
    % 먼저 모든 데이터를 수집하여 y축 범위 계산
    allValsCharge = [];
    for rIdx = 1:length(rFieldsCharge)
        rField = rFieldsCharge{rIdx};
        [vals, ~] = build_boxplot_data(aggCharge.(socName), rField, cyclesWanted);
        if ~isempty(vals)
            allValsCharge = [allValsCharge; vals];
        end
    end
    
    % y축 범위 계산 (약간의 여백 포함)
    if ~isempty(allValsCharge)
        yMinCharge = min(allValsCharge);
        yMaxCharge = max(allValsCharge);
        yRangeCharge = yMaxCharge - yMinCharge;
        yLimCharge = [yMinCharge - 0.05*yRangeCharge, yMaxCharge + 0.05*yRangeCharge];
    else
        yLimCharge = [0, 1];
    end
    
    fig = figure('Name', sprintf('Fancy Boxplot - Charge - %s', socName), ...
                 'Position', [100, 100, 1800, 1000], 'Visible', currentVisible);
    
    for rIdx = 1:length(rFieldsCharge)
        rField = rFieldsCharge{rIdx};
        subplot(2, 3, rIdx);
        
        % 데이터 준비
        [vals, groups] = build_boxplot_data(aggCharge.(socName), rField, cyclesWanted);
        
        if isempty(vals)
            text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 14, 'Color', [0.5 0.5 0.5]);
            title(rField, 'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
            axis off;
        else
            % Fancy 박스플롯 생성
            create_fancy_boxplot(vals, groups, rField, 'Charge', yLimCharge);
        end
    end
    
    sgtitle(sprintf('%s | %s | Charge Resistance', paramLabel, socName), ...
            'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 저장
    outDir = fullfile(figuresBaseDir, paramLabel);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    saveName = sprintf('FancyBoxplot_%s_%s_Charge', paramLabel, socName);
    saveName = regexprep(saveName, '[^a-zA-Z0-9_]', '_');
    saveas(fig, fullfile(outDir, [saveName '.fig']));
    saveas(fig, fullfile(outDir, [saveName '.png']), 'png');
    fprintf('  Saved: %s\n', saveName);
    
    if strcmp(currentVisible, 'off')
        close(fig);
    end
end

% 충전 저항 차이 시각화 (R5-R1, R10-R1, R30-R1, R60-R1) - R3-R1 제외
diffFieldNames = {'R5_R1', 'R10_R1', 'R30_R1', 'R60_R1'};
socNamesChargeDiff = fieldnames(aggChargeDiff);
for s = 1:length(socNamesChargeDiff)
    socName = socNamesChargeDiff{s};
    
    % SOC 70만 visible on, 나머지는 off
    if contains(socName, '70') || strcmp(socName, 'SOC70')
        currentVisible = 'on';
    else
        currentVisible = 'off';
    end
    
    % 먼저 모든 데이터를 수집하여 y축 범위 계산
    allValsChargeDiff = [];
    for dIdx = 1:length(diffFieldNames)
        diffFieldName = diffFieldNames{dIdx};
        [vals, ~] = build_boxplot_data(aggChargeDiff.(socName), diffFieldName, cyclesWanted);
        if ~isempty(vals)
            allValsChargeDiff = [allValsChargeDiff; vals];
        end
    end
    
    % y축 범위 계산
    if ~isempty(allValsChargeDiff)
        yMinChargeDiff = min(allValsChargeDiff);
        yMaxChargeDiff = max(allValsChargeDiff);
        yRangeChargeDiff = yMaxChargeDiff - yMinChargeDiff;
        yLimChargeDiff = [yMinChargeDiff - 0.05*yRangeChargeDiff, yMaxChargeDiff + 0.05*yRangeChargeDiff];
    else
        yLimChargeDiff = [0, 1];
    end
    
    fig = figure('Name', sprintf('Fancy Boxplot - Charge Diff - %s', socName), ...
                 'Position', [100, 100, 1800, 1000], 'Visible', currentVisible);
    
    for dIdx = 1:length(diffFieldNames)
        diffFieldName = diffFieldNames{dIdx};
        subplot(2, 2, dIdx);  % 2x2 subplot (R3-R1 제외로 4개만)
        
        % 데이터 준비
        [vals, groups] = build_boxplot_data(aggChargeDiff.(socName), diffFieldName, cyclesWanted);
        
        if isempty(vals)
            text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 14, 'Color', [0.5 0.5 0.5]);
            title(diffFieldName, 'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
            axis off;
        else
            % Fancy 박스플롯 생성
            create_fancy_boxplot(vals, groups, diffFieldName, 'Charge', yLimChargeDiff);
        end
    end
    
    sgtitle(sprintf('%s | %s | Charge Resistance Difference (R-R1)', paramLabel, socName), ...
            'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 저장
    outDir = fullfile(figuresBaseDir, paramLabel);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    saveName = sprintf('FancyBoxplot_%s_%s_Charge_Diff', paramLabel, socName);
    saveName = regexprep(saveName, '[^a-zA-Z0-9_]', '_');
    saveas(fig, fullfile(outDir, [saveName '.fig']));
    saveas(fig, fullfile(outDir, [saveName '.png']), 'png');
    fprintf('  Saved: %s\n', saveName);
    
    if strcmp(currentVisible, 'off')
        close(fig);
    end
end

% 방전 저항 시각화
socNamesDischarge = fieldnames(aggDischarge);
for s = 1:length(socNamesDischarge)
    socName = socNamesDischarge{s};
    
    % SOC 70만 visible on, 나머지는 off
    if contains(socName, '70') || strcmp(socName, 'SOC70')
        currentVisible = 'on';
    else
        currentVisible = 'off';
    end
    
    % 먼저 모든 데이터를 수집하여 y축 범위 계산
    allValsDischarge = [];
    for rIdx = 1:length(rFieldsDischarge)
        rField = rFieldsDischarge{rIdx};
        [vals, ~] = build_boxplot_data(aggDischarge.(socName), rField, cyclesWanted);
        if ~isempty(vals)
            allValsDischarge = [allValsDischarge; vals];
        end
    end
    
    % y축 범위 계산 (약간의 여백 포함)
    if ~isempty(allValsDischarge)
        yMinDischarge = min(allValsDischarge);
        yMaxDischarge = max(allValsDischarge);
        yRangeDischarge = yMaxDischarge - yMinDischarge;
        yLimDischarge = [yMinDischarge - 0.05*yRangeDischarge, yMaxDischarge + 0.05*yRangeDischarge];
    else
        yLimDischarge = [0, 1];
    end
    
    fig = figure('Name', sprintf('Fancy Boxplot - Discharge - %s', socName), ...
                 'Position', [100, 100, 1800, 1000], 'Visible', currentVisible);
    
    for rIdx = 1:length(rFieldsDischarge)
        rField = rFieldsDischarge{rIdx};
        subplot(2, 3, rIdx);
        
        % 데이터 준비
        [vals, groups] = build_boxplot_data(aggDischarge.(socName), rField, cyclesWanted);
        
        if isempty(vals)
            text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 14, 'Color', [0.5 0.5 0.5]);
            title(rField, 'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
            axis off;
        else
            % Fancy 박스플롯 생성
            create_fancy_boxplot(vals, groups, rField, 'Discharge', yLimDischarge);
        end
    end
    
    sgtitle(sprintf('%s | %s | Discharge Resistance', paramLabel, socName), ...
            'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 저장
    outDir = fullfile(figuresBaseDir, paramLabel);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    saveName = sprintf('FancyBoxplot_%s_%s_Discharge', paramLabel, socName);
    saveName = regexprep(saveName, '[^a-zA-Z0-9_]', '_');
    saveas(fig, fullfile(outDir, [saveName '.fig']));
    saveas(fig, fullfile(outDir, [saveName '.png']), 'png');
    fprintf('  Saved: %s\n', saveName);
    
    if strcmp(currentVisible, 'off')
        close(fig);
    end
end

% 방전 저항 차이 시각화 (R3-R1, R5-R1, R10-R1, R30-R1, R60-R1)
socNamesDischargeDiff = fieldnames(aggDischargeDiff);
for s = 1:length(socNamesDischargeDiff)
    socName = socNamesDischargeDiff{s};
    
    % SOC 70만 visible on, 나머지는 off
    if contains(socName, '70') || strcmp(socName, 'SOC70')
        currentVisible = 'on';
    else
        currentVisible = 'off';
    end
    
    % 먼저 모든 데이터를 수집하여 y축 범위 계산
    allValsDischargeDiff = [];
    for dIdx = 1:length(diffFieldNames)
        diffFieldName = diffFieldNames{dIdx};
        [vals, ~] = build_boxplot_data(aggDischargeDiff.(socName), diffFieldName, cyclesWanted);
        if ~isempty(vals)
            allValsDischargeDiff = [allValsDischargeDiff; vals];
        end
    end
    
    % y축 범위 계산
    if ~isempty(allValsDischargeDiff)
        yMinDischargeDiff = min(allValsDischargeDiff);
        yMaxDischargeDiff = max(allValsDischargeDiff);
        yRangeDischargeDiff = yMaxDischargeDiff - yMinDischargeDiff;
        yLimDischargeDiff = [yMinDischargeDiff - 0.05*yRangeDischargeDiff, yMaxDischargeDiff + 0.05*yRangeDischargeDiff];
    else
        yLimDischargeDiff = [0, 1];
    end
    
    fig = figure('Name', sprintf('Fancy Boxplot - Discharge Diff - %s', socName), ...
                 'Position', [100, 100, 1800, 1000], 'Visible', currentVisible);
    
    for dIdx = 1:length(diffFieldNames)
        diffFieldName = diffFieldNames{dIdx};
        subplot(2, 2, dIdx);  % 2x2 subplot (R3-R1 제외로 4개만)
        
        % 데이터 준비
        [vals, groups] = build_boxplot_data(aggDischargeDiff.(socName), diffFieldName, cyclesWanted);
        
        if isempty(vals)
            text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 14, 'Color', [0.5 0.5 0.5]);
            title(diffFieldName, 'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
            axis off;
        else
            % Fancy 박스플롯 생성
            create_fancy_boxplot(vals, groups, diffFieldName, 'Discharge', yLimDischargeDiff);
        end
    end
    
    sgtitle(sprintf('%s | %s | Discharge Resistance Difference (R-R1)', paramLabel, socName), ...
            'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 저장
    outDir = fullfile(figuresBaseDir, paramLabel);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    saveName = sprintf('FancyBoxplot_%s_%s_Discharge_Diff', paramLabel, socName);
    saveName = regexprep(saveName, '[^a-zA-Z0-9_]', '_');
    saveas(fig, fullfile(outDir, [saveName '.fig']));
    saveas(fig, fullfile(outDir, [saveName '.png']), 'png');
    fprintf('  Saved: %s\n', saveName);
    
    if strcmp(currentVisible, 'off')
        close(fig);
    end
end

fprintf('\n=== Fancy Boxplot Visualization Complete ===\n');

%% 로컬 함수들

function [vals, groups] = build_boxplot_data(dataStruct, rField, cyclesWanted)
    vals = [];
    groups = categorical(cell(0, 1), cyclesWanted, 'Ordinal', true);
    for c = 1:length(cyclesWanted)
        cycleName = cyclesWanted{c};
        cycleField = make_valid_cycle_field(cycleName);
        if isfield(dataStruct, cycleField) && isfield(dataStruct.(cycleField), rField)
            rVals = dataStruct.(cycleField).(rField);
            rVals = rVals(:);
            rVals = rVals(isfinite(rVals) & rVals ~= 0);  % 0 값 제외
            if ~isempty(rVals)
                vals = [vals; rVals];
                groups = [groups; repmat(categorical({cycleName}, cyclesWanted, 'Ordinal', true), numel(rVals), 1)];
            end
        end
    end
end

function cycleField = make_valid_cycle_field(cycleName)
    if ~isempty(regexp(cycleName, '^\d', 'once'))
        cycleField = ['cyc_' cycleName];
    else
        cycleField = cycleName;
    end
end

function create_fancy_boxplot(vals, groups, rField, eventType, yLimits)
    % daboxplot 사용 여부 확인
    useDaboxplot = exist('daboxplot', 'file') == 2;
    
    if useDaboxplot
        % daboxplot 사용
        % 데이터를 행렬로 변환: 각 열이 하나의 사이클(조건)
        uniqueGroups = unique(groups);
        nCycles = length(uniqueGroups);
        
        % 각 사이클별 데이터를 열로 구성
        dataMatrix = [];
        cycleLabels = cell(nCycles, 1);
        maxLen = 0;
        
        % 먼저 최대 길이 확인
        for c = 1:nCycles
            groupVals = vals(groups == uniqueGroups(c));
            if length(groupVals) > maxLen
                maxLen = length(groupVals);
            end
        end
        
        % 데이터 행렬 생성 (NaN으로 패딩)
        dataMatrix = NaN(maxLen, nCycles);
        for c = 1:nCycles
            groupVals = vals(groups == uniqueGroups(c));
            cycleLabels{c} = char(uniqueGroups(c));
            dataMatrix(1:length(groupVals), c) = groupVals(:);
        end
        
        % 색상 팔레트 설정 (그라데이션)
        % 충전: 빨강, 방전: 파랑
        if strcmp(eventType, 'Charge')
            baseColor = [0.8 0.2 0.2];  % 빨간색 계열
            colors = [linspace(baseColor(1), 0.9, nCycles)', ...
                      linspace(baseColor(2), 0.4, nCycles)', ...
                      linspace(baseColor(3), 0.4, nCycles)'];
        else  % Discharge
            baseColor = [0.2 0.4 0.8];  % 파란색 계열
            colors = [linspace(baseColor(1), 0.4, nCycles)', ...
                      linspace(baseColor(2), 0.7, nCycles)', ...
                      linspace(baseColor(3), 0.9, nCycles)'];
        end
        
        % daboxplot 생성 (여러 조건, 1개 그룹)
        % 이미지 스타일: 테두리만, 위스커 표시, 아웃라이어 표시
        % try-catch로 에러 처리
        daboxplotSuccess = false;
        try
            h = daboxplot(dataMatrix, ...
                'fill', 0, ...                      % 테두리만 (채워지지 않음)
                'colors', colors, ...                % 색상
                'whiskers', 1, ...                   % 위스커 표시
                'scatter', 0, ...                    % 스캐터 없음
                'mean', 0, ...                       % 평균 표시 안 함
                'outliers', 0, ...                   % 아웃라이어 표시 안 함
                'linkline', 1, ...                   % 그룹 내 박스 연결선 (interaction 효과 강조)
                'xtlabels', cycleLabels, ...         % x축 레이블
                'boxspacing', 0.8, ...               % 박스 간격
                'boxwidth', 0.7);                    % 박스 너비
            daboxplotSuccess = true;
        catch ME
            warning('daboxplot failed: %s. Using alternative options.', ME.message);
            % 대안 옵션으로 재시도 (더 단순한 옵션)
            try
                h = daboxplot(dataMatrix, ...
                    'fill', 0, ...
                    'colors', colors, ...
                    'whiskers', 0, ...               % 위스커 없이 시도
                    'scatter', 0, ...
                    'mean', 0, ...
                    'outliers', 0, ...              % 아웃라이어 표시 안 함
                    'linkline', 1, ...              % 그룹 내 박스 연결선
                    'xtlabels', cycleLabels, ...
                    'boxspacing', 0.8, ...
                    'boxwidth', 0.7);
                daboxplotSuccess = true;
            catch ME2
                warning('daboxplot failed again: %s. Falling back to default boxplot.', ME2.message);
                daboxplotSuccess = false;
            end
        end
        
        if daboxplotSuccess
            % 박스 및 선 두께 조정 (박스만 1.5배, 위스커, 연결선은 반으로)
            if isfield(h, 'bx')
                for j = 1:numel(h.bx)
                    if isvalid(h.bx(j))
                        currentWidth = get(h.bx(j), 'LineWidth');
                        set(h.bx(j), 'LineWidth', currentWidth * 1.5);  % 박스 두께 1.5배
                    end
                end
            end
            if isfield(h, 'md')
                for j = 1:numel(h.md)
                    if isvalid(h.md(j))
                        % 중앙선은 유지 (변경 없음)
                    end
                end
            end
            if isfield(h, 'wh')
                for j = 1:numel(h.wh)
                    if isvalid(h.wh(j))
                        currentWidth = get(h.wh(j), 'LineWidth');
                        set(h.wh(j), 'LineWidth', currentWidth * 0.5);  % 위스커 두께 반으로
                    end
                end
            end
            if isfield(h, 'ln')
                for j = 1:numel(h.ln)
                    if isvalid(h.ln(j))
                        currentWidth = get(h.ln(j), 'LineWidth');
                        set(h.ln(j), 'LineWidth', currentWidth * 0.5);  % 연결선 두께 반으로
                    end
                end
            end
            
            % y축 범위 설정
            ylim(yLimits);
            
            % 제목 및 레이블
            title(rField, 'Interpreter', 'none', 'FontSize', 13, 'FontWeight', 'bold');
            xlabel('Cycle', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Resistance (mΩ)', 'FontSize', 12, 'FontWeight', 'bold');
            
            % 그리드 (minor grid 제거하여 덜 촘촘하게)
            grid on;
            set(gca, 'GridAlpha', 0.3);
            
            % 축 스타일
            set(gca, 'FontSize', 11, 'LineWidth', 1);
            
            % 배경색
            set(gca, 'Color', [0.98 0.98 0.98]);
        else
            useDaboxplot = false;  % 기본 boxplot 사용하도록 설정
        end
    end
    
    if ~useDaboxplot
        % 기본 boxplot 사용
        % 기본 boxplot 사용 (daboxplot이 없는 경우)
        % 색상: 충전=빨강, 방전=파랑
        nGroups = length(unique(groups));
        if strcmp(eventType, 'Charge')
            baseColor = [0.8 0.2 0.2];  % 빨간색 계열
            colors = [linspace(baseColor(1), 0.9, nGroups)', ...
                      linspace(baseColor(2), 0.4, nGroups)', ...
                      linspace(baseColor(3), 0.4, nGroups)'];
        else  % Discharge
            baseColor = [0.2 0.4 0.8];  % 파란색 계열
            colors = [linspace(baseColor(1), 0.4, nGroups)', ...
                      linspace(baseColor(2), 0.7, nGroups)', ...
                      linspace(baseColor(3), 0.9, nGroups)'];
        end
        
        % 박스플롯 생성
        bp = boxplot(vals, groups, ...
            'Symbol', 'o', ...
            'OutlierSize', 4, ...
            'Colors', colors, ...
            'Widths', 0.6, ...
            'Positions', 1:nGroups);
        
        % 박스플롯 스타일 개선 (박스만 1.5배, 위스커는 반으로)
        h = findobj(gca, 'Tag', 'Box');
        for j = 1:length(h)
            patch(get(h(j), 'XData'), get(h(j), 'YData'), ...
                  colors(j,:), 'FaceAlpha', 0.7, 'EdgeColor', colors(j,:)*0.6, 'LineWidth', 2.8125);  % 박스 두께 1.5배 (1.875 -> 2.8125)
        end
        
        % 중앙선 스타일 (유지)
        h = findobj(gca, 'Tag', 'Median');
        for j = 1:length(h)
            set(h(j), 'Color', 'k', 'LineWidth', 5);
        end
        
        % 위스커 스타일 (반으로)
        h = findobj(gca, 'Tag', 'Upper Whisker');
        for j = 1:length(h)
            set(h(j), 'LineStyle', '-', 'LineWidth', 1.875, 'Color', colors(j,:)*0.8);  % 위스커 두께 반으로 (3.75 -> 1.875)
        end
        h = findobj(gca, 'Tag', 'Lower Whisker');
        for j = 1:length(h)
            set(h(j), 'LineStyle', '-', 'LineWidth', 1.875, 'Color', colors(j,:)*0.8);  % 위스커 두께 반으로 (3.75 -> 1.875)
        end
        
        % 아웃라이어 스타일
        h = findobj(gca, 'Tag', 'Outliers');
        for j = 1:length(h)
            set(h(j), 'MarkerEdgeColor', colors(j,:)*0.6, 'MarkerFaceColor', colors(j,:)*0.3, ...
                      'MarkerSize', 5);
        end
        
        % y축 범위 설정
        ylim(yLimits);
        
        % 제목 및 레이블
        title(rField, 'Interpreter', 'none', 'FontSize', 13, 'FontWeight', 'bold');
        xlabel('Cycle', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Resistance (mΩ)', 'FontSize', 12, 'FontWeight', 'bold');
        
        % 그리드
        grid on;
        grid minor;
        set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
        
        % 축 스타일
        set(gca, 'FontSize', 11, 'LineWidth', 1);
        set(gca, 'Box', 'on');
        
        % 통계 정보 표시 (평균값)
        uniqueGroups = unique(groups);
        for j = 1:length(uniqueGroups)
            groupVals = vals(groups == uniqueGroups(j));
            if ~isempty(groupVals)
                meanVal = mean(groupVals);
                text(j, meanVal, sprintf('μ=%.2f', meanVal), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                     'FontSize', 9, 'Color', [0.3 0.3 0.3], 'FontWeight', 'bold');
            end
        end
        
        % 배경색
        set(gca, 'Color', [0.98 0.98 0.98]);
    end
end
