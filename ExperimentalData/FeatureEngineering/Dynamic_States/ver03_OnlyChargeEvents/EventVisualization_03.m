%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 03_EventVisualization.m
% 목적: 조합/채널/SOC별 시간저항(Rchg) 박스플롯 시각화
%
% 입력:
%   - Results/<paramLabel>/Lab_DC_Events_Features_*cyc.mat (02번 스크립트 출력)
%
% 출력:
%   - Results/figures/EventVisualization_Boxplots/<paramLabel>/*.fig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Event Visualization (Boxplots) ===\n');

%% 설정
baseResultsDir = fullfile(pwd, 'Results');
figuresBaseDir = fullfile(baseResultsDir, 'figures', 'EventVisualization_Boxplots');
if ~exist(figuresBaseDir, 'dir')
    mkdir(figuresBaseDir);
end

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

% 결과 폴더 이름에서 사용할 이벤트 타입 (폴더명과 맞춰야 함)
event_type_selection = 'Both';

cyclesWanted = {'0cyc','200cyc','400cyc','600cyc','800cyc'};
timePoints = [1, 3, 5, 10, 30, 60];
rFields = arrayfun(@(t) sprintf('Rchg_%ds', t), timePoints, 'UniformOutput', false);

figVisible = 'off';

fprintf('Base results directory: %s\n', baseResultsDir);
fprintf('Event type selection: %s\n', event_type_selection);

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
    matFiles = dir(fullfile(inputDir, 'Lab_DC_Events_Features_*cyc.mat'));
    if isempty(matFiles)
        fprintf('  WARNING: No feature files found in %s\n', inputDir);
        continue;
    end
    
    fprintf('\n=== ParamSet: %s ===\n', paramLabel);
    fprintf('Input directory: %s\n', inputDir);
    
    % 집계 구조체
    agg = struct();           % charge only: agg.(soc).(channel).(cycle).(rField)
    aggAllCharge = struct();  % aggAllCharge.(soc).(cycle).(rField)
    aggAllDischarge = struct(); % aggAllDischarge.(soc).(cycle).(rField)
    
    for i = 1:length(matFiles)
        fileName = matFiles(i).name;
        filePath = fullfile(inputDir, fileName);
        token = regexp(fileName, 'Lab_DC_Events_Features_(\d+cyc)\.mat', 'tokens', 'once');
        if isempty(token)
            continue;
        end
        cycleName = token{1};
        cycleField = make_valid_cycle_field(cycleName);
        
        dataStruct = load(filePath);
        vars = fieldnames(dataStruct);
        dataVarName = '';
        for v = 1:length(vars)
            if contains(vars{v}, 'Lab_DC_DCIR_')
                dataVarName = vars{v};
                break;
            end
        end
        if isempty(dataVarName)
            dataVarName = vars{1};
        end
        data = dataStruct.(dataVarName);
        
        chFields = fieldnames(data);
        chargeCh = chFields(contains(chFields, '_Charge'));
        dischargeCh = chFields(contains(chFields, '_Discharge'));
        
        for c = 1:length(chargeCh)
            chName = chargeCh{c};
            socs = fieldnames(data.(chName));
            for s = 1:length(socs)
                socName = socs{s};
                profs = fieldnames(data.(chName).(socName));
                for pIdx = 1:length(profs)
                    profName = profs{pIdx};
                    evts = fieldnames(data.(chName).(socName).(profName));
                    for e = 1:length(evts)
                        evt = data.(chName).(socName).(profName).(evts{e});
                        for rIdx = 1:length(rFields)
                            rField = rFields{rIdx};
                            if ~isfield(evt, rField)
                                continue;
                            end
                            rVal = evt.(rField);
                            if ~isfinite(rVal)
                                continue;
                            end
                            
                            if ~isfield(agg, socName)
                                agg.(socName) = struct();
                            end
                            if ~isfield(agg.(socName), chName)
                                agg.(socName).(chName) = struct();
                            end
                            if ~isfield(agg.(socName).(chName), cycleField)
                                agg.(socName).(chName).(cycleField) = struct();
                            end
                            if ~isfield(agg.(socName).(chName).(cycleField), rField)
                                agg.(socName).(chName).(cycleField).(rField) = [];
                            end
                            agg.(socName).(chName).(cycleField).(rField) = ...
                                [agg.(socName).(chName).(cycleField).(rField); rVal];
                            
                            if ~isfield(aggAllCharge, socName)
                                aggAllCharge.(socName) = struct();
                            end
                            if ~isfield(aggAllCharge.(socName), cycleField)
                                aggAllCharge.(socName).(cycleField) = struct();
                            end
                            if ~isfield(aggAllCharge.(socName).(cycleField), rField)
                                aggAllCharge.(socName).(cycleField).(rField) = [];
                            end
                            aggAllCharge.(socName).(cycleField).(rField) = ...
                                [aggAllCharge.(socName).(cycleField).(rField); rVal];
                        end
                    end
                end
            end
        end

        rFieldsDis = strrep(rFields, 'Rchg_', 'Rdchg_');
        for c = 1:length(dischargeCh)
            chName = dischargeCh{c};
            socs = fieldnames(data.(chName));
            for s = 1:length(socs)
                socName = socs{s};
                profs = fieldnames(data.(chName).(socName));
                for pIdx = 1:length(profs)
                    profName = profs{pIdx};
                    evts = fieldnames(data.(chName).(socName).(profName));
                    for e = 1:length(evts)
                        evt = data.(chName).(socName).(profName).(evts{e});
                        for rIdx = 1:length(rFieldsDis)
                            rField = rFieldsDis{rIdx};
                            if ~isfield(evt, rField)
                                continue;
                            end
                            rVal = evt.(rField);
                            if ~isfinite(rVal)
                                continue;
                            end
                            
                            if ~isfield(aggAllDischarge, socName)
                                aggAllDischarge.(socName) = struct();
                            end
                            if ~isfield(aggAllDischarge.(socName), cycleField)
                                aggAllDischarge.(socName).(cycleField) = struct();
                            end
                            if ~isfield(aggAllDischarge.(socName).(cycleField), rField)
                                aggAllDischarge.(socName).(cycleField).(rField) = [];
                            end
                            aggAllDischarge.(socName).(cycleField).(rField) = ...
                                [aggAllDischarge.(socName).(cycleField).(rField); rVal];
                        end
                    end
                end
            end
        end
    end
    
    if isempty(fieldnames(agg)) && isempty(fieldnames(aggAllCharge)) && isempty(fieldnames(aggAllDischarge))
        fprintf('  WARNING: No Rchg data found for %s\n', paramLabel);
        continue;
    end
    
    % 출력 폴더
    outDir = fullfile(figuresBaseDir, paramLabel);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    
    % SOC별 채널 시각화
    socNames = fieldnames(agg);
    for s = 1:length(socNames)
        socName = socNames{s};
        chNames = fieldnames(agg.(socName));
        for c = 1:length(chNames)
            chName = chNames{c};
            
            fig = figure('Name', sprintf('Boxplot %s %s %s', paramLabel, socName, chName), ...
                         'Position', [100, 100, 1600, 900], 'Visible', figVisible);
            
            for rIdx = 1:length(rFields)
                rField = rFields{rIdx};
                subplot(2, 3, rIdx);
                [vals, groups] = build_boxplot_data(agg.(socName).(chName), rField, cyclesWanted);
                if isempty(vals)
                    text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
                         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                    title(rField, 'Interpreter', 'none');
                    axis off;
                else
                    boxplot(vals, groups, 'Symbol', '');
                    title(rField, 'Interpreter', 'none');
                    xlabel('Cycle');
                    ylabel('Resistance (mOhm)');
                    grid on;
                end
            end
            
            sgtitle(sprintf('%s | %s | %s', paramLabel, socName, chName), ...
                    'Interpreter', 'none', 'FontWeight', 'bold');
            
            chDir = fullfile(outDir, 'channels', socName);
            if ~exist(chDir, 'dir')
                mkdir(chDir);
            end
            saveName = sprintf('Boxplot_%s_%s_%s.fig', paramLabel, socName, chName);
            saveName = regexprep(saveName, '[^a-zA-Z0-9_\.]', '_');
            saveas(fig, fullfile(chDir, saveName));
            close(fig);
        end
    end
    
    % SOC별 전체 채널 통합 시각화
    socNamesAll = fieldnames(aggAllCharge);
    for s = 1:length(socNamesAll)
        socName = socNamesAll{s};
        
        fig = figure('Name', sprintf('Boxplot %s %s AllChannels Charge', paramLabel, socName), ...
                     'Position', [100, 100, 1600, 900], 'Visible', figVisible);
        
        for rIdx = 1:length(rFields)
            rField = rFields{rIdx};
            subplot(2, 3, rIdx);
            [vals, groups] = build_boxplot_data(aggAllCharge.(socName), rField, cyclesWanted);
            pVal = compute_pvalue(vals, groups);
            if isempty(vals)
                text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                title(rField, 'Interpreter', 'none');
                axis off;
            else
                boxplot(vals, groups, 'Symbol', '');
                title(rField, 'Interpreter', 'none');
                xlabel('Cycle');
                ylabel('Resistance (mOhm)');
                grid on;
            end
            fprintf('[PVAL][%s][%s][AllChannels][Charge] %s: %s\n', ...
                paramLabel, socName, rField, pvalue_str(pVal));
        end
        
        sgtitle(sprintf('%s | %s | AllChannels | Charge', paramLabel, socName), ...
                'Interpreter', 'none', 'FontWeight', 'bold');
        
        allDir = fullfile(outDir, 'all_channels', socName);
        if ~exist(allDir, 'dir')
            mkdir(allDir);
        end
        saveName = sprintf('Boxplot_%s_%s_AllChannels_Charge.fig', paramLabel, socName);
        saveName = regexprep(saveName, '[^a-zA-Z0-9_\.]', '_');
        saveas(fig, fullfile(allDir, saveName));
        close(fig);
    end
    
    socNamesAllDis = fieldnames(aggAllDischarge);
    for s = 1:length(socNamesAllDis)
        socName = socNamesAllDis{s};
        
        fig = figure('Name', sprintf('Boxplot %s %s AllChannels Discharge', paramLabel, socName), ...
                     'Position', [100, 100, 1600, 900], 'Visible', figVisible);
        
        for rIdx = 1:length(rFields)
            rField = strrep(rFields{rIdx}, 'Rchg_', 'Rdchg_');
            subplot(2, 3, rIdx);
            [vals, groups] = build_boxplot_data(aggAllDischarge.(socName), rField, cyclesWanted);
            pVal = compute_pvalue(vals, groups);
            if isempty(vals)
                text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                title(rField, 'Interpreter', 'none');
                axis off;
            else
                boxplot(vals, groups, 'Symbol', '');
                title(rField, 'Interpreter', 'none');
                xlabel('Cycle');
                ylabel('Resistance (mOhm)');
                grid on;
            end
            fprintf('[PVAL][%s][%s][AllChannels][Discharge] %s: %s\n', ...
                paramLabel, socName, rField, pvalue_str(pVal));
        end
        
        sgtitle(sprintf('%s | %s | AllChannels | Discharge', paramLabel, socName), ...
                'Interpreter', 'none', 'FontWeight', 'bold');
        
        allDir = fullfile(outDir, 'all_channels', socName);
        if ~exist(allDir, 'dir')
            mkdir(allDir);
        end
        saveName = sprintf('Boxplot_%s_%s_AllChannels_Discharge.fig', paramLabel, socName);
        saveName = regexprep(saveName, '[^a-zA-Z0-9_\.]', '_');
        saveas(fig, fullfile(allDir, saveName));
        close(fig);
    end
end

fprintf('\n=== Boxplot Visualization Complete ===\n');

%% 로컬 함수
function [vals, groups] = build_boxplot_data(dataStruct, rField, cyclesWanted)
    vals = [];
    groups = categorical(cell(0, 1), cyclesWanted, 'Ordinal', true);
    for c = 1:length(cyclesWanted)
        cycleName = cyclesWanted{c};
        cycleField = make_valid_cycle_field(cycleName);
        if isfield(dataStruct, cycleField) && isfield(dataStruct.(cycleField), rField)
            rVals = dataStruct.(cycleField).(rField);
            rVals = rVals(:);
            if ~isempty(rVals)
                vals = [vals; rVals];
                groups = [groups; repmat(categorical({cycleName}, cyclesWanted, 'Ordinal', true), numel(rVals), 1)];
            end
        end
    end
    
    if isempty(vals)
        return;
    end
    for c = 1:length(cyclesWanted)
        vals = [vals; NaN];
        groups = [groups; categorical({cyclesWanted{c}}, cyclesWanted, 'Ordinal', true)];
    end
end

function cycleField = make_valid_cycle_field(cycleName)
    if ~isempty(regexp(cycleName, '^\d', 'once'))
        cycleField = ['cyc_' cycleName];
    else
        cycleField = cycleName;
    end
end

function pVal = compute_pvalue(vals, groups)
    pVal = NaN;
    if isempty(vals) || isempty(groups)
        return;
    end
    mask = ~isnan(vals);
    if ~any(mask)
        return;
    end
    vals = vals(mask);
    groups = groups(mask);
    if numel(unique(groups)) < 2
        return;
    end
    pVal = kruskalwallis(vals, groups, 'off');
end

function outStr = pvalue_str(pVal)
    if isnan(pVal)
        outStr = 'NaN';
    else
        outStr = sprintf('%.4g', pVal);
    end
end
