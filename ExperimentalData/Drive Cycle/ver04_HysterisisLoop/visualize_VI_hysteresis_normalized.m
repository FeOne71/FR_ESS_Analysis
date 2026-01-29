%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-I 궤적(Hysteresis Loop) 정규화 시각화 스크립트
% 목적: SOC 90, 70, 50 데이터에 대해 전압 정규화(ΔV)를 적용하여 비교
%       - ΔV = V_terminal - V_OCV(SOC)
%       - SOC별 Dynamics 변화 비교
%       - 채널별 SOH 비교
% 
% 분석 항목:
%   1. 기울기 (Slope): Ohmic Resistance
%   2. 루프 면적 (Area): Diffusion 성능
%   3. 비대칭성 (Asymmetry): 충전/방전 Kinetics 차이
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 설정
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
ocvDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';
ocvDataFile = fullfile(ocvDataDir, 'OCV_integrated.mat');
targetSOCs = {'SOC90', 'SOC70', 'SOC50'};  % 분석할 SOC 레벨
outputDir = fullfile(pwd, 'VI_Hysteresis_Plots_Normalized');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% OCV 데이터 로드
fprintf('=== Loading OCV Data ===\n');
if ~exist(ocvDataFile, 'file')
    error('OCV data file not found: %s\nPlease run RPT_Postprocessing.m first.', ocvDataFile);
end

ocvData = load(ocvDataFile, 'OCV_data');
OCV_data = ocvData.OCV_data;
fprintf('OCV data loaded from: %s\n', ocvDataFile);

%% OCV 데이터 시각화
fprintf('\n=== Visualizing OCV Data ===\n');
ocvFields = fieldnames(OCV_data);
ocvCycleFields = {};
for f = 1:length(ocvFields)
    fieldName = ocvFields{f};
    if contains(fieldName, 'OCV_integrated')
        ocvCycleFields{end+1} = fieldName;
    end
end

if ~isempty(ocvCycleFields)
    % 사이클 번호 추출 및 정렬
    cycleNums = zeros(length(ocvCycleFields), 1);
    for f = 1:length(ocvCycleFields)
        fieldName = ocvCycleFields{f};
        numMatch = regexp(fieldName, '(\d+)cyc', 'tokens');
        if ~isempty(numMatch)
            cycleNums(f) = str2double(numMatch{1}{1});
        else
            cycleNums(f) = 999;
        end
    end
    [~, sortIdx] = sort(cycleNums);
    ocvCycleFields = ocvCycleFields(sortIdx);
    
    % OCV 곡선 플롯
    fig_ocv = figure('Position', [100, 100, 1400, 800]);
    hold on;
    
    colors_ocv = lines(length(ocvCycleFields));
    legendEntries_ocv = {};
    
    for f = 1:length(ocvCycleFields)
        fieldName = ocvCycleFields{f};
        ocvStruct = OCV_data.(fieldName);
        
        % 사이클 번호 추출
        numMatch = regexp(fieldName, '(\d+)cyc', 'tokens');
        if ~isempty(numMatch)
            cycleNum = numMatch{1}{1};
            cycleLabel = sprintf('%scyc', cycleNum);
        else
            cycleLabel = fieldName;
        end
        
        % SOC_grid와 V_avg_SOC 확인
        if isfield(ocvStruct, 'SOC_grid') && isfield(ocvStruct, 'V_avg_SOC')
            soc_grid = ocvStruct.SOC_grid;
            v_avg_soc = ocvStruct.V_avg_SOC;
            
            % 유효한 데이터만 플롯
            validIdx = ~isnan(soc_grid) & ~isnan(v_avg_soc);
            if sum(validIdx) > 0
                plot(soc_grid(validIdx), v_avg_soc(validIdx), 'o-', ...
                     'Color', colors_ocv(f, :), 'LineWidth', 2, ...
                     'MarkerSize', 6, 'DisplayName', cycleLabel);
                legendEntries_ocv{end+1} = cycleLabel;
            end
        elseif isfield(ocvStruct, 'OCV_SOC_func')
            % 함수가 있는 경우 SOC 범위에서 샘플링
            soc_range = 0:5:100;
            v_ocv = zeros(size(soc_range));
            for s = 1:length(soc_range)
                try
                    v_ocv(s) = ocvStruct.OCV_SOC_func(soc_range(s));
                catch
                    v_ocv(s) = NaN;
                end
            end
            validIdx = ~isnan(v_ocv);
            if sum(validIdx) > 0
                plot(soc_range(validIdx), v_ocv(validIdx), 'o-', ...
                     'Color', colors_ocv(f, :), 'LineWidth', 2, ...
                     'MarkerSize', 6, 'DisplayName', cycleLabel);
                legendEntries_ocv{end+1} = cycleLabel;
            end
        end
    end
    
    xlabel('SOC (%)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('OCV (V)', 'FontSize', 12, 'FontWeight', 'bold');
    title('OCV Curves by Cycle', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    if ~isempty(legendEntries_ocv)
        legend(legendEntries_ocv, 'Location', 'best', 'FontSize', 10);
    end
    hold off;
    
    fprintf('OCV visualization completed (figure displayed, not saved)\n');
else
    fprintf('No OCV cycle data found to visualize\n');
end


%% 데이터 파일 찾기
matFiles = dir(fullfile(dataDir, 'parsedDriveCycle_*cyc_filtered.mat'));
if isempty(matFiles)
    error('No data files found in %s', dataDir);
end

fprintf('=== V-I Hysteresis Loop Normalized Visualization ===\n');
fprintf('Data directory: %s\n', dataDir);
fprintf('Target SOCs: %s\n', strjoin(targetSOCs, ', '));
fprintf('Found %d cycle files\n', length(matFiles));

%% 모든 사이클 파일에서 데이터 수집
allData = struct();
cycleNameMap = containers.Map();
cycleNameReverseMap = containers.Map();

for i = 1:length(matFiles)
    fileName = matFiles(i).name;
    token = regexp(fileName, 'parsedDriveCycle_(\d+cyc)_filtered', 'tokens');
    if isempty(token), continue; end
    cycleType = token{1}{1};
    
    % 유효한 필드명으로 변환
    if ~isempty(regexp(cycleType, '^\d', 'once'))
        cycleType_valid = ['cyc_' cycleType];
    else
        cycleType_valid = cycleType;
    end
    cycleType_valid = regexprep(cycleType_valid, '[^a-zA-Z0-9_]', '_');
    
    cycleNameMap(cycleType) = cycleType_valid;
    cycleNameReverseMap(cycleType_valid) = cycleType;
    
    fprintf('\n--- Loading %s (field: %s) ---\n', cycleType, cycleType_valid);
    dataStruct = load(fullfile(dataDir, fileName));
    
    varName = sprintf('parsedDriveCycle_%s', cycleType);
    if ~isfield(dataStruct, varName)
        fprintf('  WARNING: Variable %s not found. Skipping.\n', varName);
        continue;
    end
    data_var = dataStruct.(varName);
    
    channels = fieldnames(data_var);
    fprintf('  Channels found: %s\n', strjoin(channels, ', '));
    
    % 각 채널별로 모든 SOC 데이터 수집
    for chIdx = 1:length(channels)
        chName = channels{chIdx};
        
        if ~isfield(data_var, chName) || ~isstruct(data_var.(chName))
            continue;
        end
        
        % 모든 SOC 필드 확인
        socFields = fieldnames(data_var.(chName));
        
        for sIdx = 1:length(socFields)
            socField = socFields{sIdx};
            socField_upper = upper(socField);
            
            % 타겟 SOC 중 하나인지 확인
            targetSOC_found = false;
            targetSOC_name = '';
            for tIdx = 1:length(targetSOCs)
                targetSOC_upper = upper(targetSOCs{tIdx});
                if contains(socField_upper, targetSOC_upper) || ...
                   (contains(socField_upper, 'SOC') && ...
                    (contains(socField_upper, '90') || contains(socField_upper, '70') || contains(socField_upper, '50')))
                    targetSOC_found = true;
                    % SOC90, SOC70, SOC50 형식으로 정규화
                    if contains(socField_upper, '90')
                        targetSOC_name = 'SOC90';
                    elseif contains(socField_upper, '70')
                        targetSOC_name = 'SOC70';
                    elseif contains(socField_upper, '50')
                        targetSOC_name = 'SOC50';
                    else
                        targetSOC_name = socField;
                    end
                    break;
                end
            end
            
            if ~targetSOC_found
                continue;
            end
            
            socData = data_var.(chName).(socField);
            if ~isstruct(socData)
                continue;
            end
            
            % Profile 필드 찾기
            if isfield(socData, 'Profile') && isstruct(socData.Profile)
                profs = fieldnames(socData.Profile);
                useProfileStruct = true;
            else
                allFields = fieldnames(socData);
                profs = {};
                for fIdx = 1:length(allFields)
                    fieldName = allFields{fIdx};
                    if startsWith(fieldName, 'DC') || startsWith(fieldName, 'Profile')
                        profs{end+1} = fieldName;
                    end
                end
                useProfileStruct = false;
            end
            
            % 각 Profile별 데이터 저장
            for pIdx = 1:length(profs)
                profName = profs{pIdx};
                
                if useProfileStruct
                    if ~isfield(socData.Profile, profName)
                        continue;
                    end
                    profData = socData.Profile.(profName);
                else
                    if ~isfield(socData, profName)
                        continue;
                    end
                    profData = socData.(profName);
                end
                
                if ~isfield(profData, 'I') || ~isfield(profData, 'V') || ~isfield(profData, 't')
                    continue;
                end
                
                I = profData.I;
                V = profData.V;
                t = profData.t;
                
                % 시간 변환
                if isa(t, 'duration')
                    t_seconds = seconds(t);
                else
                    t_seconds = t;
                end
                
                % NaN 제거
                validIdx = ~isnan(I) & ~isnan(V) & ~isnan(t_seconds);
                I = I(validIdx);
                V = V(validIdx);
                t_seconds = t_seconds(validIdx);
                
                if isempty(I) || isempty(V)
                    continue;
                end
                
                % OCV 계산: RPT Postprocessing에서 생성한 통합 OCV 데이터 사용
                % 사이클명에서 숫자 추출: '0cyc' -> '0', '200cyc' -> '200'
                cycleNumToken = regexp(cycleType, '(\d+)cyc', 'tokens');
                if ~isempty(cycleNumToken)
                    cycleNum = cycleNumToken{1}{1};
                else
                    cycleNum = '';
                end
                ocvFieldName = sprintf('OCV_{integrated}_{%s}', cycleNum);
                
                if isfield(OCV_data, ocvFieldName)
                    ocvStruct = OCV_data.(ocvFieldName);
                    if isfield(ocvStruct, 'OCV_SOC_func')
                        % SOC 값 추출 (데이터에서 SOC 정보가 있는지 확인)
                        % SOC는 targetSOC_name에서 추출 (SOC90 -> 90, SOC70 -> 70, SOC50 -> 50)
                        socValue = str2double(regexp(targetSOC_name, '\d+', 'match', 'once'));
                        if isnan(socValue)
                            % SOC 값을 찾을 수 없으면 Idle 상태의 평균 전압 사용
                            Cnom = 64;
                            idle_current_threshold = Cnom * 0.01;
                            idle_idx = abs(I) < idle_current_threshold;
                            if sum(idle_idx) > 10
                                OCV = mean(V(idle_idx));
                            else
                                OCV = (max(V) + min(V)) / 2;
                            end
                        else
                            % OCV 함수를 사용하여 해당 SOC의 OCV 계산
                            OCV = ocvStruct.OCV_SOC_func(socValue);
                        end
                    else
                        % OCV 함수가 없으면 SOC_grid와 V_avg_SOC를 사용하여 보간
                        if isfield(ocvStruct, 'SOC_grid') && isfield(ocvStruct, 'V_avg_SOC')
                            socValue = str2double(regexp(targetSOC_name, '\d+', 'match', 'once'));
                            if ~isnan(socValue)
                                OCV = interp1(ocvStruct.SOC_grid, ocvStruct.V_avg_SOC, socValue, 'linear');
                            else
                                OCV = mean(ocvStruct.V_avg_SOC);
                            end
                        else
                            % 대체 방법: Idle 상태의 평균 전압
                            Cnom = 64;
                            idle_current_threshold = Cnom * 0.01;
                            idle_idx = abs(I) < idle_current_threshold;
                            if sum(idle_idx) > 10
                                OCV = mean(V(idle_idx));
                            else
                                OCV = (max(V) + min(V)) / 2;
                            end
                        end
                    end
                else
                    % 해당 사이클의 OCV 데이터가 없으면 Idle 상태의 평균 전압 사용
                    fprintf('  WARNING: OCV data not found for cycle %s. Using idle voltage.\n', cycleType);
                    Cnom = 64;
                    idle_current_threshold = Cnom * 0.01;
                    idle_idx = abs(I) < idle_current_threshold;
                    if sum(idle_idx) > 10
                        OCV = mean(V(idle_idx));
                    else
                        OCV = (max(V) + min(V)) / 2;
                    end
                end
                
                % ΔV 계산
                deltaV = V - OCV;
                
                % 데이터 저장
                if ~isfield(allData, chName)
                    allData.(chName) = struct();
                end
                
                if ~isfield(allData.(chName), cycleType_valid)
                    allData.(chName).(cycleType_valid) = struct();
                end
                
                if ~isfield(allData.(chName).(cycleType_valid), targetSOC_name)
                    allData.(chName).(cycleType_valid).(targetSOC_name) = struct();
                end
                
                allData.(chName).(cycleType_valid).(targetSOC_name).(profName) = struct();
                allData.(chName).(cycleType_valid).(targetSOC_name).(profName).I = I;
                allData.(chName).(cycleType_valid).(targetSOC_name).(profName).V = V;
                allData.(chName).(cycleType_valid).(targetSOC_name).(profName).deltaV = deltaV;
                allData.(chName).(cycleType_valid).(targetSOC_name).(profName).t = t_seconds;
                allData.(chName).(cycleType_valid).(targetSOC_name).(profName).OCV = OCV;
                
                fprintf('    Saved: %s - %s - %s - %s (OCV=%.3fV)\n', ...
                        chName, cycleType, targetSOC_name, profName, OCV);
            end
        end
    end
end

%% 채널 번호별로 그룹화 및 시각화
channels_all = fieldnames(allData);
fprintf('\n=== Visualizing Normalized V-I Hysteresis Loops ===\n');
fprintf('Total channel entries: %d\n', length(channels_all));

% 채널 번호별로 그룹화 (ch9, ch10, ... ch16)
channelGroups = containers.Map();
for c = 1:length(channels_all)
    chName = channels_all{c};
    % ch9_Drive_0cyc -> ch9
    num_match = regexp(chName, 'ch(\d+)', 'tokens');
    if ~isempty(num_match)
        chNum = num_match{1}{1};
        chGroupKey = sprintf('ch%s', chNum);
        
        if ~isKey(channelGroups, chGroupKey)
            channelGroups(chGroupKey) = {};
        end
        channelList = channelGroups(chGroupKey);
        channelList{end+1} = chName;
        channelGroups(chGroupKey) = channelList;
    end
end

% 채널 그룹 목록 정렬
channelGroupKeys = keys(channelGroups);
channelGroupNums = zeros(length(channelGroupKeys), 1);
for k = 1:length(channelGroupKeys)
    key = channelGroupKeys{k};
    num_match = regexp(key, 'ch(\d+)', 'tokens');
    if ~isempty(num_match)
        channelGroupNums(k) = str2double(num_match{1}{1});
    else
        channelGroupNums(k) = 999;
    end
end
[~, sort_idx] = sort(channelGroupNums);
channelGroupKeys = channelGroupKeys(sort_idx);

fprintf('Found %d channel groups: %s\n', length(channelGroupKeys), strjoin(channelGroupKeys, ', '));

%% 8개 채널 비교 플롯 (특정 사이클, 특정 SOC, 대표 DC)
% 같은 사이클의 채널들을 찾아서 비교
targetSOC_multi = 'SOC70';
representativeDC_multi = 'DC2';  % 대표 DC

% 첫 번째 사이클 찾기 (채널 이름에서 사이클 추출)
firstCycle_display = '';
firstCycle_valid = '';
if ~isempty(channels)
    % 첫 번째 채널에서 사이클 정보 추출
    firstChName = channels{1};
    cycle_match = regexp(firstChName, '_(\d+cyc)', 'tokens');
    if ~isempty(cycle_match)
        firstCycle_display = cycle_match{1}{1};
        % 유효한 필드명으로 변환
        if ~isempty(regexp(firstCycle_display, '^\d', 'once'))
            firstCycle_valid = ['cyc_' firstCycle_display];
        else
            firstCycle_valid = firstCycle_display;
        end
        firstCycle_valid = regexprep(firstCycle_valid, '[^a-zA-Z0-9_]', '_');
    else
        % 채널 이름에서 사이클을 찾을 수 없으면 첫 번째 채널의 사이클 사용
        if isfield(allData, firstChName)
            cycles = fieldnames(allData.(firstChName));
            if ~isempty(cycles)
                firstCycle_valid = cycles{1};
                if isKey(cycleNameReverseMap, firstCycle_valid)
                    firstCycle_display = cycleNameReverseMap(firstCycle_valid);
                else
                    firstCycle_display = firstCycle_valid;
                end
            end
        end
    end
end

if ~isempty(firstCycle_valid) && ~isempty(firstCycle_display)
    % 같은 사이클의 채널들 찾기 (ch9~ch16, 같은 사이클)
    sameCycleChannels = {};
    for chIdx = 1:length(channels)
        chName = channels{chIdx};
        % 채널 이름에서 사이클 확인
        chCycle_match = regexp(chName, '_(\d+cyc)', 'tokens');
        chCycle_display = '';
        if ~isempty(chCycle_match)
            chCycle_display = chCycle_match{1}{1};
        end
        
        % 같은 사이클이고 데이터가 있는 채널만 추가
        if (isempty(chCycle_display) || strcmp(chCycle_display, firstCycle_display)) && ...
           isfield(allData, chName) && ...
           isfield(allData.(chName), firstCycle_valid) && ...
           isfield(allData.(chName).(firstCycle_valid), targetSOC_multi) && ...
           isfield(allData.(chName).(firstCycle_valid).(targetSOC_multi), representativeDC_multi)
            sameCycleChannels{end+1} = chName;
        end
        
        % 8개 채널을 찾으면 중단
        if length(sameCycleChannels) >= 8
            break;
        end
    end
    
    if length(sameCycleChannels) >= 8
        % 8개 채널 비교 플롯
        fig_multi = figure('Position', [100, 100, 2000, 1200]);
        sgtitle(sprintf('8-Channel Comparison: %s - %s (%s)', firstCycle_display, targetSOC_multi, representativeDC_multi), ...
                'FontSize', 16, 'FontWeight', 'bold');
        
        % 8개 채널을 2x4 서브플롯으로 배치
        for chIdx = 1:min(8, length(sameCycleChannels))
            chName = sameCycleChannels{chIdx};
            
            subplot(2, 4, chIdx);
            hold on;
            
            % 채널 번호 추출
            num_match = regexp(chName, 'ch(\d+)', 'tokens');
            if ~isempty(num_match)
                chNum = num_match{1}{1};
            else
                chNum = chName;
            end
            
            data = allData.(chName).(firstCycle_valid).(targetSOC_multi).(representativeDC_multi);
            I = data.I;
            deltaV = data.deltaV;
            t_seconds = data.t;
            
            % 기울기 계산 및 가이드라인
            valid_idx = abs(I) > 0.1;
            slope = NaN;
            area = NaN;
            if sum(valid_idx) > 10
                p = polyfit(I(valid_idx), deltaV(valid_idx), 1);
                slope = p(1);
                intercept = p(2);
                I_range = [min(I), max(I)];
                deltaV_fit = I_range * slope + intercept;
                plot(I_range, deltaV_fit, 'r--', 'LineWidth', 1.5);
            end
            
            % 루프 면적 계산
            if length(I) > 2
                area = abs(polyarea(I, deltaV));
            end
            
            % 데이터가 너무 많으면 샘플링 (성능 개선)
            max_points = 5000;
            if length(I) > max_points
                sample_idx = round(linspace(1, length(I), max_points));
                I_plot = I(sample_idx);
                deltaV_plot = deltaV(sample_idx);
                t_plot = t_seconds(sample_idx);
            else
                I_plot = I;
                deltaV_plot = deltaV;
                t_plot = t_seconds;
            end
            
            % 시간 색상 매핑
            scatter(I_plot, deltaV_plot, 20, t_plot, 'filled');
            plot(I_plot, deltaV_plot, 'k-', 'LineWidth', 0.8);
            
            xlabel('Current (A)', 'FontSize', 10);
            ylabel('\Delta V (V)', 'FontSize', 10);
            title(sprintf('Ch%s', chNum), 'FontSize', 11, 'FontWeight', 'bold');
            grid on;
            
            % 분석 통계 표시
            if ~isnan(slope) && ~isnan(area)
                text(0.02, 0.98, sprintf('R=%.3f\nA=%.4f', slope, area), ...
                     'Units', 'normalized', 'VerticalAlignment', 'top', ...
                     'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
            end
            hold off;
        end
        
        % 색상바 추가
        colorbar;
        colormap('jet');
        c = colorbar;
        c.Label.String = 'Time (s)';
        
        savePath_multi = fullfile(outputDir, sprintf('8Channel_Comparison_%s_%s_%s.fig', firstCycle_display, targetSOC_multi, representativeDC_multi));
        saveas(fig_multi, savePath_multi);
        fprintf('  Saved 8-Channel Comparison: %s\n', savePath_multi);
        close(fig_multi);
    end
end

% SOC별 색상 (90: 파랑, 70: 초록, 50: 빨강)
socColors = containers.Map({'SOC90', 'SOC70', 'SOC50'}, ...
                           {[0 0 1], [0 0.8 0], [1 0 0]});
socLabels = containers.Map({'SOC90', 'SOC70', 'SOC50'}, ...
                          {'SOC 90', 'SOC 70', 'SOC 50'});

% 사이클별 색상 (lines 함수 사용)
cycleColors = lines(10);  % 최대 10개 사이클 지원

%% 채널 그룹별로 처리 (ch9, ch10, ... ch16)
for chGroupIdx = 1:length(channelGroupKeys)
    chGroupKey = channelGroupKeys{chGroupIdx};
    channelList = channelGroups(chGroupKey);
    
    fprintf('\n=== Processing Channel Group: %s (%d entries) ===\n', chGroupKey, length(channelList));
    
    % 해당 채널 그룹의 모든 사이클 수집
    allCyclesForGroup = {};
    for chIdx = 1:length(channelList)
        chName = channelList{chIdx};
        if isfield(allData, chName)
            cycles = fieldnames(allData.(chName));
            for c = 1:length(cycles)
                cycle_valid = cycles{c};
                if isKey(cycleNameReverseMap, cycle_valid)
                    cycle_display = cycleNameReverseMap(cycle_valid);
                else
                    cycle_display = cycle_valid;
                end
                if ~ismember(cycle_display, allCyclesForGroup)
                    allCyclesForGroup{end+1} = cycle_display;
                end
            end
        end
    end
    
    % 사이클 정렬
    cycle_nums = zeros(length(allCyclesForGroup), 1);
    for c = 1:length(allCyclesForGroup)
        cyc_str = allCyclesForGroup{c};
        num_match = regexp(cyc_str, '(\d+)cyc', 'tokens');
        if ~isempty(num_match)
            cycle_nums(c) = str2double(num_match{1}{1});
        else
            cycle_nums(c) = 999;
        end
    end
    [~, cyc_sort_idx] = sort(cycle_nums);
    allCyclesForGroup = allCyclesForGroup(cyc_sort_idx);
    
    fprintf('  Found cycles: %s\n', strjoin(allCyclesForGroup, ', '));
    
    %% 방법 A: SOC별 Dynamics 비교 (모든 사이클을 하나의 플롯에)
    % 모든 사이클의 데이터를 수집하여 하나의 플롯에 표시
    if ~isempty(allCyclesForGroup)
        fprintf('  Creating SOC Comparison plot (all cycles)...\n');
        
        % 사용 가능한 SOC 확인 (첫 번째 채널에서)
        availableSOCs = {};
        for chIdx = 1:length(channelList)
            chName = channelList{chIdx};
            if isfield(allData, chName)
                cycles = fieldnames(allData.(chName));
                if ~isempty(cycles)
                    cycle_valid = cycles{1};
                    if isKey(cycleNameReverseMap, cycle_valid)
                        cycle_display = cycleNameReverseMap(cycle_valid);
                    else
                        cycle_display = cycle_valid;
                    end
                    socs = fieldnames(allData.(chName).(cycle_valid));
                    availableSOCs = intersect(targetSOCs, socs);
                    break;
                end
            end
        end
        
        if ~isempty(availableSOCs)
            fprintf('    Found %d SOC levels: %s\n', length(availableSOCs), strjoin(availableSOCs, ', '));
            
            % 대표 주행 프로파일만 시각화 (DC2)
            representativeDCs = {'DC2'};
            profileOrder = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};
            
            % 각 SOC별로 플롯 생성 (SOC별 서브플롯, 각 서브플롯에 모든 사이클 표시)
            numSOCs = length(availableSOCs);
            numCols = min(3, numSOCs);
            numRows = ceil(numSOCs / numCols);
            
            fig1 = figure('Position', [100, 100, 1800, 1200]);
            sgtitle(sprintf('V-I Hysteresis Loop (Normalized): %s - All Cycles (DC2)', chGroupKey), ...
                    'FontSize', 16, 'FontWeight', 'bold');
            
            % 수치 데이터 수집용 (모든 사이클, 모든 SOC, 모든 DC)
            summaryData = struct();
            analysisResults = struct();
            
            % 각 SOC별 서브플롯
            for socIdx = 1:length(availableSOCs)
                socName = availableSOCs{socIdx};
                
                subplot(numRows, numCols, socIdx);
                hold on;
                
                legendEntries = {};
                
                % 모든 사이클에 대해 플롯
                for cycIdx = 1:length(allCyclesForGroup)
                    cycle_display = allCyclesForGroup{cycIdx};
                    
                    % 사이클의 유효한 필드명 찾기
                    if isKey(cycleNameMap, cycle_display)
                        cycle_valid = cycleNameMap(cycle_display);
                    else
                        if ~isempty(regexp(cycle_display, '^\d', 'once'))
                            cycle_valid = ['cyc_' cycle_display];
                        else
                            cycle_valid = cycle_display;
                        end
                        cycle_valid = regexprep(cycle_valid, '[^a-zA-Z0-9_]', '_');
                    end
                    
                    % 해당 사이클의 채널 찾기
                    foundData = false;
                    for chIdx = 1:length(channelList)
                        chName = channelList{chIdx};
                        if isfield(allData, chName) && ...
                           isfield(allData.(chName), cycle_valid) && ...
                           isfield(allData.(chName).(cycle_valid), socName) && ...
                           isfield(allData.(chName).(cycle_valid).(socName), 'DC2')
                            
                            data = allData.(chName).(cycle_valid).(socName).DC2;
                            I = data.I;
                            deltaV = data.deltaV;
                            t_seconds = data.t;
                            foundData = true;
                            break;
                        end
                    end
                    
                    if ~foundData
                        continue;
                    end
                    
                    % 분석 통계 계산 (모든 DC에 대해)
                    valid_idx = abs(I) > 0.1;
                    if sum(valid_idx) > 10
                        p = polyfit(I(valid_idx), deltaV(valid_idx), 1);
                        slope = p(1);
                    else
                        slope = NaN;
                    end
                    
                    if length(I) > 2
                        area = abs(polyarea(I, deltaV));
                    else
                        area = NaN;
                    end
                    
                    charge_idx = I > 0.1;
                    discharge_idx = I < -0.1;
                    if sum(charge_idx) > 5 && sum(discharge_idx) > 5
                        charge_deltaV_range = max(deltaV(charge_idx)) - min(deltaV(charge_idx));
                        discharge_deltaV_range = max(deltaV(discharge_idx)) - min(deltaV(discharge_idx));
                        asymmetry = abs(charge_deltaV_range - discharge_deltaV_range) / ...
                                    max(charge_deltaV_range, discharge_deltaV_range);
                    else
                        asymmetry = NaN;
                    end
                    
                    % 모든 DC에 대해 수치 데이터 저장 (사이클별, SOC별)
                    % cycle_display를 유효한 필드명으로 변환
                    cycle_field = cycle_display;
                    if ~isempty(regexp(cycle_field, '^\d', 'once'))
                        cycle_field = ['cyc_' cycle_field];
                    end
                    cycle_field = regexprep(cycle_field, '[^a-zA-Z0-9_]', '_');
                    
                    for dcIdx = 1:length(profileOrder)
                        dcName = profileOrder{dcIdx};
                        if ~isfield(summaryData, dcName)
                            summaryData.(dcName) = struct();
                        end
                        if ~isfield(summaryData.(dcName), socName)
                            summaryData.(dcName).(socName) = struct();
                        end
                        if ~isfield(summaryData.(dcName).(socName), cycle_field)
                            summaryData.(dcName).(socName).(cycle_field) = struct();
                        end
                        summaryData.(dcName).(socName).(cycle_field).slope = slope;
                        summaryData.(dcName).(socName).(cycle_field).area = area;
                        summaryData.(dcName).(socName).(cycle_field).asymmetry = asymmetry;
                        % 원본 cycle_display도 저장 (CSV 출력용)
                        summaryData.(dcName).(socName).(cycle_field).cycle_display = cycle_display;
                    end
                    
                    % 분석 결과 저장 (MAT 파일용)
                    for dcIdx = 1:length(profileOrder)
                        dcName = profileOrder{dcIdx};
                        resultKey = sprintf('%s_%s_%s', cycle_display, socName, dcName);
                        if ~isempty(regexp(resultKey, '^\d', 'once'))
                            resultKey = ['cyc_' resultKey];
                        end
                        resultKey = regexprep(resultKey, '[^a-zA-Z0-9_]', '_');
                        analysisResults.(resultKey) = struct();
                        analysisResults.(resultKey).slope = slope;
                        analysisResults.(resultKey).area = area;
                        analysisResults.(resultKey).asymmetry = asymmetry;
                        analysisResults.(resultKey).channel = chGroupKey;
                        analysisResults.(resultKey).cycle = cycle_display;
                        analysisResults.(resultKey).SOC = socName;
                        analysisResults.(resultKey).profile = dcName;
                    end
                    
                    % DC2만 시각화 (모든 사이클을 색상 다르게)
                    % 데이터 샘플링
                    max_points = 5000;
                    if length(I) > max_points
                        sample_idx = round(linspace(1, length(I), max_points));
                        I_plot = I(sample_idx);
                        deltaV_plot = deltaV(sample_idx);
                        t_plot = t_seconds(sample_idx);
                    else
                        I_plot = I;
                        deltaV_plot = deltaV;
                        t_plot = t_seconds;
                    end
                    
                    % 사이클별 색상 사용
                    color_cycle = cycleColors(mod(cycIdx-1, size(cycleColors,1))+1, :);
                    
                    % 시간에 따른 색상 매핑으로 산점도 플롯
                    scatter(I_plot, deltaV_plot, 20, t_plot, 'filled');
                    
                    % 시간 순서로 선 연결 (사이클별 색상)
                    plot(I_plot, deltaV_plot, 'Color', color_cycle, 'LineWidth', 1.5, ...
                         'DisplayName', cycle_display);
                    
                    legendEntries{end+1} = cycle_display;
                end
                
                % 서브플롯 완성
                xlabel('Current (A)', 'FontSize', 11);
                ylabel('\Delta V (V - V_{OCV})', 'FontSize', 11);
                title(sprintf('%s - %s', chGroupKey, socName), 'FontSize', 12, 'FontWeight', 'bold');
                grid on;
                if ~isempty(legendEntries)
                    legend(legendEntries, 'Location', 'best', 'FontSize', 8);
                end
                
                % 색상바 추가 (시간 정보)
                colorbar;
                colormap('jet');
                c = colorbar;
                c.Label.String = 'Time (s)';
                
                hold off;
            end
            
            fprintf('    Saving SOC Comparison plot...\n');
            savePath1 = fullfile(outputDir, sprintf('VI_Hysteresis_SOC_Comparison_%s_AllCycles.fig', chGroupKey));
            saveas(fig1, savePath1);
            fprintf('  Saved SOC Comparison: %s\n', savePath1);
            close(fig1);
            
            % CSV 파일 생성 (모든 사이클, 모든 SOC, 모든 DC)
            fprintf('    Creating summary statistics CSV...\n');
            csvPath = fullfile(outputDir, sprintf('Summary_Statistics_%s_AllCycles.csv', chGroupKey));
            fid = fopen(csvPath, 'w');
            fprintf(fid, 'Cycle,DC,SOC,Slope_Ohm,Area_VA,Asymmetry\n');
            for dcIdx = 1:length(profileOrder)
                dcName = profileOrder{dcIdx};
                if isfield(summaryData, dcName)
                    for socIdx = 1:length(availableSOCs)
                        socName = availableSOCs{socIdx};
                        if isfield(summaryData.(dcName), socName)
                            for cycIdx = 1:length(allCyclesForGroup)
                                cycle_display = allCyclesForGroup{cycIdx};
                                % 유효한 필드명으로 변환
                                cycle_field = cycle_display;
                                if ~isempty(regexp(cycle_field, '^\d', 'once'))
                                    cycle_field = ['cyc_' cycle_field];
                                end
                                cycle_field = regexprep(cycle_field, '[^a-zA-Z0-9_]', '_');
                                
                                if isfield(summaryData.(dcName).(socName), cycle_field)
                                    s = summaryData.(dcName).(socName).(cycle_field);
                                    fprintf(fid, '%s,%s,%s,%.6f,%.6f,%.6f\n', ...
                                        cycle_display, dcName, socName, s.slope, s.area, s.asymmetry);
                                end
                            end
                        end
                    end
                end
            end
            fclose(fid);
            fprintf('  Saved Summary Statistics: %s\n', csvPath);
            
            % MAT 파일 저장 (모든 사이클, 모든 SOC, 모든 DC)
            fprintf('    Saving analysis results MAT...\n');
            savePath_analysis = fullfile(outputDir, sprintf('Analysis_Results_%s_AllCycles.mat', chGroupKey));
            save(savePath_analysis, 'analysisResults');
            fprintf('  Saved Analysis Results: %s\n', savePath_analysis);
        end
    end
    
    fprintf('  Completed processing channel group: %s\n', chGroupKey);
end

fprintf('\n=== Visualization Complete ===\n');
fprintf('Output directory: %s\n', outputDir);
