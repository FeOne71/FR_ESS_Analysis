clear; clc; close all;

% 1. 데이터 로드
% 현재 스크립트가 있는 폴더의 Results 폴더에서 데이터 로드
scriptDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(scriptDir, 'Results');
dataFile = fullfile(resultsDir, 'DriveCycle_Summary_Table.mat');

% 데이터 파일 존재 확인
if ~exist(dataFile, 'file')
    error('데이터 파일이 없습니다: %s\n\n먼저 다음 스크립트를 순서대로 실행하세요:\n1. DriveCycle_EventAnalysis_01.m (Feature 추출)\n2. DriveCycle_DataAggregation_02.m (테이블 생성)', dataFile);
end

% EventAnalysis 결과 파일 확인
eventFiles = dir(fullfile(resultsDir, 'Lab_DC_DCIR_*cyc_Events.mat'));
if ~isempty(eventFiles)
    % EventAnalysis 파일이 있으면 Feature 필드가 있는지 빠르게 확인
    fprintf('EventAnalysis 결과 파일 발견: %d개\n', length(eventFiles));
    fprintf('Feature 필드 존재 여부 확인 중...\n');
    
    % 첫 번째 파일에서 Feature 필드 확인
    try
        firstFile = fullfile(resultsDir, eventFiles(1).name);
        load(firstFile);
        structName = who('-file', firstFile);
        if ~isempty(structName)
            eval(sprintf('tempData = %s;', structName{1}));
            chFields = fieldnames(tempData);
            if ~isempty(chFields)
                % 첫 번째 채널, 첫 번째 SOC, 첫 번째 프로파일, 첫 번째 이벤트 확인
                firstCh = chFields{1};
                if isfield(tempData, firstCh)
                    chData = tempData.(firstCh);
                    socFields = fieldnames(chData);
                    if ~isempty(socFields)
                        socData = chData.(socFields{1});
                        profFields = fieldnames(socData);
                        if ~isempty(profFields)
                            profData = socData.(profFields{1});
                            evtFields = fieldnames(profData);
                            if ~isempty(evtFields)
                                evtData = profData.(evtFields{1});
                                allFields = fieldnames(evtData);
                                hasFeat = any(startsWith(allFields, 'Feat_'));
                                if ~hasFeat
                                    fprintf('\n⚠️  확인: EventAnalysis 결과 파일에 Feature 필드가 없습니다!\n');
                                    fprintf('   → DriveCycle_EventAnalysis_01.m을 다시 실행해야 합니다.\n\n');
                                else
                                    fprintf('✓ EventAnalysis 파일에 Feature 필드가 있습니다.\n');
                                    fprintf('   → DataAggregation만 다시 실행하면 됩니다.\n\n');
                                end
                            end
                        end
                    end
                end
            end
        end
        clear tempData;
    catch
        fprintf('  (Feature 확인 중 오류 발생, 무시하고 계속...)\n\n');
    end
end

fprintf('Loading data from: %s\n', dataFile);
load(dataFile);

%% ============================================================
%% [NEW] Advanced Outlier Removal (모든 Feature 적용)
%% ============================================================

% 1. dQ 변수 매핑 (Feat_Cap을 dQ로 사용)
if ismember('Feat_Cap', summaryTable.Properties.VariableNames)
    summaryTable.Feat_dQ = summaryTable.Feat_Cap; 
    fprintf('Mapped "Feat_Cap" to "Feat_dQ"\n');
end

% -------------------------------------------------------------
% Feature 1: dQ/dV Variance (가장 중요)
% -------------------------------------------------------------
% 현상: 전압 평탄 구간(Plateau)에서 분모가 0이 되어 무한대(Inf)나 수천 단위로 튐.
% 조치: 0~100 범위를 벗어나면 노이즈로 간주.
if ismember('Feat_Var_dQdV', summaryTable.Properties.VariableNames)
    raw = summaryTable.Feat_Var_dQdV;
    % Threshold: 100 (경험적 수치, 배터리 특성에 따라 조절 가능)
    mask = raw > 100 | isinf(raw) | isnan(raw); 
    summaryTable.Feat_Var_dQdV(mask) = NaN;
    fprintf('Cleaned Feat_Var_dQdV: Removed %d outliers\n', sum(mask));
end

% -------------------------------------------------------------
% Feature 2: Voltage Skewness (전압 왜도)
% -------------------------------------------------------------
% 현상: 정규분포(-3 ~ +3)를 크게 벗어나는 값은 계산 오류일 가능성 높음.
% 조치: 절댓값이 5 이상이면 제거.
if ismember('Feat_Skew_V', summaryTable.Properties.VariableNames)
    raw = summaryTable.Feat_Skew_V;
    mask = abs(raw) > 5 | isinf(raw);
    summaryTable.Feat_Skew_V(mask) = NaN;
    fprintf('Cleaned Feat_Skew_V: Removed %d outliers\n', sum(mask));
end

% -------------------------------------------------------------
% Feature 3: Voltage Slope (전압 기울기)
% -------------------------------------------------------------
% 현상: 충전 중 전압 상승(양수), 방전 중 전압 하강(음수)
%       순식간에 튀는 경우(너무 큰 절댓값)는 측정 오류.
% 조치: 절댓값이 0.1 V/s 이상이면 제거 (충전/방전 모두 지원)
if ismember('Feat_Slope_V', summaryTable.Properties.VariableNames)
    raw = summaryTable.Feat_Slope_V;
    % [FIX] 방전(Discharge)일 경우 기울기가 음수이므로, 절댓값으로 크기만 제한
    % 기울기가 너무 가파른 것(> 0.1 V/s)만 제거
    mask = abs(raw) > 0.1 | isinf(raw); 
    summaryTable.Feat_Slope_V(mask) = NaN;
    fprintf('Cleaned Feat_Slope_V: Removed %d outliers (using |slope| > 0.1)\n', sum(mask));
end

% -------------------------------------------------------------
% Feature 4: dQ (Capacity Throughput)
% -------------------------------------------------------------
% 현상: 노이즈로 인해 용량이 0에 가깝거나, 배터리 용량을 초과하는 경우.
% 조치: 최소값(0.001 Ah), 최대값(배터리 정격 용량의 20% 등) 설정.
if ismember('Feat_dQ', summaryTable.Properties.VariableNames)
    raw = summaryTable.Feat_dQ;
    % 세그먼트가 너무 짧아서 용량이 거의 없거나(0.001 미만), 
    % 부분 충전인데 너무 큰 값(예: 10Ah 이상)은 제거
    mask = raw < 0.001 | raw > 10.0; 
    summaryTable.Feat_dQ(mask) = NaN;
    fprintf('Cleaned Feat_dQ: Removed %d outliers\n', sum(mask));
end

% -------------------------------------------------------------
% Feature 5: dV/dQ Variance (Differential Voltage - 방전 데이터 맞춤)
% -------------------------------------------------------------
% 현상: dV/dQ는 저항과 유사한 특성을 가짐. Variance는 변동성을 나타냄.
% 조치: Variance가 비정상적으로 큰 값(예: 1000 이상)은 노이즈로 간주.
if ismember('Feat_Var_dVdQ', summaryTable.Properties.VariableNames)
    raw = summaryTable.Feat_Var_dVdQ;
    % Variance는 제곱값이므로 범위가 큼. 1000 이상이면 제거
    mask = raw > 1000 | isinf(raw) | isnan(raw);
    summaryTable.Feat_Var_dVdQ(mask) = NaN;
    fprintf('Cleaned Feat_Var_dVdQ: Removed %d outliers\n', sum(mask));
end

% -------------------------------------------------------------
% Feature 6: dV/dQ Mean (Differential Voltage Mean)
% -------------------------------------------------------------
% 현상: dV/dQ 평균값이 비정상적으로 큰 값은 계산 오류일 가능성.
% 조치: 사용자 코드에서 limit_val = 50을 사용했으므로, 비슷한 범위 적용.
if ismember('Feat_Mean_dVdQ', summaryTable.Properties.VariableNames)
    raw = summaryTable.Feat_Mean_dVdQ;
    mask = abs(raw) > 100 | isinf(raw) | isnan(raw);
    summaryTable.Feat_Mean_dVdQ(mask) = NaN;
    fprintf('Cleaned Feat_Mean_dVdQ: Removed %d outliers\n', sum(mask));
end

fprintf('\n');

% 1.5. 테이블에 있는 Feature 필드 확인
allColumns = summaryTable.Properties.VariableNames;
featColumns = allColumns(startsWith(allColumns, 'Feat_'));
fprintf('Available Feature columns (%d):\n', length(featColumns));
for i = 1:min(20, length(featColumns))
    fprintf('  %s\n', featColumns{i});
end
if length(featColumns) > 20
    fprintf('  ... and %d more\n', length(featColumns) - 20);
end
fprintf('\n');

% 2. 보고 싶은 Feature 정의 (논문 핵심 Feature)
targetFeatures = {'Feat_Var_dQdV', 'Feat_Var_dVdQ', 'Feat_Skew_V', 'Feat_Mean_I', 'Feat_Slope_V'};
targetLabels = {'dQ/dV Variance', 'dV/dQ Variance', 'Voltage Skewness', 'Mean Current', 'Voltage Slope'};

% 사용 가능한 Feature만 필터링
availableFeatures = {};
availableLabels = {};
for i = 1:length(targetFeatures)
    if ismember(targetFeatures{i}, allColumns)
        availableFeatures{end+1} = targetFeatures{i};
        availableLabels{end+1} = targetLabels{i};
    else
        fprintf('WARNING: Feature ''%s'' not found in table. Skipping...\n', targetFeatures{i});
    end
end

if isempty(availableFeatures)
    fprintf('\n========================================\n');
    fprintf('ERROR: Feature 필드가 테이블에 없습니다!\n');
    fprintf('========================================\n\n');
    fprintf('원인: DriveCycle_Summary_Table.mat가 Feature 추출 이전에 생성되었습니다.\n\n');
    fprintf('해결 방법:\n');
    fprintf('1. DriveCycle_EventAnalysis_01.m을 실행하여 Feature를 추출하세요\n');
    fprintf('   (Feature 추출 코드가 추가된 후로는 처음 실행하는 것임)\n');
    fprintf('2. DriveCycle_DataAggregation_02.m을 다시 실행하여\n');
    fprintf('   Feature가 포함된 새로운 테이블을 생성하세요\n\n');
    fprintf('참고: EventAnalysis 결과 파일은 Results 폴더에 있습니다.\n');
    fprintf('      하지만 Feature 추출 코드 추가 이전에 생성된 파일들입니다.\n\n');
    error('Feature 필드가 없습니다. 위 안내를 따라 EventAnalysis와 DataAggregation을 다시 실행하세요.');
end

fprintf('Using %d out of %d requested features:\n', length(availableFeatures), length(targetFeatures));
for i = 1:length(availableFeatures)
    fprintf('  - %s\n', availableLabels{i});
end
fprintf('\n');

targetFeatures = availableFeatures;
targetLabels = availableLabels;

% === 설정: 분석할 이벤트 타입 선택 ===
% 'Charge' 또는 'Discharge' 또는 {'Charge', 'Discharge'} (둘 다)
eventTypesToPlot = {'Discharge'};  % 원하는 이벤트 타입으로 변경하세요
% eventTypesToPlot = {'Discharge'};  % 방전만 보고 싶으면 이 줄의 주석을 제거하고 위 줄을 주석 처리
% eventTypesToPlot = {'Charge', 'Discharge'};  % 둘 다 보고 싶으면 이 줄의 주석을 제거

fprintf('========================================\n');
fprintf('Feature Statistics & Visualization\n');
fprintf('========================================\n');
fprintf('Event Types to Plot: %s\n\n', strjoin(eventTypesToPlot, ', '));

% 3. 채널별로 나누어 그리기 (각 채널마다 figure 생성)
uniqueSOCs = unique(summaryTable.SOC(~isnan(summaryTable.SOC)));
uniqueSOCs = sort(uniqueSOCs, 'descend'); % SOC90, SOC70, SOC50 순서로 정렬
uniqueChannels = sort(unique(summaryTable.Channel(~isnan(summaryTable.Channel))));

for eventTypeIdx = 1:length(eventTypesToPlot)
    eventType = eventTypesToPlot{eventTypeIdx};
    
    fprintf('========================================\n');
    fprintf('Analyzing %s Events\n', eventType);
    fprintf('========================================\n\n');
    
    % 채널별로 figure 생성
    for chIdx = 1:length(uniqueChannels)
        channel = uniqueChannels(chIdx);
        
        fprintf('========================================\n');
        fprintf('Channel %d Analysis (%s)\n', channel, eventType);
        fprintf('========================================\n\n');
        
        % 채널별 figure 생성 (SOC별 subplot: 3행, Feature별 subplot: 2x2)
        figure('Name', sprintf('Channel %d - Features vs Cycle (%s)', channel, eventType), ...
            'Position', [100 + (chIdx-1)*50, 100 + (chIdx-1)*50, 1600, 1200]);
        
        % 각 SOC별 통계 정보 저장
        allStatsData = struct();
        
        % SOC별로 subplot 생성 (3개: SOC90, SOC70, SOC50)
        for s = 1:length(uniqueSOCs)
            soc = uniqueSOCs(s);
            
            % 해당 채널, SOC 및 이벤트 타입 데이터 필터링
            idx = summaryTable.Channel == channel & summaryTable.SOC == soc & strcmp(summaryTable.EventType, eventType);
            dataSOC = summaryTable(idx, :);
            
            % Feature별 subplot 레이아웃 (동적 계산)
            numFeatures = length(targetFeatures);
            numSOCs = length(uniqueSOCs);
            featCols = numFeatures; % Feature 개수만큼 열 생성
            
            % SOC subplot 위치 계산 (3행 1열 레이아웃)
            socSubplotIdx = s; % 1, 2, 3 (SOC90, SOC70, SOC50)
            
            % 통계 정보 저장을 위한 구조체 초기화
            statsData = struct();
            
            % 각 Feature별 subplot 생성
            for f = 1:length(targetFeatures)
                featName = targetFeatures{f};
                
                % 전체 subplot 위치: (SOC 행, Feature 열) -> numSOCs행 numFeatures열 그리드
                % SOC는 행(1-3), Feature는 열(1-numFeatures)
                subplotIdx = (s - 1) * numFeatures + f; % SOC별로 numFeatures개씩 배치
                subplot(numSOCs, numFeatures, subplotIdx);
                hold on;
        
                % NaN 값을 제거
                validIdx = ~isnan(dataSOC.(featName)) & ~isnan(dataSOC.Cycle);
                if sum(validIdx) < 5
                    text(0.5, 0.5, sprintf('Insufficient data\n(N=%d)', sum(validIdx)), ...
                        'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 10);
                    xlabel('Cycle');
                    ylabel(targetLabels{f});
                    title(sprintf('%s (SOC %d%%)', targetLabels{f}, soc));
                    grid on;
                    continue;
                end
        
                cycles_valid = dataSOC.Cycle(validIdx);
                featValues_valid = dataSOC.(featName)(validIdx);
                
                % 해당 채널 데이터만 플로팅 (채널별 figure이므로 하나의 선/점만)
                try
                    plot(cycles_valid, featValues_valid, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
                catch ME
                    fprintf('Error plotting %s: %s\n', featName, ME.message);
                    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
                        'HorizontalAlignment', 'center', 'Units', 'normalized');
                    continue;
                end
        
        % === 통계 분석 ===
        % 1. 기본 통계
        feat_mean = mean(featValues_valid);
        feat_std = std(featValues_valid);
        feat_median = median(featValues_valid);
        feat_min = min(featValues_valid);
        feat_max = max(featValues_valid);
        feat_cv = feat_std / abs(feat_mean) * 100; % Coefficient of Variation (%)
        
        % 2. Cycle과 Feature 간 상관관계 분석
        [R, P] = corrcoef(cycles_valid, featValues_valid);
        corr_coef = R(1,2);
        p_value = P(1,2);
        
        % 3. 선형 추세 분석
        p = polyfit(cycles_valid, featValues_valid, 1);
        slope = p(1);
        intercept = p(2);
        y_fit = polyval(p, cycles_valid);
        ss_res = sum((featValues_valid - y_fit).^2);
        ss_tot = sum((featValues_valid - mean(featValues_valid)).^2);
        r_squared = 1 - (ss_res / ss_tot);
        
        % 4. 추세선 그리기
        cycle_sorted = sort(unique(cycles_valid));
        if length(cycle_sorted) > 1
            y_trend = polyval(p, cycle_sorted);
            plot(cycle_sorted, y_trend, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Linear Trend');
        end
        
        % 5. Cycle별 평균 및 표준편차
        uniqueCycles = sort(unique(cycles_valid));
        cycle_means = zeros(size(uniqueCycles));
        cycle_stds = zeros(size(uniqueCycles));
        for c = 1:length(uniqueCycles)
            cycleIdx = cycles_valid == uniqueCycles(c);
            cycle_means(c) = mean(featValues_valid(cycleIdx));
            cycle_stds(c) = std(featValues_valid(cycleIdx));
        end
        errorbar(uniqueCycles, cycle_means, cycle_stds, 'ko-', 'LineWidth', 1.5, ...
            'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Cycle Mean±SD');
        
        % 6. 통계 정보 텍스트로 표시
        stats_text = sprintf(['Mean: %.3g (SD: %.3g, CV: %.1f%%)\n' ...
            'Range: [%.3g, %.3g], Median: %.3g\n' ...
            'Correlation: r=%.3f (p=%.3f)\n' ...
            'Trend: slope=%.3g/cycle, R²=%.3f'], ...
            feat_mean, feat_std, feat_cv, feat_min, feat_max, feat_median, ...
            corr_coef, p_value, slope, r_squared);
        
        % 텍스트 박스 위치 (그래프 우측 상단)
        xlim_current = xlim;
        ylim_current = ylim;
        text_x = xlim_current(1) + 0.02 * (xlim_current(2) - xlim_current(1));
        text_y = ylim_current(2) - 0.05 * (ylim_current(2) - ylim_current(1));
        text(text_x, text_y, stats_text, ...
            'FontSize', 8, 'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', ...
            'Margin', 3);
        
        xlabel('Cycle', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel(targetLabels{f}, 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s vs Cycle (SOC %d%%)', targetLabels{f}, soc), ...
            'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        legend('Location', 'best', 'FontSize', 8);
        
        % 통계 정보 저장 (콘솔 출력용)
        statsData.(featName).mean = feat_mean;
        statsData.(featName).std = feat_std;
        statsData.(featName).median = feat_median;
        statsData.(featName).min = feat_min;
        statsData.(featName).max = feat_max;
        statsData.(featName).cv = feat_cv;
        statsData.(featName).correlation = corr_coef;
        statsData.(featName).p_value = p_value;
        statsData.(featName).slope = slope;
        statsData.(featName).r_squared = r_squared;
        statsData.(featName).N = sum(validIdx);
    end
    
        sgtitle(sprintf('Channel %d - Feature Evolution Analysis (SOC %d%%, %s)', channel, soc, eventType), ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        % === 콘솔에 상세 통계 출력 ===
        fprintf('\n--- Feature Statistics (Channel %d, SOC %d%%, %s) ---\n', channel, soc, eventType);
        fprintf('%-20s | %8s | %8s | %8s | %10s | %10s | %8s | %8s | %10s\n', ...
            'Feature', 'Mean', 'SD', 'CV(%)', 'Corr(r)', 'p-value', 'Slope', 'R²', 'N');
        fprintf('%s\n', repmat('-', 1, 120));
        
        for f = 1:length(targetFeatures)
            featName = targetFeatures{f};
            if isfield(statsData, featName)
                s = statsData.(featName);
                fprintf('%-20s | %8.3g | %8.3g | %8.1f | %10.3f | %10.3f | %8.3g | %8.3f | %10d\n', ...
                    targetLabels{f}, s.mean, s.std, s.cv, s.correlation, s.p_value, s.slope, s.r_squared, s.N);
            end
        end
        fprintf('\n');
        
        % === 해석 가이드 출력 ===
        fprintf('--- Interpretation Guide ---\n');
        fprintf('Correlation (r): -1 to +1, closer to ±1 means stronger linear relationship\n');
        fprintf('p-value: < 0.05 indicates statistically significant correlation\n');
        fprintf('Slope: Change in feature per cycle (positive = increasing, negative = decreasing)\n');
        fprintf('R²: Proportion of variance explained by linear trend (0-1, higher is better)\n');
        fprintf('CV%%: Coefficient of Variation = (SD/Mean)×100, measures variability\n');
        fprintf('\n');
        end  % End of SOC loop
    end  % End of channel loop
end  % End of eventType loop

fprintf('========================================\n');
fprintf('Visualization Complete\n');
fprintf('========================================\n\n');
