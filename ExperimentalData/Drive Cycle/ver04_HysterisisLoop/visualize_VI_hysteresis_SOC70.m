%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-I 궤적(Hysteresis Loop) 시각화 스크립트
% 목적: SOC70 데이터에 대해 채널별로 V-I 궤적을 시각화
%       가속/감속이 빈번한 구간에서 타원의 두께를 비교
% 
% x축: 전류(I), y축: 전압(V)
% Charge Transfer 지배적: 직선에 가까운 궤적
% Diffusion 심함: 동일 전류값에서 전압이 이력(Hysteresis)을 그리며 타원형
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 설정
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
targetSOC = 'SOC70';  % SOC70만 분석
outputDir = fullfile(pwd, 'VI_Hysteresis_Plots');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% 데이터 파일 찾기
matFiles = dir(fullfile(dataDir, 'parsedDriveCycle_*cyc_filtered.mat'));
if isempty(matFiles)
    error('No data files found in %s', dataDir);
end

fprintf('=== V-I Hysteresis Loop Visualization (SOC70) ===\n');
fprintf('Data directory: %s\n', dataDir);
fprintf('Target SOC: %s\n', targetSOC);
fprintf('Found %d cycle files\n', length(matFiles));

%% 모든 사이클 파일에서 데이터 수집
allData = struct();
cycleNameMap = containers.Map();  % 원본 사이클명 -> 유효한 필드명 매핑
cycleNameReverseMap = containers.Map();  % 유효한 필드명 -> 원본 사이클명 매핑

for i = 1:length(matFiles)
    fileName = matFiles(i).name;
    token = regexp(fileName, 'parsedDriveCycle_(\d+cyc)_filtered', 'tokens');
    if isempty(token), continue; end
    cycleType = token{1}{1};
    
    % 유효한 필드명으로 변환 (숫자로 시작하면 'cyc_' 접두사 추가)
    if ~isempty(regexp(cycleType, '^\d', 'once'))
        cycleType_valid = ['cyc_' cycleType];
    else
        cycleType_valid = cycleType;
    end
    % 특수문자를 언더스코어로 변환
    cycleType_valid = regexprep(cycleType_valid, '[^a-zA-Z0-9_]', '_');
    
    cycleNameMap(cycleType) = cycleType_valid;
    cycleNameReverseMap(cycleType_valid) = cycleType;
    
    fprintf('\n--- Loading %s (field: %s) ---\n', cycleType, cycleType_valid);
    dataStruct = load(fullfile(dataDir, fileName));
    
    % 데이터 변수명 찾기
    varName = sprintf('parsedDriveCycle_%s', cycleType);
    if ~isfield(dataStruct, varName)
        fprintf('  WARNING: Variable %s not found. Skipping.\n', varName);
        continue;
    end
    data_var = dataStruct.(varName);
    
    % 채널 목록
    channels = fieldnames(data_var);
    fprintf('  Channels found: %s\n', strjoin(channels, ', '));
    
    % 각 채널별로 SOC70 데이터 수집
    for chIdx = 1:length(channels)
        chName = channels{chIdx};
        
        if ~isfield(data_var, chName) || ~isstruct(data_var.(chName))
            continue;
        end
        
        % SOC70 확인 (다양한 형식 지원: SOC70, SOC_70 등)
        socFields = fieldnames(data_var.(chName));
        targetSOC_found = false;
        targetSOC_field = '';
        
        % SOC 필드명 찾기 (대소문자 무시, 언더스코어 유무 무시)
        for sIdx = 1:length(socFields)
            socField = socFields{sIdx};
            socField_upper = upper(socField);
            targetSOC_upper = upper(targetSOC);
            % SOC70, SOC_70, SOC70_ 등 다양한 형식 매칭
            if contains(socField_upper, targetSOC_upper) || ...
               (contains(socField_upper, 'SOC') && contains(socField_upper, '70'))
                targetSOC_found = true;
                targetSOC_field = socField;
                break;
            end
        end
        
        if ~targetSOC_found
            continue;
        end
        
        % Profile 필드 찾기
        socData = data_var.(chName).(targetSOC_field);
        if ~isstruct(socData)
            continue;
        end
        
        % Profile이 별도 구조체인지, 직접 필드인지 확인
        if isfield(socData, 'Profile') && isstruct(socData.Profile)
            profs = fieldnames(socData.Profile);
            useProfileStruct = true;
        else
            % DC1, DC2 등으로 시작하는 필드 찾기
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
            
            % Profile 데이터 가져오기
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
            
            % 데이터 저장 (유효한 필드명 사용)
            if ~isfield(allData, chName)
                allData.(chName) = struct();
            end
            
            if ~isfield(allData.(chName), cycleType_valid)
                allData.(chName).(cycleType_valid) = struct();
            end
            
            allData.(chName).(cycleType_valid).(profName) = struct();
            allData.(chName).(cycleType_valid).(profName).I = profData.I;
            allData.(chName).(cycleType_valid).(profName).V = profData.V;
            allData.(chName).(cycleType_valid).(profName).t = profData.t;
            
            fprintf('    Saved: %s - %s - %s\n', chName, cycleType, profName);
        end
    end
end

%% 채널별로 V-I 궤적 시각화
channels = fieldnames(allData);
fprintf('\n=== Visualizing V-I Hysteresis Loops ===\n');
fprintf('Total channels: %d\n', length(channels));

% Profile 이름 정렬 (일관된 순서를 위해)
profileOrder = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

for chIdx = 1:length(channels)
    chName = channels{chIdx};
    fprintf('\n--- Processing Channel: %s ---\n', chName);
    
    % 해당 채널의 사이클 목록
    cycles = fieldnames(allData.(chName));
    
    % 사이클 정렬 (숫자 순서대로)
    cycle_nums = zeros(length(cycles), 1);
    for c = 1:length(cycles)
        cyc_str = cycles{c};
        % cyc_ 접두사 제거하여 숫자 추출
        if startsWith(cyc_str, 'cyc_')
            cyc_str_clean = cyc_str(5:end);
        else
            cyc_str_clean = cyc_str;
        end
        num_match = regexp(cyc_str_clean, '(\d+)cyc', 'tokens');
        if ~isempty(num_match)
            cycle_nums(c) = str2double(num_match{1}{1});
        else
            cycle_nums(c) = 999;
        end
    end
    [~, cyc_sort_idx] = sort(cycle_nums);
    cycles = cycles(cyc_sort_idx);
    
    % Figure 생성 (각 채널당 하나)
    fig = figure('Position', [100, 100, 1600, 1200]);
    sgtitle(sprintf('V-I Hysteresis Loop: %s (%s)', chName, targetSOC), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % 각 사이클별로 서브플롯 생성
    numCycles = length(cycles);
    numCols = min(3, numCycles);
    numRows = ceil(numCycles / numCols);
    
    for cycIdx = 1:numCycles
        cycleType_valid = cycles{cycIdx};
        % 원본 사이클명 가져오기 (표시용)
        if isKey(cycleNameReverseMap, cycleType_valid)
            cycleType_display = cycleNameReverseMap(cycleType_valid);
        else
            cycleType_display = cycleType_valid;
        end
        
        subplot(numRows, numCols, cycIdx);
        hold on;
        
        % Profile별로 플롯
        profs = fieldnames(allData.(chName).(cycleType_valid));
            
            % Profile 순서 정렬
            profs_sorted = {};
            for p = 1:length(profileOrder)
                if ismember(profileOrder{p}, profs)
                    profs_sorted{end+1} = profileOrder{p};
                end
            end
            % 정렬되지 않은 Profile 추가
            for p = 1:length(profs)
                if ~ismember(profs{p}, profs_sorted)
                    profs_sorted{end+1} = profs{p};
                end
            end
            
            colors = lines(length(profs_sorted));
            legendEntries = {};
            
            for pIdx = 1:length(profs_sorted)
                profName = profs_sorted{pIdx};
                
                if ~isfield(allData.(chName).(cycleType_valid), profName)
                    continue;
                end
                
                data = allData.(chName).(cycleType_valid).(profName);
                I = data.I;
                V = data.V;
                t = data.t;
                
                % 시간을 duration에서 초로 변환
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
                
                % 전류 변화율 계산 (가속/감속 구간 감지)
                if length(I) > 1 && length(t_seconds) > 1
                    dt_vec = [1; diff(t_seconds)];
                    dt_vec(dt_vec == 0) = 1;  % 0으로 나누기 방지
                    dI_dt = [0; diff(I) ./ dt_vec(2:end)];
                    if std(dI_dt) > 0
                        dI_dt(abs(dI_dt) > 10*std(dI_dt)) = 0;  % 이상치 제거
                        accel_decel_idx = abs(dI_dt) > 0.5 * std(abs(dI_dt));
                    else
                        accel_decel_idx = false(size(I));
                    end
                else
                    accel_decel_idx = false(size(I));
                end
                
                % 전체 데이터 플롯 (연한 색상)
                scatter(I, V, 15, colors(pIdx, :), 'filled', ...
                        'DisplayName', [profName ' (all)']);
                
                % 가속/감속 구간 강조 (진한 색상, 큰 마커)
                if any(accel_decel_idx)
                    scatter(I(accel_decel_idx), V(accel_decel_idx), 40, ...
                            colors(pIdx, :), 'filled', 'MarkerEdgeColor', 'k', ...
                            'LineWidth', 1.5, ...
                            'DisplayName', [profName ' (accel/decel)']);
                end
                
                % 시간 순서로 선 연결 (타원형 궤적 강조)
                plot(I, V, 'Color', colors(pIdx, :), 'LineWidth', 1.0, ...
                     'LineStyle', '--');
                
                legendEntries{end+1} = profName;
            end
        
        xlabel('Current (A)', 'FontSize', 12);
        ylabel('Voltage (V)', 'FontSize', 12);
        title(sprintf('%s', cycleType_display), 'FontSize', 12);
        grid on;
        legend(legendEntries, 'Location', 'best', 'FontSize', 8);
        hold off;
    end
    
    % Figure 저장
    savePath = fullfile(outputDir, sprintf('VI_Hysteresis_%s_%s.fig', chName, targetSOC));
    saveas(fig, savePath);
    fprintf('  Saved: %s\n', savePath);
    
    % 각 Profile별로 개별 플롯도 생성 (더 자세한 분석용)
    for cycIdx = 1:numCycles
        cycleType_valid = cycles{cycIdx};
        % 원본 사이클명 가져오기 (표시용)
        if isKey(cycleNameReverseMap, cycleType_valid)
            cycleType_display = cycleNameReverseMap(cycleType_valid);
        else
            cycleType_display = cycleType_valid;
        end
        
        profs = fieldnames(allData.(chName).(cycleType_valid));
        
        % Profile별 개별 플롯
        fig2 = figure('Position', [100, 100, 1600, 1200]);
        sgtitle(sprintf('V-I Hysteresis Loop: %s - %s (%s)', chName, cycleType_display, targetSOC), ...
                'FontSize', 16, 'FontWeight', 'bold');
        
        numProfs = length(profs);
        numCols2 = min(4, numProfs);
        numRows2 = ceil(numProfs / numCols2);
        
        for pIdx = 1:numProfs
            profName = profs{pIdx};
            
            if ~isfield(allData.(chName).(cycleType_valid), profName)
                continue;
            end
            
            data = allData.(chName).(cycleType_valid).(profName);
            I = data.I;
            V = data.V;
            t = data.t;
            
            % 시간을 duration에서 초로 변환
            if isa(t, 'duration')
                t_seconds = seconds(t);
            else
                t_seconds = t;
            end
            
            % NaN 제거
            validIdx = ~isnan(I) & ~isnan(V);
            I = I(validIdx);
            V = V(validIdx);
            t_seconds = t_seconds(validIdx);
            
            if isempty(I) || isempty(V)
                continue;
            end
            
            subplot(numRows2, numCols2, pIdx);
            hold on;
            
            % 전류 변화율 계산 (가속/감속 구간 감지)
            if length(I) > 1 && length(t_seconds) > 1
                dt_vec = [1; diff(t_seconds)];
                dt_vec(dt_vec == 0) = 1;  % 0으로 나누기 방지
                dI_dt = [0; diff(I) ./ dt_vec(2:end)];
                if std(dI_dt) > 0
                    dI_dt(abs(dI_dt) > 10*std(dI_dt)) = 0;  % 이상치 제거
                    accel_decel_idx = abs(dI_dt) > 0.5 * std(abs(dI_dt));
                else
                    accel_decel_idx = false(size(I));
                end
            else
                accel_decel_idx = false(size(I));
            end
            
            % 전체 데이터 플롯 (시간에 따른 색상)
            scatter(I, V, 20, t_seconds, 'filled');
            
            % 가속/감속 구간 강조
            if any(accel_decel_idx)
                scatter(I(accel_decel_idx), V(accel_decel_idx), 50, ...
                        t_seconds(accel_decel_idx), 'filled', ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            end
            
            % 시간 순서로 선 연결 (타원형 궤적 강조)
            plot(I, V, 'k-', 'LineWidth', 0.5);
            
            colorbar;
            colormap('jet');
            c = colorbar;
            c.Label.String = 'Time (s)';
            
            xlabel('Current (A)', 'FontSize', 11);
            ylabel('Voltage (V)', 'FontSize', 11);
            title(sprintf('%s', profName), 'FontSize', 11);
            grid on;
            hold off;
            
            % 통계 정보 추가 (Hysteresis 두께 계산)
            I_range = max(I) - min(I);
            V_range = max(V) - min(V);
            
            % 동일 전류값에서의 전압 범위 계산 (Hysteresis 두께)
            I_bins = linspace(min(I), max(I), 50);
            V_hysteresis = zeros(size(I_bins));
            for binIdx = 1:length(I_bins)-1
                bin_mask = I >= I_bins(binIdx) & I < I_bins(binIdx+1);
                if sum(bin_mask) > 1
                    V_hysteresis(binIdx) = max(V(bin_mask)) - min(V(bin_mask));
                end
            end
            avg_hysteresis_thickness = mean(V_hysteresis(V_hysteresis > 0));
            max_hysteresis_thickness = max(V_hysteresis);
            
            text(0.05, 0.95, ...
                 sprintf('I range: %.2f A\nV range: %.3f V\nAvg Hyst: %.4f V\nMax Hyst: %.4f V', ...
                         I_range, V_range, avg_hysteresis_thickness, max_hysteresis_thickness), ...
                 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');
        end
        
        savePath2 = fullfile(outputDir, sprintf('VI_Hysteresis_%s_%s_%s_detail.fig', chName, cycleType_display, targetSOC));
        saveas(fig2, savePath2);
        close(fig2);
    end
    
    close(fig);
    
    % 8개 주행 부하(DC1-DC8)를 한눈에 비교하는 플롯 생성
    fig3 = figure('Position', [100, 100, 2000, 1000]);
    sgtitle(sprintf('V-I Hysteresis Loop Comparison: %s (%s) - All Drive Cycles', chName, targetSOC), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % 모든 사이클의 데이터를 수집하여 DC별로 비교
    allCycles = fieldnames(allData.(chName));
    
    % 사이클 정렬 (숫자 순서대로)
    cycle_nums = zeros(length(allCycles), 1);
    for c = 1:length(allCycles)
        cyc_str = allCycles{c};
        if startsWith(cyc_str, 'cyc_')
            cyc_str_clean = cyc_str(5:end);
        else
            cyc_str_clean = cyc_str;
        end
        num_match = regexp(cyc_str_clean, '(\d+)cyc', 'tokens');
        if ~isempty(num_match)
            cycle_nums(c) = str2double(num_match{1}{1});
        else
            cycle_nums(c) = 999;
        end
    end
    [~, cyc_sort_idx] = sort(cycle_nums);
    allCycles = allCycles(cyc_sort_idx);
    
    % DC1-DC8에 대해 서브플롯 생성 (2행 4열)
    for dcIdx = 1:length(profileOrder)
        dcName = profileOrder{dcIdx};
        
        subplot(2, 4, dcIdx);
        hold on;
        
        colors_cycle = lines(length(allCycles));
        legendEntries_dc = {};
        
        for cycIdx = 1:length(allCycles)
            cycleType_valid = allCycles{cycIdx};
            % 원본 사이클명 가져오기 (표시용)
            if isKey(cycleNameReverseMap, cycleType_valid)
                cycleType_display = cycleNameReverseMap(cycleType_valid);
            else
                cycleType_display = cycleType_valid;
            end
            
            if ~isfield(allData.(chName).(cycleType_valid), dcName)
                continue;
            end
            
            data = allData.(chName).(cycleType_valid).(dcName);
            I = data.I;
            V = data.V;
            t = data.t;
            
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
            
            % 전류 변화율 계산
            if length(I) > 1 && length(t_seconds) > 1
                dt_vec = [1; diff(t_seconds)];
                dt_vec(dt_vec == 0) = 1;  % 0으로 나누기 방지
                dI_dt = [0; diff(I) ./ dt_vec(2:end)];
                if std(dI_dt) > 0
                    dI_dt(abs(dI_dt) > 10*std(dI_dt)) = 0;
                    accel_decel_idx = abs(dI_dt) > 0.5 * std(abs(dI_dt));
                else
                    accel_decel_idx = false(size(I));
                end
            else
                accel_decel_idx = false(size(I));
            end
            
            % 전체 데이터 플롯
            scatter(I, V, 10, colors_cycle(cycIdx, :), 'filled');
            
            % 가속/감속 구간 강조
            if any(accel_decel_idx)
                scatter(I(accel_decel_idx), V(accel_decel_idx), 30, ...
                        colors_cycle(cycIdx, :), 'filled', ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
            end
            
            % 시간 순서로 선 연결
            plot(I, V, 'Color', colors_cycle(cycIdx, :), 'LineWidth', 1.0, ...
                 'LineStyle', '-');
            
            legendEntries_dc{end+1} = cycleType_display;
        end
        
        xlabel('Current (A)', 'FontSize', 11);
        ylabel('Voltage (V)', 'FontSize', 11);
        title(sprintf('%s', dcName), 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        if ~isempty(legendEntries_dc)
            legend(legendEntries_dc, 'Location', 'best', 'FontSize', 7);
        end
        hold off;
    end
    
    % DC 비교 플롯 저장
    savePath3 = fullfile(outputDir, sprintf('VI_Hysteresis_%s_%s_DC_Comparison.fig', chName, targetSOC));
    saveas(fig3, savePath3);
    fprintf('  Saved DC Comparison: %s\n', savePath3);
    close(fig3);
end

fprintf('\n=== Visualization Complete ===\n');
fprintf('Output directory: %s\n', outputDir);
