%% Field Module R50s Visualization - Generate Visualizations from Saved Data
% This script loads saved resistance calculation results and generates visualizations
% 
% Process:
% 1. Load saved resistance data from MAT file
% 2. Generate 3D scatter plots
% 3. Generate 3D bar plots
% 4. Generate monthly dV/dI scatter plots
% 5. Generate module-wise boxplot
% 6. Save figures to Module Resistance folder

clear; clc; close all;

%% Configuration
fontSize = 12;
saveDir = 'Module Resistance';

% Topology
Ns_cells_per_module = 14;  % 14 cells in series per module
numModules = 17;  % Module01 to Module17

% Target Rack and Module for visualization (can be modified)
targetRack = 'Rack01';
targetModule = 1:17;  % Module number (1-17) - Visualize all modules

%% Load saved calculation results
fprintf('=== Loading Saved Resistance Data ===\n');
dataFile = fullfile(saveDir, sprintf('%s_Resistance_Summary.mat', targetRack));

if ~exist(dataFile, 'file')
    error('Calculation results file not found: %s\nPlease run Field_Module_R_Calculation.m first.', dataFile);
end

load(dataFile, 'module_data');
fprintf('Data loaded from: %s\n', dataFile);

%% Create Rack-specific folder for figures
rackSaveDir = fullfile(saveDir, targetRack);
if ~exist(rackSaveDir, 'dir')
    mkdir(rackSaveDir);
end

%% Process each target module for visualization
if isscalar(targetModule)
    targetModules = targetModule;
else
    targetModules = targetModule;
end

%% Step 1: 모든 모듈의 저항값 수집하여 최대/최소 찾기
fprintf('=== Collecting all resistance values to determine global Y-axis range ===\n');
all_resistance_values = [];

for modIdx = 1:length(targetModules)
    modNum = targetModules(modIdx);
    moduleName = sprintf('Module%02d', modNum);
    
    % Check if module has any data
    hasData = false;
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if isfield(module_data, moduleName) && isfield(module_data.(moduleName), cellName) && ...
           ~isempty(module_data.(moduleName).(cellName).dates)
            hasData = true;
            break;
        end
    end
    
    if ~hasData
        continue;
    end
    
    % 전압 구간별 박스플롯용 데이터 수집
    if isfield(module_data, moduleName)
        v_bins = 3.2:0.05:4.2;
        num_bins = length(v_bins) - 1;
        years_list = [2021, 2022, 2023];
        
        for cellNum = 1:Ns_cells_per_module
            cellName = sprintf('Cell%02d', cellNum);
            if isfield(module_data.(moduleName), cellName)
                c_data = module_data.(moduleName).(cellName);
                
                if isfield(c_data, 'V_avg') && isfield(c_data, 'R50s_daily') && ...
                   ~isempty(c_data.V_avg) && ~isempty(c_data.R50s_daily)
                    % 유효한 저항값만 수집
                    valid_R = c_data.R50s_daily(~isnan(c_data.R50s_daily)) * 1000; % mΩ
                    all_resistance_values = [all_resistance_values; valid_R(:)];
                end
            end
        end
    end
end

% 최대/최소 저항값 계산
if ~isempty(all_resistance_values)
    global_R_min = min(all_resistance_values);
    global_R_max = max(all_resistance_values);
    % 여유를 위해 5% 여백 추가
    global_R_range = global_R_max - global_R_min;
    global_ylim = [max(0, global_R_min - global_R_range * 0.05), global_R_max + global_R_range * 0.05];
    fprintf('  Global resistance range: [%.4f, %.4f] mΩ\n', global_ylim(1), global_ylim(2));
else
    % 데이터가 없으면 기본값 사용
    global_ylim = [0, 1.5];
    fprintf('  No resistance data found, using default range: [0, 1.5] mΩ\n');
end

%% Step 2: 각 모듈별로 시각화 생성
for modIdx = 1:length(targetModules)
    modNum = targetModules(modIdx);
    moduleName = sprintf('Module%02d', modNum);
    
    fprintf('\n=== Creating Visualizations for %s %s ===\n', targetRack, moduleName);
    
    % Check if module has any data
    hasData = false;
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if isfield(module_data, moduleName) && isfield(module_data.(moduleName), cellName) && ...
           ~isempty(module_data.(moduleName).(cellName).dates)
            hasData = true;
            break;
        end
    end
    
    if ~hasData
        fprintf('  No data available for %s %s\n', targetRack, moduleName);
        continue;  % Skip this module, continue to next
    end
    
    fprintf('  Creating visualizations for %s %s...\n', targetRack, moduleName);
    
    % Collect all cell data
    all_cell_nums = [];
    all_dates = [];
    all_R50s = [];
    all_T = [];
    
    % Get all unique dates across all cells
    all_unique_dates = [];
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if isfield(module_data, moduleName) && isfield(module_data.(moduleName), cellName) && ...
           ~isempty(module_data.(moduleName).(cellName).dates)
            if isempty(all_unique_dates)
                % First datetime array - use it directly
                all_unique_dates = module_data.(moduleName).(cellName).dates(:);
            else
                % Concatenate with existing datetime array
                all_unique_dates = [all_unique_dates; module_data.(moduleName).(cellName).dates(:)];
            end
        end
    end
    
    % Get unique dates and sort if data exists
    if ~isempty(all_unique_dates)
        all_unique_dates = unique(all_unique_dates);
        all_unique_dates = sort(all_unique_dates);
    else
        fprintf('  No date data available for visualization\n');
        continue;  % Skip this module
    end
    
    % Create date index mapping
    date_indices = (1:length(all_unique_dates))';
    
    % Prepare matrices for bar3
    R50s_matrix = NaN(length(all_unique_dates), Ns_cells_per_module);
    T_matrix = NaN(length(all_unique_dates), Ns_cells_per_module);
    
    % Fill matrices and collect scatter data
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if ~isfield(module_data, moduleName) || ~isfield(module_data.(moduleName), cellName) || ...
           isempty(module_data.(moduleName).(cellName).dates)
            continue;
        end
        
        dates = module_data.(moduleName).(cellName).dates;
        R50s = module_data.(moduleName).(cellName).R50s_daily;
        T = module_data.(moduleName).(cellName).T_avg;
        
        % Find date indices using ismember (much faster than find in loop)
        [~, date_indices] = ismember(dates, all_unique_dates);
        valid_date_idx = date_indices > 0;
        
        if any(valid_date_idx)
            date_indices_valid = date_indices(valid_date_idx);
            R50s_valid = R50s(valid_date_idx);
            T_valid = T(valid_date_idx);
            
            % Fill matrices
            R50s_matrix(date_indices_valid, cellNum) = R50s_valid * 1000;  % mΩ
            T_matrix(date_indices_valid, cellNum) = T_valid;
            
            % Collect for scatter
            all_cell_nums = [all_cell_nums; repmat(cellNum, length(date_indices_valid), 1)];
            all_dates = [all_dates; date_indices_valid];
            all_R50s = [all_R50s; R50s_valid * 1000];
            all_T = [all_T; T_valid];
        end
    end
    
    if isempty(all_cell_nums)
        fprintf('  No valid data points for visualization\n');
        continue;  % Skip this module
    end
    
    % Prepare month labels for y-axis
    dateYears = zeros(size(all_unique_dates));
    dateMonths = zeros(size(all_unique_dates));
    for i = 1:length(all_unique_dates)
        dateYears(i) = all_unique_dates(i).Year;
        dateMonths(i) = all_unique_dates(i).Month;
    end
    
    uniqueMonths = unique(dateYears*100 + dateMonths);
    monthLabels = {};
    monthPositions = [];
    
    for i = 1:length(uniqueMonths)
        monthVal = uniqueMonths(i);
        yearVal = floor(monthVal/100);
        monthVal = mod(monthVal, 100);
        
        monthIdx = find(dateYears == yearVal & dateMonths == monthVal, 1);
        if ~isempty(monthIdx)
            monthLabels{end+1} = sprintf('%d-%02d', yearVal, monthVal);
            monthPositions(end+1) = monthIdx;
        end
    end
    
    %% Figure 1: 3D Scatter Plot - COMMENTED OUT
    %{
    figure('Position', [100, 100, 1000, 800]);
    
    scatter3(all_cell_nums, all_dates, all_R50s, 60, all_T, 'filled', 'MarkerEdgeColor', 'k');
    
    colormap(flipud(autumn));
    
    xlabel('Cell Number', 'FontSize', fontSize, 'FontWeight', 'bold');
    ylabel('Date Index', 'FontSize', fontSize, 'FontWeight', 'bold');
    zlabel('R50s (mΩ)', 'FontSize', fontSize, 'FontWeight', 'bold');
    title(sprintf('3D Scatter: %s %s', targetRack, moduleName), 'FontSize', fontSize+2);
    
    % --- [정육면체 형태 및 축 설정 핵심] ---
    pbaspect([1 1 1]); % 시각적 박스 비율을 1:1:1 정육면체로 강제 고정
    grid on;
    set(gca, 'BoxStyle', 'full', 'Box', 'on'); % 박스 형태 강조
    % ---------------------------------------
    
    % X-axis: Cell numbers 1~14
    xticks(1:14);
    xticklabels(1:14);
    xlim([0.5, 14.5]);
    
    % Y-axis: Date labels
    set(gca, 'YTick', monthPositions);
    set(gca, 'YTickLabel', monthLabels);
    ylim([1, length(all_unique_dates)]);
    
    ax3d = gca;
    c = colorbar(ax3d, 'eastoutside');
    c.Label.String = 'Temp (°C)';
    
    view(45, 30); % 보는 각도 최적화
    
    % Save figure
    saveas(gcf, fullfile(rackSaveDir, sprintf('%s_%s_R50s_3D_Scatter.fig', targetRack, moduleName)));
    close(gcf);
    %}
    
    %% Figure 2: 개선된 3D Bar Plot (Cell vs Date vs Resistance + Temperature Color) - COMMENTED OUT
    %{
    figure('Position', [100, 100, 1200, 850], 'Color', 'w');
    
    % 바 그래프 생성 (width=0.8로 설정하여 바 사이 간격 확보)
    h = bar3(R50s_matrix, 0.8);
    
    % 온도 기반 색상 맵핑 및 스타일 설정
    for i = 1:length(h)
        zData = get(h(i), 'ZData');
        xData = get(h(i), 'XData');
        yData = get(h(i), 'YData');
        
        % 현재 바(Cell i)의 데이터가 비어있지 않은지 확인
        if ~isempty(xData)
            cellNum = i; % bar3의 각 핸들은 하나의 열(Cell)을 담당함
            
            % 온도 매트릭스(T_matrix)에서 해당 셀의 온도 추출하여 색상 입히기
            % bar3는 하나의 바당 6개의 면을 생성하므로 CData 크기를 맞춰야 함
            cData = NaN(size(zData));
            for row = 1:6:size(zData, 1)
                dateIdx = floor((row-1)/6) + 1;
                if dateIdx <= size(T_matrix, 1)
                    tempVal = T_matrix(dateIdx, cellNum);
                    cData(row:row+5, :) = tempVal;
                end
            end
            
            set(h(i), 'CData', cData, ...
                'FaceColor', 'flat', ...   % 면 전체에 일정한 온도색 입힘
                'EdgeColor', [0.2 0.2 0.2], ... % 바 테두리를 진한 회색으로
                'FaceAlpha', 0.85);        % 약간 투명하게 하여 가독성 향상
        end
    end
    
    % 컬러바 설정 (온도 표시)
    colormap(flipud(autumn)); % 온도가 높을수록 붉은색
    c = colorbar;
    c.Label.String = 'Average Battery Temperature (°C)';
    c.Label.FontWeight = 'bold';
    if ~all(isnan(T_matrix), 'all')
        caxis([nanmin(T_matrix,[],'all')-1, nanmax(T_matrix,[],'all')+1]);
    end
    
    % 축 및 레이블 설정
    xlabel('Cell Number', 'FontSize', fontSize, 'FontWeight', 'bold');
    ylabel('Time (Date Index)', 'FontSize', fontSize, 'FontWeight', 'bold');
    zlabel('Resistance R50s (mΩ)', 'FontSize', fontSize, 'FontWeight', 'bold');
    title({sprintf('%s %s Resistance Analysis', targetRack, moduleName), ...
           'Z-Axis: Resistance | Color: Temperature'}, 'FontSize', fontSize+2);

    % --- [핵심: 정육면체 모양으로 비율 고정] ---
    % X:Y:Z의 시각적 비율을 1 : 1.5 : 1로 설정 (날짜가 많으므로 Y를 조금 더 길게)
    pbaspect([1 1.5 1]);
    
    % X축: 1~14 셀 번호
    xticks(1:14);
    xticklabels(1:14);
    xlim([0.5, 14.5]);
    
    % Y축: 날짜 레이블 (기존 monthLabels 활용)
    if ~isempty(monthPositions)
        set(gca, 'YTick', monthPositions);
        set(gca, 'YTickLabel', monthLabels);
    end
    ylim([0.5, length(all_unique_dates) + 0.5]);
    
    % Z축: 0부터 시작하여 근거 제시
    z_vals = R50s_matrix(~isnan(R50s_matrix));
    if ~isempty(z_vals)
        zlim([0, max(z_vals)*1.1]);
    end
    
    grid on;
    set(gca, 'BoxStyle', 'full', 'Box', 'on');
    view(-35, 30); % 뒤쪽 바가 더 잘 보이도록 각도 조정
    
    % Save figure
    saveas(gcf, fullfile(rackSaveDir, sprintf('%s_%s_R50s_3D_Bar.fig', targetRack, moduleName)));
    close(gcf);
    %}
    %}
    
    %% Figure 3: Monthly R=dV/dI Plot (separate figure for each month) - COMMENTED OUT
    %{
    % Get all unique dates for this module
    all_dates_for_plot = [];
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if isfield(module_data, moduleName) && isfield(module_data.(moduleName), cellName) && ...
           ~isempty(module_data.(moduleName).(cellName).date_for_materials)
            % Ensure date_for_materials is a datetime array
            dates_cell = module_data.(moduleName).(cellName).date_for_materials(:);
            if ~isdatetime(dates_cell)
                dates_cell = datetime(dates_cell);
            end
            
            if isempty(all_dates_for_plot)
                % First datetime array - use it directly
                all_dates_for_plot = dates_cell;
            else
                % Concatenate with existing datetime array
                all_dates_for_plot = [all_dates_for_plot; dates_cell];
            end
        end
    end
    
    % Get unique dates and sort if data exists
    if ~isempty(all_dates_for_plot)
        % Ensure it's a datetime array before unique/sort
        if ~isdatetime(all_dates_for_plot)
            all_dates_for_plot = datetime(all_dates_for_plot);
        end
        all_dates_for_plot = unique(all_dates_for_plot);
        all_dates_for_plot = sort(all_dates_for_plot);
        
        % unique() and sort() may change the type, so check again
        if ~isdatetime(all_dates_for_plot)
            all_dates_for_plot = datetime(all_dates_for_plot);
        end
    end
    
    if ~isempty(all_dates_for_plot)
        % Final check: Ensure it's a datetime array before using year/month functions
        if ~isdatetime(all_dates_for_plot)
            try
                all_dates_for_plot = datetime(all_dates_for_plot);
            catch
                fprintf('  ERROR: Cannot convert all_dates_for_plot to datetime. Type: %s\n', class(all_dates_for_plot));
                continue;  % Skip this module
            end
        end
        
        % Group dates by month
        try
            % Extract year and month from datetime array using array indexing
            dateYears = zeros(size(all_dates_for_plot));
            dateMonths = zeros(size(all_dates_for_plot));
            for i = 1:length(all_dates_for_plot)
                dateYears(i) = all_dates_for_plot(i).Year;
                dateMonths(i) = all_dates_for_plot(i).Month;
            end
        catch ME
            fprintf('  ERROR in year/month functions: %s\n', ME.message);
            fprintf('  all_dates_for_plot type: %s, size: %s\n', class(all_dates_for_plot), mat2str(size(all_dates_for_plot)));
            continue;  % Skip this module
        end
        
        uniqueMonthKeys = unique(dateYears*100 + dateMonths);
        
        % Process each month separately
        for monthIdx = 1:length(uniqueMonthKeys)
            monthKey = uniqueMonthKeys(monthIdx);
            yearVal = floor(monthKey/100);
            monthVal = mod(monthKey, 100);
            
            % Find dates in this month
            month_mask = (dateYears*100 + dateMonths) == monthKey;
            dates_in_month = all_dates_for_plot(month_mask);
            
            if isempty(dates_in_month)
                continue;
            end
            
            % Calculate subplot grid size for this month
            numDates = length(dates_in_month);
            cols = ceil(sqrt(numDates));
            rows = ceil(numDates / cols);
            
            figure('Position', [100, 100, 1600, 1200]);
            
            for d = 1:numDates
                date = dates_in_month(d);
                subplot(rows, cols, d);
                
                % Collect all dV and dI for this date across all cells
                all_dV_date = [];
                all_dI_date = [];
                
                for cellNum = 1:Ns_cells_per_module
                    cellName = sprintf('Cell%02d', cellNum);
                    if isfield(module_data, moduleName) && isfield(module_data.(moduleName), cellName)
                        % Use ismember instead of find for better performance
                        date_match = ismember(module_data.(moduleName).(cellName).date_for_materials, date);
                        if any(date_match)
                            date_idx = find(date_match, 1);
                            if date_idx <= length(module_data.(moduleName).(cellName).dI_materials)
                                dI_cell = module_data.(moduleName).(cellName).dI_materials{date_idx};
                                dV_cell = module_data.(moduleName).(cellName).dV_materials{date_idx};
                                
                                if ~isempty(dI_cell) && ~isempty(dV_cell)
                                    all_dI_date = [all_dI_date; dI_cell(:)];
                                    all_dV_date = [all_dV_date; dV_cell(:)];
                                end
                            end
                        end
                    end
                end
                
                if ~isempty(all_dI_date) && length(all_dI_date) >= 10
                    % Scatter plot
                    scatter(all_dI_date, all_dV_date, 30, 'filled', 'MarkerFaceAlpha', 0.6);
                    hold on;
                    
                    % Linear fitting using matrix operation (faster than polyfit)
                    X = [all_dI_date, ones(size(all_dI_date))];
                    coeffs = X \ all_dV_date;
                    R_fit = coeffs(1);
                    dI_fit = linspace(min(all_dI_date), max(all_dI_date), 100);
                    dV_fit = coeffs(1) * dI_fit + coeffs(2);
                    
                    % Plot fitted line
                    plot(dI_fit, dV_fit, 'r-', 'LineWidth', 2);
                    
                    % Calculate R-squared
                    dV_predicted = X * coeffs;
                    ss_res = sum((all_dV_date - dV_predicted).^2);
                    ss_tot = sum((all_dV_date - mean(all_dV_date)).^2);
                    r_squared = 1 - (ss_res / ss_tot);
                    
                    % Title with date and R value
                    title(sprintf('%s\nR=%.4f mΩ, R²=%.4f', datestr(date, 'yyyy-mm-dd'), R_fit*1000, r_squared), ...
                        'FontSize', fontSize-2);
                    xlabel('dI (A)', 'FontSize', fontSize-2);
                    ylabel('dV (V)', 'FontSize', fontSize-2);
                    grid on;
                else
                    title(sprintf('%s\nNo data', datestr(date, 'yyyy-mm-dd')), 'FontSize', fontSize-2);
                end
            end
            
            sgtitle(sprintf('%s %s - Monthly R=dV/dI Plots (%d-%02d)', targetRack, moduleName, yearVal, monthVal), ...
                'FontSize', fontSize+2, 'FontWeight', 'bold');
            
            % Save figure with month identifier
            saveas(gcf, fullfile(rackSaveDir, sprintf('%s_%s_Daily_R_dVdI_%d%02d.fig', targetRack, moduleName, yearVal, monthVal)));
            close(gcf);
        end  % End of monthIdx loop
    end  % End of if ~isempty(all_dates_for_plot)
    %}
    
    %% Figure 4: 전압 구간별 저항 분포 박스플롯 (공통 전압 bin마다 별도 figure, 셀별×연도별)
    % Check if V_avg field exists
    has_V_avg = false;
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if isfield(module_data, moduleName) && isfield(module_data.(moduleName), cellName) && ...
           isfield(module_data.(moduleName).(cellName), 'V_avg') && ...
           isfield(module_data.(moduleName).(cellName), 'dates')
            has_V_avg = true;
            break;
        end
    end
    
    if has_V_avg
        % 1. 전압 구간 설정 (3.2V ~ 4.2V까지 0.05V 간격)
        v_bins = 3.2:0.05:4.2;
        num_bins = length(v_bins) - 1;
        
        % 2. 각 전압 bin별로 3개 연도(2021, 2022, 2023) 모두에 데이터가 있는지 확인
        common_bins = [];
        years_list = [2021, 2022, 2023];  % 확인할 연도
        
        for b = 1:num_bins
            v_start = v_bins(b);
            v_end = v_bins(b+1);
            
            % 각 연도별로 데이터 존재 여부 확인
            has_data_by_year = false(1, length(years_list));
            
            for cellNum = 1:Ns_cells_per_module
                cellName = sprintf('Cell%02d', cellNum);
                if isfield(module_data, moduleName) && isfield(module_data.(moduleName), cellName)
                    c_data = module_data.(moduleName).(cellName);
                    
                    if isfield(c_data, 'V_avg') && isfield(c_data, 'dates') && ...
                       ~isempty(c_data.V_avg) && ~isempty(c_data.dates)
                        
                        % 전압 조건에 맞는 인덱스 찾기
                        voltage_idx = (c_data.V_avg >= v_start) & (c_data.V_avg < v_end) & ~isnan(c_data.R50s_daily);
                        
                        if any(voltage_idx)
                            % 해당 인덱스의 날짜에서 연도 추출
                            dates_in_bin = c_data.dates(voltage_idx);
                            years_in_bin = year(dates_in_bin);
                            
                            % 각 연도별로 데이터 존재 확인
                            for y = 1:length(years_list)
                                if any(years_in_bin == years_list(y))
                                    has_data_by_year(y) = true;
                                end
                            end
                        end
                    end
                end
            end
            
            % 3개 연도 모두에 데이터가 있으면 공통 bin으로 추가
            if all(has_data_by_year)
                common_bins = [common_bins; b];
            end
        end
        
        % 3. 공통 전압 bin마다 별도의 figure 생성
        for bin_idx = 1:length(common_bins)
            b = common_bins(bin_idx);
            v_start = v_bins(b);
            v_end = v_bins(b+1);
            bin_label = sprintf('%.2f-%.2fV', v_start, v_end);
            
            figure('Position', [100, 100, 1600, 800], 'Color', 'w');
            
            % 4. 각 셀별, 연도별로 데이터 수집
            all_R_combined = [];
            cell_group = [];  % 셀 번호 그룹 (1-14)
            year_group = []; % 연도 그룹 (2021, 2022, 2023)
            
            for cellNum = 1:Ns_cells_per_module
                cellName = sprintf('Cell%02d', cellNum);
                if isfield(module_data, moduleName) && isfield(module_data.(moduleName), cellName)
                    c_data = module_data.(moduleName).(cellName);
                    
                    if isfield(c_data, 'V_avg') && isfield(c_data, 'dates') && ...
                       ~isempty(c_data.V_avg) && ~isempty(c_data.dates)
                        
                        % 전압 조건에 맞는 인덱스 찾기
                        voltage_idx = (c_data.V_avg >= v_start) & (c_data.V_avg < v_end) & ~isnan(c_data.R50s_daily);
                        
                        if any(voltage_idx)
                            dates_in_bin = c_data.dates(voltage_idx);
                            R_in_bin = c_data.R50s_daily(voltage_idx) * 1000; % mΩ 단위
                            years_in_bin = year(dates_in_bin);
                            
                            % 연도별로 데이터 수집
                            for y = 1:length(years_list)
                                year_val = years_list(y);
                                year_idx = (years_in_bin == year_val);
                                
                                if any(year_idx)
                                    R_year = R_in_bin(year_idx);
                                    all_R_combined = [all_R_combined; R_year(:)];
                                    cell_group = [cell_group; repmat(cellNum, length(R_year), 1)];
                                    year_group = [year_group; repmat(y, length(R_year), 1)];
                                end
                            end
                        end
                    end
                end
            end
            
            if ~isempty(all_R_combined)
                % 연도별 색상 정의: 2021=초록, 2022=파랑, 2023=빨강
                year_colors_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
                year_colors_map(2021) = [0, 0.8, 0];      % 초록
                year_colors_map(2022) = [0, 0.4, 1];      % 파랑
                year_colors_map(2023) = [1, 0, 0];        % 빨강
                
                % boxplot with two grouping variables: cell number and year
                % X축: 셀 번호 (1-14), 각 셀 위치마다 연도별 박스가 나란히
                bp = boxplot(all_R_combined, {cell_group, year_group}, ...
                    'LabelVerbosity', 'minor', 'PlotStyle', 'traditional');
                
                % X축 레이블 설정: 셀 번호만 표시
                ax = gca;
                unique_cells = unique(cell_group);
                num_years = length(years_list);
                
                % 박스플롯 내부 선(중앙값선, 평균선 등) 색상 및 두께 변경
                % boxplot의 구조: 각 셀별로 연도 순서대로 박스가 생성됨
                % {cell_group, year_group} 순서이므로 같은 셀 내에서 연도 순서대로
                
                % 모든 박스 관련 객체 찾기
                h_boxes = findobj(gca, 'Tag', 'Box');
                h_medians = findobj(gca, 'Tag', 'Median');
                h_whiskers = findobj(gca, 'Tag', 'Whisker');
                h_outliers = findobj(gca, 'Tag', 'Outliers');
                
                % 박스 순서: 각 셀별로 연도 순서대로 (2021, 2022, 2023)
                box_idx = 0;
                for c = 1:length(unique_cells)
                    cell_num = unique_cells(c);
                    for y = 1:length(years_list)
                        year_val = years_list(y);
                        box_idx = box_idx + 1;
                        
                        % 해당 연도의 색상 가져오기
                        if isKey(year_colors_map, year_val)
                            box_color = year_colors_map(year_val);
                        else
                            box_color = [0.5, 0.5, 0.5];  % 기본 회색
                        end
                        
                        % 박스 색상 및 두께 변경 (반값: 4.5 -> 2.25)
                        if box_idx <= length(h_boxes)
                            set(h_boxes(box_idx), 'LineWidth', 2.25, 'Color', box_color);
                        end
                        
                        % 중앙값선 색상 및 두께 변경 (반값: 4.5 -> 2.25)
                        if box_idx <= length(h_medians)
                            set(h_medians(box_idx), 'LineWidth', 2.25, 'Color', box_color);
                        end
                        
                        % Whisker 색상 변경
                        if box_idx <= length(h_whiskers)
                            set(h_whiskers(box_idx), 'LineWidth', 2, 'Color', box_color);
                        end
                    end
                end
                
                % 모든 선 객체 찾아서 색상 및 두께 변경 (추가 선들)
                h_all_lines = findobj(gca, 'Type', 'line');
                box_idx = 0;
                for c = 1:length(unique_cells)
                    for y = 1:length(years_list)
                        year_val = years_list(y);
                        box_idx = box_idx + 1;
                        
                        if isKey(year_colors_map, year_val)
                            box_color = year_colors_map(year_val);
                        else
                            box_color = [0.5, 0.5, 0.5];
                        end
                        
                        % 각 박스와 관련된 선들 찾기 (위치 기반)
                        for line_idx = 1:length(h_all_lines)
                            line_obj = h_all_lines(line_idx);
                            % 이미 처리한 객체는 제외
                            if ~ismember(line_obj, [h_boxes; h_medians; h_whiskers; h_outliers])
                                % 박스 근처의 선인지 확인 (X 위치 기반)
                                x_data = get(line_obj, 'XData');
                                if ~isempty(x_data)
                                    expected_x = box_idx;
                                    if any(abs(x_data - expected_x) < 0.3)
                                        set(line_obj, 'Color', box_color);
                                        current_width = get(line_obj, 'LineWidth');
                                        if current_width > 0 && current_width < 2
                                            set(line_obj, 'LineWidth', current_width * 1.5);  % 반값: 3배 -> 1.5배
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                % 각 셀의 중앙 위치 계산
                tick_positions = [];
                tick_labels = {};
                for c = 1:length(unique_cells)
                    cell_num = unique_cells(c);
                    % 각 셀의 중앙 위치: (셀번호-1) * 연도수 + (연도수+1)/2
                    center_pos = (c - 1) * num_years + (num_years + 1) / 2;
                    tick_positions = [tick_positions; center_pos];
                    tick_labels{end+1} = sprintf('Cell%02d', cell_num);
                end
                
                set(gca, 'XTick', tick_positions);
                set(gca, 'XTickLabel', tick_labels);
                ax.XAxis.TickLabelRotation = 45;
                
                grid on;
                xlabel('Cell Number', 'FontSize', fontSize, 'FontWeight', 'bold');
                ylabel('Resistance R50s (mΩ)', 'FontSize', fontSize, 'FontWeight', 'bold');
                ylim(global_ylim);  % 모든 박스플롯에 동일한 Y축 범위 적용
                
                % 셀별 박스플롯 사이에 경계선 추가 (y축 범위에 맞게)
                hold on;
                for c = 1:length(unique_cells)-1
                    % 각 셀 그룹의 끝 위치 계산
                    boundary_pos = c * num_years + 0.5;
                    % 수직선 그리기 (y축 범위 전체)
                    plot([boundary_pos, boundary_pos], global_ylim, 'k-', 'LineWidth', 1.5);
                end
                title(sprintf('%s %s: Resistance Distribution by Cell (Voltage: %s)', targetRack, moduleName, bin_label), ...
                    'FontSize', fontSize+2, 'FontWeight', 'bold');
                
                % 범례 추가 (연도별) - 색상 지정: 2021=초록, 2022=파랑, 2023=빨강
                hold on;
                legend_handles = [];
                legend_labels = {};
                % 연도별 색상 정의
                year_colors = containers.Map('KeyType', 'double', 'ValueType', 'any');
                year_colors(2021) = [0, 0.8, 0];      % 초록
                year_colors(2022) = [0, 0.4, 1];      % 파랑
                year_colors(2023) = [1, 0, 0];        % 빨강
                
                colors = zeros(length(years_list), 3);
                for y = 1:length(years_list)
                    year_val = years_list(y);
                    if isKey(year_colors, year_val)
                        colors(y,:) = year_colors(year_val);
                    else
                        colors(y,:) = [0.5, 0.5, 0.5];  % 기본 회색
                    end
                    h = plot(NaN, NaN, 's', 'Color', colors(y,:), 'MarkerFaceColor', colors(y,:), 'MarkerSize', 10);
                    legend_handles = [legend_handles; h];
                    legend_labels{end+1} = sprintf('%d', year_val);
                end
                legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', fontSize-2);
                
                % 전체 평균선 추가
                yline(nanmean(all_R_combined), 'r--', sprintf('Total Average: %.4f mΩ', nanmean(all_R_combined)), ...
                    'LineWidth', 1.5, 'FontSize', fontSize);
                
                % 통계 분석 수행 (셀별로만)
                fprintf('\n  === Statistical Analysis for %s (Voltage: %s) - Cell-wise ===\n', moduleName, bin_label);
                
                % 셀별 연도 간 차이 검정
                cell_stats = struct();
                for c = 1:length(unique_cells)
                    cell_num = unique_cells(c);
                    cell_idx = (cell_group == cell_num);
                    cell_R = all_R_combined(cell_idx);
                    cell_years = year_group(cell_idx);
                    
                    if length(unique(cell_years)) >= 2 && length(cell_R) >= 6
                        try
                            [p_cell, ~, stats_cell] = kruskalwallis(cell_R, cell_years, 'off');
                            cell_stats.(sprintf('Cell%02d', cell_num)).p_value = p_cell;
                            cell_stats.(sprintf('Cell%02d', cell_num)).mean_R = nanmean(cell_R);
                            cell_stats.(sprintf('Cell%02d', cell_num)).std_R = nanstd(cell_R);
                            
                            % 연도별 평균값
                            for y = 1:length(years_list)
                                year_R = cell_R(cell_years == y);
                                if ~isempty(year_R)
                                    cell_stats.(sprintf('Cell%02d', cell_num)).(sprintf('Y%d_mean', years_list(y))) = nanmean(year_R);
                                    cell_stats.(sprintf('Cell%02d', cell_num)).(sprintf('Y%d_std', years_list(y))) = nanstd(year_R);
                                    cell_stats.(sprintf('Cell%02d', cell_num)).(sprintf('Y%d_N', years_list(y))) = length(year_R);
                                end
                            end
                        catch
                            cell_stats.(sprintf('Cell%02d', cell_num)).p_value = NaN;
                        end
                    end
                end
                
                % 모듈별 통계표 형식으로 출력
                fprintf('\n    === Module Statistics Table ===\n');
                fprintf('    Module: %s | Voltage Bin: %s\n', moduleName, bin_label);
                fprintf('    %s\n', repmat('=', 1, 120));
                fprintf('    %-8s | %-12s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s\n', ...
                    'Cell', 'p-value', 'Mean', 'Std', 'Y2021_M', 'Y2021_S', 'Y2022_M', 'Y2022_S', 'Y2023_M', 'Y2023_S');
                fprintf('    %s\n', repmat('-', 1, 120));
                
                significant_cells = [];
                for c = 1:length(unique_cells)
                    cell_num = unique_cells(c);
                    cell_name = sprintf('Cell%02d', cell_num);
                    if isfield(cell_stats, cell_name) && isfield(cell_stats.(cell_name), 'p_value')
                        p_cell = cell_stats.(cell_name).p_value;
                        mean_R = cell_stats.(cell_name).mean_R;
                        std_R = cell_stats.(cell_name).std_R;
                        
                        if ~isnan(p_cell)
                            sig_marker = '';
                            if p_cell < 0.001
                                sig_marker = '***';
                                significant_cells = [significant_cells; cell_num];
                            elseif p_cell < 0.01
                                sig_marker = '**';
                            elseif p_cell < 0.05
                                sig_marker = '*';
                            end
                            
                            % 연도별 데이터 추출
                            y2021_mean = NaN; y2021_std = NaN;
                            y2022_mean = NaN; y2022_std = NaN;
                            y2023_mean = NaN; y2023_std = NaN;
                            
                            if isfield(cell_stats.(cell_name), 'Y2021_mean')
                                y2021_mean = cell_stats.(cell_name).Y2021_mean;
                                y2021_std = cell_stats.(cell_name).Y2021_std;
                            end
                            if isfield(cell_stats.(cell_name), 'Y2022_mean')
                                y2022_mean = cell_stats.(cell_name).Y2022_mean;
                                y2022_std = cell_stats.(cell_name).Y2022_std;
                            end
                            if isfield(cell_stats.(cell_name), 'Y2023_mean')
                                y2023_mean = cell_stats.(cell_name).Y2023_mean;
                                y2023_std = cell_stats.(cell_name).Y2023_std;
                            end
                            
                            fprintf('    %-8s | %-10.4e%s | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f\n', ...
                                cell_name, p_cell, sig_marker, mean_R, std_R, ...
                                y2021_mean, y2021_std, y2022_mean, y2022_std, y2023_mean, y2023_std);
                        else
                            fprintf('    %-8s | %-12s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s\n', ...
                                cell_name, 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A');
                        end
                    else
                        fprintf('    %-8s | %-12s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s\n', ...
                            cell_name, 'No Data', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A');
                    end
                end
                fprintf('    %s\n', repmat('=', 1, 120));
                
                % 유의한 차이를 보이는 셀 요약
                if ~isempty(significant_cells)
                    fprintf('    Significant cells (p < 0.05): %s\n', mat2str(significant_cells));
                else
                    fprintf('    No cells with significant year-to-year difference (p >= 0.05)\n');
                end
                
            else
                text(0.5, 0.5, 'No Data in This Voltage Range', 'HorizontalAlignment', 'center');
            end
            
            % 파일명에 전압 구간 포함 (예: 3.6-3.7V -> 3_6-3_7V)
            voltage_str = strrep(bin_label, '.', '_');
            saveas(gcf, fullfile(rackSaveDir, sprintf('%s_%s_Voltage_%s_Boxplot.fig', targetRack, moduleName, voltage_str)));
            close(gcf);
            
            % 통계 결과를 구조체로 저장 (셀별 통계만)
            stats_result = struct();
            stats_result.voltage_bin = bin_label;
            stats_result.module_name = moduleName;
            stats_result.cell_stats = cell_stats;
            stats_result.significant_cells = significant_cells;
            
            % 통계 결과를 MAT 파일로 저장
            stats_filename = fullfile(rackSaveDir, sprintf('%s_%s_Voltage_%s_Statistics.mat', targetRack, moduleName, voltage_str));
            save(stats_filename, 'stats_result', '-v7.3');
            
            fprintf('  %s: Voltage binning boxplot and statistics saved for %s\n', moduleName, bin_label);
        end
        
        if isempty(common_bins)
            fprintf('  %s: No common voltage bins found across all years\n', moduleName);
        end
    else
        fprintf('  %s: V_avg field not found, skipping voltage binning boxplot\n', moduleName);
    end
    
end  % End of modIdx loop

%% Module-wise Cell Resistance Distribution Boxplot (Latest Date)
fprintf('\n=== Creating Module-wise Boxplot (Latest Date) ===\n');

% Collect latest resistance data for all modules
latest_R_matrix = NaN(Ns_cells_per_module, numModules);  % [14 cells x 17 modules]

for modNum = 1:numModules
    moduleName_temp = sprintf('Module%02d', modNum);
    
    % Check if this module has data
    if ~isfield(module_data, moduleName_temp)
        continue;
    end
    
    % Get latest date for this module
    all_dates_module = [];
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if isfield(module_data.(moduleName_temp), cellName) && ...
           ~isempty(module_data.(moduleName_temp).(cellName).dates)
            all_dates_module = [all_dates_module; module_data.(moduleName_temp).(cellName).dates];
        end
    end
    
    if isempty(all_dates_module)
        continue;
    end
    
    % Get unique dates and find latest
    all_dates_module = unique(all_dates_module);
    all_dates_module = sort(all_dates_module);
    latest_date = all_dates_module(end);
    
    % Extract resistance values for latest date
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if isfield(module_data.(moduleName_temp), cellName) && ...
           ~isempty(module_data.(moduleName_temp).(cellName).dates)
            dates = module_data.(moduleName_temp).(cellName).dates;
            R50s = module_data.(moduleName_temp).(cellName).R50s_daily;
            
            % Find latest date index
            date_match = (dates == latest_date);
            if any(date_match)
                date_idx = find(date_match, 1);
                latest_R_matrix(cellNum, modNum) = R50s(date_idx) * 1000;  % Convert to mΩ
            end
        end
    end
end

% Create boxplot only if we have data
valid_modules = ~all(isnan(latest_R_matrix), 1);
if any(valid_modules)
    % Filter out modules with no data
    latest_R_matrix_filtered = latest_R_matrix(:, valid_modules);
    module_labels = cell(1, sum(valid_modules));
    module_idx = 1;
    for modNum = 1:numModules
        if valid_modules(modNum)
            module_labels{module_idx} = sprintf('Mod%02d', modNum);
            module_idx = module_idx + 1;
        end
    end
    
    figure('Position', [100, 100, 1200, 600]);
    boxplot(latest_R_matrix_filtered, 'Labels', module_labels);
    
    grid on;
    xlabel('Module Number (Rack01)', 'FontSize', fontSize, 'FontWeight', 'bold');
    ylabel('Resistance R50s (mΩ)', 'FontSize', fontSize, 'FontWeight', 'bold');
    ylim(global_ylim);  % 모든 박스플롯에 동일한 Y축 범위 적용
    title(sprintf('Rack01: Cell Resistance Distribution per Module (Latest Date: %s)', ...
        datestr(latest_date, 'yyyy-mm-dd')), 'FontSize', fontSize+2, 'FontWeight', 'bold');
    
    % Add total average line
    hold on;
    total_avg = nanmean(latest_R_matrix_filtered(:));
    if ~isnan(total_avg)
        yline(total_avg, 'r--', sprintf('Total Average: %.4f mΩ', total_avg), ...
            'LineWidth', 1.5, 'FontSize', fontSize);
    end
    
    % 박스플롯 내부 선(중앙값선, 평균선 등) 두께를 반값으로 (4.5 -> 2.25)
    h = findobj(gca, 'Tag', 'Box');
    for i = 1:length(h)
        set(h(i), 'LineWidth', 2.25);  % 반값: 4.5 -> 2.25
    end
    
    % 중앙값선 (Median line)
    h_median = findobj(gca, 'Tag', 'Median');
    for i = 1:length(h_median)
        set(h_median(i), 'LineWidth', 2.25);
    end
    
    % 평균선 (Mean line, if shown)
    h_mean = findobj(gca, 'Tag', 'Mean');
    for i = 1:length(h_mean)
        set(h_mean(i), 'LineWidth', 2.25);
    end
    
    % 박스 외곽선
    h_outline = findobj(gca, 'Type', 'line');
    for i = 1:length(h_outline)
        current_width = get(h_outline(i), 'LineWidth');
        if current_width > 0 && current_width < 2  % 박스플롯 내부 선들만
            set(h_outline(i), 'LineWidth', current_width * 1.5);  % 반값: 3배 -> 1.5배
        end
    end
    
    % Save figure
    saveas(gcf, fullfile(rackSaveDir, sprintf('%s_Module_Resistance_Boxplot_Latest.fig', targetRack)));
    close(gcf);
    
    fprintf('  Boxplot saved: %d modules with data\n', sum(valid_modules));
else
    fprintf('  No data available for boxplot\n');
end

fprintf('\n=== Visualization Complete ===\n');
fprintf('All figures saved to: %s\n', rackSaveDir);

