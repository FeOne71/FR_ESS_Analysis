%% Field Module R50s Analysis - Module-Level Resistance Calculation
% This script calculates and visualizes module-level R50s resistance
% from RARD synchronized data (40-50 second intervals)
% 
% Process:
% 1. Load module data from RARDsync folder
% 2. Calculate module voltage from cell voltages (14 cells in series)
% 3. Calculate R50s using dV/dI linear fitting (no moving average)
% 4. Generate module-specific visualizations
% 5. Save figures to Module Resistance folder

clear; clc; close all;

%% Configuration
fontSize = 12;
saveDir = 'Module Resistance';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Data directory
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\RARDsync';

% Topology
Np = 2;  % 2 parallel
Ns_cells_per_module = 14;  % 14 cells in series per module

% Years to process
years = {'2021', '2022', '2023'};

% Target Rack and Module for visualization (can be modified)
targetRack = 'Rack01';
targetModule = 1:17;  % Module number (1-17)

% Voltage binning configuration (for SOC-based resistance analysis)
enable_voltage_binning = true;  % Set to false to disable voltage filtering
V_min = 3.0;  % Minimum voltage (V)
V_max = 4.2;  % Maximum voltage (V)
V_bin_size = 0.1;  % Voltage bin size (V)
V_bins = V_min:V_bin_size:V_max;  % Voltage bins [3.0, 3.1, ..., 4.2]

rackName = targetRack;  % For data processing
numModules = 17;  % Module01 to Module17 (for data loading)

%% Initialize data storage for all target modules
module_data = struct();
% Initialize structure for all modules in targetModule array
if isscalar(targetModule)
    targetModules = targetModule;
else
    targetModules = targetModule;
end
for modIdx = 1:length(targetModules)
    modNum = targetModules(modIdx);
    moduleName = sprintf('Module%02d', modNum);
    module_data.(moduleName) = struct();
end
% Note: Pre-allocation will be done after file list is collected

%% Load and process data
fprintf('=== Loading Module Data from RARDsync ===\n');
debug_stats = struct();
debug_stats.files_loaded = 0;
debug_stats.files_failed = 0;
debug_stats.total_days_processed = 0;
debug_stats.start_time = tic;

% Performance optimization: Collect all file paths first
fprintf('Collecting all mat file paths...\n');
allMatFiles = [];
for year_idx = 1:length(years)
    year = years{year_idx};
    yearPath = fullfile(dataDir, year);
    if ~exist(yearPath, 'dir')
        continue;
    end
    
    % Use recursive dir to find all Raw_*.mat files
    matFiles = dir(fullfile(yearPath, '**', 'Raw_*.mat'));
    for f = 1:length(matFiles)
        allMatFiles = [allMatFiles; struct('folder', matFiles(f).folder, 'name', matFiles(f).name)];
    end
end

fprintf('Found %d mat files to process\n', length(allMatFiles));

% Pre-allocate cell arrays based on file count (critical for performance)
numFiles = length(allMatFiles);
% Pre-allocate for all target modules
if isscalar(targetModule)
    targetModules = targetModule;
else
    targetModules = targetModule;
end
for modIdx = 1:length(targetModules)
    modNum = targetModules(modIdx);
    moduleName = sprintf('Module%02d', modNum);
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        module_data.(moduleName).(cellName) = struct();
        % Pre-allocate cell arrays with known size
        module_data.(moduleName).(cellName).dates_cell = cell(numFiles, 1);
        module_data.(moduleName).(cellName).R50s_daily_cell = cell(numFiles, 1);
        module_data.(moduleName).(cellName).T_avg_cell = cell(numFiles, 1);
        module_data.(moduleName).(cellName).n_points_cell = cell(numFiles, 1);
        module_data.(moduleName).(cellName).r_squared_cell = cell(numFiles, 1);
        module_data.(moduleName).(cellName).dV_materials = cell(numFiles, 1);
        module_data.(moduleName).(cellName).dI_materials = cell(numFiles, 1);
        module_data.(moduleName).(cellName).date_for_materials_cell = cell(numFiles, 1);
    end
end

% Process files using parallel processing
% Pre-allocate results cell array for parfor compatibility
results_cell = cell(numFiles, 1);
files_loaded_count = zeros(numFiles, 1);
files_failed_count = zeros(numFiles, 1);
days_processed_count = zeros(numFiles, 1);

parfor f = 1:numFiles
    matFilePath = fullfile(allMatFiles(f).folder, allMatFiles(f).name);
    dayKey = allMatFiles(f).name(1:end-4);  % Remove .mat extension
    
    % Extract date from filename (Raw_YYYYMMDD)
    if length(dayKey) >= 12 && strcmp(dayKey(1:4), 'Raw_')
        dateStr = dayKey(5:12);
        year_num = str2double(dateStr(1:4));
        month = str2double(dateStr(5:6));
        day = str2double(dateStr(7:8));
        date = datetime(year_num, month, day);
    else
        continue;
    end
    
    if mod(f, 10) == 1 || f == numFiles
        fprintf('  Processing file %d/%d: %s\n', f, numFiles, allMatFiles(f).name);
    end
    
    % Initialize result structure for this file (for all modules)
    file_result = struct();
    file_result.date = date;
    file_result.module_data = struct();  % Store data for each module
    file_result.files_loaded = 0;
    file_result.files_failed = 0;
    file_result.days_processed = 0;
            
    try
        % Use matfile for better performance (loads only what's needed)
        m = matfile(matFilePath);
        if ~isprop(m, 'Raw')
            continue;
        end
        
        Raw_all = m.Raw;  % Load Raw variable
        file_result.files_loaded = 1;
        
        if ~isfield(Raw_all, rackName)
            continue;
        end
        
        rackData = Raw_all.(rackName);
        
        % Process all target modules
        if isscalar(targetModule)
            targetModules = targetModule;
        else
            targetModules = targetModule;
        end
        
        % Process each target module
        for modIdx = 1:length(targetModules)
            modNum = targetModules(modIdx);
            moduleName_local = sprintf('Module%02d', modNum);
            
            if ~isfield(rackData, moduleName_local)
                continue;
            end
            
            moduleData = rackData.(moduleName_local);
            
            % Check for required fields
            if ~isfield(moduleData, 'Time') || isempty(moduleData.Time)
                continue;
            end
            
            % Extract current from RBMS synchronized data
            if ~isfield(moduleData, 'RBMS_DCCurrent_A')
                continue;
            end
            
            I_rack = moduleData.RBMS_DCCurrent_A;
            
            % Extract temperature
            if ~isfield(moduleData, 'RBMS_AverageMT_degC')
                continue;
            end
            
            T_batt = moduleData.RBMS_AverageMT_degC;
            
            % Convert to cell current
            I_cell = I_rack / Np;
            
            % Process each cell individually
            for cellNum = 1:Ns_cells_per_module
            cellName = sprintf('Cell%02d', cellNum);
            cellField = sprintf('M%d_Cell%d', modNum, cellNum);
            
            if ~isfield(moduleData, cellField)
                continue;
            end
            
            % Extract individual cell voltage
            V_cell = moduleData.(cellField);
            
            % Remove NaN values
            valid_idx = ~isnan(V_cell) & ~isnan(I_cell) & ~isnan(T_batt);
            V_cell = V_cell(valid_idx);
            I_cell_valid = I_cell(valid_idx);
            T_batt_valid = T_batt(valid_idx);
            
            % Apply voltage binning filter if enabled
            if enable_voltage_binning
                voltage_mask = (V_cell >= V_min) & (V_cell <= V_max);
                if any(voltage_mask)
                    V_cell = V_cell(voltage_mask);
                    I_cell_valid = I_cell_valid(voltage_mask);
                    T_batt_valid = T_batt_valid(voltage_mask);
                else
                    continue;  % No data in voltage range
                end
            end
            
            if length(V_cell) < 10
                continue;
            end
            
            % Calculate dV and dI (NO moving average)
            dI = diff(I_cell_valid);
            dV = diff(V_cell);
            T_inst = T_batt_valid(1:end-1);
            
            % Filter valid data points
            valid_idx = ~isnan(dI) & ~isnan(dV) & dI ~= 0 & dV ~= 0;
            dI_valid = dI(valid_idx);
            dV_valid = dV(valid_idx);
            T_inst_valid = T_inst(valid_idx);
            
            % Calculate R50s using matrix operation (faster than polyfit): dV = a * dI
            if length(dI_valid) >= 10
                % Use matrix operation instead of polyfit for better performance
                % dV = R * dI, so R = dI \ dV (left division)
                X = [dI_valid, ones(size(dI_valid))];  % Include intercept for robustness
                coeffs = X \ dV_valid;
                R50s_daily = coeffs(1);  % Slope
                
                % Calculate R-squared
                dV_predicted = X * coeffs;
                ss_res = sum((dV_valid - dV_predicted).^2);
                ss_tot = sum((dV_valid - mean(dV_valid)).^2);
                r_squared = 1 - (ss_res / ss_tot);
                
                % Store results in file_result structure (module-specific)
                if ~isfield(file_result.module_data, moduleName_local)
                    file_result.module_data.(moduleName_local) = struct();
                    file_result.module_data.(moduleName_local).cell_data = cell(Ns_cells_per_module, 1);
                end
                
                cell_result = struct();
                cell_result.cellNum = cellNum;
                cell_result.date = date;
                cell_result.R50s_daily = R50s_daily;
                cell_result.T_avg = mean(T_inst_valid);
                cell_result.n_points = length(dI_valid);
                cell_result.r_squared = r_squared;
                cell_result.dV_valid = dV_valid;
                cell_result.dI_valid = dI_valid;
                
                file_result.module_data.(moduleName_local).cell_data{cellNum} = cell_result;
                file_result.days_processed = file_result.days_processed + 1;
            end
            end  % End of cellNum loop
        end  % End of modIdx loop
        
    catch ME
        fprintf('      ERROR processing %s: %s\n', allMatFiles(f).name, ME.message);
        file_result.files_failed = 1;
    end
    
    % Store result for this file
    results_cell{f} = file_result;
    files_loaded_count(f) = file_result.files_loaded;
    files_failed_count(f) = file_result.files_failed;
    days_processed_count(f) = file_result.days_processed;
end

% Aggregate results from parfor into module_data structure
fprintf('\n=== Aggregating Parallel Processing Results ===\n');
if isscalar(targetModule)
    targetModules = targetModule;
else
    targetModules = targetModule;
end

for f = 1:numFiles
    file_result = results_cell{f};
    if isempty(file_result) || ~isfield(file_result, 'module_data')
        continue;
    end
    
    % Process each module in the result
    module_names = fieldnames(file_result.module_data);
    for modIdx = 1:length(module_names)
        moduleName_result = module_names{modIdx};
        
        if ~isfield(module_data, moduleName_result)
            continue;
        end
        
        module_result = file_result.module_data.(moduleName_result);
        if ~isfield(module_result, 'cell_data')
            continue;
        end
        
        for cellNum = 1:Ns_cells_per_module
            if isempty(module_result.cell_data{cellNum})
                continue;
            end
            
            cell_result = module_result.cell_data{cellNum};
            cellName = sprintf('Cell%02d', cellNum);
            
            % Store results using direct index (pre-allocated cell arrays)
            module_data.(moduleName_result).(cellName).dates_cell{f} = cell_result.date;
            module_data.(moduleName_result).(cellName).R50s_daily_cell{f} = cell_result.R50s_daily;
            module_data.(moduleName_result).(cellName).T_avg_cell{f} = cell_result.T_avg;
            module_data.(moduleName_result).(cellName).n_points_cell{f} = cell_result.n_points;
            module_data.(moduleName_result).(cellName).r_squared_cell{f} = cell_result.r_squared;
            module_data.(moduleName_result).(cellName).dV_materials{f} = cell_result.dV_valid;
            module_data.(moduleName_result).(cellName).dI_materials{f} = cell_result.dI_valid;
            module_data.(moduleName_result).(cellName).date_for_materials_cell{f} = cell_result.date;
        end
    end
end

% Update debug stats
debug_stats.files_loaded = sum(files_loaded_count);
debug_stats.files_failed = sum(files_failed_count);
debug_stats.total_days_processed = sum(days_processed_count);

fprintf('\n=== DATA LOADING SUMMARY ===\n');
fprintf('Files loaded: %d, Failed: %d, Days processed: %d, Loading time: %.1f sec\n', ...
    debug_stats.files_loaded, debug_stats.files_failed, debug_stats.total_days_processed, toc(debug_stats.start_time));

%% Convert cell arrays to regular arrays and sort by date
fprintf('\n=== Converting Cell Arrays and Sorting Data by Date ===\n');
if isscalar(targetModule)
    targetModules = targetModule;
else
    targetModules = targetModule;
end

for modIdx = 1:length(targetModules)
    modNum = targetModules(modIdx);
    moduleName = sprintf('Module%02d', modNum);
    
    if ~isfield(module_data, moduleName)
        continue;
    end
    
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        
        % Remove empty cells and convert to regular arrays
        keepIdx = ~cellfun(@isempty, module_data.(moduleName).(cellName).dates_cell);
        
        if any(keepIdx)
        % Convert dates cell array to datetime array (only non-empty cells)
        dates_temp = module_data.(moduleName).(cellName).dates_cell(keepIdx);
        module_data.(moduleName).(cellName).dates = [dates_temp{:}]';
        
        % Convert numeric cell arrays (only non-empty cells)
        module_data.(moduleName).(cellName).R50s_daily = [module_data.(moduleName).(cellName).R50s_daily_cell{keepIdx}]';
        module_data.(moduleName).(cellName).T_avg = [module_data.(moduleName).(cellName).T_avg_cell{keepIdx}]';
        module_data.(moduleName).(cellName).n_points = [module_data.(moduleName).(cellName).n_points_cell{keepIdx}]';
        module_data.(moduleName).(cellName).r_squared = [module_data.(moduleName).(cellName).r_squared_cell{keepIdx}]';
        
        % Sort by date
        [module_data.(moduleName).(cellName).dates, sortIdx] = sort(module_data.(moduleName).(cellName).dates);
        module_data.(moduleName).(cellName).R50s_daily = module_data.(moduleName).(cellName).R50s_daily(sortIdx);
        module_data.(moduleName).(cellName).T_avg = module_data.(moduleName).(cellName).T_avg(sortIdx);
        module_data.(moduleName).(cellName).n_points = module_data.(moduleName).(cellName).n_points(sortIdx);
        module_data.(moduleName).(cellName).r_squared = module_data.(moduleName).(cellName).r_squared(sortIdx);
        
        % Sort dV/dI materials accordingly
        dV_temp = module_data.(moduleName).(cellName).dV_materials(keepIdx);
        dI_temp = module_data.(moduleName).(cellName).dI_materials(keepIdx);
        dates_materials_temp = module_data.(moduleName).(cellName).date_for_materials_cell(keepIdx);
        
        module_data.(moduleName).(cellName).dV_materials = dV_temp(sortIdx);
        module_data.(moduleName).(cellName).dI_materials = dI_temp(sortIdx);
        % Convert cell array to datetime array
        dates_materials_sorted = dates_materials_temp(sortIdx);
        if ~isempty(dates_materials_sorted)
            % Convert cell array to datetime array safely
            try
                % Check if all elements are datetime
                all_datetime = true;
                for idx = 1:length(dates_materials_sorted)
                    if ~isdatetime(dates_materials_sorted{idx})
                        all_datetime = false;
                        break;
                    end
                end
                
                if all_datetime
                    module_data.(moduleName).(cellName).date_for_materials = [dates_materials_sorted{:}]';
                else
                    % Convert each element to datetime if needed
                    datetime_array = [];
                    for idx = 1:length(dates_materials_sorted)
                        if isdatetime(dates_materials_sorted{idx})
                            datetime_array = [datetime_array; dates_materials_sorted{idx}];
                        else
                            datetime_array = [datetime_array; datetime(dates_materials_sorted{idx})];
                        end
                    end
                    module_data.(moduleName).(cellName).date_for_materials = datetime_array;
                end
            catch ME
                fprintf('  WARNING: Error converting date_for_materials for %s: %s\n', cellName, ME.message);
                module_data.(moduleName).(cellName).date_for_materials = [];
            end
        else
            module_data.(moduleName).(cellName).date_for_materials = [];
        end
        
        % Clear cell arrays to save memory
        module_data.(moduleName).(cellName).dates_cell = {};
        module_data.(moduleName).(cellName).R50s_daily_cell = {};
        module_data.(moduleName).(cellName).T_avg_cell = {};
        module_data.(moduleName).(cellName).n_points_cell = {};
        module_data.(moduleName).(cellName).r_squared_cell = {};
        module_data.(moduleName).(cellName).date_for_materials_cell = {};
        end
    end  % End of cellNum loop
    
    % Print summary for this module
    total_points = 0;
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if isfield(module_data.(moduleName), cellName) && ...
           ~isempty(module_data.(moduleName).(cellName).dates)
            total_points = total_points + length(module_data.(moduleName).(cellName).dates);
        end
    end
    if total_points > 0
        fprintf('  %s: %d total data points across all cells\n', moduleName, total_points);
    end
end  % End of modIdx loop

%% ========================================================================
%% STEP 1: FAST CALCULATION - Save results only (NO visualization)
%% ========================================================================
% Visualization is separated to improve performance
% Results are saved and can be visualized later using saved data

fprintf('\n=== Saving Calculation Results (No Visualization) ===\n');
save(fullfile(saveDir, sprintf('%s_Resistance_Summary.mat', targetRack)), 'module_data', '-v7.3');
fprintf('Calculation results saved to: %s_Resistance_Summary.mat\n', targetRack);

%% ========================================================================
%% STEP 2: VISUALIZATION (Optional - Comment out for faster processing)
%% ========================================================================
% Uncomment the section below to generate visualizations
% This can be run separately after calculation is complete

% enable_visualization = false;  % Set to true to enable visualization
enable_visualization = true;  % Set to false to skip visualization for speed

if enable_visualization
    fprintf('\n=== Creating Visualizations ===\n');
    
    if isscalar(targetModule)
        targetModules = targetModule;
    else
        targetModules = targetModule;
    end
    
    % Create Rack-specific folder
    rackSaveDir = fullfile(saveDir, targetRack);
    if ~exist(rackSaveDir, 'dir')
        mkdir(rackSaveDir);
    end
    
    % Process each target module for visualization
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
    
    %% Figure 1: 3D Scatter Plot
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
    
    %% Figure 2: 3D Bar Plot
    figure('Position', [100, 100, 1000, 800]);
    
    h = bar3(R50s_matrix, 0.8); % 바 두께 조절 (0.8)
    
    % 온도 기반 색상 입히기 로직
    for i = 1:length(h)
        xData = get(h(i), 'XData');
        yData = get(h(i), 'YData');
        zData = get(h(i), 'ZData');
        if ~isempty(xData) && ~isempty(yData)
            cellNum = round(nanmean(xData(:)));
            if cellNum >= 1 && cellNum <= 14
                tempData = NaN(size(zData));
                for row = 1:size(zData, 1)
                    dateIdx = round(nanmean(yData(row, :)));
                    if dateIdx >= 1 && dateIdx <= size(T_matrix, 1)
                        tempData(row, :) = T_matrix(dateIdx, cellNum);
                    end
                end
                set(h(i), 'CData', tempData, 'FaceColor', 'interp');
            end
        end
    end
    
    colormap(flipud(autumn));
    
    xlabel('Cell Number', 'FontSize', fontSize, 'FontWeight', 'bold');
    ylabel('Date Index', 'FontSize', fontSize, 'FontWeight', 'bold');
    zlabel('R50s (mΩ)', 'FontSize', fontSize, 'FontWeight', 'bold');
    title(sprintf('3D Bar: %s %s', targetRack, moduleName), 'FontSize', fontSize+2);
    
    % --- [정육면체 형태 및 축 설정 핵심] ---
    pbaspect([1 1 1]); % 시각적 박스 비율을 1:1:1로 고정
    grid on;
    set(gca, 'BoxStyle', 'full', 'Box', 'on');
    % ---------------------------------------
    
    xticks(1:14);
    xticklabels(1:14);
    xlim([0.5, 14.5]);
    
    set(gca, 'YTick', monthPositions);
    set(gca, 'YTickLabel', monthLabels);
    ylim([0.5, length(all_unique_dates)+0.5]);
    
    ax3d = gca;
    c = colorbar(ax3d, 'eastoutside');
    c.Label.String = 'Temp (°C)';
    
    view(45, 30);
    
    % Save figure
    saveas(gcf, fullfile(rackSaveDir, sprintf('%s_%s_R50s_3D_Bar.fig', targetRack, moduleName)));
    close(gcf);
    
    %% Figure 3: Monthly R=dV/dI Plot (separate figure for each month)
    % Get all unique dates for this module
    all_dates_for_plot = [];
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if ~isempty(module_data.(moduleName).(cellName).date_for_materials)
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
    end  % End of modIdx loop (visualization)
end  % End of enable_visualization

%% Summary Statistics
fprintf('\n=== Summary Statistics ===\n');
if isscalar(targetModule)
    targetModules = targetModule;
else
    targetModules = targetModule;
end

for modIdx = 1:length(targetModules)
    modNum = targetModules(modIdx);
    moduleName = sprintf('Module%02d', modNum);
    
    if ~isfield(module_data, moduleName)
        continue;
    end
    
    all_R50s_module = [];
    total_points = 0;
    
    for cellNum = 1:Ns_cells_per_module
        cellName = sprintf('Cell%02d', cellNum);
        if isfield(module_data.(moduleName), cellName) && ...
           ~isempty(module_data.(moduleName).(cellName).R50s_daily)
            R50s_data = module_data.(moduleName).(cellName).R50s_daily;
            all_R50s_module = [all_R50s_module; R50s_data];
            total_points = total_points + length(R50s_data);
        end
    end
    
    if ~isempty(all_R50s_module)
        fprintf('%s: %d total points across all cells, Mean: %.4f mΩ, Std: %.4f mΩ, Range: %.4f-%.4f mΩ\n', ...
            moduleName, total_points, mean(all_R50s_module)*1000, std(all_R50s_module)*1000, ...
            min(all_R50s_module)*1000, max(all_R50s_module)*1000);
    end
end  % End of module statistics loop

%% Module-wise Cell Resistance Distribution Boxplot (Latest Date)
fprintf('\n=== Creating Module-wise Boxplot (Latest Date) ===\n');

% Collect latest resistance data for all modules
% Note: Currently only targetModule is processed, so boxplot will show only one module
% To show all 17 modules, need to process all modules first
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
    title(sprintf('Rack01: Cell Resistance Distribution per Module (Latest Date: %s)', ...
        datestr(latest_date, 'yyyy-mm-dd')), 'FontSize', fontSize+2, 'FontWeight', 'bold');
    
    % Add total average line
    hold on;
    total_avg = nanmean(latest_R_matrix_filtered(:));
    if ~isnan(total_avg)
        yline(total_avg, 'r--', sprintf('Total Average: %.4f mΩ', total_avg), ...
            'LineWidth', 2, 'FontSize', fontSize);
    end
    
    % Save figure
    saveas(gcf, fullfile(rackSaveDir, sprintf('%s_Module_Resistance_Boxplot_Latest.fig', targetRack)));
    close(gcf);
    
    fprintf('  Boxplot saved: %d modules with data\n', sum(valid_modules));
else
    fprintf('  No data available for boxplot\n');
end

fprintf('\nModule R50s analysis complete! Results saved to: %s\n', saveDir);

