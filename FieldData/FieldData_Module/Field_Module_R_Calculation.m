%% Field Module R50s Calculation - Module-Level Resistance Calculation
% This script calculates module-level R50s resistance from RARD synchronized data
% and saves the results for later visualization
% 
% Process:
% 1. Load module data from RARDsync folder
% 2. Calculate R50s using dV/dI linear fitting (no moving average)
% 3. Save results to MAT file for visualization

clear; clc; close all;

%% Configuration
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

% Target Rack and Module for calculation
targetRack = 'Rack01';
targetModule = 1:17;  % Module number (1-17) - Process all modules

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
        module_data.(moduleName).(cellName).V_avg_cell = cell(numFiles, 1);  % Average voltage for each day
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
                    cell_result.V_avg = mean(V_cell);  % Average voltage for this day
                    
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
            module_data.(moduleName_result).(cellName).V_avg_cell{f} = cell_result.V_avg;
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
            module_data.(moduleName).(cellName).V_avg = [module_data.(moduleName).(cellName).V_avg_cell{keepIdx}]';
            
            % Sort by date
            [module_data.(moduleName).(cellName).dates, sortIdx] = sort(module_data.(moduleName).(cellName).dates);
            module_data.(moduleName).(cellName).R50s_daily = module_data.(moduleName).(cellName).R50s_daily(sortIdx);
            module_data.(moduleName).(cellName).T_avg = module_data.(moduleName).(cellName).T_avg(sortIdx);
            module_data.(moduleName).(cellName).n_points = module_data.(moduleName).(cellName).n_points(sortIdx);
            module_data.(moduleName).(cellName).r_squared = module_data.(moduleName).(cellName).r_squared(sortIdx);
            module_data.(moduleName).(cellName).V_avg = module_data.(moduleName).(cellName).V_avg(sortIdx);
            
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
            module_data.(moduleName).(cellName).V_avg_cell = {};
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

%% Save calculation results
fprintf('\n=== Saving Calculation Results ===\n');
save(fullfile(saveDir, sprintf('%s_Resistance_Summary.mat', targetRack)), 'module_data', '-v7.3');
fprintf('Calculation results saved to: %s_Resistance_Summary.mat\n', targetRack);

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

fprintf('\n=== Calculation Complete ===\n');
fprintf('Results saved to: %s\n', fullfile(saveDir, sprintf('%s_Resistance_Summary.mat', targetRack)));
fprintf('Run Field_Module_R_Visualization.m to generate visualizations.\n');

