%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Ambient Temperature to Existing MAT Files
% This script reads existing MAT files and adds ambient temperature data
% from corresponding PLC files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% Main script starts here
    % Base paths
    mat_base_path = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
    plc_base_path = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS';
    
    % Process years
    years = {'2021'}; %, '2022', '2023'};
    
    for year_idx = 1:length(years)
        year = years{year_idx};
        year_path = fullfile(mat_base_path, year);
        
        if ~exist(year_path, 'dir')
            continue;
        end
        
        % Get month folders
        month_folders = dir(year_path);
        month_folders = month_folders([month_folders.isdir] & ~strcmp({month_folders.name}, '.') & ~strcmp({month_folders.name}, '..'));
        
        for month_idx = 1:length(month_folders)
            month_folder = month_folders(month_idx).name;
            month_path = fullfile(year_path, month_folder);
            
            % Get MAT files in month folder
            mat_files = dir(fullfile(month_path, 'Raw_*.mat'));
            
            for mat_idx = 1:length(mat_files)
                mat_file = mat_files(mat_idx).name;
                mat_path = fullfile(month_path, mat_file);
                
                % Extract date from filename
                date_str = extractBetween(mat_file, 'Raw_', '.mat');
                if isempty(date_str)
                    continue;
                end
                date_str = date_str{1};
                
                fprintf('Processing: %s\n', mat_file);
                
                % Find corresponding PLC file
                plc_file = findPLCFile(plc_base_path, year, date_str);
                if isempty(plc_file)
                    fprintf('  No PLC file found for %s\n', date_str);
                    continue;
                end
                
                % Read ambient temperature from PLC file
                ambient_temp = readAmbientTemperature(plc_file);
                if isempty(ambient_temp)
                    fprintf('  No ambient temperature data found in %s\n', plc_file);
                    continue;
                end
                
                % Load existing MAT file
                try
                    % Debug: Check what's in the MAT file
                    fprintf('  Debug: Loading MAT file: %s\n', mat_path);
                    mat_info = whos('-file', mat_path);
                    fprintf('  Debug: Variables in MAT file: %s\n', strjoin({mat_info.name}, ', '));
                    
                    load(mat_path, 'Raw');
                    
                    % Add ambient temperature to all modules
                    Raw = addAmbientToAllModules(Raw, ambient_temp);
                    
                    % Save updated MAT file
                    save(mat_path, 'Raw');
                    fprintf('  Ambient temperature added successfully\n');
                    
                catch ME
                    fprintf('  Error processing %s: %s\n', mat_file, ME.message);
                end
            end
        end
    end
    
    fprintf('Ambient temperature addition completed!\n');

function plc_file = findPLCFile(plc_base_path, year, date_str)
    % Find PLC file for given date
    year_month = date_str(1:6); % YYYYMM
    day = date_str(7:8); % DD
    
    % Construct PLC file path
    plc_folder = fullfile(plc_base_path, sprintf('%s_KIMJ', year_month), date_str);
    
    if ~exist(plc_folder, 'dir')
        plc_file = [];
        return;
    end
    
    % Find PLC file based on year format
    if strcmp(year, '2021')
        % 21년도 형식: 20210601_LGCHEM_PLC#1
        plc_pattern = sprintf('%s_LGCHEM_PLC#*.csv', date_str);
    else
        % 22,23년도 형식: JXR_BSC_PLC1_20220601
        plc_pattern = sprintf('JXR_BSC_PLC*_%s*.csv', date_str);
    end
    
    % Find PLC file
    plc_files = dir(fullfile(plc_folder, plc_pattern));
    
    % Filter out Config files
    valid_files = [];
    for i = 1:length(plc_files)
        if ~contains(plc_files(i).name, 'Config')
            valid_files(end+1) = i;
        end
    end
    
    if isempty(valid_files)
        plc_file = [];
        return;
    end
    
    plc_files = plc_files(valid_files);
    
    % Use all found files (combine multiple PLC files for one day)
    plc_files_list = {};
    for i = 1:length(plc_files)
        plc_files_list{end+1} = fullfile(plc_folder, plc_files(i).name);
    end
    plc_file = plc_files_list;  % Return list of all PLC files
end

function ambient_temp = readAmbientTemperature(plc_files)
    % Read ambient temperature from multiple PLC files
    try
        all_plc_data = [];
        
        % Process each PLC file
        for file_idx = 1:length(plc_files)
            plc_file = plc_files{file_idx};
            fprintf('    Debug: Reading PLC file %d/%d: %s\n', file_idx, length(plc_files), plc_file);
            
            % Read PLC data with original column names preserved
            plc_data = readtable(plc_file, 'VariableNamingRule', 'preserve');
            
            % Find ambient temperature column
            ambient_col = find(contains(plc_data.Properties.VariableNames, 'Ambient Temperature', 'IgnoreCase', true));
            if isempty(ambient_col)
                fprintf('    Debug: No Ambient Temperature column in %s\n', plc_file);
                continue;
            end
            
            % Find time column
            time_col = find(contains(plc_data.Properties.VariableNames, 'Time', 'IgnoreCase', true));
            if isempty(time_col)
                fprintf('    Debug: No Time column in %s\n', plc_file);
                continue;
            end
            
            % Extract time and ambient temperature data
            time_data = plc_data{:, time_col};
            ambient_data = plc_data{:, ambient_col};
            
            % Convert duration to datetime with date
            if isduration(time_data)
                % Get date from filename (e.g., 20210601_LGCHEM_PLC#1.csv -> 2021-06-01)
                [~, filename, ~] = fileparts(plc_file);
                date_str = extractBetween(filename, 1, 8);  % YYYYMMDD
                if ~isempty(date_str)
                    base_date = datetime(date_str{1}, 'InputFormat', 'yyyyMMdd');
                    time_data = base_date + time_data;  % Add date to duration
                end
            end
            
            % Combine data
            file_data = table(time_data, ambient_data, 'VariableNames', {'Time', 'Ambient_Temperature'});
            
            if isempty(all_plc_data)
                all_plc_data = file_data;
            else
                all_plc_data = vertcat(all_plc_data, file_data);
            end
            
            fprintf('    Debug: Added %d data points from %s\n', height(file_data), plc_file);
        end
        
        ambient_temp = all_plc_data;
        fprintf('    Debug: Total PLC data points: %d\n', height(ambient_temp));
        
    catch ME
        fprintf('  Error reading PLC files: %s\n', ME.message);
        ambient_temp = [];
    end
end

function Raw = addAmbientToAllModules(Raw, ambient_temp)
    % Add ambient temperature to all racks in Raw structure (ver04 format)
    % Synchronize PLC time with Raw time (RBMS time)
    
    % Get all rack fields
    rack_fields = fieldnames(Raw);
    
    for rack_idx = 1:length(rack_fields)
        rack_field = rack_fields{rack_idx};
        if ~startsWith(rack_field, 'Rack')
            continue;
        end
        
        % Get Raw time (RBMS time from ver04)
        raw_time = Raw.(rack_field).Time;
        
        % DEBUG: Check time formats
        fprintf('  Debug: %s - Raw time type: %s, sample: %s\n', rack_field, class(raw_time), string(raw_time(1)));
        fprintf('  Debug: %s - PLC time type: %s, sample: %s\n', rack_field, class(ambient_temp.Time), string(ambient_temp.Time(1)));
        fprintf('  Debug: %s - Raw time length: %d, PLC time length: %d\n', rack_field, length(raw_time), height(ambient_temp));
        
        % DEBUG: Check date mismatch
        raw_date = extractBetween(string(raw_time(1)), 1, 10);  % YYYY-MM-DD
        plc_date = extractBetween(string(ambient_temp.Time(1)), 1, 10);  % YYYY-MM-DD
        if ~strcmp(raw_date, plc_date)
            fprintf('  Debug: DATE MISMATCH! Raw: %s, PLC: %s\n', raw_date, plc_date);
        end
        
        % Synchronize ambient temperature with Raw time
        synced_ambient = syncAmbientWithRawTime(ambient_temp, raw_time);
        
        % Add synchronized ambient temperature to this rack
        Raw.(rack_field).Ambient_Temperature = synced_ambient;
    end
end

function synced_ambient = syncAmbientWithRawTime(ambient_temp, raw_time)
    % Synchronize PLC ambient temperature with Raw time
    % Only use exact time matches, no interpolation - OPTIMIZED VERSION
    
    % Extract PLC time and temperature
    if istable(ambient_temp) && height(ambient_temp) > 0
        plc_time = ambient_temp.Time;
        plc_temp = ambient_temp.Ambient_Temperature;
    else
        synced_ambient = [];
        return;
    end
    
    % Convert PLC time to match Raw time type (cached conversion)
    if isdatetime(raw_time)
        % Raw is datetime, convert PLC to datetime
        if ~isdatetime(plc_time)
            plc_time = datetime(plc_time);
        end
    elseif isstring(raw_time) || ischar(raw_time)
        % Raw is string, convert PLC to string
        if isdatetime(plc_time)
            plc_time = string(plc_time);
        elseif ~isstring(plc_time) && ~ischar(plc_time)
            plc_time = string(plc_time);
        end
    end
    
    % OPTIMIZATION 1: Use ismember for O(n+m) instead of O(n×m)
    [~, matched_indices] = ismember(raw_time, plc_time);
    
    % DEBUG: Check matching results
    fprintf('    Debug: Total raw times: %d, Total PLC times: %d\n', length(raw_time), length(plc_time));
    fprintf('    Debug: Matched indices range: %d to %d\n', min(matched_indices), max(matched_indices));
    fprintf('    Debug: Non-zero matches: %d\n', sum(matched_indices > 0));
    
    % DEBUG: Check which time points are missing
    if length(raw_time) ~= length(plc_time)
        missing_count = length(raw_time) - length(plc_time);
        fprintf('    Debug: Missing %d time points in PLC data\n', missing_count);
        
        % Find missing time points
        missing_indices = find(matched_indices == 0);
        if ~isempty(missing_indices)
            fprintf('    Debug: Missing time points (first 5): ');
            for i = 1:min(5, length(missing_indices))
                fprintf('%s ', string(raw_time(missing_indices(i))));
            end
            if length(missing_indices) > 5
                fprintf('... (and %d more)', length(missing_indices) - 5);
            end
            fprintf('\n');
        end
    end
    
    % OPTIMIZATION 2: Pre-allocate arrays
    valid_matches = matched_indices > 0;
    num_matches = sum(valid_matches);
    
    if num_matches == 0
        fprintf('    Debug: No time matches found!\n');
        synced_ambient = [];
        return;
    end
    
    fprintf('    Debug: Found %d matching times\n', num_matches);
    
    % OPTIMIZATION 3: Vectorized extraction
    matched_temps = plc_temp(matched_indices(valid_matches));
    
    % Return as double array (matching ver04 format)
    synced_ambient = matched_temps;
end
