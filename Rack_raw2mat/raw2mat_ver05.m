%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert RBMS CSV files to MAT file with Ambient Temperature - Integrated Version
% Combines ver04 RBMS processing with add_ambient_temperature PLC processing
% Date: 2025-01-XX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. File Directory
clc; clear; close all;

data_dir = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS';
save_dir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
years    =  {'202106_KIMJ'}; %,'202306_KIMJ'};    %   {'202106_KIMJ'}; %,'202206_KIMJ','202306_KIMJ'};

%% 2. Folder Traversal
for y = 1:length(years)
    year_folder = fullfile(data_dir, years{y});
    if ~exist(year_folder, 'dir')
        continue;
    end

    year_str = extractBefore(years{y}, '_');
    year_only = str2double(year_str(1:4));
    month_only = str2double(year_str(5:6));

    output_year_folder = fullfile(save_dir, year_str(1:4));
    if ~exist(output_year_folder, 'dir')
        mkdir(output_year_folder);
    end

    output_month_folder = fullfile(output_year_folder, year_str(1:6));
    if ~exist(output_month_folder, 'dir')
        mkdir(output_month_folder);
    end

    last_day = eomday(year_only, month_only);

    %% 3. Daily CSV files Traversal (RBMS + PLC)
    for day = 1:last_day
        day_str = sprintf('%s%02d', year_str(1:6), day);
        day_folder = fullfile(year_folder, day_str);
        mat_file = fullfile(output_month_folder, sprintf('Raw_%s.mat', day_str));

        fprintf('Processing day: %s\n', day_str);

        if ~exist(day_folder, 'dir')
            fprintf('No files found for %s, skipping...\n', day_str);
            continue;
        end

        %% RBMS File Processing
        Raw = struct();
        for rackIdx = 1:8
            % 2022/2023년도 형식: JXR_BSC_Rack1_20220601.csv
            if str2double(year_str(1:4)) >= 2022
                file_pattern = fullfile(day_folder, sprintf('JXR_BSC_Rack%d_%s*.csv', rackIdx, day_str));
            else
                % 2021년도 형식: *RBMS[01]*.csv
                file_pattern = fullfile(day_folder, sprintf('*RBMS[%02d]*.csv', rackIdx));
            end
            file_list = dir(file_pattern);
            if isempty(file_list)
                fprintf('  No RBMS file: %s\n', file_pattern);
                continue;
            end
            T_cells = cell(1, length(file_list));
            cell_idx = 1;
            for f = 1:length(file_list)
                rbms_file = fullfile(day_folder, file_list(f).name);
                fid = fopen(rbms_file, 'r');
                for i = 1:11
                    fgetl(fid);
                end
                header_line = fgetl(fid);   % 12th line: variable names
                fclose(fid);
                T = readtable(rbms_file, 'HeaderLines', 11, 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
                if isempty(T)
                    fprintf('  Warning: Empty table for file %s\n', rbms_file);
                    continue;
                end
                header_row = T.Properties.VariableNames;
                wanted_vars_raw = { ...
                    'Counter', 'Time', 'Charge', 'Status', 'Active Mdls(ea)', ...
                    'SOC(%)', 'SOH(%)', 'DC Current(A)', 'DC Chg. Current Limit(A)', ...
                    'DC Dchg. Current Limit(A)', 'DC Power(kW)', 'DC Chg. Power Limit(kW)', ...
                    'DC Dchg. Power Limit(kW)', 'Sum. C.V.(V)', 'Average C.V.(V)', ...
                    'Highest C.V.(V)', 'Lowest C.V.(V)', 'Highest C.V. Pos(MBMS)', ...
                    'Highest C.V. Pos(Cell)', 'Lowest C.V. Pos(MBMS)', 'Lowest C.V. Pos(Cell)', ...
                    'Average M.T.(oC)', 'Highest M.T.(oC)', 'Lowest M.T.(oC)'};
                selected_cols = [];
                for k = 1:length(wanted_vars_raw)
                    idx = find(strcmp(strtrim(header_row), wanted_vars_raw{k}), 1);
                    if ~isempty(idx)
                        selected_cols(end+1) = idx;
                    end
                end
                if isempty(selected_cols)
                    fprintf('  Warning: No matching columns in file %s\n', rbms_file);
                    continue;
                end
                T_sel = T(:, selected_cols);
                clean_var_names = cell(size(selected_cols));
                for i = 1:length(selected_cols)
                    clean_var_names{i} = convertVariableName(header_row{selected_cols(i)}, false);
                end
                T_sel.Properties.VariableNames = clean_var_names;
                for j = 1:length(clean_var_names)
                    var_name = clean_var_names{j};
                    if strcmpi(var_name, 'Time')
                        T_sel.(var_name) = string(T_sel.(var_name));
                    elseif ismember(var_name, {'Charge', 'Status'})
                        T_sel.(var_name) = string(T_sel.(var_name));
                    else
                        data_col = T_sel.(var_name);
                        if iscell(data_col)
                            numeric_data = nan(size(data_col));
                            for k = 1:length(data_col)
                                if ~isempty(data_col{k}) && ~isnan(str2double(data_col{k}))
                                    numeric_data(k) = str2double(data_col{k});
                                end
                            end
                            T_sel.(var_name) = numeric_data;
                        else
                            T_sel.(var_name) = double(T_sel.(var_name));
                        end
                    end
                end
                T_cells{cell_idx} = T_sel;
                cell_idx = cell_idx + 1;
            end
            % cell array에서 빈 셀 제거
            T_cells = T_cells(~cellfun('isempty',T_cells));
            if isempty(T_cells)
                continue;
            end
            T_all = vertcat(T_cells{:});
            rack_field = sprintf('Rack%02d', rackIdx);
            S = struct();
            field_names = T_all.Properties.VariableNames;
            for i = 1:length(field_names)
                field = convertVariableName(field_names{i}, false);
                S.(field) = T_all.(field_names{i});
            end
            Raw.(rack_field) = S;
        end

        %% PLC Temperatures Processing
        fprintf('  Processing PLC Temperatures...\n');
        
        % Find PLC files
        plc_files = findPLCFiles(data_dir, year_str, day_str);
        if ~isempty(plc_files)
            % Read temperatures from PLC files
            plc_temp_data = readPLCTemperatures(plc_files);
            if ~isempty(plc_temp_data)
                % Add temperatures to all racks
                Raw = addPLCTempsToAllRacks(Raw, plc_temp_data);
                fprintf('  PLC temperatures added successfully\n');
            else
                fprintf('  Warning: No PLC temperature data found\n');
            end
        else
            fprintf('  Warning: No PLC files found for %s\n', day_str);
        end

        % 저장
        save(mat_file, 'Raw', '-v7.3');
        fprintf('Raw.mat saved: %s\n', mat_file);
    end
end

%% Helper Functions

function clean_name = convertVariableName(original_name, is_total)
    clean_name = original_name;
    clean_name = strrep(clean_name, ' ', '');
    clean_name = strrep(clean_name, '.', '');
    clean_name = strrep(clean_name, '(%)', 'Pct');
    clean_name = strrep(clean_name, '(A)', '_A');
    clean_name = strrep(clean_name, '(V)', '_V');
    clean_name = strrep(clean_name, '(kW)', '_kW');
    clean_name = strrep(clean_name, '(Wh)', '_Wh');
    clean_name = strrep(clean_name, '(mOhm)', '_mOhm');
    clean_name = strrep(clean_name, '(mA)', '_mA');
    clean_name = strrep(clean_name, '(oC)', '_degC');
    clean_name = regexprep(clean_name, '[^a-zA-Z0-9_]', '');
    if is_total
        clean_name = ['Total_' clean_name];
    end
    if isempty(clean_name) || isstrprop(clean_name(1), 'digit')
        clean_name = ['Var_' clean_name];
    end
end

function plc_files = findPLCFiles(data_dir, year_str, day_str)
    % Find PLC files for given date
    year_month = year_str(1:6); % YYYYMM
    
    % Construct PLC file path
    plc_folder = fullfile(data_dir, sprintf('%s_KIMJ', year_month), day_str);
    
    if ~exist(plc_folder, 'dir')
        plc_files = {};
        return;
    end
    
    % Find PLC file based on year format (exclude Config files)
    if str2double(year_str(1:4)) == 2021
        % 21년도 형식: 20210601_LGCHEM_PLC#1
        plc_pattern = sprintf('%s_LGCHEM_PLC#*.csv', day_str);
    else
        % 22,23년도 형식: JXR_BSC_PLC1_20220601
        plc_pattern = sprintf('JXR_BSC_PLC*_%s*.csv', day_str);
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
        plc_files = {};
        return;
    end
    
    plc_files = plc_files(valid_files);
    
    % Convert to full paths
    plc_files_list = {};
    for i = 1:length(plc_files)
        plc_files_list{end+1} = fullfile(plc_folder, plc_files(i).name);
    end
    plc_files = plc_files_list;
end

function temp_data = readPLCTemperatures(plc_files)
    % Read ambient and battery temperature from multiple PLC files
    try
        all_plc_data = [];
        
        % Process each PLC file
        for file_idx = 1:length(plc_files)
            plc_file = plc_files{file_idx};
            fprintf('    Reading PLC file %d/%d: %s\n', file_idx, length(plc_files), plc_file);
            
            % Read PLC data with original column names preserved
            plc_data = readtable(plc_file, 'VariableNamingRule', 'preserve');
            
            % Find ambient temperature column
            ambient_col = find(contains(plc_data.Properties.VariableNames, 'Ambient Temperature', 'IgnoreCase', true));
            if isempty(ambient_col)
                fprintf('    Warning: No Ambient Temperature column in %s\n', plc_file);
                continue;
            end
            
            battemp_col = find(contains(plc_data.Properties.VariableNames, 'Battery Temperature', 'IgnoreCase', true));
            if isempty(battemp_col)
                fprintf('    Warning: No Battery Temperature column in %s\n', plc_file);
                continue;
            end 

            % Find time column
            time_col = find(contains(plc_data.Properties.VariableNames, 'Time', 'IgnoreCase', true));
            if isempty(time_col)
                fprintf('    Warning: No Time column in %s\n', plc_file);
                continue;
            end
            
            % Extract time and ambient temperature data
            time_data = plc_data{:, time_col};
            ambient_data = plc_data{:, ambient_col};
            battemp_data = plc_data{:, battemp_col};        
            
            % Convert duration to datetime with date
            if isduration(time_data)
                % Get date from filename
                [~, filename, ~] = fileparts(plc_file);
                
                % Extract date based on file naming convention
                if contains(filename, 'JXR_BSC_')
                    % 2022/2023년도 형식: JXR_BSC_PLC1_20220601 -> 20220601
                    date_match = regexp(filename, '(\d{8})', 'tokens');
                    if ~isempty(date_match)
                        date_str = date_match{1}{1};
                        base_date = datetime(date_str, 'InputFormat', 'yyyyMMdd');
                        time_data = base_date + time_data;  % Add date to duration
                    end
                else
                    % 2021년도 형식: 20210601_LGCHEM_PLC#1 -> 20210601
                    date_str = extractBetween(filename, 1, 8);  % YYYYMMDD
                    if ~isempty(date_str)
                        base_date = datetime(date_str{1}, 'InputFormat', 'yyyyMMdd');
                        time_data = base_date + time_data;  % Add date to duration
                    end
                end
            end
            
            % Combine data
            file_data = table(time_data, ambient_data, battemp_data, 'VariableNames', {'Time', 'Ambient_Temperature', 'Battery_Temperature'});
            
            if isempty(all_plc_data)
                all_plc_data = file_data;
            else
                all_plc_data = vertcat(all_plc_data, file_data);
            end
            
            fprintf('    Added %d data points from %s\n', height(file_data), plc_file);
        end
        
        temp_data = all_plc_data;
        if ~isempty(temp_data)
            fprintf('    Total PLC data points: %d\n', height(temp_data));
        else
            fprintf('    Total PLC data points: 0\n');
        end
        
    catch ME
        fprintf('  Error reading PLC files: %s\n', ME.message);
        temp_data = [];
    end
end

function Raw = addPLCTempsToAllRacks(Raw, plc_temp_data)
    % Add PLC temperatures to all racks in Raw structure
    % Synchronize PLC time with Raw time (RBMS time)
    
    % Get all rack fields
    rack_fields = fieldnames(Raw);
    
    for rack_idx = 1:length(rack_fields)
        rack_field = rack_fields{rack_idx};
        if ~startsWith(rack_field, 'Rack')
            continue;
        end
        
        % Get Raw time (RBMS time)
        raw_time = Raw.(rack_field).Time;
        
        % Synchronize temperatures with Raw time
        [synced_ambient, synced_batt] = syncPLCTempsWithRawTime(plc_temp_data, raw_time);
        
        % Add synchronized temperatures to this rack
        Raw.(rack_field).Ambient_Temperature = synced_ambient;
        Raw.(rack_field).Battery_Temperature = synced_batt;
    end
end

function [synced_ambient, synced_batt] = syncPLCTempsWithRawTime(plc_temp_data, raw_time)
    % Synchronize PLC temperatures with Raw time
    % Only use exact time matches, no interpolation - OPTIMIZED VERSION
    
    % Initialize output arrays with NaN so they match raw_time length
    synced_ambient = nan(length(raw_time), 1);
    synced_batt = nan(length(raw_time), 1);
    
    % Extract PLC time and temperature
    if istable(plc_temp_data) && height(plc_temp_data) > 0
        plc_time = plc_temp_data.Time;
        plc_ambient = plc_temp_data.Ambient_Temperature;
        plc_batt = plc_temp_data.Battery_Temperature;
    else
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
    
    % OPTIMIZATION: Use ismember for O(n+m) instead of O(n×m)
    [~, matched_indices] = ismember(raw_time, plc_time);
    
    % Process matches
    valid_matches = matched_indices > 0;
    num_matches = sum(valid_matches);
    
    if num_matches > 0
        % Vectorized assignment to the correct positions matching raw_time
        match_idx = matched_indices(valid_matches);
        synced_ambient(valid_matches) = plc_ambient(match_idx);
        synced_batt(valid_matches) = plc_batt(match_idx);
    end
end
    