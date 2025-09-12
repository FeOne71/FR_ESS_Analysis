%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert RBMS and RARD CSV files to MAT file - Each Rack as a sub-struct
% RBMS & RARD time synchronization
% Features:
% - RBMS data processing (Rack01-08)
% - RARD data processing (Module-level data)
% - Time synchronization between RBMS and RARD data
% - RARD time as reference for synchronization
% - 8 Rack / 17 Module / 14 Cell
% Date: 2025-09-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. File Directory
clc; clear; close all;

dataFile = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS';
years    = {'202106_KIMJ', '202206_KIMJ', '202306_KIMJ'}; %,'202206_KIMJ','202306_KIMJ'};
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\RARDsync';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end


%% 2. Folder Traversal
for y = 1:length(years)
    year_folder = fullfile(dataFile, years{y});
    if ~exist(year_folder, 'dir')
        continue;
    end

    year_str = extractBefore(years{y}, '_');
    year_only = str2double(year_str(1:4));
    month_only = str2double(year_str(5:6));

    output_year_folder = fullfile(saveDir, year_str(1:4));
    if ~exist(output_year_folder, 'dir')
        mkdir(output_year_folder);
    end

    output_month_folder = fullfile(output_year_folder, year_str(1:6));
    if ~exist(output_month_folder, 'dir')
        mkdir(output_month_folder);
    end

    last_day = eomday(year_only, month_only);

    %% 3. Daily CSV files Traversal (RBMS & RARD)
    for day = 1:last_day
        day_str = sprintf('%s%02d', year_str(1:6), day);
        day_folder = fullfile(year_folder, day_str);
        mat_file = fullfile(output_month_folder, sprintf('Raw_%s.mat', day_str));

        fprintf('Processing day: %s\n', day_str);
        tic;

        if ~exist(day_folder, 'dir')
            fprintf('No files found for %s, skipping...\n', day_str);
            continue;
        end

        %% RARD File Processing
        rard_data = struct();
        rard_file = findRARDFile(day_folder, day_str);
        
        if ~isempty(rard_file)
            fprintf('  Found RARD file: %s\n', rard_file);
            rard_data = processRARDFile(rard_file, day_str);
        else
            fprintf('  No RARD file found for %s\n', day_str);
        end

        %% RBMS File Processing (for synchronization only)
        rbms_data = struct();
        for rackIdx = 1:8
            % 연도에 따른 파일 패턴 결정
            if year_only == 2021
            file_pattern = fullfile(day_folder, sprintf('*RBMS[%02d]*.csv', rackIdx));
            else
                % 2022년 이후: JXR_BSC_Rack1_YYYYMMDD.csv, JXR_BSC_Rack1_YYYYMMDD_1.csv, ...
                file_pattern = fullfile(day_folder, sprintf('*_Rack%d_%s*.csv', rackIdx, day_str));
            end
            
            file_list = dir(file_pattern);
            if isempty(file_list)
                fprintf('  No file: %s\n', file_pattern);
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
                
                % % 디버깅: 모든 헤더 출력
                % fprintf('  Debug: All headers in %s:\n', file_list(f).name);
                % for h = 1:length(header_row)
                %     fprintf('    [%d] %s\n', h, header_row{h});
                % end
                % 
                wanted_vars_raw = { ...
                    'Counter', 'Time', 'Charge', 'Status', 'Active Mdls(ea)', ...
                    'SOC(%)', 'SOH(%)', 'DC Current(A)', 'DC Chg. Current Limit(A)', ...
                    'DC Dchg. Current Limit(A)', 'DC Power(kW)', 'DC Chg. Power Limit(kW)', ...
                    'DC Dchg. Power Limit(kW)', 'Sum. C.V.(V)', 'Average C.V.(V)', ...
                    'Highest C.V.(V)', 'Lowest C.V.(V)', 'Highest C.V. Pos.(MBMS)', ...
                    'Highest C.V. Pos.(Cell)', 'Lowest C.V. Pos.(MBMS)', 'Lowest C.V. Pos.(Cell)', ...
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
                
                % % 디버깅: 찾은 변수들 출력
                % found_vars = header_row(selected_cols);
                % fprintf('  Debug: Found %d variables in %s\n', length(found_vars), file_list(f).name);
                % for v = 1:length(found_vars)
                %     fprintf('    %s\n', found_vars{v});
                % end
                T_sel = T(:, selected_cols);
                clean_var_names = cell(size(selected_cols));
                for i = 1:length(selected_cols)
                    clean_var_names{i} = convertVariableName(header_row{selected_cols(i)});
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
                field = convertVariableName(field_names{i});
                S.(field) = T_all.(field_names{i});
            end
            rbms_data.(rack_field) = S;
        end

        %% RBMS와 RARD 데이터 시간 동기화 및 통합 저장
        Raw = struct();
        if ~isempty(fieldnames(rard_data)) && ~isempty(fieldnames(rbms_data))
            fprintf('  Synchronizing RBMS and RARD data...\n');
            Raw = synchronizeRBMSRARD(rbms_data, rard_data);
        else
            fprintf('  Both RARD and RBMS data required for synchronization\n');
        end

        % 저장
        save(mat_file, 'Raw', '-v7.3');
        elapsed_time = toc;
        fprintf('Raw.mat saved: %s (%.2f seconds)\n', mat_file, elapsed_time);
    end
end

function rard_file = findRARDFile(day_folder, day_str)
    % Function to find RARD file
    % Pattern: *_RARD.csv or *_RARD_YYYYMMDD.csv
    
    % day_folder에서 찾기
    rard_pattern1 = fullfile(day_folder, '*_RARD.csv');
    rard_files1 = dir(rard_pattern1);
    
    rard_pattern2 = fullfile(day_folder, sprintf('*_RARD_%s.csv', day_str));
    rard_files2 = dir(rard_pattern2);
    
    if ~isempty(rard_files1)
        rard_file = fullfile(day_folder, rard_files1(1).name);
    elseif ~isempty(rard_files2)
        rard_file = fullfile(day_folder, rard_files2(1).name);
    else
        rard_file = '';
    end
end

function rard_data = processRARDFile(rard_file, day_str)
    % Function to process RARD file
    [T, header_line, cell_header_line] = readRARDTable(rard_file);
    if isempty(T)
        rard_data = struct();
        return;
    end
    T = processTimeData(T, day_str);
    rard_data = splitByModule(T, rard_file, header_line, cell_header_line);
    fprintf('    RARD data processed: %d Module found\n', length(fieldnames(rard_data)));
end

function [T, header_line, cell_header_line] = readRARDTable(rard_file)
    % Read RARD table
    fid = fopen(rard_file, 'r');
    if fid == -1
        fprintf('    Error: Cannot open RARD file %s\n', rard_file);
        T = [];
        header_line = '';
        return;
    end
    fclose(fid);
    
    % ver04 method: directly check header lines before using readtable
    fid = fopen(rard_file, 'r');
    for i = 1:3
        fgetl(fid);  % Skip first 3 rows
    end
    header_line = fgetl(fid);   % 4th line: module numbers
    cell_header_line = fgetl(fid);   % 5th line: cell headers
    fclose(fid);
    
    T = readtable(rard_file, 'HeaderLines', 4, 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
    
    if isempty(T)
        fprintf('    Warning: Empty RARD table\n');
        return;
    end
    
    
    % Check required columns
    rack_col = find(strcmp(T.Properties.VariableNames, 'Rack No.'));
    time_col = find(strcmp(T.Properties.VariableNames, 'Time'));
    
    if isempty(rack_col) || isempty(time_col)
        fprintf('    Warning: Required columns not found in RARD file\n');
        T = [];
    end
end

function T = processTimeData(T, day_str)
    % Process time data
    base_date = datetime(day_str, 'InputFormat', 'yyyyMMdd');
    time_col = find(strcmp(T.Properties.VariableNames, 'Time'));
    time_data = T.(T.Properties.VariableNames{time_col});
    
    if iscell(time_data) || isstring(time_data)
        time_str = string(time_data);
        if contains(time_str(1), '-') || contains(time_str(1), '/')
            T.Time = datetime(time_str, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        else
            T.Time = base_date + duration(time_str, 'InputFormat', 'HH:mm:ss');
        end
    elseif isnumeric(time_data)
        T.Time = base_date + seconds(time_data);
    end
end

function rard_data = splitByModule(T, rard_file, header_line, cell_header_line)
    % Split data by Rack and Module (optimized version)
    rard_data = struct();
    rack_col = find(strcmp(T.Properties.VariableNames, 'Rack No.'));
    time_col = find(strcmp(T.Properties.VariableNames, 'Time'));
    
    rack_numbers = unique(T.(T.Properties.VariableNames{rack_col}));
    rack_numbers = rack_numbers(~isnan(rack_numbers));
    
    fprintf('    Found %d racks in RARD data: %s\n', length(rack_numbers), mat2str(rack_numbers));
    
    % Process all module data at once (vectorized)
    for rack_num = rack_numbers'
        rack_idx = round(rack_num);
        if rack_idx >= 1 && rack_idx <= 8
            rack_mask = T.(T.Properties.VariableNames{rack_col}) == rack_num;
            rack_data = T(rack_mask, :);
            
            fprintf('    Processing Rack%02d: %d data points, %d columns\n', rack_idx, height(rack_data), width(rack_data));
            
            if ~isempty(rack_data)
                % Process all modules at once
                module_data = extractAllModules(rack_data, rack_idx, time_col, rack_col, header_line, cell_header_line);
                if ~isempty(module_data)
                    % Merge module data directly
                    module_fields = fieldnames(module_data);
                    for m = 1:length(module_fields)
                        rard_data.(module_fields{m}) = module_data.(module_fields{m});
                    end
                end
            end
        end
    end
end

function module_data = extractAllModules(rack_data, rack_idx, time_col, rack_col, header_line, cell_header_line)
    module_data = struct();
    
    time_data = rack_data{:, time_col};
    rack_data_array = rack_data{:, rack_col};
    counter_data = rack_data{:, 1};
    
    module_positions = parseModulePositions(header_line);
    
    cell_positions = parseCellPositions(cell_header_line);
    
    rack_field_name = sprintf('Rack%02d', rack_idx);
    module_data.(rack_field_name) = struct();
    
    for module_num = 1:17
        if isfield(module_positions, sprintf('Module%02d', module_num))
            start_col = module_positions.(sprintf('Module%02d', module_num));
            
            module_field_name = sprintf('Module%02d', module_num);
            
            module_data.(rack_field_name).(module_field_name) = struct();
            module_data.(rack_field_name).(module_field_name).Time = time_data;
            module_data.(rack_field_name).(module_field_name).Rack_No_ = rack_data_array;
            module_data.(rack_field_name).(module_field_name).Counter = counter_data;
            
            for cell_idx = 1:14
                cell_field_name = sprintf('M%d_Cell%d', module_num, cell_idx);
                if isfield(cell_positions, cell_field_name)
                    col_idx = cell_positions.(cell_field_name);
                    if col_idx <= width(rack_data)
                        cell_data = rack_data{:, col_idx};
                        if istable(cell_data)
                            module_data.(rack_field_name).(module_field_name).(cell_field_name) = table2array(cell_data);
                        else
                            module_data.(rack_field_name).(module_field_name).(cell_field_name) = cell_data;
                        end
                    end
                end
            end
        end
    end
end

function module_positions = parseModulePositions(header_line)
    module_positions = struct();
    
    if isempty(header_line)
        fprintf('    Error: Empty header line\n');
        return;
    end
    
    % fprintf('    Debug: Row 4 content: %s\n', header_line);
    
    parts = strsplit(header_line, ',');
    
    for col = 4:length(parts)  % 4열부터 시작
        cell_content = strtrim(parts{col});
        if ~isempty(cell_content)
            % [Module#01] 형태에서 숫자 추출
            match = regexp(cell_content, 'Module#(\d+)', 'tokens');
            if ~isempty(match)
                module_num = str2double(match{1}{1});
                if module_num >= 1 && module_num <= 17
                    module_name = sprintf('Module%02d', module_num);
                    if ~isfield(module_positions, module_name)
                        module_positions.(module_name) = col;
                        % fprintf('    Debug: Found %s at column %d\n', module_name, col);
                    end
                end
            end
        end
    end
    fprintf('    Debug: Total modules found: %d\n', length(fieldnames(module_positions)));
end

function cell_positions = parseCellPositions(cell_header_line)
    cell_positions = struct();
    
    if isempty(cell_header_line)
        fprintf('    Error: Empty cell header line\n');
        return;
    end
    
    % fprintf('    Debug: Row 5 content: %s\n', cell_header_line);
    
    parts = strsplit(cell_header_line, ',');
    
    for col = 4:length(parts)  % 4열부터 시작
        cell_content = strtrim(parts{col});
        if ~isempty(cell_content)
            % Cell#01(V) 형태에서 셀 번호 추출
            match = regexp(cell_content, 'Cell#(\d+)\(V\)', 'tokens');
            if ~isempty(match)
                cell_num = str2double(match{1}{1});
                if cell_num >= 1 && cell_num <= 14
                    % 모듈 번호 계산: 각 모듈은 20개 컬럼을 차지
                    % Module01: 4-23, Module02: 24-43, Module03: 44-63, ...
                    module_num = floor((col - 4) / 20) + 1;
                    if module_num >= 1 && module_num <= 17
                        cell_field_name = sprintf('M%d_Cell%d', module_num, cell_num);
                        cell_positions.(cell_field_name) = col;
                        % fprintf('    Debug: Found %s at column %d\n', cell_field_name, col);
                    end
                end
            end
        end
    end
    fprintf('    Debug: Total cells found: %d\n', length(fieldnames(cell_positions)));
end

function Raw = synchronizeRBMSRARD(rbms_data, rard_data)
    Raw = struct();
    rbms_time_cache = struct();
    for rack_idx = 1:8
        rack_name = sprintf('Rack%02d', rack_idx);
        if isfield(rbms_data, rack_name)
            rbms_rack_data = rbms_data.(rack_name);
            if isfield(rbms_rack_data, 'Time')
                rbms_time_cache.(rack_name) = string(rbms_rack_data.Time);
            end
        end
    end
    
    % Process all modules for each rack
    rard_racks = fieldnames(rard_data);
    
    for rack_idx = 1:length(rard_racks)
        rack_name = rard_racks{rack_idx};  % e.g., Rack01
        rack_data = rard_data.(rack_name);
        
        % Process all modules for the corresponding rack
        if isfield(rbms_data, rack_name) && isfield(rbms_time_cache, rack_name)
            rbms_rack_data = rbms_data.(rack_name);
            rbms_time = rbms_time_cache.(rack_name);
            
            % Create rack structure
            Raw.(rack_name) = struct();
            
            % Process each module
            for module_num = 1:17
                module_name = sprintf('Module%02d', module_num);
                if isfield(rack_data, module_name)
                    rard_module_data = rack_data.(module_name);
                    combined_table = syncModuleDataUltraFast(rard_module_data, rbms_rack_data, rbms_time);
                    Raw.(rack_name).(module_name) = combined_table;
                end
            end
        else
            % Store only RARD data if RBMS data is not available
            Raw.(rack_name) = rack_data;
        end
    end
    
    fprintf('    All modules synchronized!\n');
end

function combined_table = syncModuleDataUltraFast(rard_data, rbms_data, rbms_time)
    if isempty(rbms_time)
        combined_table = [];
        return;
    end
    
    rard_time = string(rard_data.Time);
    
    rard_datetime = datetime(rard_time);
    rbms_datetime = datetime(rbms_time);
    
    start_time = max(min(rbms_datetime), min(rard_datetime));
    end_time = min(max(rbms_datetime), max(rard_datetime));
    
    rard_mask = (rard_datetime >= start_time) & (rard_datetime <= end_time);
    
    if ~any(rard_mask)
        combined_table = rard_data;  % RARD 데이터만이라도 반환
        return;
    end
    
    rard_data_sync = struct();
    rard_data_sync.Time = rard_data.Time(rard_mask);
    rard_data_sync.Rack_No_ = rard_data.Rack_No_(rard_mask);
    rard_data_sync.Counter = rard_data.Counter(rard_mask);
    
    rard_fields = fieldnames(rard_data);
    for f = 1:length(rard_fields)
        field_name = rard_fields{f};
        if startsWith(field_name, 'M') && contains(field_name, '_Cell')
            rard_data_sync.(field_name) = rard_data.(field_name)(rard_mask);
        end
    end
    
    rard_time_sync = string(rard_data_sync.Time);
    [~, matched_indices] = ismember(rard_time_sync, rbms_time);
    matched_indices(matched_indices == 0) = NaN;
    
    matched_count = sum(~isnan(matched_indices));
    % fprintf('    Debug: Matched %d out of %d RARD times with RBMS\n', matched_count, length(rard_time_sync));
    
    rbms_sync = struct();
    rbms_fields = fieldnames(rbms_data);
    valid_idx = ~isnan(matched_indices);
    valid_indices = matched_indices(valid_idx);
    
    for j = 1:length(rbms_fields)
        field_name = rbms_fields{j};
        field_data = rbms_data.(field_name);
        
        if ~strcmp(field_name, 'Time')
            if isnumeric(field_data)
                matched_data = nan(size(matched_indices));
                if any(valid_idx)
                    matched_data(valid_idx) = field_data(valid_indices);
                end
                rbms_sync.(field_name) = matched_data;
            else
                rbms_sync.(field_name) = field_data;
            end
        end
    end
    
    combined_table = struct();
    combined_table.Time = rard_data_sync.Time;
    combined_table.Rack_No_ = rard_data_sync.Rack_No_;
    combined_table.Counter = rard_data_sync.Counter;
    
    rard_fields = fieldnames(rard_data_sync);
    for f = 1:length(rard_fields)
        field_name = rard_fields{f};
        if startsWith(field_name, 'M') && contains(field_name, '_Cell')
            combined_table.(field_name) = rard_data_sync.(field_name);
        end
    end
    
    rbms_sync_fields = fieldnames(rbms_sync);
    % fprintf('    Debug: RBMS sync fields found: %d\n', length(rbms_sync_fields));
    for k = 1:length(rbms_sync_fields)
        field_name = rbms_sync_fields{k};
        rbms_field_name = sprintf('RBMS_%s', field_name);
        combined_table.(rbms_field_name) = rbms_sync.(field_name);
        % fprintf('    Debug: Added %s to combined_table\n', rbms_field_name);
    end
end


function clean_name = convertVariableName(original_name)
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
    clean_name = strrep(clean_name, 'Pos.(MBMS)', '_PosMBMS');
    clean_name = strrep(clean_name, 'Pos.(Cell)', '_PosCell');
    clean_name = regexprep(clean_name, '[^a-zA-Z0-9_]', '');
end

