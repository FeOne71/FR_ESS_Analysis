%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New Rack Data 
% Convert Rack Data CSV files to MAT file
% Multi-year processing: 2023, 2024, 2025
% Rack structure: 17 modules in sereies
% Modu structure: 2P 14S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% File Directory
clc; clear; close all;

data_dir  = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\김제 배터리 데이터';
years     = {'2023', '2024', '2025'};
save_dir = 'New';

%% File Traversal

for y = 1:length(years)
    year_str = years{y};
    year_only = str2double(year_str);
    
    % data_dir에서 해당 연도의 파일들을 찾아서 처리
    file_pattern = fullfile(data_dir, sprintf('%s*.csv', year_str));
    file_list = dir(file_pattern);
    
    fprintf('Found %d files for year %s\n', length(file_list), year_str);
    
    if isempty(file_list)
        fprintf('No files found for year %s, skipping...\n', year_str);
        continue;
    end
    
    % 연도별 폴더 생성
    output_year_folder = fullfile(save_dir, year_str);
    fprintf('Creating folder: %s\n', save_dir);
    
    % save_dir 경로를 단계별로 생성
    path_parts = strsplit(save_dir, '\');
    current_path = '';
    for p = 1:length(path_parts)
        if ~isempty(path_parts{p})
            if isempty(current_path)
                current_path = path_parts{p};
            else
                current_path = fullfile(current_path, path_parts{p});
            end
            if ~exist(current_path, 'dir')
                mkdir(current_path);
                fprintf('Created: %s\n', current_path);
            end
        end
    end
    
    fprintf('Creating year folder: %s\n', output_year_folder);
    if ~exist(output_year_folder, 'dir')
        mkdir(output_year_folder);
        fprintf('Created year folder: %s\n', output_year_folder);
    end
    
    % 각 파일을 처리
    for f = 1:length(file_list)
        file_info = file_list(f);
        input_file = fullfile(data_dir, file_info.name);
        
        % 파일명에서 날짜 추출 (YYYYMMDD_KIMJ_LGE_01_01_01 형식)
        date_part = file_info.name(1:8); % YYYYMMDD 부분
        year_month = date_part(1:6); % YYYYMM 부분
        day_str = date_part;
        
        % 월별 폴더 생성
        output_month_folder = fullfile(output_year_folder, year_month);
        if ~exist(output_month_folder, 'dir')
            mkdir(output_month_folder);
        end
        
        % 출력 MAT 파일명
        output_mat_file = fullfile(output_month_folder, sprintf('Raw_%s.mat', day_str));
        
        fprintf('Processing day: %s\n', day_str);
        
        %% CSV File Processing
        try
            % Read CSV file
            T = readtable(input_file, 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
            if isempty(T)
                fprintf('  Warning: Empty table for file %s\n', input_file);
                continue;
            end
            
            header_row = T.Properties.VariableNames;
            wanted_vars_raw = { ...
                'Date_Time', 'Rack_Status', 'Charge_Mode', 'Discharge_Mode', 'CB_Feedback', ...
                'Contactor_Feedback', 'CellBalancing', 'SOC_BMS', 'SOH_BMS', 'DCCurrent', ...
                'DCPower', 'DCchgPowerLimit', 'DCdchgPowerLimit', 'CVavg', 'CVmax', 'CVmin', ...
                'CVmaxNoMBMS', 'CVmaxNoCell', 'CVminNoMBMS', 'CVminNoCell', ...
                'MTavg', 'MTmax', 'MTmin', 'MTmaxNoMBMS', 'MTminNoMBMS'};
            
            selected_cols = [];
            for k = 1:length(wanted_vars_raw)
                idx = find(strcmp(strtrim(header_row), wanted_vars_raw{k}), 1);
                if ~isempty(idx)
                    selected_cols(end+1) = idx;
                end
            end
            
            if isempty(selected_cols)
                fprintf('  Warning: No matching columns in file %s\n', input_file);
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
                data_col = T_sel.(var_name);
                if iscell(data_col)
                    numeric_data = nan(size(data_col));
                    for k = 1:length(data_col)
                        if ~isempty(data_col{k}) && ~isnan(str2double(data_col{k}))
                            numeric_data(k) = str2double(data_col{k});
                        end
                    end
                    T_sel.(var_name) = numeric_data;
                elseif isa(data_col, 'duration')
                    % duration 타입을 초 단위로 변환하되 원본도 저장
                    if strcmp(var_name, 'DateTime') || strcmp(var_name, 'Date_Time')
                        % 원본 Date_Time 저장
                        T_sel.Date_Time = data_col;
                        % 초 단위 변환 데이터 추가
                        T_sel.Date_Time_seconds = seconds(data_col);
                    else
                        T_sel.(var_name) = seconds(data_col);
                    end
                else
                    T_sel.(var_name) = double(T_sel.(var_name));
                end
            end
            
            % Raw 구조체 생성 
            Raw = struct();
            field_names = T_sel.Properties.VariableNames;
            for i = 1:length(field_names)
                field = convertVariableName(field_names{i}, false);
                Raw.(field) = T_sel.(field_names{i});
            end
            
            % MAT 파일로 저장
            save(output_mat_file, 'Raw', '-v7.3');
            fprintf('Raw.mat saved: %s\n', output_mat_file);
            
        catch ME
            fprintf('Error processing %s: %s\n', file_info.name, ME.message);
        end
    end
end

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