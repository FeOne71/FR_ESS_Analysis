%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert RBMS CSV files to MAT file - Each Rack as a sub-struct
% Date: 2025-07-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. File Directory
clc; clear; close all;

data_dir = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS';
years    = {'202106_KIMJ'}; %,'202206_KIMJ','202306_KIMJ'};
% save_dir = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\raw2mat_ver04\Rack_raw2mat';
save_dir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';

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

    %% 3. Daily CSV files Traversal (RBMS)
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
            file_pattern = fullfile(day_folder, sprintf('*RBMS[%02d]*.csv', rackIdx));
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

        % 저장
        save(mat_file, 'Raw', '-v7.3');
        fprintf('Raw.mat saved: %s\n', mat_file);
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

