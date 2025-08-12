%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert KIMJ ESS CSV file to MAT file - All variablles included
% Time sync: BSC, PLC interp1
% Date: 2025-07-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. File Directory
clc; clear; close all;

data_dir = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS';
years    = {'202106_KIMJ','202206_KIMJ','202306_KIMJ'};  
save_dir = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\raw2mat_ver04';

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

    %% 3. Daily
    for day = 1:last_day
        day_str = sprintf('%s%02d', year_str(1:6), day);  
        day_folder = fullfile(year_folder, day_str);  
        mat_file = fullfile(output_month_folder, sprintf('Raw_%s.mat', day_str));

        fprintf('Processing day: %s\n', day_str);
        
        % 날짜 폴더 확인
        if ~exist(day_folder, 'dir')
            fprintf('No files found for %s, skipping...\n', day_str);
            continue;
        end

        %% 2-1. 파일 찾기
        % BSC 파일 찾기
        bsc_files = find_files(day_folder, day_str, year_only, 'BSC');
        
        % Total BSC 파일 찾기
        total_files = find_files(day_folder, day_str, year_only, 'TOTAL');
        
        % Online BSC 파일 찾기
        online_files = find_files(day_folder, day_str, year_only, 'ONLINE');
        
        % PLC 파일 찾기
        plc_files = find_files(day_folder, day_str, year_only, 'PLC');
        
        % 파일이 모두 없으면 다음 날짜로 넘어감
        if isempty(bsc_files) && isempty(total_files) && isempty(online_files) && isempty(plc_files)
            fprintf('No files found for %s, skipping...\n', day_str);
            continue;
        end

        %% 3. 구조체 초기화
        Raw = struct();
        
        %% 4. BSC 파일 처리 - 새 함수 사용
        [Raw.BSC_Time, Raw.BSC_Charge, Raw.BSC_Status] = processBSCFiles(bsc_files);
        
        %% 5. Total BSC 데이터 처리 (readtable 사용)
        
        % 변수 초기화
        Raw.Total_Average_SOC    = [];
        Raw.Total_Highest_SOC    = [];
        Raw.Total_Lowest_SOC     = [];
        Raw.Total_Average_CV_Sum = [];
        Raw.Total_Average_SOH    = [];
        Raw.Total_Average_CV     = [];
        Raw.Total_Average_MT     = [];
        Raw.Total_Highest_MT     = [];
        Raw.Total_Lowest_MT      = [];
        
        for i = 1:length(total_files)
            if ~isempty(total_files)
                fprintf('  Processing Total file %d/%d\n', i, length(total_files));
            end
            try
                % 파일 불러오기
                opts = detectImportOptions(total_files{i});
                opts.VariableNamingRule = 'preserve';  % 원본 컬럼명 유지
                
                % 필요한 변수 목록
                vars_needed = {'AverageSOC____1', 'HighestSOC____1', 'LowestSOC____1', ...
                              'AverageC_V_Sum_V__1', 'AverageSOH____1', 'AverageC_V__V__1', ...
                              'AverageM_T__oC__1', 'HighestM_T__oC__1', 'LowestM_T__oC__1'};
                
                % 변수 이름이 존재하는지 확인
                vars_found = intersect(opts.VariableNames, vars_needed);
                
                if ~isempty(vars_found)
                    opts.SelectedVariableNames = vars_found;
                    T = readtable(total_files{i}, opts);
                    
                    % 각 변수 추출 및 배열에 추가
                    for j = 1:length(vars_found)
                        var_name = vars_found{j};
                        switch var_name
                            case 'AverageSOC____1'
                                Raw.Total_Average_SOC = [Raw.Total_Average_SOC; double(T.AverageSOC____1)];
                            case 'HighestSOC____1'
                                Raw.Total_Highest_SOC = [Raw.Total_Highest_SOC; double(T.HighestSOC____1)];
                            case 'LowestSOC____1'
                                Raw.Total_Lowest_SOC = [Raw.Total_Lowest_SOC; double(T.LowestSOC____1)];
                            case 'AverageC_V_Sum_V__1'
                                Raw.Total_Average_CV_Sum = [Raw.Total_Average_CV_Sum; double(T.AverageC_V_Sum_V__1)];
                            case 'AverageSOH____1'
                                Raw.Total_Average_SOH = [Raw.Total_Average_SOH; double(T.AverageSOH____1)];
                            case 'AverageC_V__V__1'
                                Raw.Total_Average_CV = [Raw.Total_Average_CV; double(T.AverageC_V__V__1)];
                            case 'AverageM_T__oC__1'
                                Raw.Total_Average_MT = [Raw.Total_Average_MT; double(T.AverageM_T__oC__1)];
                            case 'HighestM_T__oC__1'
                                Raw.Total_Highest_MT = [Raw.Total_Highest_MT; double(T.HighestM_T__oC__1)];
                            case 'LowestM_T__oC__1'
                                Raw.Total_Lowest_MT = [Raw.Total_Lowest_MT; double(T.LowestM_T__oC__1)];
                        end
                    end
                end
            catch ME
                fprintf('Warning: Error processing Total file %s: %s\n', total_files{i}, ME.message);
            end
        end
        
        %% 6. Online BSC 데이터 처리 (readtable 사용)
        
        % 변수 초기화
        Raw.Online_DC_Current = [];
        Raw.Online_Highest_DC = [];
        Raw.Online_Lowest_DC = [];
        Raw.Online_DC_Power = [];
        
        for i = 1:length(online_files)
            if ~isempty(online_files)
                fprintf('  Processing Online file %d/%d\n', i, length(online_files));
            end
            try
                % 파일 불러오기
                opts = detectImportOptions(online_files{i});
                opts.VariableNamingRule = 'preserve';  % 원본 컬럼명 유지
                
                % 필요한 변수 목록
                vars_needed = {'DCCurrent_A_', 'HighestDCCurrent_A_', 'LowestDCCurrent_A_', 'DC Power(kW)'};
                
                % 변수 이름이 존재하는지 확인
                vars_found = intersect(opts.VariableNames, vars_needed);
                
                if ~isempty(vars_found)
                    opts.SelectedVariableNames = vars_found;
                    T = readtable(online_files{i}, opts);
                    
                    % 각 변수 추출 및 배열에 추가
                    for j = 1:length(vars_found)
                        var_name = vars_found{j};
                        switch var_name
                            case 'DCCurrent_A_'
                                Raw.Online_DC_Current = [Raw.Online_DC_Current; double(T.DCCurrent_A_)];
                            case 'HighestDCCurrent_A_'
                                Raw.Online_Highest_DC = [Raw.Online_Highest_DC; double(T.HighestDCCurrent_A_)];
                            case 'LowestDCCurrent_A_'
                                Raw.Online_Lowest_DC = [Raw.Online_Lowest_DC; double(T.LowestDCCurrent_A_)];
                            case 'DC Power(kW)'
                                Raw.Online_DC_Power = [Raw.Online_DC_Power; double(T.("DC Power(kW)"))];
                        end
                    end
                end
            catch ME
                fprintf('Warning: Error processing Online file %s: %s\n', online_files{i}, ME.message);
            end
        end
        
        %% 7. PLC 데이터 처리 (readtable 사용)
        
        % 변수 초기화
        Raw.Plc_Time = datetime.empty;
        Raw.Plc_Battery_Temperature = [];
        Raw.Plc_Humidity = [];
        Raw.Plc_Ambient_Temperature = [];
        
        % 기준 날짜 설정
        baseDate = datetime(day_str, 'InputFormat', 'yyyyMMdd');
        
        for i = 1:length(plc_files)
            if ~isempty(plc_files)
                fprintf('  Processing PLC file %d/%d\n', i, length(plc_files));
            end
            try
                % 파일 불러오기
                opts = detectImportOptions(plc_files{i});
                opts.VariableNamingRule = 'preserve';  % 원본 컬럼명 유지
                
                % 열 이름 소문자로 변환하여 검색 (대소문자 구분 없이)
                varLower = lower(opts.VariableNames);
                
                % 필요한 열 찾기
                time_idx = find(contains(varLower, 'time'));
                batt_temp_idx = find(contains(varLower, 'battery') & contains(varLower, 'temp'));
                humid_idx = find(contains(varLower, 'humid'));
                amb_temp_idx = find(contains(varLower, 'ambient') & contains(varLower, 'temp'));
                
                % 최소한 시간 컬럼은 있어야 함
                if isempty(time_idx)
                    continue;
                end
                
                % 필요한 컬럼만 선택
                selected_cols = unique([time_idx(1), ...
                                      if_not_empty(batt_temp_idx, 1), ...
                                      if_not_empty(humid_idx, 1), ...
                                      if_not_empty(amb_temp_idx, 1)]);
                opts.SelectedVariableNames = opts.VariableNames(selected_cols);
                
                % 테이블 읽기
                T = readtable(plc_files{i}, opts);
                
                % 시간 데이터 처리
                time_var = opts.VariableNames{time_idx(1)};
                time_data = T.(time_var);
                
                % 시간 변환
                try
                    % 시간 데이터 형식에 따른 처리
                    if isduration(time_data)
                        file_time = baseDate + time_data;
                    elseif iscell(time_data) || isstring(time_data)
                        % 전체 날짜+시간 형식인지 확인
                        sample = string(time_data(1));
                        if contains(sample, '-') || contains(sample, '/')
                            file_time = datetime(time_data, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
                        else
                            % 시간만 있는 형식
                            file_time = baseDate + duration(time_data, 'InputFormat', 'HH:mm:ss');
                        end
                    elseif isnumeric(time_data)
                        % 숫자 형식인 경우
                        file_time = baseDate + seconds(time_data);
                    else
                        continue;
                    end
                    
                    % 배터리 온도 추출
                    temp_data = [];
                    if ~isempty(batt_temp_idx)
                        temp_var = opts.VariableNames{batt_temp_idx(1)};
                        temp_data = double(T.(temp_var));
                    else
                        % 비어있는 경우 NaN으로 채움
                        temp_data = nan(height(T), 1);
                    end
                    
                    % 습도 추출
                    humid_data = [];
                    if ~isempty(humid_idx)
                        humid_var = opts.VariableNames{humid_idx(1)};
                        humid_data = double(T.(humid_var));
                    else
                        humid_data = nan(height(T), 1);
                    end
                    
                    % 주변 온도 추출
                    amb_temp_data = [];
                    if ~isempty(amb_temp_idx)
                        amb_temp_var = opts.VariableNames{amb_temp_idx(1)};
                        amb_temp_data = double(T.(amb_temp_var));
                    else
                        amb_temp_data = nan(height(T), 1);
                    end
                    
                    % 배열에 추가
                    Raw.Plc_Time = [Raw.Plc_Time; file_time];
                    Raw.Plc_Battery_Temperature = [Raw.Plc_Battery_Temperature; temp_data];
                    Raw.Plc_Humidity = [Raw.Plc_Humidity; humid_data];
                    Raw.Plc_Ambient_Temperature = [Raw.Plc_Ambient_Temperature; amb_temp_data];
                catch ME
                    fprintf('Warning: Time conversion error in PLC file %s: %s\n', plc_files{i}, ME.message);
                end
            catch ME
                fprintf('Warning: Error processing PLC file %s: %s\n', plc_files{i}, ME.message);
            end
        end
        
        % PLC 데이터 정렬 및 중복 제거
        if ~isempty(Raw.Plc_Time)
            % 정렬
            [Raw.Plc_Time, idx] = sort(Raw.Plc_Time);
            
            if ~isempty(Raw.Plc_Battery_Temperature)
                Raw.Plc_Battery_Temperature = Raw.Plc_Battery_Temperature(idx);
            end
            
            if ~isempty(Raw.Plc_Humidity)
                Raw.Plc_Humidity = Raw.Plc_Humidity(idx);
            end
            
            if ~isempty(Raw.Plc_Ambient_Temperature)
                Raw.Plc_Ambient_Temperature = Raw.Plc_Ambient_Temperature(idx);
            end
            
            % 중복 제거
            [Raw.Plc_Time, ia, ~] = unique(Raw.Plc_Time);
            
            if ~isempty(Raw.Plc_Battery_Temperature)
                Raw.Plc_Battery_Temperature = Raw.Plc_Battery_Temperature(ia);
            end
            
            if ~isempty(Raw.Plc_Humidity)
                Raw.Plc_Humidity = Raw.Plc_Humidity(ia);
            end
            
            if ~isempty(Raw.Plc_Ambient_Temperature)
                Raw.Plc_Ambient_Temperature = Raw.Plc_Ambient_Temperature(ia);
            end
        end
        
        %% 8. 시간 동기화
        
        % 기본값 초기화
        Raw.sync_Time = datetime.empty;
        Raw.sync_Time_num = [];
        Raw.Plc_Battery_Temperature_sync = [];
        Raw.Plc_Humidity_sync = [];
        Raw.Plc_Ambient_Temperature_sync = [];
        
        % BSC와 PLC 데이터가 모두 있는 경우에만 동기화
        if ~isempty(Raw.BSC_Time) && ~isempty(Raw.Plc_Time)
            % 시간 범위 결정
            start_idx = max(Raw.BSC_Time(1), Raw.Plc_Time(1));
            end_idx = min(Raw.BSC_Time(end), Raw.Plc_Time(end));
            
            % 범위 내 BSC 시간 선택
            bsc_mask = (Raw.BSC_Time >= start_idx) & (Raw.BSC_Time <= end_idx);
            
            if sum(bsc_mask) > 0
                Raw.sync_Time = Raw.BSC_Time(bsc_mask);
                Raw.sync_Time_num = datenum(Raw.sync_Time);
                
                % 보간 수행 (extrapolation 제거)
                if length(Raw.sync_Time) >= 2 && length(Raw.Plc_Time) >= 2
                    try
                        plc_time_num = datenum(Raw.Plc_Time);
                        
                        % 배터리 온도 보간
                        if ~isempty(Raw.Plc_Battery_Temperature)
                            Raw.Plc_Battery_Temperature_sync = interp1(plc_time_num, Raw.Plc_Battery_Temperature, Raw.sync_Time_num, 'linear');
                        end
                        
                        % 습도 보간
                        if ~isempty(Raw.Plc_Humidity)
                            Raw.Plc_Humidity_sync = interp1(plc_time_num, Raw.Plc_Humidity, Raw.sync_Time_num, 'linear');
                        end
                        
                        % 주변 온도 보간
                        if ~isempty(Raw.Plc_Ambient_Temperature)
                            Raw.Plc_Ambient_Temperature_sync = interp1(plc_time_num, Raw.Plc_Ambient_Temperature, Raw.sync_Time_num, 'linear');
                        end
                    catch
                        % 오류 무시하고 계속 진행
                    end
                end
            end
        end
        
        %% 9. 저장
        save(mat_file, 'Raw', '-v7.3');
        fprintf('Raw.mat saved: %s\n', mat_file);
    end
end

%% BSC 파일 처리 함수 
function [BSC_Time, BSC_Charge, BSC_Status] = processBSCFiles(bsc_files)
    % 초기화
    allBSC_data = table();
    selectedCols_bsc = [2,4,5];
    
    if isempty(bsc_files)
        BSC_Time = datetime.empty;
        BSC_Charge = string.empty;
        BSC_Status = string.empty;
        return;
    end
    
    % 첫 번째 파일에서 헤더 읽기
    C = readcell(bsc_files{1}, 'Delimiter', ',');
    header = C(5,:);
    header_sel = strtrim(header(selectedCols_bsc));
    varNames = matlab.lang.makeValidName(header_sel);
    
    % 모든 BSC 파일 처리
    for i = 1:length(bsc_files)
        try
            C = readcell(bsc_files{i}, 'Delimiter', ',');
            
            % 데이터 행 추출 (6행부터)
            dataCells = C(6:end,:);
            
            % 테이블로 변환
            T_part = cell2table(dataCells(:, selectedCols_bsc), 'VariableNames', varNames);
            
            % 데이터 타입 변환
            T_part.Charge = string(T_part.Charge);
            T_part.Status = string(T_part.Status);
            T_part.Time = datetime(T_part.Time(:), 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
            
            % 전체 데이터에 추가
            allBSC_data = [allBSC_data; T_part];
        catch e
            fprintf('Error processing BSC file %s: %s\n', bsc_files{i}, e.message);
        end
    end
    
    % 출력 변수 설정
    if ~isempty(allBSC_data)
        BSC_Time = allBSC_data.Time;
        BSC_Charge = allBSC_data.Charge;
        BSC_Status = allBSC_data.Status;
    else
        BSC_Time = datetime.empty;
        BSC_Charge = string.empty;
        BSC_Status = string.empty;
    end
end

%% 파일 패턴에 따라 파일 목록 찾기
function files = find_files(folder, day_str, year, file_type)
    files = {};
    
    if ~exist(folder, 'dir')
        return;
    end
    
    all_csv_files = dir(fullfile(folder, '*.csv'));
    if isempty(all_csv_files)
        return;
    end
    
    all_filenames = {all_csv_files.name};
    
    switch upper(file_type)
        case 'BSC'
            if year == 2021
                % 첫 번째 BSC 파일 (접미사 없음)
                pattern1 = sprintf('%s_LGCHEM_BSC.csv', day_str);
                idx = find(strcmpi(all_filenames, pattern1));
                if ~isempty(idx)
                    files{end+1} = fullfile(folder, all_filenames{idx(1)});
                end
                
                % 이후 파일: _1, _2, ... _7
                for i = 1:7
                    pattern = sprintf('%s_LGCHEM_BSC_%d.csv', day_str, i);
                    idx = find(strcmpi(all_filenames, pattern));
                    if ~isempty(idx)
                        files{end+1} = fullfile(folder, all_filenames{idx(1)});
                    end
                end
            else                
                % 2022년 이후 
                possibleSuffixes = ["", "_1", "_2"];
                for s = 1:length(possibleSuffixes)
                    pattern = sprintf('JXR_BSC_Section_%s%s.csv', day_str, possibleSuffixes(s));
                    idx = find(strcmpi(all_filenames, pattern));
                    if ~isempty(idx)
                        files{end+1} = fullfile(folder, all_filenames{idx(1)});
                    end
                end
            end
            
        case 'TOTAL'
            if year == 2021
                pattern1 = sprintf('Total_%s_LGCHEM_BSC.csv', day_str);
                idx = find(strcmpi(all_filenames, pattern1));
                if ~isempty(idx)
                    files{end+1} = fullfile(folder, all_filenames{idx(1)});
                end
                
                for i = 1:10
                    pattern = sprintf('Total_%s_LGCHEM_BSC_%d.csv', day_str, i);
                    idx = find(strcmpi(all_filenames, pattern));
                    if ~isempty(idx)
                        files{end+1} = fullfile(folder, all_filenames{idx(1)});
                    end
                end
            else
                % 2022년 이후 
                possibleSuffixes = ["", "_1", "_2"];
                for s = 1:length(possibleSuffixes)
                    pattern = sprintf('Total_JXR_BSC_Section_%s%s.csv', day_str, possibleSuffixes(s));
                    idx = find(strcmpi(all_filenames, pattern));
                    if ~isempty(idx)
                        files{end+1} = fullfile(folder, all_filenames{idx(1)});
                    end
                end
            end
            
        case 'ONLINE'
            if year == 2021
                pattern1 = sprintf('Online_%s_LGCHEM_BSC.csv', day_str);
                idx = find(strcmpi(all_filenames, pattern1));
                if ~isempty(idx)
                    files{end+1} = fullfile(folder, all_filenames{idx(1)});
                end
                
                for i = 1:10
                    pattern = sprintf('Online_%s_LGCHEM_BSC_%d.csv', day_str, i);
                    idx = find(strcmpi(all_filenames, pattern));
                    if ~isempty(idx)
                        files{end+1} = fullfile(folder, all_filenames{idx(1)});
                    end
                end
            else
                % 2022년 이후 
                possibleSuffixes = ["", "_1", "_2"];
                for s = 1:length(possibleSuffixes)
                    pattern = sprintf('Online_JXR_BSC_Section_%s%s.csv', day_str, possibleSuffixes(s));
                    idx = find(strcmpi(all_filenames, pattern));
                    if ~isempty(idx)
                        files{end+1} = fullfile(folder, all_filenames{idx(1)});
                    end
                end
            end
            
        case 'PLC'
            if year == 2021
                pattern1 = sprintf('%s_LGCHEM_PLC#1.csv', day_str);
                idx = find(strcmpi(all_filenames, pattern1));
                if ~isempty(idx)
                    files{end+1} = fullfile(folder, all_filenames{idx(1)});
                end
                
                for i = 1:5
                    pattern = sprintf('%s_LGCHEM_PLC#1_%d.csv', day_str, i);
                    idx = find(strcmpi(all_filenames, pattern));
                    if ~isempty(idx)
                        files{end+1} = fullfile(folder, all_filenames{idx(1)});
                    end
                end
            else
                % 2022년 이후 
                possibleSuffixes = ["", "_1", "_2"];
                for s = 1:length(possibleSuffixes)
                    pattern = sprintf('JXR_BSC_PLC1_%s%s.csv', day_str, possibleSuffixes(s));
                    idx = find(strcmpi(all_filenames, pattern));
                    if ~isempty(idx)
                        files{end+1} = fullfile(folder, all_filenames{idx(1)});
                    end
                end
            end
    end
end

function result = if_not_empty(arr, default_value)
    if isempty(arr)
        result = [];
    else
        result = arr(default_value);
    end
end                