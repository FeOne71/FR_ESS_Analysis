%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find missing dates in CSV files
% 폴더에서 CSV 파일명의 날짜를 추출하여 누락된 날짜를 찾는 스크립트
% 파일명 형식: YYYYMMDD_KIMJ_LGE_01_01_01.csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 경로 설정
data_dir = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\김제 배터리 데이터';

fprintf('Checking CSV files in: %s\n', data_dir);
fprintf('%s\n', repmat('=', 1, 80));

%% CSV 파일 목록 읽기
csv_files = dir(fullfile(data_dir, '*.csv'));
fprintf('Total CSV files found: %d\n', length(csv_files));

if isempty(csv_files)
    error('No CSV files found in the directory');
end

%% 파일명에서 날짜 추출
file_dates = [];
file_names = {csv_files.name};

for i = 1:length(file_names)
    filename = file_names{i};
    if length(filename) >= 8
        date_str = filename(1:8); % YYYYMMDD
        try
            date_num = datetime(date_str, 'InputFormat', 'yyyyMMdd');
            file_dates = [file_dates; date_num];
        catch
            fprintf('Warning: Could not parse date from file: %s\n', filename);
        end
    end
end

if isempty(file_dates)
    error('No valid dates found in file names');
end

% 중복 제거 및 정렬
file_dates = unique(file_dates);
file_dates = sort(file_dates);

fprintf('Unique dates found: %d\n', length(file_dates));
fprintf('Date range: %s to %s\n', datestr(min(file_dates)), datestr(max(file_dates)));

%% 연도/월별로 그룹화
years = year(file_dates);
months = month(file_dates);
year_month_pairs = unique([years, months], 'rows');

%% 누락된 날짜 찾기
missing_dates_all = [];

for ym_idx = 1:size(year_month_pairs, 1)
    y = year_month_pairs(ym_idx, 1);
    m = year_month_pairs(ym_idx, 2);
    
    % 해당 월의 모든 날짜 생성
    last_day = eomday(y, m);
    all_dates_in_month = datetime(y, m, 1:last_day)';
    
    % 해당 월의 실제 파일 날짜
    mask = (years == y) & (months == m);
    actual_dates = file_dates(mask);
    
    % 누락된 날짜 찾기
    missing_dates = setdiff(all_dates_in_month, actual_dates);
    
    if ~isempty(missing_dates)
        for d = 1:length(missing_dates)
            missing_dates_all = [missing_dates_all; struct('Year', y, 'Month', m, 'Date', missing_dates(d))];
        end
    end
end

%% 결과 출력
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('Missing Dates Summary\n');
fprintf('%s\n', repmat('=', 1, 80));

if isempty(missing_dates_all)
    fprintf('No missing dates found!\n');
else
    fprintf('Total missing dates: %d\n\n', length(missing_dates_all));
    
    % 표 형식으로 출력 (탭 구분 - Excel에 복사 가능)
    % 각 날짜를 개별 행으로 출력
    fprintf('Year\tMonth\tDate\n');
    fprintf('%s\n', repmat('-', 1, 80));
    
    % 연도/월/날짜 순으로 정렬
    [~, sort_idx] = sort([missing_dates_all.Date]);
    missing_dates_sorted = missing_dates_all(sort_idx);
    
    for i = 1:length(missing_dates_sorted)
        date_str = datestr(missing_dates_sorted(i).Date, 'yyyymmdd');
        fprintf('%d\t%d\t%s\n', missing_dates_sorted(i).Year, missing_dates_sorted(i).Month, date_str);
    end
    
    fprintf('%s\n', repmat('-', 1, 80));
    
    % 연도/월별 요약 (참고용)
    fprintf('\n%s\n', repmat('=', 1, 80));
    fprintf('Summary by Year/Month (for reference)\n');
    fprintf('%s\n', repmat('=', 1, 80));
    
    % 연도/월별로 누락된 날짜 그룹화
    missing_by_ym = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:length(missing_dates_all)
        key = sprintf('%d_%02d', missing_dates_all(i).Year, missing_dates_all(i).Month);
        if isKey(missing_by_ym, key)
            missing_by_ym(key) = [missing_by_ym(key); missing_dates_all(i).Date];
        else
            missing_by_ym(key) = missing_dates_all(i).Date;
        end
    end
    
    % 키를 정렬
    keys = keys(missing_by_ym);
    keys_sorted = sort(keys);
    
    fprintf('Year\tMonth\tCount\n');
    fprintf('%s\n', repmat('-', 1, 80));
    
    for k = 1:length(keys_sorted)
        key = keys_sorted{k};
        parts = strsplit(key, '_');
        y = str2double(parts{1});
        m = str2double(parts{2});
        dates = missing_by_ym(key);
        fprintf('%d\t%d\t%d\n', y, m, length(dates));
    end
    
    fprintf('%s\n', repmat('-', 1, 80));
end

%% 통계 정보
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('Statistics\n');
fprintf('%s\n', repmat('=', 1, 80));

% 전체 기간 계산
start_date = min(file_dates);
end_date = max(file_dates);
total_days = days(end_date - start_date) + 1;
expected_files = total_days;
actual_files = length(file_dates);
missing_count = length(missing_dates_all);

fprintf('Date range: %s to %s\n', datestr(start_date), datestr(end_date));
fprintf('Total days in range: %d\n', total_days);
fprintf('Files found: %d\n', actual_files);
fprintf('Missing dates: %d\n', missing_count);
fprintf('Coverage: %.2f%%\n', (actual_files / total_days) * 100);

fprintf('\nScript completed!\n');
