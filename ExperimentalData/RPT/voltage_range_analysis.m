%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3년치 Raw MAT 파일 전압 범위 분석
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% 파일 경로 설정
folderPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\raw2mat_ver04';
savePath = fullfile(pwd, 'Voltage_Analysis_Results');

if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%% MAT 파일 찾기 (연도별 폴더 구조)
all_mat_files = [];
year_list = {'2021', '2022', '2023'};

for year_idx = 1:length(year_list)
    current_year = year_list{year_idx};
    year_path = fullfile(folderPath, current_year);
    
    if exist(year_path, 'dir')
        months = dir(year_path);
        months = months([months.isdir]);
        months = months(~ismember({months.name}, {'.', '..'}));
        
        for month_idx = 1:length(months)
            current_month = months(month_idx).name;
            month_path = fullfile(year_path, current_month);
            
            % 해당 월의 Raw_YYYYMMDD.mat 파일들 찾기
            mat_files = dir(fullfile(month_path, 'Raw_*.mat'));
            
            for file_idx = 1:length(mat_files)
                file_info = mat_files(file_idx);
                file_info.fullpath = fullfile(month_path, file_info.name);
                file_info.year = current_year;
                file_info.month = current_month;
                file_info.date = file_info.name(5:12); % Raw_YYYYMMDD.mat에서 날짜 추출
                all_mat_files = [all_mat_files; file_info];
            end
        end
    end
end

fprintf('Found %d MAT files\n', length(all_mat_files));

%% 전압 및 SOC 데이터 수집
all_voltage_data = [];
all_soc_data = [];
all_date_data = [];

for file_idx = 1:length(all_mat_files)
    file_info = all_mat_files(file_idx);
    filepath = file_info.fullpath;
    
    fprintf('Processing: %s (%s-%s)\n', file_info.name, file_info.year, file_info.month);
    
    try
        % MAT 파일 로드
        data = load(filepath);
        
        % Raw 구조체에서 전압 및 SOC 데이터 추출
        if isfield(data, 'Raw')
            raw_data = data.Raw;
            
            % Total_AverageCVSum_V (전압) 데이터 추출
            if isfield(raw_data, 'Total_AverageCVSum_V')
                voltage_data = raw_data.Total_AverageCVSum_V;
                
                % Total_AverageSOC 데이터 추출
                if isfield(raw_data, 'Total_AverageSOC')
                    soc_data = raw_data.Total_AverageSOC;
                    
                    % 데이터 수집
                    all_voltage_data = [all_voltage_data; voltage_data(:)];
                    all_soc_data = [all_soc_data; soc_data(:)];
                    all_date_data = [all_date_data; repmat(datetime(file_info.date, 'InputFormat', 'yyyyMMdd'), length(voltage_data), 1)];
                    
                    fprintf('Extracted %d points: V=%.3f-%.3f V, SOC=%.1f-%.1f%%\n', ...
                           length(voltage_data), min(voltage_data), max(voltage_data), ...
                           min(soc_data), max(soc_data));
                else
                    fprintf('No Total_AverageSOC field found in %s\n', file_info.name);
                end
            else
                fprintf('No Total_AverageCVSum_V field found in %s\n', file_info.name);
            end
        else
            fprintf('No Raw structure found in %s\n', file_info.name);
        end
        
    catch ME
        fprintf('Error processing %s: %s\n', file_info.name, ME.message);
    end
end

%% 일자별 전압/SOC 범위 및 월별 평균 분석
if ~isempty(all_voltage_data)
    fprintf('\n=== Daily Voltage/SOC Range and Monthly Average Analysis ===\n');
    fprintf('Total data points: %d\n', length(all_voltage_data));
    fprintf('Date range: %s to %s\n', datestr(min(all_date_data)), datestr(max(all_date_data)));
    
    % 일자별 전압/SOC 범위 계산
    unique_dates = unique(all_date_data);
    daily_voltage_ranges = zeros(length(unique_dates), 2); % min, max
    daily_soc_ranges = zeros(length(unique_dates), 2); % min, max
    
    for date_idx = 1:length(unique_dates)
        date_val = unique_dates(date_idx);
        date_mask = all_date_data == date_val;
        date_voltage = all_voltage_data(date_mask);
        date_soc = all_soc_data(date_mask);
        
        daily_voltage_ranges(date_idx, 1) = min(date_voltage);
        daily_voltage_ranges(date_idx, 2) = max(date_voltage);
        daily_soc_ranges(date_idx, 1) = min(date_soc);
        daily_soc_ranges(date_idx, 2) = max(date_soc);
    end
    
    % 월별 평균 범위 계산
    years_data = year(all_date_data);
    months_data = month(all_date_data);
    unique_years_months = unique([years_data, months_data], 'rows');
    monthly_voltage_ranges = zeros(size(unique_years_months, 1), 2); % min, max
    monthly_soc_ranges = zeros(size(unique_years_months, 1), 2); % min, max
    
    for m_idx = 1:size(unique_years_months, 1)
        year_val = unique_years_months(m_idx, 1);
        month_val = unique_years_months(m_idx, 2);
        
        % 해당 월의 일자별 범위들 찾기
        month_daily_voltage_ranges = [];
        month_daily_soc_ranges = [];
        
        for date_idx = 1:length(unique_dates)
            date_val = unique_dates(date_idx);
            if year(date_val) == year_val && month(date_val) == month_val
                month_daily_voltage_ranges = [month_daily_voltage_ranges; daily_voltage_ranges(date_idx, :)];
                month_daily_soc_ranges = [month_daily_soc_ranges; daily_soc_ranges(date_idx, :)];
            end
        end
        
        % 월별 평균 범위 (일자별 범위들의 평균)
        if ~isempty(month_daily_voltage_ranges)
            monthly_voltage_ranges(m_idx, 1) = mean(month_daily_voltage_ranges(:, 1)); % 평균 최소값
            monthly_voltage_ranges(m_idx, 2) = mean(month_daily_voltage_ranges(:, 2)); % 평균 최대값
        end
        if ~isempty(month_daily_soc_ranges)
            monthly_soc_ranges(m_idx, 1) = mean(month_daily_soc_ranges(:, 1)); % 평균 최소값
            monthly_soc_ranges(m_idx, 2) = mean(month_daily_soc_ranges(:, 2)); % 평균 최대값
        end
    end
    
    % 연도별 평균 범위 계산
    unique_years = unique(years_data);
    year_voltage_ranges = zeros(length(unique_years), 2); % min, max
    year_soc_ranges = zeros(length(unique_years), 2); % min, max
    
    for year_idx = 1:length(unique_years)
        year_val = unique_years(year_idx);
        
        % 해당 연도의 월별 범위들
        year_month_mask = unique_years_months(:, 1) == year_val;
        year_monthly_voltage_ranges = monthly_voltage_ranges(year_month_mask, :);
        year_monthly_soc_ranges = monthly_soc_ranges(year_month_mask, :);
        
        % 연도별 평균 범위 (월별 범위들의 평균)
        year_voltage_ranges(year_idx, 1) = mean(year_monthly_voltage_ranges(:, 1)); % 평균 최소값
        year_voltage_ranges(year_idx, 2) = mean(year_monthly_voltage_ranges(:, 2)); % 평균 최대값
        year_soc_ranges(year_idx, 1) = mean(year_monthly_soc_ranges(:, 1)); % 평균 최소값
        year_soc_ranges(year_idx, 2) = mean(year_monthly_soc_ranges(:, 2)); % 평균 최대값
    end
    
    % 연도별 범위 요약 그래프
    figure('Name', 'Year-wise Voltage and SOC Ranges', 'Position', [100 100 1200 800]);
    
    % Subplot 1: 연도별 전압 범위
    subplot(1,2,1);
    hold on;
    for year_idx = 1:length(unique_years)
        year_val = unique_years(year_idx);
        min_v = year_voltage_ranges(year_idx, 1);
        max_v = year_voltage_ranges(year_idx, 2);
        
        % 전압 범위를 선으로 표시
        plot([year_val, year_val], [min_v, max_v], 'b-o', 'LineWidth', 4, 'MarkerSize', 10);
        
        % 범위 값을 텍스트로 표시
        text(year_val, max_v + (max_v - min_v) * 0.1, sprintf('%.0f-%.0fV', min_v, max_v), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
    xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Voltage Range [V]', 'FontSize', 12, 'FontWeight', 'bold');
    title('Voltage Range by Year', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    xticks(unique_years); % 연도만 표시
    
    % Subplot 2: 연도별 SOC 범위
    subplot(1,2,2);
    hold on;
    for year_idx = 1:length(unique_years)
        year_val = unique_years(year_idx);
        min_soc = year_soc_ranges(year_idx, 1);
        max_soc = year_soc_ranges(year_idx, 2);
        
        % SOC 범위를 선으로 표시
        plot([year_val, year_val], [min_soc, max_soc], 'r-s', 'LineWidth', 4, 'MarkerSize', 10);
        
        % 범위 값을 텍스트로 표시
        text(year_val, max_soc + (max_soc - min_soc) * 0.1, sprintf('%.1f-%.1f%%', min_soc, max_soc), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
    xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('SOC Range [%]', 'FontSize', 12, 'FontWeight', 'bold');
    title('SOC Range by Year', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    xticks(unique_years); % 연도만 표시
    
    % 저장
    figName = fullfile(savePath, 'year_voltage_soc_ranges.fig');
    savefig(gcf, figName);
    fprintf('Saved: %s\n', figName);
    
    % 데이터 저장
    voltage_analysis = struct();
    voltage_analysis.all_voltage_data = all_voltage_data;
    voltage_analysis.all_soc_data = all_soc_data;
    voltage_analysis.all_date_data = all_date_data;
    voltage_analysis.daily_voltage_ranges = daily_voltage_ranges;
    voltage_analysis.daily_soc_ranges = daily_soc_ranges;
    voltage_analysis.monthly_voltage_ranges = monthly_voltage_ranges;
    voltage_analysis.monthly_soc_ranges = monthly_soc_ranges;
    voltage_analysis.year_voltage_ranges = year_voltage_ranges;
    voltage_analysis.year_soc_ranges = year_soc_ranges;
    
    matName = fullfile(savePath, 'voltage_soc_analysis.mat');
    save(matName, 'voltage_analysis');
    fprintf('Saved: %s\n', matName);
    
else
    fprintf('No voltage data found!\n');
end

fprintf('\n=== Voltage Analysis Completed ===\n'); 