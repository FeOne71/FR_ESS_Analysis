%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rack Data DCIR Charge - New Logic (idle -> charging) - Multi Year Analysis
% Peak detection: NewLogic (idle -> charging)
% 1초 샘플링 데이터, 새로운 피크 검출 로직 적용
% 2021, 2022, 2023년도 분석

clc; clear; close all;

%% Directory
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
yearList = {'2021', '2022', '2023'};
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\Charge\DCIR_Charge_Onori_NewLogic');

if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

%% Parameters
C_nom_cell = 128;
thr = C_nom_cell * 0.1;   % Initial current threshold (A) 
dI  = C_nom_cell * 0.2;    % Current change threshold (A)
ddI = 1;                   % Continuous current increase threshold (A)

%% Visualization variables
Markersize = 15;
Fontsize = 10;
LineWidth = 3;
BarLineWidth = 1.5;

rackNames_all = {'Rack01'}; %, 'Rack02', 'Rack03', 'Rack04', 'Rack05', 'Rack06', 'Rack07', 'Rack08'};
global_eventStruct = struct();
for r = 1:length(rackNames_all)
    global_eventStruct.(rackNames_all{r}) = struct();
end

% Visualization variables
all_peak_times = {};
all_peak_currents = {};
all_peak_voltages = {};
all_raw_times = {};
all_raw_currents = {};
all_raw_voltages = {};
all_original_data = {};  % 원본 데이터 저장용
all_peak_counts = [];
all_filenames = {};
all_dcir_values = [];  % DCIR 값 수집용

for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    fprintf('Processing year: %s\n', year);
    eventStruct = struct();
    eventCount = 0;
         yearPath = fullfile(dataDir, year);
     fprintf('Year path: %s\n', yearPath);
     fprintf('Year path exists: %d\n', exist(yearPath, 'dir'));
     monthDirs = dir(fullfile(yearPath, '20*'));
     fprintf('Found %d month directories\n', length(monthDirs));
     for m = 1:length(monthDirs)
             if ~monthDirs(m).isdir, continue; end
             monthPath = fullfile(yearPath, monthDirs(m).name);
             fprintf('Processing month: %s\n', monthDirs(m).name);
             matFiles = dir(fullfile(monthPath, '*.mat'));
             fprintf('Found %d mat files in %s\n', length(matFiles), monthDirs(m).name);
             for f = 1:length(matFiles)
                         matFilePath = fullfile(monthPath, matFiles(f).name);
                           fprintf('Loading file: %s\n', matFiles(f).name);
              load(matFilePath);
              
              % Check if Raw_Rack variable exists
              if ~exist('Raw_Rack', 'var')
                  fprintf('Warning: Raw_Rack variable not found in file %s\n', matFiles(f).name);
                  continue;
              end
              
              fprintf('Raw_Rack structure fields: %s\n', strjoin(fieldnames(Raw_Rack), ', '));
              rackNames = rackNames_all;
            
            % Initialize daily peak tracking for this file
            daily_peak_times = {};
            daily_peak_currents = {};
            daily_peak_voltages = {};
            daily_raw_times = {};
            daily_raw_currents = {};
            daily_raw_voltages = {};
            daily_original_data = {};  % 원본 데이터 저장용
            daily_total_peaks = 0;
            
            for rack_idx = 1:length(rackNames)
                rackName = rackNames{rack_idx};
                if ~isfield(Raw_Rack, rackName), continue; end
                rackData = Raw_Rack.(rackName);
                t = rackData.Time;
                I = rackData.DCCurrent_A;
                V = rackData.AverageCV_V;
                soc = rackData.SOCPct;
                T_batt = rackData.AverageMT_degC;
                dc_power = rackData.DCPower_kW;
                
                % Convert time data to numeric
                if iscell(t) || isstring(t)
                    t = datetime(t);
                end
                if isdatetime(t)
                    t = seconds(t - t(1));
                end
                N = length(I);
                
                % Derivative of Current (원본 전류의 미분값)
                for i = 1:length(I)-1
                    dI_dt(i) = (I(i+1) - I(i)) / (t(i+1) - t(i));
                end
                
                % New Peak Detection Logic (idle -> charging)
                PeakTime = {};
                PeakCurrent = {};
                PeakVoltage = {};
                z = 1;
                prev_peak_end_idx = 0;  % 이전 피크 끝점 인덱스
                
                for i = 1:(length(I)-1)
                    % 1. idle 구간에서 charging으로 전환되는 지점 찾기
                    if I(i) > -thr && I(i) < thr  % idle 구간
                        if i+1 <= length(I) && I(i+1) >= thr  % 다음 시점이 charging
                            % 2. 피크 시작점 확인 (idle의 마지막 시점)
                            peak_start_idx = i;
                            
                            % 3. 피크 끝점 찾기 (전류가 더 이상 증가하지 않는 첫 시점)
                            peak_end_idx = peak_start_idx;
                            for j = peak_start_idx+1:length(I)-1
                                if I(j+1) <= I(j)  % 전류 증가가 멈춘 시점
                                    peak_end_idx = j;
                                    break;                                
                                end
                            end
                            
                            % 4. 피크 길이 확인 (최소 2개 포인트 이상)
                            if peak_end_idx > peak_start_idx
                                % 5. 중복 피크 방지 (겹침 확인)
                                if peak_start_idx > prev_peak_end_idx
                                    % 6. 피크 검출 조건 확인 (NotFiltered와 동일)
                                    if (I(peak_end_idx) - I(peak_start_idx)) >= dI
                                        if (I(peak_start_idx+1) - I(peak_start_idx)) > ddI
                                            flag = 1;
                                            % 구간 내 전류/미분 체크
                                            for zi = peak_start_idx:peak_end_idx-1
                                                if zi <= length(dI_dt) && (dI_dt(zi) < 0 || I(zi+1) < 0)
                                                    flag = 0;
                                                    break;
                                                end
                                            end
                                            
                                            if flag == 1
                                                % 피크 저장
                                                peak_length = peak_end_idx - peak_start_idx + 1;
                                                PeakTime{z} = t(peak_start_idx:peak_end_idx);
                                                PeakCurrent{z} = I(peak_start_idx:peak_end_idx);
                                                PeakVoltage{z} = V(peak_start_idx:peak_end_idx);
                                                z = z + 1;
                                                prev_peak_end_idx = peak_end_idx;  % 끝점 업데이트
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Accumulate daily peaks from all racks
                if ~isempty(PeakTime)
                    % 각 랙별로 cell 배열로 저장
                    if isempty(daily_peak_times)
                        daily_peak_times = {PeakTime};
                        daily_peak_currents = {PeakCurrent};
                        daily_peak_voltages = {PeakVoltage};
                        daily_raw_times = {t};
                        daily_raw_currents = {I};
                        daily_raw_voltages = {V};
                    else
                        daily_peak_times{end+1} = PeakTime;
                        daily_peak_currents{end+1} = PeakCurrent;
                        daily_peak_voltages{end+1} = PeakVoltage;
                        daily_raw_times{end+1} = t;
                        daily_raw_currents{end+1} = I;
                        daily_raw_voltages{end+1} = V;
                    end
                    daily_total_peaks = daily_total_peaks + length(PeakTime);
                    fprintf('Rack %s: %d peaks detected\n', rackName, length(PeakTime));
                else
                    fprintf('Rack %s: No peaks detected\n', rackName);
                end
                
                % 원본 데이터 저장 (시간, 전압, 전류, 온도, SOC)
                if isempty(daily_original_data)
                    daily_original_data = {struct('time', rackData.Time, 'voltage', rackData.AverageCV_V, ...
                        'current', rackData.DCCurrent_A, 'temperature', rackData.AverageMT_degC, 'soc', rackData.SOCPct)};
                else
                    daily_original_data{end+1} = struct('time', rackData.Time, 'voltage', rackData.AverageCV_V, ...
                        'current', rackData.DCCurrent_A, 'temperature', rackData.AverageMT_degC, 'soc', rackData.SOCPct);
                end
                
                % Resistance Computation
                for i = 1:length(PeakTime)
                    DV = (PeakVoltage{i}(end) - PeakVoltage{i}(1));
                    DI = PeakCurrent{i}(end) - PeakCurrent{i}(1);
                    
                    if DI > 0 && PeakCurrent{i}(end) > 0 && DV > 0
                        dcir_val = (DV / DI) * 1000;
                        all_dcir_values = [all_dcir_values; dcir_val];  % DCIR 값 수집
                    else
                        dcir_val = NaN;
                    end
                    
                    eventCount = eventCount + 1;
                    evtName = sprintf('event%d', eventCount);
                    
                    % Extract date from filename
                    date_str = matFiles(f).name(5:12);  % Extract YYYYMMDD from Raw_YYYYMMDD.mat
                    formatted_date = sprintf('%s-%s-%s', date_str(1:4), date_str(5:6), date_str(7:8));
                    
                    eventStruct.(evtName).rack_name = rackName;
                    eventStruct.(evtName).year = year;
                    eventStruct.(evtName).date = formatted_date;
                    eventStruct.(evtName).start_idx = find(t >= PeakTime{i}(1), 1);
                    eventStruct.(evtName).end_idx = find(t >= PeakTime{i}(end), 1);
                    eventStruct.(evtName).charge_duration = length(PeakTime{i});
                    eventStruct.(evtName).avg_current = mean(PeakCurrent{i});
                    eventStruct.(evtName).t_seq = rackData.Time(eventStruct.(evtName).start_idx:eventStruct.(evtName).end_idx);
                    eventStruct.(evtName).I_seq = PeakCurrent{i};
                    eventStruct.(evtName).V_seq = PeakVoltage{i};
                    eventStruct.(evtName).PeakChgR = dcir_val;
                    eventStruct.(evtName).PeakChgR_DV = DV;
                    eventStruct.(evtName).PeakChgR_DI = DI;
                    eventStruct.(evtName).PeakChgR_V1 = PeakVoltage{i}(1);
                    eventStruct.(evtName).PeakChgR_V2 = PeakVoltage{i}(end);
                    eventStruct.(evtName).PeakChgR_I1 = PeakCurrent{i}(1);
                    eventStruct.(evtName).PeakChgR_I2 = PeakCurrent{i}(end);
                end
            end
            
            % Store daily data if any peaks were found
            if daily_total_peaks > 0
                all_peak_times{end+1} = daily_peak_times;
                all_peak_currents{end+1} = daily_peak_currents;
                all_peak_voltages{end+1} = daily_peak_voltages;
                all_raw_times{end+1} = daily_raw_times;
                all_raw_currents{end+1} = daily_raw_currents;
                all_raw_voltages{end+1} = daily_raw_voltages;
                all_original_data{end+1} = daily_original_data;  % 원본 데이터 저장
                all_peak_counts(end+1) = daily_total_peaks;
                all_filenames{end+1} = matFiles(f).name;
                fprintf('Date %s: Total %d peaks detected\n', matFiles(f).name, daily_total_peaks);
            end
        end
    end
    
    % Save structure (by rack/year/date)
    for rack_idx = 1:length(rackNames_all)
        rackName = rackNames_all{rack_idx};
        eventNames = fieldnames(eventStruct);
        rack_events = {};
        for i = 1:length(eventNames)
            evt = eventStruct.(eventNames{i});
            if strcmp(evt.rack_name, rackName)
                rack_events{end+1} = eventNames{i};
            end
        end
        if isempty(rack_events), continue; end
        
        % Initialize rack structure if not exists
        if ~isfield(global_eventStruct, rackName)
            global_eventStruct.(rackName) = struct();
        end
        
        % Initialize year structure if not exists
        year_str = sprintf('year_%s', year);
        if ~isfield(global_eventStruct.(rackName), year_str)
            global_eventStruct.(rackName).(year_str) = struct();
        end
        
        % Group events by date
        date_events = struct();
        for i = 1:length(rack_events)
            evt = eventStruct.(rack_events{i});
            date_key = evt.date;
            
            % Convert date to valid field name (add E prefix and replace - with _)
            valid_date_key = ['E' strrep(date_key, '-', '_')];
            
            if ~isfield(date_events, valid_date_key)
                date_events.(valid_date_key) = struct();
            end
            
            date_events.(valid_date_key).(rack_events{i}) = evt;
        end
        
        % Save events by date
        date_names = fieldnames(date_events);
        for d = 1:length(date_names)
            date_key = date_names{d};
            global_eventStruct.(rackName).(year_str).(date_key) = date_events.(date_key);
        end
    end
end


%% DCIR Histogram
if ~isempty(all_dcir_values)
    % 유효한 DCIR 값만 필터링 (NaN 제거)
    valid_dcir = all_dcir_values(~isnan(all_dcir_values));
    
    if ~isempty(valid_dcir)
        % 1σ 이상치 제거
        mean_val = mean(valid_dcir);
        std_val = std(valid_dcir);
        lower_bound = mean_val - std_val;
        upper_bound = mean_val + std_val;
        
        % 이상치 제거 전 데이터 수
        original_count = length(valid_dcir);
        
        % 1σ 범위 내 데이터만 유지
        valid_dcir_filtered = valid_dcir(valid_dcir >= lower_bound & valid_dcir <= upper_bound);
        
        % 이상치 제거 후 데이터 수
        filtered_count = length(valid_dcir_filtered);
        removed_count = original_count - filtered_count;
        
        fprintf('Original %d values, Removed %d outliers (%.1f%%), Final %d values\n', ...
            original_count, removed_count, (removed_count/original_count)*100, filtered_count);
        
        valid_dcir = valid_dcir_filtered; % Use filtered data
        
        figure(2);
        hold on;
        grid on;
        box on;
        
        % 히스토그램 생성 (20개 bin)
        h = histogram(valid_dcir, 20, 'FaceColor', '#0073C2', 'EdgeColor', 'black', 'LineWidth', 1);
        
        % 통계 정보 계산
        mean_dcir = mean(valid_dcir);
        std_dcir = std(valid_dcir);
        median_dcir = median(valid_dcir);
        
        % 정규분포 곡선 추가
        x_curve = linspace(0, 1, 1000);
        y_curve = normpdf(x_curve, mean_dcir, std_dcir);
        scale_factor = max(h.Values) / max(y_curve);
        y_curve_scaled = y_curve * scale_factor;
        plot(x_curve, y_curve_scaled, 'r-', 'LineWidth', 3, ...
            'DisplayName', sprintf('Normal (μ=%.3f, σ=%.3f)', mean_dcir, std_dcir));
        
        % 평균선 추가
        xline(mean_dcir, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Mean: %.3f mΩ', mean_dcir));
        
        xlabel('Resistance [mΩ]', 'interpreter', 'tex');
        xlim([0 1]);
        ylabel('Frequency', 'interpreter', 'tex');
        title('Peak_{Chg} R Distribution 2021-2023 (1σ Outliers Removed)', 'FontSize', 14);
        legend('Location', 'best');
        set(gca, 'fontsize', Fontsize);
        set(gca, 'ticklabelinterpreter', 'tex');
        
        % 통계 정보 텍스트 박스 추가
        stats_text = sprintf(['Statistics (1σ outliers removed):\n' ...
            'Total Events: %d\n' ...
            'Mean: %.3f mΩ\n' ...
            'Std: %.3f mΩ\n' ...
            'Median: %.3f mΩ\n' ...
            'Min: %.3f mΩ\n' ...
            'Max: %.3f mΩ'], ...
            length(valid_dcir), mean_dcir, std_dcir, median_dcir, ...
            min(valid_dcir), max(valid_dcir));
        
        text(0.02, 0.98, stats_text, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 10, ...
            'BackgroundColor', 'white', 'EdgeColor', 'black');
        
        % 히스토그램 저장
        hist_filename = fullfile(saveDir, 'DCIR_Histogram_NewLogic_2021_2022_2023.fig');
        saveas(gcf, hist_filename);
        fprintf('DCIR Histogram saved to: %s\n', hist_filename);
        
        fprintf('DCIR Histogram complete: %d valid DCIR values\n', length(valid_dcir));
    else
        fprintf('No valid DCIR values found for histogram.\n');
    end
else
    fprintf('No DCIR values collected for histogram.\n');
end

%% Visualization
if ~isempty(all_peak_times)
    % rack01의 피크가 가장 많은 날짜 찾기
    rack01_peak_counts = zeros(length(all_peak_times), 1);
    for i = 1:length(all_peak_times)
        if ~isempty(all_peak_times{i}) && length(all_peak_times{i}) >= 1
            rack01_peak_counts(i) = length(all_peak_times{i}{1});  % rack01은 첫 번째 랙
        end
    end
    
    [~, max_idx] = max(rack01_peak_counts);
    max_peaks_filename = all_filenames{max_idx};
    
    % 최다 검출 일자의 데이터 (rack01)
    t = all_raw_times{max_idx}{1};  % rack01의 정규화된 시간 데이터
    I = all_raw_currents{max_idx}{1};  % rack01의 전류 데이터
    V = all_raw_voltages{max_idx}{1};  % rack01의 전압 데이터
    
    % 저장된 원본 데이터 사용
    original_data = all_original_data{max_idx}{1};  % rack01의 원본 데이터
    original_time = datetime(original_data.time);  % string을 datetime으로 변환
    original_V = original_data.voltage;  % 원본 전압 데이터
    original_I = original_data.current;  % 원본 전류 데이터
    
    % Extract date from filename (Raw_YYYYMMDD.mat)
    date_str = max_peaks_filename(5:12);  % Extract YYYYMMDD
    formatted_date = sprintf('%s-%s-%s', date_str(1:4), date_str(5:6), date_str(7:8));
    
    % 해당 일자의 피크 데이터 (rack01)
    PeakTime = all_peak_times{max_idx}{1};
    PeakCurrent = all_peak_currents{max_idx}{1};
    PeakVoltage = all_peak_voltages{max_idx}{1};
    
    figure(1);
    
    % Voltage subplot
    subplot(2,1,1);
    hold on; 
    grid on;
    box on;
    plot(original_time, original_V, 'Color', '#C0C0C0', 'Linewidth', LineWidth-1, 'HandleVisibility', 'off');
    
    % 피크 구간 오버랩 (t_seq 사용)
    eventNames = fieldnames(eventStruct);
    for i = 1:length(eventNames)
        evt = eventStruct.(eventNames{i});
        if strcmp(evt.rack_name, 'Rack01') && strcmp(evt.date, formatted_date)
            % 디버깅: 피크 데이터 확인
            fprintf('Debug: evt.t_seq length = %d, evt.V_seq length = %d\n', length(evt.t_seq), length(evt.V_seq));
            if ~isempty(evt.t_seq) && ~isempty(evt.V_seq)
                plot(datetime(evt.t_seq), evt.V_seq, 'r-', 'Linewidth', LineWidth+1, 'Color', '#0073C2');
            end
        end
    end
    ylabel('Voltage [V]', 'interpreter', 'tex');
    set(gca, 'fontsize', Fontsize);
    set(gca, 'ticklabelinterpreter', 'tex');
    title('Rack01 Voltage Data with Detected Peaks');
    
    % Current subplot
    subplot(2,1,2);
    hold on;
    grid on;
    box on;
    plot(original_time, original_I, 'Color', '#C0C0C0', 'Linewidth', LineWidth-1);
    
    % 피크 구간 오버랩 (t_seq 사용)
    for i = 1:length(eventNames)
        evt = eventStruct.(eventNames{i});
        if strcmp(evt.rack_name, 'Rack01') && strcmp(evt.date, formatted_date)
            % 디버깅: 피크 데이터 확인
            fprintf('Debug: evt.t_seq length = %d, evt.I_seq length = %d\n', length(evt.t_seq), length(evt.I_seq));
            if ~isempty(evt.t_seq) && ~isempty(evt.I_seq)
                plot(datetime(evt.t_seq), evt.I_seq, 'r-', 'Linewidth', LineWidth+1, 'Color', '#0073C2');
            end
        end
    end
    xlabel('Time', 'interpreter', 'tex');
    ylabel('Current [A]', 'interpreter', 'tex');
    set(gca, 'fontsize', Fontsize);
    set(gca, 'ticklabelinterpreter', 'tex');
    title('Rack01 Current Data with Detected Peaks');
    
    % Calculate rack01 peaks for the selected day
    rack01_peaks_for_day = length(PeakTime);
    
    sgtitle(sprintf('Date: %s, Peaks: %d 2021-2023 Data', formatted_date, rack01_peaks_for_day), 'FontSize', 16);
    
    fig_filename = fullfile(saveDir, 'Onori_Method_NewLogic_Visualization_2021_2022_2023.fig');
    saveas(gcf, fig_filename);
    fprintf('Figure saved to: %s\n', fig_filename);
    
    fprintf('Visualization complete: %d peaks detected (New Logic).\n', length(PeakTime));
else
    fprintf('No peaks detected.\n');
end

save(fullfile(saveDir, 'all_chg_events_onori_newlogic_all_years_2021_2022_2023.mat'), 'global_eventStruct');
fprintf('Processing complete\n');
fprintf('Results saved to: %s\n', saveDir);
