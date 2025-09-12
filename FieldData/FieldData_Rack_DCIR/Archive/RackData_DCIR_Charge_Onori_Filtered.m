%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rack Data DCIR Charge - Onori Sliding Window Method (fig_3.m 스타일)
% Peak detection: Onori Sliding Window Method (fig_3.m 스타일)
% 1초 샘플링 데이터, 30초 슬라이딩 윈도우, 피크 감지 및 DCIR 계산
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir  = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Rack_raw2mat';
yearList = {'2021'}; %, '2022', '2023'}; 
saveDir  = fullfile('G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_Rack_DCIR\DCIR_Charge_Onori');
if ~exist(saveDir, 'dir')
    mkdir(saveDir); 
end

%% Parameters
C_nom_cell = 128;
dt  = 4;                   % Current monotonic increase interval
thr = C_nom_cell * 0.02;   % Initial current threshold (A)
dI  = C_nom_cell * 0.2;    % Current change threshold after 3s (A)
ddI = 1;                   % Continuous current increase threshold (A)
te  = 4;                   % Filter window size

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

for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    fprintf('Processing year: %s\n', year);
    eventStruct = struct();
    eventCount = 0;
    yearPath = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));
    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, 'Raw_*.mat'));
        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            load(matFilePath); 
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
                
                % Filter on Current
                I_filt = zeros(N,1);
                Pre = 0*ones(te/2,1);
                Post = 0*ones(te/2,1);
                I_calc = [Pre; I; Post];
                for i = 1:N
                    for m = 0:te
                        I_filt(i) = I_filt(i) + I_calc(i+m);
                    end
                    I_filt(i) = I_filt(i)/(te+1);
                end
                
                % Derivative of Current
                for i = 1:length(I)-1
                    dI_dt(i) = (I(i+1) - I(i)) / (t(i+1) - t(i));
                end
                
                % MA on Current Derivative
                N_dI = length(dI_dt);
                filt_dI_dt = zeros(N_dI,1);
                Pre = 0*ones(te/2,1);
                Post = 0*ones(te/2,1);
                calc_dI_dt = [Pre; dI_dt'; Post];
                for i = 1:N_dI
                    for m = 0:te
                        filt_dI_dt(i) = filt_dI_dt(i) + calc_dI_dt(i+m);
                    end
                    filt_dI_dt(i) = filt_dI_dt(i)/(te+1);
                end
                
                % Driving Peaks Identification
                PeakTime = {};
                PeakCurrent = {};
                PeakVoltage = {};
                z = 1;
                
                for i = 1:(length(I_filt)-dt)
                    if (I_filt(i+dt) - I_filt(i)) > dI
                        if I(i) > -thr && I(i) < thr
                            if (I_filt(i+1) - I_filt(i)) > ddI
                                flag = 1;
                                for zi = 1:dt
                                    if filt_dI_dt(i+zi-1) < 0 || I(i+zi) < 0
                                        flag = 0;
                                    end
                                end
                                if flag == 1
                                    if z == 1
                                        for j = 1:dt
                                            PeakTime{z}(j) = t(i+j-1);
                                            PeakCurrent{z}(j) = I(i+j-1);
                                            PeakVoltage{z}(j) = V(i+j-1);
                                        end
                                        z = z + 1;
                                    else
                                        if (t(i) - PeakTime{z-1}(1)) > dt
                                            for j = 1:dt
                                                PeakTime{z}(j) = t(i+j-1);
                                                PeakCurrent{z}(j) = I(i+j-1);
                                                PeakVoltage{z}(j) = V(i+j-1);
                                            end
                                            z = z + 1;
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
                    
                    if DI > 0 && PeakCurrent{i}(end) > 0
                        dcir_val = (DV / DI) * 1000;
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
                    eventStruct.(evtName).end_idx = eventStruct.(evtName).start_idx + dt - 1;
                    eventStruct.(evtName).charge_duration = dt;
                    eventStruct.(evtName).avg_current = mean(PeakCurrent{i});
                    eventStruct.(evtName).t_seq = rackData.Time(eventStruct.(evtName).start_idx:eventStruct.(evtName).end_idx);
                    eventStruct.(evtName).I_seq = PeakCurrent{i};
                    eventStruct.(evtName).V_seq = PeakVoltage{i};
                    eventStruct.(evtName).DCIR_Onori = dcir_val;
                    eventStruct.(evtName).DCIR_Onori_DV = DV;
                    eventStruct.(evtName).DCIR_Onori_DI = DI;
                    eventStruct.(evtName).DCIR_Onori_V1 = PeakVoltage{i}(1);
                    eventStruct.(evtName).DCIR_Onori_V2 = PeakVoltage{i}(end);
                    eventStruct.(evtName).DCIR_Onori_I1 = PeakCurrent{i}(1);
                    eventStruct.(evtName).DCIR_Onori_I2 = PeakCurrent{i}(end);
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
    
    sgtitle(sprintf('Date: %s, Peaks: %d', formatted_date, rack01_peaks_for_day), 'FontSize', 16);
    
    fig_filename = fullfile(saveDir, 'Onori_Method_Visualization.fig');
    saveas(gcf, fig_filename);
    fprintf('Figure saved to: %s\n', fig_filename);
    
    fprintf('Visualization complete: %d peaks detected.\n', length(PeakTime));
else
    fprintf('No peaks detected.\n');
end

save(fullfile(saveDir, 'all_chg_events_onori_all_years.mat'), 'global_eventStruct');
fprintf('Processing complete!\n');
fprintf('Results saved to: %s\n', saveDir);
