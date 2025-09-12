%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RackData_DCIR_Charge_Onori_OldData_v2.m
% Battery Rack DCIR Analysis for Old Data v2
% 1-second peak detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

%% Parameters
% Peak detection parameters
Cnom = 128;                   % Rack nominal Capacity (Ah)
thr = [Cnom*0.05, Cnom*0.1]; % Idle threshold [charge, discharge] (A)
dI_chg = Cnom * 0.15;         % Charge current change threshold (A)        I(i+1) - I(i) > 16.64A
dI_dch = Cnom * 0.15;         % Discharge current change threshold (A)     I(i+1) - I(i) < -25.6A
dt = 1;   % Peak duration (1 sample for 1Hz data = 1 second)
% Visualization parameters
Fontsize = 12;
LineWidth = 2;

%% Data paths
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR','Total_Peak_OldData_v2');

% Create save directory if it doesn't exist
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Initialize data structures
yearList = {'2022', '2023'};
rackNames_all = {'Rack01'}; % , 'Rack02', 'Rack03', 'Rack04', 'Rack05', 'Rack06', 'Rack07', 'Rack08'};

% Global event structure for all years
global_eventStruct.ChgPeak = struct();
global_eventStruct.DchPeak = struct();
for r = 1:length(rackNames_all)
    global_eventStruct.ChgPeak.(rackNames_all{r}) = struct();
    global_eventStruct.DchPeak.(rackNames_all{r}) = struct();
end

% Visualization variables
all_peak_times = {};
all_peak_currents = {};
all_peak_voltages = {};
all_raw_times = {};
all_raw_currents = {};
all_raw_voltages = {};

all_original_data = {};
all_peak_counts = [];
all_filenames = {};

all_dcir_values = [];


%% Process each year
for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    fprintf('==================================================\n');
    fprintf('Processing year: %s\n', year);

    % 연도별 SOC 통계 초기화
    year_soc_all = [];
    year_files_count = 0;

    % Initialize eventStruct for current year only
    eventStruct = struct();
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

            rackNames = rackNames_all;

            daily_original_data = {};
            daily_eventCount = 0;

            for rack_idx = 1:length(rackNames)
                rackName = rackNames{rack_idx};
                fprintf('  %s: ', rackName);

                if ~isfield(Raw_Rack, rackName)
                    fprintf('0 chgPeaks 0 dchPeaks (no data)\n');
                    continue;
                end

                rackData = Raw_Rack.(rackName);
                t = rackData.Time;
                I = rackData.DCCurrent_A;
                V = rackData.AverageCV_V;
                soc = rackData.SOCPct;
                T_batt = rackData.AverageMT_degC;
                P = rackData.DCPower_kW;

                % SOC 범위 디버깅
                fprintf('    SOC Range: %.1f%% - %.1f%% (Mean: %.1f%%)\n', ...
                    min(soc), max(soc), mean(soc));

                % 연도별 SOC 데이터 수집
                year_soc_all = [year_soc_all; soc(:)];
                year_files_count = year_files_count + 1;

                % Convert time data to seconds-from-start
                t = datetime(t);
                t = seconds(t - t(1));

                %% Charge Peak Detection (1 second peak resistance logic)
                PeakTime = {};
                PeakCurrent = {};
                PeakVoltage = {};
                PeakPower_chg = {};
                PeakSoC = [];
                PeakTemp = [];
                z = 1;

                % SOC 필터링 범위 디버깅
                soc_in_chg_range = sum(soc >= 65 & soc <= 70);
                fprintf('    Charge SOC Filter (65-70%%): %d/%d samples (%.1f%%)\n', ...
                    soc_in_chg_range, length(soc), 100*soc_in_chg_range/length(soc));

                for i = 1:(length(I)-dt)
                    % Check if total current change from idle to dt later exceeds threshold
                    if I(i+dt) - I(i) > dI_chg
                        % Check if current is in idle state (using raw current)
                        if I(i) > -thr(1) && I(i) < thr(1)
                            % The variable flag is introduced to check if the peak is always increasing and the current never changes sign
                            flag = 1;
                            for zi = 1:dt
                                if I(i+zi-1) < 0 || I(i+zi) < 0
                                    flag = 0;
                                end
                            end
                            if flag == 1
                                % Store peak data for two time points (start and 1 second later)
                                if z == 1
                                    % Store peak data for two time points
                                    PeakTime{z}(1) = t(i);      % Start time
                                    PeakTime{z}(2) = t(i+1);    % 1 second later
                                    PeakCurrent{z}(1) = I(i);   % Start current
                                    PeakCurrent{z}(2) = I(i+1); % 1 second later current
                                    PeakVoltage{z}(1) = V(i);   % Start voltage
                                    PeakVoltage{z}(2) = V(i+1); % 1 second later voltage
                                    PeakPower_chg{z}(1) = P(i); % Start power
                                    PeakPower_chg{z}(2) = P(i+1); % 1 second later power
                                    PeakSoC(z) = soc(i);
                                    PeakTemp(z) = T_batt(i);
                                    z = z + 1;
                                else
                                    % Additional condition that discards a new peak if a sufficient amount of time (dt=1s) is not passed from the beginning of the last peak
                                    if (t(i) - PeakTime{z-1}(1)) > dt
                                        % Store peak data for two time points
                                        PeakTime{z}(1) = t(i);      % Start time
                                        PeakTime{z}(2) = t(i+1);    % 1 second later
                                        PeakCurrent{z}(1) = I(i);   % Start current
                                        PeakCurrent{z}(2) = I(i+1); % 1 second later current
                                        PeakVoltage{z}(1) = V(i);   % Start voltage
                                        PeakVoltage{z}(2) = V(i+1); % 1 second later voltage
                                        PeakPower_chg{z}(1) = P(i); % Start power
                                        PeakPower_chg{z}(2) = P(i+1); % 1 second later power
                                        PeakSoC(z) = soc(i);
                                        PeakTemp(z) = T_batt(i);
                                        z = z + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % Save charge events (all events like reference)
            valid_chg_count = 0;  % Initialize counter
            for i = 1:(z-1)
                % Calculate DCIR for all peaks (like reference)
                dV = PeakVoltage{i}(end) - PeakVoltage{i}(1);
                dI = PeakCurrent{i}(end) - PeakCurrent{i}(1);
                dcir_val = (dV / dI) * 1000; % mΩ

                % Save all events regardless of DCIR value (like reference)
                valid_chg_count = valid_chg_count + 1;

                % Use file date YYYYMMDD for per-day keys
                date_str = matFiles(f).name;
                if length(date_str) >= 12
                    date_str = date_str(5:12);
                else
                    date_str = sprintf('%s%s01', year, monthDirs(m).name); % fallback YYYYMM01
                end
                evtName = sprintf('ChgEvent_%s_%s_%d', rackName, date_str, valid_chg_count);
                eventStruct.(evtName).rack_name = rackName;
                eventStruct.(evtName).date = date_str;
                eventStruct.(evtName).type = 'charge';
                eventStruct.(evtName).start_idx = i;  % Peak detection index
                eventStruct.(evtName).end_idx = i+dt-1;  % Peak end index (dt time steps later)
                eventStruct.(evtName).peak_duration = PeakTime{i}(end) - PeakTime{i}(1);
                eventStruct.(evtName).avg_current = mean(PeakCurrent{i});

                % Store peak data for dt time steps
                eventStruct.(evtName).t_seq = PeakTime{i};
                eventStruct.(evtName).I_seq = PeakCurrent{i};
                eventStruct.(evtName).V_seq = PeakVoltage{i};
                eventStruct.(evtName).P_seq = PeakPower_chg{i};

                % Store DCIR values
                eventStruct.(evtName).dV = dV;
                eventStruct.(evtName).dI = dI;
                eventStruct.(evtName).PeakChgR = dcir_val;
            end

            %% Discharge Peak Detection (1 second peak resistance logic)
            PeakTimeDis = {};
            PeakCurrentDis = {};
            PeakVoltageDis = {};
            PeakPower_dch = {};
            PeakSoCDis = [];
            PeakTempDis = [];
            zd = 1;

            % SOC 필터링 범위 디버깅
            soc_in_dch_range = sum(soc >= 60 & soc <= 70);
            fprintf('    Discharge SOC Filter (60-70%%): %d/%d samples (%.1f%%)\n', ...
                soc_in_dch_range, length(soc), 100*soc_in_dch_range/length(soc));

            for i = 1:(length(I)-dt)
                % Check if total current change from idle to dt later exceeds threshold
                dI_actual = I(i+dt) - I(i);
                if dI_actual < -dI_dch
                    % Check if current is in idle state (using raw current)
                    if I(i) > -thr(2) && I(i) < thr(2)
                        % The variable flag is introduced to check if the peak is always decreasing and the current never changes sign
                        flag = 1;
                        for zi = 1:dt
                            if I(i+zi-1) > 0 || I(i+zi) > 0
                                flag = 0;
                            end
                        end
                        if flag == 1
                            % Store peak data for two time points (start and 1 second later)
                            if zd == 1
                                % Store peak data for two time points
                                PeakTimeDis{zd}(1) = t(i);      % Start time
                                PeakTimeDis{zd}(2) = t(i+1);    % 1 second later
                                PeakCurrentDis{zd}(1) = I(i);   % Start current
                                PeakCurrentDis{zd}(2) = I(i+1); % 1 second later current
                                PeakVoltageDis{zd}(1) = V(i);   % Start voltage
                                PeakVoltageDis{zd}(2) = V(i+1); % 1 second later voltage
                                PeakPower_dch{zd}(1) = P(i);    % Start power
                                PeakPower_dch{zd}(2) = P(i+1);  % 1 second later power
                                PeakSoCDis(zd) = soc(i);
                                PeakTempDis(zd) = T_batt(i);
                                zd = zd + 1;
                            else
                                % Additional condition that discards a new peak if a sufficient amount of time (dt=1s) is not passed from the beginning of the last peak
                                if (t(i) - PeakTimeDis{zd-1}(1)) > dt
                                    % Store peak data for two time points
                                    PeakTimeDis{zd}(1) = t(i);      % Start time
                                    PeakTimeDis{zd}(2) = t(i+1);    % 1 second later
                                    PeakCurrentDis{zd}(1) = I(i);   % Start current
                                    PeakCurrentDis{zd}(2) = I(i+1); % 1 second later current
                                    PeakVoltageDis{zd}(1) = V(i);   % Start voltage
                                    PeakVoltageDis{zd}(2) = V(i+1); % 1 second later voltage
                                    PeakPower_dch{zd}(1) = P(i);    % Start power
                                    PeakPower_dch{zd}(2) = P(i+1);  % 1 second later power
                                    PeakSoCDis(zd) = soc(i);
                                    PeakTempDis(zd) = T_batt(i);
                                    zd = zd + 1;
                                end
                            end
                        end
                    end
                end
            end
        end

        % Save discharge events (all events like reference)
        valid_dch_count = 0;  % Initialize counter
        for i = 1:(zd-1)
            % Calculate DCIR for all peaks (like reference)
            dV = PeakVoltageDis{i}(end) - PeakVoltageDis{i}(1);
            dI = PeakCurrentDis{i}(end) - PeakCurrentDis{i}(1);
            dcir_val = (dV / dI) * 1000; % mΩ

            % Save all events regardless of DCIR value (like reference)
            valid_dch_count = valid_dch_count + 1;

            % Use file date YYYYMMDD for per-day keys
            date_str = matFiles(f).name;
            if length(date_str) >= 12
                date_str = date_str(5:12);
            else
                date_str = sprintf('%s%s01', year, monthDirs(m).name); % fallback YYYYMM01
            end
            evtName = sprintf('DisEvent_%s_%s_%d', rackName, date_str, valid_dch_count);
            eventStruct.(evtName).rack_name = rackName;
            eventStruct.(evtName).date = date_str;
            eventStruct.(evtName).type = 'discharge';
            eventStruct.(evtName).start_idx = i;  % Peak detection index
            eventStruct.(evtName).end_idx = i+dt-1;  % Peak end index (dt time steps later)
            eventStruct.(evtName).peak_duration = PeakTimeDis{i}(end) - PeakTimeDis{i}(1);
            eventStruct.(evtName).avg_current = mean(PeakCurrentDis{i});

            % Store peak data for dt time steps
            eventStruct.(evtName).t_seq = PeakTimeDis{i};
            eventStruct.(evtName).I_seq = PeakCurrentDis{i};
            eventStruct.(evtName).V_seq = PeakVoltageDis{i};
            eventStruct.(evtName).P_seq = PeakPower_dch{i};

            % Store DCIR values
            eventStruct.(evtName).dV = dV;
            eventStruct.(evtName).dI = dI;
            eventStruct.(evtName).PeakDisR = dcir_val;
        end

        % Debug output for each rack
        fprintf('%d chgPeaks %d dchPeaks\n', z-1, zd-1);

        % Store original data for visualization
        if ~isempty(PeakTime) || ~isempty(PeakTimeDis)
            daily_original_data{end+1} = struct('time', t, 'voltage', V, 'current', I);
        end
    end


    % Store daily data
    if ~isempty(daily_original_data)
        all_original_data{end+1} = daily_original_data;
        all_filenames{end+1} = matFiles(f).name;
    end

    % Save structure by rack/year/date
    year_str = sprintf('year_%s', year);

    % Save charge events (Rack -> ChgPeak -> EYYYYMM -> EYYYYMMDD -> EventN)
    for rack_idx = 1:length(rackNames_all)
    rackName = rackNames_all{rack_idx};
    eventNames = fieldnames(eventStruct);
    rack_events = {};

    for i = 1:length(eventNames)
        evt = eventStruct.(eventNames{i});

        if isfield(evt,'type') && strcmp(evt.type,'charge') && strcmp(evt.rack_name, rackName)
            rack_events{end+1} = eventNames{i};
        end
    end

    if isempty(rack_events), continue; end
    if ~isfield(global_eventStruct.ChgPeak, rackName)
        global_eventStruct.ChgPeak.(rackName) = struct();
    end
    for i = 1:length(rack_events)
        evt = eventStruct.(rack_events{i});
        % Keys
        month_key = ['E' evt.date(1:6)];
        date_key  = ['E' evt.date];
        if ~isfield(global_eventStruct.ChgPeak.(rackName), month_key)
            global_eventStruct.ChgPeak.(rackName).(month_key) = struct();
        end
        if ~isfield(global_eventStruct.ChgPeak.(rackName).(month_key), date_key)
            global_eventStruct.ChgPeak.(rackName).(month_key).(date_key) = struct();
        end
        existing = fieldnames(global_eventStruct.ChgPeak.(rackName).(month_key).(date_key));
        nextIdx = numel(existing) + 1;
        evKey = sprintf('Event%d', nextIdx);
        global_eventStruct.ChgPeak.(rackName).(month_key).(date_key).(evKey) = evt;
    end
end

% Save discharge events (Rack -> DchPeak -> EYYYYMM -> EYYYYMMDD -> EventN)
for rack_idx = 1:length(rackNames_all)
    rackName = rackNames_all{rack_idx};
    eventNames = fieldnames(eventStruct);
    rack_events = {};

    for i = 1:length(eventNames)
        evt = eventStruct.(eventNames{i});
        if isfield(evt,'type') && strcmp(evt.type,'discharge') && strcmp(evt.rack_name, rackName)
            rack_events{end+1} = eventNames{i};
        end
    end

    if isempty(rack_events), continue; end
    if ~isfield(global_eventStruct.DchPeak, rackName)
        global_eventStruct.DchPeak.(rackName) = struct();
    end
    for i = 1:length(rack_events)
        evt = eventStruct.(rack_events{i});
        month_key = ['E' evt.date(1:6)];
        date_key  = ['E' evt.date];
        if ~isfield(global_eventStruct.DchPeak.(rackName), month_key)
            global_eventStruct.DchPeak.(rackName).(month_key) = struct();
        end
        if ~isfield(global_eventStruct.DchPeak.(rackName).(month_key), date_key)
            global_eventStruct.DchPeak.(rackName).(month_key).(date_key) = struct();
        end
        existing = fieldnames(global_eventStruct.DchPeak.(rackName).(month_key).(date_key));
        nextIdx = numel(existing) + 1;
        evKey = sprintf('Event%d', nextIdx);
        global_eventStruct.DchPeak.(rackName).(month_key).(date_key).(evKey) = evt;
    end
end

% 연도별 SOC 통계 요약
if ~isempty(year_soc_all)
    fprintf('\n--- Year %s SOC Summary ---\n', year);
    fprintf('Files processed: %d\n', year_files_count);
    fprintf('Total SOC samples: %d\n', length(year_soc_all));
    fprintf('SOC Range: %.1f%% - %.1f%%\n', min(year_soc_all), max(year_soc_all));
    fprintf('SOC Mean: %.1f%% (±%.1f%%)\n', mean(year_soc_all), std(year_soc_all));

    % SOC 구간별 분포
    soc_ranges = [0 20; 20 40; 40 60; 60 80; 80 100];
    fprintf('SOC Distribution:\n');
    for sr = 1:size(soc_ranges, 1)
        count = sum(year_soc_all >= soc_ranges(sr,1) & year_soc_all < soc_ranges(sr,2));
        fprintf('  %d-%d%%: %d samples (%.1f%%)\n', ...
            soc_ranges(sr,1), soc_ranges(sr,2), count, 100*count/length(year_soc_all));
    end

    % 현재 필터링 범위 적중률
    chg_range_count = sum(year_soc_all >= 50 & year_soc_all <= 70);
    dch_range_count = sum(year_soc_all >= 60 & year_soc_all <= 90);
    fprintf('Current Filter Coverage:\n');
    fprintf('  Charge (50-70%%): %d samples (%.1f%%)\n', ...
        chg_range_count, 100*chg_range_count/length(year_soc_all));
    fprintf('  Discharge (60-90%%): %d samples (%.1f%%)\n', ...
        dch_range_count, 100*dch_range_count/length(year_soc_all));
end

% Save current year data
save(fullfile(saveDir, sprintf('Peak_OldData_v2_%s.mat', year)), 'global_eventStruct');
fprintf('\nYear %s processing complete\n', year);
fprintf('==================================================\n');

% Clear year-specific eventStruct to free memory, but keep global_eventStruct
clear eventStruct;
end

% Save total combined data from all years
save(fullfile(saveDir, sprintf('Peak_OldData_v2_%s_total.mat', strjoin(yearList,'_'))), 'global_eventStruct');
fprintf('Total data processing complete\n');

%% Visualization with Year-specific 3σ Outlier Removal (Reference Style)
% Reference-style histfit visualization for each rack with year-specific outlier removal
for r = 1:length(rackNames_all)
    rackName = rackNames_all{r};

    % Collect charge peak data by year
    chg_vals_by_year = cell(1, length(yearList));
    for yi = 1:length(yearList)
        current_year = yearList{yi};
        chg_vals = [];
        if isfield(global_eventStruct.ChgPeak, rackName)
            month_names = fieldnames(global_eventStruct.ChgPeak.(rackName));
            fprintf('Debug: Year %s, Found months: %s\n', current_year, strjoin(month_names, ', '));
            for m = 1:length(month_names)
                month_key = month_names{m};
                % Check if this month belongs to current year
                if length(month_key) >= 7 && strcmp(month_key(2:5), current_year)
                    date_names = fieldnames(global_eventStruct.ChgPeak.(rackName).(month_key));
                    for d = 1:length(date_names)
                        dayEvents = global_eventStruct.ChgPeak.(rackName).(month_key).(date_names{d});
                        evtNamesDay = fieldnames(dayEvents);
                        for e = 1:length(evtNamesDay)
                            evt = dayEvents.(evtNamesDay{e});
                            if isfield(evt, 'PeakChgR') && ~isnan(evt.PeakChgR) && ~isinf(evt.PeakChgR) && evt.PeakChgR > 0
                                chg_vals = [chg_vals, evt.PeakChgR];
                            end
                        end
                    end
                end
            end
        end
        chg_vals_by_year{yi} = chg_vals;
    end

    % Collect discharge peak data by year
    dch_vals_by_year = cell(1, length(yearList));
    for yi = 1:length(yearList)
        current_year = yearList{yi};
        dch_vals = [];
        if isfield(global_eventStruct.DchPeak, rackName)
            month_names = fieldnames(global_eventStruct.DchPeak.(rackName));
            fprintf('Debug: Year %s, Found months: %s\n', current_year, strjoin(month_names, ', '));
            for m = 1:length(month_names)
                month_key = month_names{m};
                % Check if this month belongs to current year
                if length(month_key) >= 7 && strcmp(month_key(2:5), current_year)
                    date_names = fieldnames(global_eventStruct.DchPeak.(rackName).(month_key));
                    for d = 1:length(date_names)
                        dayEvents = global_eventStruct.DchPeak.(rackName).(month_key).(date_names{d});
                        evtNamesDay = fieldnames(dayEvents);
                        for e = 1:length(evtNamesDay)
                            evt = dayEvents.(evtNamesDay{e});
                            if isfield(evt, 'PeakDisR') && ~isnan(evt.PeakDisR) && ~isinf(evt.PeakDisR) && evt.PeakDisR > 0
                                dch_vals = [dch_vals, evt.PeakDisR];
                            end
                        end
                    end
                end
            end
        end
        dch_vals_by_year{yi} = dch_vals;
    end

    % Apply 3σ outlier removal for each year (Reference style with while loop)
    chg_vals_cleaned = [];
    dch_vals_cleaned = [];
    chg_vals_cleaned_by_year = cell(1, length(yearList));
    dch_vals_cleaned_by_year = cell(1, length(yearList));

    % Charge peaks - year-specific outlier removal (Reference style)
    for yi = 1:length(yearList)
        year_vals = chg_vals_by_year{yi};
        if ~isempty(year_vals)
            % 필터링 전 통계 (1단계에서 이미 NaN/Inf/음수 제거됨)
            original_count = length(year_vals);
            nan_count = sum(isnan(year_vals));
            inf_count = sum(isinf(year_vals));
            neg_count = sum(year_vals < 0);

            % Remove NaN and Inf values first (allow negative values like reference)
            year_vals = year_vals(~isnan(year_vals) & ~isinf(year_vals));

            fprintf('    Debug Charge %s: NaN=%d, Inf=%d, Neg=%d (should be 0)\n', ...
                yearList{yi}, nan_count, inf_count, neg_count);

            if ~isempty(year_vals) && length(year_vals) >= 3
                % Fit normal distribution and calculate 3σ bounds (Reference style)
                pd = fitdist(year_vals', 'Normal');
                out_max = pd.mu + 3*pd.sigma;
                out_min = pd.mu - 3*pd.sigma;  % Allow negative values like reference

                % Remove outliers using while loop (Reference style)
                r = 1;
                while r <= length(year_vals)
                    if(year_vals(r) > out_max || year_vals(r) < out_min)
                        year_vals(r) = [];
                    else
                        r = r + 1;
                    end
                end
                chg_vals_cleaned = [chg_vals_cleaned, year_vals];
                chg_vals_cleaned_by_year{yi} = year_vals;

                fprintf('  %s Charge %s: Original %d, Removed %d outliers (3σ), Final %d\n', ...
                    rackName, yearList{yi}, length(chg_vals_by_year{yi}), length(chg_vals_by_year{yi}) - length(year_vals), length(year_vals));
            else
                % 데이터가 부족하면 3σ 제거 없이 원본 데이터 사용
                chg_vals_cleaned = [chg_vals_cleaned, year_vals];
                chg_vals_cleaned_by_year{yi} = year_vals;
                fprintf('  %s Charge %s: Not enough data (%d samples) for outlier removal, using original data\n', ...
                    rackName, yearList{yi}, length(year_vals));
            end
        end
    end

    % Discharge peaks - year-specific outlier removal (Reference style)
    for yi = 1:length(yearList)
        year_vals = dch_vals_by_year{yi};
        if ~isempty(year_vals)
            % Remove NaN and Inf values first (allow negative values like reference)
            year_vals = year_vals(~isnan(year_vals) & ~isinf(year_vals));

            if ~isempty(year_vals) && length(year_vals) >= 3
                % Fit normal distribution and calculate 3σ bounds (Reference style)
                pd = fitdist(year_vals', 'Normal');
                out_max = pd.mu + 3*pd.sigma;
                out_min = pd.mu - 3*pd.sigma;  % Allow negative values like reference

                % Remove outliers using while loop (Reference style)
                r = 1;
                while r <= length(year_vals)
                    if(year_vals(r) > out_max || year_vals(r) < out_min)
                        year_vals(r) = [];
                    else
                        r = r + 1;
                    end
                end
                dch_vals_cleaned = [dch_vals_cleaned, year_vals];
                dch_vals_cleaned_by_year{yi} = year_vals;

                fprintf('  %s Discharge %s: Original %d, Removed %d outliers (3σ), Final %d\n', ...
                    rackName, yearList{yi}, length(dch_vals_by_year{yi}), length(dch_vals_by_year{yi}) - length(year_vals), length(year_vals));
            else
                % 데이터가 부족하면 3σ 제거 없이 원본 데이터 사용
                dch_vals_cleaned = [dch_vals_cleaned, year_vals];
                dch_vals_cleaned_by_year{yi} = year_vals;
                fprintf('  %s Discharge %s: Not enough data (%d samples) for outlier removal, using original data\n', ...
                    rackName, yearList{yi}, length(year_vals));
            end
        end
    end

    % Create histfit plots with min-max/30 binning for cleaned data
    % Create year-overlay histfit plots with year-specific min-max/30
    if ~isempty(chg_vals_cleaned)
        figure('Name', sprintf('Hist_Chg_%s_3sigma', rackName), 'NumberTitle', 'off', 'Visible', 'on');
        fprintf('Debug: Created charge histogram figure\n');
        hold on;

        % Colors for years - Charge (Reference style)
        yearColors = {'#0073C2', '#FF8C00'};  % Blue for 2021, Orange for 2022

        % Plot cleaned data with year-specific colors
        legend_handles = [];

        % Use existing cleaned data for plotting
        % chg_vals_cleaned_by_year is already populated from the outlier removal section

        % Plot each year's cleaned data with overlay (Reference style)
        legend_handles = [];
        for yi = 1:length(yearList)
            if ~isempty(chg_vals_cleaned_by_year{yi})
                % Use histfit with 30 bins
                n_bins = 30;  % 30개 bin
                h = histfit(chg_vals_cleaned_by_year{yi}, n_bins);
                h(1).EdgeColor = yearColors{yi};
                h(1).FaceColor = yearColors{yi};
                h(1).FaceAlpha = 0.7;
                h(2).Color = yearColors{yi};
                h(2).LineWidth = 4;

                % Store legend handle
                legend_handles = [legend_handles, h(1)];

                % Calculate and display mean and std
                mean_val = mean(chg_vals_cleaned_by_year{yi});
                std_val = std(chg_vals_cleaned_by_year{yi});

                % Add text annotation
                text(0.02, 0.98-yi*0.1, sprintf('%s: Mean=%.4f, STD=%.4f', yearList{yi}, mean_val, std_val), ...
                    'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
                    'Color', 'black', 'VerticalAlignment', 'top');
            end
        end

        % Add legend
        legend(legend_handles, yearList, 'Location', 'best');

        xlabel('R_{CHG} [m\Omega]'); ylabel('Frequency [-]');
        title(sprintf('Charge Peaks - %s (3σ outliers removed)', rackName));
        xlim([0 inf]);  % Reference style display range
        set(findall(gcf,'-property','FontSize'),'FontSize',20);
        set(findall(gcf,'-property','interpreter'),'interpreter','tex');
        set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex');
        saveas(gcf, fullfile(saveDir, sprintf('Hist_Chg_%s_3sigma.fig', rackName)));
    end

    if ~isempty(dch_vals_cleaned)
        figure('Name', sprintf('Hist_Dch_%s_3sigma', rackName), 'NumberTitle', 'off', 'Visible', 'on');
        hold on;

        % Colors for years - Discharge (Reference style)
        yearColors = {'#0073C2', '#FF8C00'};  % Blue for 2021, Orange for 2022

        % Plot cleaned data with year-specific colors
        legend_handles = [];

        % Use existing cleaned data for plotting
        % dch_vals_cleaned_by_year is already populated from the outlier removal section

        % Plot each year's cleaned data with overlay (Reference style)
        legend_handles = [];
        for yi = 1:length(yearList)
            if ~isempty(dch_vals_cleaned_by_year{yi})
                % Use histfit with 30 bins
                n_bins = 30;  % 30개 bin
                h = histfit(dch_vals_cleaned_by_year{yi}, n_bins);
                h(1).EdgeColor = yearColors{yi};
                h(1).FaceColor = yearColors{yi};
                h(1).FaceAlpha = 0.7;
                h(2).Color = yearColors{yi};
                h(2).LineWidth = 4;

                % Store legend handle
                legend_handles = [legend_handles, h(1)];

                % Calculate and display mean and std
                mean_val = mean(dch_vals_cleaned_by_year{yi});
                std_val = std(dch_vals_cleaned_by_year{yi});

                % Add text annotation
                text(0.02, 0.98-yi*0.1, sprintf('%s: Mean=%.4f, STD=%.4f', yearList{yi}, mean_val, std_val), ...
                    'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
                    'Color', 'black', 'VerticalAlignment', 'top');
            end
        end

        % Add legend
        legend(legend_handles, yearList, 'Location', 'best');

        xlabel('R_{DCH} [m\Omega]'); ylabel('Frequency [-]');
        title(sprintf('Discharge Peaks - %s (3σ outliers removed)', rackName));
        xlim([0 inf]);  % Reference style display range
        set(findall(gcf,'-property','FontSize'),'FontSize',20);
        set(findall(gcf,'-property','interpreter'),'interpreter','tex');
        set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex');
        saveas(gcf, fullfile(saveDir, sprintf('Hist_Dch_%s_3sigma.fig', rackName)));
    end

    %% Boxplot Visualization with Year-specific 3σ Outlier Removal
    % Create boxplot figures for each rack with year-specific outlier removal

    % Charge peaks boxplot
    if ~isempty(chg_vals_cleaned)
        figure('Name', sprintf('Boxplot_Chg_%s_3sigma', rackName), 'NumberTitle', 'off', 'Visible', 'on');
        hold on;

        % Use existing cleaned data for boxplot
        if ~isempty(chg_vals_cleaned)
            % Prepare data for daboxplot
            daboxplot_data = [];
            daboxplot_labels = {};

            for yi = 1:length(yearList)
                if ~isempty(chg_vals_cleaned_by_year{yi})
                    daboxplot_data = [daboxplot_data, chg_vals_cleaned_by_year{yi}];
                    daboxplot_labels = [daboxplot_labels, repmat({yearList{yi}}, 1, length(chg_vals_cleaned_by_year{yi}))];
                end
            end

            % Create daboxplot with existing data
            h = daboxplot(daboxplot_data, daboxplot_labels);

            % Set transparency for daboxplot
            for i = 1:length(h)
                if isfield(h(i), 'FaceAlpha')
                    h(i).FaceAlpha = 0.5;
                end
            end

            % Set x-axis labels to actual years
            unique_labels = unique(daboxplot_labels);
            set(gca, 'XTick', 1:length(unique_labels), 'XTickLabel', unique_labels);

            % Calculate and plot mean values with std
            mean_values = [];
            std_values = [];
            for yi = 1:length(yearList)
                if ~isempty(chg_vals_cleaned_by_year{yi})
                    mean_values = [mean_values, mean(chg_vals_cleaned_by_year{yi})];
                    std_values = [std_values, std(chg_vals_cleaned_by_year{yi})];
                end
            end

            % Plot mean line
            hold on;
            plot(1:length(mean_values), mean_values, 'r-', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'red');

            % Add mean and std values as text (same position as histogram)
            for i = 1:length(mean_values)
                text(0.02, 0.98-i*0.1, sprintf('%s: Mean=%.4f, STD=%.4f', yearList{i}, mean_values(i), std_values(i)), ...
                    'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
                    'Color', 'black', 'VerticalAlignment', 'top');
            end

            xlabel('Year'); ylabel('R_{CHG} [m\Omega]');
            title(sprintf('Charge Peaks Boxplot - %s (3σ outliers removed)', rackName));
            % legend('Mean', 'Location', 'best');
            set(findall(gcf,'-property','FontSize'),'FontSize',20);
            set(findall(gcf,'-property','interpreter'),'interpreter','tex');
            set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex');
            saveas(gcf, fullfile(saveDir, sprintf('Boxplot_Chg_%s_3sigma.fig', rackName)));
        end
    end

    % Discharge peaks boxplot
    if ~isempty(dch_vals_cleaned)
        figure('Name', sprintf('Boxplot_Dch_%s_3sigma', rackName), 'NumberTitle', 'off', 'Visible', 'on');
        hold on;

        % Use existing cleaned data for boxplot
        if ~isempty(dch_vals_cleaned)
            % Prepare data for daboxplot
            daboxplot_data = [];
            daboxplot_labels = {};

            for yi = 1:length(yearList)
                if ~isempty(dch_vals_cleaned_by_year{yi})
                    daboxplot_data = [daboxplot_data, dch_vals_cleaned_by_year{yi}];
                    daboxplot_labels = [daboxplot_labels, repmat({yearList{yi}}, 1, length(dch_vals_cleaned_by_year{yi}))];
                end
            end

            % Create daboxplot with existing data
            h = daboxplot(daboxplot_data, daboxplot_labels);

            % Set transparency for daboxplot
            for i = 1:length(h)
                if isfield(h(i), 'FaceAlpha')
                    h(i).FaceAlpha = 0.5;
                end
            end

            % Set x-axis labels to actual years
            unique_labels = unique(daboxplot_labels);
            set(gca, 'XTick', 1:length(unique_labels), 'XTickLabel', unique_labels);

            % Calculate and plot mean values with std
            mean_values = [];
            std_values = [];
            for yi = 1:length(yearList)
                if ~isempty(dch_vals_cleaned_by_year{yi})
                    mean_values = [mean_values, mean(dch_vals_cleaned_by_year{yi})];
                    std_values = [std_values, std(dch_vals_cleaned_by_year{yi})];
                end
            end

            % Plot mean line
            hold on;
            plot(1:length(mean_values), mean_values, 'r-', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'red');

            % Add mean and std values as text (same position as histogram)
            for i = 1:length(mean_values)
                text(0.02, 0.98-i*0.1, sprintf('%s: Mean=%.4f, STD=%.4f', yearList{i}, mean_values(i), std_values(i)), ...
                    'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
                    'Color', 'black', 'VerticalAlignment', 'top');
            end

            xlabel('Year'); ylabel('R_{DCH} [m\Omega]');
            title(sprintf('Discharge Peaks Boxplot - %s (3σ outliers removed)', rackName));
            % legend('Mean', 'Location', 'best');
            set(findall(gcf,'-property','FontSize'),'FontSize',20);
            set(findall(gcf,'-property','interpreter'),'interpreter','tex');
            set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex');
            saveas(gcf, fullfile(saveDir, sprintf('Boxplot_Dch_%s_3sigma.fig', rackName)));
        end
    end
end
