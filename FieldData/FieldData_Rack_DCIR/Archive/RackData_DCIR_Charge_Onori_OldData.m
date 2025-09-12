%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RackData_DCIR_Charge_Onori_OldData.m
% Battery Rack DCIR Analysis for Old Data
% Based on fig_4_5.m methodology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

%% Parameters
% Peak detection parameters
Cnom = 128;                   % Rack nominal Capacity (Ah)
thr = [Cnom*0.03, Cnom*0.05]; % Idle threshold [charge, discharge] (A)
dI_chg = Cnom * 0.10;          % Charge current change threshold (A)
dI_dch = Cnom * 0.10;          % Discharge current change threshold (A)
ddI = 1;                       % Instantaneous current change threshold (A)
te = 8;   % Moving average window (6 samples for 1Hz data = 6 seconds)
dt = 1;   % Peak duration (1 sample for 1Hz data = 1 second)

% Visualization parameters
Fontsize = 12;
LineWidth = 2;

%% Data paths
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';

saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR','TotalPeak_OldData');


% Create save directory if it doesn't exist
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Initialize data structures
yearList = {'2021', '2022'};
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
    fprintf('Processing year: %s\n', year);

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

                % Convert time data to seconds-from-start
                t = datetime(t);
                t = seconds(t - t(1));

                % Filter on Current (Moving average of 1s)
                N = length(I);
                I_filt = zeros(N,1);
                Pre = 0*ones(te/2,1);
                Post = 0*ones(te/2,1);
                I_calc = [Pre;I;Post];
                for i = 1:N
                    for m = 0:te
                        I_filt(i) = I_filt(i) + I_calc(i+m);
                    end
                    I_filt(i) = I_filt(i)/(te+1);
                end

                % Derivative of Current
                for i = 1:length(I) - 1
                    dI_dt(i) = (I(i+1) - I(i)) / (t(i+1) - t(i));
                end

                % MA on Current Derivative
                N = length(dI_dt);
                filt_dI_dt = zeros(N,1);
                Pre = 0*ones(te/2,1);
                Post = 0*ones(te/2,1);
                calc_dI_dt  = [Pre;dI_dt';Post];
                for i = 1:N
                    for m = 0:te
                        filt_dI_dt(i) = filt_dI_dt(i) + calc_dI_dt(i+m);
                    end
                    filt_dI_dt(i) = filt_dI_dt(i)/(te+1);
                end

                %% Charge Peak Detection (1 second peak resistance logic)
                PeakTime = {};
                PeakCurrent = {};
                PeakVoltage = {};
                PeakPower_chg = {};
                PeakSoC = [];
                PeakTemp = [];
                z = 1;

                for i = 1:(length(I)-dt)
                    % Check if total current change from idle to 1 second later exceeds threshold (using moving average)
                    if (I_filt(i+dt) - I_filt(i)) > dI_chg
                        % Check if current is in idle state (using raw current)
                        if I(i) >-thr(1) && I(i)<thr(1)
                            % The variable ddI is introduced to ensure that a peak has an initial increase of at least 1A (peaks with initial flat plateaus are discarded)
                            if (I_filt(i+1) - I_filt(i)) > ddI
                                % The variable flag is introduced to check if the peak is always increasing and the current never changes sign
                                flag = 1;
                                for zi = 1:dt
                                    if filt_dI_dt(i+zi-1) < 0 || I(i+zi) < 0
                                        flag = 0;
                                    end
                                end
                                if flag == 1
                                    if z == 1
                                        % Store peak data for dt time steps
                                        for j=1:dt
                                            PeakTime{z}(j) = t(i+j-1);
                                            PeakCurrent{z}(j) = I(i+j-1);
                                            PeakVoltage{z}(j) = V(i+j-1);
                                            PeakPower_chg{z}(j) = P(i+j-1);
                                            
                                        end
                                        PeakSoC(z) = soc(i);
                                        PeakTemp(z) = T_batt(i);
                                        z = z + 1;
                                    else
                                        % Additional condition that discards a new peak if a sufficient amount of time (dt=1s) is not passed from the beginning of the last peak
                                        if (t(i) - PeakTime{z-1}(1)) > dt
                                            % Store peak data for dt time steps
                                            for j=1:dt
                                                PeakTime{z}(j) = t(i+j-1);
                                                PeakCurrent{z}(j) = I(i+j-1);
                                                PeakVoltage{z}(j) = V(i+j-1);
                                                PeakPower_chg{z}(j) = P(i+j-1);
                                            end
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

                % Save charge events
                for i = 1:(z-1)
                    % Use file date YYYYMMDD for per-day keys
                    date_str = matFiles(f).name;
                    if length(date_str) >= 12
                        date_str = date_str(5:12);
                    else
                        date_str = sprintf('%s%s01', year, monthDirs(m).name); % fallback YYYYMM01
                    end
                    evtName = sprintf('ChgEvent_%s_%s_%d', rackName, date_str, i);
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

                    % Store power and filtered current data for the same time points
                    start_time = PeakTime{i}(1);
                    end_time = PeakTime{i}(end);
                    time_indices = find(t >= start_time & t <= end_time);
                    eventStruct.(evtName).I_filt_seq = I_filt(time_indices);  % Moving average current data

                    % Calculate DCIR
                    dV = PeakVoltage{i}(end) - PeakVoltage{i}(1);
                    dI = PeakCurrent{i}(end) - PeakCurrent{i}(1);
                    eventStruct.(evtName).dV = dV;
                    eventStruct.(evtName).dI = dI;
                    if dI ~= 0
                        dcir_val = (dV / dI) * 1000; % mΩ
                    else
                        dcir_val = NaN;
                    end
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

                for i = 1:(length(I)-dt)
                    % Check if total current change from idle to 1 second later exceeds threshold (using moving average)
                    if (I_filt(i+dt) - I_filt(i)) < -dI_dch
                        % Check if current is in idle state (using raw current)
                        if I(i) >-thr(2) && I(i)<thr(2)
                            % The variable ddI is introduced to ensure that a peak has an initial decrease of at least -1A (peaks with initial flat plateaus are discarded)
                            if (I_filt(i+1) - I_filt(i)) < -ddI
                                % The variable flag is introduced to check if the peak is always decreasing and the current never changes sign
                                flag = 1;
                                for zi = 1:dt
                                    if filt_dI_dt(i+zi-1) > 0 || I(i+zi) > 0
                                        flag = 0;
                                    end
                                end
                                if flag == 1
                                    if zd == 1
                                        % Store peak data for dt time steps
                                        for j=1:dt
                                            PeakTimeDis{zd}(j) = t(i+j-1);
                                            PeakCurrentDis{zd}(j) = I(i+j-1);
                                            PeakVoltageDis{zd}(j) = V(i+j-1);
                                            PeakPower_dch{zd}(j) = P(i+j-1);

                                        end
                                        PeakSoCDis(zd) = soc(i);
                                        PeakTempDis(zd) = T_batt(i);
                                        zd = zd + 1;
                                    else
                                        % Additional condition that discards a new peak if a sufficient amount of time (dt=1s) is not passed from the beginning of the last peak
                                        if (t(i) - PeakTimeDis{zd-1}(1)) > dt
                                            % Store peak data for dt time steps
                                            for j=1:dt
                                                PeakTimeDis{zd}(j) = t(i+j-1);
                                                PeakCurrentDis{zd}(j) = I(i+j-1);
                                                PeakVoltageDis{zd}(j) = V(i+j-1);
                                                PeakPower_dch{zd}(j) = P(i+j-1);
                                            end
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
            end

            % Save discharge events
            for i = 1:(zd-1)
                % Use file date YYYYMMDD for per-day keys
                date_str = matFiles(f).name;
                if length(date_str) >= 12
                    date_str = date_str(5:12);
                else
                    date_str = sprintf('%s%s01', year, monthDirs(m).name); % fallback YYYYMM01
                end
                evtName = sprintf('DisEvent_%s_%s_%d', rackName, date_str, i);
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

                % Store power and filtered current data for the same time points
                start_time = PeakTimeDis{i}(1);
                end_time = PeakTimeDis{i}(end);
                time_indices = find(t >= start_time & t <= end_time);
                eventStruct.(evtName).I_filt_seq = I_filt(time_indices);  % Moving average current data

                % Calculate DCIR
                dV = PeakVoltageDis{i}(end) - PeakVoltageDis{i}(1);
                dI = PeakCurrentDis{i}(end) - PeakCurrentDis{i}(1);
                eventStruct.(evtName).dV = dV;
                eventStruct.(evtName).dI = dI;
                if dI ~= 0
                    dcir_val = (dV / dI) * 1000; % mΩ
                else
                    dcir_val = NaN;

                end
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

    % User override: do not remove outliers; keep all events in structure

    % Save current year data
    save(fullfile(saveDir, sprintf('Peak_OldData_%s.mat', year)), 'global_eventStruct');
    fprintf('Year %s processing complete\n', year);

    % Clear year-specific eventStruct to free memory, but keep global_eventStruct
    clear eventStruct;
end

% Save total combined data from all years
save(fullfile(saveDir, sprintf('Peak_OldData_%s_total.mat', strjoin(yearList,'_'))), 'global_eventStruct');
fprintf('Total data processing complete\n');

%% Visualization with Year-specific 3σ Outlier Removal (Reference Style)
% Reference-style histfit visualization for each rack with year-specific outlier removal
for r = 1:length(rackNames_all)
    rackName = rackNames_all{r};

    % Collect charge peak data by year
    chg_vals_by_year = cell(1, length(yearList));
    for yi = 1:length(yearList)
        year = yearList{yi};
        chg_vals = [];
        if isfield(global_eventStruct.ChgPeak, rackName)
            month_names = fieldnames(global_eventStruct.ChgPeak.(rackName));
            for m = 1:length(month_names)
                month_key = month_names{m};
                % Check if this month belongs to current year
                if length(month_key) >= 7 && strcmp(month_key(2:5), year)
                    date_names = fieldnames(global_eventStruct.ChgPeak.(rackName).(month_key));
                    for d = 1:length(date_names)
                        dayEvents = global_eventStruct.ChgPeak.(rackName).(month_key).(date_names{d});
                        evtNamesDay = fieldnames(dayEvents);
                        for e = 1:length(evtNamesDay)
                            evt = dayEvents.(evtNamesDay{e});
                            if isfield(evt, 'PeakChgR') && ~isnan(evt.PeakChgR)
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
        year = yearList{yi};
        dch_vals = [];
        if isfield(global_eventStruct.DchPeak, rackName)
            month_names = fieldnames(global_eventStruct.DchPeak.(rackName));
            for m = 1:length(month_names)
                month_key = month_names{m};
                % Check if this month belongs to current year
                if length(month_key) >= 7 && strcmp(month_key(2:5), year)
                    date_names = fieldnames(global_eventStruct.DchPeak.(rackName).(month_key));
                    for d = 1:length(date_names)
                        dayEvents = global_eventStruct.DchPeak.(rackName).(month_key).(date_names{d});
                        evtNamesDay = fieldnames(dayEvents);
                        for e = 1:length(evtNamesDay)
                            evt = dayEvents.(evtNamesDay{e});
                            if isfield(evt, 'PeakDisR') && ~isnan(evt.PeakDisR)
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
    
    % Charge peaks - year-specific outlier removal (Reference style)
    for yi = 1:length(yearList)
        year_vals = chg_vals_by_year{yi};
        if ~isempty(year_vals)
            % Fit normal distribution and calculate 3σ bounds (Reference style)
            pd = fitdist(year_vals', 'Normal');
            out_max = pd.mu + 3*pd.sigma;
            out_min = pd.mu - 3*pd.sigma;
            
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
            
            fprintf('  %s Charge %s: Original %d, Removed %d outliers (3σ), Final %d\n', ...
                rackName, yearList{yi}, length(chg_vals_by_year{yi}), length(chg_vals_by_year{yi}) - length(year_vals), length(year_vals));
        end
    end
    
    % Discharge peaks - year-specific outlier removal (Reference style)
    for yi = 1:length(yearList)
        year_vals = dch_vals_by_year{yi};
        if ~isempty(year_vals)
            % Fit normal distribution and calculate 3σ bounds (Reference style)
            pd = fitdist(year_vals', 'Normal');
            out_max = pd.mu + 3*pd.sigma;
            out_min = pd.mu - 3*pd.sigma;
            
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
            
            fprintf('  %s Discharge %s: Original %d, Removed %d outliers (3σ), Final %d\n', ...
                rackName, yearList{yi}, length(dch_vals_by_year{yi}), length(dch_vals_by_year{yi}) - length(year_vals), length(year_vals));
        end
    end

        % Create histfit plots with min-max/30 binning for cleaned data
    % Create year-overlay histfit plots with year-specific min-max/30
    if ~isempty(chg_vals_cleaned)
        figure('Name', sprintf('Hist_Chg_%s_3sigma', rackName), 'NumberTitle', 'off', 'Visible', 'on');
        hold on;
        
        % Colors for years - Charge (Blue tones)
        yearColors = {'#0073C2', '#FF8C00'};  % Blue for 2021, Orange for 2022
        
        % Plot each year's data with year-specific min-max/30
        legend_handles = [];
        for yi = 1:length(yearList)
            year_vals = chg_vals_by_year{yi};
            if ~isempty(year_vals)
                % Apply 3σ outlier removal for this year
                pd = fitdist(year_vals', 'Normal');
                out_max = pd.mu + 3*pd.sigma;
                out_min = pd.mu - 3*pd.sigma;
                
                r = 1;
                while r <= length(year_vals)
                    if(year_vals(r) > out_max || year_vals(r) < out_min)
                        year_vals(r) = [];
                    else
                        r = r + 1;
                    end
                end
                
                if ~isempty(year_vals)
                    % Year-specific min-max/30 binning
                    min_val = min(year_vals);
                    max_val = max(year_vals);
                    bin_width = (max_val - min_val) / 30;
                    edges = min_val:bin_width:max_val;
                    
                    h = histfit(year_vals, 30, 'normal');  % Use 30 bins instead of edges
                    h(1).EdgeColor = yearColors{yi};
                    h(1).FaceColor = yearColors{yi};
                    h(1).FaceAlpha = 0.7;
                    h(2).Color = yearColors{yi};
                    h(2).LineWidth = 3;
                    
                    % Store only histogram bars for legend
                    legend_handles = [legend_handles, h(1)];
                end
            end
        end
        
        xlabel('R_{CHG} [m\Omega]'); ylabel('Frequency [-]'); 
        title(sprintf('Charge Peaks - %s (3σ outliers removed)', rackName));
        legend(legend_handles, yearList, 'Location', 'best');
        set(findall(gcf,'-property','FontSize'),'FontSize',20);
        set(findall(gcf,'-property','interpreter'),'interpreter','tex');
        set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex');
        saveas(gcf, fullfile(saveDir, sprintf('Hist_Chg_%s_3sigma.fig', rackName)));
    end
    
    if ~isempty(dch_vals_cleaned)
        figure('Name', sprintf('Hist_Dch_%s_3sigma', rackName), 'NumberTitle', 'off', 'Visible', 'on');
        hold on;
        
        % Colors for years - Discharge (Same as Charge)
        yearColors = {'#0073C2', '#FF8C00'};  % Blue for 2021, Orange for 2022
        
        % Plot each year's data with year-specific min-max/30
        legend_handles = [];
        for yi = 1:length(yearList)
            year_vals = dch_vals_by_year{yi};
            if ~isempty(year_vals)
                % Apply 3σ outlier removal for this year
                pd = fitdist(year_vals', 'Normal');
                out_max = pd.mu + 3*pd.sigma;
                out_min = pd.mu - 3*pd.sigma;
                
                r = 1;
                while r <= length(year_vals)
                    if(year_vals(r) > out_max || year_vals(r) < out_min)
                        year_vals(r) = [];
                    else
                        r = r + 1;
                    end
                end
                
                if ~isempty(year_vals)
                    % Year-specific min-max/30 binning
                    min_val = min(year_vals);
                    max_val = max(year_vals);
                    bin_width = (max_val - min_val) / 30;
                    edges = min_val:bin_width:max_val;
                    
                    h = histfit(year_vals, 30, 'normal');  % Use 30 bins instead of edges
                    h(1).EdgeColor = yearColors{yi};
                    h(1).FaceColor = yearColors{yi};
                    h(1).FaceAlpha = 0.7;
                    h(2).Color = yearColors{yi};
                    h(2).LineWidth = 3;
                    
                    % Store only histogram bars for legend
                    legend_handles = [legend_handles, h(1)];
                end
            end
        end
        
        xlabel('R_{DCH} [m\Omega]'); ylabel('Frequency [-]'); 
        title(sprintf('Discharge Peaks - %s (3σ outliers removed)', rackName));
        legend(legend_handles, yearList, 'Location', 'best');
        set(findall(gcf,'-property','FontSize'),'FontSize',20);
        set(findall(gcf,'-property','interpreter'),'interpreter','tex');
        set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex');
        saveas(gcf, fullfile(saveDir, sprintf('Hist_Dch_%s_3sigma.fig', rackName)));
    end
end


