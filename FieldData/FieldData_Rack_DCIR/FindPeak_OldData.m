%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FindPeak_OldData.m
% Battery Rack DCIR Analysis using MATLAB findpeaks function
% - Uses findpeaks to detect current peaks during charge/discharge transitions
% - Calculates resistance from voltage change at current peaks
% - Creates year-wise histograms for charge and discharge peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

%% Parameters
% Peak detection parameters
Cnom = 128;                   % Rack nominal Capacity (Ah)
thr_idle = Cnom*0.02;         % Idle threshold (A) 
min_peak_height = Cnom*0.07;   % Minimum peak height (A)
min_peak_distance = 1;        % Minimum distance between peaks (samples)
min_peak_prominence = Cnom*0.1; % Minimum peak prominence (A) - PeakDetection.txt 스타일로 완화 

% Cell conversion parameters
Ns = 17*14;    % 238s (series cells)
Np = 2;        % 2p (parallel strings)

% Visualization parameters
Fontsize = 12;
LineWidth = 2;

%% Data paths
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR','FindPeak_Results');

% Create save directory if it doesn't exist
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Initialize data structures
yearList = {'Year2021'}; %, 'Year2023', 'Year2025'};
yearLabels = {'2021'}; %, '2023', '2025'};  % For display purposes
rackNames_all = {'Rack01'};

% Results storage
Results = struct();
for y = 1:length(yearList)
    Results.(yearList{y}).ChgPeaks = [];
    Results.(yearList{y}).DchPeaks = [];
    Results.(yearList{y}).ChgResistances = [];
    Results.(yearList{y}).DchResistances = [];
end

%% Process each year
for y = 1:length(yearList)
    year = yearList{y};
    yearLabel = yearLabels{y};
    fprintf('\n=== Processing Year: %s ===\n', yearLabel);
    
    % Set data path based on year
    if strcmp(year, 'Year2021')
        dataPath = fullfile(dataDir, '2021', '202106');
        dataFiles = dir(fullfile(dataPath, 'Raw_*.mat'));
    elseif strcmp(year, 'Year2022')
        dataPath = fullfile(dataDir, '2022', '202206');  % 2022년 경로로 수정
        dataFiles = dir(fullfile(dataPath, 'Raw_*.mat'));
    elseif strcmp(year, 'Year2023')
        dataPath = fullfile(dataDir, '2023', '202310');  % 2022년 경로로 수정
        dataFiles = dir(fullfile(dataPath, 'Raw_*.mat'));        
    elseif strcmp(year, 'Year2023')
        dataPath = fullfile(dataDir, 'New', '2023', '202310');
        dataFiles = dir(fullfile(dataPath, 'Raw_*.mat'));
    elseif strcmp(year, 'Year2024')
        dataPath = fullfile(dataDir, 'New', '2024', '202401');
        dataFiles = dir(fullfile(dataPath, 'Raw_*.mat'));
    elseif strcmp(year, 'Year2025')
        dataPath = fullfile(dataDir, 'New', '2025', '202501');
        dataFiles = dir(fullfile(dataPath, 'Raw_*.mat'));
    end
    
    if isempty(dataFiles)
        fprintf('No data files found for year %s\n', year);
        continue;
    end
    
    year_chg_peaks = [];
    year_dch_peaks = [];
    year_chg_resistances = [];
    year_dch_resistances = [];
    year_chg_peak_times = [];
    year_dch_peak_times = [];
    
    %% Process each file in the year
    for f = 1:length(dataFiles)
        dataFile = fullfile(dataPath, dataFiles(f).name);
        fprintf('Processing: %s\n', dataFiles(f).name);
        
        try
            % Load data
            S = load(dataFile);
            if isfield(S, 'Raw_Rack')
                Raw_Rack = S.Raw_Rack;
            else
                % Find the struct containing rack data
                fieldNames = fieldnames(S);
                for fn = 1:length(fieldNames)
                    if isstruct(S.(fieldNames{fn}))
                        Raw_Rack = S.(fieldNames{fn});
                        break;
                    end
                end
            end
            
            % Process each rack
            for r = 1:length(rackNames_all)
                rackName = rackNames_all{r};
                
                if ~isfield(Raw_Rack, rackName)
                    continue;
                end
                
                rackData = Raw_Rack.(rackName);
                
                % Extract signals
                if strcmp(year, 'Year2021')
                    t = rackData.Time;
                    I = rackData.DCCurrent_A;
                    V = rackData.AverageCV_V;
                elseif strcmp(year, 'Year2022')
                    % 2022 data format
                    if isfield(rackData, 'Date_Time')
                        if isduration(rackData.Date_Time)
                            t = datetime(2022,6,12) + rackData.Date_Time;
                        else
                            t = datetime(rackData.Date_Time);
                        end
                    else
                        % 2022년 Time 필드가 이미 datetime 형식
                        if isdatetime(rackData.Time)
                            t = rackData.Time;
                        else
                            t = datetime(rackData.Time);
                        end
                    end
                    I = rackData.DCCurrent_A;
                    V = rackData.AverageCV_V;
                else
                    % New data format (2023, 2024, 2025)
                    if isfield(rackData, 'Date_Time')
                        if isduration(rackData.Date_Time)
                            if strcmp(year, 'Year2023')
                                t = datetime(2023,10,16) + rackData.Date_Time;
                            elseif strcmp(year, 'Year2024')
                                t = datetime(2024,1,1) + rackData.Date_Time;
                            else  % Year2025
                                t = datetime(2025,1,1) + rackData.Date_Time;
                            end
                        else
                            t = datetime(rackData.Date_Time);
                        end
                    else
                        if strcmp(year, 'Year2023')
                            t = datetime(2023,10,16) + duration(string(rackData.Time),'InputFormat','hh:mm:ss');
                        elseif strcmp(year, 'Year2024')
                            t = datetime(2024,1,1) + duration(string(rackData.Time),'InputFormat','hh:mm:ss');
                        else  % Year2025
                            t = datetime(2025,1,1) + duration(string(rackData.Time),'InputFormat','hh:mm:ss');
                        end
                    end
                    I = rackData.DCCurrent;
                    V = rackData.CVavg;
                end
                
                % Convert time to seconds
                t = datetime(t);
                t_sec = seconds(t - t(1));
                
                % Convert to cell units
                I_cell = I / Np;  % A per cell
                
                %% Find current peaks using findpeaks
                
                % Charge peaks (positive current changes)
                [chg_peaks, chg_locs] = findpeaks(I_cell, ...
                    'MinPeakProminence', min_peak_prominence/Np);
                
                % Discharge peaks (negative current changes)
                [dch_peaks, dch_locs] = findpeaks(-I_cell, ...
                    'MinPeakProminence', min_peak_prominence/Np);
                dch_peaks = -dch_peaks;  % Convert back to negative
                
                %% Filter peaks based on idle condition
                % Check if current was in idle state before the peak
                valid_chg_peaks = [];
                valid_chg_locs = [];
                valid_dch_peaks = [];
                valid_dch_locs = [];
                
                % Check charge peaks
                for p = 1:length(chg_peaks)
                    peak_loc = chg_locs(p);
                    if peak_loc > 3  % Ensure we have enough history (완화: 5 -> 3)
                        % Check if current was in idle state 3 samples before (완화: 5 -> 3)
                        idle_before = all(abs(I_cell(peak_loc-3:peak_loc-1)) < thr_idle/Np);
                        if idle_before
                            valid_chg_peaks = [valid_chg_peaks, chg_peaks(p)];
                            valid_chg_locs = [valid_chg_locs, peak_loc];
                        end
                    end
                end
                
                % Check discharge peaks
                for p = 1:length(dch_peaks)
                    peak_loc = dch_locs(p);
                    if peak_loc > 3  % Ensure we have enough history (완화: 5 -> 3)
                        % Check if current was in idle state 3 samples before (완화: 5 -> 3)
                        idle_before = all(abs(I_cell(peak_loc-3:peak_loc-1)) < thr_idle/Np);
                        if idle_before
                            valid_dch_peaks = [valid_dch_peaks, dch_peaks(p)];
                            valid_dch_locs = [valid_dch_locs, peak_loc];
                        end
                    end
                end
                
                %% Calculate resistance from voltage change
                chg_resistances = [];
                dch_resistances = [];
                chg_peak_times = [];
                dch_peak_times = [];
                
                % Calculate resistance for charge peaks
                for p = 1:length(valid_chg_peaks)
                    peak_loc = valid_chg_locs(p);
                    if peak_loc > 5 && peak_loc < length(V) - 5  % 완화: 10 -> 5
                        % Voltage change from 3 samples before to 3 samples after (완화: 5 -> 3)
                        V_before = mean(V(peak_loc-3:peak_loc-1));
                        V_after = mean(V(peak_loc+1:peak_loc+3));
                        dV = V_after - V_before;
                        dI = valid_chg_peaks(p);
                        
                        if dI ~= 0
                            R = abs(dV / dI) * 1000;  % mOhm
                            chg_resistances = [chg_resistances, R];
                            % Store absolute time instead of relative time
                            chg_peak_times = [chg_peak_times, t(peak_loc)];
                        end
                    end
                end
                
                % Calculate resistance for discharge peaks
                for p = 1:length(valid_dch_peaks)
                    peak_loc = valid_dch_locs(p);
                    if peak_loc > 5 && peak_loc < length(V) - 5  % 완화: 10 -> 5
                        % Voltage change from 3 samples before to 3 samples after (완화: 5 -> 3)
                        V_before = mean(V(peak_loc-3:peak_loc-1));
                        V_after = mean(V(peak_loc+1:peak_loc+3));
                        dV = V_after - V_before;
                        dI = valid_dch_peaks(p);
                        
                        if dI ~= 0
                            R = abs(dV / dI) * 1000;  % mOhm
                            dch_resistances = [dch_resistances, R];
                            % Store absolute time instead of relative time
                            dch_peak_times = [dch_peak_times, t(peak_loc)];
                        end
                    end
                end
                
                % Store results
                year_chg_peaks = [year_chg_peaks, valid_chg_peaks];
                year_dch_peaks = [year_dch_peaks, valid_dch_peaks];
                year_chg_resistances = [year_chg_resistances, chg_resistances];
                year_dch_resistances = [year_dch_resistances, dch_resistances];
                year_chg_peak_times = [year_chg_peak_times, chg_peak_times];
                year_dch_peak_times = [year_dch_peak_times, dch_peak_times];
                
                fprintf('  %s: %d chg peaks, %d dch peaks\n', ...
                    rackName, length(valid_chg_peaks), length(valid_dch_peaks));
            end
            
        catch ME
            fprintf('Error processing %s: %s\n', dataFiles(f).name, ME.message);
            continue;
        end
    end
    
    % Store year results
    Results.(year).ChgPeaks = year_chg_peaks;
    Results.(year).DchPeaks = year_dch_peaks;
    Results.(year).ChgResistances = year_chg_resistances;
    Results.(year).DchResistances = year_dch_resistances;
    Results.(year).ChgPeakTimes = year_chg_peak_times;
    Results.(year).DchPeakTimes = year_dch_peak_times;
    
    fprintf('Year %s Summary: %d chg peaks, %d dch peaks\n', ...
        yearLabel, length(year_chg_peaks), length(year_dch_peaks));
end

%% Create histograms
figure('Name', 'Peak Resistance Histograms by Year', 'NumberTitle', 'off');
set(gcf, 'Position', [100, 100, 1200, 800]);

% Define colors for each year
yearColors = [
    0.0000 0.4470 0.7410;  % 2021 - blue
    0.8500 0.3250 0.0980;  % 2023 - orange
    0.4940 0.1840 0.5560   % 2025 - purple
];

% Charge resistance histograms
subplot(2,1,1);
hold on; grid on;
for y = 1:length(yearList)
    year = yearList{y};
    yearLabel = yearLabels{y};
    resistances = Results.(year).ChgResistances;
    
    if ~isempty(resistances)
        histogram(resistances, 'BinWidth', 0.1, 'FaceAlpha', 0.6, ...
            'FaceColor', yearColors(y,:), 'DisplayName', yearLabel);
    end
end
xlabel('Resistance [mΩ]');
ylabel('Count');
title('Charge Peak Resistances');
legend('Location', 'best');
set(gca, 'FontSize', Fontsize);

% Discharge resistance histograms
subplot(2,1,2);
hold on; grid on;
for y = 1:length(yearList)
    year = yearList{y};
    yearLabel = yearLabels{y};
    resistances = Results.(year).DchResistances;
    
    if ~isempty(resistances)
        histogram(resistances, 'BinWidth', 0.1, 'FaceAlpha', 0.6, ...
            'FaceColor', yearColors(y,:), 'DisplayName', yearLabel);
    end
end
xlabel('Resistance [mΩ]');
ylabel('Count');
title('Discharge Peak Resistances');
legend('Location', 'best');
set(gca, 'FontSize', Fontsize);

%% Create time-current plots with peak markers
fprintf('\n=== Creating Time-Current Plots ===\n');

% Check if Results structure exists and has data
if ~exist('Results', 'var')
    fprintf('ERROR: Results structure not found!\n');
    return;
end

figure('Name', 'Time-Current Plots with Peak Markers', 'NumberTitle', 'off');
set(gcf, 'Position', [100, 100, 1200, 800]);

% Create subplots for each year
numYears = length(yearList);
for y = 1:numYears
    year = yearList{y};
    yearLabel = yearLabels{y};
    
    fprintf('Processing year: %s\n', yearLabel);
    
    % Check if this year has data in Results
    if ~isfield(Results, year)
        fprintf('WARNING: No data found for year %s\n', year);
        continue;
    end
    
    subplot(numYears, 1, y);
    
    % Set data path based on year
    dataPath = '';
    if strcmp(year, 'Year2021')
        dataPath = fullfile(dataDir, '2021', '202106');
    elseif strcmp(year, 'Year2022')
        dataPath = fullfile(dataDir, '2022', '202206');
    elseif strcmp(year, 'Year2023')
        dataPath = fullfile(dataDir, 'New', '2023', '202310');
    elseif strcmp(year, 'Year2024')
        dataPath = fullfile(dataDir, 'New', '2024', '202401');
    elseif strcmp(year, 'Year2025')
        dataPath = fullfile(dataDir, 'New', '2025', '202501');
    end
    
    dataFiles = dir(fullfile(dataPath, 'Raw_*.mat'));
    fprintf('Data path: %s\n', dataPath);
    fprintf('Found %d data files\n', length(dataFiles));
    
    if ~isempty(dataFiles)
        % Find the file with the most peaks
        max_peaks = 0;
        best_file_idx = 1;
        
        for f = 1:length(dataFiles)
            dataFile = fullfile(dataPath, dataFiles(f).name);
            
            try
                S = load(dataFile);
                
                if isfield(S, 'Raw_Rack')
                    Raw_Rack = S.Raw_Rack;
                else
                    fieldNames = fieldnames(S);
                    for fn = 1:length(fieldNames)
                        if isstruct(S.(fieldNames{fn}))
                            Raw_Rack = S.(fieldNames{fn});
                            break;
                        end
                    end
                end
                
                if isfield(Raw_Rack, 'Rack01')
                    rackData = Raw_Rack.Rack01;
                    
                    % Extract time data to check peak count
                    if strcmp(year, 'Year2021')
                        t = rackData.Time;
                    elseif strcmp(year, 'Year2022')
                        if isfield(rackData, 'Date_Time')
                            if isduration(rackData.Date_Time)
                                t = datetime(2022,6,12) + rackData.Date_Time;
                            else
                                t = datetime(rackData.Date_Time);
                            end
                        else
                            if isdatetime(rackData.Time)
                                t = rackData.Time;
                            else
                                t = datetime(rackData.Time);
                            end
                        end
                    else
                        if isfield(rackData, 'Date_Time')
                            if isduration(rackData.Date_Time)
                                if strcmp(year, 'Year2023')
                                    t = datetime(2023,10,16) + rackData.Date_Time;
                                elseif strcmp(year, 'Year2024')
                                    t = datetime(2024,1,1) + rackData.Date_Time;
                                else
                                    t = datetime(2025,1,1) + rackData.Date_Time;
                                end
                            else
                                t = datetime(rackData.Date_Time);
                            end
                        else
                            if strcmp(year, 'Year2023')
                                t = datetime(2023,10,16) + duration(string(rackData.Time),'InputFormat','hh:mm:ss');
                            elseif strcmp(year, 'Year2024')
                                t = datetime(2024,1,1) + duration(string(rackData.Time),'InputFormat','hh:mm:ss');
                            else
                                t = datetime(2025,1,1) + duration(string(rackData.Time),'InputFormat','hh:mm:ss');
                            end
                        end
                    end
                    
                    % Count peaks in this file
                    file_start_time = t(1);
                    file_end_time = t(end);
                    
                    chg_peak_times_abs = Results.(year).ChgPeakTimes;
                    dch_peak_times_abs = Results.(year).DchPeakTimes;
                    
                    chg_mask = (chg_peak_times_abs >= file_start_time) & (chg_peak_times_abs <= file_end_time);
                    dch_mask = (dch_peak_times_abs >= file_start_time) & (dch_peak_times_abs <= file_end_time);
                    
                    total_peaks = sum(chg_mask) + sum(dch_mask);
                    
                    if total_peaks > max_peaks
                        max_peaks = total_peaks;
                        best_file_idx = f;
                    end
                end
            catch ME
                continue;
            end
        end
        
        % Plot only the file with the most peaks
        dataFile = fullfile(dataPath, dataFiles(best_file_idx).name);
        fprintf('Plotting file with most peaks: %s (%d peaks)\n', dataFiles(best_file_idx).name, max_peaks);
        
        try
            S = load(dataFile);
            
            if isfield(S, 'Raw_Rack')
                Raw_Rack = S.Raw_Rack;
            else
                fieldNames = fieldnames(S);
                for fn = 1:length(fieldNames)
                    if isstruct(S.(fieldNames{fn}))
                        Raw_Rack = S.(fieldNames{fn});
                        break;
                    end
                end
            end
            
            if isfield(Raw_Rack, 'Rack01')
                rackData = Raw_Rack.Rack01;
                
                % Extract time and current data
                if strcmp(year, 'Year2021')
                    t = rackData.Time;
                    I = rackData.DCCurrent_A;
                elseif strcmp(year, 'Year2022')
                    if isfield(rackData, 'Date_Time')
                        if isduration(rackData.Date_Time)
                            t = datetime(2022,6,12) + rackData.Date_Time;
                        else
                            t = datetime(rackData.Date_Time);
                        end
                    else
                        if isdatetime(rackData.Time)
                            t = rackData.Time;
                        else
                            t = datetime(rackData.Time);
                        end
                    end
                    I = rackData.DCCurrent_A;
                else
                    if isfield(rackData, 'Date_Time')
                        if isduration(rackData.Date_Time)
                            if strcmp(year, 'Year2023')
                                t = datetime(2023,10,16) + rackData.Date_Time;
                            elseif strcmp(year, 'Year2024')
                                t = datetime(2024,1,1) + rackData.Date_Time;
                            else
                                t = datetime(2025,1,1) + rackData.Date_Time;
                            end
                        else
                            t = datetime(rackData.Date_Time);
                        end
                    else
                        if strcmp(year, 'Year2023')
                            t = datetime(2023,10,16) + duration(string(rackData.Time),'InputFormat','hh:mm:ss');
                        elseif strcmp(year, 'Year2024')
                            t = datetime(2024,1,1) + duration(string(rackData.Time),'InputFormat','hh:mm:ss');
                        else
                            t = datetime(2025,1,1) + duration(string(rackData.Time),'InputFormat','hh:mm:ss');
                        end
                    end
                    I = rackData.DCCurrent;
                end
                
                % Convert to cell current
                I_cell = I / Np;
                
                % Convert time to seconds for plotting
                if isdatetime(t)
                    t_sec = seconds(t - t(1));
                else
                    t_sec = seconds(t);
                end
                
                % Plot current vs time
                plot(t_sec, I_cell, 'k-', 'LineWidth', 0.5);
                hold on; grid on;
                
                % Find peaks that belong to this specific file
                file_start_time = t(1);
                file_end_time = t(end);
                
                % Filter peaks that fall within this file's time range
                chg_peak_times_abs = Results.(year).ChgPeakTimes;
                chg_peak_currents = Results.(year).ChgPeaks;
                dch_peak_times_abs = Results.(year).DchPeakTimes;
                dch_peak_currents = Results.(year).DchPeaks;
                
                % Find peaks within this file's time range
                chg_mask = (chg_peak_times_abs >= file_start_time) & (chg_peak_times_abs <= file_end_time);
                dch_mask = (dch_peak_times_abs >= file_start_time) & (dch_peak_times_abs <= file_end_time);
                
                % Mark charge peaks for this file
                if any(chg_mask)
                    chg_peak_times_rel = seconds(chg_peak_times_abs(chg_mask) - t(1));
                    plot(chg_peak_times_rel, chg_peak_currents(chg_mask), 'ro', 'MarkerSize', 8, ...
                        'MarkerFaceColor', 'red', 'DisplayName', 'Charge Peaks');
                end
                
                % Mark discharge peaks for this file
                if any(dch_mask)
                    dch_peak_times_rel = seconds(dch_peak_times_abs(dch_mask) - t(1));
                    plot(dch_peak_times_rel, dch_peak_currents(dch_mask), 'bo', 'MarkerSize', 8, ...
                        'MarkerFaceColor', 'blue', 'DisplayName', 'Discharge Peaks');
                end
                
                xlabel('Time [s]');
                ylabel('Current [A]');
                title(sprintf('Year %s - %s (Most Peaks: %d)', yearLabel, dataFiles(best_file_idx).name, max_peaks));
                legend('Location', 'best');
                set(gca, 'FontSize', 12);
            end
        catch ME
            fprintf('Error loading best file %s: %s\n', dataFiles(best_file_idx).name, ME.message);
        end
    end
end

%% Print summary statistics
fprintf('\n=== Peak Resistance Summary ===\n');
fprintf('%-8s %-12s %-12s %-12s %-12s\n', 'Year', 'Chg_Count', 'Chg_Mean', 'Dch_Count', 'Dch_Mean');
fprintf('%s\n', repmat('-', 1, 60));

for y = 1:length(yearList)
    year = yearList{y};
    yearLabel = yearLabels{y};
    chg_res = Results.(year).ChgResistances;
    dch_res = Results.(year).DchResistances;
    
    chg_mean = mean(chg_res);
    dch_mean = mean(dch_res);
    
    fprintf('%-8s %-12d %-12.3f %-12d %-12.3f\n', ...
        yearLabel, length(chg_res), chg_mean, length(dch_res), dch_mean);
end

%% Save results
save(fullfile(saveDir, 'FindPeak_Results.mat'), 'Results', '-v7.3');
fprintf('\nResults saved to: %s\n', fullfile(saveDir, 'FindPeak_Results.mat'));
