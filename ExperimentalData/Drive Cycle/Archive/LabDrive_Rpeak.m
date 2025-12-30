%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab Drive Cycle Peak Resistance Analysis
% OldData v1 logic applied to parsed Drive Cycle data
% Ch9 only, 0cyc/200cyc/400cyc comparison with histogram visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Parameters (OldData v1 logic)
Cnom = 64;                    % Cell nominal capacity [Ah] (Drive Cycle is cell level)
thr = [Cnom*0.05, Cnom*0.07];  % Idle threshold [charge, discharge] [A]
dI_chg = Cnom * 0.12;         % Charge current change threshold [A] 
dI_dch = Cnom * 0.2;         % Discharge current change threshold [A] 
ddI = 1;                      % Instantaneous current change threshold [A]
te = 10;                      % Moving average window (1 second = 10 samples)
dt = 10;                      % Peak duration (1 second = 10 samples)

fprintf('=== Lab Drive Cycle Peak Resistance Analysis ===\n');
fprintf('Parameters: Cnom=%.0fAh, thr=[%.1f, %.1f]A, dI_chg=%.1fA, dI_dch=%.1fA\n', ...
    Cnom, thr(1), thr(2), dI_chg, dI_dch);

%% Load data for all cycles - Auto-detect available cycles
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
allCycleData = struct();

fprintf('\n=== Loading Data ===\n');

% Auto-detect available cycle files
matFiles = dir(fullfile(dataDir, 'parsedDriveCycle_*cyc_filtered.mat'));
cycleTypes = {};

for i = 1:length(matFiles)
    filename = matFiles(i).name;
    % Extract cycle type from filename (e.g., '0cyc', '200cyc', '400cyc', '600cyc')
    match = regexp(filename, 'parsedDriveCycle_(\d+cyc)_filtered\.mat', 'tokens');
    if ~isempty(match)
        cycleTypes{end+1} = match{1}{1};
    end
end

% Sort cycle types numerically
cycleNumbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), cycleTypes);
[sortedNumbers, sortIdx] = sort(cycleNumbers);
cycleTypes = cycleTypes(sortIdx);

fprintf('Found cycles: %s\n', strjoin(cycleTypes, ', '));

% Load data for each detected cycle
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    dataFile = fullfile(dataDir, sprintf('parsedDriveCycle_%s_filtered.mat', cycleType));
    
    if exist(dataFile, 'file')
        fprintf('Loading: %s\n', dataFile);
        load(dataFile);
        eval(sprintf('allCycleData.(''cycle_%s'') = parsedDriveCycle_%s;', cycleType, cycleType));
    else
        fprintf('File not found: %s\n', dataFile);
    end
end

%% Initialize results structure
results = struct();
results.Ch9_Charging = struct();
results.Ch9_Discharging = struct();

% Initialize for each detected cycle
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    cycleFieldName = sprintf('cycle_%s', cycleType);
    results.Ch9_Charging.(cycleFieldName) = struct();
    results.Ch9_Discharging.(cycleFieldName) = struct();
    
    % Initialize for each SOC level
    socLevels = {'SOC90', 'SOC70', 'SOC50'};
    for socIdx = 1:length(socLevels)
        socLevel = socLevels{socIdx};
        results.Ch9_Charging.(cycleFieldName).(socLevel) = struct();
        results.Ch9_Discharging.(cycleFieldName).(socLevel) = struct();
        
        % Initialize for each profile (DC1-DC8)
        for profIdx = 1:8
            profileName = sprintf('DC%d', profIdx);
            results.Ch9_Charging.(cycleFieldName).(socLevel).(profileName) = [];
            results.Ch9_Discharging.(cycleFieldName).(socLevel).(profileName) = [];
        end
    end
end

%% Process Ch9 data for all cycles
fprintf('\n=== Processing Ch9 Data ===\n');

for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    
    cycleFieldName = sprintf('cycle_%s', cycleType);
    if ~isfield(allCycleData, cycleFieldName)
        fprintf('Skipping %s: No data loaded\n', cycleType);
        continue;
    end
    
    fprintf('\n--- Processing %s ---\n', cycleType);
    cycleData = allCycleData.(cycleFieldName);
    
    % Find Ch9 data
    ch9FieldName = '';
    fields = fieldnames(cycleData);
    for i = 1:length(fields)
        if contains(fields{i}, 'ch9')
            ch9FieldName = fields{i};
            break;
        end
    end
    
    if isempty(ch9FieldName)
        fprintf('Ch9 data not found in %s\n', cycleType);
        continue;
    end
    
    fprintf('Found Ch9 field: %s\n', ch9FieldName);
    ch9Data = cycleData.(ch9FieldName);
    
    % Process each SOC level
    socLevels = fieldnames(ch9Data);
    for socIdx = 1:length(socLevels)
        socLevel = socLevels{socIdx};
        fprintf('  Processing %s...\n', socLevel);
        
        socData = ch9Data.(socLevel);
        profiles = fieldnames(socData);
        
        for profIdx = 1:length(profiles)
            profileName = profiles{profIdx};
            profileData = socData.(profileName);
            
            % Extract time series data
            V = profileData.V;
            I = profileData.I;
            t = profileData.t;
            
            % Convert time to seconds if needed
            if isa(t, 'duration')
                t_seconds = seconds(t);
            else
                t_seconds = t;
            end
            
            % Normalize time to start from 0
            t_normalized = t_seconds - t_seconds(1);
            
            % Apply OldData v1 peak detection logic - EXACT COPY from fig4_5
            [PeakTime_Chg, PeakCurrent_Chg, PeakVoltage_Chg, PeakTime_Dchg, PeakCurrent_Dchg, PeakVoltage_Dchg] = ...
                detectPeaks_OldData_v1(V, I, t_normalized, thr, dI_chg, dI_dch, ddI, te, dt);
            
            % Calculate DCIR for charging peaks - EXACT COPY from fig4_5
            if ~isempty(PeakVoltage_Chg)
                R_Chg = calculateDCIR_peaks(PeakVoltage_Chg, PeakCurrent_Chg);
                results.Ch9_Charging.(cycleFieldName).(socLevel).(profileName) = ...
                    [results.Ch9_Charging.(cycleFieldName).(socLevel).(profileName), R_Chg];
            end
            
            % Calculate DCIR for discharging peaks - EXACT COPY from fig4_5
            if ~isempty(PeakVoltage_Dchg)
                R_Dchg = calculateDCIR_peaks(PeakVoltage_Dchg, PeakCurrent_Dchg);
                results.Ch9_Discharging.(cycleFieldName).(socLevel).(profileName) = ...
                    [results.Ch9_Discharging.(cycleFieldName).(socLevel).(profileName), R_Dchg];
            end
            
            % Count peaks
            chg_count = length(PeakVoltage_Chg);
            dchg_count = length(PeakVoltage_Dchg);
            fprintf('    %s: %d charging, %d discharging peaks\n', ...
                profileName, chg_count, dchg_count);
        end
    end
end

%% Generate Statistics
fprintf('\n=== Statistics Summary ===\n');
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    cycleFieldName = sprintf('cycle_%s', cycleType);
    
    total_chg_peaks = 0;
    total_dchg_peaks = 0;
    
    % Sum across all SOC levels and profiles
    socLevels = {'SOC90', 'SOC70', 'SOC50'};
    for socIdx = 1:length(socLevels)
        socLevel = socLevels{socIdx};
        for profIdx = 1:8
            profileName = sprintf('DC%d', profIdx);
            chg_peaks = length(results.Ch9_Charging.(cycleFieldName).(socLevel).(profileName));
            dchg_peaks = length(results.Ch9_Discharging.(cycleFieldName).(socLevel).(profileName));
            total_chg_peaks = total_chg_peaks + chg_peaks;
            total_dchg_peaks = total_dchg_peaks + dchg_peaks;
        end
    end
    
    fprintf('%s: %d charging peaks, %d discharging peaks\n', ...
        cycleType, total_chg_peaks, total_dchg_peaks);
end

%% Create Visualization
fprintf('\n=== Creating Visualizations ===\n');

% Figure 1: Charging DCIR by Peak Order
fprintf('\n=== Creating Chg Rpeak by Peak Order ===\n');
figure('Name', 'Ch9 Chg R_{Peak} by Peak Order', 'Position', [100, 100, 1200, 800]);
colors = {'b', 'r', 'g'}; % Blue for 0cyc, Red for 200cyc, Green for 400cyc

for profIdx = 1:8
    profileName = sprintf('DC%d', profIdx);
    subplot(2, 4, profIdx);
    
    hold on;
    for cycleIdx = 1:length(cycleTypes)
        cycleType = cycleTypes{cycleIdx};
        cycleFieldName = sprintf('cycle_%s', cycleType);
        
        % Get DCIR values for this profile and cycle (SOC50 only)
        dcir_data = results.Ch9_Charging.(cycleFieldName).SOC50.(profileName);
        
        if ~isempty(dcir_data)
            peak_order = 1:length(dcir_data);
            plot(peak_order, dcir_data, 'o-', 'Color', colors{cycleIdx}, ...
                'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', cycleType);
        end
    end
    
    title(sprintf('%s Charging', profileName));
    xlabel('Peak Order');
    ylabel('DCIR [mΩ]');
    legend('Location', 'best');
    grid on;
    hold off;
end

sgtitle('Ch9 Chg Rpeak by Peak Order (SOC50)');

% Figure 2: Discharging DCIR by Peak Order
fprintf('\n=== Creating Dchg Rpeak by Peak Order ===\n');
figure('Name', 'Ch9 Dchg R_{Peak} by Peak Order', 'Position', [150, 150, 1200, 800]);

for profIdx = 1:8
    profileName = sprintf('DC%d', profIdx);
    subplot(2, 4, profIdx);
    
    hold on;
    for cycleIdx = 1:length(cycleTypes)
        cycleType = cycleTypes{cycleIdx};
        cycleFieldName = sprintf('cycle_%s', cycleType);
        
        % Get DCIR values for this profile and cycle (SOC50 only)
        dcir_data = results.Ch9_Discharging.(cycleFieldName).SOC50.(profileName);
        
        if ~isempty(dcir_data)
            peak_order = 1:length(dcir_data);
            plot(peak_order, dcir_data, 'o-', 'Color', colors{cycleIdx}, ...
                'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', cycleType);
        end
    end
    
    title(sprintf('%s Discharging', profileName));
    xlabel('Peak Order');
    ylabel('DCIR [mΩ]');
    legend('Location', 'best');
    grid on;
    hold off;
end

sgtitle('Ch9 Dchg Rpeak by Peak Order (SOC50)');

% Figure 3: Peak Detection Visualization
fprintf('\n=== Creating Peak Detection Visualization ===\n');
figure('Name', 'Ch9 Peak Detection Overlay', 'Position', [200, 200, 1600, 800]);

% Process one cycle for visualization (use 0cyc as example)
cycleType = '0cyc';
cycleFieldName = sprintf('cycle_%s', cycleType);

if isfield(allCycleData, cycleFieldName)
    cycleData = allCycleData.(cycleFieldName);
    
    % Find Ch9 data
    ch9FieldName = '';
    fields = fieldnames(cycleData);
    for i = 1:length(fields)
        if contains(fields{i}, 'ch9')
            ch9FieldName = fields{i};
            break;
        end
    end
    
    if ~isempty(ch9FieldName)
        ch9Data = cycleData.(ch9FieldName);
        
        % Process only SOC50
        if isfield(ch9Data, 'SOC50')
            socData = ch9Data.SOC50;
            profiles = fieldnames(socData);
            
            for profIdx = 1:length(profiles)
                profileName = profiles{profIdx};
                profileData = socData.(profileName);
                
                % Extract time series data
                V = profileData.V;
                I = profileData.I;
                t = profileData.t;
                
                % Convert time to seconds if needed
                if isa(t, 'duration')
                    t_seconds = seconds(t);
                else
                    t_seconds = t;
                end
                
                % Normalize time to start from 0
                t_normalized = t_seconds - t_seconds(1);
                
                % Apply peak detection - EXACT COPY from fig4_5
                [PeakTime_Chg, PeakCurrent_Chg, PeakVoltage_Chg, PeakTime_Dchg, PeakCurrent_Dchg, PeakVoltage_Dchg] = ...
                    detectPeaks_OldData_v1(V, I, t_normalized, thr, dI_chg, dI_dch, ddI, te, dt);
                
                % Create subplot
                subplot(2, 4, profIdx);
                
                % Plot current vs time
                plot(t_normalized, I, 'Color', [0.69, 0.69, 0.69], 'LineWidth', 1, 'DisplayName', 'Current');
                hold on;
                
                % Overlay charging peaks - EXACT COPY from fig4_5
                chg_plotted = 0;
                if ~isempty(PeakVoltage_Chg)
                    for i = 1:length(PeakVoltage_Chg)
                        peak_time = PeakTime_Chg{i};
                        peak_current = PeakCurrent_Chg{i};
                        plot(peak_time, peak_current, 'r-', 'LineWidth', 3, 'DisplayName', 'Charging Peak');
                        chg_plotted = chg_plotted + 1;
                    end
                end
                
                % Overlay discharging peaks - EXACT COPY from fig4_5
                dchg_plotted = 0;
                if ~isempty(PeakVoltage_Dchg)
                    for i = 1:length(PeakVoltage_Dchg)
                        peak_time = PeakTime_Dchg{i};
                        peak_current = PeakCurrent_Dchg{i};
                        plot(peak_time, peak_current, 'b-', 'LineWidth', 3, 'DisplayName', 'Discharging Peak');
                        dchg_plotted = dchg_plotted + 1;
                    end
                end
                
                % Formatting
                xlabel('Time [s]');
                ylabel('Current [A]');
                grid on;
                
                % Add peak count to title (actual plotted peaks)
                title(sprintf('%s - Chg:%d, Dchg:%d', profileName, chg_plotted, dchg_plotted));
                                
                hold off;
            end
        end
    end
end

sgtitle(sprintf('Ch9 Peak Detection Overlay (%s) - Red: Charging, Blue: Discharging', cycleType));

%% Save results
savePath = 'LabDrive_Rpeak_Results.mat';
save(savePath, 'results');
fprintf('\nResults saved to: %s\n', savePath);

fprintf('\n=== Analysis Complete ===\n');

%% Sub-functions

function [PeakTime_Chg, PeakCurrent_Chg, PeakVoltage_Chg, PeakTime_Dchg, PeakCurrent_Dchg, PeakVoltage_Dchg] = detectPeaks_OldData_v1(V, I, t, thr, dI_chg, dI_dch, ddI, te, dt)
    % EXACT COPY of fig4_5 logic - Acceleration and Braking peak detection
    % Parameters: thr=[3.2, 3.2]A, dI_chg=9.6A, dI_dch=9.6A, ddI=1A, te=10, dt=10
    
    % Initialize peak structures - EXACT COPY from fig4_5
    clear PeakTime PeakVoltage PeakCurrent PeakT PeakEpoch R TimeR TimeEpoch PeakVoltageIn PeakCurrentIn PeakCurrentMax Capacity;
    
    % Filter on Current - EXACT COPY from fig4_5
    % Moving average of 1s --> 10 sampling points (our data: 0.1s sampling)
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
    
    % Derivative of Current - EXACT COPY from fig4_5
    dI_dt = zeros(length(I)-1, 1);
    for i = 1:length(I)-1
        dI_dt(i) = (I(i+1) - I(i)) / (t(i+1) - t(i));
    end
    
    % MA on Current Derivative - EXACT COPY from fig4_5
    N = length(dI_dt);
    filt_dI_dt = zeros(N,1);
    Pre = 0*ones(te/2,1);
    Post = 0*ones(te/2,1);
    calc_dI_dt  = [Pre;dI_dt;Post];
    for i = 1:N
        for m = 0:te
            filt_dI_dt(i) = filt_dI_dt(i) + calc_dI_dt(i+m);
        end
        filt_dI_dt(i) = filt_dI_dt(i)/(te+1);
    end
    
    % CHARGING Peak Detection - Using BRAKING logic from fig4_5 (lines 156-197)
    % In fig4_5: Braking = Charging, but current is negative
    % In our data: Charging = positive current, so we invert the current conditions
    clear PeakTime PeakVoltage PeakCurrent;
    PeakTime = {}; PeakVoltage = {}; PeakCurrent = {};  % Initialize empty cell arrays
    z = 1;
    for i = (te+1):(length(I)-dt-te)
        if (I(i+dt)-I(i)) > dI_chg  % Positive change for charging (using raw current)
            if abs(I(i)) < thr(1)
                % The variable ddI is introduced to ensure that a peak has an initial increase 
                % of at least 1A (peaks with initial flat plateaus are discarded)                
                if (I(i+1)-I(i)) > ddI
                    % The variable flag is introduced to check if the peak is always increasing and the current never changes sign 
                    flag = 1;
                    for zi = 1:dt
                        if filt_dI_dt(i+zi-1) < 0 || I(i+zi) < 0
                            flag = 0;
                        end
                    end
                    if flag == 1
                        if z == 1
                            % Debug: Print values for detected charging peaks
                            fprintf('DEBUG CHARGING PEAK DETECTED: Peak#%d, i=%d, t=%.1f, I_start=%.3f, I_end=%.3f, I_filt_change=%.3f, ddI_check=%.3f\n', ...
                                z, i, t(i), I(i), I(i+dt), I_filt(i+dt)-I_filt(i), I_filt(i+1)-I_filt(i));
                            for j = 1:dt
                                PeakTime{z}(j) = t(i+j-1);
                                PeakCurrent{z}(j) = I(i+j-1);
                                PeakVoltage{z}(j) = V(i+j-1);
                            end
                            z = z + 1;
                        else
                            if (t(i) - PeakTime{z-1}(1)) > dt
                                % Debug: Print values for detected charging peaks
                                fprintf('DEBUG CHARGING PEAK DETECTED: Peak#%d, i=%d, t=%.1f, I_start=%.3f, I_end=%.3f, I_filt_change=%.3f, ddI_check=%.3f\n', ...
                                    z, i, t(i), I(i), I(i+dt), I_filt(i+dt)-I_filt(i), I_filt(i+1)-I_filt(i));
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
    
    % Store charging peaks
    PeakTime_Chg = PeakTime;
    PeakCurrent_Chg = PeakCurrent;
    PeakVoltage_Chg = PeakVoltage;
    
    % DISCHARGING Peak Detection - Using ACCELERATION logic from fig4_5 (lines 70-114)
    % In fig4_5: Acceleration = Discharging, but current is positive
    % In our data: Discharging = negative current, so we invert the current conditions
    clear PeakTime PeakVoltage PeakCurrent PeakT PeakEpoch R TimeR TimeEpoch PeakVoltageIn PeakCurrentIn PeakCurrentMax Capacity;
    PeakTime = {}; PeakVoltage = {}; PeakCurrent = {};  % Initialize empty cell arrays
    z = 1;
    for i = (te+1):(length(I_filt)-dt-te)
        if (I(i+dt)-I(i)) < -dI_dch  % Negative change for discharging (using raw current)
            if I(i) > -thr(2) && I(i) < thr(2)
                % The variable ddI is introduced to ensure that a peak has an initial decrease 
                % of at least -1A (peaks with initial flat plateaus are discarded)                
                if (I(i+1)-I(i)) < -ddI
                    % The variable flag is introduced to check if the peak is always decreasing and the current never changes sign 
                    flag = 1;
                    for zi = 1:dt
                        if filt_dI_dt(i+zi-1) > 0 || I(i+zi) > 0
                            flag = 0;
                        end
                    end
                    if flag == 1
                        if z == 1
                            % Debug: Print values for detected discharging peaks
                            fprintf('DEBUG DISCHARGING PEAK DETECTED: Peak#%d, i=%d, t=%.1f, I_start=%.3f, I_end=%.3f, I_filt_change=%.3f, ddI_check=%.3f\n', ...
                                z, i, t(i), I(i), I(i+dt), I_filt(i+dt)-I_filt(i), I_filt(i+1)-I_filt(i));
                            % The current, voltage, time signals are saved for dt time steps once a peak has been detected
                            for j = 1:dt
                                PeakTime{z}(j) = t(i+j-1);
                                PeakCurrent{z}(j) = I(i+j-1);
                                PeakVoltage{z}(j) = V(i+j-1);
                            end
                            z = z + 1;
                        else
                            % Additional condition that discards a new peak if a sufficient amount of time (dt=1s) is not passed from
                            % the beginning of the last peak (in this way the same peak cannot be saved twice in case of numerical errors) 
                            if (t(i) - PeakTime{z-1}(1)) > dt
                                % Debug: Print values for detected discharging peaks
                                fprintf('DEBUG DISCHARGING PEAK DETECTED: Peak#%d, i=%d, t=%.1f, I_start=%.3f, I_end=%.3f, I_filt_change=%.3f, ddI_check=%.3f\n', ...
                                    z, i, t(i), I(i), I(i+dt), I_filt(i+dt)-I_filt(i), I_filt(i+1)-I_filt(i));
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
    
    % Store discharging peaks
    PeakTime_Dchg = PeakTime;
    PeakCurrent_Dchg = PeakCurrent;
    PeakVoltage_Dchg = PeakVoltage;
end

function R = calculateDCIR_peaks(PeakVoltage, PeakCurrent)
    % EXACT COPY of fig4_5 resistance computation (lines 117-121)
    % for i = 1:length(PeakTime)
    %     DV = -(PeakVoltage{i}(end)-PeakVoltage{i}(1));
    %     DI = PeakCurrent{i}(end)-PeakCurrent{i}(1);
    %     R(i) = DV/DI*1000;
    % end
    
    R = [];
    for i = 1:length(PeakVoltage)
        DV = (PeakVoltage{i}(end) - PeakVoltage{i}(1));  % No negative sign for our data
        DI = PeakCurrent{i}(end) - PeakCurrent{i}(1);
        if DI ~= 0
            R(i) = (DV / DI) * 1000; % Convert to mΩ
        else
            R(i) = NaN;
        end
    end
end
