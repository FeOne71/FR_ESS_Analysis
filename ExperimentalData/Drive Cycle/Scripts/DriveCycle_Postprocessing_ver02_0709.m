%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drive Cycle Data Parser (Modified)
% Parse 8 channels of real-load profile data by 3 SOC levels
% Remove charge/discharge data after final rest period
% File name: DriveCycle_Postprocessing_ver02_0709.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% Path configuration - use current directory
folderPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\csv';
saveFolder = fullfile(folderPath, 'parsed_data');
if ~exist(saveFolder, 'dir'); mkdir(saveFolder); end

% File list (8 channels)
fileNames = {
    'Ch9_Drive_200cyc.csv';
    'Ch10_Drive_200cyc.csv';
    'Ch11_Drive_200cyc.csv';
    'Ch12_Drive_200cyc.csv';
    'Ch13_Drive_200cyc.csv';
    'Ch14_Drive_200cyc.csv';
    'Ch15_Drive_200cyc.csv';
    'Ch16_Drive_200cyc.csv'
};

% Step Index definition by SOC
SOC90_stepIndex = [5, 7, 9, 11, 13, 15, 17, 19];
SOC70_stepIndex = [23, 25, 27, 29, 31, 33, 35, 37];
SOC50_stepIndex = [41, 43, 45, 47, 49, 51, 53, 55];

% Real-load profile names (8 profiles)
profileNames = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

% Rest time configuration (seconds)
initialRestTime = 8 * 60;  % Initial 8-minute rest period
finalRestTime = 8 * 60;    % Final 8-minute rest period (minimum)

% Initialize data structure for each channel
parsedDriveCycle_200cyc = struct();

fprintf('Starting real-load profile data parsing...\n');

for i = 1:length(fileNames)
    filename = fileNames{i};
    filepath = fullfile(folderPath, filename);
    
    if exist(filepath, 'file') ~= 2
        fprintf('File not found: %s\n', filepath);
        continue;
    end
    
    fprintf('Processing: %s\n', filename);
    
    % Read CSV file
    T = readtable(filepath);
    
    % Extract data
    stepIndex = T{:,2};     % Step Index
    stepType = T{:,3};      % Step Type
    time = T{:,5};          % Time [s]
    totalTime = T{:,6};     % Total Time [s]
    current = T{:,7};       % Current [A]
    voltage = T{:,8};       % Voltage [V]
    
    % Extract channel name
    [~, baseName, ~] = fileparts(filename);
    channelName = extractBetween(baseName, 'Ch', '_');
    channelFieldName = sprintf('ch%s_Drive_200cyc', channelName{1});
    
    % Initialize structure for each channel
    parsedDriveCycle_200cyc.(channelFieldName) = struct();
    
    % Extract and process SOC90 data
    fprintf('  Extracting SOC90 data...\n');
    for j = 1:length(SOC90_stepIndex)
        stepIdx = SOC90_stepIndex(j);
        profileName = profileNames{j};
        
        % Find data with matching Step Index and Step Type="SIM"
        mask = (stepIndex == stepIdx) & strcmp(stepType, 'SIM');
        
        if any(mask)
            % Extract original data
            step_voltage = voltage(mask);
            step_current = current(mask);
            step_time = time(mask);
            step_totalTime = totalTime(mask);
            
            % Remove data after final rest period
            [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
                removeFinalRestData(step_voltage, step_current, step_time, step_totalTime, finalRestTime);
            
            % Store in structure
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).V = filtered_V;
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).I = filtered_I;
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).t = filtered_t;
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).totalTime = filtered_totalTime;
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).stepIndex = stepIdx;
            
            fprintf('    %s: %d -> %d data points (removed after final rest)\n', ...
                    profileName, length(step_voltage), length(filtered_V));
        else
            fprintf('    %s: No data found (Step Index %d)\n', profileName, stepIdx);
        end
    end
    
    % Extract and process SOC70 data
    fprintf('  Extracting SOC70 data...\n');
    for j = 1:length(SOC70_stepIndex)
        stepIdx = SOC70_stepIndex(j);
        profileName = profileNames{j};
        
        % Find data with matching Step Index and Step Type="SIM"
        mask = (stepIndex == stepIdx) & strcmp(stepType, 'SIM');
        
        if any(mask)
            % Extract original data
            step_voltage = voltage(mask);
            step_current = current(mask);
            step_time = time(mask);
            step_totalTime = totalTime(mask);
            
            % Remove data after final rest period
            [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
                removeFinalRestData(step_voltage, step_current, step_time, step_totalTime, finalRestTime);
            
            % Store in structure
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).V = filtered_V;
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).I = filtered_I;
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).t = filtered_t;
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).totalTime = filtered_totalTime;
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).stepIndex = stepIdx;
            
            fprintf('    %s: %d -> %d data points (removed after final rest)\n', ...
                    profileName, length(step_voltage), length(filtered_V));
        else
            fprintf('    %s: No data found (Step Index %d)\n', profileName, stepIdx);
        end
    end
    
    % Extract and process SOC50 data
    fprintf('  Extracting SOC50 data...\n');
    for j = 1:length(SOC50_stepIndex)
        stepIdx = SOC50_stepIndex(j);
        profileName = profileNames{j};
        
        % Find data with matching Step Index and Step Type="SIM"
        mask = (stepIndex == stepIdx) & strcmp(stepType, 'SIM');
        
        if any(mask)
            % Extract original data
            step_voltage = voltage(mask);
            step_current = current(mask);
            step_time = time(mask);
            step_totalTime = totalTime(mask);
            
            % Remove data after final rest period
            [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
                removeFinalRestData(step_voltage, step_current, step_time, step_totalTime, finalRestTime);
            
            % Store in structure
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).V = filtered_V;
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).I = filtered_I;
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).t = filtered_t;
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).totalTime = filtered_totalTime;
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).stepIndex = stepIdx;
            
            fprintf('    %s: %d -> %d data points (removed after final rest)\n', ...
                    profileName, length(step_voltage), length(filtered_V));
        else
            fprintf('    %s: No data found (Step Index %d)\n', profileName, stepIdx);
        end
    end
    
    fprintf('  %s completed\n\n', channelName{1});
end

% Save results
savePath = fullfile(saveFolder, 'parsedDriveCycle_200cyc_filtered.mat');
save(savePath, 'parsedDriveCycle_200cyc');

fprintf('Parsing completed! Results saved to: %s\n', savePath);

% Print structure summary
fprintf('\n=== Parsing Results Summary ===\n');
channels = fieldnames(parsedDriveCycle_200cyc);
for i = 1:length(channels)
    channelName = channels{i};
    fprintf('Channel: %s\n', channelName);
    
    % Check number of data for each SOC
    if isfield(parsedDriveCycle_200cyc.(channelName), 'SOC90')
        soc90_fields = fieldnames(parsedDriveCycle_200cyc.(channelName).SOC90);
        fprintf('  SOC90: %d profiles\n', length(soc90_fields));
    end
    
    if isfield(parsedDriveCycle_200cyc.(channelName), 'SOC70')
        soc70_fields = fieldnames(parsedDriveCycle_200cyc.(channelName).SOC70);
        fprintf('  SOC70: %d profiles\n', length(soc70_fields));
    end
    
    if isfield(parsedDriveCycle_200cyc.(channelName), 'SOC50')
        soc50_fields = fieldnames(parsedDriveCycle_200cyc.(channelName).SOC50);
        fprintf('  SOC50: %d profiles\n', length(soc50_fields));
    end
end

fprintf('\nTotal parsed data: %d channels × 3 SOCs × 8 profiles = %d datasets\n', ...
        length(channels), length(channels) * 3 * 8);

%% Sub-function: Remove data after final rest period
function [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
    removeFinalRestData(voltage, current, time, totalTime, finalRestTime)
    
    % Convert time to seconds and normalize to start from 0
    if isduration(time)
        time_sec = seconds(time);
    else
        time_sec = time;
    end
    time_normalized = time_sec - time_sec(1);
    
    % Find points where current is below 2A (rest period)
    current_threshold = 2;  % Consider current below 2A as rest
    rest_mask = abs(current) < current_threshold;
    
    % Total data length
    totalLength = length(current);
    
    % Debug: Print basic information
    fprintf('      [DEBUG] Total data length: %d, Last current: %.2fA\n', totalLength, current(end));
    
    % Find all rest periods lasting 8 minutes or more in the entire dataset
    longRestPeriods = [];
    
    i = 1;
    while i <= totalLength
        if rest_mask(i)
            % Rest period start detected
            restStart = i;
            
            % Find the end of continuous rest period
            while i <= totalLength && rest_mask(i)
                i = i + 1;
            end
            restEnd = i - 1;
            
            % Calculate rest period duration
            restDuration = time_normalized(restEnd) - time_normalized(restStart);
            
            % Debug: Print found rest period information
            fprintf('      [DEBUG] Rest period found: %d~%d (%.0f sec, %.1f min)\n', ...
                    restStart, restEnd, restDuration, restDuration/60);
            
            % Save only rest periods lasting 8 minutes or more
            if restDuration >= finalRestTime
                longRestPeriods = [longRestPeriods; restStart, restEnd, restDuration];
                fprintf('      [DEBUG] → Selected as 8+ minute rest period\n');
            end
        else
            i = i + 1;
        end
    end
    
    % Check if there are rest periods lasting 8 minutes or more
    if ~isempty(longRestPeriods)
        % Select the last long rest period
        lastRestPeriod = longRestPeriods(end, :);
        restStartIdx = lastRestPeriod(1);
        restEndIdx = lastRestPeriod(2);
        restDuration = lastRestPeriod(3);
        
        % Save data only up to the end of final rest period
        cutoffIdx = restEndIdx;
        
        filtered_V = voltage(1:cutoffIdx);
        filtered_I = current(1:cutoffIdx);
        filtered_t = time(1:cutoffIdx);
        filtered_totalTime = totalTime(1:cutoffIdx);
        
        fprintf('      Final rest period found. data point starts from %d to %d (%.0f sec, %.1f min)\n', ...
                restStartIdx, restEndIdx, restDuration, restDuration/60);
        fprintf('      Data trimmed: %d -> %d points\n', totalLength, cutoffIdx);
        
    else
        % Keep original data if no rest periods lasting 8 minutes or more
        filtered_V = voltage;
        filtered_I = current;
        filtered_t = time;
        filtered_totalTime = totalTime;
        
        fprintf('      No rest periods lasting 8+ minutes found - keeping original data\n');
    end
end 