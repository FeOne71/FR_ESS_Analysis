%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab Data - Drive Cycle DCIR Analysis (RPT0)
% Event detection and time-based DCIR calculation from real load profile data
% Logic adapted from FieldData_DCIR_Charge_CurrentClustering_Auto.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Load Data
load('G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\parsedDriveCycle\parsedDriveCycle_0cyc_filtered.mat');

% Check what variables are loaded
fprintf('=== Loaded variables ===\n');
vars = who;
for i = 1:length(vars)
    fprintf('%s\n', vars{i});
end

% Check if parsedDriveCycle_0cyc exists, if not try parsedDriveCycle_0cyc
if exist('parsedDriveCycle_0cyc', 'var')
    fprintf('\nUsing parsedDriveCycle_0cyc\n');
    data_var = parsedDriveCycle_0cyc;
elseif exist('parsedDriveCycle_0cyc', 'var')
    fprintf('\nUsing parsedDriveCycle_0cyc (0cyc file contains 0cyc data)\n');
    data_var = parsedDriveCycle_0cyc;
else
    error('No parsedDriveCycle variable found!');
end

% Check data structure
fprintf('\n=== Data structure check ===\n');
fprintf('Data variable class: %s\n', class(data_var));
fprintf('Data variable size: ');
disp(size(data_var));

if isstruct(data_var)
    fprintf('Data variable fields:\n');
    fields = fieldnames(data_var);
    for i = 1:length(fields)
        fprintf('  %s\n', fields{i});
    end
end

%% Settings (64Ah cell)
dt_list = [1, 3, 5, 10, 30, 50];
Cnom = 64;    % Battery capacity [Ah]
Pnom = 0.235; % Nominal Power    [kW] 3.68V * 64Ah = 235.52W
current_threshold = Cnom * 0.02;  % Idle current threshold [A]
min_duration      = 10;           % Minimum driving duration [s]
max_P_std         = Pnom * 0.03;  % Maximum power standard deviation [kW] (cell level)
max_I_std         = Cnom * 0.025; % Maximum current standard deviation [A]

fprintf('=== 0cyc DCIR Analysis ===\n');
fprintf('Current threshold: %.2f A\n', current_threshold);
fprintf('Min duration: %d s\n', min_duration);
fprintf('Max P std: %.3f kW\n', max_P_std);
fprintf('Max I std: %.2f A\n', max_I_std);

%% Pulse Detection & DCIR calculation
Lab_DC_DCIR_0cyc = struct();

channels = fieldnames(data_var);
fprintf('\nFound %d channels: ', length(channels));
for i = 1:length(channels)
    fprintf('%s ', channels{i});
end
fprintf('\n');

for ch_idx = 1:length(channels)
    channelName = channels{ch_idx};
    channel_data = data_var.(channelName);
    soc_levels = fieldnames(channel_data);
    
    fprintf('\n=== Processing %s ===\n', channelName);
    fprintf('SOC levels: ');
    for i = 1:length(soc_levels)
        fprintf('%s ', soc_levels{i});
    end
    fprintf('\n');
    
    chg_struct_name = [channelName '_ChgEvent'];
    dchg_struct_name = [channelName '_DchEvent'];
    if ~isfield(Lab_DC_DCIR_0cyc, chg_struct_name)
        Lab_DC_DCIR_0cyc.(chg_struct_name) = struct();
    end
    if ~isfield(Lab_DC_DCIR_0cyc, dchg_struct_name)
        Lab_DC_DCIR_0cyc.(dchg_struct_name) = struct();
    end
    
    for soc_idx = 1:length(soc_levels)
        soc_level = soc_levels{soc_idx};
        soc_data = channel_data.(soc_level);
        profiles = fieldnames(soc_data);
        
        fprintf('  %s: ', soc_level);
        for i = 1:length(profiles)
            fprintf('%s ', profiles{i});
        end
        fprintf('\n');
        
        if ~isfield(Lab_DC_DCIR_0cyc.(chg_struct_name), soc_level)
            Lab_DC_DCIR_0cyc.(chg_struct_name).(soc_level) = struct();
        end
        if ~isfield(Lab_DC_DCIR_0cyc.(dchg_struct_name), soc_level)
            Lab_DC_DCIR_0cyc.(dchg_struct_name).(soc_level) = struct();
        end
        
        for prof_idx = 1:length(profiles)
            profile_name = profiles{prof_idx};
            profile_data = soc_data.(profile_name);
            
            if ~isfield(Lab_DC_DCIR_0cyc.(chg_struct_name).(soc_level), profile_name)
                Lab_DC_DCIR_0cyc.(chg_struct_name).(soc_level).(profile_name) = struct();
            end
            if ~isfield(Lab_DC_DCIR_0cyc.(dchg_struct_name).(soc_level), profile_name)
                Lab_DC_DCIR_0cyc.(dchg_struct_name).(soc_level).(profile_name) = struct();
            end
            
            V = profile_data.V;
            I = profile_data.I;
            t = profile_data.t;
            totalTime = profile_data.totalTime;
            stepIndex = profile_data.stepIndex;
            
            % Debug: Print data info
            fprintf('    %s: %d points, I range [%.2f, %.2f] A, V range [%.2f, %.2f] V\n', ...
                profile_name, length(I), min(I), max(I), min(V), max(V));
            
            if isa(t, 'duration')
                t_seconds = seconds(t);
            else
                t_seconds = t;
            end
            
            %% Step 1: Idle -> Load transition detection
            is_idle = abs(I) < current_threshold;
            is_driving = abs(I) >= current_threshold;
            
            idle_to_driving = find(is_idle(1:end-1) & is_driving(2:end));
            
            fprintf('      Idle points: %d/%d (%.1f%%)\n', sum(is_idle), length(is_idle), sum(is_idle)/length(is_idle)*100);
            fprintf('      Driving points: %d/%d (%.1f%%)\n', sum(is_driving), length(is_driving), sum(is_driving)/length(is_driving)*100);
            fprintf('      Transitions found: %d\n', length(idle_to_driving));
            
            if isempty(idle_to_driving)
                fprintf('      → SKIP: No transitions detected\n');
                continue;  % Skip if no transition
            end
            
            %% Step 2: Event Detection (charging/discharging split)
            chg_event_count = 0;
            dchg_event_count = 0;
            for i = 1:length(idle_to_driving)
                idx1 = idle_to_driving(i);
                start_driving_idx = idx1 + 1;
                
                % Find end index of driving
                driving_end_idx = start_driving_idx;
                while driving_end_idx <= length(is_driving) && is_driving(driving_end_idx)
                    driving_end_idx = driving_end_idx + 1;
                end
                driving_end_idx = driving_end_idx - 1;
                
                start_idx = idx1;
                end_idx = driving_end_idx;
                
                % Check driving duration (actual time)
                driving_time = t_seconds(driving_end_idx) - t_seconds(start_driving_idx);
                fprintf('      Event %d: duration = %.1f s', i, driving_time);
                
                if driving_time < min_duration
                    fprintf(' → SKIP: Too short (< %d s)\n', min_duration);
                    continue;
                end
                
                % Extract segment
                t_seg = t_seconds(start_idx:end_idx);
                I_seg = I(start_idx:end_idx);
                V_seg = V(start_idx:end_idx);
                
                P_seg = V_seg .* I_seg / 1000;  % [kW]
                
                % idx2: end of driving segment
                idx2 = driving_end_idx - start_idx + 1;
                
                %% Step 3: Stability check
                power_std   = std(P_seg(3:idx2-2));
                current_std = std(I_seg(3:idx2-2));
                
                fprintf(', P_std=%.3f kW, I_std=%.2f A', power_std, current_std);
                
                if power_std >= max_P_std || current_std >= max_I_std
                    fprintf(' → FILTERED: Unstable\n');
                    continue;
                end                   
                
                % Charging/discharging split
                mean_current = mean(I_seg);
                if mean_current > 0
                    chg_event_count = chg_event_count + 1;
                    evtName = sprintf('event%d', chg_event_count);
                    target_struct = Lab_DC_DCIR_0cyc.(chg_struct_name);
                    fprintf(' → CHARGING (mean I=%.2f A)\n', mean_current);
                elseif mean_current < 0
                    dchg_event_count = dchg_event_count + 1;
                    evtName = sprintf('event%d', dchg_event_count);
                    target_struct = Lab_DC_DCIR_0cyc.(dchg_struct_name);
                    fprintf(' → DISCHARGING (mean I=%.2f A)\n', mean_current);
                else
                    fprintf(' → SKIP: Zero current\n');
                    continue;
                end
                
                % Store event in the appropriate structure
                target_struct.(soc_level).(profile_name).(evtName).channel = channelName;
                target_struct.(soc_level).(profile_name).(evtName).soc_level = soc_level;
                target_struct.(soc_level).(profile_name).(evtName).profile_name = profile_name;
                target_struct.(soc_level).(profile_name).(evtName).stepIndex = stepIndex;
                target_struct.(soc_level).(profile_name).(evtName).event_number = str2double(evtName(6:end));
                target_struct.(soc_level).(profile_name).(evtName).transition_idx = idx1;
                target_struct.(soc_level).(profile_name).(evtName).driving_duration = driving_time;
                target_struct.(soc_level).(profile_name).(evtName).t = t_seg;
                target_struct.(soc_level).(profile_name).(evtName).I = I_seg;
                target_struct.(soc_level).(profile_name).(evtName).V = V_seg;
                target_struct.(soc_level).(profile_name).(evtName).P = P_seg;
                target_struct.(soc_level).(profile_name).(evtName).I_std = current_std;
                target_struct.(soc_level).(profile_name).(evtName).P_std = power_std;
                for dt_idx = 1:length(dt_list)
                    dt_sec = dt_list(dt_idx);
                    if length(t_seg) > 1
                        dt = t_seg(2) - t_seg(1); % Sample interval
                        idx_dt = round(dt_sec / dt + 1); % Integer index
                    else
                        idx_dt = 1;
                    end
                    if idx_dt <= length(I_seg)
                        V1 = V_seg(1);
                        V2 = V_seg(idx_dt);
                        I1 = I_seg(1);
                        I2 = I_seg(idx_dt);
                        dV = V2 - V1;
                        dI = I2 - I1;
                        if dI > 0 && dV > 0
                            dcir_val = (dV / dI) * 1000;
                        elseif dI < 0 && dV < 0
                            dcir_val = (dV / dI) * 1000;
                        else
                            dcir_val = NaN;
                        end
                    else
                        V1 = NaN; V2 = NaN; I1 = NaN; I2 = NaN;
                        dV = NaN; dI = NaN;
                        dcir_val = NaN;
                    end
                    fieldName = sprintf('DCIR_%ds', dt_sec);
                    target_struct.(soc_level).(profile_name).(evtName).(fieldName).val = dcir_val;
                    target_struct.(soc_level).(profile_name).(evtName).(fieldName).V1 = V1;
                    target_struct.(soc_level).(profile_name).(evtName).(fieldName).V2 = V2;
                    target_struct.(soc_level).(profile_name).(evtName).(fieldName).I1 = I1;
                    target_struct.(soc_level).(profile_name).(evtName).(fieldName).I2 = I2;
                    target_struct.(soc_level).(profile_name).(evtName).(fieldName).dV = dV;
                    target_struct.(soc_level).(profile_name).(evtName).(fieldName).dI = dI;
                end
                % Save back to main structure
                if mean(I_seg) > 0
                    Lab_DC_DCIR_0cyc.(chg_struct_name) = target_struct;
                else
                    Lab_DC_DCIR_0cyc.(dchg_struct_name) = target_struct;
                end
            end
            
            fprintf('      → Total events: %d charging, %d discharging\n', chg_event_count, dchg_event_count);
        end
    end
end

%% Save results
save('Lab_DC_DCIR_0cyc_Events.mat', 'Lab_DC_DCIR_0cyc');

%% Summary
fprintf('\n=== Analysis Summary ===\n');
total_events = 0;
for ch_idx = 1:length(channels)
    channelName = channels{ch_idx};
    chg_struct_name = [channelName '_ChgEvent'];
    dchg_struct_name = [channelName '_DchEvent'];
    
    chg_events = 0;
    dchg_events = 0;
    
    if isfield(Lab_DC_DCIR_0cyc, chg_struct_name)
        chg_fields = fieldnames(Lab_DC_DCIR_0cyc.(chg_struct_name));
        for i = 1:length(chg_fields)
            if isfield(Lab_DC_DCIR_0cyc.(chg_struct_name).(chg_fields{i}), 'SOC90')
                events = fieldnames(Lab_DC_DCIR_0cyc.(chg_struct_name).(chg_fields{i}).SOC90);
                chg_events = chg_events + length(events);
            end
        end
    end
    
    if isfield(Lab_DC_DCIR_0cyc, dchg_struct_name)
        dchg_fields = fieldnames(Lab_DC_DCIR_0cyc.(dchg_struct_name));
        for i = 1:length(dchg_fields)
            if isfield(Lab_DC_DCIR_0cyc.(dchg_struct_name).(dchg_fields{i}), 'SOC90')
                events = fieldnames(Lab_DC_DCIR_0cyc.(dchg_struct_name).(dchg_fields{i}).SOC90);
                dchg_events = dchg_events + length(events);
            end
        end
    end
    
    fprintf('%s: %d charging, %d discharging events\n', channelName, chg_events, dchg_events);
    total_events = total_events + chg_events + dchg_events;
end

fprintf('Total events detected: %d\n', total_events);

%% Start Visualization
fprintf('Starting visualization script...\n');
DriveCycle_DCIR_Visualization; 