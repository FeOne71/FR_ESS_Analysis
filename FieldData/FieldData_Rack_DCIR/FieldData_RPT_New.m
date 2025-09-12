%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldData_RPT_New.m
% RPT-style pulse resistance extraction for New data (single day)
% - Input: D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New\2023\Raw_20231016.mat
% - Variables used: SOC_BMS, CVavg, DCPower, DCCurrent
% - Convert to CELL units: I_cell = I/2; P_cell = P*1e3/(17*14*2)
% - Detect pulses in time window 11:27:00–11:31:30 (idle -> negative -> idle)
% - Resistance computed at 10 s from idle-last index; clamp to pulse end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%% Parameters
Cnom = 128;                    % Rack nominal Capacity (Ah)
thr_A = Cnom*0.02;             % Idle threshold (A) for |I| < thr
R_after_sec = 10;              % 10-second evaluation

% winStartTOD = duration(10,46,0); % 2023
% winEndTOD   = duration(10,52,0); % 2023

winStartTOD = duration(10,46,0); % 2025
winEndTOD   = duration(10,52,0); % 2025

% Topology
Ns = 17*14;    % 238s
Np = 2;        % 2p

%% Paths
% dataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New\2023\202310\Raw_20231016.mat';
dataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New\2025\202507\Raw_20250711.mat';
saveDir  = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\FieldData_RPT_results';
if ~exist(saveDir,'dir'); mkdir(saveDir); end

%% Load
if ~exist(dataFile,'file')
    % Try to find a file matching the date in the same folder
    [dpath,~,~] = fileparts(dataFile);
    cand = dir(fullfile(dpath, 'Raw_20250711*.mat'));
    if isempty(cand)
        error('Data file not found: %s', dataFile);
    else
        dataFile = fullfile(dpath, cand(1).name);
    end
end

S = load(dataFile);            % loads variables

% Determine the data struct variable
if isfield(S,'Raw_Rack')
    Raw_Rack = S.Raw_Rack;
else
    % Find first struct containing expected fields
    vars = fieldnames(S);
    Raw_Rack = [];
    for vi = 1:numel(vars)
        val = S.(vars{vi});
        if isstruct(val)
            f = fieldnames(val);
            if any(strcmpi(f,'CVavg')) && any(strcmpi(f,'DCCurrent'))
                Raw_Rack = val; break;
            end
        end
    end
    if isempty(Raw_Rack)
        error('No suitable struct (with CVavg/DCCurrent) found in %s', dataFile);
    end
end

%% Racks
rackNames = {'Rack01'}; % single rack in New data

%% Result
Result = struct();

for r = 1:numel(rackNames)
    rack = rackNames{r};
    % New data may be stored directly as a struct (no Rack01 field)
    if isfield(Raw_Rack, rack)
        D = Raw_Rack.(rack);
    else
        D = Raw_Rack; % direct struct (single rack)
    end

    % Time
    % Date_Time may be duration or string; anchor to known date
    if isfield(D,'Date_Time')
        if isduration(D.Date_Time)
            % t = datetime(2023,10,16) + D.Date_Time;
            t = datetime(2025,07,11) + D.Date_Time;
        else
            t = datetime(D.Date_Time);
        end
    elseif isfield(D,'Time')
        % fallback
        t = datetime(2025,07,11) + duration(string(D.Time),'InputFormat','hh:mm:ss');
    else
        error('No Date_Time/Time found.');
    end
    t0 = t(1);
    tsec = seconds(t - t0);

    % Signals
    I_rack = D.DCCurrent(:);
    Vcell_avg = D.CVavg(:);
    P_rack_kW = D.DCPower(:);   % assumed in kW
    % Optional: soc = D.SOC_BMS(:);

    % To cell units
    I_cell = I_rack / Np;
    P_cell_W = (P_rack_kW*1e3) / (Ns*Np);

    % Window
    tod = timeofday(t);
    mask = (tod >= winStartTOD) & (tod <= winEndTOD);
    idx_map = find(mask);
    if isempty(idx_map)
        Result.(rack) = struct();
        Result.(rack).R_10s_mOhm = [];
        continue;
    end

    % Arrays in window
    t_det = t(mask); tsec_det = tsec(mask);
    I_det = I_rack(mask); I_cell_det = I_cell(mask);
    Vcell_det = Vcell_avg(mask); P_cell_W_det = P_cell_W(mask);

    % Pulse detection: idle -> negative -> idle
    isIdle = abs(I_det) < thr_A;   isNeg = I_det < -thr_A;
    pulses = struct('start_idle_idx',{}, 'neg_start_idx',{}, 'end_idx',{}, 'next_idle_idx',{});
    N = numel(I_det); i = 1;
    while i < N
        if ~isIdle(i), i = i + 1; continue; end
        idle_end = i; while idle_end < N && isIdle(idle_end+1), idle_end = idle_end + 1; end
        j = idle_end + 1; if j > N, break; end
        if ~isNeg(j), i = j; continue; end
        k = j; while k < N && ~isIdle(k+1), k = k + 1; end
        next_idle = k + 1;
        pulses(end+1) = struct('start_idle_idx',idle_end,'neg_start_idx',j,'end_idx',k,'next_idle_idx',next_idle); %#ok<AGROW>
        if numel(pulses) >= 10, break; end
        i = next_idle;
    end

    % Evaluate resistance at 1s, 5s, 10s after start for first up to 10 pulses
    numPulses = min(10, numel(pulses));
    R_1s_mOhm = nan(1, numPulses);
    R_5s_mOhm = nan(1, numPulses);
    R_10s_mOhm = nan(1, numPulses);
    dV_1s_list = nan(1, numPulses);
    dV_5s_list = nan(1, numPulses);
    dV_10s_list = nan(1, numPulses);
    I_1s_list = nan(1, numPulses);
    I_avg_5s_list = nan(1, numPulses);
    I_avg_10s_list = nan(1, numPulses);

    rackStruct = struct();

    % Precompute display names for user readability
    displayNames = {'1stPulse','2ndPulse','3rdPulse','4thPulse','5thPulse', ...
                    '6thPulse','7thPulse','8thPulse','9thPulse','10thPulse'};

    for p = 1:numPulses
        % Use LAST idle sample as the pulse reference (start index)
        idx0 = pulses(p).start_idle_idx;
        t0s  = tsec_det(idx0);
        t1s_1 = t0s + 1;   % 1 second
        t1s_5 = t0s + 5;   % 5 seconds
        t1s_10 = t0s + R_after_sec;  % 10 seconds
        
        % Find indices at different time points
        idx1_1s = find(tsec_det >= t1s_1, 1, 'first');
        idx1_5s = find(tsec_det >= t1s_5, 1, 'first');
        idx1 = find(tsec_det >= t1s_10, 1, 'first');
        
        if isempty(idx1_1s), idx1_1s = numel(tsec_det); end
        if isempty(idx1_5s), idx1_5s = numel(tsec_det); end
        if isempty(idx1), idx1 = numel(tsec_det); end
        
        % ensure evaluation stays within the pulse segment
        pulse_end_idx = pulses(p).end_idx;
        if idx1_1s > pulse_end_idx, idx1_1s = pulse_end_idx; end
        if idx1_5s > pulse_end_idx, idx1_5s = pulse_end_idx; end
        if idx1 > pulse_end_idx, idx1 = pulse_end_idx; end

        % Compute dV for different time points
        dV_1s = Vcell_det(idx1_1s) - Vcell_det(idx0);      % V change at 1s
        dV_5s = Vcell_det(idx1_5s) - Vcell_det(idx0);      % V change at 5s  
        dV_10s = Vcell_det(idx1) - Vcell_det(idx0);        % V change at 10s
        
        % Get current values for different resistance calculations
        pulse_start = pulses(p).neg_start_idx;
        pulse_end = pulses(p).end_idx;
        
        % R1s: use current at 1s point
        I_1s = I_cell_det(idx1_1s);
        
        % R5s: use average current from I(1) to I(5)
        if pulse_end - pulse_start >= 4  % ensure we have at least 5 samples
            I_avg_5s = mean(I_cell_det(pulse_start:pulse_start+4));  % I(1) to I(5)
        else
            I_avg_5s = mean(I_cell_det(pulse_start:pulse_end));  % fallback
        end
        
        % R10s: use average current from I(3) to I(10)
        if pulse_end - pulse_start >= 7  % ensure we have at least 8 samples (3 to 10)
            I_avg_10s = mean(I_cell_det(pulse_start+2:pulse_start+9));  % I(3) to I(10)
        else
            I_avg_10s = mean(I_cell_det(pulse_start:pulse_end));  % fallback
        end
        
        % Calculate resistances
        if I_1s == 0
            R_1s = NaN;
        else
            R_1s = (dV_1s / I_1s) * 1000;
        end
        
        if I_avg_5s == 0
            R_5s = NaN;
        else
            R_5s = (dV_5s / I_avg_5s) * 1000;
        end
        
        if I_avg_10s == 0
            R_10s = NaN;
        else
            R_10s = (dV_10s / I_avg_10s) * 1000;
        end

        R_1s_mOhm(p) = R_1s;
        R_5s_mOhm(p) = R_5s;
        R_10s_mOhm(p) = R_10s;
        
        dV_1s_list(p) = dV_1s;
        dV_5s_list(p) = dV_5s;
        dV_10s_list(p) = dV_10s;
        I_1s_list(p) = I_1s;
        I_avg_5s_list(p) = I_avg_5s;
        I_avg_10s_list(p) = I_avg_10s;

        % Store sequences around the pulse start up to t1s
        % store the entire pulse segment from last idle to just before next idle
        seg_start = pulses(p).start_idle_idx;
        seg_end   = pulses(p).end_idx;
        seq_idx_det = seg_start:seg_end;
        seq_idx = idx_map(seq_idx_det);

        pulseStruct = struct();
        pulseStruct.t_datetime = t(seq_idx);
        pulseStruct.t_sec = tsec(seq_idx);
        pulseStruct.Vcell_avg = Vcell_avg(seq_idx);
        pulseStruct.I_cell = I_cell(seq_idx);
        pulseStruct.P_cell_W = P_cell_W(seq_idx);
        pulseStruct.dV_1s = dV_1s; pulseStruct.I_1s = I_1s; pulseStruct.R_1s_mOhm = R_1s;
        pulseStruct.dV_5s = dV_5s; pulseStruct.I_avg_5s = I_avg_5s; pulseStruct.R_5s_mOhm = R_5s;
        pulseStruct.dV_10s = dV_10s; pulseStruct.I_avg_10s = I_avg_10s; pulseStruct.R_10s_mOhm = R_10s;

        % Valid MATLAB field name (cannot start with digit)
        pulseName = sprintf('Pulse%02d', p);
        % pulseStruct.displayName = displayNames{p};
        rackStruct.(pulseName) = pulseStruct;
    end

    Result.(rack) = rackStruct;

    % Keep for plotting and table
    Result.(rack).R_1s_mOhm = R_1s_mOhm;
    Result.(rack).R_5s_mOhm = R_5s_mOhm;
    Result.(rack).R_10s_mOhm = R_10s_mOhm;
end

%% Create resistance table
fprintf('\n=== Resistance Table (Cell Units) ===\n');
fprintf('%-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'Rack', 'Pulse', 'R1s[mΩ]', 'R5s[mΩ]', 'R10s[mΩ]', 'dV1s[V]', 'dV5s[V]', 'dV10s[V]', 'I_avg[A]');
fprintf('%s\n', repmat('-', 1, 80));

for r = 1:numel(rackNames)
    rack = rackNames{r};
    if ~isfield(Result, rack) || ~isfield(Result.(rack),'R_1s_mOhm')
        continue;
    end
    
    R_1s = Result.(rack).R_1s_mOhm;
    R_5s = Result.(rack).R_5s_mOhm;
    R_10s = Result.(rack).R_10s_mOhm;
    
    for p = 1:min(10, numel(R_1s))
        if ~isnan(R_1s(p))
            % Get dV and I values from first pulse for display
            pulseName = sprintf('Pulse%02d', p);
            if isfield(Result.(rack), pulseName)
                pulse = Result.(rack).(pulseName);
                dV1s = pulse.dV_1s;
                dV5s = pulse.dV_5s;
                dV10s = pulse.dV_10s;
                I_avg = pulse.I_avg_10s;  % use 10s average for display
            else
                dV1s = NaN; dV5s = NaN; dV10s = NaN; I_avg = NaN;
            end
            
            fprintf('%-8s %-8d %-8.2f %-8.2f %-8.2f %-8.3f %-8.3f %-8.3f %-8.2f\n', ...
                rack, p, R_1s(p), R_5s(p), R_10s(p), dV1s, dV5s, dV10s, I_avg);
        end
    end
end

%% Plot
figure('Name','RPT 10s Resistance per Pulse (New, Discharge)','NumberTitle','off'); hold on; grid on;
rackColors = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840; 0 0 0];
rackMarkers = {'o','s','^','v','d','>','<','p'};
for r = 1:numel(rackNames)
    rack = rackNames{r};
    if ~isfield(Result,rack) || ~isfield(Result.(rack),'R_10s_mOhm'), continue; end
    Rv = Result.(rack).R_10s_mOhm; col = rackColors(min(r, size(rackColors,1)), :); mk = rackMarkers{mod(r-1,numel(rackMarkers))+1};
    plot(1:numel(Rv), Rv, '-', 'Color', col, 'LineWidth', 1.8, 'Marker', mk, 'MarkerSize', 6, 'DisplayName', rack);
end
xlabel('Pulse index'); ylabel('R_{10s} [m\Omega]'); 
title('Per-rack 10s Resistance of first 10 negative-current pulses');
xticks(1:10); legend('Location','best');
set(findall(gcf,'-property','FontSize'),'FontSize',14);

%% Save
% save(fullfile(saveDir,'FieldData_RPT_New_20231016.mat'), 'Result', '-v7.3');
save(fullfile(saveDir,'FieldData_RPT_New_20250711.mat'), 'Result', '-v7.3');
fprintf('Completed RPT (New) for %s\n', dataFile);


