%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldData_RPT_Old.m
% RPT-style pulse resistance extraction for Old data (single day)
% Requirements from user:
%  - Load: Rack_raw2mat\2021\202106\Raw_20210607.mat
%  - Data has 8 racks (Rack01..Rack08) in struct Raw_Rack
%  - Convert rack signals to cell units:
%       I_cell = DCCurrent_A / 2
%       P_cell = DCPower_kW * 1e3 / (17*14*2)
%  - Idle if |I| < threshold (use OldData_v1 parameters)
%  - Extract pulses that leave idle into negative current (discharge)
%  - For each of the first 10 pulses, compute resistance 10 s after pulse start
%  - Plot per-rack pulse resistances (1..10) in one figure
%  - Save structure: RackXX with 1stPulse..10thPulse, store time series and dV,dI and R metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% Parameters (from OldData_v1.m)
Cnom = 128;                    % Rack nominal Capacity (Ah)
thr_A = Cnom*0.02;              % Idle threshold (A) for |I| < thr
R_after_sec = 10;              % Resistance evaluation time after pulse start (s)
useWindow = true;              % detect only inside a time-of-day window
winStartTOD = duration(14,45,0); % inclusive start (time-of-day)
winEndTOD   = duration(14,50,0); % inclusive end (time-of-day)

% Topology
Ns = 17*14;                    % series strings per rack (238s)
Np = 2;                        % parallel strings per rack (2p)

%% Paths
dataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\2021\202106\Raw_20210607.mat';
saveDir  = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\FieldData_RPT_results';
if ~exist(saveDir,'dir'); mkdir(saveDir); end

%% Load
S = load(dataFile);            % loads Raw_Rack
if ~isfield(S,'Raw_Rack')
    error('Raw_Rack not found in %s', dataFile);
end
Raw_Rack = S.Raw_Rack;

%% Rack list
rackNames = {'Rack01','Rack02','Rack03','Rack04','Rack05','Rack06','Rack07','Rack08'};

%% Helpers
% (No smoothing per user's request)

%% Result structure
Result = struct();

%% Process each rack
for r = 1:numel(rackNames)
    rack = rackNames{r};
    if ~isfield(Raw_Rack, rack); continue; end
    D = Raw_Rack.(rack);

    % Time handling
    t = datetime(D.Time);          % string -> datetime
    t0 = t(1);
    tsec = seconds(t - t0);        % seconds from start

    % Signals (rack units)
    I_rack = D.DCCurrent_A(:);                     % A
    Vcell_avg = D.AverageCV_V(:);                  % V (cell average)
    P_rack_kW = D.DCPower_kW(:);                   % kW (rack)

    % Convert to cell units
    I_cell = I_rack / Np;                          % A per cell
    P_cell_W = (P_rack_kW*1e3) / (Ns*Np);          % W per cell

    % Use cell-level signals only (per requirement)

    % Restrict to time-of-day window if enabled
    if useWindow
        tod = timeofday(t);
        mask = (tod >= winStartTOD) & (tod <= winEndTOD);
    else
        mask = true(size(t));
    end
    idx_map = find(mask);
    if isempty(idx_map)
        Result.(rack) = struct();
        Result.(rack).R_10s_mOhm = [];
        continue;
    end

    % Detection arrays within window
    t_det = t(mask);
    tsec_det = tsec(mask);
    
    I_det = I_rack(mask);
    I_cell_det = I_cell(mask);    
    
    Vcell_det = Vcell_avg(mask);
    
    P_rack_kW_det = P_rack_kW(mask);
    P_cell_W_det = P_cell_W(mask);

    % Find pulses as segments within window: idle (|I|<thr) -> negative (I<-thr) -> idle
    isIdle = abs(I_det) < thr_A;     % idle mask
    isNeg  = I_det < -thr_A;         % negative mask
    pulses = struct('start_idle_idx',{}, 'neg_start_idx',{}, 'end_idx',{}, 'next_idle_idx',{});

    N = numel(I_det);
    i = 1;
    while i < N
        if ~isIdle(i)
            i = i + 1; continue; % need to start from idle
        end

        % extend through idle region
        idle_end = i;
        while idle_end < N && isIdle(idle_end+1)
            idle_end = idle_end + 1;
        end

        % first sample after idle
        j = idle_end + 1;
        if j > N
            break;
        end

        % require negative region to form a pulse
        if ~isNeg(j)
            i = j; % move forward to seek next idle
            continue;
        end

        % walk through negative region until next idle begins
        k = j;
        while k < N && ~isIdle(k+1)
            k = k + 1;
        end

        % now k is the last non-idle sample within the pulse
        next_idle = k + 1;            % first idle sample after pulse
        pulses(end+1) = struct( ...
            'start_idle_idx', idle_end, ...
            'neg_start_idx',  j, ...
            'end_idx',        k, ...
            'next_idle_idx',  next_idle); %#ok<AGROW>

        % limit to first 10 pulses
        if numel(pulses) >= 10
            break;
        end

        % continue scanning from next idle region
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
        seq_idx = idx_map(seq_idx_det); % map back to original indices
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

%% Visualization: Time series plots for winStartTOD to winEndTOD
figure('Name','Time Series: Current, Voltage, Power (14:45-14:50)','NumberTitle','off');
tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% Get data from first rack for time series plot
rack = rackNames{1};
if isfield(Raw_Rack, rack)
    D = Raw_Rack.(rack);
    t = datetime(D.Time);
    I_rack = D.DCCurrent_A(:);
    Vcell_avg = D.AverageCV_V(:);
    P_rack_kW = D.DCPower_kW(:);
    
    % Convert to cell units
    I_cell = I_rack / Np;                          % A per cell
    P_cell_W = (P_rack_kW*1e3) / (Ns*Np);          % W per cell
    
    % Apply time window
    tod = timeofday(t);
    mask = (tod >= winStartTOD) & (tod <= winEndTOD);
    t_win = t(mask);
    I_win = I_cell(mask);
    V_win = Vcell_avg(mask);
    P_win = P_cell_W(mask);
    
    % (1,1) Current
    ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');
    plot(ax1, t_win, I_win, '-', 'Color', [0.8 0.3 0.1], 'LineWidth', 3);
    title(ax1, 'Cell Current I_{cell} [A]'); xlabel(ax1,'Time'); ylabel(ax1,'I_{cell} [A]');
    xlim(ax1, [min(t_win) max(t_win)]);
    
    % (2,1) Voltage
    ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');
    plot(ax2, t_win, V_win, '-', 'Color', [0.2 0.6 0.8], 'LineWidth', 3);
    title(ax2, 'Average Cell Voltage V_{cell,avg} [V]'); xlabel(ax2,'Time'); ylabel(ax2,'V_{cell,avg} [V]');
    xlim(ax2, [min(t_win) max(t_win)]);
    
    % (3,1) Power
    ax3 = nexttile(tl,3); hold(ax3,'on'); grid(ax3,'on');
    plot(ax3, t_win, P_win, '-', 'Color', [0.5 0.2 0.6], 'LineWidth', 3);
    title(ax3, 'Cell Power P_{cell} [W]'); xlabel(ax3,'Time'); ylabel(ax3,'P_{cell} [W]');
    xlim(ax3, [min(t_win) max(t_win)]);
    
    % Format time axis
    for ax = [ax1, ax2, ax3]
        ax.XTickLabelRotation = 45;
        ax.TickDir = 'out';
    end
end

set(findall(gcf,'-property','FontSize'),'FontSize',12);

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

%% Visualization: one figure with 8 racks, pulse 1..10 resistance
figure('Name','RPT 10s Resistance per Pulse (Discharge)','NumberTitle','off');
hold on; grid on;
rackColors = [
    0.0000 0.4470 0.7410;  % blue
    0.8500 0.3250 0.0980;  % orange
    0.9290 0.6940 0.1250;  % yellow
    0.4940 0.1840 0.5560;  % purple
    0.4660 0.6740 0.1880;  % green
    0.3010 0.7450 0.9330;  % cyan
    0.6350 0.0780 0.1840;  % red-brown
    0.0000 0.0000 0.0000   % black
];
rackMarkers = {'o','s','^','v','d','>','<','p'};
for r = 1:numel(rackNames)
    rack = rackNames{r};
    if ~isfield(Result, rack) || ~isfield(Result.(rack),'R_10s_mOhm'); continue; end
    Rv = Result.(rack).R_10s_mOhm;
    col = rackColors(min(r, size(rackColors,1)), :);
    mark = rackMarkers{mod(r-1, numel(rackMarkers)) + 1};
    plot(1:numel(Rv), Rv, '-', 'Color', col, 'LineWidth', 1.8, ...
         'Marker', mark, 'MarkerSize', 6, 'DisplayName', rack);
end
xlabel('Pulse index'); ylabel('R_{10s} [m\Omega]');
title('Per-rack 10s Resistance of first 10 negative-current pulses');
xticks(1:10);
legend('Location','best');
set(findall(gcf,'-property','FontSize'),'FontSize',14);

%% Save result structure
save(fullfile(saveDir,'FieldData_RPT_Old_20210607.mat'), 'Result', '-v7.3');

fprintf('Completed RPT extraction for %s\n', dataFile);


