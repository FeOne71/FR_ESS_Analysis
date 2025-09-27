%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FielddQdV_Simple.m
% ESS Rack01 - Simple Grouping and Visualization
% - Charge event detection (>=180s)
% - 2 power groups + voltage start/end groups
% - Visualization: x-axis=time, y-axis=voltage, yearly curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_dQdV','dQdV_Simple_Results');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

%% Parameters
Cnom_cell = 64;
idle_thr = Cnom_cell*0.05;
charge_thr = idle_thr;
min_avg_current = Cnom_cell*0.15;
min_charge_duration = 180;
moving_avg_window = 30;  % Moving average window for voltage smoothing

%% Rack Configuration
Np = 2;

%% Targets
yearList = {'2021','2022','2023'};
rackNames_all = {'Rack01'};

%% Grouping parameters
voltage_bins = 3.6:0.01:4.2;
n_voltage_bins = numel(voltage_bins) - 1;

%% Colors for years
year_colors = struct();
year_colors.year_2021 = [0 0.8 0];  % Green
year_colors.year_2022 = [0 0 1];    % Blue  
year_colors.year_2023 = [1 0 0];    % Red

%% Step 1: Detect charge events
fprintf('=== Step 1: Charge event detection (>=%ds, avg_current>=%.1fA) ===\n', min_charge_duration, min_avg_current);
all_charge_events = struct();

for yi = 1:numel(yearList)
    year = yearList{yi};
    fprintf('Processing year: %s\n', year);
    
    yearPath = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));
    
    year_events = struct();
    ev_count = 0;
    
    for m = 1:numel(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, '*.mat'));
        [~, idx] = sort({matFiles.name}); matFiles = matFiles(idx);
        
        for f = 1:numel(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            S = load(matFilePath);
            if ~isfield(S, 'Raw'), continue; end
            
            for ri = 1:numel(rackNames_all)
                rackName = rackNames_all{ri};
                if ~isfield(S.Raw, rackName), continue; end
                
                R = S.Raw.(rackName);
                t = R.Time;
                I = R.DCCurrent_A / Np;
                V = R.AverageCV_V;
                P = R.DCPower_kW;
                
                t = seconds(datetime(t) - datetime(t(1)));
                
                % Find charge events
                ch_mask = (I > charge_thr) & (I > 0);
                charge_starts = find(diff([false; ch_mask]) == 1);
                charge_ends = find(diff([ch_mask; false]) == -1);
                
                if numel(charge_starts) > numel(charge_ends)
                    charge_starts = charge_starts(1:numel(charge_ends));
                elseif numel(charge_ends) > numel(charge_starts)
                    charge_ends = charge_ends(1:numel(charge_starts));
                end
                
                for si = 1:numel(charge_starts)
                    sIdx = charge_starts(si);
                    eIdx = charge_ends(si);
                    dur = t(eIdx) - t(sIdx);
                    
                    if dur >= min_charge_duration
                        avg_current = mean(I(sIdx:eIdx), 'omitnan');
                        if avg_current >= min_avg_current
                            ev_count = ev_count + 1;
                            key = sprintf('Event_%d', ev_count);
                            year_events.(key).time = t(sIdx:eIdx) - t(sIdx);
                            year_events.(key).voltage = V(sIdx:eIdx);
                            year_events.(key).power = mean(P(sIdx:eIdx), 'omitnan');
                            year_events.(key).start_voltage = V(sIdx);
                            year_events.(key).end_voltage = V(eIdx);
                        end
                    end
                end
            end
        end
    end
    
    all_charge_events.(['year_' year]) = year_events;
    fprintf('Year %s: %d charge events detected\n', year, ev_count);
end

%% Step 2: Determine power groups
fprintf('\n=== Step 2: Determine power groups ===\n');

all_powers = [];
for yi = 1:numel(yearList)
    year = yearList{yi};
    E = all_charge_events.(['year_' year]);
    names = fieldnames(E);
    for e = 1:numel(names)
        all_powers(end+1) = E.(names{e}).power;
    end
end

[~, power_centers] = kmeans(all_powers(:), 2, 'Replicates', 10);
power_centers = sort(power_centers);
power_ranges = [power_centers(1)*0.9, power_centers(1)*1.1; power_centers(2)*0.9, power_centers(2)*1.1];

fprintf('Power groups: %.1f-%.1f kW, %.1f-%.1f kW\n', power_ranges(1,1), power_ranges(1,2), power_ranges(2,1), power_ranges(2,2));

%% Step 3: Group and visualize
fprintf('\n=== Step 3: Group and visualize ===\n');

valid_groups = 0;
for pi = 1:2  % 2 power groups
    for si = 1:n_voltage_bins  % start voltage bins
        for ei = 1:n_voltage_bins  % end voltage bins
            
            % Collect events for this group
            group_events = struct();
            all_durations = [];
            
            for yi = 1:numel(yearList)
                year = yearList{yi};
                E = all_charge_events.(['year_' year]);
                names = fieldnames(E);
                
                year_events = [];
                for e = 1:numel(names)
                    ev = E.(names{e});
                    
                    % Check if event matches this group
                    if ev.power >= power_ranges(pi,1) && ev.power < power_ranges(pi,2)
                        start_bin = find(ev.start_voltage >= voltage_bins(1:end-1) & ev.start_voltage < voltage_bins(2:end), 1);
                        end_bin = find(ev.end_voltage >= voltage_bins(1:end-1) & ev.end_voltage < voltage_bins(2:end), 1);
                        
                        if ~isempty(start_bin) && ~isempty(end_bin) && start_bin == si && end_bin == ei
                            year_events = [year_events; e];
                            all_durations = [all_durations; numel(ev.time)];
                        end
                    end
                end
                
                group_events.(['year_' year]) = year_events;
            end
            
            % Check if all years have at least 5 events
            all_years_have_enough = true;
            for yi = 1:numel(yearList)
                year = yearList{yi};
                if numel(group_events.(['year_' year])) < 5
                    all_years_have_enough = false;
                    break;
                end
            end
            
            if all_years_have_enough
                valid_groups = valid_groups + 1;
                common_time = min(all_durations);
                
                fprintf('Valid group %d: P%d_S%d_E%d (total events: %d, common time: %d)\n', ...
                    valid_groups, pi, si, ei, sum(cellfun(@numel, struct2cell(group_events))), common_time);
                
                % Create individual events figure
                figure('Position',[100,100,1200,800]);
                hold on; grid on;
                
                for yi = 1:numel(yearList)
                    year = yearList{yi};
                    year_events = group_events.(['year_' year]);
                    E = all_charge_events.(['year_' year]);
                    names = fieldnames(E);
                    
                    for e = 1:numel(year_events)
                        ev = E.(names{year_events(e)});
                        t = ev.time(1:min(common_time, numel(ev.time)));
                        V = ev.voltage(1:min(common_time, numel(ev.voltage)));
                        
                        % Apply moving average to voltage for smoothing
                        if numel(V) >= moving_avg_window
                            V_smooth = movmean(V, moving_avg_window);
                        else
                            V_smooth = V;  % Use original if not enough data points
                        end
                        
                        color_key = ['year_' year];
                        plot(t, V_smooth, '-', 'Color', year_colors.(color_key), 'LineWidth', 0.5);
                    end
                end
                
                xlabel('Time (s)');
                ylabel('Voltage (V)');
                title(sprintf('Group P%d_S%d_E%d: %.1f-%.1fkW, %.2f-%.2fV to %.2f-%.2fV', ...
                    pi, si, ei, power_ranges(pi,1), power_ranges(pi,2), ...
                    voltage_bins(si), voltage_bins(si+1), voltage_bins(ei), voltage_bins(ei+1)));
                
                % Create legend
                legend_handles = [];
                legend_labels = {};
                for yi = 1:numel(yearList)
                    year = yearList{yi};
                    if numel(group_events.(['year_' year])) > 0
                        color_key = ['year_' year];
                        h = plot(NaN, NaN, '-', 'Color', year_colors.(color_key), 'LineWidth', 2);
                        legend_handles = [legend_handles, h];
                        legend_labels{end+1} = sprintf('%s (%d events)', year, numel(group_events.(['year_' year])));
                    end
                end
                legend(legend_handles, legend_labels, 'Location', 'best');
                
                % Save figure
                outFig = fullfile(saveDir, sprintf('Group_P%d_S%d_E%d.fig', pi, si, ei));
                saveas(gcf, outFig);
                
                % Create average curves figure
                figure('Position',[100,100,1200,800]);
                hold on; grid on;
                
                avg_legend_handles = [];
                avg_legend_labels = {};
                
                for yi = 1:numel(yearList)
                    year = yearList{yi};
                    year_events = group_events.(['year_' year]);
                    E = all_charge_events.(['year_' year]);
                    names = fieldnames(E);
                    
                    if numel(year_events) > 0
                        % Collect all voltage data
                        all_times = [];
                        all_voltages = [];
                        
                        for e = 1:numel(year_events)
                            ev = E.(names{year_events(e)});
                            t = ev.time(1:min(common_time, numel(ev.time)));
                            V = ev.voltage(1:min(common_time, numel(ev.voltage)));
                            all_times = [all_times; t];
                            all_voltages = [all_voltages; V];
                        end
                        
                        % Calculate average curve
                        time_bins = 0:1:common_time;
                        avg_voltages = zeros(size(time_bins));
                        
                        for ti = 1:numel(time_bins)
                            bin_mask = abs(all_times - time_bins(ti)) <= 0.5;
                            if sum(bin_mask) > 0
                                avg_voltages(ti) = mean(all_voltages(bin_mask), 'omitnan');
                            end
                        end
                        
                        % Plot average curve
                        valid_mask = ~isnan(avg_voltages) & avg_voltages > 0;
                        if sum(valid_mask) > 1
                            color_key = ['year_' year];
                            h = plot(time_bins(valid_mask), avg_voltages(valid_mask), '-', 'Color', year_colors.(color_key), 'LineWidth', 3);
                            avg_legend_handles = [avg_legend_handles, h];
                            avg_legend_labels{end+1} = sprintf('%s (avg, %d events)', year, numel(year_events));
                        end
                    end
                end
                
                xlabel('Time (s)');
                ylabel('Voltage (V)');
                title(sprintf('Group P%d_S%d_E%d - Average Curves: %.1f-%.1fkW, %.2f-%.2fV to %.2f-%.2fV', ...
                    pi, si, ei, power_ranges(pi,1), power_ranges(pi,2), ...
                    voltage_bins(si), voltage_bins(si+1), voltage_bins(ei), voltage_bins(ei+1)));
                legend(avg_legend_handles, avg_legend_labels, 'Location', 'best');
                
                % Save average figure
                outFig_avg = fullfile(saveDir, sprintf('Group_P%d_S%d_E%d_Average.fig', pi, si, ei));
                saveas(gcf, outFig_avg);
                
                % Create dQ/dV and dV/dQ figure
                figure('Position',[100,100,1200,800]);
                subplot(2,1,1); % dQ/dV plot
                hold on; grid on;
                
                dQdV_legend_handles = [];
                dQdV_legend_labels = {};
                
                for yi = 1:numel(yearList)
                    year = yearList{yi};
                    year_events = group_events.(['year_' year]);
                    E = all_charge_events.(['year_' year]);
                    names = fieldnames(E);
                    
                    if numel(year_events) > 0
                        % Collect all voltage and current data for dQ/dV calculation
                        all_voltages_dQdV = [];
                        all_dQdV = [];
                        
                        for e = 1:numel(year_events)
                            ev = E.(names{year_events(e)});
                            t = ev.time(1:min(common_time, numel(ev.time)));
                            V = ev.voltage(1:min(common_time, numel(ev.voltage)));
                            
                            % Calculate dQ/dV (current / dV/dt)
                            if numel(V) > 1
                                dV = diff(V);
                                dt = diff(t);
                                dV_dt = dV ./ dt;
                                
                                % Avoid division by zero
                                valid_idx = abs(dV_dt) > 1e-6;
                                if sum(valid_idx) > 0
                                    % For dQ/dV, we need current data - approximate from power/voltage
                                    % Since we don't have current stored, we'll use a simplified approach
                                    % dQ/dV ≈ (Power/Voltage) / (dV/dt)
                                    P_avg = ev.power * 1000; % Convert to W
                                    V_mid = (V(1:end-1) + V(2:end)) / 2;
                                    I_approx = P_avg ./ V_mid; % Approximate current
                                    
                                    dQ_dV = I_approx(valid_idx) ./ dV_dt(valid_idx);
                                    V_plot = V_mid(valid_idx);
                                    
                                    all_voltages_dQdV = [all_voltages_dQdV; V_plot];
                                    all_dQdV = [all_dQdV; dQ_dV];
                                end
                            end
                        end
                        
                        % Calculate average dQ/dV curve
                        if ~isempty(all_voltages_dQdV)
                            voltage_bins_dQdV = voltage_bins;
                            avg_dQdV = zeros(size(voltage_bins_dQdV));
                            
                            for vi = 1:numel(voltage_bins_dQdV)
                                bin_mask = abs(all_voltages_dQdV - voltage_bins_dQdV(vi)) <= 0.01;
                                if sum(bin_mask) > 0
                                    avg_dQdV(vi) = mean(all_dQdV(bin_mask), 'omitnan');
                                end
                            end
                            
                            % Plot average dQ/dV curve
                            valid_mask = ~isnan(avg_dQdV) & avg_dQdV ~= 0;
                            if sum(valid_mask) > 1
                                color_key = ['year_' year];
                                h = plot(voltage_bins_dQdV(valid_mask), avg_dQdV(valid_mask), '-', 'Color', year_colors.(color_key), 'LineWidth', 2);
                                dQdV_legend_handles = [dQdV_legend_handles, h];
                                dQdV_legend_labels{end+1} = sprintf('%s (avg, %d events)', year, numel(year_events));
                            end
                        end
                    end
                end
                
                xlabel('Voltage (V)');
                ylabel('dQ/dV (Ah/V)');
                title(sprintf('Group P%d_S%d_E%d - dQ/dV Curves: %.1f-%.1fkW, %.2f-%.2fV to %.2f-%.2fV', ...
                    pi, si, ei, power_ranges(pi,1), power_ranges(pi,2), ...
                    voltage_bins(si), voltage_bins(si+1), voltage_bins(ei), voltage_bins(ei+1)));
                legend(dQdV_legend_handles, dQdV_legend_labels, 'Location', 'best');
                
                subplot(2,1,2); % dV/dQ plot
                hold on; grid on;
                
                dVdQ_legend_handles = [];
                dVdQ_legend_labels = {};
                
                for yi = 1:numel(yearList)
                    year = yearList{yi};
                    year_events = group_events.(['year_' year]);
                    E = all_charge_events.(['year_' year]);
                    names = fieldnames(E);
                    
                    if numel(year_events) > 0
                        % Collect all voltage and current data for dV/dQ calculation
                        all_voltages_dVdQ = [];
                        all_dVdQ = [];
                        
                        for e = 1:numel(year_events)
                            ev = E.(names{year_events(e)});
                            t = ev.time(1:min(common_time, numel(ev.time)));
                            V = ev.voltage(1:min(common_time, numel(ev.voltage)));
                            
                            % Calculate dV/dQ (dV/dt / current)
                            if numel(V) > 1
                                dV = diff(V);
                                dt = diff(t);
                                dV_dt = dV ./ dt;
                                
                                % For dV/dQ, we need current data - approximate from power/voltage
                                P_avg = ev.power * 1000; % Convert to W
                                V_mid = (V(1:end-1) + V(2:end)) / 2;
                                I_approx = P_avg ./ V_mid; % Approximate current
                                
                                % Avoid division by zero
                                valid_idx = abs(I_approx) > 1e-6;
                                if sum(valid_idx) > 0
                                    dV_dQ = dV_dt(valid_idx) ./ I_approx(valid_idx);
                                    V_plot = V_mid(valid_idx);
                                    
                                    all_voltages_dVdQ = [all_voltages_dVdQ; V_plot];
                                    all_dVdQ = [all_dVdQ; dV_dQ];
                                end
                            end
                        end
                        
                        % Calculate average dV/dQ curve
                        if ~isempty(all_voltages_dVdQ)
                            voltage_bins_dVdQ = voltage_bins;
                            avg_dVdQ = zeros(size(voltage_bins_dVdQ));
                            
                            for vi = 1:numel(voltage_bins_dVdQ)
                                bin_mask = abs(all_voltages_dVdQ - voltage_bins_dVdQ(vi)) <= 0.01;
                                if sum(bin_mask) > 0
                                    avg_dVdQ(vi) = mean(all_dVdQ(bin_mask), 'omitnan');
                                end
                            end
                            
                            % Plot average dV/dQ curve
                            valid_mask = ~isnan(avg_dVdQ) & avg_dVdQ ~= 0;
                            if sum(valid_mask) > 1
                                color_key = ['year_' year];
                                h = plot(voltage_bins_dVdQ(valid_mask), avg_dVdQ(valid_mask), '-', 'Color', year_colors.(color_key), 'LineWidth', 2);
                                dVdQ_legend_handles = [dVdQ_legend_handles, h];
                                dVdQ_legend_labels{end+1} = sprintf('%s (avg, %d events)', year, numel(year_events));
                            end
                        end
                    end
                end
                
                xlabel('Voltage (V)');
                ylabel('dV/dQ (V/Ah)');
                title(sprintf('Group P%d_S%d_E%d - dV/dQ Curves: %.1f-%.1fkW, %.2f-%.2fV to %.2f-%.2fV', ...
                    pi, si, ei, power_ranges(pi,1), power_ranges(pi,2), ...
                    voltage_bins(si), voltage_bins(si+1), voltage_bins(ei), voltage_bins(ei+1)));
                legend(dVdQ_legend_handles, dVdQ_legend_labels, 'Location', 'best');
                
                % Save dQ/dV and dV/dQ figure
                outFig_dQdV = fullfile(saveDir, sprintf('Group_P%d_S%d_E%d_dQdV_dVdQ.fig', pi, si, ei));
                saveas(gcf, outFig_dQdV);
            end
        end
    end
end

fprintf('\nTotal valid groups: %d\n', valid_groups);
fprintf('All done. Results saved to: %s\n', saveDir);
