%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FielddQdV_OldData.m
% ESS Rack01 - HI(Partial Capacity) by Current Groups and Start Voltage
% - Charge event detection (>=30s, avg_current>=9.6A)
% - 2D grouping: (average current range) × (start voltage range)
% - Voltage bins: 3.4~4.2V, 0.01V intervals (80 bins)
% - HI(V_i) = sum(I*dt) within bin [mAh]
% - Visualization: per group, yearly curves with different colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_dQdV','dQdV_Results_OldData');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

%% Parameters
Cnom_cell = 64;                 % nominal per-cell capacity [Ah]
idle_thr   = Cnom_cell*0.05;    % [A] idle threshold (charge/discharge)
charge_thr = idle_thr;          % [A] charge detection threshold
min_avg_current = Cnom_cell*0.15; % [A] minimum average current for charge events (9.6A)

%% Rack Configuration
Np = 2;       % 2p

%% Visualization parameters
Fontsize  = 12;
LineWidth = 2;

%% Targets
yearList       = {'2021','2022','2023'};
rackNames_all  = {'Rack01'};

% Fixed voltage bins: 3.6~4.0V, 0.01V intervals
edges_V   = 3.6:0.01:4.0;
centers_V = (edges_V(1:end-1)+edges_V(2:end))/2;
nVB       = numel(centers_V);

%% Grouping parameters
% Start voltage range (V) - must be within voltage bin range (3.6-4.0V)
start_voltage_range = [3.6, 4.0];  % [min, max] - single range
% Current groups will be determined from actual data distribution

%% Containers
all_charge_events = struct();

%% Step 1: Detect charge events (>=30s, avg_current>=9.6A)
fprintf('=== Step 1: Charge event detection (>=30s, avg_current>=%.1fA) ===\n', min_avg_current);
for yi = 1:numel(yearList)
    year = yearList{yi};
    fprintf('Processing year: %s\n', year);

    yearPath  = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));

    year_events = struct();
    ev_count = 0;

    for m = 1:numel(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles  = dir(fullfile(monthPath, '*.mat'));
        [~, idx]  = sort({matFiles.name}); matFiles = matFiles(idx);

        for f = 1:numel(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            S = load(matFilePath);
            if ~isfield(S, 'Raw'), continue; end

            for ri = 1:numel(rackNames_all)
                rackName = rackNames_all{ri};
                if ~isfield(S.Raw, rackName), continue; end

                R = S.Raw.(rackName);
                t   = R.Time;
                I   = R.DCCurrent_A;
                V   = R.AverageCV_V;
                soc = R.SOCPct;
                P   = R.DCPower_kW;  % Use existing DC Power variable

                % Per-cell current
                I = I / Np;

                % time in seconds-from-start
                t = seconds(datetime(t) - datetime(t(1)));

                % Detect charge (I > 0 and above threshold)
                ch_mask = (I > charge_thr) & (I > 0);

                % Continuous segments of charge
                segs = find_continuous_segments(ch_mask);

                for si = 1:numel(segs)
                    sIdx = segs(si).start;
                    eIdx = segs(si).end;
                    dur  = t(eIdx) - t(sIdx);
                    if dur >= 30
                        % Calculate average current and power for this segment
                        avg_current = mean(I(sIdx:eIdx), 'omitnan');
                        avg_power = mean(P(sIdx:eIdx), 'omitnan');  % Use existing DC Power
                        
                        % Filter by minimum average current (9.6A)
                        if avg_current >= min_avg_current
                            ev_count = ev_count + 1;
                            key = sprintf('Event_%d', ev_count);
                            year_events.(key).time     = t(sIdx:eIdx);
                            year_events.(key).current  = I(sIdx:eIdx);
                            year_events.(key).voltage  = V(sIdx:eIdx);
                            year_events.(key).power    = P(sIdx:eIdx);
                            year_events.(key).soc      = soc(sIdx:eIdx);
                            year_events.(key).date     = matFiles(f).name(1:end-4);
                            year_events.(key).duration = dur;
                            year_events.(key).avg_current = avg_current;
                            year_events.(key).avg_voltage = mean(V(sIdx:eIdx), 'omitnan');
                            year_events.(key).DCPower_kW = avg_power;  % Use existing DC Power
                            year_events.(key).start_voltage = V(sIdx);
                        end
                    end
                end
            end
        end
    end
                
    all_charge_events.(['year_' year]) = year_events;
    fprintf('Year %s: %d charge events detected (>=30s, avg_current>=%.1fA)\n', year, ev_count, min_avg_current);
end

%% Step 2: Determine power groups from actual data
fprintf('\n=== Step 2: Determine power groups from actual data ===\n');

% Collect all average powers from all years
all_DCPower_kW = [];
for yi = 1:numel(yearList)
    year = yearList{yi};
    E = all_charge_events.(['year_' year]);
    names = fieldnames(E);
    for e = 1:numel(names)
        ev = E.(names{e});
        all_DCPower_kW(end+1) = ev.DCPower_kW;
    end
end

fprintf('Total events: %d\n', numel(all_DCPower_kW));
fprintf('DC Power range: %.3f - %.3f kW\n', min(all_DCPower_kW), max(all_DCPower_kW));
fprintf('DC Power percentiles: 25%%=%.3f, 50%%=%.3f, 75%%=%.3f kW\n', ...
    prctile(all_DCPower_kW, 25), prctile(all_DCPower_kW, 50), prctile(all_DCPower_kW, 75));

% Determine 2 power groups using k-means
K = 2;
[group_assignments, power_centers] = kmeans(all_DCPower_kW(:), K, 'Replicates', 10);
power_centers = sort(power_centers); % Sort from low to high

% Create power ranges based on centers (±10% tolerance)
tol_ratio = 0.10;
power_ranges = zeros(K, 2);
for k = 1:K
    center = power_centers(k);
    power_ranges(k, 1) = center * (1 - tol_ratio);
    power_ranges(k, 2) = center * (1 + tol_ratio);
end

fprintf('DC Power groups determined:\n');
for k = 1:K
    fprintf('  Group %d: center=%.3f kW, range=%.3f-%.3f kW\n', ...
        k, power_centers(k), power_ranges(k,1), power_ranges(k,2));
end

% Create group combinations
n_power_ranges = size(power_ranges, 1);
total_groups = n_power_ranges;  % Only power groups, no voltage sub-groups

fprintf('Start voltage range: %.1f-%.1f V (single range)\n', start_voltage_range(1), start_voltage_range(2));
fprintf('Total groups: %d (power groups only)\n', total_groups);

% Filter events that don't belong to any group
fprintf('\nFiltering events that do not belong to any group...\n');
for yi = 1:numel(yearList)
    year = yearList{yi};
    E = all_charge_events.(['year_' year]);
    names = fieldnames(E);
    nE = numel(names);
    
    if nE == 0, continue; end
    
    filtered_events = struct();
    filtered_count = 0;
    excluded_count = 0;
    
    for e = 1:nE
        ev = E.(names{e});
        avg_curr = ev.avg_current;
        start_volt = ev.start_voltage;
        
        % Check if event belongs to any group
        belongs_to_group = false;
        avg_power = ev.DCPower_kW;
        for pi = 1:n_power_ranges
            if avg_power >= power_ranges(pi,1) && avg_power < power_ranges(pi,2) && ...
               start_volt >= start_voltage_range(1) && start_volt < start_voltage_range(2)
                belongs_to_group = true;
                break;
            end
        end
        
        if belongs_to_group
            filtered_count = filtered_count + 1;
            key = sprintf('Event_%d', filtered_count);
            filtered_events.(key) = ev;
        else
            excluded_count = excluded_count + 1;
        end
    end
    
    all_charge_events.(['year_' year]) = filtered_events;
    fprintf('Year %s: %d events kept, %d events excluded (not in any group)\n', ...
        year, filtered_count, excluded_count);
end

%% Step 3: HI extraction per group with common time range
fprintf('\n=== Step 3: HI extraction per group with common time range ===\n');

yearly_HI = struct();

% First pass: collect events for each group across all years
group_events_all_years = struct();

for pi = 1:n_power_ranges
    group_id = pi;  % Direct mapping: group_id = pi
    group_events_all_years.(['group_' num2str(group_id)]) = struct();
    
    for yi = 1:numel(yearList)
        year = yearList{yi};
        E = all_charge_events.(['year_' year]);
        names = fieldnames(E);
        nE = numel(names);
        
        if nE == 0, continue; end
        
        % Filter events for this power group
        group_events = [];
        for e = 1:nE
            ev = E.(names{e});
            avg_power = ev.DCPower_kW;
            start_volt = ev.start_voltage;
            
            % Check if event matches this group
            if avg_power >= power_ranges(pi,1) && avg_power < power_ranges(pi,2) && ...
               start_volt >= start_voltage_range(1) && start_volt < start_voltage_range(2)
                group_events = [group_events; e];
            end
        end
        
        group_events_all_years.(['group_' num2str(group_id)]).(['year_' year]) = group_events;
    end
end

% Second pass: find common time range and extract HI
for pi = 1:n_power_ranges
    group_id = pi;  % Direct mapping: group_id = pi
        
        % Check if all years have at least 3 events
        all_years_have_events = true;
        min_events = inf;
        for yi = 1:numel(yearList)
            year = yearList{yi};
            if ~isfield(group_events_all_years.(['group_' num2str(group_id)]), ['year_' year])
                all_years_have_events = false;
                break;
            end
            n_events = numel(group_events_all_years.(['group_' num2str(group_id)]).(['year_' year]));
            if n_events < 3
                all_years_have_events = false;
                break;
            end
            min_events = min(min_events, n_events);
        end
        
        if ~all_years_have_events
            fprintf('  Group %d (P%.1f-%.1fW, V%.1f-%.1fV): Not all years have >=3 events, skipping\n', ...
                group_id, power_ranges(pi,1), power_ranges(pi,2), ...
                start_voltage_range(1), start_voltage_range(2));
            continue;
        end
        
        % Find common time range across all years
        common_time_range = inf;
        for yi = 1:numel(yearList)
            year = yearList{yi};
            E = all_charge_events.(['year_' year]);
            names = fieldnames(E);
            group_events = group_events_all_years.(['group_' num2str(group_id)]).(['year_' year]);
            
            for e = 1:numel(group_events)
                ev = E.(names{group_events(e)});
                common_time_range = min(common_time_range, ev.duration);
            end
        end
        
        fprintf('  Group %d (P%.1f-%.1fW, V%.1f-%.1fV): Common time range = %.1f s\n', ...
            group_id, power_ranges(pi,1), power_ranges(pi,2), ...
            start_voltage_range(1), start_voltage_range(2), common_time_range);
        
        % Extract HI for each year using common time range
        for yi = 1:numel(yearList)
            year = yearList{yi};
            E = all_charge_events.(['year_' year]);
            names = fieldnames(E);
            group_events = group_events_all_years.(['group_' num2str(group_id)]).(['year_' year]);
            
            HI_perEvent = [];
            for e = 1:numel(group_events)
                ev = E.(names{group_events(e)});
                
                % Restrict to common time range
                t = ev.time;
                I = ev.current;
                V = ev.voltage;
                
                % Filter to common time range
                time_mask = t <= common_time_range;
                t = t(time_mask);
                I = I(time_mask);
                V = V(time_mask);
                
                if numel(t) < 2, continue; end
                
                % Irregular sampling guard
                dt_vec = [0; diff(t)];
                valid  = isfinite(I) & isfinite(V) & isfinite(dt_vec) & (dt_vec > 0);
                if nnz(valid) < 2, continue; end
                
                % Assign each sample's contribution I*dt to the bin of its V
                binIdx = discretize(V, edges_V);
                sel = valid & ~isnan(binIdx) & (binIdx >= 1) & (binIdx <= nVB);
                if ~any(sel), continue; end
                
                Asec = accumarray(binIdx(sel), I(sel).*dt_vec(sel), [nVB 1], @sum, 0);  % [A*s]
                QmAh = (Asec.'/3600)*1000;  % row vector [1 x nVB], mAh
                
                HI_perEvent = [HI_perEvent; QmAh];
            end
            
            % Aggregate (median) per bin across events
            if ~isempty(HI_perEvent)
                med = median(HI_perEvent, 1, 'omitnan');
                fprintf('    Year %s: %d events, HI extracted\n', year, size(HI_perEvent,1));
                
                yearly_HI.(['year_' year]).(['group_' num2str(group_id)]).HI_med_mAh = med;
                yearly_HI.(['year_' year]).(['group_' num2str(group_id)]).HI_perEvent = HI_perEvent;
                yearly_HI.(['year_' year]).(['group_' num2str(group_id)]).n_events = size(HI_perEvent,1);
                yearly_HI.(['year_' year]).(['group_' num2str(group_id)]).power_range = power_ranges(pi,:);
                yearly_HI.(['year_' year]).(['group_' num2str(group_id)]).voltage_range = start_voltage_range;
                yearly_HI.(['year_' year]).(['group_' num2str(group_id)]).common_time_s = common_time_range;
            end
        end
end

%% Step 4: Visualization
fprintf('\n=== Step 4: Visualization ===\n');

% Figure 1: All individual charge events by year
figure('Position',[100,100,1400,900]);
hold on; grid on; set(gca,'FontSize',Fontsize);
colors = [1 0 0; 0 1 0; 0 0 1];  % Red, Green, Blue for 2021, 2022, 2023
leg = {};

for yi = 1:numel(yearList)
    year = yearList{yi};
    if ~isfield(yearly_HI, ['year_' year]), continue; end
    
    % Collect all individual events for this year
    year_all_events = [];
    for pi = 1:n_power_ranges
        group_id = pi;  % Direct mapping: group_id = pi
            
        if isfield(yearly_HI.(['year_' year]), ['group_' num2str(group_id)])
            G = yearly_HI.(['year_' year]).(['group_' num2str(group_id)]);
            if isfield(G, 'HI_perEvent') && ~isempty(G.HI_perEvent)
                year_all_events = [year_all_events; G.HI_perEvent];
            end
        end
    end
    
    if isempty(year_all_events), continue; end
    
    % Plot each individual event
    x = centers_V;
    for e = 1:size(year_all_events, 1)
        event_data = year_all_events(e, :);
        valid_idx = ~isnan(event_data) & event_data > 0;
        if any(valid_idx)
            plot(x(valid_idx), event_data(valid_idx), '-', 'Color', colors(yi,:), ...
                'LineWidth', 0.5);
        end
    end
    
    leg{end+1} = year;
end

xlabel('Partial voltage range (V)');
ylabel('Partial capacity (mAh)');
title('All individual charge events by year');
legend(leg, 'Location','best');
xlim([edges_V(1), edges_V(end)]);

outFig1 = fullfile(saveDir, 'HI_allIndividualEvents.fig');
saveas(gcf, outFig1);

% Figure 2: Median curves by year (all groups combined)
figure('Position',[200,200,1400,900]);
hold on; grid on; set(gca,'FontSize',Fontsize);
colors = [1 0 0; 0 1 0; 0 0 1];  % Red, Green, Blue for 2021, 2022, 2023
leg = {};

for yi = 1:numel(yearList)
    year = yearList{yi};
    if ~isfield(yearly_HI, ['year_' year]), continue; end
    
    % Collect all groups for this year
    year_data = [];
    for pi = 1:n_power_ranges
        group_id = pi;  % Direct mapping: group_id = pi
            
        if isfield(yearly_HI.(['year_' year]), ['group_' num2str(group_id)])
            G = yearly_HI.(['year_' year]).(['group_' num2str(group_id)]);
            if ~isempty(G.HI_med_mAh) && ~all(isnan(G.HI_med_mAh))
                year_data = [year_data; G.HI_med_mAh];
            end
        end
    end
    
    if isempty(year_data), continue; end
    
    % Calculate median across all groups for this year
    med_all_groups = median(year_data, 1, 'omitnan');
    
    x = centers_V;
    valid_idx = ~isnan(med_all_groups) & med_all_groups > 0;
    if any(valid_idx)
        plot(x(valid_idx), med_all_groups(valid_idx), 'o-', 'Color', colors(yi,:), ...
            'MarkerSize', 6, 'LineWidth', 2, 'MarkerFaceColor', colors(yi,:));
    end
    
    leg{end+1} = year;
end

xlabel('Partial voltage range (V)');
ylabel('Partial capacity (mAh)');
title('Partial capacity by partial voltage range - Median curves');
legend(leg, 'Location','best');
xlim([edges_V(1), edges_V(end)]);

outFig2 = fullfile(saveDir, 'HI_medianCurves.fig');
saveas(gcf, outFig2);

% Save results
save(fullfile(saveDir,'HI_byGroups_Results.mat'), 'yearly_HI','edges_V','centers_V','yearList','power_ranges','start_voltage_range');

fprintf('\nAll done. Results saved to: %s\n', saveDir);

%% --------- Utility: find continuous true segments
function segments = find_continuous_segments(mask)
segments = struct('start',{},'end',{});
if isempty(mask), return; end
dm = diff([false; mask(:); false]);
s  = find(dm == 1);
e  = find(dm == -1) - 1;
for i = 1:numel(s)
    segments(i).start = s(i);
    segments(i).end   = e(i);
end
end