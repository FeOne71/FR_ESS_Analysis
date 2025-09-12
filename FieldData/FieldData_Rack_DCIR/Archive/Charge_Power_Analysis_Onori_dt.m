%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rack Data DCIR Charge - Onori Method with dt=5s Time Window
% Charging event extraction and impedance, power analysis using Onori logic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory - Multiple years
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New';
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\Charge_Power_Analysis_Onori_dt5');
%%
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters
% Battery pack capacity [Ah]
C_nom_rack = 128;               % Ah
P_nom_rack = 112;               % 3.68*64/1000 [kW] | 0.23552 * 2P * 14S * 17S = 112
idle_thr = C_nom_rack * 0.1;   % Initial current threshold (A)
min_duration = 300;             % [s] - Charging duration (300 seconds or more)
max_P_std = P_nom_rack * 0.0001; % Max power standard deviation [kW] 0.0112kW
max_I_std = C_nom_rack * 0.005; % Max current standard deviation [A] 1.28A

% Event similarity filtering parameters
% Fixed voltage bins for charging start voltage
V_bins = 3.5:0.05:4.2;          % [V] - Fixed voltage bins (3.5V ~ 4.2V, 0.1V intervals)
P_bins = 0.05:0.05:10;               % [kW] - Power bins (0.05kW ~ 1kW, 0.05kW intervals)
P_similarity_threshold = 0.003;  % [%] - Power similarity threshold (5%)

% Reference method: Linear fitting of voltage signal (no voltage range analysis)
% We'll use SOC-based analysis like the reference

% Onori parameters for resistance calculation (like reference)
dt = [100 200];                  % Time window: [0.1s 100s] (like reference)

% Sampling time [s]
Ts = 1;

% Plot thinning (reference-like) - downsample plotted points
% Reference used step=100 with Ts=0.01s (~1s 간격). 우리 데이터는 Ts=1s이므로
% 1~5초 간격이 적절. 기본 5초로 설정
plot_step_seconds = 5;   % [s]

%% Initialize data structures
% Year-based data structure: year -> month -> events
yearly_data = struct();
event_counter = 0;

% Store actual event data for plotting
all_events_data = struct();
all_events_data.t = {};
all_events_data.I = {};
all_events_data.V = {};
all_events_data.year = [];
all_events_data.month = [];
all_events_data.R0_instant = {};  % Charging Impedance [mΩ]
all_events_data.DV_DQ_instant = {};  % Voltage-Charge Characteristic [V/Ah]
all_events_data.P_loss_instant = {};
all_events_data.eta_loss_instant = {};
all_events_data.V_start = [];
all_events_data.V_fin = [];
all_events_data.P_avg = [];

% Add SOC-based analysis data (like reference)
all_events_data.SOC_seq = {};              % SOC time series
all_events_data.SOC_R0_avg = {};           % Average R0 for SOC analysis
all_events_data.SOC_P_loss_avg = {};       % Average P_loss for SOC analysis
all_events_data.SOC_eta_loss_avg = {};     % Average eta_loss for SOC analysis

%% Folder traversal - Multiple years
yearDirs = dir(fullfile(dataDir, '20*'));
for y = 1:length(yearDirs)
    if ~yearDirs(y).isdir
        continue;
    end
    yearPath = fullfile(dataDir, yearDirs(y).name);

    % Extract year number
    year_name = yearDirs(y).name;
    yearnum = str2num(year_name);
    fprintf('Processing year: %d\n', yearnum);

    % Initialize yearly data structure
    year_key = sprintf('year_%04d', yearnum);
    yearly_data.(year_key) = struct();
    yearly_data.(year_key).events_count = 0;

    % Get month directories for this year
    monthDirs = dir(fullfile(yearPath, '20*'));
    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir
            continue;
        end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, '*.mat'));

        % Extract month number from folder name (handle different naming formats)
        folder_name = monthDirs(m).name;
        if length(folder_name) >= 6
            monthnum = str2num(folder_name(5:6));
        elseif length(folder_name) >= 2
            monthnum = str2num(folder_name(end-1:end));
        else
            fprintf('Warning: Cannot parse month from folder name: %s\n', folder_name);
            continue;
        end
        fprintf('Processing year %d, month: %d (%d files)\n', yearnum, monthnum, length(matFiles));

        % Initialize monthly data structure within year
        month_key = sprintf('month_%02d', monthnum);
        yearly_data.(year_key).(month_key) = struct();
        yearly_data.(year_key).(month_key).events = struct();
        yearly_data.(year_key).(month_key).events_count = 0;

        for f = 1:length(matFiles)
            fprintf('Processing daily file: %s (%d/%d)\n', matFiles(f).name, f, length(matFiles));
            matFilePath = fullfile(monthPath, matFiles(f).name);
            load(matFilePath);

            % Extract data from Raw structure
            t = Raw.Date_Time_seconds;
            I = Raw.DCCurrent;
            V_raw = Raw.CVavg;
            T_batt = Raw.MTavg;
            soc = Raw.SOC_BMS;

            % Apply moving average to voltage to reduce noise (window=50)
            V = movmean(V_raw, 10);

            % Calculate power using smoothed voltage
            P = V .* I / 1000; % Convert to kW

            %% Event extraction: Idle to Charging transition
            % Find idle periods (low current)
            is_idle = abs(I) < idle_thr;
            is_charging = I > idle_thr;

            % Find transitions from idle to charging
            idle_to_charge = find(is_idle(1:end-1) & is_charging(2:end));

            %% Process each charging event
            for i = 1:length(idle_to_charge)
                start_idx = idle_to_charge(i);
                start_charge_idx = start_idx + 1;

                % Find end of charging period
                end_idx = start_charge_idx;
                while end_idx <= length(is_charging) && is_charging(end_idx)
                    end_idx = end_idx + 1;
                end
                end_idx = end_idx - 3; % 1 second before charging end

                % Condition 1 : Check charging duration
                charge_duration = end_idx - start_charge_idx + 1;
                if charge_duration < min_duration
                    continue;
                end

                % Extract event data (full charging event)
                t_seg_full = t(start_idx:end_idx);
                I_seg_full = I(start_idx:end_idx);
                V_seg_full = V(start_idx:end_idx);
                P_seg_full = P(start_idx:end_idx);
                T_seg_full = T_batt(start_idx:end_idx);
                soc_seg_full = soc(start_idx:end_idx);

                % Remove initial and final transient phases (similar to reference)
                % For 300s events, remove 30s from start and 30s from end
                remove_initial = 30;  % 30 seconds from start
                remove_final = 30;    % 30 seconds from end
                
                if length(t_seg_full) > (remove_initial + remove_final)
                    % Extract stable charging period (excluding transients)
                    t_seg = t_seg_full(remove_initial+1:end-remove_final);
                    I_seg = I_seg_full(remove_initial+1:end-remove_final);
                    V_seg = V_seg_full(remove_initial+1:end-remove_final);
                    P_seg = P_seg_full(remove_initial+1:end-remove_final);
                    T_seg = T_seg_full(remove_initial+1:end-remove_final);
                    soc_seg = soc_seg_full(remove_initial+1:end-remove_final);
                else
                    % If event is too short, use full event
                    t_seg = t_seg_full;
                    I_seg = I_seg_full;
                    V_seg = V_seg_full;
                    P_seg = P_seg_full;
                    T_seg = T_seg_full;
                    soc_seg = soc_seg_full;
                end

                % Linear fitting of voltage signal (like reference)
                % Select voltage points that are monotonically increasing
                clear x_mark;
                clear y_mark;
                k = 1;
                
                for j = 1:length(V_seg)
                    if j == 1
                        x_mark(k) = t_seg(j);
                        y_mark(k) = V_seg(j);
                        k = k + 1;
                    else
                        flag = 1;
                        for xi = 1:k-1
                            if V_seg(j) <= y_mark(xi)
                                flag = 0;
                            end
                        end
                        
                        if flag == 1
                            x_mark(k) = t_seg(j);
                            y_mark(k) = V_seg(j);
                            k = k + 1;
                        end
                    end
                end
                
                % Linear interpolation of voltage signal
                V_fitted = interp1(x_mark, y_mark, t_seg, 'linear', 'extrap');

                % Check stability
                power_std = std(P_seg(3:end-3));
                current_std = std(I_seg(3:end-3));

                % Condition 2 : Power Current std
                if power_std > max_P_std || current_std > max_I_std
                    continue;
                end

                % % Condition 3: SOC diff
                % if abs(max(soc_seg) - min(soc_seg)) < 1 
                %     continue;
                % end                

                %% Calculate R₀ using Onori method (like reference)
                % Method 1: Charging Impedance (Z_CHG) - exactly like reference
                % Use fitted voltage signal like reference
                R0_instant = cell(1, length(dt));

                for k_dt = 1:length(dt)
                    if length(I_seg) > dt(k_dt)
                        R0_instant{k_dt} = zeros(1, length(I_seg)-dt(k_dt));
                        
                        for j = 1:(length(I_seg)-dt(k_dt))
                            V1 = V_fitted(j);
                            V2 = V_fitted(j+dt(k_dt));
                            I_now = I_seg(j);

                            dV = V2 - V1;

                            % Charging Impedance calculation: exactly like reference
                            if I_now ~= 0
                                R0_instant{k_dt}(j) = abs(dV / I_now) * (1000);  
                            else
                                R0_instant{k_dt}(j) = NaN;
                            end
                        end
                    else
                        R0_instant{k_dt} = [];
                    end
                end



                % Method 2: Voltage-Charge Characteristic (DV/DQ) - exactly like reference
                % DV/DQ = R_pseudo / (1000 × dt × Ts) × 3600 [V/Ah]
                DV_DQ_instant = cell(1, length(dt));
                
                for k_dt = 1:length(dt)
                    if ~isempty(R0_instant{k_dt})
                        DV_DQ_instant{k_dt} = R0_instant{k_dt} / (1000 * dt(k_dt) * Ts) * 3600;  % Exactly like reference
                    else
                        DV_DQ_instant{k_dt} = [];
                    end
                end

                % Apply moving average to smooth R₀ values (like reference)
                R0_smoothed = cell(1, length(dt));
                for k_dt = 1:length(dt)
                    if ~isempty(R0_instant{k_dt})
                        R0_smoothed{k_dt} = movmean(R0_instant{k_dt}, 10);  % Reduced window size
                    else
                        R0_smoothed{k_dt} = [];
                    end
                end

                % Check if we have valid R₀ values
                valid_R0 = [];
                for k_dt = 1:length(dt)
                    if ~isempty(R0_smoothed{k_dt})
                        valid_R0 = [valid_R0, R0_smoothed{k_dt}(~isnan(R0_smoothed{k_dt}))];
                    end
                end
                
                if isempty(valid_R0)
                    continue; % Skip events with no valid R₀
                end

                % Use first valid R₀ as representative value for monthly averages
                R0_rep = valid_R0(1);

                %% Calculate power loss and efficiency metrics (time-varying) - like reference
                % Calculate P_loss and η_loss for each time point
                % P_loss(t) = I(t)² × R₀(t), η_loss(t) = R₀(t) / Z(t)
                P_loss_instant = cell(1, length(dt));
                eta_loss_instant = cell(1, length(dt));

                for k_dt = 1:length(dt)
                    if ~isempty(R0_smoothed{k_dt})
                        P_loss_instant{k_dt} = zeros(size(R0_smoothed{k_dt}));
                        eta_loss_instant{k_dt} = zeros(size(R0_smoothed{k_dt}));
                        
                        for j = 1:length(R0_smoothed{k_dt})
                            % Calculate impedance: Z(t) = V_fitted(t)/I(t) (use fitted voltage)
                            Z_instant_k = V_fitted(j) / I_seg(j);

                            % Calculate power loss: P_loss(t) = I(t)² × R₀(t)
                            if ~isnan(R0_smoothed{k_dt}(j))
                                P_loss_instant{k_dt}(j) = (I_seg(j)^2) * R0_smoothed{k_dt}(j) / 1000; % Convert to kW
                            else
                                P_loss_instant{k_dt}(j) = NaN;
                            end

                            % Calculate loss ratio: η_loss(t) = R₀(t) / Z(t)
                            if Z_instant_k > 0 && ~isnan(R0_smoothed{k_dt}(j))
                                eta_loss_instant{k_dt}(j) = R0_smoothed{k_dt}(j) / Z_instant_k;
                            else
                                eta_loss_instant{k_dt}(j) = NaN;
                            end
                        end
                    else
                        P_loss_instant{k_dt} = [];
                        eta_loss_instant{k_dt} = [];
                    end
                end

                % SOC-based analysis (like reference)
                % Calculate average values for SOC analysis
                SOC_R0_avg = [];
                SOC_P_loss_avg = [];
                SOC_eta_loss_avg = [];
                
                for k_dt = 1:length(dt)
                    if ~isempty(R0_smoothed{k_dt})
                        SOC_R0_avg = [SOC_R0_avg, mean(R0_smoothed{k_dt}(~isnan(R0_smoothed{k_dt})))];
                        SOC_P_loss_avg = [SOC_P_loss_avg, mean(P_loss_instant{k_dt}(~isnan(P_loss_instant{k_dt})))];
                        SOC_eta_loss_avg = [SOC_eta_loss_avg, mean(eta_loss_instant{k_dt}(~isnan(eta_loss_instant{k_dt})))];
                    end
                end

                % Get representative values (average of stable charging period)
                % Since we already removed transients, use the full stable period
                V_charging = V_seg;
                I_charging = I_seg;
                P_charging = P_seg;
                T_charging = T_seg;

                V_rep = mean(V_charging);
                I_rep = mean(I_charging);
                T_rep = mean(T_charging);
                P_rep = mean(P_charging);
                
                % Calculate eta_loss_rep from cell structure
                eta_loss_all = [];
                for k_dt = 1:length(dt)
                    if ~isempty(eta_loss_instant{k_dt})
                        eta_loss_all = [eta_loss_all, eta_loss_instant{k_dt}(~isnan(eta_loss_instant{k_dt}))];
                    end
                end
                eta_loss_rep = mean(eta_loss_all);

                % Get stable charging period start voltage, end voltage and average power for similarity filtering
                V1 = V_seg(1);                   % Stable charging period start voltage
                V_fin = V_seg(end);              % Stable charging period end voltage
                P_avg = mean(P_charging(3:end-3));    % Average power (excluding first/last 3 points)

                %% Store results
                event_counter = event_counter + 1;

                % Store event data in yearly_data structure: year -> month -> events
                event_key = sprintf('event_%03d', yearly_data.(year_key).(month_key).events_count + 1);
                yearly_data.(year_key).(month_key).events.(event_key) = struct();

                % Store reference values for R₀ calculation
                yearly_data.(year_key).(month_key).events.(event_key).V1 = V1;  % V(start)
                yearly_data.(year_key).(month_key).events.(event_key).V2 = V_fin;  % V(end)
                % yearly_data.(year_key).(month_key).events.(event_key).I1 = I1;  % I(start)
                yearly_data.(year_key).(month_key).events.(event_key).I2 = I_seg(end);  % I(end)

                % Store time series data
                yearly_data.(year_key).(month_key).events.(event_key).T_seq = T_seg;
                yearly_data.(year_key).(month_key).events.(event_key).V_seq = V_seg;
                yearly_data.(year_key).(month_key).events.(event_key).I_seq = I_seg;
                yearly_data.(year_key).(month_key).events.(event_key).P_seq = P_seg;
                yearly_data.(year_key).(month_key).events.(event_key).R0_seq = R0_smoothed;  % Charging Impedance [mΩ] (cell)
                yearly_data.(year_key).(month_key).events.(event_key).DV_DQ_seq = DV_DQ_instant;  % Voltage-Charge Characteristic [V/Ah] (cell)
                yearly_data.(year_key).(month_key).events.(event_key).P_loss_seq = P_loss_instant;  % Power loss (cell)
                yearly_data.(year_key).(month_key).events.(event_key).eta_loss_seq = eta_loss_instant;  % Loss ratio (cell)

                % Store condition values for similarity filtering
                yearly_data.(year_key).(month_key).events.(event_key).V_start = V1;
                yearly_data.(year_key).(month_key).events.(event_key).V_fin = V_fin;
                yearly_data.(year_key).(month_key).events.(event_key).P_avg_condition = P_avg;

                % Store SOC-based analysis data (like reference)
                yearly_data.(year_key).(month_key).events.(event_key).SOC_seq = soc_seg;
                yearly_data.(year_key).(month_key).events.(event_key).SOC_R0_avg = SOC_R0_avg;
                yearly_data.(year_key).(month_key).events.(event_key).SOC_P_loss_avg = SOC_P_loss_avg;
                yearly_data.(year_key).(month_key).events.(event_key).SOC_eta_loss_avg = SOC_eta_loss_avg;

                % Update event counts
                yearly_data.(year_key).(month_key).events_count = yearly_data.(year_key).(month_key).events_count + 1;
                yearly_data.(year_key).events_count = yearly_data.(year_key).events_count + 1;

                % Store actual event data for plotting (normalize time to start from 0)
                t_normalized = t_seg - t_seg(1); % Start from 0 seconds
                all_events_data.t{end+1} = t_normalized;
                all_events_data.I{end+1} = I_seg;
                all_events_data.V{end+1} = V_seg;
                all_events_data.year = [all_events_data.year, yearnum];
                all_events_data.month = [all_events_data.month, monthnum];
                all_events_data.R0_instant{end+1} = R0_smoothed;  % Charging Impedance [mΩ] (cell-of-cells per dt)
                all_events_data.DV_DQ_instant{end+1} = DV_DQ_instant;  % Voltage-Charge Characteristic [V/Ah] (cell-of-cells per dt)
                all_events_data.P_loss_instant{end+1} = P_loss_instant;  % Power loss (cell-of-cells per dt)
                all_events_data.eta_loss_instant{end+1} = eta_loss_instant;  % Loss ratio (cell-of-cells per dt)
                all_events_data.V_start = [all_events_data.V_start, V1];
                all_events_data.V_fin = [all_events_data.V_fin, V_fin];
                all_events_data.P_avg = [all_events_data.P_avg, P_avg];

                % Store SOC-based analysis data (like reference)
                all_events_data.SOC_seq{end+1} = soc_seg;
                all_events_data.SOC_R0_avg{end+1} = SOC_R0_avg;
                all_events_data.SOC_P_loss_avg{end+1} = SOC_P_loss_avg;
                all_events_data.SOC_eta_loss_avg{end+1} = SOC_eta_loss_avg;


                fprintf('Event %d: R₀=%.4f mΩ, η_loss=%.4f, V_avg=%.2f V\n', ...
                    event_counter, R0_rep, eta_loss_rep, V_rep);
            end
        end



        fprintf('Year %d, Month %d: %d events processed\n', yearnum, monthnum, yearly_data.(year_key).(month_key).events_count);
    end
    fprintf('Year %d: Total %d events processed\n', yearnum, yearly_data.(year_key).events_count);
end

%% Event similarity filtering function - Find common conditions across all years using V_start and V_fin ranges
function [common_conditions, filtered_events_by_year] = find_common_similar_conditions(all_events_data, V_bins, P_bins)
% Find charging conditions that exist across all years with same V_start, V_fin, and P_avg ranges
year_list = unique(all_events_data.year);
n_years = length(year_list);

% Get all V_start, V_fin, and P_avg values
all_V_starts = all_events_data.V_start;
all_V_fins = all_events_data.V_fin;
all_P_avgs = all_events_data.P_avg;

% Find unique condition clusters using V_start and V_fin bins
common_conditions = [];
filtered_events_by_year = struct();

% For each V_start bin, check if similar conditions exist across all years
for v_start_bin = 1:(length(V_bins)-1)
    V_start_bin_start = V_bins(v_start_bin);
    V_start_bin_end = V_bins(v_start_bin + 1);
    V_start_bin_center = (V_start_bin_start + V_start_bin_end) / 2;

    % For each V_fin bin
    for v_fin_bin = 1:(length(V_bins)-1)
        V_fin_bin_start = V_bins(v_fin_bin);
        V_fin_bin_end = V_bins(v_fin_bin + 1);
        V_fin_bin_center = (V_fin_bin_start + V_fin_bin_end) / 2;

        % For each P_avg bin
        for p_bin = 1:(length(P_bins)-1)
            P_bin_start = P_bins(p_bin);
            P_bin_end = P_bins(p_bin + 1);
            P_bin_center = (P_bin_start + P_bin_end) / 2;

            % Find events in this V_start, V_fin, and P_avg bin combination
            events_in_condition = find(all_V_starts >= V_start_bin_start & all_V_starts < V_start_bin_end & ...
                all_V_fins >= V_fin_bin_start & all_V_fins < V_fin_bin_end & ...
                all_P_avgs >= P_bin_start & all_P_avgs < P_bin_end);

            if isempty(events_in_condition)
                continue; % Skip if no events in this condition
            end

            % Check if this condition exists in all years
            condition_exists_in_all_years = true;
            yearly_event_counts = zeros(1, n_years);

            for y = 1:n_years
                yearnum = year_list(y);
                year_events_idx = find(all_events_data.year == yearnum);

                if isempty(year_events_idx)
                    condition_exists_in_all_years = false;
                    break;
                end

                % Get events for this year
                V_starts_year = all_events_data.V_start(year_events_idx);
                V_fins_year = all_events_data.V_fin(year_events_idx);
                P_avgs_year = all_events_data.P_avg(year_events_idx);

                % Find events in this year that fall within this condition
                similar_in_year = [];
                for j = 1:length(year_events_idx)
                    V_start_in_bin = V_starts_year(j) >= V_start_bin_start && V_starts_year(j) < V_start_bin_end;
                    V_fin_in_bin = V_fins_year(j) >= V_fin_bin_start && V_fins_year(j) < V_fin_bin_end;
                    P_avg_in_bin = P_avgs_year(j) >= P_bin_start && P_avgs_year(j) < P_bin_end;

                    if V_start_in_bin && V_fin_in_bin && P_avg_in_bin
                        similar_in_year = [similar_in_year, year_events_idx(j)];
                    end
                end

                yearly_event_counts(y) = length(similar_in_year);

                % Check if this year has at least one event (no minimum requirement)
                if length(similar_in_year) == 0
                    condition_exists_in_all_years = false;
                    break;
                end
            end

            % If condition exists in all years, add to common conditions
            if condition_exists_in_all_years
                condition_idx = length(common_conditions) + 1;
                common_conditions(condition_idx).V_start_bin_start = V_start_bin_start;
                common_conditions(condition_idx).V_start_bin_end = V_start_bin_end;
                common_conditions(condition_idx).V_start_bin_center = V_start_bin_center;
                common_conditions(condition_idx).V_fin_bin_start = V_fin_bin_start;
                common_conditions(condition_idx).V_fin_bin_end = V_fin_bin_end;
                common_conditions(condition_idx).V_fin_bin_center = V_fin_bin_center;
                common_conditions(condition_idx).P_bin_center = P_bin_center;
                common_conditions(condition_idx).events = events_in_condition;
                common_conditions(condition_idx).yearly_event_counts = yearly_event_counts;
                common_conditions(condition_idx).total_events = sum(yearly_event_counts);

                % Organize events by year
                for y = 1:n_years
                    yearnum = year_list(y);
                    year_field = sprintf('year_%d', yearnum);
                    year_events = events_in_condition(all_events_data.year(events_in_condition) == yearnum);
                    filtered_events_by_year(condition_idx).(year_field) = year_events;
                end

            end
        end
    end

    % Sort conditions by total number of events (descending)
    if ~isempty(common_conditions)
        [~, sort_idx] = sort([common_conditions.total_events], 'descend');
        common_conditions = common_conditions(sort_idx);
        filtered_events_by_year = filtered_events_by_year(sort_idx);
    end
end
end

%% Find common similar conditions across all years
[common_conditions, filtered_events_by_year] = find_common_similar_conditions(all_events_data, V_bins, P_bins);

fprintf('\nFound %d common charging conditions across all years:\n', length(common_conditions));
for i = 1:length(common_conditions)
    year_list = unique(all_events_data.year);
    fprintf('Condition %d: V_{start}=%.3f-%.3fV (center=%.3fV), V_{fin}=%.3f-%.3fV (center=%.3fV), P_{avg}=%.1fkW (Total: %d events)\n', ...
        i, common_conditions(i).V_start_bin_start, common_conditions(i).V_start_bin_end, common_conditions(i).V_start_bin_center, ...
        common_conditions(i).V_fin_bin_start, common_conditions(i).V_fin_bin_end, common_conditions(i).V_fin_bin_center, ...
        common_conditions(i).P_bin_center, ...
        common_conditions(i).total_events);
    for y = 1:length(year_list)
        yearnum = year_list(y);
        fprintf('  Year %d: %d events\n', yearnum, common_conditions(i).yearly_event_counts(y));
    end
end

%% Figure 1: Best common charging condition - Time vs Current/Voltage (only if common conditions exist)
if ~isempty(common_conditions)
    % Use only the first (best) common condition (most events)
    cond_idx = 1;
    condition = common_conditions(cond_idx);

    year_list = unique(all_events_data.year);
    n_years = length(year_list);

    % Create separate figures for each year
    for y = 1:n_years
        yearnum = year_list(y);
        year_field = sprintf('year_%d', yearnum);
        year_events_idx = filtered_events_by_year(cond_idx).(year_field);

        if ~isempty(year_events_idx)
            figure('Position', [100, 100, 1400, 800]);

            % Subplot 1: Current vs Time for this year
            subplot(2, 1, 1);
            hold on;
            colors = lines(length(year_events_idx));

            for j = 1:length(year_events_idx)
                event_idx = year_events_idx(j);
                plot(all_events_data.t{event_idx}, all_events_data.I{event_idx}, ...
                    'Color', colors(j, :), 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('Event %d', j));
            end

            xlabel('Time (s)');
            ylabel('Current (A)');
            title(sprintf('Year %d - Common Condition Events: V_{start}=%.3f-%.3fV, V_{fin}=%.3f-%.3fV, P_{avg}=%.3fkW - Current vs Time (n=%d events)', ...
                yearnum, condition.V_start_bin_start, condition.V_start_bin_end, condition.V_fin_bin_start, condition.V_fin_bin_end, condition.P_bin_center, length(year_events_idx)));
            grid on;

            % Subplot 2: Voltage vs Time for this year
            subplot(2, 1, 2);
            hold on;

            for j = 1:length(year_events_idx)
                event_idx = year_events_idx(j);
                plot(all_events_data.t{event_idx}, all_events_data.V{event_idx}, ...
                    'Color', colors(j, :), 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('Event %d', j));
            end

            xlabel('Time (s)');
            ylabel('Voltage (V)');
            title(sprintf('Year %d - Common Condition Events: V_{start}=%.3f-%.3fV, V_{fin}=%.3f-%.3fV, P_{avg}=%.3fkW - Voltage vs Time (n=%d events)', ...
                yearnum, condition.V_start_bin_start, condition.V_start_bin_end, condition.V_fin_bin_start, condition.V_fin_bin_end, condition.P_bin_center, length(year_events_idx)));
            grid on;

            % Save figure for this year
            saveas(gcf, fullfile(saveDir, sprintf('best_common_condition_time_series_year_%d_Onori_dt5.fig', yearnum)));
        end
    end
else
    fprintf('No common charging conditions found across all years. Skipping Figure 1.\n');
end

% %% Figure 2: All common charging conditions - R₀, DV/DQ, P_loss, η_loss vs Time (Yearly basis)
% if ~isempty(common_conditions)
%     % Process all common conditions
%     for cond_idx = 1:length(common_conditions)
%         condition = common_conditions(cond_idx);
% 
%         year_list = unique(all_events_data.year);
%         n_years = length(year_list);
% 
%         figure('Position', [100, 100, 1400, 800]);
% 
%         for y = 1:n_years
%             yearnum = year_list(y);
%             year_field = sprintf('year_%d', yearnum);
%             year_events_idx = filtered_events_by_year(cond_idx).(year_field);
% 
%             if ~isempty(year_events_idx)
%                 % Subplot 1: R₀ (Charging Impedance) vs Time for this year
%                 subplot(4, n_years, y);
%                 hold on;
%                 colors = lines(length(year_events_idx));
%                 sz = 20;  % Marker size
% 
%                 for j = 1:length(year_events_idx)
%                     event_idx = year_events_idx(j);
%                     % R₀ varies with time - plot only valid values (exclude NaN)
%                     t_event = all_events_data.t{event_idx};
%                     R0_event = all_events_data.R0_instant{event_idx};
% 
%                     % Find valid (non-NaN) values
%                     valid_idx = ~isnan(R0_event);
%                     if any(valid_idx)
%                         scatter(t_event(valid_idx), R0_event(valid_idx), sz, 'o', 'filled', ...
%                             'Color', colors(j, :), 'DisplayName', sprintf('Event %d', j));
%                     end
%                 end
% 
%                 xlabel('Time (s)');
%                 ylabel('R₀ (mΩ)');
%                 title(sprintf('Year %d - Charging Impedance vs Time (n=%d)', yearnum, length(year_events_idx)));
%                 grid on;
% 
%                 % Subplot 2: DV/DQ (Voltage-Charge Characteristic) vs Time for this year
%                 subplot(4, n_years, y + n_years);
%                 hold on;
% 
%                 for j = 1:length(year_events_idx)
%                     event_idx = year_events_idx(j);
%                     % DV/DQ varies with time - plot only valid values (exclude NaN)
%                     t_event = all_events_data.t{event_idx};
%                     DV_DQ_event = all_events_data.DV_DQ_instant{event_idx};
% 
%                     % Find valid (non-NaN) values
%                     valid_idx = ~isnan(DV_DQ_event);
%                     if any(valid_idx)
%                         scatter(t_event(valid_idx), DV_DQ_event(valid_idx), sz, 's', 'filled', ...
%                             'Color', colors(j, :), 'DisplayName', sprintf('Event %d', j));
%                     end
%                 end
% 
%                 xlabel('Time (s)');
%                 ylabel('DV/DQ (V/Ah)');
%                 title(sprintf('Year %d - Voltage-Charge Characteristic vs Time (n=%d)', yearnum, length(year_events_idx)));
%                 grid on;
% 
%                 % Subplot 3: P_loss vs Time for this year
%                 subplot(4, n_years, y + 2*n_years);
%                 hold on;
% 
%                 for j = 1:length(year_events_idx)
%                     event_idx = year_events_idx(j);
%                     % P_loss varies with time - plot only valid values (exclude NaN)
%                     t_event = all_events_data.t{event_idx};
%                     P_loss_event = all_events_data.P_loss_instant{event_idx};
% 
%                     % Find valid (non-NaN) values
%                     valid_idx = ~isnan(P_loss_event);
%                     if any(valid_idx)
%                         scatter(t_event(valid_idx), P_loss_event(valid_idx), sz, 'd', 'filled', ...
%                             'Color', colors(j, :), 'DisplayName', sprintf('Event %d', j));
%                     end
%                 end
% 
%                 xlabel('Time (s)');
%                 ylabel('P_{loss} (kW)');
%                 title(sprintf('Year %d - P_{loss} vs Time (n=%d)', yearnum, length(year_events_idx)));
%                 grid on;
% 
%                 % Subplot 4: η_loss vs Time for this year
%                 subplot(4, n_years, y + 3*n_years);
%                 hold on;
% 
%                 for j = 1:length(year_events_idx)
%                     event_idx = year_events_idx(j);
%                     % η_loss varies with time - plot only valid values (exclude NaN)
%                     t_event = all_events_data.t{event_idx};
%                     eta_loss_event = all_events_data.eta_loss_instant{event_idx};
% 
%                     % Find valid (non-NaN) values
%                     valid_idx = ~isnan(eta_loss_event);
%                     if any(valid_idx)
%                         scatter(t_event(valid_idx), eta_loss_event(valid_idx), sz, '^', 'filled', ...
%                             'Color', colors(j, :), 'DisplayName', sprintf('Event %d', j));
%                     end
%                 end
% 
%                 xlabel('Time (s)');
%                 ylabel('η_{loss}');
%                 title(sprintf('Year %d - η_{loss} vs Time (n=%d)', yearnum, length(year_events_idx)));
%                 grid on;
%             end
%         end
% 
%         sgtitle(sprintf('Common Condition %d: V_{start}=%.3f-%.3fV, V_{fin}=%.3f-%.3fV (Total: %d events)', ...
%             cond_idx, condition.V_start_bin_start, condition.V_start_bin_end, condition.V_fin_bin_start, condition.V_fin_bin_end, condition.total_events));
%         saveas(gcf, fullfile(saveDir, sprintf('common_condition_%d_R0_DVDQ_Ploss_etaloss_timeseries_Onori_dt5.fig', cond_idx)));
%     end
% else
%     fprintf('No common charging conditions found across all years. Skipping Figure 2.\n');
% end

%% Figure 3: All common charging conditions - Monthly comparison of R₀, P_loss, η_loss (only if common conditions exist)
if ~isempty(common_conditions)
    % Process all common conditions
    for cond_idx = 1:length(common_conditions)
        condition = common_conditions(cond_idx);

        figure('Position', [100, 100, 1400, 800]);

        % Get year-month combinations for this specific condition only
        year_month_combinations = [];
        for y = 1:length(year_list)
            yearnum = year_list(y);
            year_field = sprintf('year_%d', yearnum);
            if isfield(filtered_events_by_year(cond_idx), year_field)
                year_events_idx = filtered_events_by_year(cond_idx).(year_field);
                for event_idx = year_events_idx
                    year_month_combinations = [year_month_combinations; all_events_data.year(event_idx), all_events_data.month(event_idx)];
                end
            end
        end
        unique_year_months = unique(year_month_combinations, 'rows');
        unique_year_months = sortrows(unique_year_months, [1, 2]); % Sort by year, then month

        n_periods = size(unique_year_months, 1);

        % Initialize arrays to store period averages for this condition
        R0_period_avg = zeros(1, n_periods);
        P_loss_period_avg = zeros(1, n_periods);
        eta_loss_period_avg = zeros(1, n_periods);
        period_event_counts = zeros(1, n_periods);

        % Calculate period averages for this condition
        for p = 1:n_periods
            year_num = unique_year_months(p, 1);
            month_num = unique_year_months(p, 2);
            
            % Get events for this condition and period
            year_field = sprintf('year_%d', year_num);
            if isfield(filtered_events_by_year(cond_idx), year_field)
                period_events_idx = filtered_events_by_year(cond_idx).(year_field);
                % Filter by month
                period_events_idx = period_events_idx(all_events_data.month(period_events_idx) == month_num);
                fprintf('Condition %d, Period %d-%02d: Found %d events\n', cond_idx, year_num, month_num, length(period_events_idx));
            else
                period_events_idx = [];
                fprintf('Condition %d, Period %d-%02d: No events found\n', cond_idx, year_num, month_num);
            end

            if ~isempty(period_events_idx)
                % Calculate representative values for each event in this condition
                R0_this_period = [];
                P_loss_this_period = [];
                eta_loss_this_period = [];

                for event_idx = period_events_idx
                    % Handle cell structure for R0_instant
                    R0_all = [];
                    for k_dt = 1:length(dt)
                        if iscell(all_events_data.R0_instant{event_idx}) && k_dt <= length(all_events_data.R0_instant{event_idx})
                            R0_data = all_events_data.R0_instant{event_idx}{k_dt};
                            if ~isempty(R0_data)
                                R0_all = [R0_all, R0_data(~isnan(R0_data))];
                            end
                        end
                    end
                    R0_rep = mean(R0_all);
                    fprintf('  Event %d: R0_rep = %.4f (from %d valid values)\n', event_idx, R0_rep, length(R0_all));
                    
                    % Handle cell structure for P_loss_instant
                    P_loss_all = [];
                    for k_dt = 1:length(dt)
                        if iscell(all_events_data.P_loss_instant{event_idx}) && k_dt <= length(all_events_data.P_loss_instant{event_idx})
                            P_loss_data = all_events_data.P_loss_instant{event_idx}{k_dt};
                            if ~isempty(P_loss_data)
                                P_loss_all = [P_loss_all, P_loss_data(~isnan(P_loss_data))];
                            end
                        end
                    end
                    P_loss_rep = mean(P_loss_all);
                    
                    % Handle cell structure for eta_loss_instant
                    eta_loss_all = [];
                    for k_dt = 1:length(dt)
                        if iscell(all_events_data.eta_loss_instant{event_idx}) && k_dt <= length(all_events_data.eta_loss_instant{event_idx})
                            eta_loss_data = all_events_data.eta_loss_instant{event_idx}{k_dt};
                            if ~isempty(eta_loss_data)
                                eta_loss_all = [eta_loss_all, eta_loss_data(~isnan(eta_loss_data))];
                            end
                        end
                    end
                    eta_loss_rep = mean(eta_loss_all);

                    if ~isempty(R0_all) && ~isempty(P_loss_all) && ~isempty(eta_loss_all)
                        R0_this_period = [R0_this_period, R0_rep];
                        P_loss_this_period = [P_loss_this_period, P_loss_rep];
                        eta_loss_this_period = [eta_loss_this_period, eta_loss_rep * 100];
                        fprintf('    Added to period averages\n');
                    else
                        fprintf('    Skipped - empty data\n');
                    end
                end

                % Calculate period averages
                if ~isempty(R0_this_period)
                    R0_period_avg(p) = mean(R0_this_period);
                    P_loss_period_avg(p) = mean(P_loss_this_period);
                    eta_loss_period_avg(p) = mean(eta_loss_this_period);
                    period_event_counts(p) = length(R0_this_period);
                    fprintf('  Period %d-%02d: R0_{avg} = %.4f, Ploss_{avg} = %.4f, eta_loss_avg = %.4f (n=%d)\n', ...
                        year_num, month_num, R0_period_avg(p), P_loss_period_avg(p), eta_loss_period_avg(p), period_event_counts(p));
                else
                    R0_period_avg(p) = NaN;
                    P_loss_period_avg(p) = NaN;
                    eta_loss_period_avg(p) = NaN;
                    period_event_counts(p) = 0;
                    fprintf('  Period %d-%02d: No valid data\n', year_num, month_num);
                end
            else
                R0_period_avg(p) = NaN;
                P_loss_period_avg(p) = NaN;
                eta_loss_period_avg(p) = NaN;
                period_event_counts(p) = 0;
            end
        end

        % Create x-axis labels for periods
        period_labels = cell(1, n_periods);
        for p = 1:n_periods
            year_num = unique_year_months(p, 1);
            month_num = unique_year_months(p, 2);
            period_labels{p} = sprintf('%d-%02d', year_num, month_num);
        end

        % Subplot 1: R₀ vs Period
        subplot(1, 3, 1);
        hold on;
        colors = lines(n_periods);

        for p = 1:n_periods
            year_num = unique_year_months(p, 1);
            month_num = unique_year_months(p, 2);
            if ~isnan(R0_period_avg(p))
                plot(p, R0_period_avg(p), 'o-', ...
                    'Color', colors(p, :), 'MarkerSize', 8, 'LineWidth', 0.8, ...
                    'MarkerFaceColor', colors(p, :), ...
                    'DisplayName', sprintf('%d-%02d (n=%d)', year_num, month_num, period_event_counts(p)));
            end
        end

        xlabel('Period (YYYY-MM)');
        ylabel('Average R₀ (mΩ)');
        title(sprintf('Common Condition %d: V_{start}=%.3f-%.3fV, V_{fin}=%.3f-%.3fV (Total: %d events)', ...
            cond_idx, condition.V_start_bin_start, condition.V_start_bin_end, condition.V_fin_bin_start, condition.V_fin_bin_end, condition.total_events));
        grid on;
        xticks(1:n_periods);
        xticklabels(period_labels);
        legend('Location', 'best');

        % Subplot 2: P_loss vs Period
        subplot(1, 3, 2);
        hold on;

        for p = 1:n_periods
            year_num = unique_year_months(p, 1);
            month_num = unique_year_months(p, 2);
            if ~isnan(P_loss_period_avg(p))
                plot(p, P_loss_period_avg(p), 's', ...
                    'Color', colors(p, :), 'MarkerSize', 8, 'LineWidth', 0.8, ...
                    'MarkerFaceColor', colors(p, :), ...
                    'DisplayName', sprintf('%d-%02d (n=%d)', year_num, month_num, period_event_counts(p)));
            end
        end

        xlabel('Period (YYYY-MM)');
        ylabel('Average P_{loss} (kW)');
        title('Period Average P_{loss}');
        grid on;
        xticks(1:n_periods);
        xticklabels(period_labels);
        legend('Location', 'best');

        % Subplot 3: η_loss vs Period
        subplot(1, 3, 3);
        hold on;

        for p = 1:n_periods
            year_num = unique_year_months(p, 1);
            month_num = unique_year_months(p, 2);
            if ~isnan(eta_loss_period_avg(p))
                plot(p, eta_loss_period_avg(p), 'd', ...
                    'Color', colors(p, :), 'MarkerSize', 8, 'LineWidth', 0.8, ...
                    'MarkerFaceColor', colors(p, :), ...
                    'DisplayName', sprintf('%d-%02d (n=%d)', year_num, month_num, period_event_counts(p)));
            end
        end

        xlabel('Period (YYYY-MM)');
        ylabel('Average η_{loss} (%)');
        title('Period Average η_{loss} [Onori dt=5s]');
        grid on;
        xticks(1:n_periods);
        xticklabels(period_labels);
        legend('Location', 'best');

        saveas(gcf, fullfile(saveDir, sprintf('common_condition_%d_period_comparison_Onori_dt5.fig', cond_idx)));
    end
else
    fprintf('No common charging conditions found across all years. Skipping Figure 3.\n');
end

%% Figure 4: Voltage-based Analysis (reference-like) - Z_CHG vs Voltage
if ~isempty(common_conditions)
    % Process all common conditions
    for cond_idx = 1:length(common_conditions)
        condition = common_conditions(cond_idx);

        year_list = unique(all_events_data.year);
        n_years = length(year_list);

        figure('Position', [100, 100, 1600, 800]);

        for y = 1:n_years
            yearnum = year_list(y);
            year_field = sprintf('year_%d', yearnum);
            year_events_idx = filtered_events_by_year(cond_idx).(year_field);

            if ~isempty(year_events_idx)
                fprintf('Figure 4 - Condition %d, Year %d: Found %d events\n', cond_idx, yearnum, length(year_events_idx));
                % Subplot 1: Z_CHG vs Voltage for this year (reference-like)
                subplot(2, n_years, y);
                hold on;
                colors = lines(length(year_events_idx));
                sz = 20;  % Marker size like reference

                % Plot step for thinning
                plot_step = max(1, round(plot_step_seconds / Ts));

                for j = 1:length(year_events_idx)
                    event_idx = year_events_idx(j);
                    V_event = all_events_data.V{event_idx};
                    R0_event = all_events_data.R0_instant{event_idx};  % This is now a cell
                    
                    % Plot only valid (non-NaN) values for each dt
                    if iscell(R0_event)  % Check if it's a cell
                        fprintf('  Event %d: R0_event is cell with %d elements\n', j, length(R0_event));
                        for k_dt = 1:length(dt)
                            if k_dt <= length(R0_event) && ~isempty(R0_event{k_dt})
                                % downsample index mask
                                ds_idx = 1:plot_step:length(R0_event{k_dt});
                                valid_idx = ~isnan(R0_event{k_dt});
                                valid_idx = valid_idx & ismember(1:length(R0_event{k_dt}), ds_idx);
                                if any(valid_idx)
                                    % Use corresponding Voltage values
                                    V_valid = V_event(1:length(R0_event{k_dt}));
                                    scatter(V_valid(valid_idx), R0_event{k_dt}(valid_idx), sz, 'o', 'filled', ...
                                        'Color', colors(j, :), 'DisplayName', sprintf('Event %d (dt=%d)', j, dt(k_dt)));
                                    fprintf('    dt=%d: Plotted %d points\n', dt(k_dt), sum(valid_idx));
                                else
                                    fprintf('    dt=%d: No valid points\n', dt(k_dt));
                                end
                            else
                                fprintf('    dt=%d: Empty or missing data\n', dt(k_dt));
                            end
                        end
                    else
                        fprintf('  Event %d: R0_event is not a cell\n', j);
                    end
                end

                xlabel('Voltage (V)');
                ylabel('Z_{CHG} [mΩ]');
                title(sprintf('Year %d - Z_{CHG} vs Voltage (n=%d)', yearnum, length(year_events_idx)));
                grid on;
                xlim([3.5, 4.2]);

                % Subplot 2: DV/DQ vs Voltage for this year (reference-like)
                subplot(2, n_years, y + n_years);
                hold on;

                for j = 1:length(year_events_idx)
                    event_idx = year_events_idx(j);
                    V_event = all_events_data.V{event_idx};
                    DV_DQ_event = all_events_data.DV_DQ_instant{event_idx};  % This is now a cell
                    
                    % Plot only valid (non-NaN) values for each dt
                    if iscell(DV_DQ_event)  % Check if it's a cell
                        for k_dt = 1:length(dt)
                            if k_dt <= length(DV_DQ_event) && ~isempty(DV_DQ_event{k_dt})
                                ds_idx = 1:plot_step:length(DV_DQ_event{k_dt});
                                valid_idx = ~isnan(DV_DQ_event{k_dt});
                                valid_idx = valid_idx & ismember(1:length(DV_DQ_event{k_dt}), ds_idx);
                                if any(valid_idx)
                                    % Use corresponding Voltage values
                                    V_valid = V_event(1:length(DV_DQ_event{k_dt}));
                                    scatter(V_valid(valid_idx), DV_DQ_event{k_dt}(valid_idx), sz, 's', 'filled', ...
                                        'Color', colors(j, :), 'DisplayName', sprintf('Event %d (dt=%d)', j, dt(k_dt)));
                                end
                            end
                        end
                    end
                end

                xlabel('Voltage (V)');
                ylabel('DV/DQ [V/Ah]');
                title(sprintf('Year %d - DV/DQ vs Voltage (n=%d)', yearnum, length(year_events_idx)));
                grid on;
                xlim([3.5, 4.2]);
            end
        end

        sgtitle(sprintf('Common Condition %d: V_{start}=%.3f-%.3fV, V_{fin}=%.3f-%.3fV, P_{avg}=%.1fkW - Voltage-based Analysis (Total: %d events)', ...
            cond_idx, condition.V_start_bin_start, condition.V_start_bin_end, condition.V_fin_bin_start, condition.V_fin_bin_end, condition.P_bin_center, condition.total_events));
        saveas(gcf, fullfile(saveDir, sprintf('common_condition_%d_VOLTAGE_analysis_Onori_dt5.fig', cond_idx)));
    end
else
    fprintf('No common charging conditions found across all years. Skipping Figure 4.\n');
end

%% Save results
save(fullfile(saveDir, 'power_loss_analysis_results_multi_year_Onori_dt5.mat'), 'yearly_data', 'all_events_data', 'common_conditions', 'filtered_events_by_year');

fprintf('\nAnalysis completed. Total events processed: %d\n', event_counter);
fprintf('Results saved to: %s\n', saveDir);
fprintf('Figures saved:\n');
if ~isempty(common_conditions)
    year_list = unique(all_events_data.year);
    for y = 1:length(year_list)
        yearnum = year_list(y);
        fprintf('  - best_common_condition_time_series_year_%d_Onori_dt5.fig\n', yearnum);
    end
    for cond_idx = 1:length(common_conditions)
        fprintf('  - common_condition_%d_R0_Ploss_etaloss_timeseries_Onori_dt5.fig\n', cond_idx);
        fprintf('  - common_condition_%d_period_comparison_Onori_dt5.fig\n', cond_idx);
        fprintf('  - common_condition_%d_SOC_analysis_Onori_dt5.fig\n', cond_idx);
    end
else
    fprintf('  - No common conditions found, no figures generated\n');
end
