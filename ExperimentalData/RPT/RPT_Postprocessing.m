%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Postprocessing - Generate and Save Figures by Category
% Chg/Dch OCV > Avg OCV > Iterate 8Ch 
% Static Cpacity / OCV Curve integration scipt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

% Parallel processing: set true if Parallel Computing Toolbox available
useParallel = true;  % set true to use parfor over channels
if useParallel
    try
        pool = gcp('nocreate');
        if isempty(pool)
            fprintf('Starting parallel pool (may take 30-60 s)...\n');
            nWorkers = min(6, feature('numcores'));  % many clusters default to max 6
            try
                parpool('local', nWorkers);
            catch ME1
                % Try 'Threads' profile if 'local' (Processes) fails
                fprintf('Local pool failed, trying Threads profile...\n');
                fprintf('  Reason: %s\n', ME1.message);
                parpool('threads', nWorkers);
            end
        end
        pool = gcp('nocreate');
        fprintf('Using parallel pool (%s, %d workers).\n', pool.Cluster.Profile, pool.NumWorkers);
        % Attach so workers can run process_one_channel_vgrid (it has local_interp_vq/vqt inside)
        scriptDir = fileparts(mfilename('fullpath'));
        addAttachedFiles(pool, {fullfile(scriptDir, 'process_one_channel_vgrid.m')});
    catch ME
        useParallel = false;
        fprintf('Parallel pool not available, using serial.\n');
        fprintf('  Reason: %s\n', ME.message);
    end
end

%% File Directory
% ExperimentalDataPath = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data\RPT';
ExperimentalDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data\RPT';
saveDir_OCV   = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';
saveDir_Crate = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Crate_integrated';
if ~exist(saveDir_OCV,'dir'); mkdir(saveDir_OCV); end
if ~exist(saveDir_Crate,'dir'); mkdir(saveDir_Crate); end

channels = {'Ch09', 'Ch10', 'Ch11', 'Ch12', 'Ch13', 'Ch14', 'Ch15', 'Ch16'};
rpt_cycles = {'0cyc', '200cyc', '400cyc','600cyc','800cyc','1000cyc'};

%% OCV Processing
fprintf('\n=== OCV Processing ===\n');

% OCV conditions: charge (8,2), discharge (10,2)
ocv_conditions = {'charge', 'discharge'};
ocv_steps = [8, 10];

% Static Capacity conditions: discharge (3,2) only
static_capacity_step = 3;  % Step 3: Discharge

% Structure to store OCV data for each channel
channel_ocv_data = struct();
% Structure to store Static Capacity data for each channel
channel_static_capacity_data = struct();
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    channel_ocv_data.(channel) = struct();
    channel_static_capacity_data.(channel) = struct();
end

% Pre-allocate for parallel gather
ch_ocv_results = cell(length(channels), 1);
ch_static_results = cell(length(channels), 1);

if useParallel
    parfor ch_idx = 1:length(channels)
        [ch_ocv_results{ch_idx}, ch_static_results{ch_idx}] = process_one_channel_ocv( ...
            channels{ch_idx}, rpt_cycles, ExperimentalDataPath, saveDir_OCV, ...
            ocv_conditions, ocv_steps, static_capacity_step);
    end
else
    for ch_idx = 1:length(channels)
        [ch_ocv_results{ch_idx}, ch_static_results{ch_idx}] = process_one_channel_ocv( ...
            channels{ch_idx}, rpt_cycles, ExperimentalDataPath, saveDir_OCV, ...
            ocv_conditions, ocv_steps, static_capacity_step);
    end
end

% Gather into structs
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    channel_ocv_data.(channel) = ch_ocv_results{ch_idx};
    channel_static_capacity_data.(channel) = ch_static_results{ch_idx};
    if ~useParallel
        fprintf('Saved: %s\n', fullfile(saveDir_OCV, sprintf('%s_OCV.fig', channel)));
    end
end
if useParallel
    fprintf('Saved OCV figures: %s\n', saveDir_OCV);
end

fprintf('Building Average OCV by Cycle...\n');
% Figure 2: Average OCV for each cycle separately (Approach A)
figure('Name', 'Average OCV by Cycle', 'Position', [100 100 1200 800]);
hold on;

% Collect per-channel averaged OCV (charge/discharge averaged per channel) - Dynamic
per_channel_avg_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    per_channel_avg_data.(field_name) = [];
end

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
        % 데이터가 있는 경우에만 추가
        if isfield(channel_ocv_data, channel) && isfield(channel_ocv_data.(channel), field_name)
            per_channel_avg_data.(field_name) = [per_channel_avg_data.(field_name); channel_ocv_data.(channel).(field_name).avg_ocv];
        end
    end
end

% Average across channels - Dynamic
avg_ocv_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    avg_ocv_data.(field_name) = mean(per_channel_avg_data.(field_name), 1, 'omitnan');
end

soc_grid = 0:0.5:100;

% Sorted grid (already ascending, but keep consistent)
[soc_grid_sorted, sort_idx] = sort(soc_grid);

% Sort OCV data and plot dynamically
colors = {'g-', 'b-', 'y-', 'r-', 'm-', 'c-', 'k-'}; % Add more colors if needed
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    avg_ocv_sorted = avg_ocv_data.(field_name)(sort_idx);
    
    color_idx = mod(rpt_idx-1, length(colors)) + 1;
    plot(soc_grid_sorted, avg_ocv_sorted, colors{color_idx}, 'LineWidth', 3, 'DisplayName', sprintf('RPT%s Average', rpt_cycle(1:end-3)));
end

xlabel('SOC [%]');
ylabel('Voltage [V]');
title('Average OCV by Cycle (Approach A: per-channel avg → channel avg)');
legend('Location', 'best');
grid on;
xlim([0 100]);

% Save average OCV (compact = faster/smaller)
figName = fullfile(saveDir_OCV, 'Average_OCV_by_Cycle.fig');
fprintf('Saving %s...\n', 'Average_OCV_by_Cycle.fig');
savefig(gcf, figName);
close(gcf);
fprintf('Saved: %s\n', figName);

% Create OCV functions with sorted data - Dynamic
OCV_func_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    avg_ocv_sorted = avg_ocv_data.(field_name)(sort_idx);
    OCV_func_data.(field_name) = @(soc_query) interp1(soc_grid_sorted, avg_ocv_sorted, soc_query, 'linear');
end


% Compute mean capacity across channels - Dynamic
cap_list_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    cap_list_data.(field_name) = [];
end

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
        % 데이터가 있는 경우에만 추가
        if isfield(channel_ocv_data, channel) && isfield(channel_ocv_data.(channel), field_name)
            cap_list_data.(field_name) = [cap_list_data.(field_name), channel_ocv_data.(channel).(field_name).capacity];
        end
    end
end

mean_capacity_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    mean_capacity_data.(field_name) = mean(cap_list_data.(field_name), 'omitnan');
end

% Compute Static Capacity statistics across channels - Dynamic (Discharge only)
static_capacity_list_dch = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    static_capacity_list_dch.(field_name) = [];
end

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
        % 데이터가 있는 경우에만 추가
        if isfield(channel_static_capacity_data, channel) && isfield(channel_static_capacity_data.(channel), field_name)
            static_capacity_list_dch.(field_name) = [static_capacity_list_dch.(field_name), channel_static_capacity_data.(channel).(field_name).discharge];
        end
    end
end

mean_static_capacity_dch = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    mean_static_capacity_dch.(field_name) = mean(static_capacity_list_dch.(field_name), 'omitnan');
end


% OCV as a function of charge amount Q (Ah) - Dynamic
q_grid_data = struct();
OCV_Q_func_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    q_grid_data.(field_name) = mean_capacity_data.(field_name) .* (soc_grid_sorted ./ 100);
    avg_ocv_sorted = avg_ocv_data.(field_name)(sort_idx);
    OCV_Q_func_data.(field_name) = @(q_query) interp1(q_grid_data.(field_name), avg_ocv_sorted, q_query, 'linear');
end

fprintf('Building OCV_data struct...\n');
% Save OCV arrays and metadata to MAT for later use - Dynamic (soc/ocv grid removed; use RPT_VQ_grid for V-grid data)
OCV_data = struct();

% Add dynamic data for each cycle
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    
    % Add OCV data
    OCV_data.(sprintf('avg_ocv_rpt%s', rpt_cycle(1:end-3))) = avg_ocv_data.(field_name)(sort_idx);
    OCV_data.(sprintf('mean_capacity_rpt%s', rpt_cycle(1:end-3))) = mean_capacity_data.(field_name);  % OCV capacity (8-channel average)
    OCV_data.(sprintf('q_grid_rpt%s', rpt_cycle(1:end-3))) = q_grid_data.(field_name);
    OCV_data.(sprintf('OCV_SOC_func_rpt%s', rpt_cycle(1:end-3))) = OCV_func_data.(field_name);
    OCV_data.(sprintf('OCV_Q_func_rpt%s', rpt_cycle(1:end-3))) = OCV_Q_func_data.(field_name);
    
    % Add Static Capacity data (8-channel average, discharge only)
    OCV_data.(sprintf('mean_static_capacity_rpt%s', rpt_cycle(1:end-3))) = mean_static_capacity_dch.(field_name);
    
    % Add individual channel Static Capacity data (discharge only)
    for ch_idx = 1:length(channels)
        channel = channels{ch_idx};
        channel_num = channel(3:end);  % 'Ch09' -> '09'
        % 데이터가 있는 경우에만 추가
        if isfield(channel_static_capacity_data, channel) && isfield(channel_static_capacity_data.(channel), field_name)
            OCV_data.(sprintf('static_capacity_ch%s_rpt%s', channel_num, rpt_cycle(1:end-3))) = channel_static_capacity_data.(channel).(field_name).discharge;
        else
            OCV_data.(sprintf('static_capacity_ch%s_rpt%s', channel_num, rpt_cycle(1:end-3))) = NaN;
        end
    end
end

%% 8개 채널 통합 OCV 함수 생성 (OCV_integrated.m 로직)
fprintf('\n=== Creating Integrated OCV Functions ===\n');

% 8개 채널의 OCV를 8x101 행렬로 쌓기 - Dynamic
all_V_OCV_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    all_V_OCV_data.(field_name) = [];
end

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
        % 데이터가 있는 경우에만 추가
        if isfield(channel_ocv_data, channel) && isfield(channel_ocv_data.(channel), field_name)
            all_V_OCV_data.(field_name) = [all_V_OCV_data.(field_name); channel_ocv_data.(channel).(field_name).avg_ocv];
        end
    end
end

% 8채널 평균 OCV 계산 - Dynamic
V_avg_SOC_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    V_avg_SOC_data.(field_name) = mean(all_V_OCV_data.(field_name), 1, 'omitnan'); % 1x101
end


% 통합 OCV 함수 생성 - Dynamic
OCV_SOC_func_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    OCV_SOC_func_data.(field_name) = @(soc_query) interp1(soc_grid_sorted, V_avg_SOC_data.(field_name), soc_query, 'linear');
end


% 통합 구조체 생성 - Dynamic
OCV_integrated_data = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    
    OCV_integrated_data.(field_name).V_avg_SOC = V_avg_SOC_data.(field_name);
    OCV_integrated_data.(field_name).OCV_SOC_func = OCV_SOC_func_data.(field_name);
    OCV_integrated_data.(field_name).mean_capacity = mean_capacity_data.(field_name);
    
    % 개별 셀 용량 저장 (8개 채널)
    individual_ocv_capacity = [];
    individual_static_capacity = [];
    for ch_idx = 1:length(channels)
        channel = channels{ch_idx};
        % 데이터가 있는 경우에만 추가
        if isfield(channel_ocv_data, channel) && isfield(channel_ocv_data.(channel), field_name)
            individual_ocv_capacity = [individual_ocv_capacity, channel_ocv_data.(channel).(field_name).capacity];
        end
        if isfield(channel_static_capacity_data, channel) && isfield(channel_static_capacity_data.(channel), field_name)
            individual_static_capacity = [individual_static_capacity, channel_static_capacity_data.(channel).(field_name).discharge];
        end
    end
    OCV_integrated_data.(field_name).individual_ocv_capacity = individual_ocv_capacity;  % 1x8 array
    OCV_integrated_data.(field_name).individual_static_capacity = individual_static_capacity;  % 1x8 array
    
    % T1(℃): 채널별 OCV 구간 시계열 (charge step + discharge step 순서)
    OCV_integrated_data.(field_name).T1 = struct();
    for ch_idx = 1:length(channels)
        channel = channels{ch_idx};
        filepath = fullfile(ExperimentalDataPath, sprintf('%s_RPT_%s.csv', channel, rpt_cycle));
        if ~isfile(filepath)
            OCV_integrated_data.(field_name).T1.(channel) = [];
            continue;
        end
        try
            Tc = readtable(filepath, 'VariableNamingRule', 'preserve');
        catch
            OCV_integrated_data.(field_name).T1.(channel) = [];
            continue;
        end
        t1Col = [];
        for j = 1:numel(Tc.Properties.VariableNames)
            if contains(lower(char(Tc.Properties.VariableNames{j})), 't1')
                t1Col = j; break;
            end
        end
        if isempty(t1Col)
            OCV_integrated_data.(field_name).T1.(channel) = [];
            continue;
        end
        ch_mask = (Tc{:,2} == ocv_steps(1)) & (Tc{:,4} == 2);
        dch_mask = (Tc{:,2} == ocv_steps(2)) & (Tc{:,4} == 2);
        t1_vals = [Tc{ch_mask, t1Col}; Tc{dch_mask, t1Col}];
        OCV_integrated_data.(field_name).T1.(channel) = double(t1_vals(:));  % 시계열 [℃]
    end
end

% 통합 OCV 시각화 - Dynamic
figure('Name', 'Integrated SOC-OCV Curve', 'Position', [100 100 1200 800]);
hold on;
colors = {'bo-', 'go-', 'yo-', 'ro-', 'mo-', 'co-', 'ko-'}; % Add more colors if needed
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    color_idx = mod(rpt_idx-1, length(colors)) + 1;
    plot(soc_grid_sorted, V_avg_SOC_data.(field_name), colors{color_idx}, 'LineWidth', 2, 'DisplayName', sprintf('RPT%s (%s)', rpt_cycle(1:end-3), rpt_cycle));
end
xlabel('SOC [%]');
ylabel('OCV [V]');
title('Integrated SOC-OCV Curve (8-cell Average)');
grid on;
xlim([0 100]);
legend('Location', 'best');

% 통합 그래프 저장 (compact = faster/smaller)
figName = fullfile(saveDir_OCV, 'Integrated_OCV_8cell_Average.fig');
fprintf('Saving Integrated_OCV_8cell_Average.fig...\n');
savefig(gcf, figName);
close(gcf);
fprintf('Saved: %s\n', figName);

% 통합 데이터를 OCV_data에 추가 - Dynamic
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    OCV_data.(sprintf('OCV_integrated_%s', rpt_cycle(1:end-3))) = OCV_integrated_data.(field_name);
end


matSaveFile = fullfile(saveDir_OCV, 'OCV_integrated.mat');
fprintf('Saving OCV_integrated.mat...\n');
save(matSaveFile, 'OCV_data', '-v7');
fprintf('Saved OCV data MAT: %s\n', matSaveFile);

%% ========================================================================
%  C-rate Processing (0.1C/0.5C/1C/2C/3C)
%  - Charge/Discharge step idx: [28/48, 32/52, 36/56, 40/60, 44/64]
%  - Save Raw Voltage, Current, Capacity data (No interpolation)
%  - Save Channel-wise data only (Removed 8-channel average)
% =========================================================================

fprintf('\n=== C-rate Processing ===\n');

crate_labels = {'c01','c05','c1','c2','c3'};
crate_steps_charge = [28 32 36 40 44];      %0.1C 0.5C 1C 2C 3C
crate_steps_discharge = [48 52 56 60 64];   %0.1C 0.5C 1C 2C 3C
% crate_grid_points = 200; % Removed

Crate_data = struct();
crate_cell = cell(length(channels), 1);

% 채널별 처리 (병렬 가능)
if useParallel
    parfor ch_idx = 1:length(channels)
        crate_cell{ch_idx} = process_one_channel_crate(ch_idx, channels, rpt_cycles, ...
            ExperimentalDataPath, crate_labels, crate_steps_charge, crate_steps_discharge);
    end
else
    for ch_idx = 1:length(channels)
        crate_cell{ch_idx} = process_one_channel_crate(ch_idx, channels, rpt_cycles, ...
            ExperimentalDataPath, crate_labels, crate_steps_charge, crate_steps_discharge);
    end
end
for ch_idx = 1:length(channels)
    Crate_data.(sprintf('Ch%s', channels{ch_idx}(3:end))) = crate_cell{ch_idx};
end

% Removed 8-channel averaging (total8) as it relied on common grid interpolation which is now removed.

% Save C-rate data
matSaveFile_Crate = fullfile(saveDir_Crate, 'Crate_integrated.mat');
fprintf('Saving Crate_integrated.mat...\n');
save(matSaveFile_Crate, 'Crate_data', '-v7');
fprintf('Saved C-rate data MAT: %s\n', matSaveFile_Crate);

%% ========================================================================
%  V-grid Interpolation (0.001 V spacing, no extrapolation)
%  Structure: cyc0 -> Ch09..Ch16, cyc200 -> Ch09..Ch16, ...
%  Each channel: Static, OCV_charge, OCV_discharge, c01_charge, c01_discharge, ... c3_discharge
%  interp1(V, Q, V_grid, 'linear'); V_grid = min(V):0.001:max(V)
% =========================================================================
fprintf('\n=== V-grid Interpolation (Static, OCV, C-rate) ===\n');
dV = 0.001;  % V-grid spacing [V]

RPT_VQ_grid = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    cycle_key = sprintf('cyc%s', rpt_cycle(1:end-3));
    RPT_VQ_grid.(cycle_key) = struct();
    ch_grid_cell = cell(length(channels), 1);
    if useParallel
        parfor ch_idx = 1:length(channels)
            ch_grid_cell{ch_idx} = process_one_channel_vgrid(channels{ch_idx}, cycle_key, rpt_cycle, ...
                ExperimentalDataPath, Crate_data, static_capacity_step, ocv_steps, crate_labels, dV);
        end
    else
        for ch_idx = 1:length(channels)
            ch_grid_cell{ch_idx} = process_one_channel_vgrid(channels{ch_idx}, cycle_key, rpt_cycle, ...
                ExperimentalDataPath, Crate_data, static_capacity_step, ocv_steps, crate_labels, dV);
        end
    end
    for ch_idx = 1:length(channels)
        RPT_VQ_grid.(cycle_key).(channels{ch_idx}) = ch_grid_cell{ch_idx};
    end
end

% T1(℃)는 V,I,Q와 동일한 곳에서 읽어 저장됨: process_one_channel_vgrid(Static/OCV), process_one_channel_crate(C-rate) → RPT_VQ_grid 각 구간에 T1_raw, T1 포함

% Save RPT_VQ_grid (V-grid interpolated, no extrapolation)
matSaveFile_VQ = fullfile(saveDir_OCV, 'RPT_VQ_grid.mat');
fprintf('Saving RPT_VQ_grid.mat...\n');
save(matSaveFile_VQ, 'RPT_VQ_grid', '-v7');
fprintf('Saved RPT_VQ_grid MAT: %s\n', matSaveFile_VQ);

%% ========================================================================
%  C-rate Visualization (Capacity vs Voltage)
%  - 채널별: 5개 C-rate를 서브플롯(2x3)으로 구성
%  - 사이클별 색상 구분
%  - Charge / Discharge 각각 별도 fig 저장
% =========================================================================

fprintf('\n=== C-rate Visualization ===\n');
rate_names = {'0.1C','0.5C','1C','2C','3C'};
rate_labels = crate_labels;  % {'c01','c05','c1','c2','c3'}

for ch_idx = 1:length(channels)
    fprintf('  C-rate figures: %s\n', channels{ch_idx});
    channel = channels{ch_idx};
    channel_key = sprintf('Ch%s', channel(3:end));
    
    % Charge figure
    fig_chg = figure('Name', sprintf('Crate Charge - %s', channel), ...
        'Position', [100 100 1400 900], 'Color', 'w');
    tl_chg = tiledlayout(fig_chg, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    colors = lines(length(rpt_cycles));
    
    for r = 1:length(rate_labels)
        label = rate_labels{r};
        nexttile(tl_chg, r); hold on; grid on;
        for cycIdx = 1:length(rpt_cycles)
            cycle = rpt_cycles{cycIdx};
            cycle_key = sprintf('cyc%s', cycle(1:end-3));
            if isfield(Crate_data, channel_key) && ...
               isfield(Crate_data.(channel_key), cycle_key) && ...
               isfield(Crate_data.(channel_key).(cycle_key), label) && ...
               isfield(Crate_data.(channel_key).(cycle_key).(label), 'charge') && ...
               ~isempty(Crate_data.(channel_key).(cycle_key).(label).charge.Q)
           
                Q = Crate_data.(channel_key).(cycle_key).(label).charge.Q;
                V = Crate_data.(channel_key).(cycle_key).(label).charge.V;
                % I = Crate_data.(channel_key).(cycle_key).(label).charge.I; 
                
                if ~isempty(Q)
                    plot(Q, V, '-', 'LineWidth', 1.5, ...
                        'Color', colors(cycIdx,:), 'DisplayName', cycle);
                end
            end
        end
        title(sprintf('%s Charge', rate_names{r}));
        xlabel('Capacity [Ah]');
        ylabel('Voltage [V]');
        legend('Location', 'best', 'FontSize', 8);
    end
    savefig(fig_chg, fullfile(saveDir_Crate, sprintf('%s_Crate_Charge.fig', channel)));
    close(fig_chg);
    
    % Discharge figure
    fig_dch = figure('Name', sprintf('Crate Discharge - %s', channel), ...
        'Position', [100 100 1400 900], 'Color', 'w');
    tl_dch = tiledlayout(fig_dch, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for r = 1:length(rate_labels)
        label = rate_labels{r};
        nexttile(tl_dch, r); hold on; grid on;
        for cycIdx = 1:length(rpt_cycles)
            cycle = rpt_cycles{cycIdx};
            cycle_key = sprintf('cyc%s', cycle(1:end-3));
            if isfield(Crate_data, channel_key) && ...
               isfield(Crate_data.(channel_key), cycle_key) && ...
               isfield(Crate_data.(channel_key).(cycle_key), label) && ...
               isfield(Crate_data.(channel_key).(cycle_key).(label), 'discharge') && ...
               ~isempty(Crate_data.(channel_key).(cycle_key).(label).discharge.Q)
           
                Q = Crate_data.(channel_key).(cycle_key).(label).discharge.Q;
                V = Crate_data.(channel_key).(cycle_key).(label).discharge.V;
                % I = Crate_data.(channel_key).(cycle_key).(label).discharge.I;
                
                if ~isempty(Q)
                    plot(Q, V, '-', 'LineWidth', 1.5, ...
                        'Color', colors(cycIdx,:), 'DisplayName', cycle);
                end
            end
        end
        title(sprintf('%s Discharge', rate_names{r}));
        xlabel('Capacity [Ah]');
        ylabel('Voltage [V]');
        legend('Location', 'best', 'FontSize', 8);
    end
    savefig(fig_dch, fullfile(saveDir_Crate, sprintf('%s_Crate_Discharge.fig', channel)));
    close(fig_dch);
end

fprintf('\n=== RPT Postprocessing Completed ===\n');
