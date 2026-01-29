%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Postprocessing - Generate and Save Figures by Category
% Chg/Dch OCV > Avg OCV > Iterate 8Ch 
% Static Cpacity / OCV Curve integration scipt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

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

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};

    % Figure 1: Individual channel OCV
    figure('Name', sprintf('OCV - %s', channel), 'Position', [100 100 1200 800]);
    set(gcf, 'Visible', 'off');

    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
        filepath = fullfile(ExperimentalDataPath, filename);

        % Ch14 RPT800cyc는 실험 중단으로 인해 방전 OCV 데이터 이상 - 건너뛰기
        if strcmp(channel, 'Ch14') && strcmp(rpt_cycle, '800cyc')
            fprintf('Skipping: Ch14 RPT800cyc (experiment interrupted, discharge OCV data invalid)\n');
            continue;
        end

        % 파일이 없거나 읽을 수 없으면 건너뛰기
        if ~isfile(filepath)
            fprintf('Warning: File not found - %s (skipping)\n', filename);
            continue;
        end
        
        try
            T = readtable(filepath);
        catch ME
            fprintf('Warning: Unable to read file - %s (skipping)\n', filename);
            fprintf('  Error: %s\n', ME.message);
            continue;
        end

        subplot(5,2,rpt_idx);
        hold on;

        charge_ocv = [];
        discharge_ocv = [];

        for ocv_idx = 1:length(ocv_conditions)
            step_idx = ocv_steps(ocv_idx);
            ocv_data = (T{:,2} == step_idx) & (T{:,4} == 2);

            capacity_data = T{ocv_data, 9};
            voltage_data = T{ocv_data, 8};

            % Reverse discharge OCV so voltage increases with SOC
            if strcmp(ocv_conditions{ocv_idx}, 'discharge')
                voltage_data = flipud(voltage_data);
            end

            % Move avg (window size = data 수 / 200)
            window_size = round(length(voltage_data) / 200);
            if window_size < 1; window_size = 1; end
            voltage_smoothed = movmean(voltage_data, window_size);

            % 200개 간격으로 추출 (x 변수 grid 201)
            soc_data = linspace(0, 100, 201);
            
            % interp로 y 추출
            voltage_interp = interp1(linspace(0, 100, length(voltage_smoothed)), voltage_smoothed, soc_data, 'linear');

            if strcmp(ocv_conditions{ocv_idx}, 'charge')
                charge_ocv = voltage_interp;
                plot(soc_data, voltage_interp, 'b-', 'LineWidth', 2, 'DisplayName', 'Charge');
            else
                discharge_ocv = voltage_interp;
                plot(soc_data, voltage_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'Discharge');
            end
        end

        % Calculate and display average charge/discharge OCV
        avg_ocv = (charge_ocv + discharge_ocv) / 2;
        plot(soc_data, avg_ocv, 'g-', 'LineWidth', 2, 'DisplayName', 'Average');

        % Dynamic Field name mapping
        field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % '0cyc' → 'cyc0', '200cyc' → 'cyc200', etc.

        % Save soc, avg_ocv, capacity (use last capacity of each step)
        charge_mask = (T{:,2} == ocv_steps(1)) & (T{:,4} == 2);
        discharge_mask = (T{:,2} == ocv_steps(2)) & (T{:,4} == 2);
        channel_ocv_data.(channel).(field_name).soc = soc_data;
        channel_ocv_data.(channel).(field_name).avg_ocv = avg_ocv;
        channel_ocv_data.(channel).(field_name).capacity = (T{find(charge_mask,1,'last'), 9} + T{find(discharge_mask,1,'last'), 9})/2;
        
        % Extract Static Capacity (Step 3: Discharge only)
        static_dch_mask = (T{:,2} == static_capacity_step) & (T{:,4} == 2);
        if any(static_dch_mask)
            static_capacity_dch = T{find(static_dch_mask,1,'last'), 9};  % Last capacity of discharge step
            channel_static_capacity_data.(channel).(field_name).discharge = static_capacity_dch;
        else
            channel_static_capacity_data.(channel).(field_name).discharge = NaN;
            fprintf('Warning: Static Capacity (discharge) data not found for %s - %s\n', channel, rpt_cycle);
        end

        xlabel('SOC [%]');
        ylabel('Voltage [V]');
        title(sprintf('%s - %s OCV', channel, rpt_cycle));
        legend('Location', 'best');
        grid on;
        xlim([0 100]);
    end

    % Save figure
    figName = fullfile(saveDir_OCV, sprintf('%s_OCV.fig', channel));
    savefig(gcf, figName);
    close(gcf);  % 그래프 창 닫기
    fprintf('Saved: %s\n', figName);
end

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

% Save average OCV
figName = fullfile(saveDir_OCV, 'Average_OCV_by_Cycle.fig');
savefig(gcf, figName);
close(gcf);  % 그래프 창 닫기
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

% Save OCV arrays and metadata to MAT for later use - Dynamic
OCV_data = struct();
OCV_data.soc_grid = soc_grid_sorted;

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
    
    OCV_integrated_data.(field_name).SOC_grid = soc_grid_sorted;
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

% 통합 그래프 저장
figName = fullfile(saveDir_OCV, 'Integrated_OCV_8cell_Average.fig');
savefig(gcf, figName);
close(gcf);  % 그래프 창 닫기
fprintf('Saved: %s\n', figName);

% 통합 데이터를 OCV_data에 추가 - Dynamic
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3)); % Remove 'cyc' suffix
    OCV_data.(sprintf('OCV_integrated_%s', rpt_cycle(1:end-3))) = OCV_integrated_data.(field_name);
end


matSaveFile = fullfile(saveDir_OCV, 'OCV_integrated.mat');
save(matSaveFile, 'OCV_data');
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

% 채널별 처리
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    channel_key = sprintf('Ch%s', channel(3:end));
    Crate_data.(channel_key) = struct();
    
    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        cycle_key = sprintf('cyc%s', rpt_cycle(1:end-3)); % 0cyc -> cyc0 등
        
        % Ch14 RPT800cyc 스킵 (실험 중단)
        if strcmp(channel, 'Ch14') && strcmp(rpt_cycle, '800cyc')
            fprintf('Skipping C-rate: Ch14 RPT800cyc\n');
            continue;
        end
        
        filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
        filepath = fullfile(ExperimentalDataPath, filename);
        if ~isfile(filepath)
            fprintf('Warning: File not found - %s (skipping C-rate)\n', filename);
            continue;
        end
        try
            % Use readtable with 'preserve' to safely access Current, Voltage, Capacity columns
            T = readtable(filepath, 'VariableNamingRule', 'preserve');
        catch ME
            fprintf('Warning: Unable to read file - %s (skipping C-rate)\n', filename);
            fprintf('  Error: %s\n', ME.message);
            continue;
        end
        
        % 구조 초기화
        Crate_data.(channel_key).(cycle_key) = struct();
        
        for r = 1:length(crate_labels)
            label = crate_labels{r};
            step_chg = crate_steps_charge(r);
            step_dch = crate_steps_discharge(r);
            
            % Charge
            mask_chg = (T.("Step Index") == step_chg);
            V_chg = T.("Voltage(V)")(mask_chg);
            I_chg = T.("Current(A)")(mask_chg);
            Q_chg = T.("Capacity(Ah)")(mask_chg);
            if any(strcmp('Time', T.Properties.VariableNames))
                t_chg = T.("Time")(mask_chg);
            else
                t_chg = [];
            end
            
            % Discharge
            mask_dch = (T.("Step Index") == step_dch);
            V_dch = T.("Voltage(V)")(mask_dch);
            I_dch = T.("Current(A)")(mask_dch);
            Q_dch = T.("Capacity(Ah)")(mask_dch);
            if any(strcmp('Time', T.Properties.VariableNames))
                t_dch = T.("Time")(mask_dch);
            else
                t_dch = [];
            end
            
            % Save Raw Data (No interpolation)
            charge_struct = struct('V', V_chg, 'I', I_chg, 'Q', Q_chg, 't', t_chg);
            discharge_struct = struct('V', V_dch, 'I', I_dch, 'Q', Q_dch, 't', t_dch);
            
            Crate_data.(channel_key).(cycle_key).(label).charge = charge_struct;
            Crate_data.(channel_key).(cycle_key).(label).discharge = discharge_struct;
        end
    end
end

% Removed 8-channel averaging (total8) as it relied on common grid interpolation which is now removed.

% Save C-rate data
matSaveFile_Crate = fullfile(saveDir_Crate, 'Crate_integrated.mat');
save(matSaveFile_Crate, 'Crate_data');
fprintf('Saved C-rate data MAT: %s\n', matSaveFile_Crate);

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
