%% FieldQmax_dQdV.m - dQ/dV calculation for charge and discharge
% Calculates dQ/dV for charge and discharge segments
% Dates: 2021-06-03, 2023-10-16, 2024-09-09, 2025-07-11

clear all; close all; clc;

%% Parameters
C_cell_Ah = 64;                     % Cell capacity (Ah)
thr_A = C_cell_Ah * 0.05;          % Idle threshold (A)
Np = 2;                             % parallel cells (2P)
dt = 1;                             % s (고정 가정)

% Minimum duration for each date (same as Qmax scripts)
% Old data (2021): min_charge_sec = 600, min_discharge_sec = 300
% New data (2023): min_charge_sec = 300, min_discharge_sec = 300
% New data (2024): min_charge_sec = 300 (assumed), min_discharge_sec = 300 (assumed)
% New data (2025): min_charge_sec = 480, min_discharge_sec = 300
min_charge_secs = [600, 300, 300, 480];      % 각 날짜별 최소 충전 구간 길이 (초)
min_discharge_secs = [300, 150, 300, 300];   % 각 날짜별 최소 방전 구간 길이 (초)

% Common Voltage Window (사용자 정의 - 논문 분석용 핵심 구간)
% 자동 계산 대신 사용자가 직접 설정
charge_V_common_min = 3.75;         % Charge 공통 전압 구간 최소값 [V]
charge_V_common_max = 3.95;        % Charge 공통 전압 구간 최대값 [V]
discharge_V_common_min = 3.75;     % Discharge 공통 전압 구간 최소값 [V]
discharge_V_common_max = 3.90;    % Discharge 공통 전압 구간 최대값 [V]

% Voltage range for capacity extraction (Ah 계산용 전압 구간)
capacity_V_min = 3.80;             % 용량 추출 전압 구간 최소값 [V]
capacity_V_max = 3.85;             % 용량 추출 전압 구간 최대값 [V]

%% Paths (필드 데이터 로드 위치)
% 모든 필드 데이터: D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat
% Old data format (2021)
dataFile_2021 = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old\2021\202106\Raw_20210603.mat';

% New data format (2023, 2024, 2025)
dataFile_2023 = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New\2023\202310\Raw_20231016.mat';
dataFile_2024 = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New\2024\202409\Raw_20240909.mat';
dataFile_2025 = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New\2025\202507\Raw_20250711.mat';

% Data structure
dates = {'2021-06-03', '2023-10-16', '2024-09-09', '2025-07-11'};
dataFiles = {dataFile_2021, dataFile_2023, dataFile_2024, dataFile_2025};
dataTypes = {'old', 'new', 'new', 'new'};  % 데이터 형식
base_dates = {datetime(2021,6,3), datetime(2023,10,16), datetime(2024,9,9), datetime(2025,7,11)};

% Storage for all results
All_Charge_dQdV = struct();
All_Discharge_dQdV = struct();

%% Process each date
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    dataFile = dataFiles{date_idx};
    dataType = dataTypes{date_idx};
    base_date = base_dates{date_idx};
    min_charge_sec = min_charge_secs(date_idx);
    min_discharge_sec = min_discharge_secs(date_idx);
    
    fprintf('\n=== Processing %s (%s format) ===\n', date_str, dataType);
    fprintf('Min charge duration: %d sec, Min discharge duration: %d sec\n', min_charge_sec, min_discharge_sec);
    
    %% Load data
    S = load(dataFile);
    
    if strcmp(dataType, 'old')
        % Old data format
        if ~isfield(S, 'Raw')
            error('Raw not found in %s', dataFile);
        end
        Raw_Rack = S.Raw;
        if ~isfield(Raw_Rack, 'Rack01')
            error('Raw_Rack.Rack01 not found in %s', dataFile);
        end
        D = Raw_Rack.Rack01;
        
        % Time vector
        if isfield(D, 'Time')
            t = datetime(D.Time);
        elseif isfield(D, 'Date_Time')
            if isduration(D.Date_Time)
                t = base_date + D.Date_Time;
            else
                t = datetime(D.Date_Time);
            end
        else
            error('No Time/Date_Time field present.');
        end
        
        % Signals
        I_rack = D.DCCurrent_A(:);
        Vcell_avg = D.AverageCV_V(:);
        SOC_raw = D.SOCPct(:);  % Raw SOC from BMS (%)
        
    else
        % New data format
        D = S.Raw;
        
        % Time
        if isduration(D.Date_Time)
            t = base_date + D.Date_Time;
        else
            t = datetime(D.Date_Time);
        end
        
        % Signals
        I_rack = D.DCCurrent(:);
        Vcell_avg = D.CVavg(:);
        SOC_raw = D.SOC_BMS(:);  % Raw SOC from BMS (%)
    end
    
    t0 = t(1);
    tsec = seconds(t - t0);
    
    % Cell-level signals
    I_cell = I_rack / Np;  % A per cell
    
    % Smoothing: Moving average with 10 second window (disabled for segment detection)
    % window_size_samples = 10;  % 10 seconds (dt = 1s)
    % I_cell_smooth = movmean(I_cell, window_size_samples);
    % Vcell_avg_smooth = movmean(Vcell_avg, window_size_samples);
    
    % Use original data for segment detection (no smoothing)
    % I_cell = I_cell_smooth;
    % Vcell_avg = Vcell_avg_smooth;
    
    % Threshold for cell level
    thr_cell = thr_A / Np;
    
    %% Find charge and discharge segments (same method as Qmax scripts)
    isIdle = abs(I_cell) < thr_cell;
    isChg = I_cell > thr_cell;
    isDischg = I_cell < -thr_cell;
    
    % Helper function to find contiguous segments
    find_segments = @(mask) local_find_segments(mask);
    
    idleSegs = find_segments(isIdle);
    chgSegs = find_segments(isChg);
    dischgSegs = find_segments(isDischg);
    
    % Filter by minimum duration (same as Qmax scripts)
    dur = @(seg) seg(:,2) - seg(:,1) + 1;
    min_samples_chg = ceil(min_charge_sec / dt);
    min_samples_dischg = ceil(min_discharge_sec / dt);
    
    chgSegs = chgSegs(dur(chgSegs) >= min_samples_chg, :);
    dischgSegs = dischgSegs(dur(dischgSegs) >= min_samples_dischg, :);
    
    fprintf('Found %d charge segments (>= %d sec)\n', size(chgSegs,1), min_charge_sec);
    fprintf('Found %d discharge segments (>= %d sec)\n', size(dischgSegs,1), min_discharge_sec);
    
    %% Calculate dQ/dV for charge segments
    Charge_dQdV_data = [];
    % 2024년도는 Chg02만 사용 (시간: 12:55-13:13)
    if strcmp(date_str, '2024-09-09')
        % Chg02 시간 범위로 매칭
        target_chg_start = datetime(2024, 9, 9, 12, 55, 0);
        target_chg_end = datetime(2024, 9, 9, 13, 13, 0);
        chg_seg_indices = [];
        for k = 1:size(chgSegs,1)
            seg_start_time = t(chgSegs(k,1));
            seg_end_time = t(chgSegs(k,2));
            % 시간 범위로 매칭 (약 1분 허용)
            if abs(seg_start_time - target_chg_start) < minutes(5) && abs(seg_end_time - target_chg_end) < minutes(5)
                chg_seg_indices = k;
                fprintf('2024-09-09: Matched charge segment %d (Chg02) at time %s to %s\n', k, datestr(seg_start_time), datestr(seg_end_time));
                break;
            end
        end
        if isempty(chg_seg_indices)
            fprintf('Warning: 2024-09-09 Chg02 segment not found\n');
        end
    else
        chg_seg_indices = 1:size(chgSegs,1);  % 모든 세그먼트
    end
    
    for k = chg_seg_indices
        chg_start = chgSegs(k,1);
        chg_end = chgSegs(k,2);
        
        % Extract segment data
        I_seg = I_cell(chg_start:chg_end);
        V_seg = Vcell_avg(chg_start:chg_end);
        t_seg = tsec(chg_start:chg_end);
        
        % Calculate cumulative charge Q (Ah)
        % Q = ∫ I dt (시간에 대한 전류 적분)
        Q_seg = cumtrapz(t_seg, I_seg) / 3600;  % Ah
        
        % Sort by voltage (ascending for charge)
        % 전압 순서로 정렬 (dQ/dV 계산을 위해)
        [V_sorted, sort_idx] = sort(V_seg);
        Q_sorted = Q_seg(sort_idx);
        
        % Remove duplicate voltage values (보간을 위해)
        [V_unique, unique_idx, ~] = unique(V_sorted, 'stable');
        Q_unique = Q_sorted(unique_idx);
        
        % Grid 보간법: 고정된 전압 간격(0.001V)으로 Grid 생성
        if length(V_unique) > 1
            V_min = min(V_unique);
            V_max = max(V_unique);
            dV_grid = 0.005;   % 10mV로 변경 (필드 데이터 노이즈 감소)
            V_grid = V_min:dV_grid:V_max;  % Grid 전압 벡터
            
            % 보간: 불규칙한 원본 데이터를 고른 Grid 위로 보간
            Q_grid = interp1(V_unique, Q_unique, V_grid, 'linear');
            
            % [핵심] Q 데이터를 부드럽게 다듬기 (전압 도메인 스무딩)
            Q_grid = smoothdata(Q_grid, 'gaussian', 15);
            
            % dQ/dV 계산: diff(Q_grid) / dV_grid (분모가 상수!)
            if length(Q_grid) > 1
                dQ = diff(Q_grid);  % ΔQ
                dV = dV_grid;       % 고정된 간격 0.001V
                
                % 중간 전압 값 사용 (V_mid = V_grid의 중간점)
                V_mid = (V_grid(1:end-1) + V_grid(2:end)) / 2;
                
                % dQ/dV 계산 (분모가 상수이므로 안정적)
                dQdV = dQ ./ dV;  % Ah/V
                
                % Store data (use sequential index for storage)
                data_idx = length(Charge_dQdV_data) + 1;
                Charge_dQdV_data(data_idx).V_dQdV = V_mid;
                Charge_dQdV_data(data_idx).dQdV = dQdV;
                Charge_dQdV_data(data_idx).V_raw = V_seg;
                Charge_dQdV_data(data_idx).I_raw = I_seg;
                Charge_dQdV_data(data_idx).Q_raw = Q_seg;  % Store Q for V-Q plot
                Charge_dQdV_data(data_idx).t_raw = t(chg_start:chg_end);
                Charge_dQdV_data(data_idx).SOC_raw = SOC_raw(chg_start:chg_end);
                Charge_dQdV_data(data_idx).date = date_str;
                Charge_dQdV_data(data_idx).V_min = V_min;  % Store V range for common window calculation
                Charge_dQdV_data(data_idx).V_max = V_max;
            end
        end
    end
    
    %% Calculate dQ/dV for discharge segments
    Discharge_dQdV_data = [];
    % 2024년도는 Dchg03만 사용 (시간: 14:14-16:45)
    if strcmp(date_str, '2024-09-09')
        % Dchg03 시간 범위로 매칭
        target_dischg_start = datetime(2024, 9, 9, 14, 14, 0);
        target_dischg_end = datetime(2024, 9, 9, 16, 45, 0);
        dischg_seg_indices = [];
        
        % 디버깅: 모든 방전 세그먼트 시간 출력
        fprintf('2024-09-09: Available discharge segments:\n');
        for k = 1:size(dischgSegs,1)
            seg_start_time = t(dischgSegs(k,1));
            seg_end_time = t(dischgSegs(k,2));
            fprintf('  Segment %d: %s to %s\n', k, datestr(seg_start_time), datestr(seg_end_time));
        end
        fprintf('  Target (Dchg03): %s to %s\n', datestr(target_dischg_start), datestr(target_dischg_end));
        
        for k = 1:size(dischgSegs,1)
            seg_start_time = t(dischgSegs(k,1));
            seg_end_time = t(dischgSegs(k,2));
            % 시작 시간만 매칭 (표의 End 시간은 휴지 기간 포함일 수 있음)
            time_diff_start = abs(seg_start_time - target_dischg_start);
            if time_diff_start < minutes(1)  % 시작 시간 1분 이내
                dischg_seg_indices = k;
                fprintf('2024-09-09: Matched discharge segment %d (Dchg03) at time %s to %s\n', k, datestr(seg_start_time), datestr(seg_end_time));
                break;
            end
        end
        if isempty(dischg_seg_indices)
            fprintf('Warning: 2024-09-09 Dchg03 segment not found\n');
        end
    else
        dischg_seg_indices = 1:size(dischgSegs,1);  % 모든 세그먼트
    end
    
    for k = dischg_seg_indices
        dischg_start = dischgSegs(k,1);
        dischg_end = dischgSegs(k,2);
        
        % 2024년도: 14:25:49 이후는 rest 처리 (dQ/dV 계산에서 제외)
        if strcmp(date_str, '2024-09-09')
            target_rest_start = datetime(2024, 9, 9, 14, 25, 49);
            % dischg_end 시간이 14:25:49 이후인 경우, 14:25:49 이전으로 제한
            if t(dischg_end) > target_rest_start
                % 14:25:49 이전의 마지막 인덱스 찾기
                idx_limit = find(t <= target_rest_start, 1, 'last');
                if ~isempty(idx_limit) && idx_limit >= dischg_start
                    dischg_end = idx_limit;
                    fprintf('2024-09-09: Limited discharge end to 14:25:49 (idx %d, time %s)\n', ...
                        dischg_end, datestr(t(dischg_end), 'HH:MM:SS'));
                end
            end
        end
        
        % Extract segment data
        I_seg = I_cell(dischg_start:dischg_end);
        V_seg = Vcell_avg(dischg_start:dischg_end);
        t_seg = tsec(dischg_start:dischg_end);
        
        % Calculate cumulative discharge Q (Ah)
        % Q = ∫ |I| dt (시간에 대한 전류 절댓값 적분, 양수로 계산)
        Q_seg = cumtrapz(t_seg, abs(I_seg)) / 3600;  % Ah (positive)
        
        % Sort by voltage (descending for discharge)
        % 전압 순서로 내림차순 정렬 (방전은 전압이 감소하므로)
        [V_sorted, sort_idx] = sort(V_seg, 'descend');
        Q_sorted = Q_seg(sort_idx);
        
        % Remove duplicate voltage values (보간을 위해)
        [V_unique, unique_idx, ~] = unique(V_sorted, 'stable');
        Q_unique = Q_sorted(unique_idx);
        
        % Grid 보간법: 고정된 전압 간격(0.001V)으로 Grid 생성
        if length(V_unique) > 1
            V_max = max(V_unique);
            V_min = min(V_unique);
            dV_grid = 0.005;   % 10mV로 변경 (필드 데이터 노이즈 감소)
            V_grid = V_max:-dV_grid:V_min;  % Grid 전압 벡터 (높은 전압부터 낮은 전압으로)
            
            % 보간: 불규칙한 원본 데이터를 고른 Grid 위로 보간
            Q_grid = interp1(V_unique, Q_unique, V_grid, 'linear');
            
            % [핵심] Q 데이터를 부드럽게 다듬기 (전압 도메인 스무딩)
            Q_grid = smoothdata(Q_grid, 'gaussian', 15);
            
            % dQ/dV 계산: diff(Q_grid) / dV_grid (분모가 상수!)
            if length(Q_grid) > 1
                dQ = diff(Q_grid);  % ΔQ
                dV = dV_grid;       % 고정된 간격 0.001V
                
                % 중간 전압 값 사용 (V_mid = V_grid의 중간점)
                V_mid = (V_grid(1:end-1) + V_grid(2:end)) / 2;
                
                % dQ/dV 계산 (분모가 상수이므로 안정적)
                dQdV = dQ ./ dV;  % Ah/V (positive for discharge)
                
                % Store data (use sequential index for storage)
                data_idx = length(Discharge_dQdV_data) + 1;
                Discharge_dQdV_data(data_idx).V_dQdV = V_mid;
                Discharge_dQdV_data(data_idx).dQdV = dQdV;
                Discharge_dQdV_data(data_idx).V_raw = V_seg;
                Discharge_dQdV_data(data_idx).I_raw = I_seg;
                Discharge_dQdV_data(data_idx).Q_raw = Q_seg;  % Store Q for V-Q plot
                Discharge_dQdV_data(data_idx).t_raw = t(dischg_start:dischg_end);
                Discharge_dQdV_data(data_idx).SOC_raw = SOC_raw(dischg_start:dischg_end);
                Discharge_dQdV_data(data_idx).date = date_str;
                Discharge_dQdV_data(data_idx).V_min = V_min;  % Store V range for common window calculation
                Discharge_dQdV_data(data_idx).V_max = V_max;
            end
        end
    end
    
    % Store in global structure
    All_Charge_dQdV.(sprintf('date_%d', date_idx)) = Charge_dQdV_data;
    All_Discharge_dQdV.(sprintf('date_%d', date_idx)) = Discharge_dQdV_data;
    
    %% Plot for individual date (Charge and Discharge with I, V, Time-SOC subplots)
    fig_date = figure('Name', sprintf('dQ/dV Analysis: %s', date_str), 'NumberTitle', 'off');
    tl_date = tiledlayout(fig_date, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Charge dQ/dV
    ax1 = nexttile(tl_date, 1); hold(ax1, 'on'); grid(ax1, 'on');
    colors = lines(length(Charge_dQdV_data));
    for k = 1:length(Charge_dQdV_data)
        if ~isempty(Charge_dQdV_data(k).V_dQdV)
            plot(ax1, Charge_dQdV_data(k).V_dQdV, Charge_dQdV_data(k).dQdV, '-o', ...
                'Color', colors(k,:), 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', sprintf('Segment %d', k));
        end
    end
    title(ax1, sprintf('Charge dQ/dV: %s', date_str));
    xlabel(ax1, 'Voltage [V]');
    ylabel(ax1, 'dQ/dV [Ah/V]');
    legend(ax1, 'Location', 'best');
    
    % Charge Current
    ax2 = nexttile(tl_date, 2); hold(ax2, 'on'); grid(ax2, 'on');
    for k = 1:length(Charge_dQdV_data)
        if ~isempty(Charge_dQdV_data(k).I_raw)
            % For time-series plots, show markers at intervals to avoid clutter
            t_data = Charge_dQdV_data(k).t_raw;
            I_data = Charge_dQdV_data(k).I_raw;
            marker_interval = max(1, floor(length(t_data) / 50));  % Show ~50 markers max
            idx_markers = 1:marker_interval:length(t_data);
            plot(ax2, t_data, I_data, '-', ...
                'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Segment %d', k));
            plot(ax2, t_data(idx_markers), I_data(idx_markers), 'o', ...
                'Color', colors(k,:), 'MarkerSize', 4, 'MarkerFaceColor', colors(k,:), 'HandleVisibility', 'off');
        end
    end
    title(ax2, 'Charge Current');
    xlabel(ax2, 'Time');
    ylabel(ax2, 'Current [A]');
    legend(ax2, 'Location', 'best');
    
    % Charge V-Q Curve
    ax3 = nexttile(tl_date, 3); hold(ax3, 'on'); grid(ax3, 'on');
    for k = 1:length(Charge_dQdV_data)
        if ~isempty(Charge_dQdV_data(k).V_raw) && ~isempty(Charge_dQdV_data(k).Q_raw)
            Q_data = Charge_dQdV_data(k).Q_raw;
            V_data = Charge_dQdV_data(k).V_raw;
            marker_interval = max(1, floor(length(Q_data) / 50));  % Show ~50 markers max
            idx_markers = 1:marker_interval:length(Q_data);
            plot(ax3, Q_data, V_data, '-', ...
                'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Segment %d', k));
            plot(ax3, Q_data(idx_markers), V_data(idx_markers), 'o', ...
                'Color', colors(k,:), 'MarkerSize', 4, 'MarkerFaceColor', colors(k,:), 'HandleVisibility', 'off');
        end
    end
    title(ax3, 'Charge V-Q Curve');
    xlabel(ax3, 'Capacity [Ah]');
    ylabel(ax3, 'Voltage [V]');
    legend(ax3, 'Location', 'best');
    
    % Charge Time-SOC
    ax4 = nexttile(tl_date, 4); hold(ax4, 'on'); grid(ax4, 'on');
    for k = 1:length(Charge_dQdV_data)
        if ~isempty(Charge_dQdV_data(k).t_raw) && ~isempty(Charge_dQdV_data(k).SOC_raw)
            t_data = Charge_dQdV_data(k).t_raw;
            SOC_data = Charge_dQdV_data(k).SOC_raw;
            marker_interval = max(1, floor(length(t_data) / 50));  % Show ~50 markers max
            idx_markers = 1:marker_interval:length(t_data);
            plot(ax4, t_data, SOC_data, '-', ...
                'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Segment %d', k));
            plot(ax4, t_data(idx_markers), SOC_data(idx_markers), 'o', ...
                'Color', colors(k,:), 'MarkerSize', 4, 'MarkerFaceColor', colors(k,:), 'HandleVisibility', 'off');
        end
    end
    title(ax4, 'Charge Time-SOC');
    xlabel(ax4, 'Time');
    ylabel(ax4, 'SOC [%]');
    legend(ax4, 'Location', 'best');
    
    % Discharge dQ/dV
    ax5 = nexttile(tl_date, 5); hold(ax5, 'on'); grid(ax5, 'on');
    colors_dischg = lines(length(Discharge_dQdV_data));
    for k = 1:length(Discharge_dQdV_data)
        if ~isempty(Discharge_dQdV_data(k).V_dQdV)
            plot(ax5, Discharge_dQdV_data(k).V_dQdV, Discharge_dQdV_data(k).dQdV, '-o', ...
                'Color', colors_dischg(k,:), 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', sprintf('Segment %d', k));
        end
    end
    title(ax5, sprintf('Discharge dQ/dV: %s', date_str));
    xlabel(ax5, 'Voltage [V]');
    ylabel(ax5, 'dQ/dV [Ah/V]');
    legend(ax5, 'Location', 'best');
    
    % Discharge Current
    ax6 = nexttile(tl_date, 6); hold(ax6, 'on'); grid(ax6, 'on');
    for k = 1:length(Discharge_dQdV_data)
        if ~isempty(Discharge_dQdV_data(k).I_raw)
            t_data = Discharge_dQdV_data(k).t_raw;
            I_data = Discharge_dQdV_data(k).I_raw;
            marker_interval = max(1, floor(length(t_data) / 50));  % Show ~50 markers max
            idx_markers = 1:marker_interval:length(t_data);
            plot(ax6, t_data, I_data, '-', ...
                'Color', colors_dischg(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Segment %d', k));
            plot(ax6, t_data(idx_markers), I_data(idx_markers), 'o', ...
                'Color', colors_dischg(k,:), 'MarkerSize', 4, 'MarkerFaceColor', colors_dischg(k,:), 'HandleVisibility', 'off');
        end
    end
    title(ax6, 'Discharge Current');
    xlabel(ax6, 'Time');
    ylabel(ax6, 'Current [A]');
    legend(ax6, 'Location', 'best');
    
    % Discharge V-Q Curve
    ax7 = nexttile(tl_date, 7); hold(ax7, 'on'); grid(ax7, 'on');
    for k = 1:length(Discharge_dQdV_data)
        if ~isempty(Discharge_dQdV_data(k).V_raw) && ~isempty(Discharge_dQdV_data(k).Q_raw)
            Q_data = Discharge_dQdV_data(k).Q_raw;
            V_data = Discharge_dQdV_data(k).V_raw;
            marker_interval = max(1, floor(length(Q_data) / 50));  % Show ~50 markers max
            idx_markers = 1:marker_interval:length(Q_data);
            plot(ax7, Q_data, V_data, '-', ...
                'Color', colors_dischg(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Segment %d', k));
            plot(ax7, Q_data(idx_markers), V_data(idx_markers), 'o', ...
                'Color', colors_dischg(k,:), 'MarkerSize', 4, 'MarkerFaceColor', colors_dischg(k,:), 'HandleVisibility', 'off');
        end
    end
    title(ax7, 'Discharge V-Q Curve');
    xlabel(ax7, 'Capacity [Ah]');
    ylabel(ax7, 'Voltage [V]');
    legend(ax7, 'Location', 'best');
    
    % Discharge Time-SOC
    ax8 = nexttile(tl_date, 8); hold(ax8, 'on'); grid(ax8, 'on');
    for k = 1:length(Discharge_dQdV_data)
        if ~isempty(Discharge_dQdV_data(k).t_raw) && ~isempty(Discharge_dQdV_data(k).SOC_raw)
            t_data = Discharge_dQdV_data(k).t_raw;
            SOC_data = Discharge_dQdV_data(k).SOC_raw;
            marker_interval = max(1, floor(length(t_data) / 50));  % Show ~50 markers max
            idx_markers = 1:marker_interval:length(t_data);
            plot(ax8, t_data, SOC_data, '-', ...
                'Color', colors_dischg(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Segment %d', k));
            plot(ax8, t_data(idx_markers), SOC_data(idx_markers), 'o', ...
                'Color', colors_dischg(k,:), 'MarkerSize', 4, 'MarkerFaceColor', colors_dischg(k,:), 'HandleVisibility', 'off');
        end
    end
    title(ax8, 'Discharge Time-SOC');
    xlabel(ax8, 'Time');
    ylabel(ax8, 'SOC [%]');
    legend(ax8, 'Location', 'best');
    
    % Save figure
    saveas(fig_date, sprintf('FieldQmax_dQdV_%s.fig', date_str));
    saveas(fig_date, sprintf('FieldQmax_dQdV_%s.png', date_str));
end

%% Step 2: Use User-Defined Common Voltage Window
fprintf('\n=== Step 2: Using User-Defined Common Voltage Window ===\n');
fprintf('Charge Common Voltage Window: %.3f V to %.3f V\n', charge_V_common_min, charge_V_common_max);
fprintf('Discharge Common Voltage Window: %.3f V to %.3f V\n', discharge_V_common_min, discharge_V_common_max);
fprintf('Capacity Extraction Range: %.3f V to %.3f V\n', capacity_V_min, capacity_V_max);

%% Step 3: Recalculate dQ/dV within Common Voltage Window
fprintf('\n=== Step 3: Recalculating dQ/dV within Common Voltage Window ===\n');

dV_grid = 0.005;  % Fixed grid spacing (5mV)

% Recalculate charge dQ/dV within common window
if charge_V_common_min < charge_V_common_max
    % Create fixed common grid (모든 연도에 동일)
    V_grid_common_chg = charge_V_common_min:dV_grid:charge_V_common_max;
    num_grid_points = length(V_grid_common_chg);
    fprintf('Charge common grid: %d points (%.3f V to %.3f V, step %.3f V)\n', ...
        num_grid_points, charge_V_common_min, charge_V_common_max, dV_grid);
    
    for date_idx = 1:length(dates)
        date_str = dates{date_idx};
        Charge_data = All_Charge_dQdV.(sprintf('date_%d', date_idx));
        
        for k = 1:length(Charge_data)
            if ~isempty(Charge_data(k).V_raw) && ~isempty(Charge_data(k).Q_raw)
                V_seg = Charge_data(k).V_raw;
                Q_seg = Charge_data(k).Q_raw;
                
                % Sort by voltage (ascending for charge)
                [V_sorted, sort_idx] = sort(V_seg);
                Q_sorted = Q_seg(sort_idx);
                
                % Remove duplicate voltage values
                [V_unique, unique_idx, ~] = unique(V_sorted, 'stable');
                Q_unique = Q_sorted(unique_idx);
                
                % Filter to common window (원본 데이터 범위 내에서만)
                valid_idx = V_unique >= charge_V_common_min & V_unique <= charge_V_common_max;
                V_unique_common = V_unique(valid_idx);
                Q_unique_common = Q_unique(valid_idx);
                
                if length(V_unique_common) > 1
                    % Interpolate to common grid (NaN 허용 - 데이터 없는 곳은 NaN)
                    Q_grid_common = interp1(V_unique_common, Q_unique_common, V_grid_common_chg, 'linear', NaN);
                    
                    % NaN을 제거하지 않고 그대로 사용 (모든 연도가 동일한 그리드 포인트 수 유지)
                    % Smoothing (NaN은 자동으로 처리됨)
                    Q_grid_common = smoothdata(Q_grid_common, 'gaussian', 15, 'omitnan');
                    
                    % Calculate dQ/dV (NaN 포함)
                    dQ = diff(Q_grid_common);
                    V_mid_common = (V_grid_common_chg(1:end-1) + V_grid_common_chg(2:end)) / 2;
                    dQdV_common = dQ ./ dV_grid;
                    
                    % Store common window results (NaN 포함, 모든 연도가 동일한 길이)
                    Charge_data(k).V_dQdV_common = V_mid_common;
                    Charge_data(k).dQdV_common = dQdV_common;
                    Charge_data(k).V_grid_common = V_grid_common_chg;  % 전체 그리드 저장
                    Charge_data(k).Q_grid_common = Q_grid_common;      % 전체 Q 그리드 저장 (NaN 포함)
                    Charge_data(k).V_common_min = charge_V_common_min;
                    Charge_data(k).V_common_max = charge_V_common_max;
                    
                    % Capacity extraction: 특정 전압 구간 내 총 충전 용량 계산
                    if capacity_V_min >= charge_V_common_min && capacity_V_max <= charge_V_common_max
                        capacity_idx = V_grid_common_chg >= capacity_V_min & V_grid_common_chg <= capacity_V_max;
                        Q_in_range = Q_grid_common(capacity_idx);
                        Q_in_range_valid = Q_in_range(~isnan(Q_in_range));
                        if ~isempty(Q_in_range_valid)
                            capacity_Ah = max(Q_in_range_valid) - min(Q_in_range_valid);
                            Charge_data(k).capacity_Ah = capacity_Ah;
                            Charge_data(k).capacity_V_range = [capacity_V_min, capacity_V_max];
                        else
                            Charge_data(k).capacity_Ah = NaN;
                            Charge_data(k).capacity_V_range = [capacity_V_min, capacity_V_max];
                        end
                    end
                end
            end
        end
        All_Charge_dQdV.(sprintf('date_%d', date_idx)) = Charge_data;
    end
end

% Recalculate discharge dQ/dV within common window
if discharge_V_common_min < discharge_V_common_max
    % Create fixed common grid (모든 연도에 동일, ascending for interpolation)
    V_grid_common_dischg = discharge_V_common_min:dV_grid:discharge_V_common_max;
    num_grid_points = length(V_grid_common_dischg);
    fprintf('Discharge common grid: %d points (%.3f V to %.3f V, step %.3f V)\n', ...
        num_grid_points, discharge_V_common_min, discharge_V_common_max, dV_grid);
    
    for date_idx = 1:length(dates)
        date_str = dates{date_idx};
        Discharge_data = All_Discharge_dQdV.(sprintf('date_%d', date_idx));
        
        for k = 1:length(Discharge_data)
            if ~isempty(Discharge_data(k).V_raw) && ~isempty(Discharge_data(k).Q_raw)
                V_seg = Discharge_data(k).V_raw;
                Q_seg = Discharge_data(k).Q_raw;
                
                % Sort by voltage (descending for discharge)
                [V_sorted, sort_idx] = sort(V_seg, 'descend');
                Q_sorted = Q_seg(sort_idx);
                
                % Remove duplicate voltage values
                [V_unique, unique_idx, ~] = unique(V_sorted, 'stable');
                Q_unique = Q_sorted(unique_idx);
                
                % Filter to common window (원본 데이터 범위 내에서만)
                valid_idx = V_unique >= discharge_V_common_min & V_unique <= discharge_V_common_max;
                V_unique_common = V_unique(valid_idx);
                Q_unique_common = Q_unique(valid_idx);
                
                if length(V_unique_common) > 1
                    % For interpolation, we need ascending order
                    V_unique_common_asc = sort(V_unique_common, 'ascend');
                    Q_unique_common_asc = interp1(V_unique_common, Q_unique_common, V_unique_common_asc, 'linear');
                    
                    % Interpolate to common grid (NaN 허용)
                    Q_grid_common = interp1(V_unique_common_asc, Q_unique_common_asc, V_grid_common_dischg, 'linear', NaN);
                    
                    % NaN을 제거하지 않고 그대로 사용 (모든 연도가 동일한 그리드 포인트 수 유지)
                    % Smoothing (NaN은 자동으로 처리됨)
                    Q_grid_common = smoothdata(Q_grid_common, 'gaussian', 15, 'omitnan');
                    
                    % For discharge, calculate dQ/dV (Q decreases as V decreases)
                    dQ = -diff(Q_grid_common);  % Negative because Q decreases as V decreases
                    V_mid_common = (V_grid_common_dischg(1:end-1) + V_grid_common_dischg(2:end)) / 2;
                    dQdV_common = dQ ./ dV_grid;  % Should be positive
                    
                    % Store common window results (NaN 포함, 모든 연도가 동일한 길이)
                    Discharge_data(k).V_dQdV_common = V_mid_common;
                    Discharge_data(k).dQdV_common = dQdV_common;
                    Discharge_data(k).V_grid_common = V_grid_common_dischg;  % 전체 그리드 저장
                    Discharge_data(k).Q_grid_common = Q_grid_common;         % 전체 Q 그리드 저장 (NaN 포함)
                    Discharge_data(k).V_common_min = discharge_V_common_min;
                    Discharge_data(k).V_common_max = discharge_V_common_max;
                    
                    % Capacity extraction: 특정 전압 구간 내 총 방전 용량 계산
                    if capacity_V_min >= discharge_V_common_min && capacity_V_max <= discharge_V_common_max
                        capacity_idx = V_grid_common_dischg >= capacity_V_min & V_grid_common_dischg <= capacity_V_max;
                        Q_in_range = Q_grid_common(capacity_idx);
                        Q_in_range_valid = Q_in_range(~isnan(Q_in_range));
                        if ~isempty(Q_in_range_valid)
                            capacity_Ah = max(Q_in_range_valid) - min(Q_in_range_valid);
                            Discharge_data(k).capacity_Ah = capacity_Ah;
                            Discharge_data(k).capacity_V_range = [capacity_V_min, capacity_V_max];
                        else
                            Discharge_data(k).capacity_Ah = NaN;
                            Discharge_data(k).capacity_V_range = [capacity_V_min, capacity_V_max];
                        end
                    end
                end
            end
        end
        All_Discharge_dQdV.(sprintf('date_%d', date_idx)) = Discharge_data;
    end
end

%% Plot all years together - Charge dQ/dV with I, V
fig_charge_all = figure('Name', 'Charge dQ/dV: All Years (2021-2025)', 'NumberTitle', 'off');
tl_charge = tiledlayout(fig_charge_all, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Year colors mapping: [2021, 2023, 2024, 2025] -> [Blue, Green, Orange, Purple]
year_colors = [0 0 1; 0 0.8 0; 1 0.5 0; 0.8 0 0.8];  % Blue, Green, Orange, Purple

% Charge dQ/dV (Common Window)
ax1 = nexttile(tl_charge, 1); hold(ax1, 'on'); grid(ax1, 'on');
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Charge_data = All_Charge_dQdV.(sprintf('date_%d', date_idx));
    
    for k = 1:length(Charge_data)
        % Use common window data if available, otherwise use original
        if isfield(Charge_data(k), 'V_dQdV_common') && ~isempty(Charge_data(k).V_dQdV_common)
            plot(ax1, Charge_data(k).V_dQdV_common, Charge_data(k).dQdV_common, '-o', ...
                'Color', year_colors(date_idx,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
                'DisplayName', sprintf('%s Seg%d', date_str, k));
        elseif ~isempty(Charge_data(k).V_dQdV)
            plot(ax1, Charge_data(k).V_dQdV, Charge_data(k).dQdV, '-o', ...
                'Color', year_colors(date_idx,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
                'DisplayName', sprintf('%s Seg%d', date_str, k));
        end
    end
end
if ~isempty(charge_V_common_min) && ~isempty(charge_V_common_max)
    title(ax1, sprintf('Charge dQ/dV: All Years (Common: %.2f-%.2f V)', charge_V_common_min, charge_V_common_max));
else
    title(ax1, 'Charge dQ/dV: All Years');
end
xlabel(ax1, 'Voltage [V]');
ylabel(ax1, 'dQ/dV [Ah/V]');
legend(ax1, 'Location', 'best', 'NumColumns', 2);

% Charge Current (all years)
ax2 = nexttile(tl_charge, 2); hold(ax2, 'on'); grid(ax2, 'on');
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Charge_data = All_Charge_dQdV.(sprintf('date_%d', date_idx));
    
    for k = 1:length(Charge_data)
        if ~isempty(Charge_data(k).I_raw)
            % 시간을 0부터 시작하도록 normalize (초 단위)
            t_normalized = seconds(Charge_data(k).t_raw - Charge_data(k).t_raw(1));
            I_data = Charge_data(k).I_raw;
            marker_interval = max(1, floor(length(t_normalized) / 50));
            idx_markers = 1:marker_interval:length(t_normalized);
            plot(ax2, t_normalized, I_data, '-', ...
                'Color', year_colors(date_idx,:), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s Seg%d', date_str, k));
            plot(ax2, t_normalized(idx_markers), I_data(idx_markers), 'o', ...
                'Color', year_colors(date_idx,:), 'MarkerSize', 4, 'MarkerFaceColor', year_colors(date_idx,:), 'HandleVisibility', 'off');
        end
    end
end
title(ax2, 'Charge Current: All Years');
xlabel(ax2, 'Time [s]');
ylabel(ax2, 'Current [A]');
legend(ax2, 'Location', 'best', 'NumColumns', 2);

% Charge V-Q Curve (all years)
ax3 = nexttile(tl_charge, 3); hold(ax3, 'on'); grid(ax3, 'on');
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Charge_data = All_Charge_dQdV.(sprintf('date_%d', date_idx));
    
    for k = 1:length(Charge_data)
        if ~isempty(Charge_data(k).V_raw) && ~isempty(Charge_data(k).Q_raw)
            Q_data = Charge_data(k).Q_raw;
            V_data = Charge_data(k).V_raw;
            marker_interval = max(1, floor(length(Q_data) / 50));
            idx_markers = 1:marker_interval:length(Q_data);
            plot(ax3, Q_data, V_data, '-', ...
                'Color', year_colors(date_idx,:), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s Seg%d', date_str, k));
            plot(ax3, Q_data(idx_markers), V_data(idx_markers), 'o', ...
                'Color', year_colors(date_idx,:), 'MarkerSize', 4, 'MarkerFaceColor', year_colors(date_idx,:), 'HandleVisibility', 'off');
        end
    end
end
title(ax3, 'Charge V-Q Curve: All Years');
xlabel(ax3, 'Capacity [Ah]');
ylabel(ax3, 'Voltage [V]');
legend(ax3, 'Location', 'best', 'NumColumns', 2);

saveas(fig_charge_all, 'FieldQmax_dQdV_Charge_AllYears.fig');
saveas(fig_charge_all, 'FieldQmax_dQdV_Charge_AllYears.png');

%% Plot all years together - Discharge dQ/dV with I, V
fig_discharge_all = figure('Name', 'Discharge dQ/dV: All Years (2021-2025)', 'NumberTitle', 'off');
tl_discharge = tiledlayout(fig_discharge_all, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Discharge dQ/dV (Common Window)
ax1 = nexttile(tl_discharge, 1); hold(ax1, 'on'); grid(ax1, 'on');
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Discharge_data = All_Discharge_dQdV.(sprintf('date_%d', date_idx));
    
    for k = 1:length(Discharge_data)
        % Use common window data if available, otherwise use original
        if isfield(Discharge_data(k), 'V_dQdV_common') && ~isempty(Discharge_data(k).V_dQdV_common)
            plot(ax1, Discharge_data(k).V_dQdV_common, Discharge_data(k).dQdV_common, '-o', ...
                'Color', year_colors(date_idx,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
                'DisplayName', sprintf('%s Seg%d', date_str, k));
        elseif ~isempty(Discharge_data(k).V_dQdV)
            plot(ax1, Discharge_data(k).V_dQdV, Discharge_data(k).dQdV, '-o', ...
                'Color', year_colors(date_idx,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
                'DisplayName', sprintf('%s Seg%d', date_str, k));
        end
    end
end
if ~isempty(discharge_V_common_min) && ~isempty(discharge_V_common_max)
    title(ax1, sprintf('Discharge dQ/dV: All Years (Common: %.2f-%.2f V)', discharge_V_common_min, discharge_V_common_max));
else
    title(ax1, 'Discharge dQ/dV: All Years');
end
xlabel(ax1, 'Voltage [V]');
ylabel(ax1, 'dQ/dV [Ah/V]');
legend(ax1, 'Location', 'best', 'NumColumns', 2);

% Discharge Current (all years)
ax2 = nexttile(tl_discharge, 2); hold(ax2, 'on'); grid(ax2, 'on');
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Discharge_data = All_Discharge_dQdV.(sprintf('date_%d', date_idx));
    
    for k = 1:length(Discharge_data)
        if ~isempty(Discharge_data(k).I_raw)
            % 시간을 0부터 시작하도록 normalize (초 단위)
            t_normalized = seconds(Discharge_data(k).t_raw - Discharge_data(k).t_raw(1));
            I_data = Discharge_data(k).I_raw;
            marker_interval = max(1, floor(length(t_normalized) / 50));
            idx_markers = 1:marker_interval:length(t_normalized);
            plot(ax2, t_normalized, I_data, '-', ...
                'Color', year_colors(date_idx,:), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s Seg%d', date_str, k));
            plot(ax2, t_normalized(idx_markers), I_data(idx_markers), 'o', ...
                'Color', year_colors(date_idx,:), 'MarkerSize', 4, 'MarkerFaceColor', year_colors(date_idx,:), 'HandleVisibility', 'off');
        end
    end
end
title(ax2, 'Discharge Current: All Years');
xlabel(ax2, 'Time [s]');
ylabel(ax2, 'Current [A]');
legend(ax2, 'Location', 'best', 'NumColumns', 2);

% Discharge V-Q Curve (all years)
ax3 = nexttile(tl_discharge, 3); hold(ax3, 'on'); grid(ax3, 'on');
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Discharge_data = All_Discharge_dQdV.(sprintf('date_%d', date_idx));
    
    for k = 1:length(Discharge_data)
        if ~isempty(Discharge_data(k).V_raw) && ~isempty(Discharge_data(k).Q_raw)
            Q_data = Discharge_data(k).Q_raw;
            V_data = Discharge_data(k).V_raw;
            marker_interval = max(1, floor(length(Q_data) / 50));
            idx_markers = 1:marker_interval:length(Q_data);
            plot(ax3, Q_data, V_data, '-', ...
                'Color', year_colors(date_idx,:), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s Seg%d', date_str, k));
            plot(ax3, Q_data(idx_markers), V_data(idx_markers), 'o', ...
                'Color', year_colors(date_idx,:), 'MarkerSize', 4, 'MarkerFaceColor', year_colors(date_idx,:), 'HandleVisibility', 'off');
        end
    end
end
title(ax3, 'Discharge V-Q Curve: All Years');
xlabel(ax3, 'Capacity [Ah]');
ylabel(ax3, 'Voltage [V]');
legend(ax3, 'Location', 'best', 'NumColumns', 2);

saveas(fig_discharge_all, 'FieldQmax_dQdV_Discharge_AllYears.fig');
saveas(fig_discharge_all, 'FieldQmax_dQdV_Discharge_AllYears.png');

%% Step 4: Extract Peak Positions, Intensities, and Capacity
fprintf('\n=== Step 4: Extracting Peak Positions, Intensities, and Capacity ===\n');

% Function to find peaks in dQ/dV data
find_peaks_in_dQdV = @(V, dQdV, min_peak_height, min_peak_distance) ...
    find_peaks_local(V, dQdV, min_peak_height, min_peak_distance);

% Peak detection parameters
min_peak_height_chg = 5;      % Minimum peak height for charge (Ah/V)
min_peak_height_dischg = 5;   % Minimum peak height for discharge (Ah/V)
min_peak_distance = 0.05;     % Minimum distance between peaks (V)

% Extract charge peaks
fprintf('\n--- Charge Peaks ---\n');
charge_peak_table = [];
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Charge_data = All_Charge_dQdV.(sprintf('date_%d', date_idx));
    
    for k = 1:length(Charge_data)
        % Use common window data if available
        if isfield(Charge_data(k), 'V_dQdV_common') && ~isempty(Charge_data(k).V_dQdV_common)
            V_plot = Charge_data(k).V_dQdV_common;
            dQdV_plot = Charge_data(k).dQdV_common;
        elseif ~isempty(Charge_data(k).V_dQdV)
            V_plot = Charge_data(k).V_dQdV;
            dQdV_plot = Charge_data(k).dQdV;
        else
            continue;
        end
        
        % Find peaks
        [peak_V, peak_dQdV] = find_peaks_in_dQdV(V_plot, dQdV_plot, min_peak_height_chg, min_peak_distance);
        
        % Store peaks
        for p = 1:length(peak_V)
            charge_peak_table = [charge_peak_table; ...
                {date_str, k, peak_V(p), peak_dQdV(p)}];
        end
    end
end

% Display charge peak table
if ~isempty(charge_peak_table)
    fprintf('\nCharge dQ/dV Peaks:\n');
    fprintf('%-12s %-6s %-10s %-12s\n', 'Date', 'Seg', 'Voltage [V]', 'dQ/dV [Ah/V]');
    fprintf('%s\n', repmat('-', 1, 45));
    for i = 1:size(charge_peak_table, 1)
        fprintf('%-12s %-6d %-10.3f %-12.3f\n', ...
            charge_peak_table{i,1}, charge_peak_table{i,2}, ...
            charge_peak_table{i,3}, charge_peak_table{i,4});
    end
else
    fprintf('No charge peaks found.\n');
end

% Extract discharge peaks
fprintf('\n--- Discharge Peaks ---\n');
discharge_peak_table = [];
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Discharge_data = All_Discharge_dQdV.(sprintf('date_%d', date_idx));
    
    for k = 1:length(Discharge_data)
        % Use common window data if available
        if isfield(Discharge_data(k), 'V_dQdV_common') && ~isempty(Discharge_data(k).V_dQdV_common)
            V_plot = Discharge_data(k).V_dQdV_common;
            dQdV_plot = Discharge_data(k).dQdV_common;
        elseif ~isempty(Discharge_data(k).V_dQdV)
            V_plot = Discharge_data(k).V_dQdV;
            dQdV_plot = Discharge_data(k).dQdV;
        else
            continue;
        end
        
        % Find peaks
        [peak_V, peak_dQdV] = find_peaks_in_dQdV(V_plot, dQdV_plot, min_peak_height_dischg, min_peak_distance);
        
        % Store peaks
        for p = 1:length(peak_V)
            discharge_peak_table = [discharge_peak_table; ...
                {date_str, k, peak_V(p), peak_dQdV(p)}];
        end
    end
end

% Display discharge peak table
if ~isempty(discharge_peak_table)
    fprintf('\nDischarge dQ/dV Peaks:\n');
    fprintf('%-12s %-6s %-10s %-12s\n', 'Date', 'Seg', 'Voltage [V]', 'dQ/dV [Ah/V]');
    fprintf('%s\n', repmat('-', 1, 45));
    for i = 1:size(discharge_peak_table, 1)
        fprintf('%-12s %-6d %-10.3f %-12.3f\n', ...
            discharge_peak_table{i,1}, discharge_peak_table{i,2}, ...
            discharge_peak_table{i,3}, discharge_peak_table{i,4});
    end
else
    fprintf('No discharge peaks found.\n');
end

% Display capacity extraction results
fprintf('\n--- Capacity Extraction (%.2f-%.2f V) ---\n', capacity_V_min, capacity_V_max);
fprintf('Charge Capacity:\n');
fprintf('%-12s %-6s %-15s\n', 'Date', 'Seg', 'Capacity [Ah]');
fprintf('%s\n', repmat('-', 1, 35));
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Charge_data = All_Charge_dQdV.(sprintf('date_%d', date_idx));
    for k = 1:length(Charge_data)
        if isfield(Charge_data(k), 'capacity_Ah') && ~isempty(Charge_data(k).capacity_Ah)
            if ~isnan(Charge_data(k).capacity_Ah)
                fprintf('%-12s %-6d %-15.3f\n', date_str, k, Charge_data(k).capacity_Ah);
            else
                fprintf('%-12s %-6d %-15s\n', date_str, k, 'NaN (no data)');
            end
        end
    end
end

fprintf('\nDischarge Capacity:\n');
fprintf('%-12s %-6s %-15s\n', 'Date', 'Seg', 'Capacity [Ah]');
fprintf('%s\n', repmat('-', 1, 35));
for date_idx = 1:length(dates)
    date_str = dates{date_idx};
    Discharge_data = All_Discharge_dQdV.(sprintf('date_%d', date_idx));
    for k = 1:length(Discharge_data)
        if isfield(Discharge_data(k), 'capacity_Ah') && ~isempty(Discharge_data(k).capacity_Ah)
            if ~isnan(Discharge_data(k).capacity_Ah)
                fprintf('%-12s %-6d %-15.3f\n', date_str, k, Discharge_data(k).capacity_Ah);
            else
                fprintf('%-12s %-6d %-15s\n', date_str, k, 'NaN (no data)');
            end
        end
    end
end

fprintf('\n=== dQ/dV Analysis Complete ===\n');

%% Local function: find contiguous true segments
function segs = local_find_segments(mask)
    mask = mask(:)';
    d = diff([false mask false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    segs = [starts(:) ends(:)];
end

%% Local function: find peaks in dQ/dV data
function [peak_V, peak_dQdV] = find_peaks_local(V, dQdV, min_peak_height, min_peak_distance)
    % Find local maxima (peaks) in dQ/dV curve
    % Inputs:
    %   V: voltage vector
    %   dQdV: dQ/dV vector
    %   min_peak_height: minimum peak height (Ah/V)
    %   min_peak_distance: minimum distance between peaks (V)
    % Outputs:
    %   peak_V: voltage positions of peaks
    %   peak_dQdV: dQ/dV values at peaks
    
    if isempty(V) || isempty(dQdV) || length(V) ~= length(dQdV)
        peak_V = [];
        peak_dQdV = [];
        return;
    end
    
    % Use MATLAB's findpeaks function if available
    if exist('findpeaks', 'file') == 2
        try
            % Calculate minimum peak distance in samples
            if length(V) > 1
                dV_avg = mean(abs(diff(V)));
                min_peak_distance_samples = max(1, round(min_peak_distance / dV_avg));
            else
                min_peak_distance_samples = 1;
            end
            
            [peak_dQdV, peak_idx] = findpeaks(dQdV, 'MinPeakHeight', min_peak_height, ...
                'MinPeakDistance', min_peak_distance_samples);
            peak_V = V(peak_idx);
        catch
            % Fallback to simple peak detection
            peak_V = [];
            peak_dQdV = [];
            for i = 2:length(dQdV)-1
                if dQdV(i) > dQdV(i-1) && dQdV(i) > dQdV(i+1) && dQdV(i) >= min_peak_height
                    % Check distance from previous peak
                    if isempty(peak_V) || abs(V(i) - peak_V(end)) >= min_peak_distance
                        peak_V = [peak_V; V(i)];
                        peak_dQdV = [peak_dQdV; dQdV(i)];
                    end
                end
            end
        end
    else
        % Fallback: simple local maximum detection
        peak_V = [];
        peak_dQdV = [];
        for i = 2:length(dQdV)-1
            if dQdV(i) > dQdV(i-1) && dQdV(i) > dQdV(i+1) && dQdV(i) >= min_peak_height
                % Check distance from previous peak
                if isempty(peak_V) || abs(V(i) - peak_V(end)) >= min_peak_distance
                    peak_V = [peak_V; V(i)];
                    peak_dQdV = [peak_dQdV; dQdV(i)];
                end
            end
        end
    end
end
