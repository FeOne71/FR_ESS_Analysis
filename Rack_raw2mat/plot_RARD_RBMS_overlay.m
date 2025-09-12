%% RARD와 RBMS 데이터 오버레이 플롯 스크립트
% RARD 데이터: RARDSync/2021/202106/Raw_20210601.mat -> Raw.Rack01_Module01
% RBMS 데이터: 2021/202106/Raw_20210601.mat -> Raw.Rack01.DCCurrent_A

clear; clc; close all;

%% 경로 설정
base_path = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
rard_sync_path = fullfile(base_path, 'RARDSync');
rbms_path = fullfile(base_path, '2021', '202106');

%% 파일 경로 설정
day_str = '20210607';
rard_file = fullfile(rard_sync_path, '2021', '202106', sprintf('Raw_%s.mat', day_str));
rbms_file = fullfile(rbms_path, sprintf('Raw_%s.mat', day_str));

fprintf('Loading RARD data from: %s\n', rard_file);
fprintf('Loading RBMS data from: %s\n', rbms_file);

%% RARD 데이터 로드
if exist(rard_file, 'file')
    rard_data = load(rard_file);
    if isfield(rard_data.Raw, 'Rack01_Module01')
        rard_table = rard_data.Raw.Rack01_Module01;
        fprintf('RARD data loaded: %d rows, %d columns\n', height(rard_table), width(rard_table));
    else
        error('RARD data: Rack01_Module01 not found');
    end
else
    error('RARD file not found: %s', rard_file);
end

%% RBMS 데이터 로드
if exist(rbms_file, 'file')
    rbms_data = load(rbms_file);
    fprintf('RBMS data structure:\n');
    disp(fieldnames(rbms_data));
    
    if isfield(rbms_data, 'Raw_Rack') && isfield(rbms_data.Raw_Rack, 'Rack01')
        rbms_table = rbms_data.Raw_Rack.Rack01;
        fprintf('RBMS data loaded: %d rows, %d columns\n', height(rbms_table), width(rbms_table));
    elseif isfield(rbms_data, 'Raw') && isfield(rbms_data.Raw, 'Rack01')
        rbms_table = rbms_data.Raw.Rack01;
        fprintf('RBMS data loaded: %d rows, %d columns\n', height(rbms_table), width(rbms_table));
    else
        error('RBMS data: Rack01 not found in Raw_Rack or Raw');
    end
else
    error('RBMS file not found: %s', rbms_file);
end

%% 시간 동기화
% RARD 데이터 컬럼명 확인
fprintf('RARD table columns:\n');
disp(rard_table.Properties.VariableNames');

% RARD 시간 (RBMS_Time 컬럼 사용)
if ismember('RBMS_Time', rard_table.Properties.VariableNames)
    rard_time = rard_table.RBMS_Time;
    fprintf('Using RBMS_Time from RARD data\n');
elseif ismember('Time', rard_table.Properties.VariableNames)
    rard_time = rard_table.Time;
    fprintf('Using Time from RARD data\n');
else
    error('No time column found in RARD data');
end

% RBMS 시간
fprintf('RBMS data type: %s\n', class(rbms_table));
if istable(rbms_table)
    fprintf('RBMS table columns:\n');
    disp(rbms_table.Properties.VariableNames');
    if ismember('Time', rbms_table.Properties.VariableNames)
        rbms_time = rbms_table.Time;
    else
        error('No Time column found in RBMS data');
    end
elseif isstruct(rbms_table)
    fprintf('RBMS struct fields:\n');
    disp(fieldnames(rbms_table));
    if isfield(rbms_table, 'Time')
        rbms_time = rbms_table.Time;
    else
        error('No Time field found in RBMS data');
    end
else
    error('RBMS data is neither table nor struct');
end

% 시간 데이터 타입 통일
fprintf('RARD time type: %s\n', class(rard_time));
fprintf('RBMS time type: %s\n', class(rbms_time));

% RBMS 시간을 datetime으로 변환
if isstring(rbms_time) || ischar(rbms_time)
    rbms_time = datetime(rbms_time);
    fprintf('Converted RBMS time to datetime\n');
end

% 전체 데이터 사용 (공통 시간 구간 찾기 제거)
fprintf('Using all data points\n');
fprintf('RARD data points: %d\n', length(rard_time));
fprintf('RBMS data points: %d\n', length(rbms_time));

%% 데이터 추출
% RARD 데이터 (Cell#01~Cell#03) - 전체 데이터
rard_cells = [];
for i = 1:3
    cell_col = sprintf('Cell#%02d(V)', i);
    if ismember(cell_col, rard_table.Properties.VariableNames)
        rard_cells(:, i) = rard_table.(cell_col);
    else
        fprintf('Warning: %s not found in RARD data\n', cell_col);
        rard_cells(:, i) = NaN(length(rard_time), 1);
    end
end

% RBMS 데이터 (전압 데이터 - AverageCV_V) - 전체 데이터
if istable(rbms_table)
    if ismember('AverageCV_V', rbms_table.Properties.VariableNames)
        rbms_voltage = rbms_table.AverageCV_V;
        voltage_name = 'AverageCV_V';
    else
        error('AverageCV_V not found in RBMS table');
    end
elseif isstruct(rbms_table)
    if isfield(rbms_table, 'AverageCV_V')
        rbms_voltage = rbms_table.AverageCV_V;
        voltage_name = 'AverageCV_V';
    else
        error('AverageCV_V not found in RBMS struct');
    end
end

%% 플롯 생성
% 데이터 추출
rbms_current = rbms_table.DCCurrent_A;
rbms_power = rbms_table.DCPower_kW;
rard_rbms_voltage = rard_table.RBMS_AverageCV_V;
rard_rbms_current = rard_table.RBMS_DCCurrent_A;
rard_rbms_power = rard_table.RBMS_DCPower_kW;

% Figure 1: 전압
figure('Position', [100, 100, 1200, 600]);
hold on;
plot(rbms_time, rbms_voltage, 'Color', '#B0B0B0', 'LineWidth', 3, 'DisplayName', sprintf('RBMS %s', voltage_name));
plot(rard_time, rard_rbms_voltage, 'o-', 'Color', '#0073C2', 'MarkerFaceColor', '#0073C2', 'MarkerSize', 4, 'DisplayName', 'RARD Avg Cell Voltage [V]');
xlabel('Time'); ylabel('Voltage (V)'); ylim([2.7, 4.5]);
title(sprintf('RBMS vs RARD Voltage - Rack01 - %s', day_str));
legend('Location', 'southeast'); grid on;

% Figure 2: 전류
figure('Position', [200, 200, 1200, 600]);
hold on;
plot(rbms_time, rbms_current, 'Color', '#B0B0B0', 'LineWidth', 3, 'DisplayName', 'RBMS DC Current [A]');
plot(rard_time, rard_rbms_current, 'o-', 'Color', '#E18727', 'MarkerFaceColor', '#E18727', 'MarkerSize', 4, 'DisplayName', 'Synced DC Current [A]');
xlabel('Time'); ylabel('Current (A)');
title(sprintf('RBMS vs RARD Current - Rack01 - %s', day_str));
legend('Location', 'southeast'); grid on;

% Figure 3: 출력
figure('Position', [300, 300, 1200, 600]);
hold on;
plot(rbms_time, rbms_power, 'Color', '#B0B0B0', 'LineWidth', 3, 'DisplayName', 'RBMS DC Power [kW]');
plot(rard_time, rard_rbms_power, 'o-', 'Color', '#925E9F', 'MarkerFaceColor', '#925E9F', 'MarkerSize', 4, 'DisplayName', 'Synced DC Power [kW]');
xlabel('Time'); ylabel('Power (kW)');
title(sprintf('RBMS vs RARD Power - Rack01 - %s', day_str));
legend('Location', 'southeast'); grid on;

%% 통계 정보 출력
fprintf('\n=== Data Statistics ===\n');
fprintf('RBMS time range: %s to %s\n', datestr(min(rbms_time)), datestr(max(rbms_time)));
fprintf('RARD time range: %s to %s\n', datestr(min(rard_time)), datestr(max(rard_time)));
fprintf('RBMS data points: %d\n', length(rbms_time));
fprintf('RARD data points: %d\n', length(rard_time));
fprintf('RARD Cell voltage range: %.3f - %.3f V\n', min(rard_cells(:)), max(rard_cells(:)));
fprintf('RBMS Voltage range: %.3f - %.3f V\n', min(rbms_voltage), max(rbms_voltage));

fprintf('\nPlots created successfully!\n');
