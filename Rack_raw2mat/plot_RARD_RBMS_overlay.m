%% RARD와 RBMS 데이터 오버레이 플롯 스크립트
% RARDsync 폴더의 mat 파일 사용 (RARD + 동기화된 RBMS 데이터 포함)
% 구조: Raw.Rack01.Module01 (RARD 데이터 + RBMS_* 필드)

clear; clc; close all;

%% 경로 설정
base_path = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
rard_sync_path = fullfile(base_path, 'RARDsync');

%% 파일 경로 설정
day_str = '20210607';
rard_file = fullfile(rard_sync_path, '2021', '202106', sprintf('Raw_%s.mat', day_str));

fprintf('Loading data from: %s\n', rard_file);

%% RARD 데이터 로드
if exist(rard_file, 'file')
    rard_data = load(rard_file);
    if isfield(rard_data.Raw, 'Rack01') && isfield(rard_data.Raw.Rack01, 'Module01')
        rard_module = rard_data.Raw.Rack01.Module01;
        fprintf('RARD data loaded: Module01 structure found\n');
        fprintf('RARD Module01 fields: %s\n', strjoin(fieldnames(rard_module), ', '));
    else
        error('RARD data: Rack01.Module01 not found');
    end
else
    error('RARD file not found: %s', rard_file);
end

%% 시간 데이터 추출
% RARD 시간 (동기화된 시간 사용)
if isfield(rard_module, 'Time')
    rard_time = rard_module.Time;
    fprintf('Using Time from RARD data\n');
else
    error('No Time field found in RARD data');
end

% 시간 데이터 타입 통일
fprintf('Time type: %s\n', class(rard_time));
if isstring(rard_time) || ischar(rard_time)
    rard_time = datetime(rard_time);
    fprintf('Converted time to datetime\n');
end

fprintf('Data points: %d\n', length(rard_time));

%% 데이터 추출
% RARD 데이터 (M1_Cell1~M1_Cell3) - 전체 데이터
rard_cells = [];
for i = 1:3
    cell_field = sprintf('M1_Cell%d', i);
    if isfield(rard_module, cell_field)
        rard_cells(:, i) = rard_module.(cell_field);
    else
        fprintf('Warning: %s not found in RARD data\n', cell_field);
        rard_cells(:, i) = NaN(length(rard_time), 1);
    end
end

%% 플롯 생성
% RARD에 동기화된 RBMS 데이터 추출
if isfield(rard_module, 'RBMS_AverageCV_V')
    rbms_voltage = rard_module.RBMS_AverageCV_V;
    voltage_name = 'AverageCV_V';
else
    error('RBMS_AverageCV_V not found in RARD module data');
end

if isfield(rard_module, 'RBMS_DCCurrent_A')
    rbms_current = rard_module.RBMS_DCCurrent_A;
else
    error('RBMS_DCCurrent_A not found in RARD module data');
end

if isfield(rard_module, 'RBMS_DCPower_kW')
    rbms_power = rard_module.RBMS_DCPower_kW;
else
    error('RBMS_DCPower_kW not found in RARD module data');
end

% Figure 1: 전압
figure('Position', [100, 100, 1200, 600]);
hold on;
plot(rard_time, rbms_voltage, 'Color', '#B0B0B0', 'LineWidth', 3, 'DisplayName', sprintf('RBMS %s', voltage_name));
% RARD 셀 전압 평균 계산
if size(rard_cells, 2) > 0
    rard_avg_voltage = mean(rard_cells, 2, 'omitnan');
    plot(rard_time, rard_avg_voltage, 'o-', 'Color', '#0073C2', 'MarkerFaceColor', '#0073C2', 'MarkerSize', 4, 'DisplayName', 'RARD Avg Cell Voltage [V]');
end
xlabel('Time'); ylabel('Voltage (V)'); ylim([2.7, 4.5]);
title(sprintf('RBMS vs RARD Voltage - Rack01 Module01 - %s', day_str));
legend('Location', 'southeast'); grid on;

% Figure 2: 전류
figure('Position', [200, 200, 1200, 600]);
hold on;
plot(rard_time, rbms_current, 'o-', 'Color', '#E18727', 'MarkerFaceColor', '#E18727', 'MarkerSize', 4, 'DisplayName', 'RBMS DC Current [A]');
xlabel('Time'); ylabel('Current (A)');
title(sprintf('RBMS Current - Rack01 Module01 - %s', day_str));
legend('Location', 'southeast'); grid on;

% Figure 3: 출력
figure('Position', [300, 300, 1200, 600]);
hold on;
plot(rard_time, rbms_power, 'o-', 'Color', '#925E9F', 'MarkerFaceColor', '#925E9F', 'MarkerSize', 4, 'DisplayName', 'RBMS DC Power [kW]');
xlabel('Time'); ylabel('Power (kW)');
title(sprintf('RBMS Power - Rack01 Module01 - %s', day_str));
legend('Location', 'southeast'); grid on;

%% 통계 정보 출력
fprintf('\n=== Data Statistics ===\n');
fprintf('Time range: %s to %s\n', datestr(min(rard_time)), datestr(max(rard_time)));
fprintf('Data points: %d\n', length(rard_time));
fprintf('RARD Cell voltage range: %.3f - %.3f V\n', min(rard_cells(:), [], 'omitnan'), max(rard_cells(:), [], 'omitnan'));
fprintf('RBMS Voltage range: %.3f - %.3f V\n', min(rbms_voltage, [], 'omitnan'), max(rbms_voltage, [], 'omitnan'));

fprintf('\nPlots created successfully!\n');
