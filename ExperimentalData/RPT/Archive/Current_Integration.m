%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current Integration Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 데이터 로드
folderPath = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data\RPT';
channel = 'Ch16';
rpt_cycle = '200cyc';

filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
filepath = fullfile(folderPath, filename);
T = readtable(filepath);

fprintf('=== Current Integration Analysis ===\n');
fprintf('Channel: %s, RPT: %s\n', channel, rpt_cycle);
fprintf('Total rows: %d\n', height(T));

%% Time을 초로 변환
fprintf('\n=== Time Conversion ===\n');
fprintf('Time column type: %s\n', class(T.Time));

% Time을 초로 변환
if isdatetime(T.Time)
    % datetime을 seconds로 변환 (시작점을 0으로)
    time_start = T.Time(1);
    time_seconds = seconds(T.Time - time_start);
    fprintf('Time converted to seconds (relative to start)\n');
elseif isduration(T.Time)
    % duration을 seconds로 변환
    time_seconds = seconds(T.Time);
    fprintf('Duration converted to seconds\n');
else
    fprintf('Unknown time format\n');
    return;
end

%% 지정된 시간 구간에서 전류 적분
time_ranges = [50810, 51084; 5248, 53603; 56375, 57367];

fprintf('\n=== Current Integration Results ===\n');
fprintf('Range\t\tStart\t\tEnd\t\tDuration\tCurrent Integral\tAverage Current\n');
fprintf('-----\t\t-----\t\t---\t\t--------\t----------------\t---------------\n');

for i = 1:size(time_ranges, 1)
    start_time = time_ranges(i, 1);
    end_time = time_ranges(i, 2);
    
    % 해당 시간 구간의 데이터 찾기
    mask = (time_seconds >= start_time) & (time_seconds <= end_time);
    
    if sum(mask) > 0
        % 시간과 전류 데이터 추출
        range_time = time_seconds(mask);
        range_current = T.Current_A_(mask);
        
        % 전류 적분
        current_integral = trapz(range_time, range_current);
        
        % 평균 전류
        avg_current = mean(range_current);
        
        % 지속 시간
        duration = end_time - start_time;
        
        fprintf('Range %d\t%.1f\t\t%.1f\t\t%.1f s\t\t%.6f Ah\t\t%.4f A\n', ...
            i, start_time, end_time, duration, current_integral, avg_current);
        
        % 상세 정보 출력
        fprintf('  Data points: %d\n', length(range_time));
        fprintf('  Current range: %.4f ~ %.4f A\n', min(range_current), max(range_current));
        fprintf('  Time step: %.2f ~ %.2f s\n', min(diff(range_time)), max(diff(range_time)));
        fprintf('\n');
    else
        fprintf('Range %d\t%.1f\t\t%.1f\t\tNo data found\n', i, start_time, end_time);
    end
end

%% 전체 데이터 통계
fprintf('\n=== Overall Statistics ===\n');
fprintf('Total time range: %.1f ~ %.1f seconds (%.1f hours)\n', ...
    min(time_seconds), max(time_seconds), (max(time_seconds) - min(time_seconds))/3600);
fprintf('Total current range: %.4f ~ %.4f A\n', min(T.Current_A_), max(T.Current_A_));
fprintf('Average current: %.4f A\n', mean(T.Current_A_));

%% 그래프 출력 (선택사항)
figure('Name', 'Current vs Time', 'Position', [100 100 1200 800]);

% 전체 데이터 플롯
subplot(2,1,1);
plot(time_seconds, T.Current_A_, 'b-', 'LineWidth', 0.5);
hold on;

% 지정된 구간들 하이라이트
colors = {'r', 'g', 'm'};
for i = 1:size(time_ranges, 1)
    start_time = time_ranges(i, 1);
    end_time = time_ranges(i, 2);
    mask = (time_seconds >= start_time) & (time_seconds <= end_time);
    
    if sum(mask) > 0
        plot(time_seconds(mask), T.Current_A_(mask), colors{i}, 'LineWidth', 2);
    end
end

xlabel('Time [s]');
ylabel('Current [A]');
title('Current vs Time (Highlighted Ranges)');
grid on;
legend('All Data', 'Range 1', 'Range 2', 'Range 3');

% 확대된 뷰
subplot(2,1,2);
for i = 1:size(time_ranges, 1)
    start_time = time_ranges(i, 1);
    end_time = time_ranges(i, 2);
    mask = (time_seconds >= start_time) & (time_seconds <= end_time);
    
    if sum(mask) > 0
        plot(time_seconds(mask), T.Current_A_(mask), colors{i}, 'LineWidth', 2);
        hold on;
    end
end

xlabel('Time [s]');
ylabel('Current [A]');
title('Current vs Time (Selected Ranges Only)');
grid on;
legend('Range 1', 'Range 2', 'Range 3');

fprintf('\n=== Analysis Complete ===\n'); 