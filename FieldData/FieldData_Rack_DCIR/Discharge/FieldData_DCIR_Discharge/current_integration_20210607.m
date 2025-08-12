% DCCurrent_A 적산 계산 코드
% 시간 범위: 2021-06-07 13:33:32 to 2021-06-07 14:07:01

% 데이터 로드
load('Raw_20210607.m');

% 시작 및 종료 시간 설정
start_time = datetime('2021-06-07 13:33:32', 'Format', 'yyyy-MM-dd HH:mm:ss');
end_time = datetime('2021-06-07 14:07:01', 'Format', 'yyyy-MM-dd HH:mm:ss');

% 시간 범위에 해당하는 데이터 인덱스 찾기 (Raw 구조체의 Time 필드 사용)
time_idx = (Raw.Time >= start_time) & (Raw.Time <= end_time);

% 해당 시간 범위의 데이터 추출
selected_time = Raw.Time(time_idx);
selected_current = Raw.DCCurrent_A(time_idx);

% 적산 계산 (trapezoidal rule 사용)
if length(selected_time) > 1
    % 시간을 초 단위로 변환
    time_seconds = seconds(selected_time - selected_time(1));
    
    % 적산 계산 (trapezoidal integration)
    current_integral = trapz(time_seconds, selected_current);
    
    % 결과 출력
    fprintf('=== DCCurrent_A 적산 결과 ===\n');
    fprintf('시작 시간: %s\n', datestr(start_time));
    fprintf('종료 시간: %s\n', datestr(end_time));
    fprintf('총 시간: %.2f 분\n', minutes(end_time - start_time));
    fprintf('적산 값: %.4f A·s\n', current_integral);
    fprintf('적산 값 (A·h): %.6f A·h\n', current_integral/3600);
    
    % 그래프 그리기
    figure('Name', 'DCCurrent_A 적산 분석', 'Position', [100, 100, 1200, 800]);
    
    % 서브플롯 1: 전류 vs 시간
    subplot(2,2,1);
    plot(selected_time, selected_current, 'b-', 'LineWidth', 1.5);
    title('DCCurrent_A vs Time');
    xlabel('Time');
    ylabel('Current (A)');
    grid on;
    
    % 서브플롯 2: 적산 누적 그래프
    subplot(2,2,2);
    cumulative_integral = cumtrapz(time_seconds, selected_current);
    plot(time_seconds/60, cumulative_integral, 'r-', 'LineWidth', 1.5);
    title('Cumulative Current Integration');
    xlabel('Time (minutes)');
    ylabel('Cumulative Integral (A·s)');
    grid on;
    
    % 서브플롯 3: 전압 vs 시간
    subplot(2,2,3);
    plot(selected_time, Raw.DCVoltage_V(time_idx), 'm-', 'LineWidth', 1.5);
    title('DCVoltage_V vs Time');
    xlabel('Time');
    ylabel('Voltage (V)');
    grid on;
    
    % 서브플롯 4: 통계 정보
    subplot(2,2,4);
    text(0.1, 0.8, sprintf('통계 정보:\n\n최대 전류: %.2f A\n최소 전류: %.2f A\n평균 전류: %.2f A\n표준편차: %.2f A\n\n적산 결과:\n%.4f A·s\n%.6f A·h', ...
        max(selected_current), min(selected_current), mean(selected_current), std(selected_current), ...
        current_integral, current_integral/3600), 'FontSize', 10, 'VerticalAlignment', 'top');
    axis off;
    
else
    fprintf('오류: 지정된 시간 범위에 데이터가 없습니다.\n');
end

% 추가 분석: 전류 변화율 계산
if length(selected_current) > 1
    current_rate_of_change = diff(selected_current) ./ diff(time_seconds);
    
    fprintf('\n=== 전류 변화율 분석 ===\n');
    fprintf('최대 변화율: %.4f A/s\n', max(current_rate_of_change));
    fprintf('최소 변화율: %.4f A/s\n', min(current_rate_of_change));
    fprintf('평균 변화율: %.4f A/s\n', mean(current_rate_of_change));
end 