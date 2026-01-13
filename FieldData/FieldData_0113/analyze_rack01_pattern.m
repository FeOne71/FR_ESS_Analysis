%% Rack01 운전 패턴 분석 (2021년 6월 전체)
% 목적: 실제 충방전 패턴 확인 및 열화 분석 가능성 검토

clear; clc; close all;

%% 경로 설정
dataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old\2021\202106';
savePath = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_0113';

fprintf('=== Rack01 운전 패턴 분석 시작 ===\n\n');

%% 1단계: 파일 목록 확인
matFiles = dir(fullfile(dataPath, '*.mat'));
fprintf('[1단계] 파일 로드 중...\n');
fprintf('총 파일 수: %d\n\n', length(matFiles));

%% 2단계: 데이터 통합
fprintf('[2단계] Rack01 데이터 통합 중...\n');

% 첫 파일로 구조 확인
sample = load(fullfile(dataPath, matFiles(1).name));
dataLength = length(sample.Raw.Rack01.Time);
totalLength = dataLength * length(matFiles);

fprintf('일별 데이터 길이: %d\n', dataLength);
fprintf('예상 총 길이: %d (%.1f일)\n\n', totalLength, totalLength/86400);

% 데이터 통합 구조체 초기화
Rack01_June = struct();
Rack01_June.Time = cell(totalLength, 1);
Rack01_June.Counter = zeros(totalLength, 1);
Rack01_June.SOC = zeros(totalLength, 1);
Rack01_June.SOH = zeros(totalLength, 1);
Rack01_June.Current = zeros(totalLength, 1);
Rack01_June.Power = zeros(totalLength, 1);
Rack01_June.Voltage_Avg = zeros(totalLength, 1);
Rack01_June.Voltage_High = zeros(totalLength, 1);
Rack01_June.Voltage_Low = zeros(totalLength, 1);
Rack01_June.Temp_Avg = zeros(totalLength, 1);
Rack01_June.Charge = cell(totalLength, 1);
Rack01_June.Status = cell(totalLength, 1);

% 파일별 데이터 로드 및 통합
idx_start = 1;
for i = 1:length(matFiles)
    if mod(i, 5) == 0
        fprintf('  진행중: %d/%d 파일\n', i, length(matFiles));
    end
    
    data = load(fullfile(dataPath, matFiles(i).name));
    rack = data.Raw.Rack01;
    
    idx_end = idx_start + length(rack.Time) - 1;
    
    Rack01_June.Time(idx_start:idx_end) = cellstr(rack.Time);
    Rack01_June.Counter(idx_start:idx_end) = rack.Counter;
    Rack01_June.SOC(idx_start:idx_end) = rack.SOCPct;
    Rack01_June.SOH(idx_start:idx_end) = rack.SOHPct;
    Rack01_June.Current(idx_start:idx_end) = rack.DCCurrent_A;
    Rack01_June.Power(idx_start:idx_end) = rack.DCPower_kW;
    Rack01_June.Voltage_Avg(idx_start:idx_end) = rack.AverageCV_V;
    Rack01_June.Voltage_High(idx_start:idx_end) = rack.HighestCV_V;
    Rack01_June.Voltage_Low(idx_start:idx_end) = rack.LowestCV_V;
    Rack01_June.Temp_Avg(idx_start:idx_end) = rack.AverageMT_degC;
    Rack01_June.Charge(idx_start:idx_end) = cellstr(rack.Charge);
    Rack01_June.Status(idx_start:idx_end) = cellstr(rack.Status);
    
    idx_start = idx_end + 1;
end

fprintf('통합 완료! 총 데이터 포인트: %d\n\n', length(Rack01_June.SOC));

%% 3단계: 기본 통계
fprintf('[3단계] 기본 통계 분석\n');
fprintf('─────────────────────────────────────\n');
fprintf('SOC:\n');
fprintf('  범위: %.2f%% ~ %.2f%%\n', min(Rack01_June.SOC), max(Rack01_June.SOC));
fprintf('  평균: %.2f%%\n', mean(Rack01_June.SOC));
fprintf('  표준편차: %.2f%%\n\n', std(Rack01_June.SOC));

fprintf('SOH:\n');
fprintf('  범위: %.2f%% ~ %.2f%%\n', min(Rack01_June.SOH), max(Rack01_June.SOH));
fprintf('  평균: %.2f%%\n\n', mean(Rack01_June.SOH));

fprintf('전류:\n');
fprintf('  범위: %.2f A ~ %.2f A\n', min(Rack01_June.Current), max(Rack01_June.Current));
fprintf('  평균: %.2f A\n', mean(Rack01_June.Current));
fprintf('  표준편차: %.2f A\n', std(Rack01_June.Current));
fprintf('  0 아닌 값 비율: %.2f%%\n\n', 100*sum(abs(Rack01_June.Current)>0.1)/length(Rack01_June.Current));

fprintf('전압:\n');
fprintf('  평균 범위: %.4f V ~ %.4f V\n', min(Rack01_June.Voltage_Avg), max(Rack01_June.Voltage_Avg));
fprintf('  평균: %.4f V\n', mean(Rack01_June.Voltage_Avg));
fprintf('  전압 편차 (High-Low) 평균: %.4f V\n\n', mean(Rack01_June.Voltage_High - Rack01_June.Voltage_Low));

fprintf('온도:\n');
fprintf('  범위: %.1f°C ~ %.1f°C\n', min(Rack01_June.Temp_Avg), max(Rack01_June.Temp_Avg));
fprintf('  평균: %.1f°C\n\n', mean(Rack01_June.Temp_Avg));

%% 4단계: 충방전 패턴 분석
fprintf('[4단계] 충방전 패턴 분석\n');
fprintf('─────────────────────────────────────\n');

% 충전/방전 구분
charge_idx = Rack01_June.Current > 0.1;  % 충전
discharge_idx = Rack01_June.Current < -0.1;  % 방전
idle_idx = abs(Rack01_June.Current) <= 0.1;  % 대기

fprintf('운전 시간 분포:\n');
fprintf('  충전: %.2f%% (%.1f시간)\n', 100*sum(charge_idx)/length(Rack01_June.Current), sum(charge_idx)/3600);
fprintf('  방전: %.2f%% (%.1f시간)\n', 100*sum(discharge_idx)/length(Rack01_June.Current), sum(discharge_idx)/3600);
fprintf('  대기: %.2f%% (%.1f시간)\n\n', 100*sum(idle_idx)/length(Rack01_June.Current), sum(idle_idx)/3600);

if sum(charge_idx) > 0
    fprintf('충전 시:\n');
    fprintf('  평균 전류: %.2f A\n', mean(Rack01_June.Current(charge_idx)));
    fprintf('  최대 전류: %.2f A\n', max(Rack01_June.Current(charge_idx)));
    fprintf('  평균 전력: %.2f kW\n\n', mean(Rack01_June.Power(charge_idx)));
end

if sum(discharge_idx) > 0
    fprintf('방전 시:\n');
    fprintf('  평균 전류: %.2f A\n', mean(Rack01_June.Current(discharge_idx)));
    fprintf('  최소 전류: %.2f A\n', min(Rack01_June.Current(discharge_idx)));
    fprintf('  평균 전력: %.2f kW\n\n', mean(Rack01_June.Power(discharge_idx)));
end

%% 5단계: 사이클 카운팅
fprintf('[5단계] 사이클 카운팅\n');
fprintf('─────────────────────────────────────\n');

% 충방전 전환 횟수
current_sign = sign(Rack01_June.Current);
transitions = sum(abs(diff(current_sign)) > 0);

fprintf('충방전 전환 횟수: %d회\n', transitions);

% SOC 변화 범위
SOC_range = max(Rack01_June.SOC) - min(Rack01_June.SOC);
fprintf('SOC 변화 범위: %.2f%%\n', SOC_range);

% 적산 전하량 (Ah)
Ah_charge = sum(Rack01_June.Current(charge_idx)) / 3600;  % 초 -> 시간
Ah_discharge = abs(sum(Rack01_June.Current(discharge_idx)) / 3600);

fprintf('적산 충전량: %.2f Ah\n', Ah_charge);
fprintf('적산 방전량: %.2f Ah\n', Ah_discharge);
fprintf('충방전 효율: %.2f%%\n\n', 100*Ah_discharge/Ah_charge);

%% 6단계: 전압-SOC 관계 분석
fprintf('[6단계] 전압-SOC 관계 분석\n');
fprintf('─────────────────────────────────────\n');

% SOC별 전압 분포
SOC_bins = 0:5:100;
V_mean_by_SOC = zeros(length(SOC_bins)-1, 1);
V_std_by_SOC = zeros(length(SOC_bins)-1, 1);

for i = 1:length(SOC_bins)-1
    idx = Rack01_June.SOC >= SOC_bins(i) & Rack01_June.SOC < SOC_bins(i+1);
    if sum(idx) > 0
        V_mean_by_SOC(i) = mean(Rack01_June.Voltage_Avg(idx));
        V_std_by_SOC(i) = std(Rack01_June.Voltage_Avg(idx));
    else
        V_mean_by_SOC(i) = NaN;
        V_std_by_SOC(i) = NaN;
    end
end

fprintf('SOC별 평균 전압:\n');
for i = 1:length(SOC_bins)-1
    if ~isnan(V_mean_by_SOC(i))
        fprintf('  SOC %d-%d%%: %.4f V (±%.4f)\n', SOC_bins(i), SOC_bins(i+1), V_mean_by_SOC(i), V_std_by_SOC(i));
    end
end
fprintf('\n');

%% 7단계: 저항 추정 가능성 검토
fprintf('[7단계] 저항 추정 가능성 검토\n');
fprintf('─────────────────────────────────────\n');

if sum(abs(Rack01_June.Current) > 10) > 100  % 10A 이상 구간이 100개 이상
    % 전류가 있는 구간만 추출
    active_idx = abs(Rack01_June.Current) > 10;
    
    I_active = Rack01_June.Current(active_idx);
    V_active = Rack01_June.Voltage_Avg(active_idx);
    
    % 간단한 선형 회귀: V = a - b*I (방전시 I<0)
    p = polyfit(I_active, V_active, 1);
    R_apparent = -p(1) * 1000;  % mOhm 단위
    
    fprintf('추정 겉보기 저항: %.2f mOhm\n', R_apparent);
    fprintf('※ 주의: OCV를 모르므로 정확한 저항이 아님\n\n');
else
    fprintf('전류가 충분히 크지 않아 저항 추정 불가능\n\n');
end

%% 8단계: 데이터 저장
fprintf('[8단계] 데이터 저장 중...\n');

save(fullfile(savePath, 'Rack01_June2021.mat'), 'Rack01_June', '-v7.3');

fprintf('저장 완료: Rack01_June2021.mat\n\n');

%% 9단계: 시각화
fprintf('[9단계] 시각화 생성 중...\n');

figure('Position', [100 100 1400 900]);

% 시간축 생성 (일 단위)
time_days = (1:length(Rack01_June.SOC)) / 86400;

% 서브플롯 1: SOC 추이
subplot(4,1,1);
plot(time_days, Rack01_June.SOC, 'b-', 'LineWidth', 1);
ylabel('SOC [%]');
title('Rack01 운전 패턴 분석 (2021년 6월)');
grid on;
xlim([0 30]);

% 서브플롯 2: 전류
subplot(4,1,2);
plot(time_days, Rack01_June.Current, 'r-', 'LineWidth', 0.5);
ylabel('전류 [A]');
grid on;
xlim([0 30]);
yline(0, 'k--', 'LineWidth', 1);

% 서브플롯 3: 전압
subplot(4,1,3);
plot(time_days, Rack01_June.Voltage_Avg, 'g-', 'LineWidth', 1);
hold on;
plot(time_days, Rack01_June.Voltage_High, 'g:', 'LineWidth', 0.5);
plot(time_days, Rack01_June.Voltage_Low, 'g:', 'LineWidth', 0.5);
ylabel('전압 [V]');
legend('평균', '최고', '최저', 'Location', 'best');
grid on;
xlim([0 30]);

% 서브플롯 4: 온도
subplot(4,1,4);
plot(time_days, Rack01_June.Temp_Avg, 'm-', 'LineWidth', 1);
xlabel('시간 [일]');
ylabel('온도 [°C]');
grid on;
xlim([0 30]);

saveas(gcf, fullfile(savePath, 'Rack01_운전패턴_202106.fig'));
fprintf('그림 저장: Rack01_운전패턴_202106.fig\n');

% 전압-SOC 관계 플롯
figure('Position', [150 150 800 600]);
scatter(Rack01_June.SOC, Rack01_June.Voltage_Avg, 1, Rack01_June.Current, 'filled');
colorbar;
xlabel('SOC [%]');
ylabel('평균 전압 [V]');
title('전압-SOC 관계 (색상: 전류)');
grid on;

saveas(gcf, fullfile(savePath, 'Rack01_전압SOC관계_202106.fig'));
fprintf('그림 저장: Rack01_전압SOC관계_202106.fig\n\n');

fprintf('=== 분석 완료 ===\n');
fprintf('워크스페이스에 Rack01_June 구조체가 로드되어 있습니다.\n');
