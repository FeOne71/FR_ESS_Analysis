%% 필드 데이터 열화 지표 추출 스크립트
% 목적: 충방전 이벤트에서 열화 모드 정량화 지표 추출
% 방법: 동적 부하 응답 기반 (OCV 사용 안함)

clear; clc; close all;

%% 경로 설정
dataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old\2021\202106';
savePath = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_0113';

fprintf('========================================\n');
fprintf('  열화 지표 추출: 이벤트 기반 분석\n');
fprintf('========================================\n\n');

%% 파라미터 설정
Cnom = 64;  % Cell 용량
rackCapacity = 128;  % Rack 용량

% 이벤트 검출 조건 (랩 스크립트와 동일)
idle_threshold = 0.01 * 128;  % 0.01C = 1.28A (Rack 기준)
min_duration = 30;              % 최소 지속시간 30초
max_I_std = 1.5;                % 전류 표준편차 1.5A
std_check_window = [3, 30];     % 3~30초 구간에서 std 계산

fprintf('\\n=== 필드 데이터 이벤트 추출 및 열화 지표 계산 ===\\n\\n');

%% 샘플 파일 로드 (첫 번째 파일로 테스트)
matFiles = dir(fullfile(dataPath, '*.mat'));
fprintf('총 파일 수: %d개\\n\\n', length(matFiles));

% 첫 3일 데이터로 테스트
numDays = min(3, length(matFiles));
fprintf('분석 대상: %d일치 데이터\\n\\n', numDays);

% 전체 데이터 통합
allData = struct();
allData.time = [];
allData.voltage = [];
allData.current = [];
allData.soc = [];
allData.temperature = [];
allData.dayIndex = [];

for day = 1:numDays
    fprintf('로딩: Day %d/%d...\\n', day, numDays);
    
    filePath = fullfile(dataPath, matFiles(day).name);
    data = load(filePath);
    rack01 = data.Raw.Rack01;
    
    % 데이터 추가
    allData.time = [allData.time; rack01.Time];
    allData.voltage = [allData.voltage; rack01.AverageCV_V];
    allData.current = [allData.current; rack01.DCCurrent_A];
    allData.soc = [allData.soc; rack01.SOCPct];
    allData.temperature = [allData.temperature; rack01.AverageMT_degC];
    allData.dayIndex = [allData.dayIndex; ones(length(rack01.SOCPct), 1) * day];
end

fprintf('\\n총 데이터 포인트: %d개\\n', length(allData.current));
fprintf('기간: %.2f일\\n\\n', length(allData.current) / 86400);

%% 이벤트 검출
fprintf('[1단계] 이벤트 검출 중...\\n');

% Idle/Active 판별
is_idle = abs(allData.current) < idle_threshold;
is_active = ~is_idle;

% Idle -> Active 전환점
idle_to_active = find(is_idle(1:end-1) & is_active(2:end));

fprintf('  Idle->Active 전환: %d개\\n', length(idle_to_active));

% 이벤트 구조체 배열
events = struct();
events.charge = [];
events.discharge = [];

eventCount_charge = 0;
eventCount_discharge = 0;
filtered_duration = 0;
filtered_std = 0;

for k = 1:length(idle_to_active)
    idx_start = idle_to_active(k) + 1;
    
    % 이벤트 타입 판별
    event_type = sign(allData.current(idx_start));
    if event_type == 0
        continue;
    end
    
    % Active 구간 끝 찾기
    idx_end = idx_start;
    while idx_end <= length(allData.current)
        if is_active(idx_end) && sign(allData.current(idx_end)) == event_type
            idx_end = idx_end + 1;
        else
            break;
        end
    end
    idx_end = idx_end - 1;
    
    % 지속시간 체크 (30초 이상)
    duration_sec = idx_end - idx_start + 1;  % 1초 샘플링 가정
    if duration_sec < min_duration
        filtered_duration = filtered_duration + 1;
        continue;
    end
    
    % 3~30초 구간에서 전류 표준편차 체크
    check_start = idx_start + std_check_window(1);
    check_end = min(idx_start + std_check_window(2), idx_end);
    
    if check_end > check_start
        I_check = allData.current(check_start:check_end);
        I_std = std(I_check);
        
        if I_std > max_I_std
            filtered_std = filtered_std + 1;
            continue;
        end
    end
    
    % 이벤트 저장
    evt = struct();
    evt.start_idx = idx_start;
    evt.end_idx = idx_end;
    evt.duration = duration_sec;
    evt.type = event_type;
    evt.current = allData.current(idx_start:idx_end);
    evt.voltage = allData.voltage(idx_start:idx_end);
    evt.soc = allData.soc(idx_start:idx_end);
    evt.temperature = allData.temperature(idx_start:idx_end);
    
    if event_type > 0
        eventCount_charge = eventCount_charge + 1;
        events.charge = [events.charge; evt];
    else
        eventCount_discharge = eventCount_discharge + 1;
        events.discharge = [events.discharge; evt];
    end
end

fprintf('\\n이벤트 검출 완료:\\n');
fprintf('  충전 이벤트: %d개\\n', eventCount_charge);
fprintf('  방전 이벤트: %d개\\n', eventCount_discharge);
fprintf('  필터링 (지속시간): %d개\\n', filtered_duration);
fprintf('  필터링 (전류 std): %d개\\n\\n', filtered_std);

%% 열화 지표 추출
fprintf('[2단계] 열화 지표 추출 중...\\n\\n');

% 지표 저장 구조체
features = struct();
features.charge = struct();
features.discharge = struct();

% 충전 이벤트 분석
if eventCount_charge > 0
    fprintf('=== 충전 이벤트 분석 (%d개) ===\\n', eventCount_charge);
    
    % 각 이벤트별 지표 계산
    for i = 1:length(events.charge)
        evt = events.charge(i);
        
        % 1. 평균값
        features.charge(i).mean_current = mean(evt.current);
        features.charge(i).mean_voltage = mean(evt.voltage);
        features.charge(i).mean_soc = mean(evt.soc);
        features.charge(i).mean_temp = mean(evt.temperature);
        
        % 2. 전압 강하 (초기 -> 종료)
        features.charge(i).voltage_drop = evt.voltage(1) - evt.voltage(end);
        
        % 3. 저항 추정 (Step response)
        % 시작 후 1초 전압 변화
        if length(evt.voltage) >= 2
            V_before = evt.voltage(1);
            V_after_1s = evt.voltage(min(2, length(evt.voltage)));
            I_avg = mean(evt.current(1:min(2, length(evt.current))));
            
            if abs(I_avg) > 0.1
                features.charge(i).R_step = abs(V_after_1s - V_before) / abs(I_avg) * 1000;  % mOhm
            else
                features.charge(i).R_step = NaN;
            end
        else
            features.charge(i).R_step = NaN;
        end
        
        % 4. 온도 상승률 (°C/분)
        if length(evt.temperature) > 1
            temp_rise = evt.temperature(end) - evt.temperature(1);
            features.charge(i).temp_rise_rate = temp_rise / (evt.duration / 60);
        else
            features.charge(i).temp_rise_rate = 0;
        end
        
        % 5. 에너지 (Wh) - Cell 기준
        cell_current = evt.current / 2;  % 2P 구조
        cell_voltage = evt.voltage;
        power_W = cell_current .* cell_voltage;
        features.charge(i).energy_Wh = sum(power_W) / 3600;  % 1초 샘플
        
        % 6. Coulomb (Ah) - Cell 기준
        features.charge(i).coulomb_Ah = sum(cell_current) / 3600;
        
        % 7. C-rate
        features.charge(i).C_rate = mean(abs(cell_current)) / Cnom;
        
        % 8. 이벤트 정보
        features.charge(i).duration = evt.duration;
        features.charge(i).start_idx = evt.start_idx;
    end
    
    % 통계 출력
    fprintf('  평균 전압: %.4f V\\n', mean([features.charge.mean_voltage]));
    fprintf('  평균 전류: %.2f A\\n', mean([features.charge.mean_current]));
    fprintf('  평균 저항: %.2f mOhm\\n', mean([features.charge.R_step], 'omitnan'));
    fprintf('  평균 온도 상승률: %.4f °C/min\\n', mean([features.charge.temp_rise_rate]));
    fprintf('  평균 C-rate: %.4f C\\n\\n', mean([features.charge.C_rate]));
end

% 방전 이벤트 분석
if eventCount_discharge > 0
    fprintf('=== 방전 이벤트 분석 (%d개) ===\\n', eventCount_discharge);
    
    for i = 1:length(events.discharge)
        evt = events.discharge(i);
        
        % 동일한 지표 계산
        features.discharge(i).mean_current = mean(evt.current);
        features.discharge(i).mean_voltage = mean(evt.voltage);
        features.discharge(i).mean_soc = mean(evt.soc);
        features.discharge(i).mean_temp = mean(evt.temperature);
        
        features.discharge(i).voltage_drop = evt.voltage(1) - evt.voltage(end);
        
        if length(evt.voltage) >= 2
            V_before = evt.voltage(1);
            V_after_1s = evt.voltage(min(2, length(evt.voltage)));
            I_avg = mean(evt.current(1:min(2, length(evt.current))));
            
            if abs(I_avg) > 0.1
                features.discharge(i).R_step = abs(V_after_1s - V_before) / abs(I_avg) * 1000;
            else
                features.discharge(i).R_step = NaN;
            end
        else
            features.discharge(i).R_step = NaN;
        end
        
        if length(evt.temperature) > 1
            temp_rise = evt.temperature(end) - evt.temperature(1);
            features.discharge(i).temp_rise_rate = temp_rise / (evt.duration / 60);
        else
            features.discharge(i).temp_rise_rate = 0;
        end
        
        cell_current = evt.current / 2;
        cell_voltage = evt.voltage;
        power_W = cell_current .* cell_voltage;
        features.discharge(i).energy_Wh = sum(power_W) / 3600;
        features.discharge(i).coulomb_Ah = sum(cell_current) / 3600;
        features.discharge(i).C_rate = mean(abs(cell_current)) / Cnom;
        features.discharge(i).duration = evt.duration;
        features.discharge(i).start_idx = evt.start_idx;
    end
    
    fprintf('  평균 전압: %.4f V\\n', mean([features.discharge.mean_voltage]));
    fprintf('  평균 전류: %.2f A\\n', mean([features.discharge.mean_current]));
    fprintf('  평균 저항: %.2f mOhm\\n', mean([features.discharge.R_step], 'omitnan'));
    fprintf('  평균 온도 상승률: %.4f °C/min\\n', mean([features.discharge.temp_rise_rate]));
    fprintf('  평균 C-rate: %.4f C\\n\\n', mean([features.discharge.C_rate]));
end

%% SOC 구간별 분석 (동일 조건 비교용)
fprintf('[3단계] SOC 구간별 저항 분석\\n\\n');

soc_bins = [30, 40, 50, 60, 70];  % 분석할 SOC 구간
soc_tolerance = 2.5;  % ±2.5%

resistance_by_soc = struct();

for soc_target = soc_bins
    fprintf('SOC %.0f%% 구간 분석...\\n', soc_target);
    
    % 충전 이벤트
    if eventCount_charge > 0
        valid_charge = [];
        for i = 1:length(features.charge)
            soc_mean = features.charge(i).mean_soc;
            if abs(soc_mean - soc_target) <= soc_tolerance
                valid_charge = [valid_charge; features.charge(i).R_step];
            end
        end
        
        if ~isempty(valid_charge)
            resistance_by_soc.charge.(sprintf('SOC%d', soc_target)) = ...
                mean(valid_charge, 'omitnan');
            fprintf('  충전: %.2f mOhm (N=%d)\\n', ...
                    resistance_by_soc.charge.(sprintf('SOC%d', soc_target)), ...
                    sum(~isnan(valid_charge)));
        else
            resistance_by_soc.charge.(sprintf('SOC%d', soc_target)) = NaN;
            fprintf('  충전: 데이터 없음\\n');
        end
    end
    
    % 방전 이벤트
    if eventCount_discharge > 0
        valid_discharge = [];
        for i = 1:length(features.discharge)
            soc_mean = features.discharge(i).mean_soc;
            if abs(soc_mean - soc_target) <= soc_tolerance
                valid_discharge = [valid_discharge; features.discharge(i).R_step];
            end
        end
        
        if ~isempty(valid_discharge)
            resistance_by_soc.discharge.(sprintf('SOC%d', soc_target)) = ...
                mean(valid_discharge, 'omitnan');
            fprintf('  방전: %.2f mOhm (N=%d)\\n', ...
                    resistance_by_soc.discharge.(sprintf('SOC%d', soc_target)), ...
                    sum(~isnan(valid_discharge)));
        else
            resistance_by_soc.discharge.(sprintf('SOC%d', soc_target)) = NaN;
            fprintf('  방전: 데이터 없음\\n');
        end
    end
    
    fprintf('\\n');
end

%% 결과 저장
fprintf('[4단계] 결과 저장 중...\\n');

results = struct();
results.events = events;
results.features = features;
results.resistance_by_soc = resistance_by_soc;
results.parameters.idle_threshold = idle_threshold;
results.parameters.min_duration = min_duration;
results.parameters.max_I_std = max_I_std;
results.parameters.numDays = numDays;

save(fullfile(savePath, 'Field_Events_Features_2021Jun.mat'), 'results');
fprintf('저장 완료: Field_Events_Features_2021Jun.mat\\n\\n');

%% 요약 리포트
fprintf('========================================\\n');
fprintf('         분석 완료 요약\\n');
fprintf('========================================\\n\\n');

fprintf('[이벤트 통계]\\n');
fprintf('  분석 기간: %d일\\n', numDays);
fprintf('  충전 이벤트: %d개\\n', eventCount_charge);
fprintf('  방전 이벤트: %d개\\n\\n', eventCount_discharge);

fprintf('[충전 이벤트 평균 지표]\\n');
if eventCount_charge > 0
    fprintf('  전압: %.4f V\\n', mean([features.charge.mean_voltage]));
    fprintf('  전류: %.2f A\\n', mean([features.charge.mean_current]));
    fprintf('  저항: %.2f mOhm\\n', mean([features.charge.R_step], 'omitnan'));
    fprintf('  온도 상승: %.4f °C/min\\n', mean([features.charge.temp_rise_rate]));
    fprintf('  C-rate: %.4f C\\n\\n', mean([features.charge.C_rate]));
end

fprintf('[방전 이벤트 평균 지표]\\n');
if eventCount_discharge > 0
    fprintf('  전압: %.4f V\\n', mean([features.discharge.mean_voltage]));
    fprintf('  전류: %.2f A\\n', mean([features.discharge.mean_current]));
    fprintf('  저항: %.2f mOhm\\n', mean([features.discharge.R_step], 'omitnan'));
    fprintf('  온도 상승: %.4f °C/min\\n', mean([features.discharge.temp_rise_rate]));
    fprintf('  C-rate: %.4f C\\n\\n', mean([features.discharge.C_rate]));
end

fprintf('========================================\\n');
