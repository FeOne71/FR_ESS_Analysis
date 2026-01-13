%% Rack01 2021년 6월 전체 데이터 분석
% 목적: 열화 모드 정량화를 위한 운전 패턴 및 특성 분석
% Rack 구성: Cell(64Ah) -> Module(2P14S, 128Ah) -> Rack(17S, 128Ah)

clear; clc; close all;

%% 경로 설정
dataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old\2021\202106';
savePath = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_0113';

fprintf('========================================\n');
fprintf('  Rack01 2021년 6월 전체 데이터 분석\n');
fprintf('========================================\n\n');

%% 1단계: 데이터 로드
fprintf('[1단계] 데이터 로딩 중...\n');

matFiles = dir(fullfile(dataPath, '*.mat'));
numDays = length(matFiles);
fprintf('총 파일 수: %d개\n', numDays);

% 데이터 구조 초기화
allData = struct();
allData.Time = [];
allData.SOCPct = [];
allData.SOHPct = [];
allData.DCCurrent_A = [];
allData.DCPower_kW = [];
allData.SumCV_V = [];
allData.AverageCV_V = [];
allData.HighestCV_V = [];
allData.LowestCV_V = [];
allData.AverageMT_degC = [];
allData.HighestMT_degC = [];
allData.LowestMT_degC = [];
allData.Charge = {};
allData.Status = {};
allData.DayIndex = [];

% 각 파일 로드 및 통합
for i = 1:numDays
    if mod(i, 5) == 0
        fprintf('  진행: %d/%d 파일\n', i, numDays);
    end
    
    filePath = fullfile(dataPath, matFiles(i).name);
    data = load(filePath);
    rack01 = data.Raw.Rack01;
    
    % 데이터 추가
    allData.Time = [allData.Time; rack01.Time];
    allData.SOCPct = [allData.SOCPct; rack01.SOCPct];
    allData.SOHPct = [allData.SOHPct; rack01.SOHPct];
    allData.DCCurrent_A = [allData.DCCurrent_A; rack01.DCCurrent_A];
    allData.DCPower_kW = [allData.DCPower_kW; rack01.DCPower_kW];
    allData.SumCV_V = [allData.SumCV_V; rack01.SumCV_V];
    allData.AverageCV_V = [allData.AverageCV_V; rack01.AverageCV_V];
    allData.HighestCV_V = [allData.HighestCV_V; rack01.HighestCV_V];
    allData.LowestCV_V = [allData.LowestCV_V; rack01.LowestCV_V];
    allData.AverageMT_degC = [allData.AverageMT_degC; rack01.AverageMT_degC];
    allData.HighestMT_degC = [allData.HighestMT_degC; rack01.HighestMT_degC];
    allData.LowestMT_degC = [allData.LowestMT_degC; rack01.LowestMT_degC];
    allData.Charge = [allData.Charge; rack01.Charge];
    allData.Status = [allData.Status; rack01.Status];
    allData.DayIndex = [allData.DayIndex; ones(length(rack01.SOCPct), 1) * i];
end

fprintf('데이터 로딩 완료!\n');
fprintf('총 데이터 포인트: %d개\n\n', length(allData.SOCPct));

%% 2단계: 기본 통계
fprintf('[2단계] 기본 통계 분석\n');

% Rack 사양
cellCapacity_Ah = 64;
moduleCapacity_Ah = 128;  % 2P
rackCapacity_Ah = 128;     % 17S

% Cell 레벨 전류 계산 (Module이 2P 구조)
allData.CellCurrent_A = allData.DCCurrent_A / 2;

% C-rate 계산
allData.CRate_Cell = allData.CellCurrent_A / cellCapacity_Ah;
allData.CRate_Rack = allData.DCCurrent_A / rackCapacity_Ah;

fprintf('\\n=== SOC 통계 ===\n');
fprintf('  평균: %.2f%%\n', mean(allData.SOCPct, 'omitnan'));
fprintf('  최소: %.2f%%\n', min(allData.SOCPct));
fprintf('  최대: %.2f%%\n', max(allData.SOCPct));
fprintf('  표준편차: %.2f%%\n', std(allData.SOCPct, 'omitnan'));

fprintf('\\n=== SOH 통계 ===\n');
fprintf('  평균: %.2f%%\n', mean(allData.SOHPct, 'omitnan'));
fprintf('  최소: %.2f%%\n', min(allData.SOHPct));
fprintf('  최대: %.2f%%\n', max(allData.SOHPct));

fprintf('\\n=== 전류 통계 (Rack 레벨) ===\n');
fprintf('  평균: %.2f A\n', mean(allData.DCCurrent_A, 'omitnan'));
fprintf('  최소: %.2f A\n', min(allData.DCCurrent_A));
fprintf('  최대: %.2f A\n', max(allData.DCCurrent_A));
fprintf('  표준편차: %.2f A\n', std(allData.DCCurrent_A, 'omitnan'));

fprintf('\\n=== 전류 통계 (Cell 레벨) ===\n');
fprintf('  평균: %.2f A\n', mean(allData.CellCurrent_A, 'omitnan'));
fprintf('  최소: %.2f A\n', min(allData.CellCurrent_A));
fprintf('  최대: %.2f A\n', max(allData.CellCurrent_A));

fprintf('\\n=== C-rate 통계 (Cell 기준) ===\n');
fprintf('  평균: %.4fC\n', mean(abs(allData.CRate_Cell), 'omitnan'));
fprintf('  최대 충전: %.4fC\n', max(allData.CRate_Cell));
fprintf('  최대 방전: %.4fC\n', abs(min(allData.CRate_Cell)));

fprintf('\\n=== 전압 통계 (Rack) ===\n');
fprintf('  평균: %.2f V\n', mean(allData.SumCV_V, 'omitnan'));
fprintf('  최소: %.2f V\n', min(allData.SumCV_V));
fprintf('  최대: %.2f V\n', max(allData.SumCV_V));

fprintf('\\n=== 전압 통계 (Cell 평균) ===\n');
fprintf('  평균: %.4f V\n', mean(allData.AverageCV_V, 'omitnan'));
fprintf('  최소: %.4f V\n', min(allData.AverageCV_V));
fprintf('  최대: %.4f V\n', max(allData.AverageCV_V));

fprintf('\\n=== 온도 통계 ===\n');
fprintf('  평균: %.2f°C\n', mean(allData.AverageMT_degC, 'omitnan'));
fprintf('  최소: %.2f°C\n', min(allData.AverageMT_degC));
fprintf('  최대: %.2f°C\n', max(allData.AverageMT_degC));

fprintf('\\n=== 전력 통계 ===\n');
fprintf('  평균: %.2f kW\n', mean(allData.DCPower_kW, 'omitnan'));
fprintf('  최소: %.2f kW (방전)\n', min(allData.DCPower_kW));
fprintf('  최대: %.2f kW (충전)\n', max(allData.DCPower_kW));

%% 3단계: 운전 패턴 분석
fprintf('\\n[3단계] 운전 패턴 분석\n');

% 전류가 0이 아닌 구간 (실제 운전 구간)
activeIdx = abs(allData.DCCurrent_A) > 0.1;  % 0.1A 이상
numActive = sum(activeIdx);
numIdle = sum(~activeIdx);

fprintf('\\n=== 운전 상태 ===\n');
fprintf('  활성 구간: %d개 (%.2f%%)\n', numActive, numActive/length(activeIdx)*100);
fprintf('  대기 구간: %d개 (%.2f%%)\n', numIdle, numIdle/length(activeIdx)*100);

% 충전/방전 구분
chargingIdx = allData.DCCurrent_A > 0.1;
dischargingIdx = allData.DCCurrent_A < -0.1;

fprintf('\\n=== 충방전 비율 ===\n');
fprintf('  충전: %d개 (%.2f%%)\n', sum(chargingIdx), sum(chargingIdx)/length(activeIdx)*100);
fprintf('  방전: %d개 (%.2f%%)\n', sum(dischargingIdx), sum(dischargingIdx)/length(activeIdx)*100);

if numActive > 0
    fprintf('\\n=== 활성 구간 전류 통계 ===\n');
    fprintf('  평균 전류: %.2f A\n', mean(allData.DCCurrent_A(activeIdx)));
    fprintf('  평균 Cell 전류: %.2f A\n', mean(allData.CellCurrent_A(activeIdx)));
    fprintf('  평균 C-rate: %.4fC\n', mean(abs(allData.CRate_Cell(activeIdx))));
    
    if sum(chargingIdx) > 0
        fprintf('\\n=== 충전 구간 ===\n');
        fprintf('  평균 충전 전류: %.2f A\n', mean(allData.DCCurrent_A(chargingIdx)));
        fprintf('  최대 충전 전류: %.2f A\n', max(allData.DCCurrent_A(chargingIdx)));
        fprintf('  평균 충전 C-rate: %.4fC\n', mean(allData.CRate_Cell(chargingIdx)));
    end
    
    if sum(dischargingIdx) > 0
        fprintf('\\n=== 방전 구간 ===\n');
        fprintf('  평균 방전 전류: %.2f A\n', abs(mean(allData.DCCurrent_A(dischargingIdx))));
        fprintf('  최대 방전 전류: %.2f A\n', abs(min(allData.DCCurrent_A(dischargingIdx))));
        fprintf('  평균 방전 C-rate: %.4fC\n', abs(mean(allData.CRate_Cell(dischargingIdx))));
    end
end

% 전력량 계산 (1초 샘플링)
chargeEnergy_kWh = sum(allData.DCPower_kW(chargingIdx)) / 3600;
dischargeEnergy_kWh = abs(sum(allData.DCPower_kW(dischargingIdx)) / 3600);

fprintf('\\n=== 에너지 처리량 (6월 전체) ===\n');
fprintf('  충전량: %.2f kWh\n', chargeEnergy_kWh);
fprintf('  방전량: %.2f kWh\n', dischargeEnergy_kWh);
fprintf('  순 에너지: %.2f kWh\n', chargeEnergy_kWh - dischargeEnergy_kWh);
fprintf('  효율: %.2f%%\n', dischargeEnergy_kWh/chargeEnergy_kWh*100);

%% 4단계: 전압-전류 관계 분석 (저항 추정)
fprintf('\\n[4단계] 전압-전류 관계 분석\n');

% 활성 구간에서 전압-전류 상관관계
if numActive > 100
    V_active = allData.AverageCV_V(activeIdx);
    I_active = allData.CellCurrent_A(activeIdx);
    
    % 전압 변화량과 전류의 관계
    corrVI = corr(V_active, I_active);
    fprintf('\\n전압-전류 상관계수: %.4f\n', corrVI);
    
    % 충전/방전별 전압 변화
    if sum(chargingIdx) > 100 && sum(dischargingIdx) > 100
        V_charge = allData.AverageCV_V(chargingIdx);
        V_discharge = allData.AverageCV_V(dischargingIdx);
        I_charge = allData.CellCurrent_A(chargingIdx);
        I_discharge = allData.CellCurrent_A(dischargingIdx);
        
        fprintf('\\n=== 충전 시 평균 전압 ===\n');
        fprintf('  %.4f V\n', mean(V_charge));
        fprintf('\\n=== 방전 시 평균 전압 ===\n');
        fprintf('  %.4f V\n', mean(V_discharge));
        fprintf('\\n평균 전압 차이: %.4f V\n', mean(V_charge) - mean(V_discharge));
        
        % 동일 SOC 구간에서 충방전 전압 차이 (저항 추정용)
        SOC_bins = 40:5:50;  % SOC 40~50% 구간
        fprintf('\\n=== SOC 구간별 전압 차이 (저항 추정) ===\n');
        
        for soc_target = SOC_bins
            soc_range = [soc_target-0.5, soc_target+0.5];
            
            idx_c = chargingIdx & (allData.SOCPct >= soc_range(1)) & (allData.SOCPct <= soc_range(2));
            idx_d = dischargingIdx & (allData.SOCPct >= soc_range(1)) & (allData.SOCPct <= soc_range(2));
            
            if sum(idx_c) > 10 && sum(idx_d) > 10
                V_c_mean = mean(allData.AverageCV_V(idx_c));
                V_d_mean = mean(allData.AverageCV_V(idx_d));
                I_c_mean = mean(allData.CellCurrent_A(idx_c));
                I_d_mean = mean(allData.CellCurrent_A(idx_d));
                
                dV = V_c_mean - V_d_mean;
                dI = I_c_mean - I_d_mean;
                
                if abs(dI) > 1
                    R_est = dV / dI * 1000;  % mOhm
                    fprintf('  SOC %.1f%%: dV=%.4fV, dI=%.2fA, R=%.2f mOhm\n', ...
                            soc_target, dV, dI, R_est);
                end
            end
        end
    end
else
    fprintf('\\n활성 구간이 부족하여 전압-전류 관계 분석 불가\n');
end

%% 5단계: 용량 분석
fprintf('\\n[5단계] 용량 분석\n');

% SOC 범위 확인
SOC_range = max(allData.SOCPct) - min(allData.SOCPct);
fprintf('\\nSOC 변동 범위: %.2f%%\n', SOC_range);

% 완전 충방전 사이클 탐지
fprintf('\\n=== 사이클 분석 ===\n');

% SOC 변화량이 큰 구간 찾기
dSOC = diff(allData.SOCPct);
large_soc_change = find(abs(dSOC) > 5);  % 5% 이상 변화

fprintf('5%% 이상 SOC 변화 구간: %d개\n', length(large_soc_change));

% 충전 사이클 카운팅 (간단한 방법)
SOC_diff = diff(allData.SOCPct);
charge_cycles = 0;
discharge_cycles = 0;
accumulated_charge = 0;
accumulated_discharge = 0;

for i = 1:length(SOC_diff)
    if SOC_diff(i) > 0
        accumulated_charge = accumulated_charge + SOC_diff(i);
    else
        accumulated_discharge = accumulated_discharge + abs(SOC_diff(i));
    end
end

equivalent_cycles = min(accumulated_charge, accumulated_discharge) / 100;

fprintf('누적 충전 SOC: %.2f%%\n', accumulated_charge);
fprintf('누적 방전 SOC: %.2f%%\n', accumulated_discharge);
fprintf('등가 사이클 수: %.2f cycles\n', equivalent_cycles);

% 전류 적산으로 용량 추정 시도
if sum(chargingIdx) > 0
    % 충전 구간 Ah 적산
    chargeAh_cell = sum(allData.CellCurrent_A(chargingIdx)) / 3600;  % 1초 샘플
    dischargeAh_cell = abs(sum(allData.CellCurrent_A(dischargingIdx)) / 3600);
    
    fprintf('\\n=== 전류 적산 (Cell 기준) ===\n');
    fprintf('  충전량: %.2f Ah\n', chargeAh_cell);
    fprintf('  방전량: %.2f Ah\n', dischargeAh_cell);
    
    % SOC 변화와 Ah의 관계
    if accumulated_charge > 0
        capacity_est = chargeAh_cell / (accumulated_charge / 100);
        fprintf('\\n추정 Cell 용량: %.2f Ah\n', capacity_est);
        fprintf('공칭 용량 대비: %.2f%%\n', capacity_est / cellCapacity_Ah * 100);
    end
end

%% 6단계: 결과 저장
fprintf('\\n[6단계] 결과 저장 중...\n');

save(fullfile(savePath, 'June2021_Rack01_Analysis.mat'), 'allData');

% 요약 리포트 작성
reportFile = fullfile(savePath, 'June2021_Rack01_Summary.txt');
fid = fopen(reportFile, 'w');

fprintf(fid, '===============================================\\n');
fprintf(fid, '  Rack01 2021년 6월 분석 결과 요약\\n');
fprintf(fid, '===============================================\\n');
fprintf(fid, '분석 일시: %s\\n\\n', datestr(now));

fprintf(fid, '[데이터 규모]\\n');
fprintf(fid, '  총 데이터: %d개\\n', length(allData.SOCPct));
fprintf(fid, '  기간: 30일\\n\\n');

fprintf(fid, '[Rack 사양]\\n');
fprintf(fid, '  Cell 용량: %d Ah\\n', cellCapacity_Ah);
fprintf(fid, '  Module 용량: %d Ah (2P14S)\\n', moduleCapacity_Ah);
fprintf(fid, '  Rack 용량: %d Ah (17S)\\n\\n', rackCapacity_Ah);

fprintf(fid, '[운전 통계]\\n');
fprintf(fid, '  활성 비율: %.2f%%\\n', numActive/length(activeIdx)*100);
fprintf(fid, '  충전 비율: %.2f%%\\n', sum(chargingIdx)/length(activeIdx)*100);
fprintf(fid, '  방전 비율: %.2f%%\\n', sum(dischargingIdx)/length(activeIdx)*100);
fprintf(fid, '  등가 사이클: %.2f\\n\\n', equivalent_cycles);

fprintf(fid, '[SOC 범위]\\n');
fprintf(fid, '  평균: %.2f%%\\n', mean(allData.SOCPct, 'omitnan'));
fprintf(fid, '  범위: %.2f%% ~ %.2f%%\\n', min(allData.SOCPct), max(allData.SOCPct));
fprintf(fid, '  변동폭: %.2f%%\\n\\n', SOC_range);

fprintf(fid, '[전류 특성 (Cell)]\\n');
fprintf(fid, '  최대 충전: %.2f A (%.4fC)\\n', max(allData.CellCurrent_A), max(allData.CRate_Cell));
fprintf(fid, '  최대 방전: %.2f A (%.4fC)\\n\\n', abs(min(allData.CellCurrent_A)), abs(min(allData.CRate_Cell)));

fprintf(fid, '[에너지 처리량]\\n');
fprintf(fid, '  충전: %.2f kWh\\n', chargeEnergy_kWh);
fprintf(fid, '  방전: %.2f kWh\\n', dischargeEnergy_kWh);
fprintf(fid, '  효율: %.2f%%\\n\\n', dischargeEnergy_kWh/chargeEnergy_kWh*100);

fclose(fid);

fprintf('\\n저장 완료!\\n');
fprintf('  - June2021_Rack01_Analysis.mat\\n');
fprintf('  - June2021_Rack01_Summary.txt\\n');

fprintf('\\n========================================\\n');
fprintf('         분석 완료\\n');
fprintf('========================================\\n');
