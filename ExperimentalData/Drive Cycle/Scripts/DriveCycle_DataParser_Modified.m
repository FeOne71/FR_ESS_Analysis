%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drive Cycle Data Parser (Modified)
% 8개 채널의 실부하 프로파일 데이터를 3개의 SOC별로 파싱
% 후기 휴지기 이후 방전/충전 데이터 제거
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% 경로 설정 - 현재 디렉토리 사용
folderPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\csv';
saveFolder = fullfile(folderPath, 'parsed_data');
if ~exist(saveFolder, 'dir'); mkdir(saveFolder); end

% 파일 목록 (8개 채널)
fileNames = {
    'Ch9_Drive_200cyc.csv';
    'Ch10_Drive_200cyc.csv';
    'Ch11_Drive_200cyc.csv';
    'Ch12_Drive_200cyc.csv';
    'Ch13_Drive_200cyc.csv';
    'Ch14_Drive_200cyc.csv';
    'Ch15_Drive_200cyc.csv';
    'Ch16_Drive_200cyc.csv'
};

% SOC별 Step Index 정의
SOC90_stepIndex = [5, 7, 9, 11, 13, 15, 17, 19];
SOC70_stepIndex = [23, 25, 27, 29, 31, 33, 35, 37];
SOC50_stepIndex = [41, 43, 45, 47, 49, 51, 53, 55];

% 실부하 프로파일 이름 정의 (8개)
profileNames = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

% 휴지기 시간 설정 (초)
initialRestTime = 8 * 60;  % 초기 8분 휴지기
finalRestTime = 8 * 60;    % 후기 8분 휴지기 (최소값)

% 각 채널별 데이터 구조체 초기화
parsedDriveCycle_200cyc = struct();

fprintf('실부하 프로파일 데이터 파싱 시작...\n');

for i = 1:length(fileNames)
    filename = fileNames{i};
    filepath = fullfile(folderPath, filename);
    
    if exist(filepath, 'file') ~= 2
        fprintf('파일 없음: %s\n', filepath);
        continue;
    end
    
    fprintf('처리 중: %s\n', filename);
    
    % CSV 파일 읽기
    T = readtable(filepath);
    
    % 데이터 추출
    stepIndex = T{:,2};     % Step Index
    stepType = T{:,3};      % Step Type
    time = T{:,5};          % Time [s]
    totalTime = T{:,6};     % Total Time [s]
    current = T{:,7};       % Current [A]
    voltage = T{:,8};       % Voltage [V]
    
    % 채널 이름 추출
    [~, baseName, ~] = fileparts(filename);
    channelName = extractBetween(baseName, 'Ch', '_');
    channelFieldName = sprintf('ch%s_Drive_200cyc', channelName{1});
    
    % 채널별 구조체 초기화
    parsedDriveCycle_200cyc.(channelFieldName) = struct();
    
    % SOC90 데이터 추출 및 처리
    fprintf('  SOC90 데이터 추출 중...\n');
    for j = 1:length(SOC90_stepIndex)
        stepIdx = SOC90_stepIndex(j);
        profileName = profileNames{j};
        
        % 해당 Step Index와 Step Type="SIM"인 데이터 찾기
        mask = (stepIndex == stepIdx) & strcmp(stepType, 'SIM');
        
        if any(mask)
            % 원본 데이터 추출
            step_voltage = voltage(mask);
            step_current = current(mask);
            step_time = time(mask);
            step_totalTime = totalTime(mask);
            
            % 후기 휴지기 이후 데이터 제거
            [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
                removeFinalRestData(step_voltage, step_current, step_time, step_totalTime, finalRestTime);
            
            % 구조체에 저장
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).V = filtered_V;
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).I = filtered_I;
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).t = filtered_t;
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).totalTime = filtered_totalTime;
            parsedDriveCycle_200cyc.(channelFieldName).SOC90.(profileName).stepIndex = stepIdx;
            
            fprintf('    %s: %d -> %d data points (후기 휴지기 이후 제거)\n', ...
                    profileName, length(step_voltage), length(filtered_V));
        else
            fprintf('    %s: 데이터 없음 (Step Index %d)\n', profileName, stepIdx);
        end
    end
    
    % SOC70 데이터 추출 및 처리
    fprintf('  SOC70 데이터 추출 중...\n');
    for j = 1:length(SOC70_stepIndex)
        stepIdx = SOC70_stepIndex(j);
        profileName = profileNames{j};
        
        % 해당 Step Index와 Step Type="SIM"인 데이터 찾기
        mask = (stepIndex == stepIdx) & strcmp(stepType, 'SIM');
        
        if any(mask)
            % 원본 데이터 추출
            step_voltage = voltage(mask);
            step_current = current(mask);
            step_time = time(mask);
            step_totalTime = totalTime(mask);
            
            % 후기 휴지기 이후 데이터 제거
            [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
                removeFinalRestData(step_voltage, step_current, step_time, step_totalTime, finalRestTime);
            
            % 구조체에 저장
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).V = filtered_V;
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).I = filtered_I;
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).t = filtered_t;
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).totalTime = filtered_totalTime;
            parsedDriveCycle_200cyc.(channelFieldName).SOC70.(profileName).stepIndex = stepIdx;
            
            fprintf('    %s: %d -> %d data points (후기 휴지기 이후 제거)\n', ...
                    profileName, length(step_voltage), length(filtered_V));
        else
            fprintf('    %s: 데이터 없음 (Step Index %d)\n', profileName, stepIdx);
        end
    end
    
    % SOC50 데이터 추출 및 처리
    fprintf('  SOC50 데이터 추출 중...\n');
    for j = 1:length(SOC50_stepIndex)
        stepIdx = SOC50_stepIndex(j);
        profileName = profileNames{j};
        
        % 해당 Step Index와 Step Type="SIM"인 데이터 찾기
        mask = (stepIndex == stepIdx) & strcmp(stepType, 'SIM');
        
        if any(mask)
            % 원본 데이터 추출
            step_voltage = voltage(mask);
            step_current = current(mask);
            step_time = time(mask);
            step_totalTime = totalTime(mask);
            
            % 후기 휴지기 이후 데이터 제거
            [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
                removeFinalRestData(step_voltage, step_current, step_time, step_totalTime, finalRestTime);
            
            % 구조체에 저장
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).V = filtered_V;
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).I = filtered_I;
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).t = filtered_t;
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).totalTime = filtered_totalTime;
            parsedDriveCycle_200cyc.(channelFieldName).SOC50.(profileName).stepIndex = stepIdx;
            
            fprintf('    %s: %d -> %d data points (후기 휴지기 이후 제거)\n', ...
                    profileName, length(step_voltage), length(filtered_V));
        else
            fprintf('    %s: 데이터 없음 (Step Index %d)\n', profileName, stepIdx);
        end
    end
    
    fprintf('  %s 완료\n\n', channelName{1});
end

% 결과 저장
savePath = fullfile(saveFolder, 'parsedDriveCycle_200cyc_filtered.mat');
save(savePath, 'parsedDriveCycle_200cyc');

fprintf('파싱 완료! 결과가 저장되었습니다: %s\n', savePath);

% 구조체 요약 출력
fprintf('\n=== 파싱 결과 요약 ===\n');
channels = fieldnames(parsedDriveCycle_200cyc);
for i = 1:length(channels)
    channelName = channels{i};
    fprintf('채널: %s\n', channelName);
    
    % 각 SOC별 데이터 개수 확인
    if isfield(parsedDriveCycle_200cyc.(channelName), 'SOC90')
        soc90_fields = fieldnames(parsedDriveCycle_200cyc.(channelName).SOC90);
        fprintf('  SOC90: %d개 프로파일\n', length(soc90_fields));
    end
    
    if isfield(parsedDriveCycle_200cyc.(channelName), 'SOC70')
        soc70_fields = fieldnames(parsedDriveCycle_200cyc.(channelName).SOC70);
        fprintf('  SOC70: %d개 프로파일\n', length(soc70_fields));
    end
    
    if isfield(parsedDriveCycle_200cyc.(channelName), 'SOC50')
        soc50_fields = fieldnames(parsedDriveCycle_200cyc.(channelName).SOC50);
        fprintf('  SOC50: %d개 프로파일\n', length(soc50_fields));
    end
end

fprintf('\n총 파싱된 데이터: %d개 채널 × 3개 SOC × 8개 프로파일 = %d개 데이터셋\n', ...
        length(channels), length(channels) * 3 * 8);

%% 서브 함수: 후기 휴지기 이후 데이터 제거
function [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
    removeFinalRestData(voltage, current, time, totalTime, finalRestTime)
    
    % 시간을 0부터 시작하도록 정규화
    time_normalized = time - time(1);
    
    % 전류가 0에 가까운 지점 찾기 (휴지기 판단)
    current_threshold = 2;  % 1mA 이하를 휴지기로 판단
    rest_mask = abs(current) < current_threshold;
    
    % 시간 순서대로 정렬되어 있다고 가정
    % 후기 휴지기 시작점 찾기
    finalRestStartIdx = [];
    
    % 실부하 프로파일이 끝난 후 지속적인 휴지기 구간 찾기
    for i = length(rest_mask):-1:1
        if ~rest_mask(i)
            % 마지막 비휴지기 지점 찾음
            % 이후 휴지기가 finalRestTime 이상 지속되는지 확인
            if i < length(rest_mask)
                potentialRestStart = i + 1;
                restDuration = time_normalized(end) - time_normalized(potentialRestStart);
                
                if restDuration >= finalRestTime
                    finalRestStartIdx = potentialRestStart;
                    break;
                end
            end
        end
    end
    
    % 후기 휴지기 시작점이 발견되면 해당 지점 이후 데이터 제거
    if ~isempty(finalRestStartIdx)
        % 후기 휴지기 시작점에서 finalRestTime만큼 유지한 후 제거
        cutoffTime = time_normalized(finalRestStartIdx) + finalRestTime;
        cutoffIdx = find(time_normalized <= cutoffTime, 1, 'last');
        
        if ~isempty(cutoffIdx)
            filtered_V = voltage(1:cutoffIdx);
            filtered_I = current(1:cutoffIdx);
            filtered_t = time(1:cutoffIdx);
            filtered_totalTime = totalTime(1:cutoffIdx);
        else
            % cutoffIdx가 없으면 원본 데이터 유지
            filtered_V = voltage;
            filtered_I = current;
            filtered_t = time;
            filtered_totalTime = totalTime;
        end
    else
        % 후기 휴지기를 찾지 못한 경우 원본 데이터 유지
        filtered_V = voltage;
        filtered_I = current;
        filtered_t = time;
        filtered_totalTime = totalTime;
    end
end 