%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drive Cycle Data Parser (Modified)
% 8개 채널의 실부하 프로파일 데이터를 3개의 SOC별로 파싱
% 후기 휴지기 이후 방전/충전 데이터 제거
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% 경로 설정 - 현재 디렉토리 사용
dataDir = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data\Drive Cycle';
saveFolder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';

if ~exist(saveFolder, 'dir'); mkdir(saveFolder); end

% 사이클별 파일 목록 (8개 채널) - Auto-detect available cycles
% cycleTypes = {'0cyc', '200cyc', '400cyc','600cyc'};  % 하드코딩 제거

% Auto-detect available cycle types from CSV files
fprintf('사용 가능한 사이클 타입 자동 감지 중...\n');
csvFiles = dir(fullfile(dataDir, 'Ch*_Drive_*cyc.csv'));
cycleTypes = {};

for i = 1:length(csvFiles)
    filename = csvFiles(i).name;
    % Extract cycle type from filename (e.g., '0cyc', '200cyc', '400cyc', '600cyc')
    match = regexp(filename, 'Ch\d+_Drive_(\d+cyc)\.csv', 'tokens');
    if ~isempty(match)
        cycleType = match{1}{1};
        if ~ismember(cycleType, cycleTypes)
            cycleTypes{end+1} = cycleType;
        end
    end
end

% Sort cycle types numerically
cycleNumbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), cycleTypes);
[sortedNumbers, sortIdx] = sort(cycleNumbers);
cycleTypes = cycleTypes(sortIdx);

fprintf('감지된 사이클: %s\n', strjoin(cycleTypes, ', '));

% Generate file names for each detected cycle
cycleFileNames = {};
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    fileNames = {};
    for ch = 9:16  % Channels 9-16
        fileNames{end+1} = sprintf('Ch%d_Drive_%s.csv', ch, cycleType);
    end
    cycleFileNames{cycleIdx} = fileNames;
end

% SOC별 Step Index 정의
SOC90_stepIndex = [5, 7, 9, 11, 13, 15, 17, 19];
SOC70_stepIndex = [23, 25, 27, 29, 31, 33, 35, 37];
SOC50_stepIndex = [41, 43, 45, 47, 49, 51, 53, 55];

% 실부하 프로파일 이름 정의 (8개)
profileNames = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

% 휴지기 시간 설정 (초)
initialRestTime = 5 * 60;  % 초기 8분 휴지기
finalRestTime = 5 * 60;    % 후기 8분 휴지기 (최소값)

% 강제 재처리 옵션 (true로 설정하면 기존 파일이 있어도 재처리)
forceReprocess = true;

fprintf('실부하 프로파일 데이터 파싱 시작...\n');

% 각 사이클별로 처리
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    fileNames = cycleFileNames{cycleIdx};
    
    % 저장 파일 경로 확인
    savePath = fullfile(saveFolder, sprintf('parsedDriveCycle_%s_filtered.mat', cycleType));
    
    % 기존 파일이 있고 강제 재처리가 아니면 스킵
    if exist(savePath, 'file') && ~forceReprocess
        fprintf('\n=== %s 사이클: 이미 처리됨 (파일 존재) - 스킵 ===\n', cycleType);
        fprintf('파일 경로: %s\n', savePath);
        fprintf('재처리를 원하면 forceReprocess = true로 설정하세요.\n');
        continue;
    end
    
    if exist(savePath, 'file') && forceReprocess
        fprintf('\n=== %s 사이클: 강제 재처리 모드 - 기존 파일 덮어쓰기 ===\n', cycleType);
    else
        fprintf('\n=== %s 사이클 처리 시작 (새 파일) ===\n', cycleType);
    end
    
    % 각 사이클별 데이터 구조체 초기화
    eval(sprintf('parsedDriveCycle_%s = struct();', cycleType));
    fprintf('데이터 디렉토리: %s\n', dataDir);
    
    for i = 1:length(fileNames)
    filename = fileNames{i};
    filepath = fullfile(dataDir, filename);
    
    fprintf('파일 확인: %s\n', filepath);
    if exist(filepath, 'file') ~= 2
        fprintf('파일 없음: %s\n', filepath);
        continue;
    end
    
    fprintf('처리 중: %s\n', filename);
    
    % CSV 파일 읽기
    T = readtable(filepath, 'VariableNamingRule', 'preserve');
    
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
    channelFieldName = sprintf('ch%s_Drive_%s', channelName{1}, cycleType);
    
    fprintf('  [디버깅] 채널 파일 로드 완료: %s\n', filename);
    
    % 채널별 구조체 초기화
    eval(sprintf('parsedDriveCycle_%s.(channelFieldName) = struct();', cycleType));
    
    % SOC90 데이터 추출 및 처리
    fprintf('  SOC90 데이터 추출 중...\n');
    fprintf('  [디버깅] SOC별 인덱스 확인 완료: SOC90 Step Indices = %s\n', mat2str(SOC90_stepIndex));
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
            fprintf('    [디버깅] SOC90 %s 휴지구간 탐지\n', profileName);
            [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
                removeFinalRestData(step_voltage, step_current, step_time, step_totalTime, finalRestTime, 'SOC90', profileName);
            
            % 구조체에 저장
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC90.(profileName).V = filtered_V;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC90.(profileName).I = filtered_I;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC90.(profileName).t = filtered_t;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC90.(profileName).totalTime = filtered_totalTime;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC90.(profileName).stepIndex = stepIdx;', cycleType));
            
            fprintf('    %s: %d -> %d data points (후기 휴지기 이후 제거)\n', ...
                    profileName, length(step_voltage), length(filtered_V));
        else
            fprintf('    %s: 데이터 없음 (Step Index %d)\n', profileName, stepIdx);
        end
    end
    
    % SOC70 데이터 추출 및 처리
    fprintf('  SOC70 데이터 추출 중...\n');
    fprintf('  [디버깅] SOC별 인덱스 확인 완료: SOC70 Step Indices = %s\n', mat2str(SOC70_stepIndex));
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
            fprintf('    [디버깅] SOC70 %s 휴지구간 탐지\n', profileName);
            [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
                removeFinalRestData(step_voltage, step_current, step_time, step_totalTime, finalRestTime, 'SOC70', profileName);
            
            % 구조체에 저장
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC70.(profileName).V = filtered_V;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC70.(profileName).I = filtered_I;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC70.(profileName).t = filtered_t;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC70.(profileName).totalTime = filtered_totalTime;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC70.(profileName).stepIndex = stepIdx;', cycleType));
            
            fprintf('    %s: %d -> %d data points (후기 휴지기 이후 제거)\n', ...
                    profileName, length(step_voltage), length(filtered_V));
        else
            fprintf('    %s: 데이터 없음 (Step Index %d)\n', profileName, stepIdx);
        end
    end
    
    % SOC50 데이터 추출 및 처리
    fprintf('  SOC50 데이터 추출 중...\n');
    fprintf('  [디버깅] SOC별 인덱스 확인 완료: SOC50 Step Indices = %s\n', mat2str(SOC50_stepIndex));
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
            fprintf('    [디버깅] SOC50 %s 휴지구간 탐지\n', profileName);
            [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
                removeFinalRestData(step_voltage, step_current, step_time, step_totalTime, finalRestTime, 'SOC50', profileName);
            
            % 구조체에 저장
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC50.(profileName).V = filtered_V;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC50.(profileName).I = filtered_I;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC50.(profileName).t = filtered_t;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC50.(profileName).totalTime = filtered_totalTime;', cycleType));
            eval(sprintf('parsedDriveCycle_%s.(channelFieldName).SOC50.(profileName).stepIndex = stepIdx;', cycleType));
            
            fprintf('    %s: %d -> %d data points (후기 휴지기 이후 제거)\n', ...
                    profileName, length(step_voltage), length(filtered_V));
        else
            fprintf('    %s: 데이터 없음 (Step Index %d)\n', profileName, stepIdx);
        end
    end
    
    fprintf('  %s 완료\n\n', channelName{1});
    end
    
    % 각 사이클별 결과 저장
    savePath = fullfile(saveFolder, sprintf('parsedDriveCycle_%s_filtered.mat', cycleType));
    eval(sprintf('save(savePath, ''parsedDriveCycle_%s'');', cycleType));
    fprintf('%s 사이클 파싱 완료! 결과가 저장되었습니다: %s\n', cycleType, savePath);
end

fprintf('\n=== 전체 파싱 완료 ===\n');

% 각 사이클별 구조체 요약 출력
for cycleIdx = 1:length(cycleTypes)
    cycleType = cycleTypes{cycleIdx};
    fprintf('\n=== %s 사이클 파싱 결과 요약 ===\n', cycleType);
    
    eval(sprintf('currentData = parsedDriveCycle_%s;', cycleType));
    channels = fieldnames(currentData);
    
    for i = 1:length(channels)
        channelName = channels{i};
        fprintf('채널: %s\n', channelName);
        
        % 각 SOC별 데이터 개수 확인
        if isfield(currentData.(channelName), 'SOC90')
            soc90_fields = fieldnames(currentData.(channelName).SOC90);
            fprintf('  SOC90: %d개 프로파일\n', length(soc90_fields));
        end
        
        if isfield(currentData.(channelName), 'SOC70')
            soc70_fields = fieldnames(currentData.(channelName).SOC70);
            fprintf('  SOC70: %d개 프로파일\n', length(soc70_fields));
        end
        
        if isfield(currentData.(channelName), 'SOC50')
            soc50_fields = fieldnames(currentData.(channelName).SOC50);
            fprintf('  SOC50: %d개 프로파일\n', length(soc50_fields));
        end
    end
    
    fprintf('총 파싱된 데이터: %d개 채널 × 3개 SOC × 8개 프로파일 = %d개 데이터셋\n', ...
            length(channels), length(channels) * 3 * 8);
end

%% 서브 함수: 후기 휴지기 이후 데이터 제거
function [filtered_V, filtered_I, filtered_t, filtered_totalTime] = ...
    removeFinalRestData(voltage, current, time, totalTime, finalRestTime, socLevel, profileName)
    
    % time 변수를 숫자(초)로 변환 (duration 타입일 수 있음)
    if isduration(time)
        time_sec = seconds(time);
    elseif isnumeric(time)
        time_sec = time;
    else
        error('time 변수 타입을 확인할 수 없습니다: %s', class(time));
    end
    
    % 시간을 0부터 시작하도록 정규화 (숫자로 변환)
    time_normalized = time_sec - time_sec(1);
    
    % 전류가 0에 가까운 지점 찾기 (휴지기 판단)
    current_threshold = 2;  % 2A 이하를 휴지기로 판단
    rest_mask = abs(current) < current_threshold;
    
    % 연속된 휴지기 구간들 찾기
    rest_periods = [];  % [start_idx, end_idx, duration_sec] 형태로 저장
    
    i = 1;
    while i <= length(rest_mask)
        if rest_mask(i)
            % 휴지기 시작 지점
            start_idx = i;
            % 휴지기 끝 지점 찾기
            while i <= length(rest_mask) && rest_mask(i)
                i = i + 1;
            end
            end_idx = i - 1;
            
            % 휴지기 지속 시간 계산 (초 단위)
            if end_idx > start_idx
                duration_sec = time_normalized(end_idx) - time_normalized(start_idx);
            else
                duration_sec = 0;
            end
            
            rest_periods = [rest_periods; start_idx, end_idx, duration_sec];
        else
            i = i + 1;
        end
    end
    
    % finalRestTime 이상 지속되는 휴지기 구간들 필터링 (숫자 비교)
    long_rest_periods = rest_periods(rest_periods(:, 3) >= finalRestTime, :);
    
    % 디버깅: finalRestTime 이상 휴지기 구간만 출력
    fprintf('      [디버깅] %s %s %.1f분 이상 휴지기 구간 탐지 결과:\n', socLevel, profileName, finalRestTime/60);
    fprintf('        %.1f분 이상 지속되는 휴지기 구간 개수: %d개\n', finalRestTime/60, size(long_rest_periods, 1));
    if size(long_rest_periods, 1) > 0
        for k = 1:size(long_rest_periods, 1)
            % duration_sec는 이미 숫자(초) 단위
            duration_sec = long_rest_periods(k, 3);
            fprintf('          구간 %d: 시작=%d, 끝=%d, 지속시간=%.2f초 (%.2f분)\n', ...
                k, long_rest_periods(k, 1), long_rest_periods(k, 2), duration_sec, duration_sec/60);
        end
    end
    
    % 두번째 휴지기 구간 찾기 (8분 이상 지속되는 휴지기 중 두번째)
    if size(long_rest_periods, 1) >= 2
        % 두번째 휴지기 구간의 끝 지점
        second_rest_end_idx = long_rest_periods(2, 2);
        
        % 두번째 휴지기 구간의 끝 지점 이후 모든 데이터 제거
        filtered_V = voltage(1:second_rest_end_idx);
        filtered_I = current(1:second_rest_end_idx);
        % time을 숫자(초)로 변환하여 반환 (일관성 유지)
        filtered_t = time_sec(1:second_rest_end_idx);
        filtered_totalTime = totalTime(1:second_rest_end_idx);
    elseif size(long_rest_periods, 1) == 1
        % 8분 이상 휴지기가 1개만 있으면 오류
        error('ERROR: 8분 이상 지속되는 휴지기 구간이 1개만 발견되었습니다. 두번째 휴지기 구간이 필요합니다.\n휴지기 구간 정보: 시작=%d, 끝=%d, 지속시간=%.2f초', ...
            long_rest_periods(1, 1), long_rest_periods(1, 2), long_rest_periods(1, 3));
    else
        % 8분 이상 지속되는 휴지기 구간이 없으면 오류
        error('ERROR: 8분 이상 지속되는 휴지기 구간이 없습니다. 데이터 구조를 확인하세요.');
    end
end