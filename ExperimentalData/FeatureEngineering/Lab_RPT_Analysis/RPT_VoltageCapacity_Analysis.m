%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_VoltageCapacity_Analysis.m
% 목적: RPT 데이터에서 전압 구간별 용량 추출 및 dQdV 피크 분석
%
% 시각화:
%   1. 전압 구간별 용량 추출 (V-Q 곡선)
%   2. dQdV 피크 높이 및 면적 분석
%
% 데이터 소스:
%   - RPT 파싱된 파일: RPT\Postprocessing\Parsed\RPT###_ch##_parsed.mat
%   - pdata 구조체의 Step별 V, Q, dQdV_AhV 데이터 사용
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Lab RPT Analysis: Voltage-Capacity & dQdV Peak Analysis ===\n');

%% 설정
rptParsedDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';
outputDir = fullfile(pwd, 'Lab_RPT_Analysis', 'Results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% RPT 사이클 및 채널
rptCycles = [0, 200, 400, 600, 800, 1000];
channels = 9;  % ch09만 처리

% 전압 구간 설정 (용량 추출용)
voltage_ranges = [
    3.0, 3.2;   % Low voltage range
    3.2, 3.4;   % Mid-low voltage range
    3.4, 3.6;   % Mid voltage range
    3.6, 3.8;   % Mid-high voltage range
    3.8, 4.0;   % High voltage range
    4.0, 4.2    % Very high voltage range
];

% dQdV 피크 검출 파라미터
% 주의: min_peak_height는 데이터에 따라 자동 조정됨
min_peak_height_charge = 0.1;   % Ah/V (충전) - 기본값, 자동 조정됨
min_peak_height_discharge = 0.1; % Ah/V (방전) - 기본값, 자동 조정됨
min_peak_distance = 0.05;        % V (피크 간 최소 거리)
peak_height_ratio = 0.1;         % 데이터 최대값의 10%를 최소 피크 높이로 사용

% 허용 C-rate 값 (평가 대상)
target_crates = [0.1, 0.5, 1, 2, 3];
crate_tolerance = 0.02;  % 허용 오차 (2%)

% C-rate 계산용 파라미터
I_1C = 64;  % 1C 전류 (A) - 표준 용량 기준

% 시각화에 사용할 Step 타입 선택
% 'Charge': 충전 데이터만, 'Discharge': 방전 데이터만, 'Both': 둘 다
use_step_type = 'Discharge';  % 'Charge', 'Discharge', 또는 'Both'

%% 데이터 수집
fprintf('\n=== Loading RPT Data ===\n');

% 전압 구간별 용량 저장
voltage_capacity_data = struct();

% dQdV 피크 데이터 저장
dQdV_peak_data = struct();

for chIdx = 1:length(channels)
    chNum = channels(chIdx);
    chName = sprintf('ch%02d', chNum);
    
    voltage_capacity_data.(chName) = struct();
    dQdV_peak_data.(chName) = struct();
    
    for cycIdx = 1:length(rptCycles)
        cycleNum = rptCycles(cycIdx);
        cycleName = sprintf('rpt%d', cycleNum);
        
        filename = sprintf('RPT%d_ch%02d_parsed.mat', cycleNum, chNum);
        filepath = fullfile(rptParsedDir, filename);
        
        if ~exist(filepath, 'file')
            fprintf('  Warning: File not found: %s\n', filename);
            continue;
        end
        
        fprintf('  Loading: %s\n', filename);
        
        try
            data_loaded = load(filepath);
            if ~isfield(data_loaded, 'pdata')
                fprintf('    Warning: pdata not found\n');
                continue;
            end
            
            pdata = data_loaded.pdata;
            
            % 디버깅: pdata 구조 확인
            fprintf('    pdata loaded: %d steps\n', length(pdata));
            
            % pdata에서 C-rate 정보 확인 (디버깅)
            valid_steps = 0;
            for s = 1:length(pdata)
                if isfield(pdata(s), 'Crate') && ~isempty(pdata(s).Crate)
                    if isnumeric(pdata(s).Crate)
                        cr = mean(abs(pdata(s).Crate), 'omitnan');
                        if ~isnan(cr) && cr > 0
                            valid_steps = valid_steps + 1;
                        end
                    end
                end
            end
            fprintf('    Valid steps with C-rate: %d/%d\n', valid_steps, length(pdata));
            
            %% 1. 전압 구간별 용량 추출 (Multi-rate CC 데이터 사용)
            % 모든 충전/방전 Step에서 V-Q 데이터 수집
            voltage_capacity_data.(chName).(cycleName) = struct();
            
            for stepIdx = 1:length(pdata)
                step_data = pdata(stepIdx);
                
                % 충전 또는 방전 Step만 사용 (휴지 제외)
                if isfield(step_data, 'type')
                    type_val = step_data.type;
                    if ischar(type_val) || isstring(type_val)
                        stepType = char(type_val(1));  % 첫 번째 문자
                    elseif iscell(type_val) && ~isempty(type_val)
                        stepType = char(type_val{1}(1));
                    elseif isnumeric(type_val) && ~isempty(type_val)
                        % 숫자로 저장된 경우 (예: 67 = 'C', 68 = 'D')
                        stepType = char(type_val(1));
                    else
                        continue;
                    end
                else
                    continue;
                end
                
                % StepType이 'C' 또는 'D'인지 확인
                if ~ischar(stepType) || (stepType ~= 'C' && stepType ~= 'D')
                    continue;
                end
                
                % Step 타입 필터링 (시각화 옵션에 따라)
                if strcmp(use_step_type, 'Charge') && stepType ~= 'C'
                    continue;  % 충전만 사용
                elseif strcmp(use_step_type, 'Discharge') && stepType ~= 'D'
                    continue;  % 방전만 사용
                end
                % 'Both'인 경우는 필터링하지 않음
                
                % C-rate 확인 (pdata에 이미 계산되어 있을 수 있음)
                Crate = NaN;
                if isfield(step_data, 'Crate') && ~isempty(step_data.Crate)
                    if isnumeric(step_data.Crate)
                        % C-rate가 배열인 경우 평균값 사용
                        Crate = mean(abs(step_data.Crate), 'omitnan');  % 절댓값 평균 C-rate
                    end
                end
                
                % C-rate가 없거나 0이면 I에서 계산
                if isnan(Crate) || Crate == 0
                    if isfield(step_data, 'I') && ~isempty(step_data.I)
                        I_mean = mean(abs(step_data.I), 'omitnan');
                        if ~isnan(I_mean) && I_mean > 0
                            Crate = I_mean / I_1C;
                        else
                            continue;  % C-rate를 계산할 수 없으면 스킵
                        end
                    else
                        continue;
                    end
                end
                
                % C-rate 값이 허용 범위 내인지 확인 (0.1, 0.5, 1, 2, 3에 가까운 값만)
                if isnan(Crate)
                    continue;
                end
                
                % target_crates 중 가장 가까운 값 찾기
                [~, closest_idx] = min(abs(target_crates - Crate));
                closest_crate = target_crates(closest_idx);
                crate_diff = abs(Crate - closest_crate);
                
                % 허용 오차 내에 있는지 확인
                if crate_diff > crate_tolerance * closest_crate
                    continue;  % 허용 범위를 벗어남
                end
                
                % 가장 가까운 C-rate로 정규화 (표시용)
                Crate_normalized = closest_crate;
                
                V = step_data.V;
                Q = step_data.Q;
                
                % 유효한 데이터만 사용
                validIdx = ~isnan(V) & ~isnan(Q) & isfinite(V) & isfinite(Q);
                V = V(validIdx);
                Q = Q(validIdx);
                
                if length(V) < 10  % 최소 데이터 포인트 필요
                    continue;
                end
                
                % V와 Q를 정렬 (전압 순서대로)
                if stepType == 'C'  % 충전: 전압 증가
                    [V_sorted, sortIdx] = sort(V);
                    Q_sorted = Q(sortIdx);
                else  % 방전: 전압 감소
                    [V_sorted, sortIdx] = sort(V, 'descend');
                    Q_sorted = Q(sortIdx);
                end
                
                % 중복 제거 (같은 전압 값)
                [V_unique, uniqueIdx] = unique(V_sorted, 'stable');
                Q_unique = Q_sorted(uniqueIdx);
                
                if length(V_unique) < 10
                    continue;
                end
                
                % 전압 구간별 용량 계산 (ΔQ = Q(V2) - Q(V1))
                capacity_by_range = zeros(size(voltage_ranges, 1), 1);
                for rangeIdx = 1:size(voltage_ranges, 1)
                    v_min = voltage_ranges(rangeIdx, 1);
                    v_max = voltage_ranges(rangeIdx, 2);
                    
                    % 전압 구간 내 데이터 찾기
                    idx_in_range = V_unique >= v_min & V_unique <= v_max;
                    
                    if sum(idx_in_range) > 1
                        V_range = V_unique(idx_in_range);
                        Q_range = Q_unique(idx_in_range);
                        
                        % 전압 구간의 시작과 끝에서 용량 차이
                        % V_min과 V_max에 가장 가까운 점 찾기
                        [~, idx_min] = min(abs(V_range - v_min));
                        [~, idx_max] = min(abs(V_range - v_max));
                        
                        if idx_min ~= idx_max
                            Q_at_vmin = Q_range(idx_min);
                            Q_at_vmax = Q_range(idx_max);
                            
                            % ΔQ = |Q(V_max) - Q(V_min)|
                            capacity_by_range(rangeIdx) = abs(Q_at_vmax - Q_at_vmin);
                        end
                    end
                end
                
                % Step 정보 저장 (필드명을 짧게)
                % 필드명: s{StepIdx}{Type}C{Crate_str} (예: s3DC1p0)
                % 정규화된 C-rate 사용 (0.1, 0.5, 1, 2, 3)
                Crate_str = sprintf('%.1f', Crate_normalized);
                Crate_str = strrep(Crate_str, '.', 'p');  % 소수점을 p로 변경
                stepName = sprintf('s%d%cC%s', step_data.StepIdx, stepType, Crate_str);
                
                % MATLAB 필드명 유효성 검사 및 정리 (최대 63자)
                stepName = matlab.lang.makeValidName(stepName);
                if length(stepName) > 63
                    stepName = stepName(1:63);
                end
                
                voltage_capacity_data.(chName).(cycleName).(stepName).V = V_unique;
                voltage_capacity_data.(chName).(cycleName).(stepName).Q = Q_unique;
                voltage_capacity_data.(chName).(cycleName).(stepName).capacity_by_range = capacity_by_range;
                voltage_capacity_data.(chName).(cycleName).(stepName).Crate = Crate;  % 실제 측정값
                voltage_capacity_data.(chName).(cycleName).(stepName).Crate_normalized = Crate_normalized;  % 정규화된 값
                voltage_capacity_data.(chName).(cycleName).(stepName).StepType = stepType;
                voltage_capacity_data.(chName).(cycleName).(stepName).StepIdx = step_data.StepIdx;
            end
            
            %% 2. dQ/dV 피크 분석 (Multi-rate CC 데이터 사용)
            % 모든 충전/방전 Step에서 dQ/dV 계산 및 피크 찾기
            for stepIdx = 1:length(pdata)
                step_data = pdata(stepIdx);
                
                % 충전 또는 방전 Step만 사용 (휴지 제외)
                if isfield(step_data, 'type')
                    type_val = step_data.type;
                    if ischar(type_val) || isstring(type_val)
                        stepType = char(type_val(1));  % 첫 번째 문자
                    elseif iscell(type_val) && ~isempty(type_val)
                        stepType = char(type_val{1}(1));
                    elseif isnumeric(type_val) && ~isempty(type_val)
                        % 숫자로 저장된 경우 (예: 67 = 'C', 68 = 'D')
                        stepType = char(type_val(1));
                    else
                        continue;
                    end
                else
                    continue;
                end
                
                % StepType이 'C' 또는 'D'인지 확인
                if ~ischar(stepType) || (stepType ~= 'C' && stepType ~= 'D')
                    continue;
                end
                
                % Step 타입 필터링 (시각화 옵션에 따라)
                if strcmp(use_step_type, 'Charge') && stepType ~= 'C'
                    continue;  % 충전만 사용
                elseif strcmp(use_step_type, 'Discharge') && stepType ~= 'D'
                    continue;  % 방전만 사용
                end
                % 'Both'인 경우는 필터링하지 않음
                
                % C-rate 확인 (pdata에 이미 계산되어 있을 수 있음)
                Crate = NaN;
                if isfield(step_data, 'Crate') && ~isempty(step_data.Crate)
                    if isnumeric(step_data.Crate)
                        % C-rate가 배열인 경우 평균값 사용
                        Crate = mean(abs(step_data.Crate), 'omitnan');  % 절댓값 평균 C-rate
                    end
                end
                
                % C-rate가 없거나 0이면 I에서 계산
                if isnan(Crate) || Crate == 0
                    if isfield(step_data, 'I') && ~isempty(step_data.I)
                        I_mean = mean(abs(step_data.I), 'omitnan');
                        if ~isnan(I_mean) && I_mean > 0
                            Crate = I_mean / I_1C;
                        else
                            continue;  % C-rate를 계산할 수 없으면 스킵
                        end
                    else
                        continue;
                    end
                end
                
                % C-rate 값이 허용 범위 내인지 확인 (0.1, 0.5, 1, 2, 3에 가까운 값만)
                if isnan(Crate)
                    continue;
                end
                
                % target_crates 중 가장 가까운 값 찾기
                [~, closest_idx] = min(abs(target_crates - Crate));
                closest_crate = target_crates(closest_idx);
                crate_diff = abs(Crate - closest_crate);
                
                % 허용 오차 내에 있는지 확인
                if crate_diff > crate_tolerance * closest_crate
                    continue;  % 허용 범위를 벗어남
                end
                
                % 가장 가까운 C-rate로 정규화 (표시용)
                Crate_normalized = closest_crate;
                
                V = step_data.V;
                Q = step_data.Q;
                
                % 유효한 데이터만 사용
                validIdx = ~isnan(V) & ~isnan(Q) & isfinite(V) & isfinite(Q);
                V = V(validIdx);
                Q = Q(validIdx);
                
                if length(V) < 10
                    continue;
                end
                
                % V와 Q를 정렬
                if stepType == 'C'
                    [V_sorted, sortIdx] = sort(V);
                    Q_sorted = Q(sortIdx);
                else
                    [V_sorted, sortIdx] = sort(V, 'descend');
                    Q_sorted = Q(sortIdx);
                end
                
                % 중복 제거 및 보간을 위한 그리드 생성
                [V_unique, uniqueIdx] = unique(V_sorted, 'stable');
                Q_unique = Q_sorted(uniqueIdx);
                
                if length(V_unique) < 10
                    continue;
                end
                
                % dQ/dV 계산: V에 대한 Q의 미분
                % 부드러운 곡선을 위해 스무딩 적용
                dV_grid = 0.005;  % 5mV 간격
                V_grid = min(V_unique):dV_grid:max(V_unique);
                
                % 보간
                Q_interp = interp1(V_unique, Q_unique, V_grid, 'linear', 'extrap');
                
                % 스무딩 (노이즈 제거)
                Q_smooth = smoothdata(Q_interp, 'gaussian', 15);
                
                % dQ/dV 계산
                dQ = diff(Q_smooth);
                dV = dV_grid;
                dQdV = dQ ./ dV;  % Ah/V
                
                % 중간 전압 값
                V_mid = (V_grid(1:end-1) + V_grid(2:end)) / 2;
                
                % Step 타입에 따라 다른 임계값 사용
                if stepType == 'C'
                    min_height_base = min_peak_height_charge;
                else
                    min_height_base = min_peak_height_discharge;
                end
                
                % 데이터의 최대값에 비례하여 min_height 자동 조정
                max_dQdV = max(abs(dQdV));
                min_height_auto = max(min_height_base, max_dQdV * peak_height_ratio);
                
                % 디버깅: dQ/dV 통계 정보
                if stepIdx == 1 || mod(stepIdx, 20) == 0
                    fprintf('      Step %d: max_dQdV=%.4f, min_height_auto=%.4f\n', ...
                        step_data.StepIdx, max_dQdV, min_height_auto);
                end
                
                % 피크 찾기 (전압 구간 정보 전달)
                [peak_V, peak_dQdV, peak_area, peak_voltage_ranges] = find_dQdV_peaks(V_mid, dQdV, min_height_auto, min_peak_distance, voltage_ranges);
                
                % C-rate 필터링은 위에서 이미 완료됨
                % 필드명 생성 (짧게, 정규화된 C-rate 사용)
                Crate_str = sprintf('%.1f', Crate_normalized);
                Crate_str = strrep(Crate_str, '.', 'p');
                stepName = sprintf('s%d%cC%s', step_data.StepIdx, stepType, Crate_str);
                stepName = matlab.lang.makeValidName(stepName);
                if length(stepName) > 63
                    stepName = stepName(1:63);
                end
                
                if ~isfield(dQdV_peak_data.(chName), cycleName)
                    dQdV_peak_data.(chName).(cycleName) = struct();
                end
                
                dQdV_peak_data.(chName).(cycleName).(stepName).V = V_mid;
                dQdV_peak_data.(chName).(cycleName).(stepName).dQdV = dQdV;
                dQdV_peak_data.(chName).(cycleName).(stepName).peak_V = peak_V;
                dQdV_peak_data.(chName).(cycleName).(stepName).peak_dQdV = peak_dQdV;
                dQdV_peak_data.(chName).(cycleName).(stepName).peak_area = peak_area;
                dQdV_peak_data.(chName).(cycleName).(stepName).peak_voltage_ranges = peak_voltage_ranges;  % 각 피크의 전압 구간 정보
                
                % 피크가 없어도 전체 dQ/dV 곡선의 면적 계산
                if isempty(peak_V) && length(V_mid) > 1 && length(dQdV) > 1
                    % 전체 면적: ∫(dQ/dV) dV = ΔQ
                    total_dQdV_area = trapz(V_mid, abs(dQdV));
                    dQdV_peak_data.(chName).(cycleName).(stepName).total_dQdV_area = total_dQdV_area;
                end
                dQdV_peak_data.(chName).(cycleName).(stepName).Crate = Crate;  % 실제 측정값
                dQdV_peak_data.(chName).(cycleName).(stepName).Crate_normalized = Crate_normalized;  % 정규화된 값
                dQdV_peak_data.(chName).(cycleName).(stepName).StepType = stepType;
                dQdV_peak_data.(chName).(cycleName).(stepName).StepIdx = step_data.StepIdx;
                
                % 디버깅: 저장된 피크 정보 출력
                if ~isempty(peak_V)
                    total_area = sum(peak_area);
                    fprintf('      Step %d (%s): C-rate=%.3f, Peaks=%d, Total Peak Area=%.4f Ah\n', ...
                        step_data.StepIdx, stepType, Crate_normalized, length(peak_V), total_area);
                elseif isfield(dQdV_peak_data.(chName).(cycleName).(stepName), 'total_dQdV_area')
                    total_area = dQdV_peak_data.(chName).(cycleName).(stepName).total_dQdV_area;
                    fprintf('      Step %d (%s): C-rate=%.3f, No peaks found, Total dQ/dV Area=%.4f Ah\n', ...
                        step_data.StepIdx, stepType, Crate_normalized, total_area);
                end
            end
            
        catch ME
            fprintf('    Error loading %s: %s\n', filename, ME.message);
            fprintf('    Error stack: %s\n', getReport(ME));
        end
    end
end

% 데이터 수집 완료 후 요약 출력
fprintf('\n=== Data Collection Summary ===\n');
for chIdx = 1:length(channels)
    chNum = channels(chIdx);
    chName = sprintf('ch%02d', chNum);
    
    fprintf('  %s:\n', chName);
    for cycIdx = 1:length(rptCycles)
        cycleNum = rptCycles(cycIdx);
        cycleName = sprintf('rpt%d', cycleNum);
        
        if isfield(dQdV_peak_data.(chName), cycleName)
            data = dQdV_peak_data.(chName).(cycleName);
            stepNames = fieldnames(data);
            peak_count = 0;
            for s = 1:length(stepNames)
                stepData = data.(stepNames{s});
                if isfield(stepData, 'peak_V') && ~isempty(stepData.peak_V)
                    peak_count = peak_count + length(stepData.peak_V);
                end
            end
            fprintf('    %s: %d steps, %d total peaks\n', cycleName, length(stepNames), peak_count);
        else
            fprintf('    %s: No data\n', cycleName);
        end
    end
end

%% 시각화 생성
fprintf('\n=== Creating Visualizations ===\n');

%% 1. 전압 구간별 용량 추출 시각화
fprintf('  Creating voltage range capacity plots...\n');

fig1 = figure('Position', [100, 100, 2000, 1200], 'Color', 'white');
set(gcf, 'Visible', 'off');

% 서브플롯: C-rate별로 V-Q 곡선 및 전압 구간 표시
numCrates = length(target_crates);
cols = 3;
rows = 2;

chNum = channels(1);  % ch09만 처리
chName = sprintf('ch%02d', chNum);
cycleColors = lines(length(rptCycles));

for crateIdx = 1:numCrates
    targetCrate = target_crates(crateIdx);
    
    subplot(rows, cols, crateIdx);
    hold on;
    
    % 각 사이클별 V-Q 곡선 플롯 (해당 C-rate만)
    for cycIdx = 1:length(rptCycles)
        cycleNum = rptCycles(cycIdx);
        cycleName = sprintf('rpt%d', cycleNum);
        
        if isfield(voltage_capacity_data.(chName), cycleName)
            data = voltage_capacity_data.(chName).(cycleName);
            stepNames = fieldnames(data);
            
            for stepIdx = 1:length(stepNames)
                stepName = stepNames{stepIdx};
                if strcmp(stepName, 'V') || strcmp(stepName, 'Q') || strcmp(stepName, 'capacity_by_range')
                    continue;  % Skip old format fields
                end
                
                stepData = data.(stepName);
                
                % StepType 필터링
                if isfield(stepData, 'StepType')
                    stepType = stepData.StepType;
                    if strcmp(use_step_type, 'Charge') && stepType ~= 'C'
                        continue;  % 충전만 사용
                    elseif strcmp(use_step_type, 'Discharge') && stepType ~= 'D'
                        continue;  % 방전만 사용
                    end
                end
                
                % C-rate 필터링
                if isfield(stepData, 'Crate_normalized')
                    if abs(stepData.Crate_normalized - targetCrate) > 0.01
                        continue;  % 해당 C-rate가 아니면 스킵
                    end
                end
                
                if isfield(stepData, 'V') && isfield(stepData, 'Q')
                    V = stepData.V;
                    Q = stepData.Q;
                    
                    % 샘플링 (너무 많은 포인트 방지)
                    if length(V) > 3000
                        idx = round(linspace(1, length(V), 3000));
                        V = V(idx);
                        Q = Q(idx);
                    end
                    
                    plot(V, Q, '-', 'LineWidth', 1.5, 'Color', cycleColors(cycIdx, :));
                end
            end
        end
    end
    
    % 전압 구간 표시
    y_range = ylim;
    for rangeIdx = 1:size(voltage_ranges, 1)
        v_min = voltage_ranges(rangeIdx, 1);
        v_max = voltage_ranges(rangeIdx, 2);
        fill([v_min, v_max, v_max, v_min], ...
            [y_range(1), y_range(1), y_range(2), y_range(2)], ...
            [0.9 0.9 0.9], 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'LineStyle', '--');
    end
    
    xlabel('Voltage [V]', 'FontSize', 10);
    ylabel('Capacity [Ah]', 'FontSize', 10);
    title(sprintf('%s: V-Q Curves (C-rate = %.1fC)', chName, targetCrate), 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    hold off;
end

sgtitle(sprintf('Voltage Range Capacity Extraction - %s (C-rate separated)', chName), 'FontSize', 14, 'FontWeight', 'bold');

savefig(fig1, fullfile(outputDir, 'VoltageRange_Capacity.fig'));
close(fig1);
fprintf('    Saved: VoltageRange_Capacity.fig\n');

%% 1.5. V(Q) 곡선에 전압 구간 표시 (C-rate별 서브플롯)
fprintf('  Creating V(Q) curves with voltage ranges...\n');

% 표시할 C-rate 필터링 (0.1, 0.5, 2, 3만)
display_crates = [0.1, 0.5, 2, 3];
numDisplayCrates = length(display_crates);

fig1_5 = figure('Position', [100, 100, 2000, 1200], 'Color', 'white');
set(gcf, 'Visible', 'off');

cols = 2;
rows = 2;
cycleColors = lines(length(rptCycles));

chNum = channels(1);  % ch09만 처리
chName = sprintf('ch%02d', chNum);

for crateIdx = 1:numDisplayCrates
    targetCrate = display_crates(crateIdx);
    
    subplot(rows, cols, crateIdx);
    hold on;
    
    % V-Q 곡선 그리기
    for cycIdx = 1:length(rptCycles)
        cycleNum = rptCycles(cycIdx);
        cycleName = sprintf('rpt%d', cycleNum);
        
        if isfield(voltage_capacity_data.(chName), cycleName)
            data = voltage_capacity_data.(chName).(cycleName);
            stepNames = fieldnames(data);
            
            for stepIdx = 1:length(stepNames)
                stepName = stepNames{stepIdx};
                if strcmp(stepName, 'V') || strcmp(stepName, 'Q')
                    continue;
                end
                
                if isfield(data, stepName)
                    stepData = data.(stepName);
                    
                    % StepType 필터링
                    if isfield(stepData, 'StepType')
                        stepType = stepData.StepType;
                        if strcmp(use_step_type, 'Charge') && stepType ~= 'C'
                            continue;  % 충전만 사용
                        elseif strcmp(use_step_type, 'Discharge') && stepType ~= 'D'
                            continue;  % 방전만 사용
                        end
                    end
                    
                    % C-rate 필터링
                    if isfield(stepData, 'Crate_normalized')
                        if abs(stepData.Crate_normalized - targetCrate) > crate_tolerance
                            continue;  % 해당 C-rate가 아니면 스킵
                        end
                    end
                    
                    if isfield(stepData, 'V') && isfield(stepData, 'Q')
                        V = stepData.V;
                        Q = stepData.Q;
                        
                        % 샘플링 (너무 많은 포인트 방지)
                        if length(V) > 3000
                            idx = round(linspace(1, length(V), 3000));
                            V = V(idx);
                            Q = Q(idx);
                        end
                        
                        plot(V, Q, '-', 'LineWidth', 1.5, 'Color', cycleColors(cycIdx, :));
                    end
                end
            end
        end
    end
    
    % 전압 구간을 점선으로 표시
    y_range = ylim;
    for rangeIdx = 1:size(voltage_ranges, 1)
        v_min = voltage_ranges(rangeIdx, 1);
        v_max = voltage_ranges(rangeIdx, 2);
        
        % 수직 점선으로 전압 구간 경계 표시
        plot([v_min, v_min], y_range, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
        plot([v_max, v_max], y_range, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    end
    
    xlabel('Voltage [V]', 'FontSize', 11);
    ylabel('Capacity [Ah]', 'FontSize', 11);
    title(sprintf('V(Q) Curve with Voltage Ranges (C-rate = %.1fC)', targetCrate), 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    hold off;
end

sgtitle(sprintf('V(Q) Curves with Voltage Range Markers - %s', chName), 'FontSize', 14, 'FontWeight', 'bold');

savefig(fig1_5, fullfile(outputDir, 'VQ_Curves_VoltageRanges.fig'));
close(fig1_5);
fprintf('    Saved: VQ_Curves_VoltageRanges.fig\n');

%% 2. 전압 구간별 용량 추세 (C-rate별)
fprintf('  Creating voltage range capacity trend plots...\n');

numCrates = length(target_crates);
chNum = channels(1);  % ch09만 처리
chName = sprintf('ch%02d', chNum);

for crateIdx = 1:numCrates
    targetCrate = target_crates(crateIdx);
    
    fig2 = figure('Position', [100, 100, 1600, 1000], 'Color', 'white');
    set(gcf, 'Visible', 'off');
    
    for rangeIdx = 1:size(voltage_ranges, 1)
        subplot(2, 3, rangeIdx);
        hold on;
        
        v_min = voltage_ranges(rangeIdx, 1);
        v_max = voltage_ranges(rangeIdx, 2);
        
        capacities = [];
        validCycles = [];
        
        for cycIdx = 1:length(rptCycles)
            cycleNum = rptCycles(cycIdx);
            cycleName = sprintf('rpt%d', cycleNum);
            
            if isfield(voltage_capacity_data.(chName), cycleName)
                data = voltage_capacity_data.(chName).(cycleName);
                stepNames = fieldnames(data);
                
                % 해당 C-rate의 용량 평균
                capacity_sum = 0;
                step_count = 0;
                
                for stepIdx = 1:length(stepNames)
                    stepName = stepNames{stepIdx};
                    if strcmp(stepName, 'V') || strcmp(stepName, 'Q')
                        continue;
                    end
                    
                    if isfield(data, stepName)
                        stepData = data.(stepName);
                        
                        % StepType 필터링
                        if isfield(stepData, 'StepType')
                            stepType = stepData.StepType;
                            if strcmp(use_step_type, 'Charge') && stepType ~= 'C'
                                continue;  % 충전만 사용
                            elseif strcmp(use_step_type, 'Discharge') && stepType ~= 'D'
                                continue;  % 방전만 사용
                            end
                        end
                        
                        % C-rate 필터링
                        if isfield(stepData, 'Crate_normalized')
                            if abs(stepData.Crate_normalized - targetCrate) > 0.01
                                continue;  % 해당 C-rate가 아니면 스킵
                            end
                        end
                        
                        if isfield(stepData, 'capacity_by_range')
                            capacity_sum = capacity_sum + stepData.capacity_by_range(rangeIdx);
                            step_count = step_count + 1;
                        end
                    end
                end
                
                if step_count > 0
                    capacities(end+1) = capacity_sum / step_count;  % 평균 사용
                    validCycles(end+1) = cycleNum;
                end
            end
        end
        
        if ~isempty(capacities)
            plot(validCycles, capacities, 'o-', 'LineWidth', 1.5, ...
                'MarkerSize', 6, 'Color', [0.2 0.6 0.8], ...
                'MarkerFaceColor', [0.2 0.6 0.8]);
        end
        
        xlabel('Cycle', 'FontSize', 10);
        ylabel('Capacity [Ah]', 'FontSize', 10);
        title(sprintf('%.1f-%.1f V', v_min, v_max), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        hold off;
    end
    
    sgtitle(sprintf('Voltage Range Capacity Trend - %s (C-rate = %.1fC)', chName, targetCrate), 'FontSize', 14, 'FontWeight', 'bold');
    
    savefig(fig2, fullfile(outputDir, sprintf('VoltageRange_CapacityTrend_C%.1f.fig', targetCrate)));
    close(fig2);
    fprintf('    Saved: VoltageRange_CapacityTrend_C%.1f.fig\n', targetCrate);
end


%% 3. dQdV 피크 분석 시각화
fprintf('  Creating dQdV peak analysis plots...\n');

fig3 = figure('Position', [100, 100, 2000, 1200], 'Color', 'white');
set(gcf, 'Visible', 'off');

% dQdV 피크 분석 (C-rate별 서브플롯)
numCrates = length(target_crates);
cols = 3;
rows = 2;
cycleColors = lines(length(rptCycles));

chNum = channels(1);  % ch09만 처리
chName = sprintf('ch%02d', chNum);

for crateIdx = 1:numCrates
    targetCrate = target_crates(crateIdx);
    
    subplot(rows, cols, crateIdx);
    hold on;
    
    for cycIdx = 1:length(rptCycles)
        cycleNum = rptCycles(cycIdx);
        cycleName = sprintf('rpt%d', cycleNum);
        
        if isfield(dQdV_peak_data.(chName), cycleName)
            data = dQdV_peak_data.(chName).(cycleName);
            stepNames = fieldnames(data);
            
            for stepIdx = 1:length(stepNames)
                stepName = stepNames{stepIdx};
                stepData = data.(stepName);
                
                % StepType 필터링
                if isfield(stepData, 'StepType')
                    stepType = stepData.StepType;
                    if strcmp(use_step_type, 'Charge') && stepType ~= 'C'
                        continue;  % 충전만 사용
                    elseif strcmp(use_step_type, 'Discharge') && stepType ~= 'D'
                        continue;  % 방전만 사용
                    end
                end
                
                % C-rate 필터링
                if isfield(stepData, 'Crate_normalized')
                    if abs(stepData.Crate_normalized - targetCrate) > 0.01
                        continue;  % 해당 C-rate가 아니면 스킵
                    end
                end
                
                if isfield(stepData, 'V') && isfield(stepData, 'dQdV')
                    V_plot = stepData.V;
                    dQdV_plot = stepData.dQdV;
                    
                    % 샘플링
                    if length(V_plot) > 3000
                        idx = round(linspace(1, length(V_plot), 3000));
                        V_plot = V_plot(idx);
                        dQdV_plot = dQdV_plot(idx);
                    end
                    
                    plot(V_plot, dQdV_plot, '-', 'LineWidth', 1.5, ...
                        'Color', cycleColors(cycIdx, :));
                    
                    % 피크 표시 및 전압 구간 표기
                    if isfield(stepData, 'peak_V') && ~isempty(stepData.peak_V)
                        peak_V = stepData.peak_V;
                        peak_dQdV = stepData.peak_dQdV;
                        scatter(peak_V, peak_dQdV, 100, 'filled', ...
                            'MarkerFaceColor', cycleColors(cycIdx, :), ...
                            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                        
                        % 각 피크의 전압 구간 표시 (수직선 및 배경)
                        if isfield(stepData, 'peak_voltage_ranges') && ~isempty(stepData.peak_voltage_ranges)
                            y_lim = ylim;
                            for p = 1:length(peak_V)
                                v_range = stepData.peak_voltage_ranges{p};
                                if ~isempty(v_range) && length(v_range) == 2
                                    % 전압 구간을 배경으로 표시
                                    v_min = v_range(1);
                                    v_max = v_range(2);
                                    fill([v_min, v_max, v_max, v_min], ...
                                        [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                                        [0.9 0.9 0.9], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                                    
                                    % 전압 구간 경계선 표시
                                    plot([v_min, v_min], y_lim, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
                                    plot([v_max, v_max], y_lim, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    xlabel('Voltage [V]', 'FontSize', 10);
    ylabel('dQ/dV [Ah/V]', 'FontSize', 10);
    title(sprintf('%s: dQdV Peaks (C-rate = %.1fC)', chName, targetCrate), 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    hold off;
end

if strcmp(use_step_type, 'Charge')
    sgtitle('dQdV Peak Analysis - Multi-rate CC (Charge Only)', 'FontSize', 14, 'FontWeight', 'bold');
elseif strcmp(use_step_type, 'Discharge')
    sgtitle('dQdV Peak Analysis - Multi-rate CC (Discharge Only)', 'FontSize', 14, 'FontWeight', 'bold');
else
    sgtitle('dQdV Peak Analysis - Multi-rate CC (Charge & Discharge)', 'FontSize', 14, 'FontWeight', 'bold');
end

savefig(fig3, fullfile(outputDir, 'dQdV_Peak_Analysis.fig'));
close(fig3);
fprintf('    Saved: dQdV_Peak_Analysis.fig\n');

%% 4. dQdV 피크 높이 및 면적 추세 (C-rate별)
fprintf('  Creating dQdV peak height and area trend plots...\n');

fig4 = figure('Position', [100, 100, 2000, 1200], 'Color', 'white');
set(gcf, 'Visible', 'off');

numCrates = length(target_crates);
chNum = channels(1);  % ch09만 처리
chName = sprintf('ch%02d', chNum);
cycleColors = lines(length(rptCycles));

% 피크 높이 추세 (C-rate별)
for crateIdx = 1:numCrates
    targetCrate = target_crates(crateIdx);
    
    subplot(2, numCrates, crateIdx);
    hold on;
    
    peak_heights = [];
    validCycles = [];
    
    for cycIdx = 1:length(rptCycles)
        cycleNum = rptCycles(cycIdx);
        cycleName = sprintf('rpt%d', cycleNum);
        
        if isfield(dQdV_peak_data.(chName), cycleName)
            data = dQdV_peak_data.(chName).(cycleName);
            stepNames = fieldnames(data);
            
            % 해당 C-rate의 최대 피크 높이 찾기
            max_peak_height = 0;
            has_peak = false;
            
            for stepIdx = 1:length(stepNames)
                stepName = stepNames{stepIdx};
                stepData = data.(stepName);
                
                % StepType 필터링
                if isfield(stepData, 'StepType')
                    stepType = stepData.StepType;
                    if strcmp(use_step_type, 'Charge') && stepType ~= 'C'
                        continue;  % 충전만 사용
                    elseif strcmp(use_step_type, 'Discharge') && stepType ~= 'D'
                        continue;  % 방전만 사용
                    end
                end
                
                % C-rate 필터링
                if isfield(stepData, 'Crate_normalized')
                    if abs(stepData.Crate_normalized - targetCrate) > 0.01
                        continue;  % 해당 C-rate가 아니면 스킵
                    end
                end
                
                if isfield(stepData, 'peak_dQdV') && ~isempty(stepData.peak_dQdV)
                    max_peak_height = max(max_peak_height, max(stepData.peak_dQdV));
                    has_peak = true;
                end
            end
            
            if has_peak
                peak_heights(end+1) = max_peak_height;
                validCycles(end+1) = cycleNum;
            end
        end
    end
    
    if ~isempty(peak_heights)
        plot(validCycles, peak_heights, 'o-', 'LineWidth', 1.5, ...
            'MarkerSize', 6, 'Color', [0.2 0.4 0.8], ...
            'MarkerFaceColor', [0.2 0.4 0.8]);
    end
    
    xlabel('Cycle', 'FontSize', 10);
    ylabel('Peak Height [Ah/V]', 'FontSize', 10);
    title(sprintf('Peak Height (C-rate = %.1fC, Voltage Range: %.1f-%.1fV)', ...
        targetCrate, min(voltage_ranges(:)), max(voltage_ranges(:))), 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    hold off;
end

% 피크 면적 추세 (C-rate별)
for crateIdx = 1:numCrates
    targetCrate = target_crates(crateIdx);
    
    subplot(2, numCrates, numCrates + crateIdx);
    hold on;
    
    peak_areas = [];
    validCycles = [];
    
    for cycIdx = 1:length(rptCycles)
        cycleNum = rptCycles(cycIdx);
        cycleName = sprintf('rpt%d', cycleNum);
        
        if isfield(dQdV_peak_data.(chName), cycleName)
            data = dQdV_peak_data.(chName).(cycleName);
            stepNames = fieldnames(data);
            
            % 해당 C-rate의 총 피크 면적 합계
            total_peak_area = 0;
            has_area = false;
            
            for stepIdx = 1:length(stepNames)
                stepName = stepNames{stepIdx};
                stepData = data.(stepName);
                
                % StepType 필터링
                if isfield(stepData, 'StepType')
                    stepType = stepData.StepType;
                    if strcmp(use_step_type, 'Charge') && stepType ~= 'C'
                        continue;  % 충전만 사용
                    elseif strcmp(use_step_type, 'Discharge') && stepType ~= 'D'
                        continue;  % 방전만 사용
                    end
                end
                
                % C-rate 필터링
                if isfield(stepData, 'Crate_normalized')
                    if abs(stepData.Crate_normalized - targetCrate) > 0.01
                        continue;  % 해당 C-rate가 아니면 스킵
                    end
                end
                
                if isfield(stepData, 'peak_area') && ~isempty(stepData.peak_area)
                    % peak_area가 0이 아닌 값이 있는지 확인
                    non_zero_area = stepData.peak_area(stepData.peak_area ~= 0);
                    if ~isempty(non_zero_area)
                        % 각 피크의 면적을 합산 (전압 구간별로 계산된 면적)
                        total_peak_area = total_peak_area + sum(non_zero_area);
                        has_area = true;
                    end
                end
            end
            
            if has_area
                peak_areas(end+1) = total_peak_area;
                validCycles(end+1) = cycleNum;
            end
        end
    end
    
    if ~isempty(peak_areas)
        plot(validCycles, peak_areas, 'o-', 'LineWidth', 1.5, ...
            'MarkerSize', 6, 'Color', [0.8 0.2 0.2], ...
            'MarkerFaceColor', [0.8 0.2 0.2]);
        fprintf('    C-rate %.1fC: %d cycles with peak area data\n', targetCrate, length(peak_areas));
    else
        fprintf('    C-rate %.1fC: No peak area data found\n', targetCrate);
    end
    
    xlabel('Cycle', 'FontSize', 10);
    ylabel('Peak Area [Ah]', 'FontSize', 10);
    title(sprintf('Peak Area (C-rate = %.1fC, Voltage Range: %.1f-%.1fV)', ...
        targetCrate, min(voltage_ranges(:)), max(voltage_ranges(:))), 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    hold off;
end

if strcmp(use_step_type, 'Charge')
    sgtitle(sprintf('dQdV Peak Height and Area Trends - %s (Charge Only, C-rate separated)', chName), 'FontSize', 14, 'FontWeight', 'bold');
elseif strcmp(use_step_type, 'Discharge')
    sgtitle(sprintf('dQdV Peak Height and Area Trends - %s (Discharge Only, C-rate separated)', chName), 'FontSize', 14, 'FontWeight', 'bold');
else
    sgtitle(sprintf('dQdV Peak Height and Area Trends - %s (Charge & Discharge, C-rate separated)', chName), 'FontSize', 14, 'FontWeight', 'bold');
end

savefig(fig4, fullfile(outputDir, 'dQdV_Peak_Trends.fig'));
close(fig4);
fprintf('    Saved: dQdV_Peak_Trends.fig\n');

%% 데이터 저장
fprintf('\n=== Saving Data ===\n');
save(fullfile(outputDir, 'RPT_Analysis_Data.mat'), ...
    'voltage_capacity_data', 'dQdV_peak_data', 'voltage_ranges', ...
    'rptCycles', 'channels', '-v7.3');
fprintf('  Saved: RPT_Analysis_Data.mat\n');

fprintf('\n=== Lab RPT Analysis Complete ===\n');
fprintf('Results saved to: %s\n', outputDir);

%% Helper Function: dQdV 피크 찾기
function [peak_V, peak_dQdV, peak_area, peak_voltage_ranges] = find_dQdV_peaks(V, dQdV, min_height, min_distance, voltage_ranges)
    peak_V = [];
    peak_dQdV = [];
    peak_area = [];
    peak_voltage_ranges = {};
    
    if isempty(V) || isempty(dQdV) || length(V) ~= length(dQdV)
        return;
    end
    
    % findpeaks 함수 사용
    if exist('findpeaks', 'file') == 2
        try
            % 데이터의 최대값 확인 (절댓값 사용)
            max_dQdV = max(abs(dQdV));
            
            % MinPeakHeight가 데이터 최대값보다 크면 조정
            % 실제 최대값의 일부를 사용하도록 자동 조정
            if min_height > max_dQdV
                % min_height를 최대값의 10%로 조정
                min_height = max(min_height * 0.1, max_dQdV * 0.05);
            end
            
            % 최소 피크 거리를 샘플 수로 변환
            if length(V) > 1
                dV_avg = mean(abs(diff(V)));
                min_distance_samples = max(1, round(min_distance / dV_avg));
            else
                min_distance_samples = 1;
            end
            
            % 경고 억제하고 findpeaks 호출
            warning('off', 'signal:findpeaks:largeMinPeakHeight');
            
            % 절댓값으로 피크 찾기 (양수/음수 피크 모두 찾기)
            [peak_vals_pos, peak_idx_pos] = findpeaks(abs(dQdV), 'MinPeakHeight', min_height, ...
                'MinPeakDistance', min_distance_samples);
            
            % 원본 dQdV에서도 피크 찾기 (양수 피크)
            [peak_vals_orig, peak_idx_orig] = findpeaks(dQdV, 'MinPeakHeight', min_height, ...
                'MinPeakDistance', min_distance_samples);
            
            warning('on', 'signal:findpeaks:largeMinPeakHeight');
            
            % 두 결과를 합치기 (중복 제거)
            all_peak_idx = unique([peak_idx_pos(:); peak_idx_orig(:)]);
            
            if isempty(all_peak_idx)
                % 피크를 찾지 못한 경우
                peak_V = [];
                peak_dQdV = [];
                peak_area = [];
                peak_voltage_ranges = {};
                return;
            end
            
            peak_V = V(all_peak_idx);
            peak_dQdV = dQdV(all_peak_idx);
            
            % 각 피크의 면적 계산: 전압 구간별로 계산
            % 면적 공식: Area = ∫(V1 to V2) (dQ/dV) dV = ΔQ
            % 피크가 속한 전압 구간을 찾아서 해당 구간의 면적 계산
            peak_area = zeros(size(peak_V));  % 사전 할당
            peak_voltage_ranges = cell(size(peak_V));  % 각 피크의 전압 구간 저장
            
            % 전압 구간이 제공되지 않으면 기본값 사용
            if nargin < 5 || isempty(voltage_ranges)
                voltage_ranges = [
                    3.0, 3.2;   % Low voltage range
                    3.2, 3.4;   % Mid-low voltage range
                    3.4, 3.6;   % Mid voltage range
                    3.6, 3.8;   % Mid-high voltage range
                    3.8, 4.0;   % High voltage range
                    4.0, 4.2    % Very high voltage range
                ];
            end
            
            for p = 1:length(peak_V)
                v_peak = peak_V(p);
                
                % 피크가 속한 전압 구간 찾기
                found_range = false;
                for rangeIdx = 1:size(voltage_ranges, 1)
                    v_min = voltage_ranges(rangeIdx, 1);
                    v_max = voltage_ranges(rangeIdx, 2);
                    
                    if v_peak >= v_min && v_peak <= v_max
                        % 해당 전압 구간에서 면적 계산
                        idx_range = V >= v_min & V <= v_max;
                        if sum(idx_range) > 1
                            V_range = V(idx_range);
                            dQdV_range = dQdV(idx_range);
                            
                            % 면적 계산: ∫(dQ/dV) dV = ΔQ
                            area = trapz(V_range, dQdV_range);
                            peak_area(p) = abs(area);  % 면적은 양수로 저장
                            peak_voltage_ranges{p} = [v_min, v_max];  % 전압 구간 저장
                            found_range = true;
                            break;
                        end
                    end
                end
                
                % 피크가 어떤 구간에도 속하지 않으면 0
                if ~found_range
                    peak_area(p) = 0;
                    peak_voltage_ranges{p} = [];
                end
            end
            
        catch
            % Fallback: 간단한 피크 검출
            peak_V = [];
            peak_dQdV = [];
            peak_area = [];
            peak_voltage_ranges = {};
            peak_voltage_ranges = {};
            
            % 전압 구간이 제공되지 않으면 기본값 사용
            if nargin < 5 || isempty(voltage_ranges)
                voltage_ranges = [
                    3.0, 3.2; 3.2, 3.4; 3.4, 3.6;
                    3.6, 3.8; 3.8, 4.0; 4.0, 4.2
                ];
            end
            
            for i = 2:length(dQdV)-1
                if dQdV(i) > dQdV(i-1) && dQdV(i) > dQdV(i+1) && dQdV(i) >= min_height
                    if isempty(peak_V) || abs(V(i) - peak_V(end)) >= min_distance
                        peak_V = [peak_V; V(i)];
                        peak_dQdV = [peak_dQdV; dQdV(i)];
                        
                        % 피크가 속한 전압 구간 찾아서 면적 계산
                        v_peak = V(i);
                        area = 0;
                        v_range_found = [];
                        for rangeIdx = 1:size(voltage_ranges, 1)
                            v_min = voltage_ranges(rangeIdx, 1);
                            v_max = voltage_ranges(rangeIdx, 2);
                            if v_peak >= v_min && v_peak <= v_max
                                idx_range = V >= v_min & V <= v_max;
                                if sum(idx_range) > 1
                                    V_range = V(idx_range);
                                    dQdV_range = dQdV(idx_range);
                                    area = abs(trapz(V_range, dQdV_range));
                                    v_range_found = [v_min, v_max];
                                    break;
                                end
                            end
                        end
                        peak_area = [peak_area; area];
                        peak_voltage_ranges{end+1} = v_range_found;
                    end
                end
            end
        end
    else
        % Fallback: 간단한 피크 검출
        % 전압 구간이 제공되지 않으면 기본값 사용
        if nargin < 5 || isempty(voltage_ranges)
            voltage_ranges = [
                3.0, 3.2; 3.2, 3.4; 3.4, 3.6;
                3.6, 3.8; 3.8, 4.0; 4.0, 4.2
            ];
        end
        
        peak_voltage_ranges = {};
        for i = 2:length(dQdV)-1
            if dQdV(i) > dQdV(i-1) && dQdV(i) > dQdV(i+1) && dQdV(i) >= min_height
                if isempty(peak_V) || abs(V(i) - peak_V(end)) >= min_distance
                    peak_V = [peak_V; V(i)];
                    peak_dQdV = [peak_dQdV; dQdV(i)];
                    
                    % 피크가 속한 전압 구간 찾아서 면적 계산
                    v_peak = V(i);
                    area = 0;
                    v_range_found = [];
                    for rangeIdx = 1:size(voltage_ranges, 1)
                        v_min = voltage_ranges(rangeIdx, 1);
                        v_max = voltage_ranges(rangeIdx, 2);
                        if v_peak >= v_min && v_peak <= v_max
                            idx_range = V >= v_min & V <= v_max;
                            if sum(idx_range) > 1
                                V_range = V(idx_range);
                                dQdV_range = dQdV(idx_range);
                                area = abs(trapz(V_range, dQdV_range));
                                v_range_found = [v_min, v_max];
                                break;
                            end
                        end
                    end
                    peak_area = [peak_area; area];
                    peak_voltage_ranges{end+1} = v_range_found;
                end
            end
        end
    end
end
