%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase1_Extract_RPT_Features.m
% 목적: RPT 데이터에서 3가지 논문 방법으로 feature 추출
% - 논문 1 (PINN): 충전 마지막 구간에서 16개 통계 feature
% - 논문 2 (Segmentation): 11개 segment 충전 용량
% - 논문 3 (Power Curve): SOC별 충방전 가능 전력
% 
% 입력: RPT0~1000_ch09~16_parsed.mat (pdata 구조)
% 출력: 채널별 × RPT 시점별 feature 시각화
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 경로 설정
% 유틸리티 함수 경로
script_dir = fileparts(mfilename('fullpath'));
utils_dir = fullfile(script_dir, 'Utils');
if ~exist(utils_dir, 'dir')
    mkdir(utils_dir);
end
addpath(utils_dir);

% RPT 데이터 경로
rpt_base_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';

% 결과 저장 경로
result_path = fullfile(script_dir, 'Results', 'Phase1_RPT');
if ~exist(result_path, 'dir')
    mkdir(result_path);
    mkdir(fullfile(result_path, 'Data'));
    mkdir(fullfile(result_path, 'Figures', 'Feature_16D'));
    mkdir(fullfile(result_path, 'Figures', 'Segment'));
    mkdir(fullfile(result_path, 'Figures', 'PowerCurve'));
    mkdir(fullfile(result_path, 'Figures', 'Capacity'));
end

%% 파라미터 설정
rpt_cycles = [0, 200, 400, 600, 800, 1000];  % RPT 시점
channels = 9:16;  % 채널 번호
n_rpt = length(rpt_cycles);
n_channels = length(channels);

% Feature 개수
n_stat_features = 16;  % 통계 feature
n_segments = 11;  % Segment 개수
n_soc_bins = 101;  % SOC 0~100%

fprintf('====================================================\n');
fprintf('Phase 1: RPT Feature 추출 시작\n');
fprintf('====================================================\n');
fprintf('RPT 시점: %d개 [%s]\n', n_rpt, num2str(rpt_cycles));
fprintf('채널: ch%02d ~ ch%02d (%d개)\n', channels(1), channels(end), n_channels);
fprintf('총 파일: %d개\n', n_rpt * n_channels);
fprintf('====================================================\n\n');

%% 데이터 구조 초기화
Capacity_Charge_RPT = nan(n_rpt, n_channels);      % 충전 용량
Capacity_Discharge_RPT = nan(n_rpt, n_channels);   % 방전 용량 (Ground Truth)
SOH_RPT = nan(n_rpt, n_channels);                  % SOH
Features_16D_RPT = nan(n_rpt, n_channels, n_stat_features);
Q_segment_RPT = nan(n_rpt, n_channels, n_segments);
PowerCurve_Discharge_RPT = nan(n_rpt, n_channels, n_soc_bins);
PowerCurve_Charge_RPT = nan(n_rpt, n_channels, n_soc_bins);

%% Segment 전압 범위 결정
fprintf('[1/4] Segment 전압 범위 결정 중...\n');
sample_file = fullfile(rpt_base_path, 'RPT0_ch09_parsed.mat');

if ~exist(sample_file, 'file')
    error('샘플 파일이 존재하지 않습니다: %s', sample_file);
end

sample = load(sample_file);
if ~isfield(sample, 'pdata')
    error('pdata 필드가 없습니다.');
end

% 충전 step에서 전압 범위 확인
V_min_global = inf;
V_max_global = -inf;

for i = 1:length(sample.pdata)
    if sample.pdata(i).type == 'C'  % 충전
        V_min_global = min(V_min_global, min(sample.pdata(i).V));
        V_max_global = max(V_max_global, max(sample.pdata(i).V));
    end
end

if isinf(V_min_global) || isinf(V_max_global)
    error('충전 구간에서 전압 범위를 찾을 수 없습니다.');
end

V_segments = linspace(V_min_global, V_max_global, n_segments + 1);
fprintf('  전압 범위: %.3f V ~ %.3f V\n', V_min_global, V_max_global);
fprintf('  Segment: %d개\n', n_segments);
for seg = 1:n_segments
    fprintf('    Seg%02d: %.3f ~ %.3f V\n', seg, V_segments(seg), V_segments(seg+1));
end
fprintf('  완료!\n\n');

%% 메인 처리 루프
fprintf('[2/4] Feature 추출 시작...\n');
total_files = n_rpt * n_channels;
file_count = 0;
success_count = 0;

for rpt_idx = 1:n_rpt
    rpt_cycle = rpt_cycles(rpt_idx);
    
    for ch_idx = 1:n_channels
        ch_num = channels(ch_idx);
        file_count = file_count + 1;
        
        % 파일명 생성
        filename = sprintf('RPT%d_ch%02d_parsed.mat', rpt_cycle, ch_num);
        filepath = fullfile(rpt_base_path, filename);
        
        fprintf('  [%d/%d] %s\n', file_count, total_files, filename);
        
        % 파일 존재 확인
        if ~exist(filepath, 'file')
            fprintf('    경고: 파일 없음 - 스킵\n\n');
            continue;
        end
        
        try
            % 데이터 로드
            data_loaded = load(filepath);
            
            if ~isfield(data_loaded, 'pdata')
                fprintf('    경고: pdata 필드 없음 - 스킵\n\n');
                continue;
            end
            
            pdata = data_loaded.pdata;
            
            %% Feature 1: Ground Truth Capacity
            [cap_charge, cap_discharge] = extract_capacity(pdata);
            Capacity_Charge_RPT(rpt_idx, ch_idx) = cap_charge;
            Capacity_Discharge_RPT(rpt_idx, ch_idx) = cap_discharge;
            
            % SOH 계산 (첫 번째 방전 용량 기준)
            if rpt_idx == 1
                SOH_RPT(rpt_idx, ch_idx) = 100;
            else
                if ~isnan(Capacity_Discharge_RPT(1, ch_idx)) && Capacity_Discharge_RPT(1, ch_idx) > 0
                    SOH_RPT(rpt_idx, ch_idx) = (cap_discharge / Capacity_Discharge_RPT(1, ch_idx)) * 100;
                end
            end
            
            fprintf('    충전 용량: %.4f Ah, 방전 용량: %.4f Ah, SOH: %.2f%%\n', ...
                cap_charge, cap_discharge, SOH_RPT(rpt_idx, ch_idx));
            
            %% Feature 2: 논문 1 방법 - 16D Statistical Features
            features_16d = extract_16d_features(pdata, V_max_global);
            if ~isempty(features_16d) && length(features_16d) == n_stat_features
                Features_16D_RPT(rpt_idx, ch_idx, :) = features_16d;
                fprintf('    16D Features 추출 완료\n');
            else
                fprintf('    경고: 16D Features 추출 실패\n');
            end
            
            %% Feature 3: 논문 2 방법 - Segment Capacity
            q_segments = extract_segment_capacity(pdata, V_segments);
            if ~isempty(q_segments) && length(q_segments) == n_segments
                Q_segment_RPT(rpt_idx, ch_idx, :) = q_segments;
                fprintf('    Segment Features 추출 완료\n');
            else
                fprintf('    경고: Segment Features 추출 실패\n');
            end
            
            %% Feature 4: 논문 3 방법 - Power Curve
            [power_discharge, power_charge] = extract_power_curve(pdata, n_soc_bins);
            if ~isempty(power_discharge) && length(power_discharge) == n_soc_bins
                PowerCurve_Discharge_RPT(rpt_idx, ch_idx, :) = power_discharge;
                PowerCurve_Charge_RPT(rpt_idx, ch_idx, :) = power_charge;
                fprintf('    Power Curve 추출 완료\n');
            else
                fprintf('    경고: Power Curve 추출 실패\n');
            end
            
            success_count = success_count + 1;
            fprintf('    완료!\n\n');
            
        catch ME
            fprintf('    오류 발생: %s\n', ME.message);
            if ~isempty(ME.stack)
                fprintf('    위치: %s (line %d)\n\n', ME.stack(1).name, ME.stack(1).line);
            end
        end
    end
end

fprintf('Feature 추출 완료: %d/%d 성공\n\n', success_count, total_files);

%% 결과 저장
fprintf('[3/4] 결과 저장 중...\n');

save(fullfile(result_path, 'Data', 'Capacity_RPT.mat'), ...
    'Capacity_Charge_RPT', 'Capacity_Discharge_RPT', 'SOH_RPT', 'rpt_cycles', 'channels');
fprintf('  Capacity_RPT.mat 저장 완료\n');

save(fullfile(result_path, 'Data', 'Features_16D_RPT.mat'), ...
    'Features_16D_RPT', 'rpt_cycles', 'channels');
fprintf('  Features_16D_RPT.mat 저장 완료\n');

save(fullfile(result_path, 'Data', 'Q_segment_RPT.mat'), ...
    'Q_segment_RPT', 'V_segments', 'rpt_cycles', 'channels');
fprintf('  Q_segment_RPT.mat 저장 완료\n');

save(fullfile(result_path, 'Data', 'PowerCurve_RPT.mat'), ...
    'PowerCurve_Discharge_RPT', 'PowerCurve_Charge_RPT', 'rpt_cycles', 'channels');
fprintf('  PowerCurve_RPT.mat 저장 완료\n\n');

%% 시각화
fprintf('[4/4] 시각화 생성 중...\n');

% Feature 이름 정의
feature_names = {
    'V_mean', 'V_std', 'V_skewness', 'V_kurtosis', ...
    'V_time', 'V_Ah', 'V_slope', 'V_entropy', ...
    'I_mean', 'I_std', 'I_skewness', 'I_kurtosis', ...
    'I_time', 'I_Ah', 'I_slope', 'I_entropy'
};

% 1. 용량 변화 시각화
plot_capacity_trends(Capacity_Discharge_RPT, SOH_RPT, rpt_cycles, channels, result_path);

% 2. 16D Features 시각화 (채널별)
plot_16d_features(Features_16D_RPT, feature_names, rpt_cycles, channels, result_path);

% 3. Segment Features 시각화 (채널별)
plot_segment_features(Q_segment_RPT, V_segments, rpt_cycles, channels, result_path);

% 4. Power Curve 시각화 (채널별)
plot_power_curves(PowerCurve_Discharge_RPT, PowerCurve_Charge_RPT, rpt_cycles, channels, result_path);

fprintf('시각화 완료!\n\n');

fprintf('====================================================\n');
fprintf('Phase 1 완료!\n');
fprintf('결과 위치: %s\n', result_path);
fprintf('====================================================\n');