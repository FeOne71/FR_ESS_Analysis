%% Generate_MasterRulers_v3.m
% 새 MasterRuler 생성 스크립트 (v3: Y2025 제외 + 물리적(용량 기준) 5분할)
%
% 변경 사유:
%   1. Y2025를 필드 공통 분석에서 제외함으로써 방전/충전 전압창 넓게 확보.
%      충전 공통: [3.674V, 3.976V]
%      방전 공통: [3.644V, 3.872V]
%   2. 기존 v2/v3 초안의 linspace(단순 전압 등분할) 방식을 폐기하고, 
%      논문 및 ver01의 물리적 근거를 반영하여 Static(0.1C) 데이터 기반
%      "동일 용량(시간) 분할" 방식으로 전압 경계선을 재설정함.
%
% 생성:
%   D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\MasterRulers_v3.mat

clear; clc; close all;

%% ========================================================================
% Section 1: 공통 전압창 정의 및 데이터 로드
% ========================================================================
NEW_CHG_MIN = 3.674;   % V
NEW_CHG_MAX = 3.976;   % V
NEW_DCH_MIN = 3.644;   % V  (방전 끝, 낮은 쪽)
NEW_DCH_MAX = 3.872;   % V  (방전 시작, 높은 쪽)
NUM_SEG = 5;           % 구간 수

baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_rpt_vq_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');

fprintf('Loading RPT_VQ_grid.mat for Static Data...\n');
load(path_rpt_vq_mat, 'RPT_VQ_grid');

% 우리는 MasterRuler_old의 채널 구조체가 필요합니다.
oldRulerPath = fullfile(baseDir, 'FeatureEngineering', 'MasterRulers_v2.mat');
loaded = load(oldRulerPath, 'MasterRulers');
MasterRulers_old = loaded.MasterRulers;
channels = fieldnames(MasterRulers_old);

%% ========================================================================
% Section 2: 평균 Static 곡선(용량-전압) 도출
% ========================================================================
% 채널들의 평균적인 형태를 얻기 위해 Interpolation을 수행합니다.
% Static 데이터는 단일 방전 곡선이므로 충전/방전 모두 이 궤적을 레퍼런스로 사용.

% 보간을 위한 절대적인 V Grid 기준점 생성 (세밀하게)
V_ref_chg = (NEW_CHG_MIN : 0.0001 : NEW_CHG_MAX)';
V_ref_dch = (NEW_DCH_MIN : 0.0001 : NEW_DCH_MAX)';

Q_chg_all = [];
Q_dch_all = [];

for i = 1:length(channels)
    ch = channels{i};
    if ~isfield(RPT_VQ_grid.cyc0, ch) || ~isfield(RPT_VQ_grid.cyc0.(ch), 'Static')
        continue;
    end
    
    s_data = RPT_VQ_grid.cyc0.(ch).Static;
    V_raw = s_data.V_raw;
    Q_raw = s_data.Q_raw;
    
    % 중복 제거 및 단조 증가 형태로 정렬 (interp1 요구조건)
    [V_u, idx_u] = unique(V_raw, 'stable');
    Q_u = Q_raw(idx_u);
    
    % 전압 오름차순
    [V_sort, idx_sort] = sort(V_u, 'ascend');
    Q_sort = Q_u(idx_sort);
    
    % Charge 윈도우 Q 보간
    q_c = interp1(V_sort, Q_sort, V_ref_chg, 'linear', 'extrap');
    Q_chg_all = [Q_chg_all, q_c];
    
    % Discharge 윈도우 Q 보간
    q_d = interp1(V_sort, Q_sort, V_ref_dch, 'linear', 'extrap');
    Q_dch_all = [Q_dch_all, q_d];
end

% 채널 평균 용량 궤적 (Q) 계산
Q_chg_avg = mean(Q_chg_all, 2, 'omitnan');
Q_dch_avg = mean(Q_dch_all, 2, 'omitnan');

% 상대 용량 크기로 변환 (시작점을 0으로)
Q_chg_avg = Q_chg_avg - min(Q_chg_avg);
Q_dch_avg = Q_dch_avg - min(Q_dch_avg);

%% ========================================================================
% Section 3: 등용량(Equal-Capacity) 기준 5등분 전압 경계점 탐색
% ========================================================================
% Charge Master Rulers (Ascending Voltage)
target_Q_chg = linspace(min(Q_chg_avg), max(Q_chg_avg), NUM_SEG + 1);
new_V_bounds_chg = interp1(Q_chg_avg, V_ref_chg, target_Q_chg, 'linear');
new_V_bounds_chg = sort(new_V_bounds_chg, 'ascend');

% Discharge Master Rulers (Descending Voltage)
target_Q_dch = linspace(min(Q_dch_avg), max(Q_dch_avg), NUM_SEG + 1);
new_V_bounds_dch = interp1(Q_dch_avg, V_ref_dch, target_Q_dch, 'linear');
new_V_bounds_dch = sort(new_V_bounds_dch, 'descend');

fprintf('\n=== New MasterRuler v3 Voltage Boundaries (Physical/Equal-Capacity) ===\n');
fprintf('Charge    (asc) :  '); fprintf('%.4fV  ', new_V_bounds_chg); fprintf('\n');
fprintf('Discharge (desc):  '); fprintf('%.4fV  ', new_V_bounds_dch); fprintf('\n');

%% ========================================================================
% Section 4: MasterRulers_v3 구조체 생성 및 저장
% ========================================================================
MasterRulers = struct();

for k = 1:length(channels)
    ch = channels{k};
    
    % 기존 메타데이터 보존
    old_fields = fieldnames(MasterRulers_old.(ch));
    for f = 1:length(old_fields)
        fn = old_fields{f};
        MasterRulers.(ch).(fn) = MasterRulers_old.(ch).(fn);
    end
    
    % 새 경계점으로 덮어쓰기
    MasterRulers.(ch).V_bounds_chg = new_V_bounds_chg;
    MasterRulers.(ch).V_bounds_dch = new_V_bounds_dch;
end

saveDir  = fullfile(baseDir, 'FeatureEngineering');
savePath = fullfile(saveDir, 'MasterRulers_v3.mat');
save(savePath, 'MasterRulers');

fprintf('\n=== Saved physically segmented MasterRulers to: %s ===\n', savePath);
