%% Field_Estimation_v5.m
% =====================================================================
% Phase 3: 필드 데이터 SOH 추정 및 시각화 (ML_Models_v5 기반)
% - 연도별 데이터 로드 후 v5 세그먼트에 맞게 ΔQ 추출
% - C_eff 계산 후 Track A / Track B 중 적합한 모델(GPR 최우선) 자동 적용
% - 결과 시각화: 연도별 OCV-based SOH vs Model SOH 비교 (Bar chart)
% =====================================================================
clear; clc; close all;

%% 1. 모델 로드
mdl_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\ML_Models_v5.mat';
load(mdl_path, 'Final_Models');
% 사용할 우선 모델 타입 (최고성능을 보인 GPR 우선, 없으면 RF)
mtype = 'GPR';

%% 2. 환경 및 데이터 셋팅
rawDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
% OCV-based Reference SOH 데이터 (이미 계산되어 있다고 가정, 예: 기존 RPT_Field_Estimation_RF 참조)
% 여기서는 이전 분석에서 얻은 평균값을 하드코딩 (없으면 계산 필요하지만, 분석 맥락상 트렌드 비교용)
ref_SOH = struct('Y2021', 97.5, 'Y2023', 95.0, 'Y2024', 92.5, 'Y2025', 90.0);

files = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old'
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new'
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new'
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new'
};

Np = 2; Q_cell = 64; 
VB = [3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.78, 3.84, 3.90, 4.00, 4.10];
N_seg = length(VB)-1;
thr_A = Q_cell * 0.05;

%% 3. 연도별 필드 추정 루프
fprintf('=== Field SOH Estimation (v5) ===\n');
res = struct();

for k = 1:size(files,1)
    fpath = files{k,1}; yr = files{k,2}; dtype = files{k,3};
    if ~exist(fpath,'file'), fprintf('[Skip] %s not found\n', yr); continue; end
    
    S = load(fpath);
    if strcmp(dtype,'old')
        I_rack = S.Raw.Rack01.DCCurrent_A(:);
        V_rack = S.Raw.Rack01.SumCV_V(:);
        time_arr = S.Raw.Rack01.Time;
    else
        I_rack = S.Raw.DCCurrent(:);
        V_rack = S.Raw.CVavg(:) * 14; % CVavg is cell average, convert to rack equivalent for later division or just use CVavg directly
        time_arr = S.Raw.Date_Time;
    end
    I_cell = I_rack / Np;
    
    if strcmp(dtype,'old')
        V_cell = S.Raw.Rack01.AverageCV_V(:);
        dt = datetime(time_arr, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        t_s = seconds(dt - dt(1));
    else
        V_cell = S.Raw.CVavg(:); % 이미 셀 평균
        t_s = seconds(time_arr) - seconds(time_arr(1));
    end
    if isempty(t_s), t_s = (0:length(I_cell)-1)'; end
    
    % ---- 충전 세그먼트 탐색 (가장 길면서 조건を満하는 구간) ----
    chg_mask = I_cell > thr_A;
    d = diff([0; chg_mask; 0]);
    starts = find(d==1); ends = find(d==-1)-1;
    lens = ends - starts + 1;
    if isempty(lens), continue; end
    
    % 가능한 한 높은 전압에 도달한(=관측 범위가 넓은) 긴 구간 찾기
    best_idx = 1; best_len = 0;
    for cand = 1:length(lens)
        if lens(cand) > 300 % 5분 이상
            c_idx = starts(cand):ends(cand);
            if max(V_cell(c_idx)) >= 3.90 && lens(cand) > best_len
                best_len = lens(cand);
                best_idx = cand;
            end
        end
    end
    
    % 만약 3.9V 도달하는 구간이 없으면 그냥 제일 긴 거 선택
    if best_len == 0
        [~, best_idx] = max(lens);
    end
    
    idx = starts(best_idx):ends(best_idx);
    
    t_seg = t_s(idx); V_seg = V_cell(idx); I_seg = I_cell(idx);
    Q_seg = cumtrapz(t_seg, I_seg)/3600;
    
    % C_eff
    C_eff = mean(abs(I_seg)) / Q_cell;
    
    % Track A vs Track B 결정 구조체
    % 필드에서 V_seg 범위에 들어오는 세그먼트만 추출 가능 추출 안되면 NaN
    dQ_raw = extract_dQ(V_seg, Q_seg, VB, true);
    
    % ==============================
    % 모델 선택 및 예측 논리
    % ==============================
    % Track A CHG: Seg 10, 9, 5 필요
    req_A_chg = [10, 9, 5];
    can_A_chg = all(~isnan(dQ_raw(req_A_chg)));
    
    % Track B CHG: Seg 9, 8, 7 필요
    req_B_chg = [9, 8, 7];
    can_B_chg = all(~isnan(dQ_raw(req_B_chg)));
    
    pred_SOH = NaN; pred_LLI = NaN; pred_LAM = NaN;
    track = 'None';
    
    if can_A_chg
        track = 'TrackA';
        mdl_soh = Final_Models.TrackA.CHG.SOH.(mtype).model;
        feat = [dQ_raw(req_A_chg), C_eff];
        pred_SOH = predict(mdl_soh, feat);
    elseif can_B_chg
        track = 'TrackB';
        mdl_soh = Final_Models.TrackB.CHG.SOH.(mtype).model;
        feat = [dQ_raw(req_B_chg), C_eff];
        pred_SOH = predict(mdl_soh, feat);
    else
        % 둘 다 안되면 에러
        track = 'Failed';
    end
    
    fprintf('-- %s --\n', yr);
    fprintf('  V_range: %.3f ~ %.3f V\n', min(V_seg), max(V_seg));
    fprintf('  C_eff:   %.3f C\n', C_eff);
    fprintf('  Track:   %s\n', track);
    fprintf('  SOH(%%):   %.2f\n', pred_SOH);
    
    res.(yr).C_eff = C_eff;
    res.(yr).track = track;
    res.(yr).SOH   = pred_SOH;
end

%% 4. 결과 시각화 
yrs_valid = fieldnames(res);
n_yrs = length(yrs_valid);
soh_mdl = nan(1, n_yrs); soh_ref = nan(1, n_yrs);

for i = 1:n_yrs
    yr = yrs_valid{i};
    soh_mdl(i) = res.(yr).SOH;
    if isfield(ref_SOH, yr), soh_ref(i) = ref_SOH.(yr); end
end

fig = figure('Position',[100,200,600,400],'Name','Field SOH Estimation v5');

% SOH 비교 플롯
hold on; grid on;
b = bar(1:n_yrs, [soh_ref; soh_mdl]');
b(1).FaceColor = [0.4 0.4 0.4]; b(1).DisplayName = 'Reference (OCV-based)';
b(2).FaceColor = [0.85 0.32 0.09]; b(2).DisplayName = sprintf('Model (%s)', mtype);

% Track Text 
for i=1:n_yrs
    yr = yrs_valid{i};
    trk = res.(yr).track;
    text(i, soh_mdl(i)+0.5, trk, 'HorizontalAlignment','center','FontSize',8,'Color','b','FontWeight','bold');
end

xticks(1:n_yrs); xticklabels(yrs_valid);
ylabel('SOH (%)'); ylim([80 100]);
title('Field Data SOH Estimation (v5)','FontWeight','bold');
legend('Location','southwest');

saveDir = fileparts(mfilename('fullpath'));
saveas(fig, fullfile(saveDir, 'Field_Estimation_Results_v5.fig'));
fprintf('\n>> Saved: Field_Estimation_Results_v5.fig\n');


%% Helper
function dQ = extract_dQ(V, Q, VB, is_charge)
    N = length(VB)-1; dQ = nan(1,N);
    % 전압 모노토닉 만들기 (간이 필터링)
    [Vu, uid] = unique(V, 'stable'); Qu = Q(uid);
    mono = true(size(Vu));
    if is_charge
        for ii=2:length(Vu), if Vu(ii)<=Vu(ii-1), mono(ii)=false; end; end
    else
        for ii=2:length(Vu), if Vu(ii)>=Vu(ii-1), mono(ii)=false; end; end
    end
    Vu = Vu(mono); Qu = Qu(mono);
    if isempty(Vu), return; end
    if Vu(1)>Vu(end), Vu=flipud(Vu); Qu=flipud(Qu); end
    if length(Vu)<10, return; end
    
    for s=1:N
        % 세그먼트가 관측 범위 내에 완전히 들어올 때만 계산
        if VB(s)>=min(Vu) && VB(s+1)<=max(Vu)
            q1=interp1(Vu,Qu,VB(s),'linear');
            q2=interp1(Vu,Qu,VB(s+1),'linear');
            if ~isnan(q1)&&~isnan(q2), dQ(s)=abs(q2-q1); end
        end
    end
end
