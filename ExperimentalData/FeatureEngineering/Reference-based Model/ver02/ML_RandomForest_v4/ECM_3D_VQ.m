%% ECM_3D_VQ.m — ECM 보정 V-Q 데이터 저장 + OCV 기준 오차 분석
%  채널별, 사이클별, C-rate별 ECM 보정된 V-Q 데이터 저장
%  공통 V 그리드 위에서 OCV(0.05C) 기준 Q 값 비교 → C-rate 간 오차 정량화
%  기준: OCV_charge / OCV_discharge (0.05C 측정)
close all; clear; clc;

%% ========================================================================
%  Section 0: Data Load
% =========================================================================
d_vq  = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat');
d_ecm = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\ECM_2RC_AllChannels_cyc0.mat');
d_cap = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Capacity_Trend_Figures\Capacity_Data_Static.mat');

saveDir = fileparts(mfilename('fullpath'));

%% ========================================================================
%  Section 1: Configuration
% =========================================================================
crates      = {'c01','c05','c1','c2','c3'};
crate_vals  = [0.1, 0.5, 1, 2, 3];
clabels     = {'0.1C','0.5C','1C','2C','3C'};
sm_win      = 61;
V_common    = (3.0 : 0.001 : 4.2)';  % 공통 V 그리드

channels   = fieldnames(d_ecm.All_ECM);
cyc_fields = fieldnames(d_vq.RPT_VQ_grid);

%% ========================================================================
%  Section 2: ECM 보정 + OCV 기준 오차 계산
% =========================================================================
VQ_ECM = struct();       % 보정된 V-Q 데이터

err_CellID    = {};
err_Cycle     = [];
err_CrateLabel = {};
err_CrateNum  = [];
err_Mode      = {};      % 'CHG' or 'DCH'
err_RMSE_Q    = [];      % Q 오차 RMSE (Ah)
err_MAE_Q     = [];      % Q 오차 MAE (Ah)
err_MaxErr_Q  = [];      % Q 최대 오차 (Ah)
err_RMSE_pct  = [];      % 오차 % (Q_0 대비)

cnt = 0;

for ch_idx = 1:length(channels)
    ch = channels{ch_idx};
    fprintf('Processing %s ...\n', ch);

    % --- Q_0 ---
    if ~isfield(d_cap.allChannelsCapacity, ch), continue; end
    cap_data = d_cap.allChannelsCapacity.(ch);
    idx_cyc0 = find(cap_data.cycles == 0, 1);
    if isempty(idx_cyc0), idx_cyc0 = 1; end
    Q_0 = max(cap_data.Q{1, idx_cyc0});
    if isnan(Q_0) || isempty(Q_0), continue; end

    ecm_chg = d_ecm.All_ECM.(ch).charge;
    ecm_dch = d_ecm.All_ECM.(ch).discharge;

    for ci = 1:length(cyc_fields)
        cyc_key = cyc_fields{ci};
        cyc_num = sscanf(cyc_key, 'cyc%d');
        if ~isfield(d_vq.RPT_VQ_grid.(cyc_key), ch), continue; end
        ch_data = d_vq.RPT_VQ_grid.(cyc_key).(ch);

        % --- 각 mode(충전/방전)별 처리 ---
        for mode = 1:2
            if mode == 1
                suffix = '_charge'; mode_str = 'CHG';
                ecm_p = ecm_chg; is_chg = true;
                ocv_field = 'OCV_charge';
            else
                suffix = '_discharge'; mode_str = 'DCH';
                ecm_p = ecm_dch; is_chg = false;
                ocv_field = 'OCV_discharge';
            end

            % === OCV 기준 데이터 (0.05C) ===
            if ~isfield(ch_data, ocv_field), continue; end
            s_ocv = ch_data.(ocv_field);

            % OCV V-Q를 공통 V 그리드에 보간
            V_ocv_ref = double(s_ocv.V_grid(:));
            Q_ocv_ref = double(s_ocv.Q(:));
            % 오름차순 정렬
            if V_ocv_ref(1) > V_ocv_ref(end)
                V_ocv_ref = flipud(V_ocv_ref);
                Q_ocv_ref = flipud(Q_ocv_ref);
            end
            % Q를 0부터 시작
            Q_ocv_ref = Q_ocv_ref - min(Q_ocv_ref);

            valid_ocv = V_common >= min(V_ocv_ref) & V_common <= max(V_ocv_ref);
            Q_ref_on_grid = nan(size(V_common));
            Q_ref_on_grid(valid_ocv) = interp1(V_ocv_ref, Q_ocv_ref, V_common(valid_ocv), 'linear');

            % OCV 데이터 저장
            VQ_ECM.(ch).(cyc_key).(mode_str).OCV.V = V_ocv_ref;
            VQ_ECM.(ch).(cyc_key).(mode_str).OCV.Q = Q_ocv_ref;
            VQ_ECM.(ch).(cyc_key).(mode_str).OCV.Q_on_Vgrid = Q_ref_on_grid;

            % === 각 C-rate: ECM 보정 후 비교 ===
            for r = 1:length(crates)
                f = [crates{r} suffix];
                if ~isfield(ch_data, f), continue; end
                s = ch_data.(f);

                V_raw = double(s.V_raw(:));
                I_raw = double(s.I_raw(:));
                t_s   = seconds(s.t_raw(:));
                Q_raw = double(s.Q_raw(:));

                % ECM correction
                V_corr = ecm_ss(V_raw, I_raw, t_s, ecm_p, Q_0, is_chg);
                V_corr = movmean(V_corr, sm_win);

                % Q를 0부터 시작
                Q_plot = Q_raw - min(Q_raw);

                % 단조(monotone) V-Q 만들기
                [V_u, uid] = unique(V_corr, 'stable');
                Q_u = Q_plot(uid);
                mono = true(size(V_u));
                if is_chg
                    for ii = 2:length(V_u)
                        if V_u(ii) <= V_u(ii-1), mono(ii) = false; end
                    end
                else
                    for ii = 2:length(V_u)
                        if V_u(ii) >= V_u(ii-1), mono(ii) = false; end
                    end
                end
                V_u = V_u(mono); Q_u = Q_u(mono);
                if V_u(1) > V_u(end)
                    V_u = flipud(V_u); Q_u = flipud(Q_u);
                end

                % 공통 V 그리드에 보간
                Q_on_grid = nan(size(V_common));
                valid_r = V_common >= min(V_u) & V_common <= max(V_u);
                Q_on_grid(valid_r) = interp1(V_u, Q_u, V_common(valid_r), 'linear');

                % VQ 데이터 저장
                VQ_ECM.(ch).(cyc_key).(mode_str).(crates{r}).V_ocv = V_corr;
                VQ_ECM.(ch).(cyc_key).(mode_str).(crates{r}).Q     = Q_plot;
                VQ_ECM.(ch).(cyc_key).(mode_str).(crates{r}).Q_on_Vgrid = Q_on_grid;

                % --- OCV 기준 오차 계산 ---
                valid = ~isnan(Q_ref_on_grid) & ~isnan(Q_on_grid);
                if sum(valid) < 10, continue; end

                dQ = Q_on_grid(valid) - Q_ref_on_grid(valid);

                cnt = cnt + 1;
                err_CellID{cnt,1}     = ch;
                err_Cycle(cnt,1)      = cyc_num;
                err_CrateLabel{cnt,1} = clabels{r};
                err_CrateNum(cnt,1)   = crate_vals(r);
                err_Mode{cnt,1}       = mode_str;
                err_RMSE_Q(cnt,1)     = sqrt(mean(dQ.^2));
                err_MAE_Q(cnt,1)      = mean(abs(dQ));
                err_MaxErr_Q(cnt,1)   = max(abs(dQ));
                err_RMSE_pct(cnt,1)   = err_RMSE_Q(cnt,1) / Q_0 * 100;
            end
        end
    end
    fprintf('  Done: %s\n', ch);
end

%% ========================================================================
%  Section 3: 결과 테이블 + 요약 출력
% =========================================================================
ErrorTable = table(err_CellID, err_Cycle, err_Mode, err_CrateLabel, err_CrateNum, ...
    err_RMSE_Q, err_MAE_Q, err_MaxErr_Q, err_RMSE_pct, ...
    'VariableNames', {'CellID','Cycle','Mode','CrateLabel','CrateNum', ...
                      'RMSE_Q_Ah','MAE_Q_Ah','MaxErr_Q_Ah','RMSE_pct'});

fprintf('\n=============================================\n');
fprintf('  ECM Correction Error Summary (vs OCV 0.05C)\n');
fprintf('=============================================\n');

modes = {'CHG','DCH'};
for m = 1:2
    fprintf('\n--- %s ---\n', modes{m});
    idx_m = strcmp(ErrorTable.Mode, modes{m});
    for r = 1:length(crates)
        idx_r = idx_m & ErrorTable.CrateNum == crate_vals(r);
        if sum(idx_r) == 0, continue; end
        sub = ErrorTable(idx_r, :);
        fprintf('  %4s: RMSE=%.3f Ah (%.2f%%), MAE=%.3f Ah, MaxErr=%.3f Ah  [N=%d]\n', ...
            clabels{r}, mean(sub.RMSE_Q_Ah), mean(sub.RMSE_pct), ...
            mean(sub.MAE_Q_Ah), mean(sub.MaxErr_Q_Ah), height(sub));
    end
end

% 저장
save(fullfile(saveDir, 'ECM_VQ_CorrectedData.mat'), 'VQ_ECM', 'V_common', '-v7.3');
save(fullfile(saveDir, 'ECM_VQ_ErrorAnalysis.mat'), 'ErrorTable');
fprintf('\nSaved:\n');
fprintf('  Data : ECM_VQ_CorrectedData.mat\n');
fprintf('  Error: ECM_VQ_ErrorAnalysis.mat\n');

%% ========================================================================
%  Helper: ECM Steady-State Correction
% =========================================================================
function V_ocv = ecm_ss(V_raw, I_raw, t_raw, ecm, Q_0, is_charge)
    I_abs = abs(I_raw);
    Q_cum = cumtrapz(t_raw, I_abs) / 3600;
    if is_charge
        SOC = Q_cum / Q_0 * 100;
    else
        SOC = 100 - Q_cum / Q_0 * 100;
    end
    SOC = max(0, min(100, SOC));

    [soc_s, si] = sort(ecm.SOC);
    R_tot_s = (ecm.R0(si) + ecm.R1(si) + ecm.R2(si)) / 1e3;  % mΩ → Ω
    valid = soc_s >= 10 & soc_s <= 95;
    soc_s = soc_s(valid);
    R_tot_s = R_tot_s(valid);

    SOC_c = max(min(soc_s), min(max(soc_s), SOC));
    R_tot = interp1(soc_s, R_tot_s, SOC_c, 'linear');
    V_ocv = V_raw - I_raw .* R_tot;
end
