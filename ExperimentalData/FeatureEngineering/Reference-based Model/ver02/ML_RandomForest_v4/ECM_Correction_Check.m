% ECM_Correction_Check.m — C-rate별 ECM 보정 전/후 V-Q 및 dQ/dV 시각화
% V_ocv에 movmean smoothing 적용 후 dQ/dV 계산
close all; clear;

d_vq  = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat');
d_ecm = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\ECM_2RC_AllChannels_cyc0.mat');
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4';

ch = 'Ch09'; Q_0 = 63.5;
ecm_chg = d_ecm.All_ECM.(ch).charge;
ecm_dch = d_ecm.All_ECM.(ch).discharge;

% 0.05C OCV reference curve
OCV_ref_key = 'OCV_charge';
if isfield(d_vq.RPT_VQ_grid.cyc0.(ch), OCV_ref_key)
    s_ocv = d_vq.RPT_VQ_grid.cyc0.(ch).(OCV_ref_key);
    [Vu_ref, u_ref] = unique(s_ocv.V_grid); Qu_ref = s_ocv.Q(u_ref);
    if Vu_ref(1)>Vu_ref(end), Vu_ref=flipud(Vu_ref); Qu_ref=flipud(Qu_ref); end
    dqv_ref = movmean(gradient(Qu_ref)./gradient(Vu_ref), 21);
    has_ref = true;
else
    has_ref = false;
end

crates  = {'c01','c05','c1','c2','c3'};
clabels = {'0.1C','0.5C','1C','2C','3C'};
sm_win  = 61;  % smoothing window for V_ocv (points)

%% ---- Charge ----
fig1 = figure('Name','ECM CHG','Position',[30,30,1600,700]);
for r = 1:5
    f = [crates{r} '_charge'];
    if ~isfield(d_vq.RPT_VQ_grid.cyc0.(ch), f), continue; end
    s = d_vq.RPT_VQ_grid.cyc0.(ch).(f);
    t_s = seconds(s.t_raw(:)); V_raw = double(s.V_raw(:));
    I_raw = double(s.I_raw(:)); Q_raw = s.Q_raw(:);

    % ECM correction
    V_ocv = ecm_correct(V_raw, I_raw, t_s, ecm_chg, Q_0, true);
    
    % Smooth V_ocv before V-grid interp
    V_ocv_sm = movmean(V_ocv, sm_win);

    % V-Q plot (윗줄)
    subplot(2,5,r);
    plot(Q_raw, V_raw, 'b-','LineWidth',1); hold on;
    plot(Q_raw, V_ocv_sm, 'r-','LineWidth',1.2);
    xlabel('Q(Ah)'); ylabel('V'); ylim([3.0 4.3]); grid on;
    title(sprintf('CHG %s',clabels{r}),'FontWeight','bold','FontSize',10);
    if r==1, legend('V_{term}','OCV_{est}(sm)','FontSize',7,'Location','northwest'); end

    % dQ/dV comparison (아랫줄)
    [Vu_r,u] = unique(s.V_grid); Qu_r = s.Q(u);
    if Vu_r(1)>Vu_r(end), Vu_r=flipud(Vu_r); Qu_r=flipud(Qu_r); end
    dqv_r = movmean(gradient(Qu_r)./gradient(Vu_r), 21);

    % V_ocv_sm → monotone → V-grid
    [Vu_o, u2] = unique(V_ocv_sm, 'stable'); Qu_o = Q_raw(u2);
    mo = true(size(Vu_o));
    for ii=2:length(Vu_o), if Vu_o(ii)<=Vu_o(ii-1), mo(ii)=false; end; end
    Vu_o=Vu_o(mo); Qu_o=Qu_o(mo);
    Vg=(min(Vu_o):0.001:max(Vu_o))'; Qg=interp1(Vu_o,Qu_o,Vg,'linear');
    ok=~isnan(Qg); Vg=Vg(ok); Qg=Qg(ok);
    dqv_o = movmean(gradient(Qg)./gradient(Vg), 21);

    subplot(2,5,5+r);
    plot(Vu_r, dqv_r, 'b-','LineWidth',1); hold on;
    plot(Vg, dqv_o, 'r-','LineWidth',1.5);
    if has_ref, plot(Vu_ref, dqv_ref, 'k--','LineWidth',1); end
    xlabel('V'); ylabel('dQ/dV'); grid on; xlim([3.3 4.2]);
    title(sprintf('dQ/dV CHG %s',clabels{r}),'FontSize',9);
    if r==1
        if has_ref, legend('V_{term}','OCV_{est}','OCV(0.05C)','FontSize',6,'Location','best');
        else, legend('V_{term}','OCV_{est}','FontSize',6,'Location','best'); end
    end
end
sgtitle(sprintf('[Ch09 cyc0] ECM Correction — Charge (smooth win=%d)',sm_win),'FontSize',12,'FontWeight','bold');
saveas(fig1, fullfile(saveDir,'ECM_Correction_Check_CHG.fig'));

%% ---- Discharge ----
fig2 = figure('Name','ECM DCH','Position',[30,30,1600,700]);

% OCV discharge reference
if isfield(d_vq.RPT_VQ_grid.cyc0.(ch), 'OCV_discharge')
    s_ocvd = d_vq.RPT_VQ_grid.cyc0.(ch).OCV_discharge;
    [Vu_refd, u_refd] = unique(s_ocvd.V_grid); Qu_refd = s_ocvd.Q(u_refd);
    if Vu_refd(1)>Vu_refd(end), Vu_refd=flipud(Vu_refd); Qu_refd=flipud(Qu_refd); end
    dqv_refd = movmean(gradient(Qu_refd)./gradient(Vu_refd), 21);
    has_refd = true;
else, has_refd = false; end

for r = 1:5
    f = [crates{r} '_discharge'];
    if ~isfield(d_vq.RPT_VQ_grid.cyc0.(ch), f), continue; end
    s = d_vq.RPT_VQ_grid.cyc0.(ch).(f);
    t_s = seconds(s.t_raw(:)); V_raw = double(s.V_raw(:));
    I_raw = double(s.I_raw(:)); Q_raw = s.Q_raw(:);

    V_ocv = ecm_correct(V_raw, I_raw, t_s, ecm_dch, Q_0, false);
    V_ocv_sm = movmean(V_ocv, sm_win);

    subplot(2,5,r);
    plot(Q_raw, V_raw, 'b-','LineWidth',1); hold on;
    plot(Q_raw, V_ocv_sm, 'r-','LineWidth',1.2);
    xlabel('Q(Ah)'); ylabel('V'); ylim([2.8 4.3]); grid on;
    title(sprintf('DCH %s',clabels{r}),'FontWeight','bold','FontSize',10);
    if r==1, legend('V_{term}','OCV_{est}(sm)','FontSize',7,'Location','northeast'); end

    [Vu_r,u] = unique(s.V_grid); Qu_r = s.Q(u);
    if Vu_r(1)>Vu_r(end), Vu_r=flipud(Vu_r); Qu_r=flipud(Qu_r); end
    dqv_r = movmean(gradient(Qu_r)./gradient(Vu_r), 21);

    [Vu_o, u2] = unique(V_ocv_sm, 'stable'); Qu_o = Q_raw(u2);
    mo = true(size(Vu_o));
    for ii=2:length(Vu_o), if Vu_o(ii)>=Vu_o(ii-1), mo(ii)=false; end; end
    Vu_o=Vu_o(mo); Qu_o=Qu_o(mo);
    Vg=(min(Vu_o):0.001:max(Vu_o))'; Qg=interp1(Vu_o,Qu_o,Vg,'linear');
    ok=~isnan(Qg); Vg=Vg(ok); Qg=Qg(ok);
    dqv_o = movmean(gradient(Qg)./gradient(Vg), 21);

    subplot(2,5,5+r);
    plot(Vu_r, dqv_r, 'b-','LineWidth',1); hold on;
    plot(Vg, dqv_o, 'r-','LineWidth',1.5);
    if has_refd, plot(Vu_refd, dqv_refd, 'k--','LineWidth',1); end
    xlabel('V'); ylabel('dQ/dV'); grid on; xlim([3.0 4.2]);
    title(sprintf('dQ/dV DCH %s',clabels{r}),'FontSize',9);
    if r==1
        if has_refd, legend('V_{term}','OCV_{est}','OCV(0.05C)','FontSize',6,'Location','best');
        else, legend('V_{term}','OCV_{est}','FontSize',6,'Location','best'); end
    end
end
sgtitle(sprintf('[Ch09 cyc0] ECM Correction — Discharge (smooth win=%d)',sm_win),'FontSize',12,'FontWeight','bold');
saveas(fig2, fullfile(saveDir,'ECM_Correction_Check_DCH.fig'));
fprintf('Done.\n');

%% --- helper: Steady-state ECM correction ---
% CC 데이터에서 RC는 이미 정상상태 → V_ocv = V - I × (R₀+R₁+R₂)
function V_ocv = ecm_correct(V_raw, I_raw, t_raw, ecm, Q_0, is_charge)
    I_abs = abs(I_raw);
    Q_cum = cumtrapz(t_raw, I_abs)/3600;
    if is_charge, SOC = Q_cum/Q_0*100;
    else, SOC = 100 - Q_cum/Q_0*100; end
    SOC = max(0, min(100, SOC));

    [soc_s, si] = sort(ecm.SOC);
    % R_total = R₀ + R₁ + R₂ (steady-state, SOC 10~95%만 사용)
    R_total_s = (ecm.R0(si) + ecm.R1(si) + ecm.R2(si)) / 1e3;  % Ω
    valid = soc_s >= 10 & soc_s <= 95;
    soc_s = soc_s(valid);
    R_total_s = R_total_s(valid);
    SOC_c = max(min(soc_s), min(max(soc_s), SOC));

    R_total = interp1(soc_s, R_total_s, SOC_c, 'linear');
    V_ocv = V_raw - I_raw .* R_total;
end
