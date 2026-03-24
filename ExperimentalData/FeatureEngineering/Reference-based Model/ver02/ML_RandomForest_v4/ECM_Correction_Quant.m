% ECM_Correction_Quant.m — C-rate별 PeakPos 정량 비교
clear; clc;

d_vq  = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat');
d_ecm = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\ECM_2RC_AllChannels_cyc0.mat');

ch = 'Ch09'; Q_0 = 63.5;
ecm_chg = d_ecm.All_ECM.(ch).charge;
ecm_dch = d_ecm.All_ECM.(ch).discharge;
crates  = {'c01','c05','c1','c2','c3'};
clabels = {'0.1C','0.5C','1C','2C','3C'};
sm_win  = 61; ma_win = 21;

%% OCV reference (0.05C)
s_ref = d_vq.RPT_VQ_grid.cyc0.(ch).OCV_charge;
[Vu,u] = unique(s_ref.V_grid); Qu = s_ref.Q(u);
if Vu(1)>Vu(end), Vu=flipud(Vu); Qu=flipud(Qu); end
dqv = movmean(gradient(Qu)./gradient(Vu), ma_win);
[~,li] = max(abs(dqv)); PeakRef_chg = Vu(li);

s_refd = d_vq.RPT_VQ_grid.cyc0.(ch).OCV_discharge;
[Vud,ud] = unique(s_refd.V_grid); Qud = s_refd.Q(ud);
if Vud(1)>Vud(end), Vud=flipud(Vud); Qud=flipud(Qud); end
dqvd = movmean(gradient(Qud)./gradient(Vud), ma_win);
[~,lid] = max(abs(dqvd)); PeakRef_dch = Vud(lid);

fprintf('0.05C OCV Reference:  CHG PeakPos=%.4fV  DCH PeakPos=%.4fV\n\n', PeakRef_chg, PeakRef_dch);

%% CHARGE comparison
fprintf('=== CHARGE PeakPos ===\n');
fprintf('%-5s | %9s | %9s | %9s | %9s | Improved?\n', 'Rate','Orig(V)','ECM_ss(V)','Δ_orig','Δ_ecm');
fprintf('%s\n', repmat('-',1,60));
for r = 1:5
    f = [crates{r} '_charge'];
    if ~isfield(d_vq.RPT_VQ_grid.cyc0.(ch),f), continue; end
    s = d_vq.RPT_VQ_grid.cyc0.(ch).(f);
    t_s = seconds(s.t_raw(:)); V_raw = double(s.V_raw(:)); I_raw = double(s.I_raw(:)); Q_raw = s.Q_raw(:);

    % Original
    [Vu2,~] = unique(s.V_grid); Qu2 = s.Q(find(unique(s.V_grid,'stable')==Vu2));
    [Vu2,u2] = unique(s.V_grid); Qu2 = s.Q(u2);
    if Vu2(1)>Vu2(end), Vu2=flipud(Vu2); Qu2=flipud(Qu2); end
    dqv2 = movmean(gradient(Qu2)./gradient(Vu2), ma_win);
    [~,li2] = max(abs(dqv2)); PkOrig = Vu2(li2);

    % ECM ss
    V_ocv = ecm_ss(V_raw, I_raw, t_s, ecm_chg, Q_0, true);
    V_ocv_sm = movmean(V_ocv, sm_win);
    [Vu3,u3] = unique(V_ocv_sm,'stable'); Qu3 = Q_raw(u3);
    mo=true(size(Vu3)); for ii=2:length(Vu3), if Vu3(ii)<=Vu3(ii-1), mo(ii)=false; end; end
    Vu3=Vu3(mo); Qu3=Qu3(mo);
    Vg=(min(Vu3):0.001:max(Vu3))'; Qg=interp1(Vu3,Qu3,Vg,'linear'); ok=~isnan(Qg); Vg=Vg(ok); Qg=Qg(ok);
    dqv3 = movmean(gradient(Qg)./gradient(Vg), ma_win);
    [~,li3] = max(abs(dqv3)); PkECM = Vg(li3);

    d_orig = PkOrig - PeakRef_chg; d_ecm = PkECM - PeakRef_chg;
    improved = abs(d_ecm) < abs(d_orig);
    fprintf('%-5s | %9.4f | %9.4f | %+9.4f | %+9.4f | %s\n', clabels{r}, PkOrig, PkECM, d_orig, d_ecm, string(improved));
end

%% DISCHARGE comparison
fprintf('\n=== DISCHARGE PeakPos ===\n');
fprintf('%-5s | %9s | %9s | %9s | %9s | Improved?\n', 'Rate','Orig(V)','ECM_ss(V)','Δ_orig','Δ_ecm');
fprintf('%s\n', repmat('-',1,60));
for r = 1:5
    f = [crates{r} '_discharge'];
    if ~isfield(d_vq.RPT_VQ_grid.cyc0.(ch),f), continue; end
    s = d_vq.RPT_VQ_grid.cyc0.(ch).(f);
    t_s = seconds(s.t_raw(:)); V_raw = double(s.V_raw(:)); I_raw = double(s.I_raw(:)); Q_raw = s.Q_raw(:);

    [Vu2,u2] = unique(s.V_grid); Qu2 = s.Q(u2);
    if Vu2(1)>Vu2(end), Vu2=flipud(Vu2); Qu2=flipud(Qu2); end
    dqv2 = movmean(gradient(Qu2)./gradient(Vu2), ma_win);
    [~,li2] = max(abs(dqv2)); PkOrig = Vu2(li2);

    V_ocv = ecm_ss(V_raw, I_raw, t_s, ecm_dch, Q_0, false);
    V_ocv_sm = movmean(V_ocv, sm_win);
    [Vu3,u3] = unique(V_ocv_sm,'stable'); Qu3 = Q_raw(u3);
    mo=true(size(Vu3)); for ii=2:length(Vu3), if Vu3(ii)>=Vu3(ii-1), mo(ii)=false; end; end
    Vu3=Vu3(mo); Qu3=Qu3(mo);
    Vg=(min(Vu3):0.001:max(Vu3))'; Qg=interp1(Vu3,Qu3,Vg,'linear'); ok=~isnan(Qg); Vg=Vg(ok); Qg=Qg(ok);
    dqv3 = movmean(gradient(Qg)./gradient(Vg), ma_win);
    [~,li3] = max(abs(dqv3)); PkECM = Vg(li3);

    d_orig = PkOrig - PeakRef_dch; d_ecm = PkECM - PeakRef_dch;
    improved = abs(d_ecm) < abs(d_orig);
    fprintf('%-5s | %9.4f | %9.4f | %+9.4f | %+9.4f | %s\n', clabels{r}, PkOrig, PkECM, d_orig, d_ecm, string(improved));
end

%% Helper
function V_ocv = ecm_ss(V_raw, I_raw, t_raw, ecm, Q_0, is_charge)
    I_abs = abs(I_raw);
    Q_cum = cumtrapz(t_raw, I_abs)/3600;
    if is_charge, SOC = Q_cum/Q_0*100;
    else, SOC = 100 - Q_cum/Q_0*100; end
    SOC = max(0,min(100,SOC));

    [soc_s,si] = sort(ecm.SOC);
    R_tot_s = (ecm.R0(si)+ecm.R1(si)+ecm.R2(si))/1e3;

    % Remove extreme SOC points (< 10%, > 95%) where R_total is anomalous
    valid = soc_s >= 10 & soc_s <= 95;
    soc_s = soc_s(valid);
    R_tot_s = R_tot_s(valid);

    % Clamp SOC to valid LUT range
    SOC_c = max(min(soc_s), min(max(soc_s), SOC));
    R_tot = interp1(soc_s, R_tot_s, SOC_c, 'linear');
    V_ocv = V_raw - I_raw .* R_tot;
end

