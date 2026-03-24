% ECM_3D_dQdV.m — 3D 시각화: V × Cycle × dQ/dV (ECM 보정 전/후)
% Ch09 충전 데이터로 시각화
close all; clear;

d_vq  = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\RPT_VQ_grid.mat');
d_ecm = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\ECM_2RC_AllChannels_cyc0.mat');
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4';

ch = 'Ch09'; Q_0 = 63.5;
ecm_chg = d_ecm.All_ECM.(ch).charge;
ecm_dch = d_ecm.All_ECM.(ch).discharge;
crates = {'c01','c05','c1','c2','c3'};
clabels = {'0.1C','0.5C','1C','2C','3C'};
colors_cr = lines(5);
ma_win = 21; sm_win = 61;

cyc_fields = fieldnames(d_vq.RPT_VQ_grid);

%% ---- Figure 1: 보정 전 (V_terminal 기준 dQ/dV) ----
fig1 = figure('Name','Before ECM','Position',[30,30,900,700]);
hold on; view([-30 30]); grid on;
xlabel('V (V)'); ylabel('Cycle'); zlabel('dQ/dV (Ah/V)');
title('[Ch09 CHG] dQ/dV — Before ECM Correction','FontSize',13,'FontWeight','bold');

for ci = 1:length(cyc_fields)
    cyc_key = cyc_fields{ci};
    cyc_num = sscanf(cyc_key, 'cyc%d');
    if ~isfield(d_vq.RPT_VQ_grid.(cyc_key), ch), continue; end
    ch_data = d_vq.RPT_VQ_grid.(cyc_key).(ch);

    for r = 1:length(crates)
        f = [crates{r} '_charge'];
        if ~isfield(ch_data, f), continue; end
        s = ch_data.(f);
        [Vu, u] = unique(s.V_grid); Qu = s.Q(u);
        if Vu(1)>Vu(end), Vu=flipud(Vu); Qu=flipud(Qu); end
        dqv = movmean(gradient(Qu)./gradient(Vu), ma_win);
        % 3.3~4.15V 범위만
        mask = Vu >= 3.3 & Vu <= 4.15;
        plot3(Vu(mask), cyc_num*ones(sum(mask),1), dqv(mask), ...
            '-','Color',colors_cr(r,:),'LineWidth',1.2);
    end
end
% Legend
for r = 1:5
    plot3(NaN,NaN,NaN,'-','Color',colors_cr(r,:),'LineWidth',2,'DisplayName',clabels{r});
end
legend('Location','northeast','FontSize',9);
saveas(fig1, fullfile(saveDir,'ECM_3D_dQdV_Before.fig'));

%% ---- Figure 2: 보정 후 (ECM Steady-State) ----
fig2 = figure('Name','After ECM','Position',[30,30,900,700]);
hold on; view([-30 30]); grid on;
xlabel('V (V)'); ylabel('Cycle'); zlabel('dQ/dV (Ah/V)');
title('[Ch09 CHG] dQ/dV — After ECM Correction','FontSize',13,'FontWeight','bold');

for ci = 1:length(cyc_fields)
    cyc_key = cyc_fields{ci};
    cyc_num = sscanf(cyc_key, 'cyc%d');
    if ~isfield(d_vq.RPT_VQ_grid.(cyc_key), ch), continue; end
    ch_data = d_vq.RPT_VQ_grid.(cyc_key).(ch);

    for r = 1:length(crates)
        f = [crates{r} '_charge'];
        if ~isfield(ch_data, f), continue; end
        s = ch_data.(f);
        t_s = seconds(s.t_raw(:)); V_raw = double(s.V_raw(:));
        I_raw = double(s.I_raw(:)); Q_raw = s.Q_raw(:);

        % ECM steady-state correction
        V_ocv = ecm_ss(V_raw, I_raw, t_s, ecm_chg, Q_0, true);
        V_ocv_sm = movmean(V_ocv, sm_win);

        % Build V-grid from corrected voltage
        [Vu, u] = unique(V_ocv_sm, 'stable'); Qu = Q_raw(u);
        mo = true(size(Vu));
        for ii=2:length(Vu), if Vu(ii)<=Vu(ii-1), mo(ii)=false; end; end
        Vu=Vu(mo); Qu=Qu(mo);
        Vg=(min(Vu):0.001:max(Vu))'; Qg=interp1(Vu,Qu,Vg,'linear');
        ok=~isnan(Qg); Vg=Vg(ok); Qg=Qg(ok);
        dqv = movmean(gradient(Qg)./gradient(Vg), ma_win);

        mask = Vg >= 3.3 & Vg <= 4.15;
        plot3(Vg(mask), cyc_num*ones(sum(mask),1), dqv(mask), ...
            '-','Color',colors_cr(r,:),'LineWidth',1.2);
    end
end
for r = 1:5
    plot3(NaN,NaN,NaN,'-','Color',colors_cr(r,:),'LineWidth',2,'DisplayName',clabels{r});
end
legend('Location','northeast','FontSize',9);
saveas(fig2, fullfile(saveDir,'ECM_3D_dQdV_After.fig'));

fprintf('Done. 2 figs saved.\n');

%% Helper
function V_ocv = ecm_ss(V_raw, I_raw, t_raw, ecm, Q_0, is_charge)
    I_abs = abs(I_raw);
    Q_cum = cumtrapz(t_raw, I_abs)/3600;
    if is_charge, SOC = Q_cum/Q_0*100;
    else, SOC = 100 - Q_cum/Q_0*100; end
    SOC = max(0,min(100,SOC));
    [soc_s,si] = sort(ecm.SOC);
    R_tot_s = (ecm.R0(si)+ecm.R1(si)+ecm.R2(si))/1e3;
    valid = soc_s >= 10 & soc_s <= 95;
    soc_s = soc_s(valid); R_tot_s = R_tot_s(valid);
    SOC_c = max(min(soc_s),min(max(soc_s),SOC));
    R_tot = interp1(soc_s, R_tot_s, SOC_c, 'linear');
    V_ocv = V_raw - I_raw .* R_tot;
end
