% DCIR_2RC_Fitting.m
% =====================================================================
% 모든 채널 (Ch09~Ch16) DCIR 2RC ECM 파라미터 피팅
% - 방전 1C 펄스 + 충전 1C 펄스 각각 별도 피팅
% - 시각화: 실측 전압 (검정 o) + 모델 전압 (빨간 o, 펄스 구간)
% - extrap 제거 → 경계 클리핑으로 교체
% =====================================================================
clear; clc; close all;

parsedDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';
dcirFile  = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_Final\DCIR_SOC_data_all_channels_final.mat';
saveDir   = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4';

d_dcir = load(dcirFile);
channels = {'Ch09','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16'};
cycle = '0';
I_1C = 64; tol = 0.05;

opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 2000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
pulse_model = @(p, t) p(1) + p(2)*(1-exp(-t/p(3))) + p(4)*(1-exp(-t/p(5)));

All_ECM = struct();

for ch_i = 1:length(channels)
    ch = channels{ch_i};
    ch_lower = lower(ch(3:end));
    filename = sprintf('RPT%s_%s_parsed.mat', cycle, ['ch' ch_lower]);
    filepath = fullfile(parsedDir, filename);
    if ~exist(filepath,'file'), fprintf('Skip: %s\n', filename); continue; end

    fprintf('\n=== %s ===\n', ch);
    dd = load(filepath);
    pdata = dd.pdata;

    % SOC references from DCIR_Final
    if isfield(d_dcir.dcir_soc_data, ch) && isfield(d_dcir.dcir_soc_data.(ch), 'cyc0')
        soc_dch = d_dcir.dcir_soc_data.(ch).cyc0.discharge_table.SOC;
        soc_chg = d_dcir.dcir_soc_data.(ch).cyc0.charge_table.SOC;
    else
        soc_dch = []; soc_chg = [];
    end

    %% ---- Find discharge 1C pulses ----
    dch_pulse_idx = [];
    for k = 1:length(pdata)
        avg_I = mean(pdata(k).I); dur = pdata(k).t(end)-pdata(k).t(1);
        if avg_I < -(1-tol)*I_1C && avg_I > -(1+tol)*I_1C && dur > 5 && dur < 70
            dch_pulse_idx(end+1) = k; %#ok<SAGROW>
        end
    end

    %% ---- Find charge 1C pulses ----
    chg_pulse_idx = [];
    for k = 1:length(pdata)
        avg_I = mean(pdata(k).I); dur = pdata(k).t(end)-pdata(k).t(1);
        if avg_I > (1-tol)*I_1C && avg_I < (1+tol)*I_1C && dur > 5 && dur < 70
            chg_pulse_idx(end+1) = k; %#ok<SAGROW>
        end
    end
    fprintf('  Discharge pulses: %d,  Charge pulses: %d\n', length(dch_pulse_idx), length(chg_pulse_idx));

    %% ---- Fit both ----
    ECM_dch = fit_pulses(pdata, dch_pulse_idx, soc_dch, false, pulse_model, opts);
    ECM_chg = fit_pulses(pdata, chg_pulse_idx, soc_chg, true,  pulse_model, opts);

    %% ---- Visualization (discharge) ----
    n_dch = length(dch_pulse_idx);
    n_chg = length(chg_pulse_idx);
    n_cols = 5;

    % Discharge figure
    if n_dch > 0
        fig = figure('Visible','off','Name',sprintf('DCH %s',ch), ...
            'Position',[30,30,300*n_cols,250*ceil(n_dch/n_cols)]);
        for p = 1:n_dch
            k = dch_pulse_idx(p);
            t_p = pdata(k).t(:)-pdata(k).t(1);
            V_p = pdata(k).V(:);
            V_OCV = pdata(max(k-1,1)).V(end);
            I_avg = abs(mean(pdata(k).I));

            subplot(ceil(n_dch/n_cols), n_cols, p); hold on;
            % Pre-rest (up to 60s)
            if k>1 && strcmp(char(pdata(k-1).type),'R')
                t_pre = pdata(k-1).t(:)-pdata(k-1).t(1);
                valid = t_pre >= (t_pre(end)-min(60,t_pre(end)));
                t_pre_p = t_pre(valid)-t_pre(find(valid,1));
                t_pre_p = t_pre_p - (t_pre_p(end)+0.1);
                plot(t_pre_p, pdata(k-1).V(valid), 'ko','MarkerSize',2,'DisplayName','Meas');
            end
            % Pulse: measured + model
            pf = [ECM_dch.R0(p), ECM_dch.R1(p), ECM_dch.tau1(p), ECM_dch.R2(p), ECM_dch.tau2(p)] / [1e3 1e3 1 1e3 1];
            pf = [ECM_dch.R0(p)/1e3, ECM_dch.R1(p)/1e3, ECM_dch.tau1(p), ECM_dch.R2(p)/1e3, ECM_dch.tau2(p)];
            V_model = V_OCV - I_avg * pulse_model(pf, t_p); % discharge: V drops
            plot(t_p, V_p,      'ko','MarkerSize',2,'HandleVisibility','off');
            plot(t_p, V_model,  'ro','MarkerSize',2,'DisplayName','Model');
            % Post-rest (up to 120s)
            if k+1<=length(pdata) && strcmp(char(pdata(k+1).type),'R')
                t_r = pdata(k+1).t(:)-pdata(k+1).t(1);
                valid_r = t_r <= min(120,t_r(end));
                plot(t_p(end)+t_r(valid_r), pdata(k+1).V(valid_r), 'ko','MarkerSize',2,'HandleVisibility','off');
            end
            title(sprintf('SOC=%.0f%%',ECM_dch.SOC(p)),'FontSize',9,'FontWeight','bold');
            xlabel('t(s)'); ylabel('V'); grid on;
            if p==1, legend('FontSize',6,'Location','best'); end
        end
        sgtitle(sprintf('[%s cyc%s] DCH 2RC — Measured(k.o) vs Model(r.o)',ch,cycle),'FontSize',12,'FontWeight','bold');
        saveas(fig, fullfile(saveDir, sprintf('DCIR_2RC_DCH_%s_cyc%s.fig',ch,cycle)));
        close(fig);
    end

    % Charge figure
    if n_chg > 0
        fig = figure('Visible','off','Name',sprintf('CHG %s',ch), ...
            'Position',[30,30,300*n_cols,250*ceil(n_chg/n_cols)]);
        for p = 1:n_chg
            k = chg_pulse_idx(p);
            t_p = pdata(k).t(:)-pdata(k).t(1);
            V_p = pdata(k).V(:);
            V_OCV = pdata(max(k-1,1)).V(end);
            I_avg = abs(mean(pdata(k).I));

            subplot(ceil(n_chg/n_cols), n_cols, p); hold on;
            if k>1 && strcmp(char(pdata(k-1).type),'R')
                t_pre = pdata(k-1).t(:)-pdata(k-1).t(1);
                valid = t_pre >= (t_pre(end)-min(60,t_pre(end)));
                t_pre_p = t_pre(valid)-t_pre(find(valid,1));
                t_pre_p = t_pre_p - (t_pre_p(end)+0.1);
                plot(t_pre_p, pdata(k-1).V(valid), 'ko','MarkerSize',2,'DisplayName','Meas');
            end
            pf = [ECM_chg.R0(p)/1e3, ECM_chg.R1(p)/1e3, ECM_chg.tau1(p), ECM_chg.R2(p)/1e3, ECM_chg.tau2(p)];
            V_model = V_OCV + I_avg * pulse_model(pf, t_p); % charge: V rises
            plot(t_p, V_p,      'ko','MarkerSize',2,'HandleVisibility','off');
            plot(t_p, V_model,  'ro','MarkerSize',2,'DisplayName','Model');
            if k+1<=length(pdata) && strcmp(char(pdata(k+1).type),'R')
                t_r = pdata(k+1).t(:)-pdata(k+1).t(1);
                valid_r = t_r <= min(120,t_r(end));
                plot(t_p(end)+t_r(valid_r), pdata(k+1).V(valid_r), 'ko','MarkerSize',2,'HandleVisibility','off');
            end
            title(sprintf('SOC=%.0f%%',ECM_chg.SOC(p)),'FontSize',9,'FontWeight','bold');
            xlabel('t(s)'); ylabel('V'); grid on;
            if p==1, legend('FontSize',6,'Location','best'); end
        end
        sgtitle(sprintf('[%s cyc%s] CHG 2RC — Measured(k.o) vs Model(r.o)',ch,cycle),'FontSize',12,'FontWeight','bold');
        saveas(fig, fullfile(saveDir, sprintf('DCIR_2RC_CHG_%s_cyc%s.fig',ch,cycle)));
        close(fig);
    end

    All_ECM.(ch).discharge = ECM_dch;
    All_ECM.(ch).charge    = ECM_chg;

    % Print summary
    fprintf('  --- Discharge ---\n');
    fprintf('  %5s | %7s %7s %5s %7s %5s\n','SOC','R0','R1','τ1','R2','τ2');
    for p=1:length(dch_pulse_idx)
        fprintf('  %5.0f%% | %6.3f %6.3f %5.1f %6.3f %5.1f\n', ...
            ECM_dch.SOC(p),ECM_dch.R0(p),ECM_dch.R1(p),ECM_dch.tau1(p),ECM_dch.R2(p),ECM_dch.tau2(p));
    end
    fprintf('  --- Charge ---\n');
    for p=1:length(chg_pulse_idx)
        fprintf('  %5.0f%% | %6.3f %6.3f %5.1f %6.3f %5.1f\n', ...
            ECM_chg.SOC(p),ECM_chg.R0(p),ECM_chg.R1(p),ECM_chg.tau1(p),ECM_chg.R2(p),ECM_chg.tau2(p));
    end
end

%% Save
save(fullfile(saveDir,'ECM_2RC_AllChannels_cyc0.mat'), 'All_ECM');
fprintf('\n>> Saved: ECM_2RC_AllChannels_cyc0.mat\n');

%% ========================================================================
function ECM = fit_pulses(pdata, pulse_idx, soc_ref, is_charge, model, opts)
    n = length(pulse_idx);
    ECM.SOC   = nan(n,1); ECM.R0 = nan(n,1); ECM.R1 = nan(n,1);
    ECM.tau1  = nan(n,1); ECM.R2 = nan(n,1); ECM.tau2 = nan(n,1);
    ECM.R_total = nan(n,1); ECM.V_OCV = nan(n,1);

    lb = [0, 0, 0.1, 0, 1]; ub = [0.01, 0.01, 20, 0.05, 300];

    for p = 1:n
        k = pulse_idx(p);
        t_p = pdata(k).t(:) - pdata(k).t(1);
        V_p = pdata(k).V(:);
        I_avg = abs(mean(pdata(k).I));

        if k>1 && strcmp(char(pdata(k-1).type),'R')
            V_OCV = pdata(k-1).V(end);
        else
            V_OCV = V_p(1);
        end
        ECM.V_OCV(p) = V_OCV;

        % R(t) = |dV(t)| / I: sign convention for charge vs discharge
        R_data = abs(V_p - V_OCV) / I_avg;
        R0_init = R_data(1); R_inf = R_data(end);
        x0 = [R0_init, (R_inf-R0_init)*0.4, 3, (R_inf-R0_init)*0.6, 30];

        try
            [pf,~] = lsqcurvefit(model, x0, t_p, R_data, lb, ub, opts);
            ECM.R0(p) = pf(1)*1e3; ECM.R1(p) = pf(2)*1e3;
            ECM.tau1(p) = pf(3); ECM.R2(p) = pf(4)*1e3;
            ECM.tau2(p) = pf(5); ECM.R_total(p) = (pf(1)+pf(2)+pf(4))*1e3;
        catch
            fprintf('  Pulse %d FIT FAILED\n', p);
        end

        if p <= length(soc_ref)
            ECM.SOC(p) = soc_ref(p);
        else
            if is_charge, ECM.SOC(p) = (p-1)/(max(n-1,1))*100;
            else, ECM.SOC(p) = 100 - (p-1)/(max(n-1,1))*100; end
        end
    end
end
