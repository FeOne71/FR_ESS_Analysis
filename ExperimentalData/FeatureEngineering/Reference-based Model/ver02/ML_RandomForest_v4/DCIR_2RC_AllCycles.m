% DCIR_2RC_AllCycles.m
% =====================================================================
% 전 사이클(cyc0~cyc1000)에 대해 채널별 2RC ECM 파라미터 피팅
% 기존 DCIR_2RC_Fitting.m을 전 사이클로 확장
% - pdata 기반 펄스 탐색 + lsqcurvefit 2RC 피팅
% - 결과: ECM_2RC_AllChannels_AllCycles.mat
% =====================================================================
clear; clc; close all;

parsedDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';
dcirFile  = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_Final\DCIR_SOC_data_all_channels_final.mat';
saveDir   = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4';

d_dcir = load(dcirFile);
channels = {'Ch09','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16'};

% Parsed 파일명의 사이클 키 매핑
% DCIR_Final 사이클 키 → Parsed 파일명 숫자
cyc_map = {
    'cyc0',   '0',    0;
    'cyc200',  '200',  200;
    'cyc400',  '400',  400;
    'cyc600',  '600',  600;
    'x800cyc', '800',  800;
};

I_1C = 64; tol = 0.05;
opts = optimoptions('lsqcurvefit','Display','off','MaxIterations',2000,...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12);
pulse_model = @(p,t) p(1) + p(2)*(1-exp(-t/p(3))) + p(4)*(1-exp(-t/p(5)));

All_ECM = struct();

for ch_i = 1:length(channels)
    ch = channels{ch_i};
    ch_lower = lower(ch(3:end));
    fprintf('\n============================\n=== %s ===\n============================\n', ch);

    for ci = 1:size(cyc_map,1)
        dcir_key = cyc_map{ci,1};   % e.g. 'cyc0'
        file_num = cyc_map{ci,2};   % e.g. '0'
        cyc_num  = cyc_map{ci,3};   % e.g. 0

        filename = sprintf('RPT%s_%s_parsed.mat', file_num, ['ch' ch_lower]);
        filepath = fullfile(parsedDir, filename);
        if ~exist(filepath,'file')
            fprintf('  [Skip] %s not found\n', filename);
            continue;
        end

        % SOC reference
        soc_dch = []; soc_chg = [];
        if isfield(d_dcir.dcir_soc_data, ch) && isfield(d_dcir.dcir_soc_data.(ch), dcir_key)
            cyc_data = d_dcir.dcir_soc_data.(ch).(dcir_key);
            if isfield(cyc_data,'discharge_table') && isfield(cyc_data.discharge_table,'SOC')
                soc_dch = cyc_data.discharge_table.SOC;
            end
            if isfield(cyc_data,'charge_table') && isfield(cyc_data.charge_table,'SOC')
                soc_chg = cyc_data.charge_table.SOC;
            end
        end

        dd = load(filepath);
        pdata = dd.pdata;
        fprintf('  %s (cyc%d): %d segments\n', filename, cyc_num, length(pdata));

        % Find discharge 1C pulses
        dch_pulse_idx = [];
        for k = 1:length(pdata)
            avg_I = mean(pdata(k).I); dur = pdata(k).t(end)-pdata(k).t(1);
            if avg_I < -(1-tol)*I_1C && avg_I > -(1+tol)*I_1C && dur > 5 && dur < 70
                dch_pulse_idx(end+1) = k; %#ok<SAGROW>
            end
        end

        % Find charge 1C pulses
        chg_pulse_idx = [];
        for k = 1:length(pdata)
            avg_I = mean(pdata(k).I); dur = pdata(k).t(end)-pdata(k).t(1);
            if avg_I > (1-tol)*I_1C && avg_I < (1+tol)*I_1C && dur > 5 && dur < 70
                chg_pulse_idx(end+1) = k; %#ok<SAGROW>
            end
        end
        fprintf('    DCH pulses: %d, CHG pulses: %d\n', length(dch_pulse_idx), length(chg_pulse_idx));

        % Fit
        ECM_dch = fit_pulses(pdata, dch_pulse_idx, soc_dch, false, pulse_model, opts);
        ECM_chg = fit_pulses(pdata, chg_pulse_idx, soc_chg, true,  pulse_model, opts);

        % Store
        cyc_key_out = sprintf('cyc%d', cyc_num);
        All_ECM.(ch).(cyc_key_out).discharge = ECM_dch;
        All_ECM.(ch).(cyc_key_out).charge    = ECM_chg;

        % Print summary
        if ~isempty(ECM_dch.SOC)
            fprintf('    DCH R_total(mΩ): ');
            fprintf('%.2f ', ECM_dch.R_total); fprintf('\n');
        end
        if ~isempty(ECM_chg.SOC)
            fprintf('    CHG R_total(mΩ): ');
            fprintf('%.2f ', ECM_chg.R_total); fprintf('\n');
        end
    end

    % cyc1000 — DCIR_Final에 없을 수 있음, Parsed만으로 처리
    filename1000 = sprintf('RPT1000_%s_parsed.mat', ['ch' ch_lower]);
    filepath1000 = fullfile(parsedDir, filename1000);
    if exist(filepath1000,'file')
        dd = load(filepath1000);
        pdata = dd.pdata;
        dch_pulse_idx = [];
        for k = 1:length(pdata)
            avg_I = mean(pdata(k).I); dur = pdata(k).t(end)-pdata(k).t(1);
            if avg_I < -(1-tol)*I_1C && avg_I > -(1+tol)*I_1C && dur > 5 && dur < 70
                dch_pulse_idx(end+1) = k; %#ok<SAGROW>
            end
        end
        chg_pulse_idx = [];
        for k = 1:length(pdata)
            avg_I = mean(pdata(k).I); dur = pdata(k).t(end)-pdata(k).t(1);
            if avg_I > (1-tol)*I_1C && avg_I < (1+tol)*I_1C && dur > 5 && dur < 70
                chg_pulse_idx(end+1) = k; %#ok<SAGROW>
            end
        end
        fprintf('  cyc1000: DCH=%d, CHG=%d pulses\n', length(dch_pulse_idx), length(chg_pulse_idx));
        % SOC 없으면 균등 배정
        ECM_dch = fit_pulses(pdata, dch_pulse_idx, [], false, pulse_model, opts);
        ECM_chg = fit_pulses(pdata, chg_pulse_idx, [], true,  pulse_model, opts);
        All_ECM.(ch).cyc1000.discharge = ECM_dch;
        All_ECM.(ch).cyc1000.charge    = ECM_chg;
    end
end

%% Save
savePath = fullfile(saveDir, 'ECM_2RC_AllChannels_AllCycles.mat');
save(savePath, 'All_ECM');
fprintf('\n>> Saved: ECM_2RC_AllChannels_AllCycles.mat\n');

%% Print R_total 변화 요약 (Ch09 방전)
fprintf('\n=== R_total 열화 추이 (Ch09 방전, SOC 평균) ===\n');
cyc_keys_out = fieldnames(All_ECM.Ch09);
for ci = 1:length(cyc_keys_out)
    ck = cyc_keys_out{ci};
    Rt = All_ECM.Ch09.(ck).discharge.R_total;
    if ~isempty(Rt)
        fprintf('  %6s: R_total avg = %.3f mΩ\n', ck, mean(Rt,'omitnan'));
    end
end

%% ========================================================================
function ECM = fit_pulses(pdata, pulse_idx, soc_ref, is_charge, model, opts)
    n = length(pulse_idx);
    ECM.SOC     = nan(n,1); ECM.R0  = nan(n,1); ECM.R1   = nan(n,1);
    ECM.tau1    = nan(n,1); ECM.R2  = nan(n,1); ECM.tau2 = nan(n,1);
    ECM.R_total = nan(n,1); ECM.V_OCV = nan(n,1);
    lb = [0, 0, 0.1, 0, 1]; ub = [0.01, 0.01, 20, 0.05, 300];

    for p = 1:n
        k = pulse_idx(p);
        t_p   = pdata(k).t(:) - pdata(k).t(1);
        V_p   = pdata(k).V(:);
        I_avg = abs(mean(pdata(k).I));
        if k>1 && strcmp(char(pdata(k-1).type),'R')
            V_OCV = pdata(k-1).V(end);
        else
            V_OCV = V_p(1);
        end
        ECM.V_OCV(p) = V_OCV;
        R_data  = abs(V_p - V_OCV) / I_avg;
        R0_init = R_data(1); R_inf = R_data(end);
        x0 = [R0_init, (R_inf-R0_init)*0.4, 3, (R_inf-R0_init)*0.6, 30];
        try
            [pf,~] = lsqcurvefit(model, x0, t_p, R_data, lb, ub, opts);
            ECM.R0(p)    = pf(1)*1e3; ECM.R1(p)   = pf(2)*1e3;
            ECM.tau1(p)  = pf(3);     ECM.R2(p)   = pf(4)*1e3;
            ECM.tau2(p)  = pf(5);
            ECM.R_total(p) = (pf(1)+pf(2)+pf(4))*1e3;
        catch
            fprintf('  Pulse %d FIT FAILED\n', p);
        end
        if p <= length(soc_ref)
            ECM.SOC(p) = soc_ref(p);
        else
            if is_charge, ECM.SOC(p) = (p-1)/max(n-1,1)*100;
            else,         ECM.SOC(p) = 100 - (p-1)/max(n-1,1)*100; end
        end
    end
end
