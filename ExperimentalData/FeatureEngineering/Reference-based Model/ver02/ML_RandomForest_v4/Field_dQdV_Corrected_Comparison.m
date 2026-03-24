% Field_dQdV_Corrected_Comparison.m
% =====================================================================
% 비교:
%   1) Raw: 스무딩 전혀 없음 (원본 V,Q → dQ/dV)
%   2) 기존 (잘못됨): V,Q 각각 movmean(60) → dQ/dV
%   3) 수정 (V-grid): 원본 V,Q → 단조 정리 → V_grid(0.001V) interp → dQ/dV → dQ/dV movmean(21)
%
% 연도별 1×2 (충전, 방전), 3 methods 오버레이
% =====================================================================
clear; clc; close all;

rawDir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'Rack_raw2mat');
saveDir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
    'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02', 'ML_RandomForest_v4');

dataFiles = {
    fullfile(rawDir, 'Old', '2021', '202106', 'Raw_20210603.mat'), 'Y2021', 'old',  datetime(2021,6,3), 600, 300;
    fullfile(rawDir, 'New', '2023', '202310', 'Raw_20231016.mat'), 'Y2023', 'new',  datetime(2023,10,16), 300, 150;
    fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024', 'new',  datetime(2024,9,9), 300, 300;
    fullfile(rawDir, 'New', '2025', '202507', 'Raw_20250711.mat'), 'Y2025', 'new',  datetime(2024,9,9), 300, 300;
};

Np = 2; C_cell_Ah = 64; thr_A = C_cell_Ah * 0.05;

for k = 1:size(dataFiles, 1)
    fpath = dataFiles{k,1}; yr = dataFiles{k,2}; dt = dataFiles{k,3}; bd = dataFiles{k,4};
    mcs = dataFiles{k,5}; mds = dataFiles{k,6};
    if ~exist(fpath,'file'), continue; end
    fprintf('Processing %s...\n', yr);

    S = load(fpath);
    if strcmp(dt,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'), t=datetime(D.Time);
        else, if isduration(D.Date_Time), t=bd+D.Date_Time; else, t=datetime(D.Date_Time); end, end
        Ir=D.DCCurrent_A(:); Va=D.AverageCV_V(:);
    else
        D=S.Raw;
        if isduration(D.Date_Time), t=bd+D.Date_Time; else, t=datetime(D.Date_Time); end
        Ir=D.DCCurrent(:); Va=D.CVavg(:);
    end
    ts = seconds(t-t(1)); Ic = Ir/Np;

    % --- Segment extraction ---
    chg_mask = Ic > thr_A; dch_mask = Ic < -thr_A;
    segs = struct('label',{'Charge','Discharge'}, 'mask',{chg_mask,dch_mask}, ...
                  'min_s',{mcs,mds}, 'asc',{true,false});

    fig = figure('Name', sprintf('dQdV Corrected - %s', yr), ...
                 'Position', [50, 100, 1400, 500], 'Visible', 'off');

    for si = 1:2
        dd = diff([0; segs(si).mask(:); 0]);
        ss = find(dd==1); ee = find(dd==-1)-1;
        if isempty(ss), continue; end
        lens = ee-ss+1; vi = find(lens >= segs(si).min_s);
        if isempty(vi), continue; end
        [~,mx] = max(lens(vi)); ii = vi(mx);
        s_i = ss(ii); e_i = ee(ii);

        V_raw = Va(s_i:e_i);
        I_raw = Ic(s_i:e_i);
        t_seg = ts(s_i:e_i);
        Q_raw = cumtrapz(t_seg, abs(I_raw)) / 3600;
        is_asc = segs(si).asc;

        subplot(1, 2, si);

        %% ============ Method 1: Raw (no smoothing) ============
        % 단조 정리만 하고 gradient
        [Vu1, Qu1] = make_monotonic(V_raw, Q_raw, is_asc);
        if ~is_asc, Vu1=flipud(Vu1); Qu1=flipud(Qu1); end
        dV1 = gradient(Vu1); dQ1 = gradient(Qu1);
        dV1(dV1==0) = NaN;
        dQdV1 = dQ1 ./ dV1;
        dQdV1(isinf(dQdV1)|isnan(dQdV1)) = 0;
        plot(Vu1, dQdV1, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.4, 'DisplayName', '1) Raw');
        hold on;

        %% ============ Method 2: 기존 (V,Q 각각 movmean60) ============
        V_sm = movmean(V_raw, 60);
        Q_sm = movmean(Q_raw, 60);
        [Vu2, Qu2] = make_monotonic(V_sm, Q_sm, is_asc);
        if ~is_asc, Vu2=flipud(Vu2); Qu2=flipud(Qu2); end
        dV2 = gradient(Vu2); dQ2 = gradient(Qu2);
        dV2(dV2==0) = NaN;
        dQdV2 = dQ2 ./ dV2;
        dQdV2(isinf(dQdV2)|isnan(dQdV2)) = 0;
        plot(Vu2, dQdV2, 'b-', 'LineWidth', 0.8, 'DisplayName', '2) V,Q ma60 (기존-잘못됨)');

        %% ============ Method 3: V-grid 방식 (수정) ============
        % Step 1: 원본 V,Q에서 단조 정리
        [Vm, Qm] = make_monotonic(V_raw, Q_raw, is_asc);
        if is_asc
            v_lo = min(Vm); v_hi = max(Vm);
            V_grid = (v_lo : 0.001 : v_hi)';
        else
            v_hi = max(Vm); v_lo = min(Vm);
            Vm = flipud(Vm); Qm = flipud(Qm); % ascending for interp1
            V_grid = (v_lo : 0.001 : v_hi)';
        end
        % Step 2: Q를 V_grid에 보간
        Q_grid = interp1(Vm, Qm, V_grid, 'linear');
        Q_grid(isnan(Q_grid)) = [];
        V_grid = V_grid(1:length(Q_grid));
        % Step 3: dQ/dV 계산
        dV3 = gradient(V_grid);
        dQ3 = gradient(Q_grid);
        dV3(dV3==0) = NaN;
        dQdV3 = dQ3 ./ dV3;
        dQdV3(isinf(dQdV3)|isnan(dQdV3)) = 0;
        % Step 4: dQ/dV에만 movmean(21)
        dQdV3_sm = movmean(dQdV3, 21);

        plot(V_grid, dQdV3_sm, 'r-', 'LineWidth', 1.5, 'DisplayName', '3) V-grid + dQdV ma21 (수정)');

        % Y-axis auto
        dc = dQdV3_sm(~isnan(dQdV3_sm) & ~isinf(dQdV3_sm));
        if ~isempty(dc)
            q2 = prctile(dc,2); q98 = prctile(dc,98);
            mg = (q98-q2)*0.2; if mg<0.5, mg=0.5; end
            ylim([q2-mg, q98+mg]);
        end

        title(sprintf('%s — %s', segs(si).label, yr), 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('V'); ylabel('dQ/dV (Ah/V)');
        legend('Location', 'best', 'FontSize', 9); grid on;
    end

    sgtitle(sprintf('[%s] dQ/dV: Raw vs 기존(잘못) vs V-grid(수정)', yr), ...
            'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig, fullfile(saveDir, sprintf('dQdV_Corrected_%s.fig', yr)));
    close(fig);
    fprintf('  Saved: dQdV_Corrected_%s.fig\n', yr);
end

fprintf('\nDone.\n');

%% Helper
function [Vm, Qm] = make_monotonic(V_in, Q_in, is_ascending)
    [Vu, uid] = unique(V_in(:), 'stable');
    Qu = Q_in(uid);
    mono = true(size(Vu));
    if is_ascending
        for ii = 2:length(Vu)
            if Vu(ii) <= Vu(ii-1), mono(ii) = false; end
        end
    else
        for ii = 2:length(Vu)
            if Vu(ii) >= Vu(ii-1), mono(ii) = false; end
        end
    end
    Vm = Vu(mono); Qm = Qu(mono);
end
