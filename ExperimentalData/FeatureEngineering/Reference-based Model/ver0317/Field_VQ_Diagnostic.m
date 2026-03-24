% Field_VQ_Diagnostic.m
% 필드 연도별 충전/방전 V-Q 곡선 + Master Ruler 세그먼트 경계 시각화
% 목적: 필드 데이터의 전압 커버리지와 dQ 피처 매핑을 시각적으로 확인

clear; clc; close all;

%% 1. Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
rawDir  = fullfile(projDir, 'Rack_raw2mat');
visDir  = fullfile(verDir, 'Visualization');
if ~exist(visDir,'dir'), mkdir(visDir); end

%% 2. Load Master Ruler
mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
V_chg = mr.MasterRuler_ver0317.V_bounds_chg;
V_dch = mr.MasterRuler_ver0317.V_bounds_dch;

%% 3. Field Data
dataFiles = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old', datetime(2021,6,3);
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new', datetime(2023,10,16);
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new', datetime(2024,9,9);
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new', datetime(2025,7,11);
};

Np = 2; Q_nom = 64; thr_A = Q_nom * 0.05 / Np;
min_chg_sec = [600, 300, 300, 300];
min_dch_sec = [300, 150, 300, 150];
ma_window = 60;

yr_colors = [0.2 0.4 0.85; 0.9 0.35 0.1; 0.1 0.7 0.3; 0.6 0.2 0.8];

%% 4. Process and Plot
fig1 = figure('Position', [50 50 1400 900]);
sgtitle('Field Charge V-Q Curves + Master Ruler Segments', 'FontSize', 16, 'FontWeight', 'bold');

fig2 = figure('Position', [100 50 1400 900]);
sgtitle('Field Discharge V-Q Curves + Master Ruler Segments', 'FontSize', 16, 'FontWeight', 'bold');

for k = 1:size(dataFiles, 1)
    fpath=dataFiles{k,1}; yr=dataFiles{k,2}; dtype=dataFiles{k,3}; base_date=dataFiles{k,4};
    if ~exist(fpath,'file'), continue; end
    
    S = load(fpath);
    if strcmp(dtype,'old')
        D=S.Raw.Rack01;
        if isfield(D,'Time'), t=datetime(D.Time);
        elseif isduration(D.Date_Time), t=base_date+D.Date_Time;
        else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent_A(:); V_avg=D.AverageCV_V(:);
    else
        D=S.Raw;
        if isduration(D.Date_Time), t=base_date+D.Date_Time; else, t=datetime(D.Date_Time); end
        I_rack=D.DCCurrent(:); V_avg=D.CVavg(:);
    end
    tsec=seconds(t-t(1)); I_cell=I_rack/Np;
    
    chgSegs=local_find_segments(I_cell>thr_A);
    dchSegs=local_find_segments(I_cell<-thr_A);
    if ~isempty(chgSegs), dur=chgSegs(:,2)-chgSegs(:,1)+1; chgSegs=chgSegs(dur>=min_chg_sec(k),:); end
    if ~isempty(dchSegs), dur=dchSegs(:,2)-dchSegs(:,1)+1; dchSegs=dchSegs(dur>=min_dch_sec(k),:); end
    
    %% --- CHARGE ---
    if ~isempty(chgSegs)
        [~,bi]=max(chgSegs(:,2)-chgSegs(:,1));
        cs=chgSegs(bi,1); ce=chgSegs(bi,2);
        v_c=V_avg(cs:ce); i_c=I_cell(cs:ce); t_c=tsec(cs:ce);
        
        % Raw V-Q
        q_c_raw = cumtrapz(t_c, abs(i_c)) / 3600;
        
        % Smoothed V-Q
        i_c_sm=movmean(abs(i_c),ma_window);
        q_c_sm=cumtrapz(t_c,i_c_sm)/3600;
        v_c_sm=movmean(v_c,ma_window);
        
        figure(fig1);
        subplot(2,2,k); hold on; grid on; box on;
        
        % Plot raw and smoothed
        plot(q_c_raw, v_c, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'DisplayName', 'Raw');
        plot(q_c_sm, v_c_sm, '-', 'Color', yr_colors(k,:), 'LineWidth', 2, 'DisplayName', 'Smoothed');
        
        % Segment boundaries
        for b = 1:length(V_chg)
            yline(V_chg(b), '--', sprintf('%.3f', V_chg(b)), 'Color', [0.5 0 0], ...
                'LineWidth', 0.8, 'FontSize', 7, 'LabelHorizontalAlignment', 'right');
        end
        
        % Shade valid segments
        xl = xlim;
        for s = 1:12
            vlo = V_chg(s); vhi = V_chg(s+1);
            if vlo >= min(v_c)-0.02 && vhi <= max(v_c)+0.02
                patch([xl(1) xl(2) xl(2) xl(1)], [vlo vlo vhi vhi], ...
                    [0.2 0.8 0.2], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end
        
        xlabel('Capacity (Ah)'); ylabel('Voltage (V)');
        title(sprintf('%s Charge [%.3f~%.3fV] C_{eff}=%.2f', yr, min(v_c), max(v_c), mean(abs(i_c))/Q_nom), 'FontSize', 11);
        legend('Location', 'southeast');
        ylim([min(V_chg)-0.05, max(V_chg)+0.05]);
    end
    
    %% --- DISCHARGE ---
    if ~isempty(dchSegs)
        [~,bi]=max(dchSegs(:,2)-dchSegs(:,1));
        ds=dchSegs(bi,1); de=dchSegs(bi,2);
        v_d=V_avg(ds:de); i_d=I_cell(ds:de); t_d=tsec(ds:de);
        
        % Raw V-Q
        q_d_raw = cumtrapz(t_d, abs(i_d)) / 3600;
        
        % Smoothed V-Q
        i_d_sm=movmean(abs(i_d),ma_window);
        q_d_sm=cumtrapz(t_d,i_d_sm)/3600;
        v_d_sm=movmean(v_d,ma_window);
        
        figure(fig2);
        subplot(2,2,k); hold on; grid on; box on;
        
        plot(q_d_raw, v_d, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'DisplayName', 'Raw');
        plot(q_d_sm, v_d_sm, '-', 'Color', yr_colors(k,:), 'LineWidth', 2, 'DisplayName', 'Smoothed');
        
        for b = 1:length(V_dch)
            yline(V_dch(b), '--', sprintf('%.3f', V_dch(b)), 'Color', [0 0 0.5], ...
                'LineWidth', 0.8, 'FontSize', 7, 'LabelHorizontalAlignment', 'right');
        end
        
        xl = xlim;
        for s = 1:12
            vlo = min(V_dch(s), V_dch(s+1)); vhi = max(V_dch(s), V_dch(s+1));
            if vlo >= min(v_d)-0.02 && vhi <= max(v_d)+0.02
                patch([xl(1) xl(2) xl(2) xl(1)], [vlo vlo vhi vhi], ...
                    [0.2 0.2 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end
        
        xlabel('Capacity (Ah)'); ylabel('Voltage (V)');
        title(sprintf('%s Discharge [%.3f~%.3fV] C_{eff}=%.2f', yr, min(v_d), max(v_d), mean(abs(i_d))/Q_nom), 'FontSize', 11);
        legend('Location', 'northeast');
        ylim([min(V_dch)-0.05, max(V_dch)+0.05]);
    end
end

%% 5. Save
figure(fig1);
saveas(fig1, fullfile(visDir, 'Field_VQ_Charge_Diagnostic.png'));
fprintf('Saved: Field_VQ_Charge_Diagnostic.png\n');

figure(fig2);
saveas(fig2, fullfile(visDir, 'Field_VQ_Discharge_Diagnostic.png'));
fprintf('Saved: Field_VQ_Discharge_Diagnostic.png\n');

fprintf('=== Diagnostic Complete! ===\n');

%% Helper
function segs = local_find_segments(mask)
    segs=[]; n=length(mask); i=1;
    while i<=n
        if mask(i), j=i;
            while j<n && mask(j+1), j=j+1; end
            segs=[segs;i,j]; i=j+1;
        else, i=i+1; end
    end
end
