% Plot_dQdV_AllSegments_4panel.m
% 4패널: Charge/Discharge × Aging/C-rate 효과 | 전체 세그먼트 경계 표시

clear; clc; close all;

%% Paths
baseDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData';
verDir  = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
vqFile  = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
visDir  = fullfile(verDir, 'RatioModel', 'Visualization');

g = load(vqFile, 'RPT_VQ_grid');
G = g.RPT_VQ_grid;
mr = load(fullfile(verDir, 'MasterRulers_PseudoOCV.mat'));
Vb_c = mr.MasterRuler_ver0317.V_bounds_chg;  % 13개 경계 → 12 세그먼트
Vb_d = mr.MasterRuler_ver0317.V_bounds_dch;
cell_id = 'Ch09';

%% Settings
% Aging: OCV 조건으로 사이클별 비교
sel_cycs   = {'cyc0','cyc200','cyc400','cyc800'};
cyc_labels = {'Cyc 0','Cyc 200','Cyc 400','Cyc 800'};
cyc_colors = [0.1 0.1 0.85; 0.05 0.70 0.70; 0.85 0.60 0.0; 0.85 0.05 0.05];

% C-rate: cyc0에서 조건별 비교
conds    = {'OCV_charge','c01_charge','c05_charge','c1_charge','c2_charge','c3_charge'};
conds_d  = {'OCV_discharge','c01_discharge','c05_discharge','c1_discharge','c2_discharge','c3_discharge'};
cr_labels = {'0.05C(OCV)','0.1C','0.5C','1C','2C','3C'};
cr_colors = [0.0 0.0 0.9; 0.0 0.75 1.0; 0.0 0.82 0.0; 0.85 0.85 0.0; 1.0 0.5 0.0; 0.9 0.0 0.0];

smfn = @(V,Q) deal(smoothdata(diff(Q)./diff(V),'sgolay',51), (V(1:end-1)+V(2:end))/2);

%% Figure
fig = figure('Position',[40 30 1500 920]);

% ===== Panel 1: Charge Aging =====
ax1 = subplot(2,2,1); hold(ax1,'on');
for ki = 1:length(sel_cycs)
    cs = sel_cycs{ki};
    if isfield(G,cs) && isfield(G.(cs),cell_id) && isfield(G.(cs).(cell_id),'OCV_charge')
        d = G.(cs).(cell_id).OCV_charge;
        [dqdv,vm] = smfn(d.V_raw, d.Q_raw);
        plot(ax1, vm, abs(dqdv), '-', 'Color', cyc_colors(ki,:), 'LineWidth', 1.8, 'DisplayName', cyc_labels{ki});
    end
end
add_seg_labels(ax1, Vb_c, 280);
title(ax1,'Charge |dQ/dV|: Aging (OCV, Cyc 0~800)','FontSize',11,'FontWeight','bold');
xlabel(ax1,'Voltage (V)','FontSize',10); ylabel(ax1,'|dQ/dV| (Ah/V)','FontSize',10);
xlim(ax1,[3.0 4.2]); ylim(ax1,[0 290]); legend(ax1,'Location','northwest','FontSize',9);

% ===== Panel 2: Discharge Aging =====
ax2 = subplot(2,2,2); hold(ax2,'on');
for ki = 1:length(sel_cycs)
    cs = sel_cycs{ki};
    if isfield(G,cs) && isfield(G.(cs),cell_id) && isfield(G.(cs).(cell_id),'OCV_discharge')
        d = G.(cs).(cell_id).OCV_discharge;
        [dqdv,vm] = smfn(d.V_raw, d.Q_raw);
        plot(ax2, vm, abs(dqdv), '-', 'Color', cyc_colors(ki,:), 'LineWidth', 1.8, 'DisplayName', cyc_labels{ki});
    end
end
add_seg_labels(ax2, Vb_d, 220);
title(ax2,'Discharge |dQ/dV|: Aging (OCV, Cyc 0~800)','FontSize',11,'FontWeight','bold');
xlabel(ax2,'Voltage (V)','FontSize',10); ylabel(ax2,'|dQ/dV| (Ah/V)','FontSize',10);
xlim(ax2,[3.0 4.2]); ylim(ax2,[0 230]); legend(ax2,'Location','northwest','FontSize',9);

% ===== Panel 3: Charge C-rate =====
ax3 = subplot(2,2,3); hold(ax3,'on');
for ki = 1:length(conds)
    cf = conds{ki};
    if isfield(G,'cyc0') && isfield(G.cyc0,cell_id) && isfield(G.cyc0.(cell_id),cf)
        d = G.cyc0.(cell_id).(cf);
        [dqdv,vm] = smfn(d.V_raw, d.Q_raw);
        plot(ax3, vm, abs(dqdv), '-', 'Color', cr_colors(ki,:), 'LineWidth', 1.8, 'DisplayName', cr_labels{ki});
    end
end
add_seg_labels(ax3, Vb_c, 260);
title(ax3,'Charge |dQ/dV|: C-rate Effect (Cyc 0)','FontSize',11,'FontWeight','bold');
xlabel(ax3,'Voltage (V)','FontSize',10); ylabel(ax3,'|dQ/dV| (Ah/V)','FontSize',10);
xlim(ax3,[3.0 4.2]); ylim(ax3,[0 270]); legend(ax3,'Location','northwest','FontSize',9);

% ===== Panel 4: Discharge C-rate =====
ax4 = subplot(2,2,4); hold(ax4,'on');
for ki = 1:length(conds_d)
    df = conds_d{ki};
    if isfield(G,'cyc0') && isfield(G.cyc0,cell_id) && isfield(G.cyc0.(cell_id),df)
        d = G.cyc0.(cell_id).(df);
        [dqdv,vm] = smfn(d.V_raw, d.Q_raw);
        plot(ax4, vm, abs(dqdv), '-', 'Color', cr_colors(ki,:), 'LineWidth', 1.8, 'DisplayName', cr_labels{ki});
    end
end
add_seg_labels(ax4, Vb_d, 215);
title(ax4,'Discharge |dQ/dV|: C-rate Effect (Cyc 0)','FontSize',11,'FontWeight','bold');
xlabel(ax4,'Voltage (V)','FontSize',10); ylabel(ax4,'|dQ/dV| (Ah/V)','FontSize',10);
xlim(ax4,[3.0 4.2]); ylim(ax4,[0 225]); legend(ax4,'Location','northwest','FontSize',9);

sgtitle(sprintf('%s | MasterRuler 전체 세그먼트 경계 (Seg01~Seg12)', cell_id), 'FontSize',13,'FontWeight','bold');

outfile = fullfile(visDir, 'dQdV_AllSegments_4panel.png');
saveas(fig, outfile);
fprintf('Saved: %s\n', outfile);

%% Helper: 세그먼트 경계 + 번호 표시
function add_seg_labels(ax, Vb, ymax)
    for bi = 1:length(Vb)
        xline(ax, Vb(bi), '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.9, 'Alpha',0.7, 'HandleVisibility','off');
    end
    for si = 1:length(Vb)-1
        Vmid = (Vb(si) + Vb(si+1)) / 2;
        text(ax, Vmid, ymax*0.96, sprintf('S%02d',si), ...
            'FontSize',7, 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
            'Color',[0.25 0.25 0.25], 'FontWeight','bold');
    end
end
