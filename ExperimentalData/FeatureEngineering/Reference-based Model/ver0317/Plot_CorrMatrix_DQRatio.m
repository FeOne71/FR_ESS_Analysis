% Plot_CorrMatrix_DQRatio.m
% Feature × Feature + SOH Correlation Matrix (Fig.4 style), per C-rate
% English titles

clear; clc; close all;

%% Paths
verDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver0317';
d  = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat')); FM = d.FM;
visDir = fullfile(verDir, 'RatioModel', 'Visualization');
Q_nom  = 64;

%% Columns
chg_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput',false);
dch_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput',false);

%% C-rate groups
uConds    = {'c01','c05','c1','c2','c3'};
cr_labels = {'0.1C','0.5C','1.0C','2.0C','3.0C'};

%% Labels
lbl_c = [arrayfun(@(i)sprintf('Seg%02d',i),3:12,'UniformOutput',false), {'SOH'}];
lbl_d = [arrayfun(@(i)sprintf('Seg%02d',i),1:11,'UniformOutput',false), {'SOH'}];

%% Colormap: blue(neg)→white→red(pos)
function cmap = rwb_cmap(n)
    n2 = floor(n/2);
    c1 = [linspace(0,1,n2)', linspace(0,1,n2)', ones(n2,1)];
    c2 = [ones(n-n2,1), linspace(1,0,n-n2)', linspace(1,0,n-n2)'];
    cmap = [c1; c2];
end

%% PCC matrix (NaN-aware)
function R = pcc_mat(X)
    n = size(X,2); R = eye(n);
    for i=1:n
        for j=i+1:n
            m = ~isnan(X(:,i)) & ~isnan(X(:,j));
            if sum(m)>3
                r = corr(X(m,i),X(m,j)); R(i,j)=r; R(j,i)=r;
            else
                R(i,j)=NaN; R(j,i)=NaN;
            end
        end
    end
end

%% Draw one correlation matrix panel
function draw_cmat(ax, R, lbls, cr_lbl, feat_type)
    n = size(R,1);
    imagesc(ax, R, [-1 1]);
    colormap(ax, rwb_cmap(256));
    for i=1:n
        for j=1:n
            if i==j, continue; end
            if ~isnan(R(i,j))
                clr='k'; if abs(R(i,j))>0.65, clr='w'; end
                text(ax,j,i,sprintf('%.2f',R(i,j)),'HorizontalAlignment','center','FontSize',7,'Color',clr);
            end
        end
    end
    hold(ax,'on');
    plot(ax,[0.5 n+0.5],[0.5 n+0.5],'k-','LineWidth',1.5);
    % SOH border
    rectangle(ax,'Position',[n-0.5,0.5,1,n],'EdgeColor',[0 0 0],'LineWidth',2);
    rectangle(ax,'Position',[0.5,n-0.5,n,1],'EdgeColor',[0 0 0],'LineWidth',2);
    set(ax,'XTick',1:n,'XTickLabel',lbls,'YTick',1:n,'YTickLabel',lbls,...
        'FontSize',8,'TickLength',[0 0]);
    xtickangle(ax,45);
    title(ax,sprintf('%s @ %s',feat_type,cr_lbl),'FontSize',10,'FontWeight','bold');
end

%% ===== CHARGE: one figure per C-rate =====
for ci = 1:5
    cond = uConds{ci};
    idx  = strcmp(FM.Condition, cond);
    dQr  = FM{idx,chg_cols} ./ sum(FM{idx,chg_cols},2,'omitnan');
    soh  = FM.Static_Capacity(idx) / Q_nom * 100;
    Xc   = [dQr, soh];
    Rc   = pcc_mat(Xc);

    fig = figure('Position',[50 50 680 620]);
    ax  = axes('Parent',fig);
    draw_cmat(ax, Rc, lbl_c, cr_labels{ci}, 'Charge dQ Ratio');
    cb  = colorbar(ax); cb.Label.String = 'Pearson Correlation Coefficient'; cb.FontSize=9;
    sgtitle(sprintf('Correlation Matrix: Charge dQ Ratio Features vs. SOH (%s)', cr_labels{ci}),...
        'FontSize',11,'FontWeight','bold');
    fname = fullfile(visDir, sprintf('CorrMatrix_Charge_%s.png', strrep(cr_labels{ci},'.','p')));
    saveas(fig, fname);
    close(fig);
    fprintf('Saved: %s\n', fname);
end

%% ===== DISCHARGE: one figure per C-rate =====
for ci = 1:5
    cond = uConds{ci};
    idx  = strcmp(FM.Condition, cond);
    dQr  = FM{idx,dch_cols} ./ sum(FM{idx,dch_cols},2,'omitnan');
    soh  = FM.Static_Capacity(idx) / Q_nom * 100;
    Xd   = [dQr, soh];
    Rd   = pcc_mat(Xd);

    fig = figure('Position',[50 50 720 660]);
    ax  = axes('Parent',fig);
    draw_cmat(ax, Rd, lbl_d, cr_labels{ci}, 'Discharge dQ Ratio');
    cb  = colorbar(ax); cb.Label.String = 'Pearson Correlation Coefficient'; cb.FontSize=9;
    sgtitle(sprintf('Correlation Matrix: Discharge dQ Ratio Features vs. SOH (%s)', cr_labels{ci}),...
        'FontSize',11,'FontWeight','bold');
    fname = fullfile(visDir, sprintf('CorrMatrix_Discharge_%s.png', strrep(cr_labels{ci},'.','p')));
    saveas(fig, fname);
    close(fig);
    fprintf('Saved: %s\n', fname);
end

fprintf('All done.\n');
