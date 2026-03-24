% RM_02_Visualization.m
% LOCO-CV 모델 비교 Bar chart + HP 설정 테이블 시각화
% RM_01 결과를 바탕으로 PPT용 시각화 생성

clear; clc; close all;

%% Paths
verDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver0317';
rmDir  = fullfile(verDir, 'RatioModel');
visDir = fullfile(rmDir, 'Visualization');
Q_nom  = 64;

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% ===== Part 1: HP Settings Table =====
% Hyperparameter 설정을 시각적 테이블로 표시
fig1 = figure('Position',[50 50 1000 360]);
ax1 = axes('Parent',fig1);
axis(ax1,'off');
title(ax1,'Hyperparameter Settings (Initial Default Values)','FontSize',13,'FontWeight','bold');

col_names = {'Model','Key Hyperparameters','NaN Handling','Notes'};
row_data = {
    'LASSO',    'α: 4-fold CV (min MSE)',                  'Impute 0 (z-space)',     'Regularization auto-selected';
    'RF',       'Trees=100, MaxSplits=20, MinLeaf=5',       'Impute 0 (z-space)',     'Bagging ensemble';
    'GBM',      'Trees=100, LR=0.1, MaxSplits=10, MinLeaf=5','Impute 0 (z-space)',   'Gradient boosting';
    'XGBoost',  'n=100, depth=5, LR=0.1, subsample=0.8',   'Default routing†',       'Python (xgboost)';
    'LightGBM', 'n=100, depth=5, LR=0.1, min_samples=5',   'Default routing†',       'Python (lightgbm)';
};

% Create table using text
nRows = 5; nCols = 4;
colW = [0.12 0.38 0.25 0.20];
rowH = 0.14;
colX = [0.02 cumsum(colW(1:end-1))+0.02];
colors = {[0.95 0.95 0.95],[1 1 1]};
header_color = [0.18 0.31 0.49];

% Header
for c = 1:nCols
    x0 = colX(c); w = colW(c);
    rectangle(ax1,'Position',[x0,0.82,w-0.01,rowH],'FaceColor',header_color,'EdgeColor','w','LineWidth',0.5);
    text(ax1, x0+w/2-0.005, 0.82+rowH/2, col_names{c},'HorizontalAlignment','center',...
        'FontSize',10,'FontWeight','bold','Color','w','Interpreter','none');
end

% Rows
for r = 1:nRows
    y0 = 0.82 - r*rowH;
    for c = 1:nCols
        x0 = colX(c); w = colW(c);
        fc = colors{mod(r,2)+1};
        if c==1, fc=[0.86 0.91 0.96]; end  % model name column highlight
        rectangle(ax1,'Position',[x0,y0,w-0.01,rowH-0.01],'FaceColor',fc,'EdgeColor',[0.8 0.8 0.8],'LineWidth',0.5);
        text(ax1, x0+w/2-0.005, y0+rowH/2, row_data{r,c},'HorizontalAlignment','center',...
            'FontSize',8.5,'Interpreter','none','Color',[0.1 0.1 0.1]);
    end
end

% Footer note
text(ax1,0.02,0.03,'* Initial default HP. Hyperparameter optimization (e.g., Bayesian optimization) is planned as future work.',...
    'FontSize',8.5,'Color',[0.5 0.5 0.5],'Interpreter','none');
text(ax1,0.02,0.01,'† XGBoost/LightGBM native NaN routing was not trained (no NaN in Lab data); default direction applied at inference.',...
    'FontSize',8.5,'Color',[0.5 0.5 0.5],'Interpreter','none');

set(ax1,'XLim',[0 1],'YLim',[0 1]);
saveas(fig1, fullfile(visDir,'RM_HP_Settings_Table.png'));
fprintf('Saved: RM_HP_Settings_Table.png\n');

%% ===== Part 2: LOCO-CV Model Comparison =====
% Run LOCO-CV with full feature set (all valid segments + C_eff)
d = load(fullfile(verDir,'FeatureMatrix_ver0317.mat')); FM = d.FM;

chg_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i),3:12,'UniformOutput',false);
dch_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i),1:11,'UniformOutput',false);
dQ_c_ratio = FM{:,chg_cols} ./ sum(FM{:,chg_cols},2,'omitnan');
dQ_d_ratio = FM{:,dch_cols} ./ sum(FM{:,dch_cols},2,'omitnan');
X_lab = [dQ_c_ratio, dQ_d_ratio, FM.C_eff_chg, FM.C_eff_dch];
y_lab = FM.Static_Capacity / Q_nom * 100;

mu_lab = mean(X_lab,1,'omitnan'); sigma_lab = std(X_lab,0,1,'omitnan');
sigma_lab(sigma_lab==0) = 1;
X_lab_s = (X_lab - mu_lab) ./ sigma_lab;
X_lab_s(isnan(X_lab_s)) = 0;  % NaN → 0 for LASSO

cellIDs = unique(FM.CellID); nCells = length(cellIDs);
modelNames = {'LASSO','RF','GBM','XGBoost','LightGBM'};
nModels = 5;
res_rmse = zeros(1,nModels); res_mape = zeros(1,nModels); res_r2 = zeros(1,nModels);
all_preds = cell(1,nModels); all_trues = cell(1,nModels);

for m = 1:nModels
    mName = modelNames{m}; fprintf('  %s ... ', mName);
    y_true_all = []; y_pred_all = [];
    for fold = 1:nCells
        testIdx = FM.CellID == cellIDs(fold);
        trainIdx = ~testIdx;
        X_tr = X_lab_s(trainIdx,:); y_tr = y_lab(trainIdx);
        X_te = X_lab_s(testIdx,:);  y_te = y_lab(testIdx);
        switch mName
            case 'LASSO'
                [B,FI] = lasso(X_tr, y_tr, 'CV', 4);
                y_p = X_te * B(:,FI.IndexMinMSE) + FI.Intercept(FI.IndexMinMSE);
            case 'RF'
                mdl = fitrensemble(X_tr,y_tr,'Method','Bag','NumLearningCycles',100,...
                    'Learners',templateTree('MaxNumSplits',20,'MinLeafSize',5,'Surrogate','on'));
                y_p = predict(mdl,X_te);
            case 'GBM'
                mdl = fitrensemble(X_tr,y_tr,'Method','LSBoost','NumLearningCycles',100,...
                    'LearnRate',0.1,'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',5,'Surrogate','on'));
                y_p = predict(mdl,X_te);
            case 'XGBoost'
                xgb = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
                    'learning_rate',0.1,'min_child_weight',int32(3),'subsample',0.8,'random_state',int32(42)));
                xgb.fit(py.numpy.array(X_tr),py.numpy.array(y_tr));
                y_p = double(xgb.predict(py.numpy.array(X_te)));
            case 'LightGBM'
                lgb = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(100),'max_depth',int32(5),...
                    'learning_rate',0.1,'min_child_samples',int32(5),'subsample',0.8,...
                    'random_state',int32(42),'verbose',int32(-1)));
                lgb.fit(py.numpy.array(X_tr),py.numpy.array(y_tr));
                y_p = double(lgb.predict(py.numpy.array(X_te)));
        end
        y_true_all = [y_true_all; y_te]; y_pred_all = [y_pred_all; y_p(:)];
    end
    err = y_true_all - y_pred_all;
    res_rmse(m) = sqrt(mean(err.^2));
    res_mape(m) = mean(abs(err./y_true_all))*100;
    res_r2(m) = 1 - sum(err.^2)/sum((y_true_all-mean(y_true_all)).^2);
    all_trues{m} = y_true_all; all_preds{m} = y_pred_all;
    fprintf('RMSE=%.3f%%, R2=%.4f\n', res_rmse(m), res_r2(m));
end

%% Bar chart: RMSE + R²
colors_bar = [0.31 0.51 0.78; 0.47 0.72 0.44; 0.91 0.49 0.31; 0.69 0.53 0.79; 0.87 0.76 0.30];
fig2 = figure('Position',[50 50 950 380]);
sgtitle('LOCO-CV Model Comparison (Full dQ Ratio Feature Set, 23 features)','FontSize',12,'FontWeight','bold');

ax2a = subplot(1,2,1);
b1 = bar(ax2a, res_rmse, 'FaceColor','flat');
b1.CData = colors_bar;
set(ax2a,'XTickLabel',modelNames,'FontSize',10,'TickLength',[0 0]);
ylabel(ax2a,'RMSE (%)','FontSize',11);
title(ax2a,'LOCO-CV RMSE','FontSize',11,'FontWeight','bold');
for i=1:nModels
    text(ax2a,i,res_rmse(i)+0.01,sprintf('%.2f%%',res_rmse(i)),'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end
ylim(ax2a,[0 max(res_rmse)*1.2]);
grid(ax2a,'on'); grid(ax2a,'minor');

ax2b = subplot(1,2,2);
b2 = bar(ax2b, res_r2, 'FaceColor','flat');
b2.CData = colors_bar;
set(ax2b,'XTickLabel',modelNames,'FontSize',10,'TickLength',[0 0]);
ylabel(ax2b,'R²','FontSize',11);
title(ax2b,'LOCO-CV R²','FontSize',11,'FontWeight','bold');
for i=1:nModels
    text(ax2b,i,res_r2(i)+0.005,sprintf('%.3f',res_r2(i)),'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end
ylim(ax2b,[0.6 max(res_r2)*1.05]);
grid(ax2b,'on'); grid(ax2b,'minor');

saveas(fig2, fullfile(visDir,'RM_LOCO_CV_Full_Comparison.png'));
fprintf('Saved: RM_LOCO_CV_Full_Comparison.png\n');

%% Scatter: Predicted vs Actual
fig3 = figure('Position',[50 50 1200 250]);
sgtitle('LOCO-CV: Predicted vs. Actual SOH (Full Feature Set)','FontSize',12,'FontWeight','bold');
for m = 1:nModels
    ax = subplot(1,5,m);
    scatter(ax, all_trues{m}, all_preds{m}, 20, colors_bar(m,:), 'filled','MarkerFaceAlpha',0.7);
    hold(ax,'on');
    ref = [min(all_trues{m}) max(all_trues{m})];
    plot(ax,ref,ref,'k--','LineWidth',1);
    title(ax,modelNames{m},'FontSize',10,'FontWeight','bold');
    xlabel(ax,'Actual SOH (%)','FontSize',9); ylabel(ax,'Predicted SOH (%)','FontSize',9);
    text(ax,0.05,0.92,sprintf('RMSE=%.2f%%\nR²=%.3f',res_rmse(m),res_r2(m)),'Units','normalized','FontSize',8,'Color','k');
    axis(ax,'equal'); grid(ax,'on');
end
saveas(fig3, fullfile(visDir,'RM_LOCO_CV_Scatter.png'));
fprintf('Saved: RM_LOCO_CV_Scatter.png\n');

fprintf('\nAll visualizations saved to: %s\n', visDir);
