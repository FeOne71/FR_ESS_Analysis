% step1_LOCO_CV.m
% RelativeDQ: Lab 피처 생성 + C-rate별 LOCO-CV + 최종 GBM 학습 + 결과 저장
% 분모 = Seg08 + Seg09 (4개년 전부 커버, 스케일 고정)

clear; clc; close all;

%% Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
rdqDir  = fullfile(verDir, 'RatioModel', 'RelativeDQ');
visDir  = fullfile(rdqDir, 'Visualization');
if ~exist(visDir,'dir'), mkdir(visDir); end

%% Load Lab Data
d  = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat')); FM = d.FM;
Q_nom = 64;

dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);
dQ_c_raw = FM{:, dQ_c_cols};  % 215×10: col6=Seg08, col7=Seg09
dQ_d_raw = FM{:, dQ_d_cols};  % 215×11: col8=Seg08, col9=Seg09

C_eff_c = FM.C_eff_chg; C_eff_d = FM.C_eff_dch;
y_lab   = FM.Static_Capacity / Q_nom * 100;
cellIDs = FM.CellID;
conditions = FM.Condition;

%% Compute Relative dQ (fixed denominator = Seg08 + Seg09)
denom_c = dQ_c_raw(:,6) + dQ_c_raw(:,7);  % Seg08 + Seg09 (charge)
denom_d = dQ_d_raw(:,8) + dQ_d_raw(:,9);  % Seg08 + Seg09 (discharge)
rdQ_c = dQ_c_raw ./ denom_c;  % 215×10
rdQ_d = dQ_d_raw ./ denom_d;  % 215×11

X_lab = [rdQ_c, rdQ_d, C_eff_c, C_eff_d];  % 215×23
feat_names = [arrayfun(@(i) sprintf('rdQ_c%02d',i),3:12,'UniformOutput',false), ...
              arrayfun(@(i) sprintf('rdQ_d%02d',i),1:11,'UniformOutput',false), ...
              {'Ceff_c','Ceff_d'}];

fprintf('Relative dQ: %d rows × %d features, NaN=%d\n', ...
    size(X_lab,1), size(X_lab,2), sum(isnan(X_lab(:))));

% Standardize (전체 기준, NaN 제외)
mu_lab  = mean(X_lab,1,'omitnan');
sig_lab = std(X_lab,0,1,'omitnan'); sig_lab(sig_lab==0)=1;
X_lab_s = (X_lab - mu_lab) ./ sig_lab;

%% C-rate Stratified LOCO-CV
uConds = unique(conditions);
fprintf('\n=== C-rate Stratified LOCO-CV (GBM Surrogate Split) ===\n');
fprintf('%-8s %-5s %-7s %-7s\n','Cond','N','RMSE%','R2');

CV_pred = nan(size(y_lab));
cond_results = table();

for ci = 1:length(uConds)
    cond = uConds(ci);
    cM = strcmp(conditions, cond);
    X_c=X_lab_s(cM,:); y_c=y_lab(cM); ids_c=cellIDs(cM);
    uCells=unique(ids_c); y_cv=nan(size(y_c));

    for fold=1:length(uCells)
        tc=uCells(fold);
        trM=~strcmp(ids_c,tc); teM=strcmp(ids_c,tc);
        mdl=fitrensemble(X_c(trM,:),y_c(trM),'Method','LSBoost',...
            'NumLearningCycles',200,'LearnRate',0.05,...
            'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on'));
        y_cv(teM)=predict(mdl,X_c(teM,:));
    end
    CV_pred(cM)=y_cv;

    err=y_c-y_cv;
    rmse=sqrt(mean(err.^2)); r2=1-sum(err.^2)/sum((y_c-mean(y_c)).^2);
    fprintf('%-8s %-5d %-7.3f %-7.4f\n',char(cond),sum(cM),rmse,r2);
    cond_results=[cond_results; table(cond,sum(cM),rmse,r2,'VariableNames',{'Condition','N','RMSE','R2'})];
end

overall_err=y_lab-CV_pred;
ov_rmse=sqrt(mean(overall_err.^2));
ov_r2=1-sum(overall_err.^2)/sum((y_lab-mean(y_lab)).^2);
fprintf('\nOverall LOCO-CV: RMSE=%.3f%%, R2=%.4f\n',ov_rmse,ov_r2);

%% Visualization: Parity Plot + CV per C-rate
fig=figure('Position',[50 50 1100 420]);

subplot(1,3,1);
scatter(y_lab,CV_pred,30,C_eff_c,'filled','MarkerEdgeColor','none'); colorbar;
hold on; lims=[85 100]; plot(lims,lims,'k--','LineWidth',1.5);
xlabel('True SOH (%)','FontSize',11); ylabel('Predicted SOH (%)','FontSize',11);
title(sprintf('LOCO-CV Parity (RMSE=%.2f%%, R²=%.3f)',ov_rmse,ov_r2),'FontSize',11);
grid on; xlim(lims); ylim(lims); axis square;
colorTitle=colorbar; colorTitle.Label.String='C_{eff,chg}';

subplot(1,3,2);
bar(cond_results.RMSE,'FaceColor',[0.3 0.6 0.9]);
set(gca,'XTick',1:height(cond_results),'XTickLabel',cond_results.Condition,'FontSize',11,'XTickLabelRotation',30);
ylabel('RMSE (%)'); title('RMSE by C-rate'); grid on;
for i=1:height(cond_results), text(i,cond_results.RMSE(i)+0.01,sprintf('%.3f',cond_results.RMSE(i)),'HorizontalAlignment','center','FontSize',9); end

subplot(1,3,3);
bar(cond_results.R2,'FaceColor',[0.8 0.4 0.2]);
set(gca,'XTick',1:height(cond_results),'XTickLabel',cond_results.Condition,'FontSize',11,'XTickLabelRotation',30);
ylabel('R²'); title('R² by C-rate'); grid on; ylim([0 1.05]);
for i=1:height(cond_results), text(i,cond_results.R2(i)+0.01,sprintf('%.3f',cond_results.R2(i)),'HorizontalAlignment','center','FontSize',9); end

sgtitle('Relative dQ (Fixed Denom: Seg8+Seg9): C-rate Stratified LOCO-CV','FontSize',13,'FontWeight','bold');
saveas(fig,fullfile(visDir,'rdQ_LOCO_CV.png'));
fprintf('Saved: rdQ_LOCO_CV.png\n');

%% Train Final GBM on ALL lab data
fprintf('\n=== Training final GBM ===\n');
mdl_GBM_final = fitrensemble(X_lab_s, y_lab, 'Method','LSBoost',...
    'NumLearningCycles',200,'LearnRate',0.05,...
    'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on'));
fprintf('  GBM trained.\n');

%% Save
save(fullfile(rdqDir,'rdQ_model.mat'),...
    'mdl_GBM_final','mu_lab','sig_lab','y_lab','CV_pred','cond_results','feat_names');
fprintf('Saved: rdQ_model.mat\n');
fprintf('=== step1 complete! ===\n');
