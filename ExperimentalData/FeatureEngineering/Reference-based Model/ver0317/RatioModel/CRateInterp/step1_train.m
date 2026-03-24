% step1_train.m  (CRateInterp)
% C-rate별 다중 모델 학습 + LOCO-CV → 최적 모델 선택 + 저장
% 모델: LASSO, RF, GBM, XGBoost, LightGBM (5종)
% 피처: dQ_c_03~12 (10) + dQ_d_01~11 (11) + C_eff_c + C_eff_d = 23개

clear; clc; close all;

%% Paths
projDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local';
verDir  = fullfile(projDir, 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
ciDir   = fullfile(verDir, 'RatioModel', 'CRateInterp');
visDir  = fullfile(ciDir, 'Visualization');
if ~exist(visDir,'dir'), mkdir(visDir); end

pyenv('Version', 'C:\Users\Chulwon Jung\AppData\Local\Programs\Python\Python311\python.exe');

%% Load Lab Data
d  = load(fullfile(verDir, 'FeatureMatrix_ver0317.mat')); FM = d.FM;
Q_nom = 64;
dQ_c_cols = arrayfun(@(i) sprintf('dQ_c_%02d',i), 3:12, 'UniformOutput', false);
dQ_d_cols = arrayfun(@(i) sprintf('dQ_d_%02d',i), 1:11, 'UniformOutput', false);
dQ_c_raw = FM{:, dQ_c_cols};
dQ_d_raw = FM{:, dQ_d_cols};

% dQ Ratio: each segment / sum(valid segments) — reduces C-rate scale dependency
rdQ_c = dQ_c_raw ./ sum(dQ_c_raw, 2, 'omitnan');  % 215×10
rdQ_d = dQ_d_raw ./ sum(dQ_d_raw, 2, 'omitnan');  % 215×11
X_lab = [rdQ_c, rdQ_d, FM.C_eff_chg, FM.C_eff_dch];  % 215×23, NaN=0
y_lab = FM.Static_Capacity / Q_nom * 100;
cellIDs = FM.CellID;
conditions = FM.Condition;

fprintf('Lab data (ratio): %d rows × %d features, NaN=%d\n', size(X_lab,1), size(X_lab,2), sum(isnan(X_lab(:))));

%% C-rate별 모델 학습 + LOCO-CV
uConds = unique(conditions);
crate_vals = [0.1, 0.5, 1.0, 2.0, 3.0];  % c01 c05 c1 c2 c3
modelNames = {'LASSO', 'RF', 'GBM', 'XGBoost', 'LightGBM'};

CRateModels = struct();  % 저장용

fprintf('\n=== C-rate별 LOCO-CV (5 models) ===\n');
fprintf('%-5s %-12s %-7s %-7s %-7s\n', 'Crate', 'BestModel', 'RMSE%', 'MAPE%', 'R2');

for ci = 1:length(uConds)
    cond = uConds(ci); cr = crate_vals(ci);
    cM   = strcmp(conditions, cond);
    X_c  = X_lab(cM,:); y_c = y_lab(cM); ids_c = cellIDs(cM);
    uCells = unique(ids_c);

    % Standardize (Train data statistics)
    mu_c  = mean(X_c,1); sig_c = std(X_c,0,1); sig_c(sig_c==0)=1;
    X_c_s = (X_c - mu_c) ./ sig_c;

    % Preallocate CV predictions
    y_cv = zeros(size(y_c,1), length(modelNames));

    for fold = 1:length(uCells)
        tc   = uCells(fold);
        trM  = ~strcmp(ids_c,tc); teM = strcmp(ids_c,tc);
        X_tr = X_c_s(trM,:); y_tr = y_c(trM);
        X_te = X_c_s(teM,:);

        % 1. LASSO
        [B,FI] = lasso(X_tr, y_tr, 'Alpha',1,'NumLambda',50,'CV',5);
        b0 = FI.Intercept(FI.IndexMinMSE);
        y_cv(teM,1) = X_te * B(:,FI.IndexMinMSE) + b0;

        % 2. RF
        mdl_rf = fitrensemble(X_tr,y_tr,'Method','Bag','NumLearningCycles',200,...
            'Learners',templateTree('MaxNumSplits',20,'MinLeafSize',3,'Surrogate','on'));
        y_cv(teM,2) = predict(mdl_rf, X_te);

        % 3. GBM
        mdl_gbm = fitrensemble(X_tr,y_tr,'Method','LSBoost','NumLearningCycles',200,...
            'LearnRate',0.05,'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on'));
        y_cv(teM,3) = predict(mdl_gbm, X_te);

        % 4. XGBoost
        mdl_xgb = py.xgboost.XGBRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(4),...
            'learning_rate',0.05,'min_child_weight',int32(3),'subsample',0.8,...
            'colsample_bytree',0.8,'random_state',int32(42)));
        mdl_xgb.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
        y_cv(teM,4) = double(mdl_xgb.predict(py.numpy.array(X_te)));

        % 5. LightGBM
        mdl_lgb = py.lightgbm.LGBMRegressor(pyargs('n_estimators',int32(200),'max_depth',int32(4),...
            'learning_rate',0.05,'min_child_samples',int32(3),'subsample',0.8,...
            'colsample_bytree',0.8,'random_state',int32(42),'verbose',int32(-1)));
        mdl_lgb.fit(py.numpy.array(X_tr), py.numpy.array(y_tr));
        y_cv(teM,5) = double(mdl_lgb.predict(py.numpy.array(X_te)));
    end

    % Compute metrics per model
    rmse_all = zeros(1,length(modelNames));
    mape_all = zeros(1,length(modelNames));
    r2_all   = zeros(1,length(modelNames));
    for m = 1:length(modelNames)
        e = y_c - y_cv(:,m);
        rmse_all(m) = sqrt(mean(e.^2));
        mape_all(m) = mean(abs(e./y_c))*100;
        r2_all(m)   = 1 - sum(e.^2)/sum((y_c-mean(y_c)).^2);
    end

    [~, best_m] = min(rmse_all);
    fprintf('%-5.1fC %-12s %-7.3f %-7.3f %-7.4f\n', ...
        cr, modelNames{best_m}, rmse_all(best_m), mape_all(best_m), r2_all(best_m));

    % Train final model (best) on ALL C-rate data
    switch best_m
        case 1  % LASSO
            [B_f,FI_f] = lasso(X_c_s, y_c,'Alpha',1,'NumLambda',50,'CV',5);
            final_mdl = struct('type','LASSO','B',B_f(:,FI_f.IndexMinMSE),'b0',FI_f.Intercept(FI_f.IndexMinMSE));
        case 2  % RF
            final_mdl = struct('type','RF','mdl',...
                fitrensemble(X_c_s,y_c,'Method','Bag','NumLearningCycles',200,...
                'Learners',templateTree('MaxNumSplits',20,'MinLeafSize',3,'Surrogate','on')));
        case 3  % GBM
            final_mdl = struct('type','GBM','mdl',...
                fitrensemble(X_c_s,y_c,'Method','LSBoost','NumLearningCycles',200,...
                'LearnRate',0.05,'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on')));
        case {4,5}  % XGB or LGB
            % Python model trained — retrain in step2 (can't save to .mat)
            % Save as GBM fallback, mark type for retraining
            final_mdl = struct('type', modelNames{best_m}, 'mdl', ...
                fitrensemble(X_c_s,y_c,'Method','LSBoost','NumLearningCycles',200,...
                'LearnRate',0.05,'Learners',templateTree('MaxNumSplits',10,'MinLeafSize',3,'Surrogate','on')));
    end

    cname = matlab.lang.makeValidName(char(cond));
    CRateModels.(cname) = struct('cond',cond,'crate',cr,...
        'best_model',modelNames{best_m},'best_idx',best_m,...
        'mu',mu_c,'sig',sig_c,'final',final_mdl,...
        'rmse_all',rmse_all,'r2_all',r2_all,'y_c',y_c,'y_cv',y_cv);
end

%% Save
save(fullfile(ciDir,'crate_models.mat'), 'CRateModels', 'crate_vals', 'uConds', 'modelNames');
fprintf('\nSaved: crate_models.mat\n=== step1 complete! ===\n');
