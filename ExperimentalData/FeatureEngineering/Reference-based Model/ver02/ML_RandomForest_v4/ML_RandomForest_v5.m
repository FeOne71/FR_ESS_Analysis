%% ML_RandomForest_v5.m
% =====================================================================
% Phase 2: 최적 세그먼트 + C_eff 기반 최종 ML 모델 훈련 및 저장
% - Track A (넓은 관측): Seg 10 + 9 + 5 + C_eff (충전), Seg 5 + 6 + C_eff (방전)
% - Track B (좁은 관측): Seg 9 + 8 + 7 + C_eff (충전), Seg 6 + 7 + 8 + C_eff (방전)
% - SOH, LLI, LAM 모델 각각 생성
% - RF, SVM, GPR 모델 생성 및 비교
% - `ML_Models_v5.mat`로 최종 저장 (필드 적용용)
% =====================================================================
clear; clc; close all;
rng(42); % 재현성

%% 1. 데이터 로드
d_pcc = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\ver02\ML_RandomForest_v4\PCC_v5b_data.mat');
d_cap = load('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Capacity_Trend_Figures\Capacity_Data_Static.mat');

% 피처
X_chg     = d_pcc.dQ_chg_all;  % [N x 11]
X_dch     = d_pcc.dQ_dch_all;  % [N x 11]
c_eff     = d_pcc.data_CrateNum; % [N x 1]
cell_ids  = d_pcc.data_CellID;   % [N x 1]
cyc_nums  = d_pcc.data_Cycle;    % [N x 1]

% 라벨 (SOH는 이미 있음, LLI/LAM 추가)
y_soh = d_pcc.data_SOH;
y_lli = nan(size(y_soh));
y_lam = nan(size(y_soh));

channels = unique(cell_ids);
for ci = 1:length(channels)
    ch = channels{ci};
    idx_ch = strcmp(cell_ids, ch);
    
    if ~isfield(d_cap.allChannelsCapacity, ch), continue; end
    cap_d = d_cap.allChannelsCapacity.(ch);
    
    cycs = cyc_nums(idx_ch);
    for i = 1:length(cycs)
        cn = cycs(i);
        id_c = find(cap_d.cycles == cn, 1);
        if ~isempty(id_c)
            % 원래 코드(v1~v4)의 LLI, LAM 할당 로직 참고
            if size(cap_d.Q, 1) >= 4 % LLI
                y_lli(idx_ch & cyc_nums==cn) = cap_d.Q{4, id_c}; 
            end
            if size(cap_d.Q, 1) >= 5 % LAM
                y_lam(idx_ch & cyc_nums==cn) = cap_d.Q{5, id_c};
            end
        end
    end
end

% 기본 유효성 검사 (SOH 기준)
valid_all = ~isnan(y_soh) & ~isnan(c_eff);
X_chg = X_chg(valid_all,:); X_dch = X_dch(valid_all,:);
y_soh = y_soh(valid_all); y_lli = y_lli(valid_all); y_lam = y_lam(valid_all);
c_eff = c_eff(valid_all); cell_ids = cell_ids(valid_all);

%% 2. 모델 구성 설정
% Track A (넓은 관측 연도)
conf.TrackA_CHG.segs = [10, 9, 5];
conf.TrackA_DCH.segs = [5, 6];
% Track B (좁은 관측 연도)
conf.TrackB_CHG.segs = [9, 8, 7];
conf.TrackB_DCH.segs = [6, 7, 8];

labels_data = {'SOH', y_soh; 'LLI', y_lli; 'LAM', y_lam};
model_types = {'RF', 'SVM', 'GPR'};
tracks = {'TrackA', 'TrackB'};
modes  = {'CHG', 'DCH'};

% 결과 저장을 위한 구조체
Final_Models = struct();
Results_CV   = struct();

% CV 그룹
u_cells = unique(cell_ids); K = 5;
cv_indices = zeros(size(y_soh)); groups = ceil(linspace(0.01, K, length(u_cells)));
for k=1:K, tc=u_cells(groups==k); for ti=1:length(tc), cv_indices(strcmp(cell_ids, tc{ti})) = k; end; end

%% 3. 최종 모델 훈련 및 CV 평가 루프
fprintf('=== Phase 2: Model Training (Final Models v5) ===\n');

for t_idx = 1:length(tracks)
    trk = tracks{t_idx};
    for m_idx = 1:length(modes)
        md = modes{m_idx};
        segs = conf.([trk '_' md]).segs;
        
        if strcmp(md,'CHG')
            X_base = X_chg(:, segs);
        else
            X_base = X_dch(:, segs);
        end
        X_comb = [X_base, c_eff]; % [세그먼트들, C_eff]
        
        valid_comb = ~any(isnan(X_comb), 2);
        Xc = X_comb(valid_comb,:);
        cvi = cv_indices(valid_comb);
        
        fprintf('\n>> [%s_%s] Features: Segs [%s] + C_eff (N=%d)\n', ...
            trk, md, num2str(segs), size(Xc,1));
        
        for l_idx = 1:size(labels_data,1)
            lbl_name = labels_data{l_idx,1};
            yc = labels_data{l_idx,2}(valid_comb);
            
            % NaN 레이블 제거 (LLI, LAM 등)
            valid_lbl = ~isnan(yc);
            Xc_l = Xc(valid_lbl,:); yc_l = yc(valid_lbl); cvi_l = cvi(valid_lbl);
            
            if isempty(Xc_l), continue; end
            
            for mt_idx = 1:length(model_types)
                mtype = model_types{mt_idx};
                y_pred = nan(size(yc_l));
                
                % CV 평가
                for k = 1:K
                    trn = (cvi_l~=k); tst = (cvi_l==k);
                    if sum(trn)==0 || sum(tst)==0, continue; end
                    
                    switch mtype
                        case 'RF'
                            model = TreeBagger(50, Xc_l(trn,:), yc_l(trn), 'Method','regression','MinLeafSize',5);
                            y_pred(tst) = predict(model, Xc_l(tst,:));
                        case 'SVM'
                            model = fitrsvm(Xc_l(trn,:), yc_l(trn), 'KernelFunction','gaussian','Standardize',true);
                            y_pred(tst) = predict(model, Xc_l(tst,:));
                        case 'GPR'
                            model = fitrgp(Xc_l(trn,:), yc_l(trn), 'KernelFunction','squaredexponential','Standardize',true);
                            y_pred(tst) = predict(model, Xc_l(tst,:));
                    end
                end
                
                % CV Metric
                rmse = sqrt(mean((yc_l - y_pred).^2, 'omitnan'));
                r2 = 1 - (sum((yc_l - y_pred).^2,'omitnan') / sum((yc_l - mean(yc_l,'omitnan')).^2,'omitnan'));
                Results_CV.(trk).(md).(lbl_name).(mtype).R2 = r2;
                Results_CV.(trk).(md).(lbl_name).(mtype).RMSE = rmse;
                Results_CV.(trk).(md).(lbl_name).(mtype).y_true = yc_l;
                Results_CV.(trk).(md).(lbl_name).(mtype).y_pred = y_pred;
                
                if strcmp(lbl_name,'SOH') % SOH만 출력
                    fprintf('  %3s | %3s | R2=%5.3f, RMSE=%5.3f\n', lbl_name, mtype, r2, rmse);
                end
                
                % === 최종 전체 데이터로 모델 재훈련 및 저장 ===
                switch mtype
                    case 'RF'
                        final_mdl = TreeBagger(50, Xc_l, yc_l, 'Method','regression');
                    case 'SVM'
                        final_mdl = fitrsvm(Xc_l, yc_l, 'KernelFunction','gaussian','Standardize',true);
                    case 'GPR'
                        final_mdl = fitrgp(Xc_l, yc_l, 'KernelFunction','squaredexponential','Standardize',true);
                end
                
                Final_Models.(trk).(md).(lbl_name).(mtype).model = final_mdl;
                Final_Models.(trk).(md).(lbl_name).(mtype).features = segs; % 어떤 세그먼트 썼는지 기록
                Final_Models.(trk).(md).(lbl_name).(mtype).feature_names = [arrayfun(@(i) sprintf('Seg%d',i), segs, 'UniformOutput',false), {'C_eff'}];
            end
        end
    end
end

%% 4. 결과 저장
saveDir = fileparts(mfilename('fullpath'));
save(fullfile(saveDir, 'ML_Models_v5.mat'), 'Final_Models', 'Results_CV');
fprintf('\n>> Saved: ML_Models_v5.mat (Contains Final Models for Field App)\n');

%% [선택] 최고 성능 모델 Parity Plot (Track A CHG SOH RF)
fig = figure('Position',[100,100,1200,400],'Name','Parity Plot: Track A vs Track B');
y_true_A = Results_CV.TrackA.CHG.SOH.RF.y_true;
y_pred_A = Results_CV.TrackA.CHG.SOH.RF.y_pred;
y_true_B = Results_CV.TrackB.CHG.SOH.RF.y_true;
y_pred_B = Results_CV.TrackB.CHG.SOH.RF.y_pred;

subplot(1,2,1); hold on; grid on;
scatter(y_true_A, y_pred_A, 50, c_eff(1:length(y_true_A)), 'filled', 'MarkerEdgeColor','k');
colormap(gca, turbo); cb=colorbar; ylabel(cb,'C-rate'); clim([0.1 3]);
plot([75 102], [75 102], 'r--', 'LineWidth',2);
xlabel('True SOH (%)'); ylabel('Predicted SOH (%)');
r2_a = Results_CV.TrackA.CHG.SOH.RF.R2;
title(sprintf('Track A (Seg 10+9+5+C_{eff})\nRF SOH R^2 = %.3f', r2_a));
axis square; xlim([75 102]); ylim([75 102]);

subplot(1,2,2); hold on; grid on;
scatter(y_true_B, y_pred_B, 50, c_eff(1:length(y_true_B)), 'filled', 'MarkerEdgeColor','k');
colormap(gca, turbo); cb=colorbar; ylabel(cb,'C-rate'); clim([0.1 3]);
plot([75 102], [75 102], 'r--', 'LineWidth',2);
xlabel('True SOH (%)'); ylabel('Predicted SOH (%)');
r2_b = Results_CV.TrackB.CHG.SOH.RF.R2;
title(sprintf('Track B (Seg 9+8+7+C_{eff})\nRF SOH R^2 = %.3f', r2_b));
axis square; xlim([75 102]); ylim([75 102]);

saveas(fig, fullfile(saveDir, 'Parity_Plot_v5.fig'));
fprintf('>> Saved: Parity_Plot_v5.fig\n');
