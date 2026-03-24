% Vis_RF_v2.m  — SOH Parity Plot + Feature Importance
% (Train_RF_v2.m 에서 호출됨)

fprintf('\n=== Vis_RF_v2: SOH Parity Plot ===\n');

save_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', ...
    'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02', 'ML_RandomForest_v2');

if ~exist('Results_RF', 'var')
    load(fullfile(save_dir, 'Result_RF_v2.mat'), 'Results_RF');
end

%% 1. Parity Plot
yt = Results_RF.SOH.Y_true;
yp = Results_RF.SOH.Y_pred;
r2  = Results_RF.SOH.R2;
rmse = Results_RF.SOH.RMSE;
mae = Results_RF.SOH.MAE;

fig1 = figure('Name','RF_v2 Parity','Position',[100 100 550 500],'Visible','off');
hold on; grid on; box on;
scatter(yt, yp, 40, [0.2 0.4 0.85], 'filled', 'MarkerFaceAlpha', 0.6);
plot([min(yt)-2, max(yt)+2], [min(yt)-2, max(yt)+2], 'k--', 'LineWidth', 1.5);
xlabel('True SOH (%)','FontSize',12); ylabel('Predicted SOH (%)','FontSize',12);
title(sprintf('RF v2 (dQ only): SOH Parity Plot\nR²=%.4f | RMSE=%.3f | MAE=%.3f', r2, rmse, mae), ...
    'FontSize',12,'FontWeight','bold');
axis equal; xlim([min(yt)-2 max(yt)+2]); ylim([min(yt)-2 max(yt)+2]);
saveas(fig1, fullfile(save_dir, 'RF_v2_ParityPlot_SOH.fig'));
close(fig1);

%% 2. Feature Importance
imp = Results_RF.SOH.Importance;
fn  = Results_RF.feature_names;

fig2 = figure('Name','RF_v2 Importance','Position',[100 100 700 400],'Visible','off');
bar(imp, 'FaceColor', [0.3 0.6 0.9]); grid on; box on;
set(gca, 'XTick', 1:length(fn), 'XTickLabel', fn, 'FontSize', 9);
xtickangle(45);
ylabel('Permutation Importance','FontSize',11);
title('RF v2 - Feature Importance (dQ only, SOH)', 'FontSize',12,'FontWeight','bold');
saveas(fig2, fullfile(save_dir, 'RF_v2_FeatureImportance.fig'));
close(fig2);

fprintf('Figures saved to: %s\n', save_dir);
