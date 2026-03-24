% Vis_RF_v4.m — 충전/방전 분리 모델 시각화
% =====================================================================
% 1. Parity Plot (충전 vs 방전 모델 비교, 각 라벨)
% 2. Feature Importance (충전 / 방전 각각)
% 3. Error Distribution (충전 vs 방전 비교)
% =====================================================================

if ~exist('Results_v4', 'var')
    load(fullfile(save_dir, 'Result_RF_v4.mat'), 'Results_v4');
end

label_names = Results_v4.label_names;
model_types = {'chg', 'dch'};
model_labels_kor = {'충전 모델', '방전 모델'};
colors = {[0.2, 0.5, 0.85], [0.85, 0.35, 0.25]}; % blue, red

%% ============================================================
% 1. Parity Plot (1×2 per label → 3 figures: SOH, LLI, LAM)
%% ============================================================
for lbl = 1:length(label_names)
    lname = label_names{lbl};
    fig = figure('Name', sprintf('RF_v4 Parity - %s', lname), ...
                 'Position', [100, 200, 900, 400], 'Visible', 'off');
    
    for m = 1:2
        mtype = model_types{m};
        R = Results_v4.(mtype).(lname);
        
        subplot(1, 2, m);
        scatter(R.Y_true, R.Y_pred, 20, colors{m}, 'filled', 'MarkerFaceAlpha', 0.5);
        hold on;
        lims = [min([R.Y_true; R.Y_pred])*0.98, max([R.Y_true; R.Y_pred])*1.02];
        plot(lims, lims, 'k--', 'LineWidth', 1.2);
        xlabel(sprintf('True %s', lname)); ylabel(sprintf('Predicted %s', lname));
        title(sprintf('%s (R²=%.3f, RMSE=%.4f)', model_labels_kor{m}, R.R2, R.RMSE), ...
              'FontSize', 11, 'FontWeight', 'bold');
        grid on; axis equal; xlim(lims); ylim(lims);
    end
    sgtitle(sprintf('[v4] %s — Charge vs Discharge Model', lname), 'FontSize', 13, 'FontWeight', 'bold');
    saveas(fig, fullfile(save_dir, sprintf('RF_v4_Parity_%s.fig', lname)));
    close(fig);
end

%% ============================================================
% 2. Feature Importance (1×2 subplot — chg / dch)
%% ============================================================
for lbl = 1:length(label_names)
    lname = label_names{lbl};
    fig = figure('Name', sprintf('RF_v4 Importance - %s', lname), ...
                 'Position', [100, 200, 1100, 400], 'Visible', 'off');
    
    for m = 1:2
        mtype = model_types{m};
        R = Results_v4.(mtype).(lname);
        fnames = Results_v4.(mtype).feature_names;
        imp = R.Importance;
        
        subplot(1, 2, m);
        barh(imp, 'FaceColor', colors{m});
        yticks(1:length(fnames)); yticklabels(fnames);
        xlabel('Permutation Importance');
        title(sprintf('%s — %s', model_labels_kor{m}, lname), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
    end
    sgtitle(sprintf('[v4] Feature Importance — %s', lname), 'FontSize', 13, 'FontWeight', 'bold');
    saveas(fig, fullfile(save_dir, sprintf('RF_v4_Importance_%s.fig', lname)));
    close(fig);
end

%% ============================================================
% 3. Error Distribution (Charge vs Discharge overlay)
%% ============================================================
fig = figure('Name', 'RF_v4 Error Distribution', ...
             'Position', [100, 200, 1200, 350], 'Visible', 'off');

for lbl = 1:length(label_names)
    lname = label_names{lbl};
    subplot(1, 3, lbl);
    for m = 1:2
        mtype = model_types{m};
        R = Results_v4.(mtype).(lname);
        err = R.Y_pred - R.Y_true;
        histogram(err, 20, 'FaceColor', colors{m}, 'FaceAlpha', 0.5, ...
                  'EdgeColor', 'none', 'DisplayName', model_labels_kor{m});
        hold on;
    end
    xlabel(sprintf('%s Error', lname)); ylabel('Count');
    title(lname, 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
end
sgtitle('[v4] CV Error Distribution — Charge vs Discharge', 'FontSize', 13, 'FontWeight', 'bold');
saveas(fig, fullfile(save_dir, 'RF_v4_ErrorDist.fig'));
close(fig);

fprintf(' >> Vis_RF_v4: All figures saved to %s\n', save_dir);
