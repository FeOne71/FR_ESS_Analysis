%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 4: Correlation and Final Report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

projectDir = pwd;
reportDir = fullfile(projectDir, 'Final_Analysis_Report');
if ~exist(reportDir, 'dir')
    mkdir(reportDir);
end

rptFile = fullfile(projectDir, 'Results', 'Phase_1', 'RPT_Feature_Matrix.mat');
driveFile = fullfile(projectDir, 'Results', 'Phase_2', 'Driving_Feature_Matrix.mat');

if ~exist(rptFile, 'file')
    error('Missing %s', rptFile);
end
if ~exist(driveFile, 'file')
    error('Missing %s', driveFile);
end

load(rptFile, 'feature_matrix', 'feat_names', 'channels_all', 'rpt_points', 'discharge_ah', 'soh_percent');
load(driveFile, 'feature_matrix', 'feat_names', 'profiles', 'socs', 'cycles');

%% RPT feature vs capacity correlation
num_feats = numel(feat_names);
num_ch = numel(channels_all);
num_rpt = numel(rpt_points);

cap_vec = discharge_ah(:);
feat_mat = reshape(feature_matrix, num_ch * num_rpt, num_feats);

rho = NaN(1, num_feats);
for f = 1:num_feats
    rho(f) = corr(feat_mat(:, f), cap_vec, 'rows', 'complete');
end

% Heatmap
fig = figure('Visible', 'off', 'Position', [100 100 1200 400]);
imagesc(rho);
colormap('jet');
colorbar;
set(gca, 'XTick', 1:num_feats, 'XTickLabel', feat_names, 'XTickLabelRotation', 90);
set(gca, 'YTick', 1, 'YTickLabel', {'RPT'});
title('RPT Feature vs Discharge Capacity Correlation');
saveas(fig, fullfile(reportDir, 'RPT_Correlation_Heatmap.png'));
close(fig);

%% RPT vs Driving feature matching at RPT cycles
rpt_cycles = zeros(1, numel(rpt_points));
for i = 1:numel(rpt_points)
    token = regexp(rpt_points{i}, '(\d+)cyc', 'tokens', 'once');
    if ~isempty(token)
        rpt_cycles(i) = str2double(token{1});
    else
        rpt_cycles(i) = i;
    end
end

match_idx = arrayfun(@(x) find(cycles == x, 1), rpt_cycles, 'UniformOutput', false);
match_idx = cellfun(@(x) ifelse(isempty(x), NaN, x), match_idx);

drive_corr = NaN(num_feats, numel(profiles), numel(socs));
for p = 1:numel(profiles)
    for s = 1:numel(socs)
        for f = 1:num_feats
            y = squeeze(feature_matrix(p, s, :, f));
            y_match = y(match_idx);
            x = discharge_ah(1, :); % use channel 1 as reference
            if numel(y_match) == numel(x)
                drive_corr(f, p, s) = corr(y_match(:), x(:), 'rows', 'complete');
            end
        end
    end
end

% Save correlation matrices
save(fullfile(reportDir, 'Correlation_Matrices.mat'), 'rho', 'drive_corr', 'feat_names', 'profiles', 'socs');

%% Top 5 features
[~, sort_idx] = sort(abs(rho), 'descend');
top5_idx = sort_idx(1:min(5, numel(sort_idx)));
top5_features = feat_names(top5_idx);
top5_rho = rho(top5_idx);

reportPath = fullfile(reportDir, 'Top5_Features.txt');
fid = fopen(reportPath, 'w');
if fid > 0
    fprintf(fid, 'Top 5 Features by |Correlation| (RPT vs Discharge Capacity)\n');
    for i = 1:numel(top5_idx)
        fprintf(fid, '%d) %s: rho = %.3f\n', i, top5_features{i}, top5_rho(i));
    end
    fclose(fid);
end

%% Additional heatmaps for driving correlation
for s = 1:numel(socs)
    fig = figure('Visible', 'off', 'Position', [100 100 1200 600]);
    imagesc(squeeze(drive_corr(:, :, s))');
    colormap('jet');
    colorbar;
    set(gca, 'YTick', 1:numel(profiles), 'YTickLabel', profiles);
    set(gca, 'XTick', 1:num_feats, 'XTickLabel', feat_names, 'XTickLabelRotation', 90);
    title(sprintf('Driving Feature vs RPT Capacity Correlation (%s)', socs{s}));
    saveas(fig, fullfile(reportDir, sprintf('Driving_Corr_%s.png', socs{s})));
    close(fig);
end

fprintf('Phase 4 complete. Report saved to %s\n', reportDir);

%% Local helper
function y = ifelse(cond, a, b)
if cond
    y = a;
else
    y = b;
end
end
