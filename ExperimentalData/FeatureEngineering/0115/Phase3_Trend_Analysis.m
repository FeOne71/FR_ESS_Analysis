%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 3: Trend Analysis and Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

projectDir = pwd;
resultsDir = fullfile(projectDir, 'Results', 'Phase_3');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

rptFile = fullfile(projectDir, 'Results', 'Phase_1', 'RPT_Feature_Matrix.mat');
driveFile = fullfile(projectDir, 'Results', 'Phase_2', 'Driving_Feature_Matrix.mat');
target_socs = {'SOC70'}; % set {'SOC70','SOC90','SOC50'} as needed

if ~exist(rptFile, 'file')
    error('Missing %s', rptFile);
end
if ~exist(driveFile, 'file')
    error('Missing %s', driveFile);
end

load(rptFile, 'feature_matrix', 'feat_names', 'channels_all', 'rpt_points', 'discharge_ah');
rpt_feature_matrix = feature_matrix;
rpt_feat_names = feat_names;
clear feature_matrix feat_names;

load(driveFile, 'feature_matrix', 'feat_names', 'profiles', 'socs', 'cycles');
if ~isempty(target_socs)
    socs = socs(ismember(socs, target_socs));
end
if isempty(socs)
    error('No matching SOCs found. Check target_socs.');
end
drive_feature_matrix = feature_matrix;
drive_feat_names = feat_names;
clear feature_matrix feat_names;

% Use cycle numbers from file names if possible
rpt_cycles = zeros(1, numel(rpt_points));
for i = 1:numel(rpt_points)
    token = regexp(rpt_points{i}, '(\d+)cyc', 'tokens', 'once');
    if ~isempty(token)
        rpt_cycles(i) = str2double(token{1});
    else
        rpt_cycles(i) = i;
    end
end

%% RPT trend plots and fitting
rptDir = fullfile(resultsDir, 'RPT_Trends');
if ~exist(rptDir, 'dir')
    mkdir(rptDir);
end

rptTrendRows = {};
rowIdx = 1;

for f = 1:numel(rpt_feat_names)
    fig = figure('Visible', 'off', 'Position', [100 100 1200 800]);
    t = tiledlayout(2, 4);
    title(t, sprintf('RPT Trend - %s', rpt_feat_names{f}), 'Interpreter', 'none');

    for c = 1:numel(channels_all)
        nexttile;
        y = squeeze(rpt_feature_matrix(c, :, f));
        plot(rpt_cycles, y, 'o-', 'LineWidth', 1.2);
        hold on;
        [slope, r2, yfit] = fit_linear(rpt_cycles, y);
        plot(rpt_cycles, yfit, '--', 'LineWidth', 1.0);
        grid on;
        xlabel('Cycle');
        ylabel('Feature');
        title(channels_all{c}, 'Interpreter', 'none');

        rptTrendRows(rowIdx, :) = {rpt_feat_names{f}, channels_all{c}, slope, r2, trend_label(slope, y)};
        rowIdx = rowIdx + 1;
    end

    saveas(fig, fullfile(rptDir, sprintf('RPT_Trend_%s.fig', safe_name(rpt_feat_names{f}))));
    close(fig);
end

%% Driving trend plots
driveDir = fullfile(resultsDir, 'Driving_Trends');
if ~exist(driveDir, 'dir')
    mkdir(driveDir);
end

driveTrendRows = {};
rowIdx = 1;

for f = 1:numel(drive_feat_names)
    for s = 1:numel(socs)
        fig = figure('Visible', 'off', 'Position', [100 100 1200 600]);
        hold on;
        for p = 1:numel(profiles)
            y = squeeze(drive_feature_matrix(p, s, :, f));
            plot(cycles, y, 'LineWidth', 1.2, 'DisplayName', profiles{p});
            [slope, r2] = fit_linear(cycles, y);
            driveTrendRows(rowIdx, :) = {drive_feat_names{f}, profiles{p}, socs{s}, slope, r2, trend_label(slope, y)};
            rowIdx = rowIdx + 1;
        end
        grid on;
        xlabel('Cycle');
        ylabel('Feature');
        title(sprintf('Driving Trend - %s (%s)', drive_feat_names{f}, socs{s}), 'Interpreter', 'none');
        legend('Location', 'bestoutside');
        saveas(fig, fullfile(driveDir, sprintf('Driving_Trend_%s_%s.fig', safe_name(drive_feat_names{f}), socs{s})));
        close(fig);
    end
end

%% Sensitivity comparison (max slope per profile)
senseDir = fullfile(resultsDir, 'Driving_Sensitivity');
if ~exist(senseDir, 'dir')
    mkdir(senseDir);
end

for f = 1:numel(drive_feat_names)
    slope_mat = NaN(numel(profiles), numel(socs));
    for p = 1:numel(profiles)
        for s = 1:numel(socs)
            y = squeeze(drive_feature_matrix(p, s, :, f));
            [slope, ~] = fit_linear(cycles, y);
            slope_mat(p, s) = slope;
        end
    end
    slope_mean = mean(abs(slope_mat), 2, 'omitnan');
    fig = figure('Visible', 'off', 'Position', [100 100 800 500]);
    bar(slope_mean);
    set(gca, 'XTickLabel', profiles);
    xlabel('Profile');
    ylabel('Mean |Slope|');
    title(sprintf('Sensitivity - %s', drive_feat_names{f}), 'Interpreter', 'none');
    grid on;
    saveas(fig, fullfile(senseDir, sprintf('Sensitivity_%s.fig', safe_name(drive_feat_names{f}))));
    close(fig);
end

%% Save summary tables
rptTable = cell2table(rptTrendRows, 'VariableNames', ...
    {'Feature', 'Channel', 'Slope', 'R2', 'Trend'});
driveTable = cell2table(driveTrendRows, 'VariableNames', ...
    {'Feature', 'Profile', 'SOC', 'Slope', 'R2', 'Trend'});

summaryPath = fullfile(resultsDir, 'Trend_Summary.xlsx');
writetable(rptTable, summaryPath, 'Sheet', 'RPT');
writetable(driveTable, summaryPath, 'Sheet', 'Driving');

fprintf('Phase 3 complete. Saved to %s\n', resultsDir);

%% Local helpers
function [slope, r2, yfit] = fit_linear(x, y)
    yfit = NaN(size(y));
    slope = NaN;
    r2 = NaN;
    x = x(:); y = y(:);
    valid = isfinite(x) & isfinite(y);
    if sum(valid) < 2
        return;
    end
    p = polyfit(x(valid), y(valid), 1);
    yfit(valid) = polyval(p, x(valid));
    slope = p(1);
    ss_res = sum((y(valid) - yfit(valid)).^2);
    ss_tot = sum((y(valid) - mean(y(valid))).^2);
    if ss_tot > 0
        r2 = 1 - ss_res / ss_tot;
    end
end

function label = trend_label(slope, y)
    label = 'stable';
    if ~isfinite(slope)
        return;
    end
    scale = max(abs(y), [], 'omitnan');
    if isempty(scale) || scale == 0
        scale = 1;
    end
    if abs(slope) < 1e-4 * scale
        label = 'stable';
    elseif slope > 0
        label = 'increasing';
    else
        label = 'decreasing';
    end
end

function name = safe_name(s)
    name = regexprep(s, '[^a-zA-Z0-9_]', '_');
end
