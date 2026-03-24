% Visualize_Common_Voltage_Range.m
% 필드 데이터 (2021~2025)의 충전/방전 구간에 대해 용량(Q) - 전압(V) 상관관계와 공통 전압 부분을 시각화합니다.
% RPT_Field_Estimation_RF.m 의 데이터 로드 로직을 참조합니다.

clear; clc; close all;
warning on;

fprintf('=== Loading Field Data (2021~2025) ===\n');
baseDir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
rawDir  = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'Rack_raw2mat');
dataFiles = {
    fullfile(rawDir, 'Old', '2021', '202106', 'Raw_20210603.mat'), 'Y2021', 'old',  datetime(2021,6,3);
    fullfile(rawDir, 'New', '2023', '202310', 'Raw_20231016.mat'), 'Y2023', 'new',  datetime(2023,10,16);
    fullfile(rawDir, 'New', '2024', '202409', 'Raw_20240909.mat'), 'Y2024', 'new',  datetime(2024,9,9);
    % fullfile(rawDir, 'New', '2025', '202507', 'Raw_20250711.mat'), 'Y2025', 'new',  datetime(2025,7,11);
};

min_charge_secs    = [600, 300, 300, 300];   % Y2021,Y2023,Y2024,Y2025
min_discharge_secs = [300, 150, 300, 150];
dt = 1;
Np = 2; C_cell_Ah = 64; thr_A = C_cell_Ah * 0.05 / Np;

FieldData = struct();
for k = 1:size(dataFiles, 1)
    fpath      = dataFiles{k, 1};
    year_label = dataFiles{k, 2};
    dataType   = dataFiles{k, 3};
    base_date  = dataFiles{k, 4};
    min_chg_sec = min_charge_secs(k);
    min_dch_sec = min_discharge_secs(k);

    if ~exist(fpath, 'file'), continue; end
    fprintf('  Processing %s...\n', year_label);

    S = load(fpath);
    if strcmp(dataType, 'old')
        D = S.Raw.Rack01;
        if isfield(D, 'Time'), t = datetime(D.Time);
        elseif isfield(D, 'Date_Time')
            if isduration(D.Date_Time), t = base_date + D.Date_Time;
            else, t = datetime(D.Date_Time); end
        else, continue; end
        I_rack = D.DCCurrent_A(:);
        V_avg  = D.AverageCV_V(:);
    else
        D = S.Raw;
        if isduration(D.Date_Time), t = base_date + D.Date_Time;
        else, t = datetime(D.Date_Time); end
        I_rack = D.DCCurrent(:);
        V_avg  = D.CVavg(:);
    end

    tsec   = seconds(t - t(1));
    I_cell = I_rack / Np;

    isChg   = I_cell >  thr_A;
    isDch   = I_cell < -thr_A;
    chgSegs = local_find_segments(isChg);
    dchSegs = local_find_segments(isDch);

    if ~isempty(chgSegs)
        dur_c = chgSegs(:,2) - chgSegs(:,1) + 1;
        chgSegs = chgSegs(dur_c >= ceil(min_chg_sec/dt), :);
    end
    if ~isempty(dchSegs)
        dur_d = dchSegs(:,2) - dchSegs(:,1) + 1;
        dchSegs = dchSegs(dur_d >= ceil(min_dch_sec/dt), :);
    end

    if strcmp(year_label, 'Y2024')
        target_chg_start = datetime(2024, 9, 9, 12, 55, 0);
        target_chg_end   = datetime(2024, 9, 9, 13, 13, 0);
        chg_seg_idx = [];
        for si = 1:size(chgSegs, 1)
            if abs(t(chgSegs(si,1)) - target_chg_start) < minutes(5) && ...
               abs(t(chgSegs(si,2)) - target_chg_end)   < minutes(5)
                chg_seg_idx = si; break;
            end
        end
        if isempty(chg_seg_idx), chgSegs = [];
        else, chgSegs = chgSegs(chg_seg_idx, :); end

        target_dch_start = datetime(2024, 9, 9, 14, 14, 0);
        target_rest_end  = datetime(2024, 9, 9, 14, 25, 49);
        dch_seg_idx = [];
        for si = 1:size(dchSegs, 1)
            if abs(t(dchSegs(si,1)) - target_dch_start) < minutes(1)
                dch_seg_idx = si; break;
            end
        end
        if ~isempty(dch_seg_idx)
            dch_end = dchSegs(dch_seg_idx, 2);
            if t(dch_end) > target_rest_end
                idx_limit = find(t <= target_rest_end, 1, 'last');
                if ~isempty(idx_limit) && idx_limit >= dchSegs(dch_seg_idx, 1)
                    dchSegs(dch_seg_idx, 2) = idx_limit;
                end
            end
            dchSegs = dchSegs(dch_seg_idx, :);
        else
            dchSegs = [];
        end
    end

    if ~isempty(chgSegs)
        [~, best] = max(chgSegs(:,2) - chgSegs(:,1));
        chg_s = chgSegs(best, 1);  chg_e = chgSegs(best, 2);
        FieldData.(year_label).Chg.V = V_avg(chg_s:chg_e);
        FieldData.(year_label).Chg.I = I_cell(chg_s:chg_e);
        FieldData.(year_label).Chg.Q = cumtrapz(abs(I_cell(chg_s:chg_e))) / 3600;
    end
    if ~isempty(dchSegs)
        [~, best] = max(dchSegs(:,2) - dchSegs(:,1));
        dch_s = dchSegs(best, 1);  dch_e = dchSegs(best, 2);
        FieldData.(year_label).Dch.V = V_avg(dch_s:dch_e);
        FieldData.(year_label).Dch.I = I_cell(dch_s:dch_e);
        FieldData.(year_label).Dch.Q = cumtrapz(abs(I_cell(dch_s:dch_e))) / 3600;
    end
end

%% Calculate and Plot Common Voltage Range (V vs Q)
years = fieldnames(FieldData);
n_yrs = length(years);

chg_mins = zeros(1, n_yrs); chg_maxs = zeros(1, n_yrs);
dch_mins = zeros(1, n_yrs); dch_maxs = zeros(1, n_yrs);

for k = 1:n_yrs
    yr = years{k};
    if isfield(FieldData.(yr), 'Chg')
        chg_mins(k) = min(FieldData.(yr).Chg.V);
        chg_maxs(k) = max(FieldData.(yr).Chg.V);
    else
        chg_mins(k) = NaN; chg_maxs(k) = NaN;
    end
    if isfield(FieldData.(yr), 'Dch')
        dch_mins(k) = min(FieldData.(yr).Dch.V);
        dch_maxs(k) = max(FieldData.(yr).Dch.V);
    else
        dch_mins(k) = NaN; dch_maxs(k) = NaN;
    end
end

% Common ranges across all valid years
common_chg_min = max(chg_mins);
common_chg_max = min(chg_maxs);
common_dch_min = max(dch_mins);
common_dch_max = min(dch_maxs);

fprintf('\n=== Common Voltage Range ===\n');
fprintf('Charge   : [%.4f, %.4f] V\n', common_chg_min, common_chg_max);
fprintf('Discharge: [%.4f, %.4f] V\n', common_dch_min, common_dch_max);

fig = figure('Name', 'Common Voltage Range (V vs Q)', 'Position', [150, 150, 1000, 500]);
sgtitle('Field Data Common Voltage Range (2021-2025)', 'FontSize', 15, 'FontWeight', 'bold');

colors = lines(n_yrs);

% Subplot 1: Charge
subplot(1,2,1); hold on; box on; grid on;
for k = 1:n_yrs
    yr = years{k};
    if isfield(FieldData.(yr), 'Chg')
        plot(FieldData.(yr).Chg.Q, FieldData.(yr).Chg.V, '-', 'LineWidth', 2, 'Color', colors(k,:), 'DisplayName', yr);
    end
end
if common_chg_max > common_chg_min
    xl = xlim; 
    patch([xl(1) xl(2) xl(2) xl(1)], ...
        [common_chg_min common_chg_min common_chg_max common_chg_max], ...
        'r', 'FaceAlpha', 0.15, 'EdgeColor', 'r', 'LineStyle', '--', 'DisplayName', 'Common Region');
    str_chg = sprintf('Common: %.3f ~ %.3f V', common_chg_min, common_chg_max);
    text(mean(xl), common_chg_max - 0.02, str_chg, 'HorizontalAlignment','center', 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');
end
xlabel('Capacity (Ah)'); ylabel('Voltage (V)'); title('Charge Profile (V vs Q)');
legend('Location', 'best');

% Subplot 2: Discharge
subplot(1,2,2); hold on; box on; grid on;
for k = 1:n_yrs
    yr = years{k};
    if isfield(FieldData.(yr), 'Dch')
        plot(FieldData.(yr).Dch.Q, FieldData.(yr).Dch.V, '-', 'LineWidth', 2, 'Color', colors(k,:), 'DisplayName', yr);
    end
end
if common_dch_max > common_dch_min
    xl = xlim;
    patch([xl(1) xl(2) xl(2) xl(1)], ...
        [common_dch_min common_dch_min common_dch_max common_dch_max], ...
        'b', 'FaceAlpha', 0.15, 'EdgeColor', 'b', 'LineStyle', '--', 'DisplayName', 'Common Region');
    str_dch = sprintf('Common: %.3f ~ %.3f V', common_dch_min, common_dch_max);
    text(mean(xl), common_dch_max - 0.02, str_dch, 'HorizontalAlignment','center', 'FontSize', 11, 'Color', 'b', 'FontWeight', 'bold');
end
xlabel('Capacity (Ah)'); ylabel('Voltage (V)'); title('Discharge Profile (V vs Q)');
legend('Location', 'best');

% Save figure
savePath = fullfile(baseDir, 'Common_Voltage_Range_Field_V_vs_Q.fig');
saveas(fig, savePath);
fprintf('Figure saved at: %s\n', savePath);

% Helper function
function segs = local_find_segments(mask)
    segs = [];
    n = length(mask);
    i = 1;
    while i <= n
        if mask(i)
            j = i;
            while j < n && mask(j+1), j = j + 1; end
            segs = [segs; i, j];
            i = j + 1;
        else
            i = i + 1;
        end
    end
end
