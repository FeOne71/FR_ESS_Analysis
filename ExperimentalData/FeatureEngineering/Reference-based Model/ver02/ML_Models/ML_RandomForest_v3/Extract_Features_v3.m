% Extract_Features_v3.m  — [3.4, 4.0]V 범위, 11 dQ 세그먼트 피처 추출
%
% [v3] 기존 RPT_FeatureExtractor.m 대비 변경사항:
%   1. 전압 경계: MasterRuler → linspace(3.4, 4.0, 12) (11 세그먼트, ~55mV/seg)
%   2. 피처: dQ_chg_1~11 + dQ_dch_1~11 + C_eff_chg + C_eff_dch = 24개
%      (Peak, Energy, T_avg 제거)
%   3. 충전/방전 경계 동일 (대칭 경계)
%   4. 각 세그먼트의 유효성(valid 여부) 기록 → NaN 허용 (RF surrogate split 활용)
%
% 데이터 소스: RPT_VQ_grid.mat (pre-gridded V-Q curves at 0.001V resolution)

clear; clc; close all; warning off;

%% Configuration
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_rpt_vq = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
path_static = fullfile(baseDir, 'Capacity_Trend_Figures', 'Capacity_Data_Static.mat');

ver02Dir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
saveDir  = fullfile(ver02Dir, 'ML_RandomForest_v3');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

fprintf('Loading Data...\n');
load(path_rpt_vq, 'RPT_VQ_grid');
load(path_static, 'allChannelsCapacity');

channels      = fieldnames(allChannelsCapacity);
target_crates = {'c01', 'c05', 'c1', 'c2', 'c3'};
target_crates_val = [0.1, 0.5, 1, 2, 3];

% === 핵심: 새로운 전압 경계 ===
V_bounds = linspace(3.4, 4.0, 12);  % 12개 경계 → 11개 세그먼트 (~54.5mV 간격)
num_segments = length(V_bounds) - 1;
fprintf('V_bounds: %s\n', mat2str(round(V_bounds, 4)));
fprintf('Segments: %d (each ~%.1f mV)\n', num_segments, diff(V_bounds(1:2))*1000);

%% Feature Extraction Loop
fprintf('\nExtracting %d dQ Features...\n', num_segments * 2 + 2);
data_CellID     = {};
data_Cycle      = [];
data_CrateLabel = {};
data_CrateNum   = [];
data_X          = [];  % [dQ_chg(11), dQ_dch(11), C_eff_chg, C_eff_dch]

cnt = 0;
for i = 1:length(channels)
    ch = channels{i};
    cap_data = allChannelsCapacity.(ch);
    idx_cyc0 = find(cap_data.cycles == 0, 1);
    if isempty(idx_cyc0), idx_cyc0 = 1; end

    Q_0 = NaN;
    item = cap_data.Q{1, idx_cyc0};
    if ~isempty(item), Q_0 = max(item); end
    if isnan(Q_0), warning('Q_0 not found for %s, skipping.', ch); continue; end

    cyc_fields = fieldnames(RPT_VQ_grid);
    for c = 1:length(cyc_fields)
        cyc_key = cyc_fields{c};
        cyc_num = sscanf(cyc_key, 'cyc%d');
        if ~isfield(RPT_VQ_grid.(cyc_key), ch), continue; end
        ch_data = RPT_VQ_grid.(cyc_key).(ch);

        for r = 1:length(target_crates)
            crate_label = target_crates{r};
            crate_val   = target_crates_val(r);
            f_chg = [crate_label '_charge'];
            f_dch = [crate_label '_discharge'];
            if ~isfield(ch_data, f_chg) || ~isfield(ch_data, f_dch), continue; end

            s_chg = ch_data.(f_chg);
            s_dch = ch_data.(f_dch);

            % --- dQ 추출: 충전 (V_bounds 기준, 대칭 경계) ---
            dQ_chg = nan(1, num_segments);
            if length(s_chg.Q) > 1
                q_min = min(s_chg.Q); q_max = max(s_chg.Q);
                q_grid = q_min : 0.01 : q_max;
                [Q_u, uid] = unique(s_chg.Q); V_u = s_chg.V_grid(uid);
                V_curve = interp1(Q_u, V_u, q_grid, 'linear', 'extrap');

                V_range_c = [min(s_chg.V_grid), max(s_chg.V_grid)];
                QR = nan(1, length(V_bounds));
                for b = 1:length(V_bounds)
                    % 경계가 실제 V 범위 내에 있을 때만 유효
                    if V_bounds(b) >= V_range_c(1) - 0.01 && V_bounds(b) <= V_range_c(2) + 0.01
                        [~, idx_v] = min(abs(V_curve - V_bounds(b)));
                        QR(b) = q_grid(idx_v);
                    end
                end
                dQ_chg = abs(diff(QR));  % NaN 전파: 양쪽 경계 중 하나라도 없으면 NaN
            end

            % --- dQ 추출: 방전 (동일 V_bounds, 역방향) ---
            dQ_dch = nan(1, num_segments);
            if length(s_dch.Q) > 1
                q_min = min(s_dch.Q); q_max = max(s_dch.Q);
                q_grid = q_min : 0.01 : q_max;
                [Q_u, uid] = unique(s_dch.Q); V_u = s_dch.V_grid(uid);
                V_curve = interp1(Q_u, V_u, q_grid, 'linear', 'extrap');

                V_range_d = [min(s_dch.V_grid), max(s_dch.V_grid)];
                QR = nan(1, length(V_bounds));
                for b = 1:length(V_bounds)
                    if V_bounds(b) >= V_range_d(1) - 0.01 && V_bounds(b) <= V_range_d(2) + 0.01
                        [~, idx_v] = min(abs(V_curve - V_bounds(b)));
                        QR(b) = q_grid(idx_v);
                    end
                end
                dQ_dch = abs(diff(QR));
            end

            % --- C_eff ---
            C_eff_chg = mean(abs(s_chg.I_raw)) / Q_0;
            C_eff_dch = mean(abs(s_dch.I_raw)) / Q_0;

            % --- Collect ---
            cnt = cnt + 1;
            data_CellID{cnt,1}     = ch;
            data_Cycle(cnt,1)      = cyc_num;
            data_CrateLabel{cnt,1} = crate_label;
            data_CrateNum(cnt,1)   = crate_val;
            data_X(cnt,:) = [dQ_chg, dQ_dch, C_eff_chg, C_eff_dch];
        end
    end
end

%% Save
% 피처 이름 생성
feature_names = {};
for s = 1:num_segments
    feature_names{end+1} = sprintf('dQ_chg_%d', s); %#ok<SAGROW>
end
for s = 1:num_segments
    feature_names{end+1} = sprintf('dQ_dch_%d', s); %#ok<SAGROW>
end
feature_names{end+1} = 'C_eff_chg';
feature_names{end+1} = 'C_eff_dch';

FeatureTable_v3 = table(data_CellID, data_Cycle, data_CrateLabel, data_CrateNum, ...
    'VariableNames', {'CellID','Cycle','CrateLabel','CrateNum'});
for f = 1:length(feature_names)
    FeatureTable_v3.(feature_names{f}) = data_X(:, f);
end

save(fullfile(saveDir, 'Feature_Matrix_v3.mat'), 'FeatureTable_v3', 'feature_names', 'V_bounds');

fprintf('\n=== Feature Extraction v3 Complete ===\n');
fprintf('  Samples: %d | Features: %d | Segments: %d\n', cnt, length(feature_names), num_segments);
fprintf('  V range: [%.2f, %.2f] V\n', V_bounds(1), V_bounds(end));
fprintf('  NaN ratio: %.1f%%\n', 100*sum(isnan(data_X(:)))/(numel(data_X)));
fprintf('  Saved to: %s\n', saveDir);
