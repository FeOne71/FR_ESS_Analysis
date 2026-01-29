%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 1: RPT Feature Extraction
% - 16D Statistical, 11 Segments, Power Curve
% - Ground truth (discharge Ah, SOH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% Config (update paths if needed)
projectDir = pwd;
resultsDir = fullfile(projectDir, 'Results', 'Phase_1');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
logPath = fullfile(resultsDir, 'phase1_log.txt');

% RPT data location (update as needed)
rptDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';
rptFilePattern = '*.mat';
rptVarName = ''; % optional: set variable name if known

% Parameters
Cnom_Ah = 64;
I_threshold = Cnom_Ah * 0.01;
voltage_window = 0.2;
seg_edges = linspace(3.0, 4.2, 12);
target_channels = {'ch09'}; % set {} or empty to use all channels

%% Load RPT files
rptFiles = dir(fullfile(rptDataDir, rptFilePattern));
if isempty(rptFiles)
    error('No RPT files found in %s', rptDataDir);
end

% Sort files by name
[~, idx] = sort({rptFiles.name});
rptFiles = rptFiles(idx);

%% Prepare containers
feat_names = [feat_stat16_names(), feat_11segments_names(seg_edges), feat_power_curve_names()];
num_feats = numel(feat_names);

% Parse file names (channel, cycle)
fileInfo = struct('path', {}, 'name', {}, 'channel', {}, 'cycle', {}, 'cycle_num', {});
channels_all = {};
cycles_all = {};
for f = 1:numel(rptFiles)
    info = parse_rpt_filename(rptFiles(f).name);
    if isempty(info.channel) || isempty(info.cycle)
        feat_log(logPath, sprintf('Filename not matched: %s', rptFiles(f).name));
        continue;
    end
    info.path = fullfile(rptDataDir, rptFiles(f).name);
    info.name = rptFiles(f).name;
    fileInfo(end+1) = info; %#ok<AGROW>
    channels_all = unique([channels_all; {info.channel}]);
    cycles_all = unique([cycles_all; {info.cycle}]);
end

if isempty(fileInfo)
    error('No RPT files matched expected naming format.');
end

channels_all = sort(channels_all);
if ~isempty(target_channels)
    channels_all = channels_all(ismember(channels_all, target_channels));
end
if isempty(channels_all)
    error('No matching channels found. Check target_channels.');
end
cycle_nums = cellfun(@cycle_to_num, cycles_all);
[~, cyc_idx] = sort(cycle_nums);
cycles_all = cycles_all(cyc_idx);

rpt_points = cycles_all;
feature_matrix = NaN(numel(channels_all), numel(cycles_all), num_feats);
discharge_ah = NaN(numel(channels_all), numel(cycles_all));

%% Main loop
for k = 1:numel(fileInfo)
    info = fileInfo(k);
    [intrim, ok] = load_rpt_intrim(info.path, logPath);
    if ~ok
        continue;
    end
    c = find(strcmp(channels_all, info.channel), 1);
    f = find(strcmp(cycles_all, info.cycle), 1);
    if isempty(c) || isempty(f)
        continue;
    end

    [V, I, t, soc] = extract_ivt_from_intrim(intrim);
    if isempty(V) || isempty(I) || isempty(t)
        feat_log(logPath, sprintf('Missing data: %s in %s', info.channel, info.name));
        continue;
    end
    t = feat_prepare_time(t);

    % Ground truth: discharge capacity from static step (Reference 기준)
    discharge_ah(c, f) = compute_discharge_capacity(intrim, I_threshold, t);
    if ~isfinite(discharge_ah(c, f))
        feat_log(logPath, sprintf('Discharge capacity missing: %s in %s', info.channel, info.name));
    end

    % Method 1: 16D Statistical (charge tail voltage range)
    [V_seg, I_seg, t_seg] = feat_select_segment(V, I, t, 'voltage', voltage_window);
    if isempty(V_seg)
        feat_log(logPath, sprintf('Charge voltage window empty: %s in %s', info.channel, info.name));
    end
    feat16 = feat_stat16(V_seg, I_seg, t_seg);

    % Method 2: 11 Segments (charge)
    feat11 = feat_11segments(V, I, t, seg_edges);

    % Method 3: Power curve (SOC 0-100)
    if isempty(soc)
        soc = feat_estimate_soc(t, I, Cnom_Ah, 0);
    end
    pcurve = feat_power_curve(V, I, t, soc);

    feature_matrix(c, f, :) = [feat16, feat11, pcurve];
end

% SOH (%)
soh_percent = NaN(size(discharge_ah));
if ~isempty(discharge_ah)
    soh_percent = discharge_ah ./ discharge_ah(:, 1) * 100;
end

%% Save
save(fullfile(resultsDir, 'RPT_Feature_Matrix.mat'), ...
    'feature_matrix', 'feat_names', 'channels_all', 'rpt_points', ...
    'discharge_ah', 'soh_percent');

fprintf('Phase 1 complete. Saved to %s\n', resultsDir);

%% Local helpers
function info = parse_rpt_filename(name)
info.channel = '';
info.cycle = '';
info.cycle_num = NaN;
token = regexp(name, 'RPT(\d+)_ch(\d+)_parsed', 'tokens', 'once');
if isempty(token)
    token = regexp(name, 'RPT(\d+)_ch(\d+)', 'tokens', 'once');
end
if ~isempty(token)
    cyc_num = str2double(token{1});
    ch_num = token{2};
    info.cycle = sprintf('%dcyc', cyc_num);
    info.cycle_num = cyc_num;
    info.channel = sprintf('ch%s', ch_num);
end
end

function n = cycle_to_num(cyc)
token = regexp(cyc, '(\d+)cyc', 'tokens', 'once');
if isempty(token)
    n = NaN;
else
    n = str2double(token{1});
end
end

function [intrim, ok] = load_rpt_intrim(filePath, logPath)
ok = false;
intrim = struct();
tmp = load(filePath);
if ~isfield(tmp, 'intrim')
    feat_log(logPath, sprintf('Missing intrim in %s', filePath));
    return;
end
intrim = tmp.intrim;
ok = true;
end

function [V, I, t, soc] = extract_ivt_from_intrim(intrim)
V = []; I = []; t = []; soc = [];
if ~isstruct(intrim)
    return;
end
if isfield(intrim, 'V') && isfield(intrim, 'I') && isfield(intrim, 't')
    V = intrim.V;
    I = intrim.I;
    t = intrim.t;
end
if isfield(intrim, 'SOCDOD')
    soc = intrim.SOCDOD;
end
end

function dch = compute_discharge_capacity(intrim, I_threshold, t)
dch = NaN;
% Reference 기준: StepIdx==3, CycleIdx==2 의 마지막 Capacity(Ah)
if isfield(intrim, 'StepIdx') && isfield(intrim, 'CycleIdx') && isfield(intrim, 'Q')
    idx = intrim.StepIdx == 3 & intrim.CycleIdx == 2;
    if any(idx)
        dch = intrim.Q(find(idx, 1, 'last'));
        return;
    end
end
% Fallback: 마지막 DchgQ (방전 구간)
if isfield(intrim, 'DchgQ') && isfield(intrim, 'type')
    idx = intrim.type == 'D';
    if any(idx)
        dch = intrim.DchgQ(find(idx, 1, 'last'));
        return;
    end
end
% Fallback: Ah-counting
if isfield(intrim, 'I') && isfield(intrim, 't')
    I = intrim.I;
    t = feat_prepare_time(intrim.t);
    idx_dis = I < -I_threshold;
    if sum(idx_dis) >= 2
        dch = abs(trapz(t(idx_dis), I(idx_dis)) / 3600);
    end
end
end
