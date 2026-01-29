%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 2: Driving Feature Extraction
% - Uses Lab_DC_Events_Raw_*.mat (charge events)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% Config
projectDir = pwd;
resultsDir = fullfile(projectDir, 'Results', 'Phase_2');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
logPath = fullfile(resultsDir, 'phase2_log.txt');

eventDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\ver03_OnlyChargeEvents\Results';
eventPattern = 'Lab_DC_Events_Raw_*.mat';

target_channel = 'ch9'; % reference channel
seg_edges = linspace(3.0, 4.2, 12);
voltage_window = 0.2;

%% Load event files
eventFiles = dir(fullfile(eventDir, eventPattern));
if isempty(eventFiles)
    error('No event files found in %s', eventDir);
end

% Sort by cycle number if possible
cycle_nums = zeros(numel(eventFiles), 1);
for f = 1:numel(eventFiles)
    token = regexp(eventFiles(f).name, 'Lab_DC_Events_Raw_(\d+)cyc', 'tokens', 'once');
    if ~isempty(token)
        cycle_nums(f) = str2double(token{1});
    else
        cycle_nums(f) = 999999;
    end
end
[~, idx] = sort(cycle_nums);
eventFiles = eventFiles(idx);
cycles = cycle_nums(idx);

profiles = {'DC1','DC2','DC3','DC4','DC5','DC6','DC7','DC8'};
socs = {'SOC90','SOC70','SOC50'};

feat_names = [feat_stat16_names(), feat_11segments_names(seg_edges), feat_11segments_power_names(seg_edges)];
num_feats = numel(feat_names);

feature_matrix = NaN(numel(profiles), numel(socs), numel(eventFiles), num_feats);

%% Main loop
for f = 1:numel(eventFiles)
    filePath = fullfile(eventDir, eventFiles(f).name);
    data = load(filePath, 'rawEvents');
    if ~isfield(data, 'rawEvents')
        feat_log(logPath, sprintf('Missing rawEvents in %s', eventFiles(f).name));
        continue;
    end
    rawEvents = data.rawEvents;

    % select channel struct name
    structName = '';
    allStructs = fieldnames(rawEvents);
    for s = 1:numel(allStructs)
        if startsWith(allStructs{s}, target_channel) && endsWith(allStructs{s}, '_Charge')
            structName = allStructs{s};
            break;
        end
    end
    if isempty(structName)
        % fallback: any charge struct
        charge_structs = allStructs(contains(allStructs, '_Charge'));
        if isempty(charge_structs)
            feat_log(logPath, sprintf('No charge events in %s', eventFiles(f).name));
            continue;
        end
        structName = charge_structs{1};
    end

    for s = 1:numel(socs)
        socName = socs{s};
        if ~isfield(rawEvents.(structName), socName)
            continue;
        end
        for p = 1:numel(profiles)
            profName = profiles{p};
            if ~isfield(rawEvents.(structName).(socName), profName)
                continue;
            end
            events = rawEvents.(structName).(socName).(profName);
            evtNames = fieldnames(events);
            if isempty(evtNames)
                continue;
            end

            % Method 1: 16D stats per event, then average
            feat16_all = NaN(numel(evtNames), 16);
            % Method 2: 11 segments accumulated
            q_acc = zeros(1, numel(seg_edges) - 1);
            q_valid = false(1, numel(seg_edges) - 1);
            % Method 3: segment power aggregated (weighted by sample count)
            p_sum = zeros(1, numel(seg_edges) - 1);
            p_count = zeros(1, numel(seg_edges) - 1);

            for e = 1:numel(evtNames)
                evt = events.(evtNames{e});
                if ~isfield(evt, 'V') || ~isfield(evt, 'I') || ~isfield(evt, 't')
                    continue;
                end
                V = evt.V; I = evt.I; t = feat_prepare_time(evt.t);
                if numel(V) < 2 || numel(I) < 2 || numel(t) < 2
                    continue;
                end

                if (max(V) - min(V)) < 0.2
                    feat_log(logPath, sprintf('Skip short event (dV<0.2): %s %s %s %s', ...
                        eventFiles(f).name, socName, profName, evtNames{e}));
                    continue;
                end

                [feat16_evt, q_evt, p_evt, p_evt_count] = extract_combined_features( ...
                    V, I, t, seg_edges, voltage_window);

                feat16_all(e, :) = feat16_evt;

                if any(~isnan(q_evt))
                    q_acc = q_acc + nan_to_zero(q_evt);
                    q_valid = q_valid | ~isnan(q_evt);
                end

                if any(p_evt_count > 0)
                    p_sum = p_sum + nan_to_zero(p_evt) .* p_evt_count;
                    p_count = p_count + p_evt_count;
                end
            end

            feat16 = mean(feat16_all, 1, 'omitnan');
            feat11 = NaN(1, numel(seg_edges) - 1);
            feat11(q_valid) = q_acc(q_valid);
            pseg = NaN(1, numel(seg_edges) - 1);
            valid_p = p_count > 0;
            pseg(valid_p) = p_sum(valid_p) ./ p_count(valid_p);

            feature_matrix(p, s, f, :) = [feat16, feat11, pseg];
        end
    end
end

%% Save
save(fullfile(resultsDir, 'Driving_Feature_Matrix.mat'), ...
    'feature_matrix', 'feat_names', 'profiles', 'socs', 'cycles');

fprintf('Phase 2 complete. Saved to %s\n', resultsDir);

%% Local helper
function x = nan_to_zero(x)
    x(isnan(x)) = 0;
end
