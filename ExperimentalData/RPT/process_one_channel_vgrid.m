function ch_vq = process_one_channel_vgrid(channel, cycle_key, rpt_cycle, ExperimentalDataPath, Crate_data, static_capacity_step, ocv_steps, crate_labels, dV)
% Build V-grid interpolated (Static, OCV, C-rate) for one channel and one cycle.
% Saves both: (1) grid-aligned V_grid, Q, t  (2) raw V_raw, Q_raw, t_raw
ch_vq = struct();
channel_key = channel;
if strcmp(channel, 'Ch14') && strcmp(rpt_cycle, '800cyc')
    return;
end
filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
filepath = fullfile(ExperimentalDataPath, filename);
if ~isfile(filepath)
    return;
end
try
    T = readtable(filepath, 'VariableNamingRule', 'preserve');
catch
    return;
end
stepCol = 2; subCol = 4; vCol = 8; qCol = 9;
if size(T,2) < 9
    return;
end
% Time, Current, T1(온도) column
tCol = []; iCol = []; t1Col = [];
vn = T.Properties.VariableNames;
for j = 1:numel(vn)
    if isempty(tCol) && contains(lower(char(vn{j})), 'time')
        tCol = j;
    end
    if isempty(iCol) && contains(lower(char(vn{j})), 'current')
        iCol = j;
    end
    if isempty(t1Col) && contains(char(vn{j}), 'T1')
        t1Col = j;
    end
    if ~isempty(tCol) && ~isempty(iCol) && ~isempty(t1Col), break; end
end

% Helper: V,I,Q,t와 동일 마스크로 T1 읽어서 세그먼트에 포함 (원시 T1_raw + V_grid 보간 T1)
save_seg = @(Vg, Qg, tg, Vr, Qr, tr, Ir, T1r) build_seg_with_T1(Vg, Qg, tg, Vr(:), Qr(:), tr, Ir(:), T1r);

% Static (discharge only) — V,I,Q,t,T1 동일 마스크로 한 번에
mask_static = (T{:,stepCol} == static_capacity_step) & (T{:,subCol} == 2);
if sum(mask_static) >= 2
    Vr = T{mask_static, vCol}; Qr = T{mask_static, qCol};
    Ir = get_col(T, mask_static, iCol);
    T1r = get_col(T, mask_static, t1Col);
    tr = get_time_col(T, mask_static, tCol);
    if ~isempty(tr) && numel(tr) == numel(Vr)
        [Vg, Qg, tg] = local_interp_vqt(Vr, Qr, tr, dV);
        ch_vq.Static = save_seg(Vg, Qg, tg, Vr, Qr, tr, Ir, T1r);
    else
        [Vg, Qg] = local_interp_vq(Vr, Qr, dV);
        ch_vq.Static = build_seg_with_T1(Vg, Qg, [], Vr(:), Qr(:), [], Ir(:), T1r);
    end
end
% OCV charge — V,I,Q,t,T1 동일 마스크
mask_chg = (T{:,stepCol} == ocv_steps(1)) & (T{:,subCol} == 2);
mask_dch = (T{:,stepCol} == ocv_steps(2)) & (T{:,subCol} == 2);
if sum(mask_chg) >= 2
    Vr = T{mask_chg, vCol}; Qr = T{mask_chg, qCol};
    Ir = get_col(T, mask_chg, iCol);
    T1r = get_col(T, mask_chg, t1Col);
    tr = get_time_col(T, mask_chg, tCol);
    if ~isempty(tr) && numel(tr) == numel(Vr)
        [Vg, Qg, tg] = local_interp_vqt(Vr, Qr, tr, dV);
        ch_vq.OCV_charge = save_seg(Vg, Qg, tg, Vr, Qr, tr, Ir, T1r);
    else
        [Vg, Qg] = local_interp_vq(Vr, Qr, dV);
        ch_vq.OCV_charge = build_seg_with_T1(Vg, Qg, [], Vr(:), Qr(:), [], Ir(:), T1r);
    end
end
% OCV discharge (flip for ascending V) — T1도 flip하여 V_raw와 순서 일치
if sum(mask_dch) >= 2
    Vr = flipud(T{mask_dch, vCol}(:)); Qr = flipud(T{mask_dch, qCol}(:));
    Ir = get_col(T, mask_dch, iCol);
    if ~isempty(Ir), Ir = flipud(Ir(:)); end
    T1r = get_col(T, mask_dch, t1Col);
    if ~isempty(T1r), T1r = flipud(T1r(:)); end
    tr = get_time_col(T, mask_dch, tCol);
    if ~isempty(tr), tr = flipud(tr(:)); end
    if ~isempty(tr) && numel(tr) == numel(Vr)
        [Vg, Qg, tg] = local_interp_vqt(Vr, Qr, tr, dV);
        ch_vq.OCV_discharge = save_seg(Vg, Qg, tg, Vr, Qr, tr, Ir, T1r);
    else
        [Vg, Qg] = local_interp_vq(Vr, Qr, dV);
        ch_vq.OCV_discharge = build_seg_with_T1(Vg, Qg, [], Vr, Qr, [], Ir, T1r);
    end
end
% C-rate from Crate_data
for r = 1:length(crate_labels)
    label = crate_labels{r};
    % charge
    if isfield(Crate_data, channel_key) && isfield(Crate_data.(channel_key), cycle_key) && ...
       isfield(Crate_data.(channel_key).(cycle_key), label) && ...
       isfield(Crate_data.(channel_key).(cycle_key).(label), 'charge')
        S = Crate_data.(channel_key).(cycle_key).(label).charge;
        V = S.V; Q = S.Q; tr = S.t;
        Ir = get_I_from_crate(S);
        T1r = get_T1_from_crate(S);
        if numel(V) >= 2 && numel(Q) >= 2
            if ~isempty(tr) && numel(tr) == numel(V)
                [Vg, Qg, tg] = local_interp_vqt(V, Q, tr, dV);
                ch_vq.([label '_charge']) = build_seg_with_T1(Vg, Qg, tg, V(:), Q(:), tr, Ir, T1r);
            else
                [Vg, Qg] = local_interp_vq(V, Q, dV);
                ch_vq.([label '_charge']) = build_seg_with_T1(Vg, Qg, [], V(:), Q(:), [], Ir, T1r);
            end
        end
    end
    % discharge
    if isfield(Crate_data, channel_key) && isfield(Crate_data.(channel_key), cycle_key) && ...
       isfield(Crate_data.(channel_key).(cycle_key), label) && ...
       isfield(Crate_data.(channel_key).(cycle_key).(label), 'discharge')
        S = Crate_data.(channel_key).(cycle_key).(label).discharge;
        V = S.V; Q = S.Q; tr = S.t;
        Ir = get_I_from_crate(S);
        T1r = get_T1_from_crate(S);
        if numel(V) >= 2 && numel(Q) >= 2
            if ~isempty(tr) && numel(tr) == numel(V)
                [Vg, Qg, tg] = local_interp_vqt(V, Q, tr, dV);
                ch_vq.([label '_discharge']) = build_seg_with_T1(Vg, Qg, tg, V(:), Q(:), tr, Ir, T1r);
            else
                [Vg, Qg] = local_interp_vq(V, Q, dV);
                ch_vq.([label '_discharge']) = build_seg_with_T1(Vg, Qg, [], V(:), Q(:), [], Ir, T1r);
            end
        end
    end
end
end

function t = get_time_col(T, mask, tCol)
if isempty(tCol), t = []; return; end
try
    t = T{mask, tCol};
catch
    t = [];
end
end

function x = get_col(T, mask, col)
if isempty(col), x = []; return; end
try
    x = T{mask, col};
    if isduration(x), x = double(seconds(x)); end
    x = double(x(:));
catch
    x = [];
end
end

function Ir = get_I_from_crate(S)
if isfield(S, 'I') && ~isempty(S.I)
    Ir = double(S.I(:));
else
    Ir = [];
end
end

function T1r = get_T1_from_crate(S)
if isfield(S, 'T1_raw') && ~isempty(S.T1_raw)
    T1r = double(S.T1_raw(:));
else
    T1r = [];
end
end

function seg = build_seg_with_T1(Vg, Qg, tg, Vr, Qr, tr, Ir, T1r)
seg = struct('V_grid', Vg, 'Q', Qg, 't', tg, ...
    'V_raw', Vr(:), 'Q_raw', Qr(:), 't_raw', tr(:), 'I_raw', Ir(:));
if isempty(T1r)
    seg.T1_raw = [];
    seg.T1 = [];
    return;
end
T1r = double(T1r(:));
seg.T1_raw = T1r;
if numel(Vr) >= 2 && numel(T1r) == numel(Vr) && numel(Vg) >= 2
    V = Vr(:); T1 = T1r;
    [V_sorted, si] = sort(V);
    T1_sorted = T1(si);
    [V_unique, ~, ic] = unique(V_sorted);
    T1_unique = accumarray(ic, T1_sorted, [], @mean);
    if numel(V_unique) >= 2
        seg.T1 = interp1(V_unique, T1_unique, Vg, 'linear');
    else
        seg.T1 = [];
    end
else
    seg.T1 = [];
end
end

%% Local: V,Q (and optional t) interpolated to 0.001V grid (no extrap). Used only inside this file so parfor workers find it.
function [V_grid, Q_interp] = local_interp_vq(V, Q, dV)
V = V(:); Q = Q(:);
if isempty(V) || numel(V) < 2
    V_grid = []; Q_interp = []; return;
end
[V_unique, Q_unique] = unique_vq(V, Q);
if numel(V_unique) < 2
    V_grid = []; Q_interp = []; return;
end
V_min = min(V_unique); V_max = max(V_unique);
V_grid = (V_min : dV : V_max)';
Q_interp = interp1(V_unique, Q_unique, V_grid, 'linear');
end

function [V_grid, Q_interp, t_interp] = local_interp_vqt(V, Q, t, dV)
V = V(:); Q = Q(:);
t = t(:);
if isduration(t)
    t = double(seconds(t));
elseif isdatetime(t)
    t = double(seconds(t - min(t)));  % elapsed seconds
else
    t = double(t);
end
if isempty(V) || numel(V) < 2 || numel(t) ~= numel(V)
    V_grid = []; Q_interp = []; t_interp = []; return;
end
[V_sorted, si] = sort(V);
Q_sorted = Q(si); t_sorted = t(si);
[V_unique, ~, ic] = unique(V_sorted);
Q_unique = accumarray(ic, Q_sorted, [], @mean);
t_unique = accumarray(ic, t_sorted, [], @mean);
if numel(V_unique) < 2
    V_grid = []; Q_interp = []; t_interp = []; return;
end
V_min = min(V_unique); V_max = max(V_unique);
V_grid = (V_min : dV : V_max)';
Q_interp = interp1(V_unique, Q_unique, V_grid, 'linear');
t_interp = interp1(V_unique, t_unique, V_grid, 'linear');
end

function [V_unique, Q_unique] = unique_vq(V, Q)
[V_sorted, si] = sort(V);
Q_sorted = Q(si);
[V_unique, ~, ic] = unique(V_sorted);
Q_unique = accumarray(ic, Q_sorted, [], @mean);
end
