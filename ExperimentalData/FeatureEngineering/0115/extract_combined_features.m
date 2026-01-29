function [feat16, q11, p11, p_count] = extract_combined_features(V, I, t, seg_edges, voltage_window)
%EXTRACT_COMBINED_FEATURES Return [16D, 11Q, 11P] for one event.
% - 16D stats from charge tail voltage window
% - 11Q: charge Ah per voltage segment
% - 11P: mean power per voltage segment

if nargin < 4 || isempty(seg_edges)
    seg_edges = linspace(3.0, 4.2, 12);
end
if nargin < 5 || isempty(voltage_window)
    voltage_window = 0.2;
end

feat16 = NaN(1, 16);
q11 = NaN(1, numel(seg_edges) - 1);
p11 = NaN(1, numel(seg_edges) - 1);
p_count = zeros(1, numel(seg_edges) - 1);

if isempty(V) || isempty(I) || isempty(t)
    return;
end
V = V(:); I = I(:); t = t(:);
if numel(V) < 2 || numel(I) < 2 || numel(t) < 2
    return;
end

% 16D stats: charge tail voltage window (fallback to tail 20%)
[V_seg, I_seg, t_seg] = feat_select_segment(V, I, t, 'voltage', voltage_window);
if isempty(V_seg)
    [V_seg, I_seg, t_seg] = feat_select_segment(V, I, t, 'tail', 0.2);
end
feat16 = feat_stat16(V_seg, I_seg, t_seg);

% 11 segments: charge Ah
q11 = feat_11segments(V, I, t, seg_edges);

% 11 segments: mean power per segment
power = V .* I;
for b = 1:numel(seg_edges) - 1
    idx = V >= seg_edges(b) & V < seg_edges(b+1) & I > 0;
    if any(idx)
        p11(b) = mean(power(idx), 'omitnan');
        p_count(b) = sum(idx);
    end
end
end
