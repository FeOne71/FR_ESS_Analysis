function [q_bins, bin_names] = feat_11segments(V, I, t, edges)
%FEAT_11SEGMENTS Compute charge Ah per voltage bin.
% edges should be 12 values for 11 bins.

if nargin < 4 || isempty(edges)
    edges = linspace(3.0, 4.2, 12);
end
bin_names = feat_11segments_names(edges);
q_bins = NaN(1, numel(edges) - 1);

if isempty(V) || isempty(I) || isempty(t)
    return;
end
V = V(:); I = I(:); t = t(:);
if numel(V) < 2 || numel(I) < 2 || numel(t) < 2
    return;
end

for b = 1:numel(edges) - 1
    idx = V >= edges(b) & V < edges(b+1) & I > 0;
    if sum(idx) >= 2
        q_bins(b) = trapz(t(idx), I(idx)) / 3600;
    else
        q_bins(b) = NaN;
    end
end
end
