function [p_curve, bin_names] = feat_power_curve(V, I, t, soc)
%FEAT_POWER_CURVE Compute mean power by SOC bin (0-100%).
bin_names = feat_power_curve_names();
p_curve = NaN(1, 101);

if isempty(V) || isempty(I) || isempty(t)
    return;
end
V = V(:); I = I(:); t = t(:);
if numel(V) < 2 || numel(I) < 2 || numel(t) < 2
    return;
end

if nargin < 4 || isempty(soc)
    soc = [];
else
    soc = soc(:);
end

if isempty(soc) || numel(soc) ~= numel(t)
    return;
end

power = V .* I;
edges = -0.5:1:100.5;
[~, ~, bin] = histcounts(soc, edges);
for k = 1:101
    idx = bin == k;
    if any(idx)
        p_curve(k) = mean(power(idx), 'omitnan');
    end
end
end
