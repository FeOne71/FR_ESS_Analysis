function H = feat_entropy(x, num_bins)
%FEAT_ENTROPY Compute histogram-based entropy.
if nargin < 2 || isempty(num_bins)
    num_bins = 20;
end
H = NaN;
if isempty(x) || all(isnan(x))
    return;
end
x = x(:);
x = x(isfinite(x));
if numel(x) < 2
    return;
end
counts = histcounts(x, num_bins);
if isempty(counts) || sum(counts) == 0
    return;
end
p = counts / sum(counts);
p = p(p > 0);
H = -sum(p .* log(p));
end
