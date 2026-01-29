function soc = feat_estimate_soc(t, I, c_nom_ah, initial_soc)
%FEAT_ESTIMATE_SOC Estimate SOC (%) using coulomb counting.
if nargin < 4 || isempty(initial_soc)
    initial_soc = 0;
end
soc = NaN(size(t));
if isempty(t) || isempty(I) || numel(t) < 2 || numel(I) ~= numel(t)
    return;
end
t = t(:); I = I(:);
dt = [0; diff(t)];
dt(dt < 0 | ~isfinite(dt)) = 0;
ah = cumsum(I .* dt) / 3600;
soc = initial_soc + (ah / c_nom_ah) * 100;
soc = max(0, min(100, soc));
end
