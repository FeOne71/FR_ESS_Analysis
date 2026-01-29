function [feat_vec, feat_names] = feat_stat16(V, I, t)
%FEAT_STAT16 Compute 16 statistical features from V/I/t.
% Features: meanV,stdV,skewV,kurtV,meanI,stdI,skewI,kurtI,
% mean_dVdt,std_dVdt,duration_s,Ah, slopeV, entropyV, entropyI, energyWh.

feat_names = feat_stat16_names();
feat_vec = NaN(1, numel(feat_names));

if isempty(V) || isempty(I) || isempty(t)
    return;
end
V = V(:); I = I(:); t = t(:);
if numel(V) < 2 || numel(I) < 2 || numel(t) < 2
    return;
end

% Basic stats
feat_vec(1) = mean(V, 'omitnan');
feat_vec(2) = std(V, 'omitnan');
feat_vec(3) = nan_skewness(V);
feat_vec(4) = nan_kurtosis(V);
feat_vec(5) = mean(I, 'omitnan');
feat_vec(6) = std(I, 'omitnan');
feat_vec(7) = nan_skewness(I);
feat_vec(8) = nan_kurtosis(I);

% Derivatives
dt = diff(t);
dt(dt <= 0 | ~isfinite(dt)) = NaN;
if all(isnan(dt))
    return;
end
dVdt = diff(V) ./ dt;
feat_vec(9) = mean(dVdt, 'omitnan');
feat_vec(10) = std(dVdt, 'omitnan');

% Duration
feat_vec(11) = max(t) - min(t);

% Ah and energy
feat_vec(12) = trapz(t, I) / 3600;
feat_vec(16) = trapz(t, V .* I) / 3600;

% Slope (linear fit of V vs t)
valid_idx = isfinite(V) & isfinite(t);
if sum(valid_idx) >= 2
    p = polyfit(t(valid_idx), V(valid_idx), 1);
    feat_vec(13) = p(1);
end

% Entropy
feat_vec(14) = feat_entropy(V, 20);
feat_vec(15) = feat_entropy(I, 20);
end

function s = nan_skewness(x)
    x = x(:);
    x = x(isfinite(x));
    if numel(x) < 3
        s = NaN;
        return;
    end
    s = skewness(x, 0);
end

function k = nan_kurtosis(x)
    x = x(:);
    x = x(isfinite(x));
    if numel(x) < 4
        k = NaN;
        return;
    end
    k = kurtosis(x, 0);
end
