function [V_grid, Q_interp] = local_interp_vq(V, Q, dV)
% interp1(V, Q, V_grid, 'linear'), V_grid = min(V):dV:max(V), no extrapolation
% Handles duplicate V by averaging Q at same V so sample points are unique.
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
Q_interp = interp1(V_unique, Q_unique, V_grid, 'linear');  % no extrap
end

function [V_grid, Q_interp, t_interp] = local_interp_vqt(V, Q, t, dV)
% Same as local_interp_vq but also interpolate time onto V_grid.
% t must be same length as V,Q. Returns (V_grid, Q_interp, t_interp).
V = V(:); Q = Q(:); t = t(:);
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
