function [V_seg, I_seg, t_seg] = feat_select_segment(V, I, t, mode, varargin)
%FEAT_SELECT_SEGMENT Select segment for feature extraction.
% mode: 'voltage' -> use [Vmax-window, Vmax], 'tail' -> last pct of time.

V_seg = [];
I_seg = [];
t_seg = [];
if isempty(V) || isempty(I) || isempty(t)
    return;
end
V = V(:); I = I(:); t = t(:);

switch lower(mode)
    case 'voltage'
        window = 0.2;
        if ~isempty(varargin)
            window = varargin{1};
        end
        idx = I > 0;
        if ~any(idx)
            return;
        end
        V_charge = V(idx);
        Vmax = max(V_charge);
        v_low = Vmax - window;
        idx2 = idx & V >= v_low & V <= Vmax;
        if sum(idx2) >= 2
            V_seg = V(idx2); I_seg = I(idx2); t_seg = t(idx2);
        end
    case 'tail'
        pct = 0.2;
        if ~isempty(varargin)
            pct = varargin{1};
        end
        n = numel(t);
        if n < 2
            return;
        end
        start_idx = max(1, floor(n * (1 - pct)) + 1);
        V_seg = V(start_idx:end);
        I_seg = I(start_idx:end);
        t_seg = t(start_idx:end);
    otherwise
        return;
end
end
