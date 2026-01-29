function names = feat_11segments_power_names(edges)
%FEAT_11SEGMENTS_POWER_NAMES Return names for 11 segment power features.
if nargin < 1 || isempty(edges)
    edges = linspace(3.0, 4.2, 12);
end
names = cell(1, numel(edges) - 1);
for b = 1:numel(edges) - 1
    names{b} = sprintf('P_%0.2f_%0.2fV', edges(b), edges(b+1));
end
end
