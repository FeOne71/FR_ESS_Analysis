function names = feat_power_curve_names()
%FEAT_POWER_CURVE_NAMES Return names for power curve features.
names = cell(1, 101);
for k = 0:100
    names{k+1} = sprintf('P_SOC_%d', k);
end
end
