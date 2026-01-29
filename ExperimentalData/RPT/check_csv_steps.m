% Check which steps exist in Ch15 RPT200cyc CSV
clear; clc;

csvFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data\RPT\Ch15_RPT_200cyc.csv';
fprintf('Checking: %s\n', csvFile);

T = readmatrix(csvFile);

% Step Index is column 2
unique_steps = unique(T(:,2));
fprintf('\nUnique Step Indices in CSV: %s\n', mat2str(sort(unique_steps)));

% Check C-rate steps
crate_steps_charge = [28 32 36 40 44];
crate_steps_discharge = [48 52 56 60 64];

fprintf('\n=== Checking C-rate Steps ===\n');
fprintf('Charge steps (28, 32, 36, 40, 44):\n');
for i = 1:length(crate_steps_charge)
    step = crate_steps_charge(i);
    mask = (T(:,2) == step);
    count = sum(mask);
    if count > 0
        fprintf('  Step %d: %d data points\n', step, count);
        fprintf('    Voltage range: [%.3f, %.3f] V\n', min(T(mask, 8)), max(T(mask, 8)));
        fprintf('    Capacity range: [%.3f, %.3f] Ah\n', min(T(mask, 9)), max(T(mask, 9)));
    else
        fprintf('  Step %d: NO DATA\n', step);
    end
end

fprintf('\nDischarge steps (48, 52, 56, 60, 64):\n');
for i = 1:length(crate_steps_discharge)
    step = crate_steps_discharge(i);
    mask = (T(:,2) == step);
    count = sum(mask);
    if count > 0
        fprintf('  Step %d: %d data points\n', step, count);
        fprintf('    Voltage range: [%.3f, %.3f] V\n', min(T(mask, 8)), max(T(mask, 8)));
        fprintf('    Capacity range: [%.3f, %.3f] Ah\n', min(T(mask, 9)), max(T(mask, 9)));
    else
        fprintf('  Step %d: NO DATA\n', step);
    end
end

fprintf('\n=== Check Complete ===\n');
