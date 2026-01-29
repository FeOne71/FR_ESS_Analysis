% Check all steps in Ch15 RPT200cyc to find C-rate data
clear; clc;

csvFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data\RPT\Ch15_RPT_200cyc.csv';
fprintf('Checking all steps in: %s\n', csvFile);

T = readmatrix(csvFile);

% Step Index (col 2), Cycle Index (col 4), Current (col 7), Voltage (col 8), Capacity (col 9)
unique_steps = unique(T(:,2));
fprintf('\nAll Step Indices: %s\n', mat2str(sort(unique_steps)));

% Check each step for current and capacity data
fprintf('\n=== Step Details (with Current and Capacity) ===\n');
for i = 1:length(unique_steps)
    step = unique_steps(i);
    mask = (T(:,2) == step);
    if sum(mask) > 0
        I_data = T(mask, 7);
        V_data = T(mask, 8);
        Q_data = T(mask, 9);
        cycle_idx = unique(T(mask, 4));
        
        fprintf('\nStep %d: %d points, CycleIdx: %s\n', step, sum(mask), mat2str(cycle_idx));
        fprintf('  Current range: [%.3f, %.3f] A\n', min(I_data), max(I_data));
        fprintf('  Voltage range: [%.3f, %.3f] V\n', min(V_data), max(V_data));
        fprintf('  Capacity range: [%.3f, %.3f] Ah\n', min(Q_data), max(Q_data));
        
        % Check if it's charge or discharge
        if all(I_data > 0)
            fprintf('  Type: CHARGE\n');
        elseif all(I_data < 0)
            fprintf('  Type: DISCHARGE\n');
        else
            fprintf('  Type: MIXED\n');
        end
    end
end

fprintf('\n=== Check Complete ===\n');
