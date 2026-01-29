% Check CycleIdx 30 for C-rate steps in Ch15 RPT200cyc
clear; clc;

csvFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data\RPT\Ch15_RPT_200cyc.csv';
fprintf('Checking CycleIdx 30 in: %s\n', csvFile);

T = readmatrix(csvFile);

% Check CycleIdx 30
mask_cycle30 = (T(:,4) == 30);
fprintf('\nCycleIdx 30 data points: %d\n', sum(mask_cycle30));

if sum(mask_cycle30) > 0
    unique_steps_cycle30 = unique(T(mask_cycle30, 2));
    fprintf('Step indices in CycleIdx 30: %s\n', mat2str(sort(unique_steps_cycle30)));
    
    % Check C-rate steps
    crate_steps_charge = [28 32 36 40 44];
    crate_steps_discharge = [48 52 56 60 64];
    
    fprintf('\n=== C-rate Steps in CycleIdx 30 ===\n');
    fprintf('Charge steps:\n');
    for i = 1:length(crate_steps_charge)
        step = crate_steps_charge(i);
        mask = (T(:,2) == step) & (T(:,4) == 30);
        count = sum(mask);
        if count > 0
            fprintf('  Step %d: %d points\n', step, count);
        else
            fprintf('  Step %d: NO DATA\n', step);
        end
    end
    
    fprintf('\nDischarge steps:\n');
    for i = 1:length(crate_steps_discharge)
        step = crate_steps_discharge(i);
        mask = (T(:,2) == step) & (T(:,4) == 30);
        count = sum(mask);
        if count > 0
            fprintf('  Step %d: %d points\n', step, count);
        else
            fprintf('  Step %d: NO DATA\n', step);
        end
    end
else
    fprintf('CycleIdx 30 has NO DATA\n');
end

% Check all CycleIdx values
fprintf('\n=== All CycleIdx values ===\n');
all_cycle_idx = unique(T(:,4));
fprintf('CycleIdx values: %s\n', mat2str(sort(all_cycle_idx)));

% For each CycleIdx, check if C-rate steps exist
fprintf('\n=== C-rate Steps by CycleIdx ===\n');
crate_steps_all = [28 32 36 40 44 48 52 56 60 64];
for cyc = sort(all_cycle_idx)'
    mask_cyc = (T(:,4) == cyc);
    steps_in_cyc = unique(T(mask_cyc, 2));
    crate_steps_found = intersect(crate_steps_all, steps_in_cyc);
    if ~isempty(crate_steps_found)
        fprintf('CycleIdx %d: Steps %s\n', cyc, mat2str(crate_steps_found));
    end
end

fprintf('\n=== Check Complete ===\n');
