%% get_steps.m
opts = detectImportOptions('D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data\RPT\Ch09_RPT_0cyc.csv');
opts.SelectedVariableNames = {'StepIndex', 'Current_A_', 'Voltage_V_', 'T1___'};
data = readtable('D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data\RPT\Ch09_RPT_0cyc.csv', opts);

steps = data.StepIndex;
I = data.Current_A_;
V = data.Voltage_V_;
T = data.T1___;

unique_steps = unique(steps);
for i = 1:length(unique_steps)
    idx = (steps == unique_steps(i));
    I_step = I(idx);
    V_step = V(idx);
    duration_sec = sum(idx); % roughly 1 sec per row
    % DCIR pulses are usually 10s-30s with high current
    if abs(mean(I_step)) > 1 && duration_sec < 600
        fprintf('Potential DCIR Step %2d: Duration=%5d, Mean I=%.2f, Max I=%.2f, Min I=%.2f\n', ...
            unique_steps(i), duration_sec, mean(I_step), max(I_step), min(I_step));
    end
end
