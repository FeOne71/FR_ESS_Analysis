function ch_struct = process_one_channel_crate(ch_idx, channels, rpt_cycles, ExperimentalDataPath, crate_labels, crate_steps_charge, crate_steps_discharge)
% Process C-rate data for one channel (for parfor/serial).
channel = channels{ch_idx};
ch_struct = struct();
for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    cycle_key = sprintf('cyc%s', rpt_cycle(1:end-3));
    if strcmp(channel, 'Ch14') && strcmp(rpt_cycle, '800cyc')
        continue;
    end
    filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
    filepath = fullfile(ExperimentalDataPath, filename);
    if ~isfile(filepath)
        continue;
    end
    try
        T = readtable(filepath, 'VariableNamingRule', 'preserve');
    catch
        continue;
    end
    ch_struct.(cycle_key) = struct();
    t1Col = find_t1_col(T);
    for r = 1:length(crate_labels)
        label = crate_labels{r};
        step_chg = crate_steps_charge(r);
        step_dch = crate_steps_discharge(r);
        mask_chg = (T.("Step Index") == step_chg);
        V_chg = T.("Voltage(V)")(mask_chg);
        I_chg = T.("Current(A)")(mask_chg);
        Q_chg = T.("Capacity(Ah)")(mask_chg);
        if any(strcmp('Time', T.Properties.VariableNames))
            t_chg = T.("Time")(mask_chg);
        else
            t_chg = [];
        end
        T1_chg = get_col_double(T, mask_chg, t1Col);
        mask_dch = (T.("Step Index") == step_dch);
        V_dch = T.("Voltage(V)")(mask_dch);
        I_dch = T.("Current(A)")(mask_dch);
        Q_dch = T.("Capacity(Ah)")(mask_dch);
        if any(strcmp('Time', T.Properties.VariableNames))
            t_dch = T.("Time")(mask_dch);
        else
            t_dch = [];
        end
        T1_dch = get_col_double(T, mask_dch, t1Col);
        charge_struct = struct('V', V_chg, 'I', I_chg, 'Q', Q_chg, 't', t_chg, 'T1_raw', T1_chg);
        discharge_struct = struct('V', V_dch, 'I', I_dch, 'Q', Q_dch, 't', t_dch, 'T1_raw', T1_dch);
        ch_struct.(cycle_key).(label).charge = charge_struct;
        ch_struct.(cycle_key).(label).discharge = discharge_struct;
    end
end
end

function j = find_t1_col(T)
% 온도 변수: T1 (°C)
j = [];
for i = 1:numel(T.Properties.VariableNames)
    if contains(char(T.Properties.VariableNames{i}), 'T1')
        j = i; return;
    end
end
end

function x = get_col_double(T, mask, col)
if isempty(col), x = []; return; end
try
    x = double(T{mask, col}(:));
catch
    x = [];
end
end
