% One-off: Load CSV T1(℃), compute 8-cell average T Range* for Lab RPT V,I,T Table.
% Charge = OCV step 8 (sub 2), Discharge = Static step 3 (sub 2).
ExperimentalDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data\RPT';
channels = {'Ch09','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16'};
rpt_cycles = {'0cyc','200cyc','400cyc','600cyc','800cyc','1000cyc'};
cycle_nums = [0, 200, 400, 600, 800, 1000];
ocv_charge_step = 8;
static_step = 3;
sub_idx = 2;  % Cycle Index == 2
stepCol = 2; subCol = 4;

T_Charge = cell(length(rpt_cycles), 7);
T_Charge(:,1) = num2cell(cycle_nums);
T_Charge(:,2) = repmat({'Charge'}, length(rpt_cycles), 1);
T_Charge(:,3:5) = {''};
T_Charge(:,7) = {''};

T_Discharge = cell(length(rpt_cycles), 7);
T_Discharge(:,1) = num2cell(cycle_nums);
T_Discharge(:,2) = repmat({'Discharge'}, length(rpt_cycles), 1);
T_Discharge(:,3:5) = {''};
T_Discharge(:,7) = {''};

for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    % Charge: OCV step 8
    t_min_ch = []; t_max_ch = [];
    for ch_idx = 1:length(channels)
        ch = channels{ch_idx};
        filepath = fullfile(ExperimentalDataPath, sprintf('%s_RPT_%s.csv', ch, rpt_cycle));
        if ~isfile(filepath), continue; end
        try
            T = readtable(filepath, 'VariableNamingRule', 'preserve');
        catch, continue; end
        t1Col = [];
        for j = 1:numel(T.Properties.VariableNames)
            if contains(lower(char(T.Properties.VariableNames{j})), 't1')
                t1Col = j; break;
            end
        end
        if isempty(t1Col), continue; end
        mask = (T{:,stepCol} == ocv_charge_step) & (T{:,subCol} == sub_idx);
        if sum(mask) < 2, continue; end
        t1_vals = double(T{mask, t1Col}(:));
        t_min_ch = [t_min_ch, min(t1_vals)];
        t_max_ch = [t_max_ch, max(t1_vals)];
    end
    if ~isempty(t_min_ch)
        T_Charge{rpt_idx, 6} = sprintf('%.2f ~ %.2f', mean(t_min_ch,'omitnan'), mean(t_max_ch,'omitnan'));
    else
        T_Charge{rpt_idx, 6} = '';
    end

    % Discharge: Static step 3
    t_min_ch = []; t_max_ch = [];
    for ch_idx = 1:length(channels)
        ch = channels{ch_idx};
        filepath = fullfile(ExperimentalDataPath, sprintf('%s_RPT_%s.csv', ch, rpt_cycle));
        if ~isfile(filepath), continue; end
        try
            T = readtable(filepath, 'VariableNamingRule', 'preserve');
        catch, continue; end
        t1Col = [];
        for j = 1:numel(T.Properties.VariableNames)
            if contains(lower(char(T.Properties.VariableNames{j})), 't1')
                t1Col = j; break;
            end
        end
        if isempty(t1Col), continue; end
        mask = (T{:,stepCol} == static_step) & (T{:,subCol} == sub_idx);
        if sum(mask) < 2, continue; end
        t1_vals = double(T{mask, t1Col}(:));
        t_min_ch = [t_min_ch, min(t1_vals)];
        t_max_ch = [t_max_ch, max(t1_vals)];
    end
    if ~isempty(t_min_ch)
        T_Discharge{rpt_idx, 6} = sprintf('%.2f ~ %.2f', mean(t_min_ch,'omitnan'), mean(t_max_ch,'omitnan'));
    else
        T_Discharge{rpt_idx, 6} = '';
    end
end

vnames = {'Cycle','Type','V_Range','I_Range','P_Range','T_Range','t_Range'};
T_charge = cell2table(T_Charge, 'VariableNames', vnames);
T_discharge = cell2table(T_Discharge, 'VariableNames', vnames);
T_charge.Properties.VariableNames = {'Cycle','Type','V Range','I Range','P Range*','T Range*','t Range*'};
T_discharge.Properties.VariableNames = {'Cycle','Type','V Range','I Range','P Range*','T Range*','t Range*'};

fprintf('=== Lab RPT V, I, T Table (*: average of 8 cells) - T Range* from T1 ===\n\n');
fprintf('--- Charge (OCV step 8, T1 8-cell avg min~max [℃]) ---\n');
disp(T_charge);
fprintf('--- Discharge (Static step 3, T1 8-cell avg min~max [℃]) ---\n');
disp(T_discharge);
