% Quick check script for Ch15 RPT200cyc C-rate data
clear; clc;

crateMatFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Crate_integrated\Crate_integrated.mat';
load(crateMatFile, 'Crate_data');

fprintf('\n=== Checking Ch15 RPT200cyc C-rate Data ===\n');

channel_key = 'Ch15';
cycle_key = 'cyc200';
crate_labels = {'c01', 'c05', 'c1', 'c2', 'c3'};
crate_names = {'0.1C', '0.5C', '1C', '2C', '3C'};

for r = 1:length(crate_labels)
    label = crate_labels{r};
    rate_name = crate_names{r};
    
    fprintf('\n--- %s (%s) ---\n', rate_name, label);
    
    % Charge
    if isfield(Crate_data, channel_key) && ...
       isfield(Crate_data.(channel_key), cycle_key) && ...
       isfield(Crate_data.(channel_key).(cycle_key), label) && ...
       isfield(Crate_data.(channel_key).(cycle_key).(label), 'charge')
        
        cap_chg = Crate_data.(channel_key).(cycle_key).(label).charge.cap_grid;
        volt_chg = Crate_data.(channel_key).(cycle_key).(label).charge.volt_grid;
        
        fprintf('  Charge: cap_grid length = %d\n', length(cap_chg));
        if ~isempty(cap_chg)
            fprintf('    cap_grid range: [%.3f, %.3f] Ah\n', min(cap_chg), max(cap_chg));
            fprintf('    volt_grid range: [%.3f, %.3f] V\n', min(volt_chg(~isnan(volt_chg))), max(volt_chg(~isnan(volt_chg))));
            fprintf('    NaN count in cap_grid: %d\n', sum(isnan(cap_chg)));
            fprintf('    NaN count in volt_grid: %d\n', sum(isnan(volt_chg)));
        else
            fprintf('    cap_grid is EMPTY!\n');
        end
    else
        fprintf('  Charge: Field not found\n');
    end
    
    % Discharge
    if isfield(Crate_data, channel_key) && ...
       isfield(Crate_data.(channel_key), cycle_key) && ...
       isfield(Crate_data.(channel_key).(cycle_key), label) && ...
       isfield(Crate_data.(channel_key).(cycle_key).(label), 'discharge')
        
        cap_dch = Crate_data.(channel_key).(cycle_key).(label).discharge.cap_grid;
        volt_dch = Crate_data.(channel_key).(cycle_key).(label).discharge.volt_grid;
        
        fprintf('  Discharge: cap_grid length = %d\n', length(cap_dch));
        if ~isempty(cap_dch)
            fprintf('    cap_grid range: [%.3f, %.3f] Ah\n', min(cap_dch), max(cap_dch));
            fprintf('    volt_grid range: [%.3f, %.3f] V\n', min(volt_dch(~isnan(volt_dch))), max(volt_dch(~isnan(volt_dch))));
            fprintf('    NaN count in cap_grid: %d\n', sum(isnan(cap_dch)));
            fprintf('    NaN count in volt_grid: %d\n', sum(isnan(volt_dch)));
        else
            fprintf('    cap_grid is EMPTY!\n');
        end
    else
        fprintf('  Discharge: Field not found\n');
    end
end

fprintf('\n=== Check Complete ===\n');
