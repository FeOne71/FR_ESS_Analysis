% Analyze_HVAC_Field.m
% Phase 2: Analyzes Field Data MAT files to determine the HVAC sensor
% characteristics and check for the presence of PLC ambient temperature data.

clear; clc; close all;

sample_file = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old\2021\202106\Raw_20210601.mat';
fprintf('Loading sample field data: %s\n', sample_file);

if ~exist(sample_file, 'file')
    error('File not found!');
end

data = load(sample_file);

disp('--- Variables at root ---');
disp(fieldnames(data));

% Assuming the data is stored in 'Rack_data' or 'merged_data' or similar
keys = fieldnames(data);
main_var = data.(keys{1});
disp('--- Variables inside main struct ---');
var_fields = fieldnames(main_var);
disp(var_fields);

% Look for HVAC or PLC keywords
found_hvac = false;
hvac_vars = {};
for i = 1:length(var_fields)
    if contains(lower(var_fields{i}), 'hvac') || contains(lower(var_fields{i}), 'plc') || contains(lower(var_fields{i}), 'temp')
        found_hvac = true;
        hvac_vars{end+1} = var_fields{i};
    end
end

if found_hvac
    fprintf('Found potential HVAC/Temp variables: \n');
    disp(hvac_vars');
    
    % Let's plot the HVAC RoomTemp if multiple exist
    fig = figure('Name', 'HVAC Comparison', 'Position', [100 100 800 400]);
    hold on;
    legend_entries = {};
    for i = 1:length(hvac_vars)
        if contains(lower(hvac_vars{i}), 'roomtemp') || contains(lower(hvac_vars{i}), 'hvac')
            % Ensure it's a numeric array
            val = main_var.(hvac_vars{i});
            if isnumeric(val) && length(val) > 1
                plot(val, 'LineWidth', 1.5);
                legend_entries{end+1} = hvac_vars{i};
            end
        end
    end
    if ~isempty(legend_entries)
        legend(legend_entries, 'Interpreter', 'none');
        title('HVAC Room Temperature Comparison');
        xlabel('Time index'); ylabel('Temperature [C/F]');
        grid on;
        
        save_path = fullfile(pwd, 'HVAC_Field_Analysis.fig');
        savefig(fig, save_path);
        fprintf('Saved HVAC plot to %s\n', save_path);
    else
        fprintf('No valid arrays to plot.\n');
    end
else
    fprintf('WARNING: No HVAC/PLC data found in this MAT file! This means the New rack data script might not have appended it.\n');
end
