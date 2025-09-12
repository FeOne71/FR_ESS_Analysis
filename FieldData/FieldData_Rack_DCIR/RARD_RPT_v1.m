%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RARD_RPT_v1.m
%
% DCIR 계산을 위한 스크립트
% - 파일: Raw_20210607.mat
% - Idle threshold: 64*0.01 = 0.64 A
% - DCIR 공식: R = (V2 - V1) / I
% - HPPC 구간: 14:15:00 ~ 14:50:00
% - 14개 셀에 대한 DCIR 히스토그램 생성
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% Parameters
idle_thresh = 64 * 0.03;  % 0.64 A
hppc_start = '14:15:00';
hppc_end = '14:50:00';

% Select racks and modules to analyze
wantedRack = {'Rack01', 'Rack02', 'Rack03', 'Rack04', 'Rack05', 'Rack06', 'Rack07', 'Rack08'};
wantedModule = {'Module01', 'Module02', 'Module03', 'Module04', 'Module05', 'Module06', 'Module07', 'Module08', 'Module09', 'Module10', 'Module11', 'Module12', 'Module13', 'Module14', 'Module15', 'Module16', 'Module17'};

%% Load data
dataFile = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\RARDsync\2021\202106\Raw_20210607.mat';
S = load(dataFile);

% Debug: Check data structure
fprintf('=== Data Structure Debug ===\n');
fprintf('Loaded variables: %s\n', strjoin(fieldnames(S), ', '));

% Check for Raw field
if isfield(S, 'Raw')
    Raw = S.Raw;
    fprintf('Raw field found\n');
else
    error('Raw not found in %s', dataFile);
end

% Debug: Check Raw structure
fprintf('Raw fields: %s\n', strjoin(fieldnames(Raw), ', '));

% Check for rack fields (Rack01, Rack02, etc.)
rack_fields = fieldnames(Raw);
rack_fields = rack_fields(startsWith(rack_fields, 'Rack'));
% fprintf('Available racks: %s\n', strjoin(rack_fields, ', '));

% Filter wanted racks
available_racks = rack_fields;
selected_racks = intersect(wantedRack, available_racks);

if isempty(selected_racks)
    error('No wanted racks found in data. Available: %s', strjoin(available_racks, ', '));
end

% fprintf('Available racks: %s\n', strjoin(available_racks, ', '));
fprintf('Selected racks: %s\n', strjoin(selected_racks, ', '));

% Get time vector and current from first rack, first module (assuming all have same structure)
first_rack = rack_fields{1};
D = Raw.(first_rack);
module_fields = fieldnames(D);
first_module = module_fields{1};
Module = D.(first_module);

% fprintf('First module (%s) fields: %s\n', first_module, strjoin(fieldnames(Module), ', '));
% fprintf('========================\n\n');

% Time vector
if isfield(Module, 'Time')
    t = datetime(Module.Time);
elseif isfield(Module, 'Date_Time')
    if isduration(Module.Date_Time)
        t = datetime(2021,6,7) + Module.Date_Time;
    else
        t = datetime(Module.Date_Time);
    end
else
    error('No Time/Date_Time field present in module.');
end

% Rack-level signals (from RBMS fields)
if isfield(Module, 'RBMS_DCCurrent_A')
    I_rack = Module.RBMS_DCCurrent_A(:);    % A (rack)
else
    error('No RBMS_DCCurrent_A field found in module.');
end

fprintf('Time range: %s to %s\n', datestr(t(1)), datestr(t(end)));
fprintf('Data length: %d points\n', length(t));

%% Filter HPPC time range
hppc_start_time = datetime(2021,6,7) + duration(hppc_start);
hppc_end_time = datetime(2021,6,7) + duration(hppc_end);

% Find indices within HPPC range
hppc_mask = (t >= hppc_start_time) & (t <= hppc_end_time);
hppc_idx = find(hppc_mask);

if isempty(hppc_idx)
    error('No data found in HPPC time range: %s to %s', hppc_start, hppc_end);
end

% Extract HPPC data
t_hppc = t(hppc_idx);
I_hppc = I_rack(hppc_idx);

fprintf('HPPC data loaded: %d points from %s to %s\n', ...
    length(hppc_idx), datestr(t_hppc(1)), datestr(t_hppc(end)));

%% Detect idle/charge/discharge segments
isIdle = abs(I_hppc) < idle_thresh;
isCharge = I_hppc > idle_thresh;
isDischarge = I_hppc < -idle_thresh;

% Helper function to find contiguous segments
find_segments = @(mask) local_find_segments(mask);

idleSegs = find_segments(isIdle);
chargeSegs = find_segments(isCharge);
dischargeSegs = find_segments(isDischarge);

fprintf('Found %d idle segments, %d charge segments, %d discharge segments\n', ...
    size(idleSegs,1), size(chargeSegs,1), size(dischargeSegs,1));

%% DCIR Calculation for all racks, modules and cells (First pulse only)
DCIR_all_cells = [];  % Store all DCIR values for histogram
DCIR_by_rack = {};    % Store DCIR values by rack
DCIR_by_module = {};  % Store DCIR values by module

% Process selected racks
for rack_idx = 1:length(selected_racks)
    rack_name = selected_racks{rack_idx};
    D = Raw.(rack_name);
    available_modules = fieldnames(D);
    
    % Filter wanted modules for this rack
    selected_modules = intersect(wantedModule, available_modules);
    
    if isempty(selected_modules)
        fprintf('Warning: No wanted modules found in rack %s. Available: %s\n', rack_name, strjoin(available_modules, ', '));
        continue;
    end
    
    % fprintf('Processing rack %s: %d selected modules (%s)\n', rack_name, length(selected_modules), strjoin(selected_modules, ', '));
    DCIR_this_rack = [];
    
    % Process selected modules in this rack
    for module_idx = 1:length(selected_modules)
        module_name = selected_modules{module_idx};
        Module = D.(module_name);
        
        % Get cell voltage data for this module
        cell_voltages = [];
        % Extract module number from module name (e.g., Module01 -> 1)
        module_num = str2double(module_name(end-1:end));
        for i = 1:14
            cell_field = sprintf('M%d_Cell%d', module_num, i);  % M1_Cell1, M2_Cell1, etc.
            if isfield(Module, cell_field)
                cell_voltages(:,i) = Module.(cell_field)(:);
            else
                error('Cell field %s not found in module %s', cell_field, module_name);
            end
        end
        
        % Get module current data
        I_module = Module.RBMS_DCCurrent_A(:);
        
        % Extract HPPC data for this module
        V_cells_hppc = cell_voltages(hppc_idx, :);
        I_hppc = I_module(hppc_idx);
        
        fprintf('  Processing module %s: %d cells\n', module_name, size(cell_voltages,2));
        DCIR_this_module = [];
        
        % Process first discharge segment only
        if size(dischargeSegs,1) >= 1
            k = 1;  % First discharge segment only
            dch_start = dischargeSegs(k,1);
            dch_end = dischargeSegs(k,2);
            
            % Find idle segment before discharge
            prev_idle = [];
            for i = 1:size(idleSegs,1)
                if idleSegs(i,2) < dch_start
                    prev_idle = idleSegs(i,:);
                end
            end
            
            if ~isempty(prev_idle)
                % V2: Find the exact point where current exceeds idle threshold
                % Look for the transition from idle to discharge
                v2_idx = prev_idle(2) + 1;  % Start from the point right after idle ends
                
                % Find the first point where |I| > idle_thresh
                while v2_idx <= length(I_hppc) && abs(I_hppc(v2_idx)) <= idle_thresh
                    v2_idx = v2_idx + 1;
                end
                
                % Make sure we found a valid transition point
                if v2_idx <= length(I_hppc)
                    % I1: current at idle (last point of idle segment)
                    I1 = I_hppc(prev_idle(2));
                    % I2: current at idle threshold crossing point (same as V2)
                    I2 = I_hppc(v2_idx);
                    
                    % Debug: Print summary for this module
                    fprintf('    Debug: I1=%.1f, I2=%.1f, ΔI=%.1f\n', I1, I2, I2-I1);
                    
                    if I2 ~= 0
                        fprintf('    Debug: Calculating DCIR for %d cells\n', 14);
                    else
                        fprintf('    Debug: I2=0, skipping DCIR calculation\n');
                    end
                    
                    if I2 ~= 0
                        % Calculate DCIR for each cell in this module
                        for cell_idx = 1:14
                            % V1: idle voltage (last point of idle segment)
                            V1 = V_cells_hppc(prev_idle(2), cell_idx);
                            % V2: voltage at idle threshold crossing point
                            V2 = V_cells_hppc(v2_idx, cell_idx);
                            
                            % DCIR calculation: R = (V2 - V1) / (I2 - I1)
                            R = (V2 - V1) / (I2 - I1);
                            DCIR_value = R * 1000;  % Convert to mOhm
                            
                            % Debug: Print summary for first cell only
                            if cell_idx == 1
                                fprintf('    Debug: V1=%.3f, V2=%.3f, ΔV=%.3f, DCIR=%.1f mOhm\n', V1, V2, V2-V1, DCIR_value);
                            end
                            
                            DCIR_all_cells(end+1) = DCIR_value;
                            DCIR_this_rack(end+1) = DCIR_value;
                            DCIR_this_module(end+1) = DCIR_value;
                        end
                    end
                end
            end
        end
        
        % Store module data
        DCIR_by_module{end+1} = struct('rack', rack_name, 'module', module_name, 'data', DCIR_this_module);
    end
    
    % Store rack data
    DCIR_by_rack{end+1} = struct('rack', rack_name, 'data', DCIR_this_rack);
    fprintf('Rack %s: %d DCIR values\n', rack_name, length(DCIR_this_rack));
end

%% Display results
if ~isempty(DCIR_all_cells)
    fprintf('\n=== DCIR Results ===\n');
    fprintf('Total DCIR measurements: %d\n', length(DCIR_all_cells));
    fprintf('DCIR Statistics (mOhm):\n');
    fprintf('Mean: %.3f mOhm\n', mean(DCIR_all_cells));
    fprintf('Std: %.3f mOhm\n', std(DCIR_all_cells));
    fprintf('Min: %.3f mOhm\n', min(DCIR_all_cells));
    fprintf('Max: %.3f mOhm\n', max(DCIR_all_cells));
    fprintf('Median: %.3f mOhm\n', median(DCIR_all_cells));
else
    fprintf('No DCIR results found.\n');
end

%% Create DCIR Histogram and Current Plot with Markers
if ~isempty(DCIR_all_cells)
    % 1. DCIR Histogram
    fig1 = figure('Name','DCIR Histogram - All Cells','NumberTitle','off');
    
    % Create histogram
    histogram(DCIR_all_cells, 30, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'black', 'FaceAlpha', 0.7);
    
    % Add statistics lines
    hold on;
    xline(mean(DCIR_all_cells), 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Mean: %.2f mOhm', mean(DCIR_all_cells)));
    xline(median(DCIR_all_cells), 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Median: %.2f mOhm', median(DCIR_all_cells)));
    
    % Formatting
    xlabel('DCIR [mOhm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Count', 'FontSize', 12, 'FontWeight', 'bold');
    title('DCIR Distribution - All Cells', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    legend('Location', 'best');
    
    % Add text box with statistics
    stats_text = sprintf('N = %d\nMean = %.2f mOhm\nStd = %.2f mOhm\nMin = %.2f mOhm\nMax = %.2f mOhm', ...
        length(DCIR_all_cells), mean(DCIR_all_cells), std(DCIR_all_cells), ...
        min(DCIR_all_cells), max(DCIR_all_cells));
    text(0.02, 0.98, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
    
    % Set axis limits for better visualization
    xlim([min(DCIR_all_cells)*0.9, max(DCIR_all_cells)*1.1]);
    
    fprintf('\nDCIR Histogram created successfully!\n');
else
    fprintf('No data available for histogram.\n');
end

% 2. Current Plot with V1, V2, I1, I2 Markers
fig2 = figure('Name','Current Plot with DCIR Calculation Points','NumberTitle','off');

% Get current data from first selected rack and module
if ~isempty(selected_racks)
    rack_name = selected_racks{1};
    D = Raw.(rack_name);
    available_modules = fieldnames(D);
    selected_modules = intersect(wantedModule, available_modules);
    
    if ~isempty(selected_modules)
        module_name = selected_modules{1};
        Module = D.(module_name);
        
        % Get module current data
        I_module = Module.RBMS_DCCurrent_A(:);
        I_hppc = I_module(hppc_idx);
        
        % Plot current
        plot(t_hppc, I_hppc, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 1);
        hold on;
        
        % Add idle threshold lines
        yline(idle_thresh, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1);
        yline(-idle_thresh, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1);
        
        % Find and mark V1, V2, I1, I2 points
        if size(dischargeSegs,1) >= 1
            k = 1;
            dch_start = dischargeSegs(k,1);
            
            % Find idle segment before discharge
            prev_idle = [];
            for i = 1:size(idleSegs,1)
                if idleSegs(i,2) < dch_start
                    prev_idle = idleSegs(i,:);
                end
            end
            
            if ~isempty(prev_idle)
                % V2: Find the exact point where current exceeds idle threshold
                v2_idx = prev_idle(2) + 1;
                while v2_idx <= length(I_hppc) && abs(I_hppc(v2_idx)) <= idle_thresh
                    v2_idx = v2_idx + 1;
                end
                
                if v2_idx <= length(I_hppc) && (v2_idx + 1) <= length(I_hppc)
                    % Mark V1, I1 point (idle end)
                    plot(t_hppc(prev_idle(2)), I_hppc(prev_idle(2)), 'bo', 'MarkerSize', 10, 'MarkerFaceColor','b', 'LineWidth', 2);
                    text(t_hppc(prev_idle(2)), I_hppc(prev_idle(2)), ' V1,I1', 'FontSize', 10, 'Color', 'blue');
                    
                    % Mark V2, I2 point (idle threshold crossing)
                    plot(t_hppc(v2_idx), I_hppc(v2_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r', 'LineWidth', 2);
                    text(t_hppc(v2_idx), I_hppc(v2_idx), ' V2', 'FontSize', 10, 'Color', 'red');
                    
                    % Mark I2 point (same as V2)
                    plot(t_hppc(v2_idx), I_hppc(v2_idx), 'go', 'MarkerSize', 10, 'MarkerFaceColor','g', 'LineWidth', 2);
                    text(t_hppc(v2_idx), I_hppc(v2_idx), ' I2', 'FontSize', 10, 'Color', 'green');
                end
            end
        end
        
        % Formatting
        xlabel('Time', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Current [A]', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Current vs Time - %s %s (HPPC Range)', rack_name, module_name), 'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        legend('Current', 'Idle Threshold', 'V1,I1', 'V2', 'I2', 'Location', 'best');
        
        fprintf('Current plot with DCIR calculation points created successfully!\n');
    end
end


%% Helper function
function segments = local_find_segments(mask)
    % Find contiguous true segments in logical array
    % Returns Nx2 array [start_idx, end_idx]
    
    if isempty(mask)
        segments = [];
        return;
    end
    
    % Find transitions
    diff_mask = diff([false; mask(:); false]);
    starts = find(diff_mask == 1);
    ends = find(diff_mask == -1) - 1;
    
    if isempty(starts)
        segments = [];
        return;
    end
    
    segments = [starts, ends];
end