
% RPT Feature Discovery Analysis (Correlation Scan)
% Goal: Identify the optimal voltage range for SOH estimation
% Method:
% 1. Load C-rate data
% 2. Define fine voltage grid (e.g., 0.02V intervals)
% 3. Extract capacity (dQ) for each interval across all cycles
% 4. Calculate correlation between interval dQ and Total Capacity (SOH)
% 5. Visualize Correlation vs Voltage to find best features

clear; clc; close all;
warning off;

%% 1. Configuration
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing';
crateFile = fullfile(baseDir, 'Crate_integrated\Crate_integrated.mat');
% Output Directory: Current Folder/Results
currentScriptPath = fileparts(mfilename('fullpath'));
if isempty(currentScriptPath), currentScriptPath = pwd; end % Fallback if run as selection
outputDir = fullfile(currentScriptPath, 'Results');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% Analysis Settings
target_crates = [0.1, 0.5, 1, 2, 3];
crate_keys = {'c01', 'c05', 'c1', 'c2', 'c3'};
channels = {'Ch09','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16'}; % Add channels as needed
use_step_type = 'Charge'; % Discharge is standard for SOH

% Voltage Grid for Scan
v_min = 3.0; 
v_max = 4.2; 
v_step = 0.02; % 20mV resolution
v_grid = v_min:v_step:v_max; 
% Intervals: [v_grid(i), v_grid(i+1)]

fprintf('Loading Data: %s...\n', crateFile);
load(crateFile, 'Crate_data');

%% 2. Feature Scan
fprintf('\n=== Starting Feature Correlation Scan ===\n');

for k = 1:length(crate_keys)
    cKey = crate_keys{k};
    targetCrate = target_crates(k);
    
    fprintf('  Analyzing C-rate: %.1fC (%s)\n', targetCrate, cKey);
    
    % Data Collectors
    cycle_dQ_matrix = []; % Rows: Cycles, Cols: Voltage Intervals
    soh_vector = [];      % Total Capacity per cycle
    cycle_list = [];
    
    % Iterate Channels & Cycles
    for chIdx = 1:length(channels)
        chName = channels{chIdx};
        if ~isfield(Crate_data, chName), continue; end
        
        cycles = fieldnames(Crate_data.(chName));
        
        for cycIdx = 1:length(cycles)
            cycleName = cycles{cycIdx}; % e.g. 'cyc0'
            cycNum = str2double(strrep(cycleName, 'cyc', ''));
            
            if isfield(Crate_data.(chName).(cycleName), cKey)
                stepData = Crate_data.(chName).(cycleName).(cKey);
                
                % Select Data
                if strcmpi(use_step_type, 'Charge')
                    rawData = stepData.charge;
                else
                    rawData = stepData.discharge;
                end
                
                if isempty(rawData.V) || isempty(rawData.Q), continue; end
                
                V = rawData.V;
                Q = rawData.Q;
                
                % unique for interpolation
                [V_unique, unique_idx] = unique(V);
                Q_unique = Q(unique_idx);
                
                % Interpolate Q onto v_grid
                % Use linear interpolation
                % Handle out of range values with NaN
                try
                    Q_interp = interp1(V_unique, Q_unique, v_grid, 'linear');
                catch
                    % If V is not monotonic or too few points
                    continue;
                end
                
                % Calculate dQ for each interval: dQ(i) = |Q(i+1) - Q(i)|
                dQ_interval = abs(diff(Q_interp));
                
                % Check validity (if interval is outside data range, it explains nothing)
                % We can treat NaN as 0 or exclude. 
                % Ideally, if the cycling range was 2.5-4.2V, all should be valid.
                % If range was 3.0-4.2V, then 2.5-3.0 values are NaN.
                % We'll replace NaN with NaN to skip correlation calculation later, or 0?
                % Correlation with NaN propagates NaN.
                
                feature_row = dQ_interval;
                
                % SOH (Total Capacity)
                % We can use max(Q) - min(Q) of this step
                total_cap = max(Q) - min(Q);
                
                cycle_dQ_matrix = [cycle_dQ_matrix; feature_row];
                soh_vector = [soh_vector; total_cap];
                cycle_list = [cycle_list; cycNum];
            end
        end
    end
    
    % Calculate Correlation per Interval
    if size(cycle_dQ_matrix, 1) < 3
        fprintf('    Not enough cycles to correlate.\n');
        continue;
    end
    
    num_intervals = size(cycle_dQ_matrix, 2);
    R_values = zeros(1, num_intervals);
    P_values = zeros(1, num_intervals);
    
    for i = 1:num_intervals
        col_data = cycle_dQ_matrix(:, i);
        
        % Remove NaNs (where voltage range was not reached)
        valid_mask = ~isnan(col_data) & ~isnan(soh_vector);
        
        if sum(valid_mask) > 2
             [r, p] = corr(col_data(valid_mask), soh_vector(valid_mask), 'Type', 'Pearson');
             R_values(i) = r;
             P_values(i) = p;
        else
            R_values(i) = NaN;
        end
    end
    
    % Visualization
    fig = figure('Visible', 'off', 'Position', [100 100 1000 600]);
    
    % Plot 1: Correlation Scan
    subplot(2,1,1);
    % X-axis: Midpoint of voltage interval
    v_mid = (v_grid(1:end-1) + v_grid(2:end)) / 2;
    bar(v_mid, R_values, 'FaceColor', [0.2 0.6 0.8]);
    hold on;
    yline(0.7, 'r--', 'Threshold (0.7)');
    xlabel('Voltage [V]');
    ylabel('Correlation (Pearson r)');
    title(sprintf('Feature Discovery: Correlation vs Voltage (%.1fC)', targetCrate));
    ylim([-1.1 1.1]);
    grid on;
    
    % Highlight Top Zones
    % Find consecutive regions with R > 0.95
    high_corr_mask = R_values > 0.95;
    % plotting logic could be complex, simple bar is good enough for now.
    
    % Plot 2: Variance/Validity Check
    subplot(2,1,2);
    % Plot mean dQ/dV profile (just to show where the peaks are)
    % Approx dQ/dV ~ dQ_interval / v_step
    mean_dQ = mean(cycle_dQ_matrix, 1, 'omitnan');
    mean_dQdV = mean_dQ / v_step;
    plot(v_mid, mean_dQdV, 'k-', 'LineWidth', 1.5);
    ylabel('Mean dQ/dV (Ah/V)');
    xlabel('Voltage [V]');
    title('Mean dQ/dV Profile (for reference)');
    grid on;
    
    saveName = fullfile(outputDir, sprintf('Discovery_Scan_%.1fC.fig', targetCrate));
    savefig(fig, saveName);
    close(fig);
    
    fprintf('    Saved Scan: %s\n', saveName);
    
    % Report High Correlation Ranges
    fprintf('    High Correlation Segments (|r| > 0.70):\n');
    indices = find(abs(R_values) > 0.70);
    if ~isempty(indices)
        % Group indices into continuous ranges
        diff_idx = diff(indices);
        breaks = [0, find(diff_idx > 1), length(indices)];
        
        for k = 1:length(breaks)-1
            range_start = indices(breaks(k)+1);
            range_end = indices(breaks(k+1));
            
            mean_R = mean(R_values(range_start:range_end));
            v_s = v_grid(range_start);
            v_e = v_grid(range_end+1);
            
            fprintf('      Range %.2fV - %.2fV (Mean r = %.4f)\n', v_s, v_e, mean_R);
        end
    else
        fprintf('      None found.\n');
    end
end

fprintf('\n=== Feature Discovery Complete ===\n');
fprintf('Check %s for correlation profiles.\n', outputDir);
