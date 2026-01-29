
% Calculate_Diff_PValues.m
% Calculate p-values for R-R1 differences for specific parameter set

clear; clc;

% --- Configuration ---
baseDir = 'd:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\ver03_OnlyChargeEvents\Results\min10_std2_rng0_Both';
cyclesWanted = {'0cyc','200cyc','400cyc','600cyc','800cyc','1000cyc'}; % Adjust as needed
socLevels = {'SOC50', 'SOC70', 'SOC90'};
timePoints = [5, 10, 30, 60]; % Compare these with 1s

% --- Helper Functions ---
function cycleField = make_valid_cycle_field(cycleName)
    if ~isempty(regexp(cycleName, '^\d', 'once'))
        cycleField = ['cyc_' cycleName];
    else
        cycleField = cycleName;
    end
end

function pVal = compute_pvalue(vals, groups)
    pVal = NaN;
    if isempty(vals) || isempty(groups) || numel(unique(groups)) < 2
        return;
    end
    pVal = kruskalwallis(vals, groups, 'off');
end

% --- Main Logic ---
fprintf('=== P-Value Calculation for Resistance Differences ===\n');

% Initialize data structures
aggCharge = struct();
aggDischarge = struct();

% Load all cycle files
files = dir(fullfile(baseDir, 'Lab_DC_Events_Features_*cyc.mat'));

for i = 1:length(files)
    filePath = fullfile(baseDir, files(i).name);
    token = regexp(files(i).name, 'Events_Features_(\d+cyc)', 'tokens', 'once');
    if isempty(token), continue; end
    cycleName = token{1};
    if ~ismember(cycleName, cyclesWanted), continue; end
    
    cycleField = make_valid_cycle_field(cycleName);
    dataStruct = load(filePath);
    
    % Find main data variable
    vars = fieldnames(dataStruct);
    dataVarName = vars{1}; % Assuming first variable
    if isfield(dataStruct, 'Lab_DC_DCIR_Feature_Table'), dataVarName = 'Lab_DC_DCIR_Feature_Table'; end 
    % Note: Based on previous file listing, variable name might vary, but usually it's the main struct
    % Let's assume the structure from the view_file output of script 03
    % Actually, script 03 logic had finding logic:
    for v = 1:length(vars)
        if contains(vars{v}, 'Lab_DC_DCIR_')
             dataVarName = vars{v};
             break;
        end
    end
    data = dataStruct.(dataVarName);
    
    % Collect Charge Data
    chFields = fieldnames(data);
    chargeChs = chFields(contains(chFields, '_Charge'));
    for cIdx = 1:length(chargeChs)
        chName = chargeChs{cIdx};
        if ~isfield(data.(chName), 'SOC50'), continue; end % Quick check
        
        for sIdx = 1:length(socLevels)
            socName = socLevels{sIdx};
            if ~isfield(data.(chName), socName), continue; end
            
            profs = fieldnames(data.(chName).(socName));
            for pIdx = 1:length(profs)
                profName = profs{pIdx};
                evts = fieldnames(data.(chName).(socName).(profName));
                
                for eIdx = 1:length(evts)
                    evt = data.(chName).(socName).(profName).(evts{eIdx});
                    
                    % Calculate R1
                    if isfield(evt, 'Rchg_1s')
                        r1 = evt.Rchg_1s;
                        if isfinite(r1)
                            for t = timePoints
                                rField = sprintf('Rchg_%ds', t);
                                if isfield(evt, rField)
                                    val = evt.(rField) - r1;
                                    if isfinite(val)
                                         fieldKey = ['T_' num2str(t)];
                                         if ~isfield(aggCharge, socName), aggCharge.(socName) = struct(); end
                                         if ~isfield(aggCharge.(socName), fieldKey), aggCharge.(socName).(fieldKey) = struct('vals', [], 'groups', []); end
                                         aggCharge.(socName).(fieldKey).vals = [aggCharge.(socName).(fieldKey).vals; val];
                                         aggCharge.(socName).(fieldKey).groups = [aggCharge.(socName).(fieldKey).groups; {cycleName}];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % Collect Discharge Data (Similar logic)
    dischargeChs = chFields(contains(chFields, '_Discharge'));
    for cIdx = 1:length(dischargeChs)
        chName = dischargeChs{cIdx};
        
        for sIdx = 1:length(socLevels)
            socName = socLevels{sIdx};
             if ~isfield(data.(chName), socName), continue; end
             
            profs = fieldnames(data.(chName).(socName));
            for pIdx = 1:length(profs)
                profName = profs{pIdx};
                evts = fieldnames(data.(chName).(socName).(profName));
                
                for eIdx = 1:length(evts)
                    evt = data.(chName).(socName).(profName).(evts{eIdx});
                    
                     % Calculate R1
                    if isfield(evt, 'Rdchg_1s')
                        r1 = evt.Rdchg_1s;
                        if isfinite(r1)
                            for t = timePoints
                                rField = sprintf('Rdchg_%ds', t);
                                if isfield(evt, rField)
                                    val = evt.(rField) - r1;
                                    if isfinite(val)
                                         fieldKey = ['T_' num2str(t)];
                                         if ~isfield(aggDischarge, socName), aggDischarge.(socName) = struct(); end
                                         if ~isfield(aggDischarge.(socName), fieldKey), aggDischarge.(socName).(fieldKey) = struct('vals', [], 'groups', []); end
                                         aggDischarge.(socName).(fieldKey).vals = [aggDischarge.(socName).(fieldKey).vals; val];
                                         aggDischarge.(socName).(fieldKey).groups = [aggDischarge.(socName).(fieldKey).groups; {cycleName}];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% Function to format p-value
fmtP = @(p) sprintf('%.3e', p);

fprintf('\n--- Charge P-Values ---\n');
fprintf('SOC\t\tR5-R1\t\tR10-R1\t\tR30-R1\t\tR60-R1\n');
for sIdx = 1:length(socLevels)
    socName = socLevels{sIdx};
    fprintf('%s\t', socName);
    for t = timePoints
        fieldKey = ['T_' num2str(t)];
        if isfield(aggCharge, socName) && isfield(aggCharge.(socName), fieldKey)
            p = compute_pvalue(aggCharge.(socName).(fieldKey).vals, aggCharge.(socName).(fieldKey).groups);
            fprintf('%s\t', fmtP(p));
        else
            fprintf('NaN\t\t');
        end
    end
    fprintf('\n');
end

fprintf('\n--- Discharge P-Values ---\n');
fprintf('SOC\t\tR5-R1\t\tR10-R1\t\tR30-R1\t\tR60-R1\n');
for sIdx = 1:length(socLevels)
    socName = socLevels{sIdx};
    fprintf('%s\t', socName);
    for t = timePoints
        fieldKey = ['T_' num2str(t)];
        if isfield(aggDischarge, socName) && isfield(aggDischarge.(socName), fieldKey)
             p = compute_pvalue(aggDischarge.(socName).(fieldKey).vals, aggDischarge.(socName).(fieldKey).groups);
            fprintf('%s\t', fmtP(p));
        else
            fprintf('NaN\t\t');
        end
    end
    fprintf('\n');
end
