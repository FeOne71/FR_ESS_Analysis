% RPT_FeatureExtractor_ver0317.m
% Generates feature matrix from Lab RPT Data based on Reference Paper 
% (Fragmented Charge Capacity mapping to Static Capacity).

clear; clc; close all; warning off;

%% 1. Setup Directories
baseDir = 'C:\Users\Chulwon Jung\내 드라이브(jcw9406@gmail.com)\KENTECH\Projects\KEPCO_ESS_Local\ExperimentalData';
vqGridFile = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
saveDir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver0317');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

%% 2. Load Master Rulers from Pseudo-OCV
mrFile = fullfile(saveDir, 'MasterRulers_PseudoOCV.mat');
if exist(mrFile, 'file')
    S_mr = load(mrFile, 'MasterRuler_ver0317');
    MR_Chg = S_mr.MasterRuler_ver0317.V_bounds_chg;
    MR_Dch = S_mr.MasterRuler_ver0317.V_bounds_dch;
else
    error('MasterRulers_PseudoOCV.mat not found. Please run Generate_MasterRulers_PseudoOCV.m first.');
end

num_seg_chg = length(MR_Chg) - 1;
num_seg_dch = length(MR_Dch) - 1;

%% 3. Load RPT VQ Grid Data
fprintf('Loading master file: RPT_VQ_grid.mat...\n');
vqData = load(vqGridFile, 'RPT_VQ_grid');
gridData = vqData.RPT_VQ_grid;
cycles = fieldnames(gridData);

%% 4. Initialize Feature Matrix Table
varNames = {'CellID', 'CycleNum', 'Condition', 'Static_Capacity'};
for i = 1:num_seg_chg, varNames{end+1} = sprintf('dQ_c_%02d', i); end
for i = 1:num_seg_dch, varNames{end+1} = sprintf('dQ_d_%02d', i); end
varNames = [varNames, {'C_eff_chg', 'C_eff_dch'}];

varTypes = [repmat({'string'}, 1, 3), repmat({'double'}, 1, 1 + num_seg_chg + num_seg_dch + 2)];
FM = table('Size', [0, length(varNames)], 'VariableTypes', varTypes, 'VariableNames', varNames);

%% 5. Extract Features and Labels
fprintf('Extracting Fragmented Features and Labels...\n');

conditions = {'c01', 'c05', 'c1', 'c2', 'c3'};
rowIdx = 1;

for cyc_idx = 1:length(cycles)
    cycStr = cycles{cyc_idx};
    cycNum = str2double(replace(cycStr, 'cyc', ''));
    cells = fieldnames(gridData.(cycStr));
    
    for cell_idx = 1:length(cells)
        chStr = cells{cell_idx};
        cellData = gridData.(cycStr).(chStr);
        
        % 1. Extract Target Label: Static Capacity
        if ~isfield(cellData, 'Static')
            continue; % Missing Static label, cannot use this row
        end
        % Max capacity swing within Static test
        static_cap = abs(max(cellData.Static.Q_raw) - min(cellData.Static.Q_raw));
        
        % 2. Extract Features per C-rate Condition
        for c = 1:length(conditions)
            cond = conditions{c};
            chg_field = [cond '_charge'];
            dch_field = [cond '_discharge'];
            
            if isfield(cellData, chg_field) && isfield(cellData, dch_field)
                data_c = cellData.(chg_field);
                data_d = cellData.(dch_field);
                
                % Meta calculation (Effective C-rate properly bridged)
                Nominal_Capacity = 64.0; % Actual Nominal Capacity of the ESS cell
                try
                    if isduration(data_c.t_raw)
                        time_c = hours(data_c.t_raw(end) - data_c.t_raw(1));
                    else
                        time_c = (data_c.t_raw(end) - data_c.t_raw(1)) / 3600;
                    end
                    cap_c = abs(data_c.Q_raw(end) - data_c.Q_raw(1));
                    C_eff_chg = (cap_c / time_c) / Nominal_Capacity;
                    
                    if isduration(data_d.t_raw)
                        time_d = hours(data_d.t_raw(end) - data_d.t_raw(1));
                    else
                        time_d = (data_d.t_raw(end) - data_d.t_raw(1)) / 3600;
                    end
                    cap_d = abs(data_d.Q_raw(end) - data_d.Q_raw(1));
                    C_eff_dch = (cap_d / time_d) / Nominal_Capacity;
                catch
                    C_eff_chg = NaN; C_eff_dch = NaN;
                end
                
                % Fragmented dQ Extraction
                dQ_c = extract_fragmented_dQ(data_c.V_raw, data_c.Q_raw, MR_Chg);
                dQ_d = extract_fragmented_dQ(data_d.V_raw, data_d.Q_raw, MR_Dch);
                
                % Append to Table
                row = {string(chStr), cycNum, string(cond), static_cap};
                row = [row, num2cell(dQ_c), num2cell(dQ_d), {C_eff_chg, C_eff_dch}];
                FM(rowIdx, :) = row;
                rowIdx = rowIdx + 1;
            end
        end
    end
end

%% 6. Save Feature Matrix
savePath = fullfile(saveDir, 'FeatureMatrix_ver0317.mat');
save(savePath, 'FM', 'MR_Chg', 'MR_Dch');
fprintf('Feature Extraction Complete. Saved to: %s\n', savePath);
fprintf('Total rows extracted: %d\n', height(FM));

%% Helper Function
function dQ = extract_fragmented_dQ(V, Q, MR)
    num_segs = length(MR) - 1;
    dQ = zeros(1, num_segs);
    for i = 1:num_segs
        v_start = MR(i);
        v_end = MR(i+1);
        
        idx = V >= v_start & V <= v_end;
        if any(idx)
            dQ(i) = abs(max(Q(idx)) - min(Q(idx)));
        else
            % 워크스루(기획) 요구사항: 0이 아닌 NaN을 반환하여 트리 모델의 Surrogate Split 유도
            dQ(i) = NaN; 
        end
    end
end
