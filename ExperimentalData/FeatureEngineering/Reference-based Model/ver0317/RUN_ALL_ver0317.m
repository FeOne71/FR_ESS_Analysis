% RUN_ALL_ver0317.m
% 3개의 스크립트를 순차 실행하고 결과를 텍스트 파일로 자동 Export합니다.

verDir0 = fileparts(mfilename('fullpath'));
cd(verDir0);

%% Step 1: Generate Master Rulers
fprintf('=== [1/3] Generate Master Rulers ===\n');
run('Generate_MasterRulers_PseudoOCV.m');

% verDir 복원 (자식 스크립트의 clear로 날아가므로)
verDir0 = fileparts(mfilename('fullpath'));
mr = load(fullfile(verDir0, 'MasterRulers_PseudoOCV.mat'));
fid = fopen(fullfile(verDir0, 'MasterRuler_Results.txt'), 'w');
fprintf(fid, '=== Master Ruler Results (ver0317) ===\n');
vbc = mr.MasterRuler_ver0317.V_bounds_chg;
vbd = mr.MasterRuler_ver0317.V_bounds_dch;
fprintf(fid, 'Charge Segments: %d\n', length(vbc)-1);
fprintf(fid, 'Discharge Segments: %d\n\n', length(vbd)-1);
fprintf(fid, 'Charge Boundaries (V): '); fprintf(fid, '%.4f ', vbc); fprintf(fid, '\n');
for i = 1:length(vbc)-1
    fprintf(fid, '  Seg%02d: %.4fV ~ %.4fV (width=%.4fV)\n', i, vbc(i), vbc(i+1), vbc(i+1)-vbc(i));
end
fprintf(fid, '\nDischarge Boundaries (V): '); fprintf(fid, '%.4f ', vbd); fprintf(fid, '\n');
for i = 1:length(vbd)-1
    fprintf(fid, '  Seg%02d: %.4fV ~ %.4fV (width=%.4fV)\n', i, vbd(i), vbd(i+1), vbd(i+1)-vbd(i));
end
fclose(fid);
fprintf('Master Ruler results exported to MasterRuler_Results.txt\n\n');

%% Step 2: Feature Extraction
fprintf('=== [2/3] RPT Feature Extractor ===\n');
run('RPT_FeatureExtractor_ver0317.m');

verDir0 = fileparts(mfilename('fullpath'));
fm = load(fullfile(verDir0, 'FeatureMatrix_ver0317.mat'));
fid = fopen(fullfile(verDir0, 'FeatureMatrix_Summary.txt'), 'w');
fprintf(fid, '=== Feature Matrix Summary ===\n');
fprintf(fid, 'Total Rows: %d\n', height(fm.FM));
fprintf(fid, 'Total Columns: %d\n', width(fm.FM));
fprintf(fid, 'Column Names:\n');
for i = 1:width(fm.FM)
    fprintf(fid, '  %d. %s\n', i, fm.FM.Properties.VariableNames{i});
end
fprintf(fid, '\nUnique Channels: ');
fprintf(fid, '%s ', unique(fm.FM.CellID));
fprintf(fid, '\nUnique Conditions: ');
conds = unique(fm.FM.Condition);
for i = 1:numel(conds)
    fprintf(fid, '%s ', string(conds(i)));
end
fprintf(fid, '\nNaN Count per Feature:\n');
for i = 1:width(fm.FM)
    if isnumeric(fm.FM{:,i})
        nan_cnt = sum(isnan(fm.FM{:,i}));
        if nan_cnt > 0
            fprintf(fid, '  %s: %d NaNs\n', fm.FM.Properties.VariableNames{i}, nan_cnt);
        end
    end
end
fclose(fid);
fprintf('Feature Matrix summary exported to FeatureMatrix_Summary.txt\n\n');

%% Step 3: Correlation Analysis
fprintf('=== [3/3] ML_00 Correlation Analysis ===\n');
run('ML_00_Correlation_Analysis.m');

fprintf('\n=== ALL DONE! ===\n');
