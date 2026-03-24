% RUN_REMAINING_ver0317.m
% Step 2 결과 Export + Step 3 (Correlation) 실행
verDir0 = fileparts(mfilename('fullpath'));
cd(verDir0);

% Feature Matrix 요약 Export
fm = load(fullfile(verDir0, 'FeatureMatrix_ver0317.mat'));
fid = fopen(fullfile(verDir0, 'FeatureMatrix_Summary.txt'), 'w');
fprintf(fid, '=== Feature Matrix Summary ===\n');
fprintf(fid, 'Total Rows: %d\n', height(fm.FM));
fprintf(fid, 'Total Columns: %d\n', width(fm.FM));
fprintf(fid, 'Column Names:\n');
for i = 1:width(fm.FM)
    fprintf(fid, '  %d. %s\n', i, fm.FM.Properties.VariableNames{i});
end
fprintf(fid, '\nUnique CellIDs: ');
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
        fprintf(fid, '  %s: %d NaNs\n', fm.FM.Properties.VariableNames{i}, nan_cnt);
    end
end
fclose(fid);
fprintf('Feature Matrix summary exported to FeatureMatrix_Summary.txt\n\n');

%% Step 3: Correlation Analysis
fprintf('=== [3/3] ML_00 Correlation Analysis ===\n');
run('ML_00_Correlation_Analysis.m');

fprintf('\n=== ALL DONE! ===\n');
