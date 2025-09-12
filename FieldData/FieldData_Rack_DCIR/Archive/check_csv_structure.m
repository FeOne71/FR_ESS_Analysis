% Check CSV structure
clc; clear;

ambTempPath = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\JeongEupSi_MonthlyAvg_AmbTemp.csv';

if exist(ambTempPath, 'file')
    fprintf('CSV file found!\n');
    ambTempData = readtable(ambTempPath);
    
    fprintf('CSV Structure:\n');
    fprintf('  Rows: %d\n', height(ambTempData));
    fprintf('  Columns: %d\n', width(ambTempData));
    fprintf('  Column names: %s\n', strjoin(ambTempData.Properties.VariableNames, ', '));
    
    fprintf('\nFirst 5 rows:\n');
    disp(ambTempData(1:5, :));
    
    fprintf('\nAll data:\n');
    disp(ambTempData);
    
    % Check if it's numeric
    fprintf('\nData types:\n');
    for i = 1:width(ambTempData)
        col_name = ambTempData.Properties.VariableNames{i};
        col_data = ambTempData.(col_name);
        fprintf('  %s: %s\n', col_name, class(col_data));
    end
    
else
    fprintf('CSV file not found at: %s\n', ambTempPath);
end 