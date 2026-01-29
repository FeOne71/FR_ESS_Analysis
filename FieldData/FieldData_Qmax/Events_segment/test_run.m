% Test run with verbose output
clear; clc;

% Test with one file
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
typePath = fullfile(dataDir, 'New');
yearPath = fullfile(typePath, '2023');
monthPath = fullfile(yearPath, '202310');
matFiles = dir(fullfile(monthPath, '*.mat'));

fprintf('Found %d files in %s\n', length(matFiles), monthPath);

if length(matFiles) > 0
    matFilePath = fullfile(monthPath, matFiles(1).name);
    fprintf('Testing file: %s\n', matFiles(1).name);
    
    % Extract date
    [~, filename_only, ~] = fileparts(matFiles(1).name);
    dateMatch = regexp(filename_only, '\d{8}', 'once', 'match');
    if ~isempty(dateMatch)
        fileDate_str = sprintf('%s-%s-%s', dateMatch(1:4), dateMatch(5:6), dateMatch(7:8));
        fprintf('Date: %s\n', fileDate_str);
        
        % Check if RPT date
        RPT_dates = {'2021-06-03', '2023-10-16', '2024-09-09', '2025-07-11'};
        RPT_dateSet = containers.Map();
        for i = 1:length(RPT_dates)
            RPT_dateSet(RPT_dates{i}) = true;
        end
        
        if isKey(RPT_dateSet, fileDate_str)
            fprintf('This is an RPT date - will be skipped\n');
        else
            fprintf('Not an RPT date - will be processed\n');
            
            % Try loading
            try
                S = load(matFilePath);
                fprintf('File loaded successfully\n');
                if isfield(S, 'Raw')
                    fprintf('Raw field exists\n');
                    rackData = S.Raw;
                    fprintf('Fields in Raw: %s\n', strjoin(fieldnames(rackData), ', '));
                end
            catch ME
                fprintf('Error loading: %s\n', ME.message);
            end
        end
    end
end
