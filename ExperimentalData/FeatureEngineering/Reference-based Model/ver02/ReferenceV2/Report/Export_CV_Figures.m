% Export_CV_Figures.m
clear; clc; close all;

base_dir = fullfile('D:', 'JCW', 'Projects', 'KEPCO_ESS_Local', 'ExperimentalData', 'FeatureEngineering', 'Lab_RPT_Analysis', 'ver02');
rf_dir = fullfile(base_dir, 'ML_RandomForest');
save_dir = fullfile(base_dir, 'Report', 'Figures');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

figs_to_export = {
    'RF_ParityPlot_Group_KFold.fig', 'Report_CV_ParityPlot.png';
    'RF_ErrorDist_Group_KFold.fig', 'Report_CV_ErrorDist.png';
    'RF_Importance_Group_KFold.fig', 'Report_CV_FeatureImportance.png';
    'RF_Hyperparams.fig', 'Report_CV_Hyperparams.png'
};

for i = 1:size(figs_to_export, 1)
    fig_path = fullfile(rf_dir, figs_to_export{i,1});
    if exist(fig_path, 'file')
        f = openfig(fig_path, 'invisible');
        % Optional: Adjust figure size for better resolution
        set(f, 'Position', [100, 100, 1000, 400]);
        
        save_path = fullfile(save_dir, figs_to_export{i,2});
        saveas(f, save_path);
        close(f);
        fprintf('Exported: %s\n', figs_to_export{i,2});
    else
        fprintf('Not found: %s\n', fig_path);
    end
end
