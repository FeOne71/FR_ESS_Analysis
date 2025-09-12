clc; close all; clear all;

basePath = 'G:\.shortcut-targets-by-id\1L2QEuUmdCOdw7T7CnHpHMGCTi8D2y_4t\정철원\Materials\ESS_Data_Preprocessing';
savePath = fullfile(basePath, 'Dataset_Overview\Monthly_Overview');

ymFolders = {'202106','202206','202306'};
yFolders = {'2021','2022','2023'};
Cnom = 1024;

if ~exist(savePath, 'dir'), mkdir(savePath); end

for y = 1:length(yFolders)
    start_date = datetime([ymFolders{y}, '01'], 'InputFormat', 'yyyyMMdd');
    n_days = 30;

    T_all = []; V_all = []; I_all = []; SOC_all = []; C_rate_all = [];

    for d = 0:n_days-1
        date = start_date + days(d);
        ymd = datestr(date, 'yyyymmdd');
        matFile = fullfile(basePath, yFolders{y}, ymFolders{y}, sprintf('Raw_%s.mat', ymd));
        if ~isfile(matFile)
            fprintf('No Files: %s\n', matFile);
            continue;
        end
        load(matFile);
        I = Raw.Online_DC_Current;

        T_all      = [T_all; Raw.sync_Time];
        V_all      = [V_all; Raw.Total_Average_CV_Sum];
        I_all      = [I_all; I];
        SOC_all    = [SOC_all; Raw.Total_Average_SOC];
        C_rate_all = [C_rate_all; I / Cnom];
    end

    % 그래프 출력
    fig = figure('Position', [100 100 1200 900]);
    subplot(4,1,1);
    plot(T_all, V_all, 'Color', '#0073C2', 'LineWidth', 1.5); grid on;
    ylabel('Voltage [V]','FontSize',14,'FontWeight','bold');

    subplot(4,1,2);
    plot(T_all, C_rate_all, 'Color', '#E18727', 'LineWidth', 1.5); grid on;
    ylabel('C-rate [1/C]','FontSize',14,'FontWeight','bold');

    subplot(4,1,3);
    plot(T_all, SOC_all, 'Color', '#7E6148', 'LineWidth', 1.5); grid on;
    ylabel('SOC [%]','FontSize',14,'FontWeight','bold');

    subplot(4,1,4);
    plot(T_all, (I_all .* V_all)/1000, 'Color','#20854E', 'LineWidth', 1.5); grid on;
    xlabel('Time','FontSize',14,'FontWeight','bold');
    ylabel('Power [kW]','FontSize',14,'FontWeight','bold');

    sgtitle(sprintf('%s Monthly Data Overview', ymFolders{y}), 'FontSize', 16, 'FontWeight', 'bold');

    savefig(fig, fullfile(savePath, sprintf('%s_overview.fig', ymFolders{y})));
end