clc; close all; clear all;

basePath = 'G:\.shortcut-targets-by-id\1L2QEuUmdCOdw7T7CnHpHMGCTi8D2y_4t\정철원\Materials\ESS_Data_Preprocessing';
savePath = fullfile(basePath, 'Dataset_Overview\Proposal_Fig');
ymFolder = '202206';
yFolder = '2022';
if ~exist(savePath, 'dir'), mkdir(savePath); end

start_date = datetime('20220601', 'InputFormat', 'yyyyMMdd');
n_days = 7;
Cnom = 1024;

T_all = []; V_all = []; I_all = []; SOC_all = []; C_rate_all = [];

for d = 0:n_days-1
    date = start_date + days(d);
    ymd = datestr(date, 'yyyymmdd');
    matFile = fullfile(basePath, yFolder, ymFolder, sprintf('Raw_%s.mat', ymd));
    if ~isfile(matFile)
        fprintf('파일 없음: %s\n', matFile);
        continue;
    end
    load(matFile);
    I = Raw.Online_DC_Current;
    T_all      = [T_all; Raw.sync_Time];
    V_all      = [V_all; Raw.Total_Average_CV_Sum];
    I_all      = [I_all; Raw.Online_DC_Current];
    SOC_all    = [SOC_all; Raw.Total_Average_SOC];
    C_rate_all = [C_rate_all; I/Cnom];
end

two_day_mask = T_all < start_date + days(2);
T_2d   = T_all(two_day_mask);
V_2d   = V_all(two_day_mask);
I_2d   = I_all(two_day_mask);
SOC_2d = SOC_all(two_day_mask);
C_rate_2d = C_rate_all(two_day_mask);

%% 1. Voltage-Current 1주일치 
fig1 = figure('Position', [100 100 1200 900]);
yyaxis left
plot(T_all, V_all, 'Color', '#0073C2'); ylabel('Voltage [V]','FontSize',14,'FontWeight','bold');
yyaxis right
plot(T_all, C_rate_all, 'Color', '#CD534C'); ylabel('C-rate [1/C]','FontSize',14,'FontWeight','bold');
xlabel('Time','FontSize',14,'FontWeight','bold'); title('Time - V, I (One Week)','FontSize',16,'FontWeight','bold'); grid on;
saveas(fig1, fullfile(savePath, 'Voltage_Current_1week.png'));
savefig(fig1, fullfile(savePath, 'Voltage_Current_1week.fig'));

%% 2. Voltage-Current 2일치 
fig2 = figure('Position', [100 100 1200 900]);
yyaxis left
plot(T_2d, V_2d, 'Color', '#0073C2'); ylabel('Voltage [V]','FontSize',14,'FontWeight','bold');
yyaxis right
plot(T_2d, C_rate_2d, 'Color', '#CD534C'); ylabel('C-rate [1/C]','FontSize',14);
xlabel('Time','FontSize',14,'FontWeight','bold'); title('Time - V, I (2 Days)','FontSize',16,'FontWeight','bold'); grid on;
saveas(fig2, fullfile(savePath, 'Voltage_Current_2days.png'));
savefig(fig2, fullfile(savePath, 'Voltage_Current_2days.fig'));

%% 3. SOC - 1주일
fig3 = figure('Position', [100 100 1200 900]);
plot(T_all, SOC_all, 'Color', '#747678','Linewidth',2); 
ylabel('SOC [%]','FontSize',14,'FontWeight','bold'); xlabel('Time','FontSize',14,'FontWeight','bold');
title('Time - SOC (One Week)','FontSize',16,'FontWeight','bold'); grid on;
saveas(fig3, fullfile(savePath, 'SOC_1week.png'));
savefig(fig3, fullfile(savePath, 'SOC_1week.fig'));

%% 4. SOC - 2일
fig4 = figure('Position', [100 100 1200 900]);
plot(T_2d, SOC_2d, 'Color', '#747678','LineWidth',2); 
ylabel('SOC [%]','FontSize',14,'FontWeight','bold'); xlabel('Time','FontSize',14,'FontWeight','bold');
title('Time - SOC (2 Days)','FontSize',16,'FontWeight','bold'); grid on;
saveas(fig4, fullfile(savePath, 'SOC_2days.png'));
savefig(fig4, fullfile(savePath, 'SOC_2days.fig'));
