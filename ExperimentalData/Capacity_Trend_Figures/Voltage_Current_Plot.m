%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental Data Voltage & Current Plot Code
% Required .csv files
% RPT / Drive Cycle / Aging 
% Updated date: 2025-07-27 
% Updated profiles: RPT - Drive Cycle - Aging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear; close all;

% Folder directory
baseFolder       = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data';
rptFolder        = fullfile(baseFolder, 'RPT\csv');
driveCycleFolder = fullfile(baseFolder, 'Drive Cycle\csv');
agingFolder      = fullfile(baseFolder, 'Aging\csv');
saveFolder       = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\DOE_Profile';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end     


% Channel
ch = '13';

% Fig color
color_voltage = '#CD534C';
color_current = '#0073C2';

% RPT0 csv
rptFile0 = fullfile(rptFolder, sprintf('Ch%s_RPT_200cyc.csv', ch));
if ~isfile(rptFile0)
    error('RPT0 파일 없음: %s', rptFile0);
end
T_rpt0 = readtable(rptFile0, 'VariableNamingRule', 'preserve');

% Drive Cycle csv
driveCyclePattern = sprintf('Ch%s*200cyc.csv', ch);
driveCycleFiles = dir(fullfile(driveCycleFolder, driveCyclePattern));
if isempty(driveCycleFiles)
    error('Drive Cycle 파일 없음: %s', driveCyclePattern);
end
driveCycleFile = fullfile(driveCycleFolder, driveCycleFiles(1).name);
T_driveCycle = readtable(driveCycleFile, 'VariableNamingRule', 'preserve');

% Aging csv
agingPattern = sprintf('Ch%s*0to200cyc*.csv', ch);
agingFiles = dir(fullfile(agingFolder, agingPattern));
if isempty(agingFiles)
    error('Aging 파일 없음: %s', agingPattern);
end
agingFile = fullfile(agingFolder, agingFiles(1).name);
T_aging = readtable(agingFile, 'VariableNamingRule', 'preserve');

% --- 데이터 추출 (Total Time: 6열, Current: 7열, Voltage: 8열)
% RPT0 데이터
time_rpt0 = T_rpt0{:,6};      % Total Time (duration)
voltage_rpt0 = T_rpt0{:,8};   % Voltage(V)
current_rpt0 = T_rpt0{:,7};   % Current(A)

% Drive Cycle 데이터
time_driveCycle = T_driveCycle{:,6};      % Total Time (duration)
voltage_driveCycle = T_driveCycle{:,8};   % Voltage(V)
current_driveCycle = T_driveCycle{:,7};   % Current(A)

% Aging 데이터
time_aging = T_aging{:,6};      % Total Time (duration)
voltage_aging = T_aging{:,8};   % Voltage(V)
current_aging = T_aging{:,7};   % Current(A)



% --- Plot

% RPT0 plot
figure;
yyaxis left
plot(time_rpt0, voltage_rpt0, 'LineWidth', 1.5);
ylabel('Voltage [V]', 'FontWeight', 'bold');
ylim([2.5 4.3]);
yyaxis right
plot(time_rpt0, current_rpt0, 'LineWidth', 1.5);
ylabel('Current [A]', 'FontWeight', 'bold');
xlabel('Time [h:mm:ss]', 'FontWeight', 'bold');
title('RPT', 'FontWeight', 'bold');
grid on;
savefig(gcf, fullfile(saveFolder, sprintf('Ch%s_RPT0_Profile.fig', ch)));
saveas(gcf, fullfile(saveFolder, sprintf('Ch%s_RPT0_Profile.png', ch)));

% Drive Cycle plot
figure;
yyaxis left
plot(time_driveCycle, voltage_driveCycle, 'LineWidth', 1.5);
ylabel('Voltage [V]', 'FontWeight', 'bold');
ylim([2.5 4.3]);
yyaxis right
plot(time_driveCycle, current_driveCycle, 'LineWidth', 1.5);
ylabel('Current [A]', 'FontWeight', 'bold');
xlabel('Time [h:mm:ss]', 'FontWeight', 'bold');
title('Drive Cycle', 'FontWeight', 'bold');
grid on;
savefig(gcf, fullfile(saveFolder, sprintf('Ch%s_DriveCycle_Profile.fig', ch)));
saveas(gcf, fullfile(saveFolder, sprintf('Ch%s_DriveCycle_Profile.png', ch)));

% Aging plot
figure;
yyaxis left
plot(time_aging, voltage_aging, 'LineWidth', 1.5);
ylabel('Voltage [V]', 'FontWeight', 'bold');
ylim([2.5 4.3]);
yyaxis right
plot(time_aging, current_aging, 'LineWidth', 1.5);
ylabel('Current [A]', 'FontWeight', 'bold');
xlabel('Time [h:mm:ss]', 'FontWeight', 'bold');
title('Aging', 'FontWeight', 'bold');
grid on;
savefig(gcf, fullfile(saveFolder, sprintf('Ch%s_Aging_Profile.fig', ch)));
saveas(gcf, fullfile(saveFolder, sprintf('Ch%s_Aging_Profile.png', ch)));
