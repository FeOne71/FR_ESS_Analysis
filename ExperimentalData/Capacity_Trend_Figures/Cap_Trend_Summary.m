%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental Data Capacity Trend Summary Code
% Required .csv files
% RPT / Aging 
% Updated date: 2025-07-27 
% Updated profiles: RPT1(BOL) - Aging (0to200) - [[RPT2]] - Aging (201 to Max)
%                                                           ^^^^^^^^^^^^^^^^^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear; close all;

% Folder directory
dataDir   = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data';
rptFolder    = fullfile(dataDir, 'RPT');
agingFolder  = fullfile(dataDir, 'Aging');
saveFolder   = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Capacity_Trend_Figures';     


% Channel
ch = '09';

% Fig color
color_static = '#CD534C';
color_ocv    = '#0073C2';
color_aging  = '#20854E';
color_rpt200 = '#EE7733';

% RPT0 csv
rptFile0 = fullfile(rptFolder, sprintf('Ch%s_RPT_0cyc.csv', ch));
if ~isfile(rptFile0)
    error('RPT0 파일 없음: %s', rptFile0);
end
T_rpt0 = readmatrix(rptFile0);

% --- RPT0 Static Capacity (StepIdx=3, CycleIdx=2)
mask_static0 = (T_rpt0(:,2) == 3) & (T_rpt0(:,4) == 2);
Q_static0 = T_rpt0(mask_static0, 9);
Q_static_end0 = Q_static0(end);

% --- RPT0 OCV (StepIdx=10, CycleIdx=2)
mask_ocv0 = (T_rpt0(:,2) == 10) & (T_rpt0(:,4) == 2);
Q_ocv0 = T_rpt0(mask_ocv0, 9);
Q_ocv_end0 = Q_ocv0(end);

% RPT200 csv
rptFile200 = fullfile(rptFolder, sprintf('Ch%s_RPT_200cyc.csv', ch));
if ~isfile(rptFile200)
    error('RPT200 파일 없음: %s', rptFile200);
end
T_rpt200 = readmatrix(rptFile200);

% --- RPT200 Static Capacity (StepIdx=3, CycleIdx=2)
mask_static200 = (T_rpt200(:,2) == 3) & (T_rpt200(:,4) == 2);
Q_static200 = T_rpt200(mask_static200, 9);
Q_static_end200 = Q_static200(end);

% --- RPT200 OCV (StepIdx=10, CycleIdx=2)
mask_ocv200 = (T_rpt200(:,2) == 10) & (T_rpt200(:,4) == 2);
Q_ocv200 = T_rpt200(mask_ocv200, 9);
Q_ocv_end200 = Q_ocv200(end);

% RPT400 csv
rptFile400 = fullfile(rptFolder, sprintf('Ch%s_RPT_400cyc.csv', ch));
if ~isfile(rptFile400)
    error('RPT400 파일 없음: %s', rptFile400);
end
T_rpt400 = readmatrix(rptFile400);

% --- RPT400 Static Capacity (StepIdx=3, CycleIdx=2)
mask_static400 = (T_rpt400(:,2) == 3) & (T_rpt400(:,4) == 2);
Q_static400 = T_rpt400(mask_static400, 9);
Q_static_end400 = Q_static400(end);

% --- RPT400 OCV (StepIdx=10, CycleIdx=2)
mask_ocv400 = (T_rpt400(:,2) == 10) & (T_rpt400(:,4) == 2);
Q_ocv400 = T_rpt400(mask_ocv400, 9);
Q_ocv_end400 = Q_ocv400(end);

% --- Aging Capacity (0to200cyc, 200to400cyc, 400to600cyc)
% Aging 0to200cyc csv
agingFolder0 = fullfile(agingFolder, '0to200cyc');
agingPattern0 = sprintf('Ch%s*0to200cyc*.csv', ch);
agingFiles0 = dir(fullfile(agingFolder0, agingPattern0));
if isempty(agingFiles0)
    error('Aging 0to200cyc 파일 없음: %s', agingPattern0);
end
agingFile0 = fullfile(agingFolder0, agingFiles0(1).name);
T_aging0 = readmatrix(agingFile0);

% Aging 200to400cyc csv
agingFolder200 = fullfile(agingFolder, '200to400cyc');
agingPattern200 = sprintf('Ch%s*200to400cyc*.csv', ch);
agingFiles200 = dir(fullfile(agingFolder200, agingPattern200));
if isempty(agingFiles200)
    error('Aging 200to400cyc 파일 없음: %s', agingPattern200);
end
agingFile200 = fullfile(agingFolder200, agingFiles200(1).name);
T_aging200 = readmatrix(agingFile200);

% Aging 400to600cyc csv
agingFolder400 = fullfile(agingFolder, '400to600cyc');
agingPattern400 = sprintf('Ch%s*400to600cyc*.csv', ch);
agingFiles400 = dir(fullfile(agingFolder400, agingPattern400));
if isempty(agingFiles400)
    error('Aging 400to600cyc 파일 없음: %s', agingPattern400);
end
agingFile400 = fullfile(agingFolder400, agingFiles400(1).name);
T_aging400 = readmatrix(agingFile400);

% --- Aging 0to200cyc 용량 추출 (StepIdx=3, CycleIdx=1~200)
Q_aging0_end = NaN(1,200);
for cyc = 1:200
    mask = (T_aging0(:,2) == 3) & (T_aging0(:,4) == cyc);
    Q_tmp = T_aging0(mask, 9);
    if ~isempty(Q_tmp)
        Q_aging0_end(cyc) = Q_tmp(end);
    end
end

% --- Aging 200to400cyc 용량 추출 (StepIdx=3, CycleIdx=1~실제마지막사이클)
% 실제 마지막 사이클 찾기
aging200_cycles = unique(T_aging200(:,4));
aging200_cycles = aging200_cycles(aging200_cycles >= 1);  % 1 이상 사이클만
max_cycle_aging200 = max(aging200_cycles);

Q_aging200_end = NaN(1,200);
for cyc = 1:max_cycle_aging200
    mask = (T_aging200(:,2) == 3) & (T_aging200(:,4) == cyc);
    Q_tmp = T_aging200(mask, 9);
    if ~isempty(Q_tmp)
        Q_aging200_end(cyc) = Q_tmp(end);
    end
end

% --- Aging 400to600cyc 용량 추출 (StepIdx=3, CycleIdx=1~실제마지막사이클)
% 실제 마지막 사이클 찾기
aging400_cycles = unique(T_aging400(:,4));
aging400_cycles = aging400_cycles(aging400_cycles >= 1);  % 1 이상 사이클만
max_cycle_aging400 = max(aging400_cycles);

Q_aging400_end = NaN(1,200);
for cyc = 1:max_cycle_aging400
    mask = (T_aging400(:,2) == 3) & (T_aging400(:,4) == cyc);
    Q_tmp = T_aging400(mask, 9);
    if ~isempty(Q_tmp)
        Q_aging400_end(cyc) = Q_tmp(end);
    end
end

%% Plot

% --- x축 & y축 값 (RPT0 + Aging0to200 + RPT200 + Aging200to400 + RPT400 + Aging400to600)
x_rpt0_static = 1;                % RPT0 Static
x_rpt0_ocv    = 2;                % RPT0 OCV
x_aging0      = linspace(3, 10, 200);  % Aging0: Cycle 1~200 압축해서 3~10 구간에 표시
x_rpt200_static = 11;             % RPT200 Static
x_rpt200_ocv  = 12;               % RPT200 OCV
x_aging200    = linspace(13, 20, 200);  % Aging200: 201~400 압축해서 13~20 구간에 표시
x_rpt400_static = 21;             % RPT400 Static
x_rpt400_ocv  = 22;               % RPT400 OCV
x_aging400    = linspace(23, 30, 200);  % Aging400: 401~600 압축해서 23~30 구간에 표시

% --- 전체 데이터 구성
x_vals = [x_rpt0_static, x_rpt0_ocv, x_aging0, x_rpt200_static, x_rpt200_ocv, x_aging200, x_rpt400_static, x_rpt400_ocv, x_aging400];
y_vals = [Q_static_end0, Q_ocv_end0, Q_aging0_end, Q_static_end200, Q_ocv_end200, Q_aging200_end, Q_static_end400, Q_ocv_end400, Q_aging400_end];





% --- Plot
fig = figure('Name', ['Channel ' ch ' - Discharge Capacity Trend'], 'Color', 'w');
hold on;
plot(x_vals(1), y_vals(1), 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'RPT0 Static', 'Color', color_static);
plot(x_vals(2), y_vals(2), 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'RPT0 OCV', 'Color', color_ocv);
plot(x_vals(3:202), y_vals(3:202), 'o', 'LineWidth', 1.2, 'MarkerSize', 3, 'DisplayName', 'Aging 0to200', 'Color', color_aging);
plot(x_vals(203), y_vals(203), 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'RPT200 Static', 'Color', color_static);
plot(x_vals(204), y_vals(204), 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'RPT200 OCV', 'Color', color_ocv);
plot(x_vals(205:404), y_vals(205:404), 'o', 'LineWidth', 1.2, 'MarkerSize', 3, 'DisplayName', sprintf('Aging 200to%d', max_cycle_aging200), 'Color', color_aging);
plot(x_vals(405), y_vals(405), 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'RPT400 Static', 'Color', color_static);
plot(x_vals(406), y_vals(406), 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'RPT400 OCV', 'Color', color_ocv);
plot(x_vals(407:end), y_vals(407:end), 'o', 'LineWidth', 1.2, 'MarkerSize', 3, 'DisplayName', sprintf('Aging 400to%d', max_cycle_aging400), 'Color', color_aging);

% --- Text 표시 (마커 위)
text(x_vals(1), y_vals(1)+1, sprintf('%.2f', y_vals(1)), 'HorizontalAlignment', 'center', 'FontSize', 25);
text(x_vals(2), y_vals(2)+1, sprintf('%.2f', y_vals(2)), 'HorizontalAlignment', 'center', 'FontSize', 25);
text(x_vals(3), y_vals(3)+1, sprintf('%.2f', y_vals(3)), 'HorizontalAlignment', 'center', 'FontSize', 25);
text(x_vals(202), y_vals(202)+1, sprintf('%.2f', y_vals(202)), 'HorizontalAlignment', 'center', 'FontSize', 25);
text(x_vals(203), y_vals(203)+1, sprintf('%.2f', y_vals(203)), 'HorizontalAlignment', 'center', 'FontSize', 25);
text(x_vals(204), y_vals(204)+1, sprintf('%.2f', y_vals(204)), 'HorizontalAlignment', 'center', 'FontSize', 25);

% --- Aging 200to400 첫 번째와 마지막 사이클 텍스트
aging200_valid = ~isnan(Q_aging200_end);
first_idx = find(aging200_valid, 1);
last_idx = find(aging200_valid, 1, 'last');
text(x_vals(205+first_idx-1), y_vals(205+first_idx-1)+1, sprintf('%.2f', y_vals(205+first_idx-1)), 'HorizontalAlignment', 'center', 'FontSize', 25);
% 마지막 사이클 - 1의 실제 X축 위치 계산
aging200_last_x_pos = 13 + (last_idx-2) * (20-13) / 199;
text(aging200_last_x_pos, y_vals(205+last_idx-2)+1, sprintf('%.2f', y_vals(205+last_idx-2)), 'HorizontalAlignment', 'center', 'FontSize', 25);

% --- RPT400 텍스트 표시
text(x_vals(405), y_vals(405)+1, sprintf('%.2f', y_vals(405)), 'HorizontalAlignment', 'center', 'FontSize', 25);
text(x_vals(406), y_vals(406)+1, sprintf('%.2f', y_vals(406)), 'HorizontalAlignment', 'center', 'FontSize', 25);

% --- Aging 400to600 첫 번째와 마지막 사이클 텍스트
aging400_valid = ~isnan(Q_aging400_end);
first_idx_400 = find(aging400_valid, 1);
last_idx_400 = find(aging400_valid, 1, 'last');
text(x_vals(407+first_idx_400-1), y_vals(407+first_idx_400-1)+1, sprintf('%.2f', y_vals(407+first_idx_400-1)), 'HorizontalAlignment', 'center', 'FontSize', 25);
% 마지막 사이클 - 1의 실제 X축 위치 계산
aging400_last_x_pos = 23 + (last_idx_400-2) * (30-23) / 199;
text(aging400_last_x_pos, y_vals(407+last_idx_400-2)+1, sprintf('%.2f', y_vals(407+last_idx_400-2)), 'HorizontalAlignment', 'center', 'FontSize', 25);

% --- X축 설정
xlim([0 31]);
% Aging 200to400의 실제 마지막 사이클 위치 계산
aging200_last_x = 13 + (max_cycle_aging200-1) * (20-13) / 199;  % 13~20 구간에서 실제 마지막 사이클 위치
% Aging 400to600의 실제 마지막 사이클 위치 계산
aging400_last_x = 23 + (max_cycle_aging400-1) * (30-23) / 199;  % 23~30 구간에서 실제 마지막 사이클 위치
xticks([1 2 3 10 11 12 13 aging200_last_x 21 22 23 aging400_last_x]);
xticklabels({'RPT0 Static', 'RPT0 OCV', 'Aging 1', 'Aging 200', 'RPT200 Static', 'RPT200 OCV', 'Aging 201', sprintf('Aging %d', max_cycle_aging200), 'RPT400 Static', 'RPT400 OCV', 'Aging 401', sprintf('Aging %d', max_cycle_aging400)});

xlabel('Test Type / Cycle Index', 'FontSize', 20);
ylabel('Discharged Capacity [Ah]', 'FontSize', 20);
ylim([40 66]);
title(sprintf('Channel %s - Capacity Trend', ch), 'FontSize', 20);
legend('Location', 'east');
grid on;

% fig = figure('Name', ['Channel ' ch ' - Discharge Capacity Trend'], 'Color', 'w');
% hold on;
% 
% % --- Plot Points
% h1 = plot(x_vals(1), y_vals(1), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5, ...
%     'DisplayName', sprintf('Static (%.2f Ah)', y_vals(1)), 'Color', color_static);
% h2 = plot(x_vals(2), y_vals(2), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5, ...
%     'DisplayName', sprintf('OCV (%.2f Ah)', y_vals(2)), 'Color', color_ocv);
% h3 = plot(x_vals(3), y_vals(3), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5, ...
%     'DisplayName', sprintf('Aging Start (%.2f Ah)', y_vals(3)), 'Color', color_aging);
% h4 = plot(x_vals(end), y_vals(end), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5, ...
%     'DisplayName', sprintf('Aging End (%.2f Ah)', y_vals(end)), 'Color', color_aging);
% 
% % --- 연결선 (Static → OCV → Aging Start → Aging End)
% x_line = [x_vals(1), x_vals(2), x_vals(3), x_vals(end)];
% y_line = [y_vals(1), y_vals(2), y_vals(3), y_vals(end)];
% plot(x_line, y_line, '-', 'Color', '#888888', 'LineStyle', '--');
% 
% % --- 전체 Aging plot
% plot(x_vals(3:end), y_vals(3:end), '-', ...
%     'Color', color_aging, 'LineWidth', 1.2, 'DisplayName', 'Aging Capacity');
% 
% % --- Text 표시 (마커 위)
% text(x_vals(1), y_vals(1)+0.1, sprintf('%.2f', y_vals(1)), 'HorizontalAlignment', 'center');
% text(x_vals(2), y_vals(2)+0.1, sprintf('%.2f', y_vals(2)), 'HorizontalAlignment', 'center');
% text(x_vals(3), y_vals(3)+0.1, sprintf('%.2f', y_vals(3)), 'HorizontalAlignment', 'center');
% text(x_vals(end), y_vals(end)+0.1, sprintf('%.2f', y_vals(end)), 'HorizontalAlignment', 'center');
% 
% % --- X축 설정
% xlim([0 11]);
% xticks([1 2 3 10]);
% xticklabels({'Static', 'OCV', 'Aging Start', 'Aging End'});
% 
% xlabel('Test Type / Cycle Index');
% ylabel('Discharge Capacity [Ah]');
% title(sprintf('Channel %s - Capacity Trend', ch));
% legend('Location', 'best');
% grid on;

% --- 저장

savefig(fig, fullfile(saveFolder, sprintf('Ch%s_Capacity_Trend_600cyc.fig', ch)));
