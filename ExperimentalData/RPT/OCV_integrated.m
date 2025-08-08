%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OCV Integration - ECM-Ready OCV-SOC Function (8셀 평균)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% File Directory
folderPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT';
summaryFolder = fullfile(folderPath, 'Summary');
ocvMatFile = fullfile(summaryFolder, 'OCV_interpolation_functions.mat');

%% Load OCV Data
fprintf('=== ECM-Ready OCV Integration Started ===\n');
fprintf('기존 OCV 파일 로드 중...\n');
load(ocvMatFile);

%% Channels (0cyc and 200cyc data for integration)
channelNames_0cyc = {'ch9_0cyc', 'ch10_0cyc', 'ch11_0cyc', 'ch12_0cyc', ...
                     'ch13_0cyc', 'ch14_0cyc', 'ch15_0cyc', 'ch16_0cyc'};
channelNames_200cyc = {'ch9_200cyc', 'ch10_200cyc', 'ch11_200cyc', 'ch12_200cyc', ...
                       'ch13_200cyc', 'ch14_200cyc', 'ch15_200cyc', 'ch16_200cyc'};

% 2. SOC grid (모든 셀 동일)
SOC_grid = OCV_struct.ch9_0cyc.soc_grid; % 1x101, 0:1:100

% 3. 0cyc: 8개 셀의 OCV를 8x101 행렬로 쌓기
all_V_OCV_0cyc = [];
for i = 1:length(channelNames_0cyc)
    vocv = OCV_struct.(channelNames_0cyc{i}).V_OCV; % 1x101
    all_V_OCV_0cyc = [all_V_OCV_0cyc; vocv];
end

% 4. 0cyc: SOC별 평균 OCV
V_avg_SOC_0cyc = mean(all_V_OCV_0cyc, 1, 'omitnan'); % 1x101

% 5. 0cyc: 함수화
OCV_SOC_func_0cyc = @(soc_query) interp1(SOC_grid, V_avg_SOC_0cyc, soc_query, 'linear');

% 6. 0cyc: 구조체 저장
OCV_integrated_struct_0cyc.SOC_grid = SOC_grid;
OCV_integrated_struct_0cyc.V_avg_SOC = V_avg_SOC_0cyc;
OCV_integrated_struct_0cyc.OCV_SOC_func = OCV_SOC_func_0cyc; 

% 0cyc: 평균 용량 계산 및 저장
cap_chg_all_0cyc = zeros(1, length(channelNames_0cyc));
cap_dch_all_0cyc = zeros(1, length(channelNames_0cyc));
for i = 1:length(channelNames_0cyc)
    cap_chg_all_0cyc(i) = OCV_struct.(channelNames_0cyc{i}).capacity_chg;
    cap_dch_all_0cyc(i) = OCV_struct.(channelNames_0cyc{i}).capacity_dch;
end
mean_cap_chg_0cyc = mean(cap_chg_all_0cyc);
mean_cap_dch_0cyc = mean(cap_dch_all_0cyc);
mean_cap_0cyc = mean([cap_chg_all_0cyc, cap_dch_all_0cyc]);

OCV_integrated_struct_0cyc.mean_capacity_chg = mean_cap_chg_0cyc;
OCV_integrated_struct_0cyc.mean_capacity_dch = mean_cap_dch_0cyc;
OCV_integrated_struct_0cyc.mean_capacity = mean_cap_0cyc;

% 7. 200cyc: 8개 셀의 OCV를 8x101 행렬로 쌓기
all_V_OCV_200cyc = [];
for i = 1:length(channelNames_200cyc)
    vocv = OCV_struct.(channelNames_200cyc{i}).V_OCV; % 1x101
    all_V_OCV_200cyc = [all_V_OCV_200cyc; vocv];
end

% 8. 200cyc: SOC별 평균 OCV
V_avg_SOC_200cyc = mean(all_V_OCV_200cyc, 1, 'omitnan'); % 1x101

% 9. 200cyc: 함수화
OCV_SOC_func_200cyc = @(soc_query) interp1(SOC_grid, V_avg_SOC_200cyc, soc_query, 'linear');

% 10. 200cyc: 구조체 저장
OCV_integrated_struct_200cyc.SOC_grid = SOC_grid;
OCV_integrated_struct_200cyc.V_avg_SOC = V_avg_SOC_200cyc;
OCV_integrated_struct_200cyc.OCV_SOC_func = OCV_SOC_func_200cyc; 

% 200cyc: 평균 용량 계산 및 저장
cap_chg_all_200cyc = zeros(1, length(channelNames_200cyc));
cap_dch_all_200cyc = zeros(1, length(channelNames_200cyc));
for i = 1:length(channelNames_200cyc)
    cap_chg_all_200cyc(i) = OCV_struct.(channelNames_200cyc{i}).capacity_chg;
    cap_dch_all_200cyc(i) = OCV_struct.(channelNames_200cyc{i}).capacity_dch;
end
mean_cap_chg_200cyc = mean(cap_chg_all_200cyc);
mean_cap_dch_200cyc = mean(cap_dch_all_200cyc);
mean_cap_200cyc = mean([cap_chg_all_200cyc, cap_dch_all_200cyc]);

OCV_integrated_struct_200cyc.mean_capacity_chg = mean_cap_chg_200cyc;
OCV_integrated_struct_200cyc.mean_capacity_dch = mean_cap_dch_200cyc;
OCV_integrated_struct_200cyc.mean_capacity = mean_cap_200cyc;

% 시각화
figure('Name', 'Integrated SOC-OCV Curve', 'Color', 'w');
plot(SOC_grid, V_avg_SOC_0cyc, 'bo-', 'LineWidth', 2, 'DisplayName', 'RPT0 (0cyc)');
hold on;
plot(SOC_grid, V_avg_SOC_200cyc, 'ro-', 'LineWidth', 2, 'DisplayName', 'RPT200 (200cyc)');
xlabel('SOC [%]');
ylabel('OCV [V]');
title('Integrated SOC-OCV Curve (8-cell Average)');
grid on;
xlim([0 100]);
legend('Location', 'best');

% 저장
OCV_integrated_RPT0cyc = OCV_integrated_struct_0cyc;
OCV_integrated_RPT200cyc = OCV_integrated_struct_200cyc;
save('G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary\OCV_interpolation_functions.mat', 'OCV_integrated_RPT0cyc', 'OCV_integrated_RPT200cyc', '-append');