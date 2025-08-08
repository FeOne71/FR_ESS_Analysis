%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

folderPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\csv';
saveFolder = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT';

fileNames = {
    % 'Ch9_RPT_0cyc.csv';
    % 'Ch10_RPT_0cyc.csv';
    % 'Ch11_RPT_0cyc.csv';
    % 'Ch12_RPT_0cyc.csv';
    % 'Ch13_RPT_0cyc.csv';
    % 'Ch14_RPT_0cyc.csv';
    % 'Ch15_RPT_0cyc.csv';
    % 'Ch16_RPT_0cyc.csv';
    'Ch9_RPT_200cyc.csv';
    'Ch10_RPT_200cyc.csv';
    'Ch11_RPT_200cyc.csv';
    'Ch12_RPT_200cyc.csv';
    'Ch13_RPT_200cyc.csv';
    'Ch14_RPT_200cyc.csv';
    'Ch15_RPT_200cyc.csv';
    'Ch16_RPT_200cyc.csv'
};

validStepIdx = [1,2,3,4];
color_voltage = '#0073C2';  % blue
color_current = '#E18727';  % orange

for i = 1:length(fileNames)
    filename = fileNames{i};
    filepath = fullfile(folderPath, filename);

    [~, baseName, ~] = fileparts(filename);  
    channelLabel = extractBetween(baseName, 'Ch', '_');  
    channelName = sprintf('Channel %s', channelLabel{1});  % → Channel 09

    if exist(filepath, 'file') ~= 2
        fprintf('파일 없음: %s\n', filepath);
        continue
    end

    T = readtable(filepath);

    stepIndex   = T{:,2};  % Step Index
    stepType    = T{:,3};  % Step Type 
    cycleIndex  = T{:,4};  % Cycle Index 
    totalTime   = T{:,6};  % Total Time [sec]
    current     = T{:,7};  % Current [A]
    voltage     = T{:,8};  % Voltage [V]
    capacity    = T{:,9};  % Capacity [Ah]
    dQdV        = T{:,16}; % dQdV [mAh/V]

    mask = ismember(stepIndex,validStepIdx);
    if sum(mask) == 0
        fprintf('Step 1 데이터 없음: %s\n', filename);
        continue
    end

    % 구조 저장
    % parsed.stepIndex  = stepIndex(mask);
    % parsed.stepType   = stepType(mask);
    % parsed.cycleIndex = cycleIndex(mask);
    % parsed.totalTime  = totalTime(mask);
    % parsed.time_cap       = seconds(totalTime(mask) - totalTime(find(mask,1))); % 상대시간
    % parsed.time       = seconds(totalTime - totalTime(1));
    % parsed.current    = current(mask);
    % parsed.voltage    = voltage(mask);
    % parsed.capacity   = capacity(mask);

    parsed.stepIndex  = stepIndex;
    parsed.stepType   = stepType;
    parsed.cycleIndex = cycleIndex;
    parsed.totalTime  = totalTime;
    parsed.time_cap       = seconds(totalTime(mask) - totalTime(find(mask,1))); % 상대시간
    parsed.time       = seconds(totalTime - totalTime(1));
    parsed.current    = current;
    parsed.voltage    = voltage;
    parsed.capacity   = capacity;


    % 저장
    [~, baseName, ~] = fileparts(filename);
    matname = sprintf('%s_StaticCap_RPT200.mat', baseName);
    save(fullfile(saveFolder, matname), 'parsed');

%% === Figure ===
    figure('Name', baseName, 'NumberTitle', 'off');
    set(gcf, 'Position', [100 100 800 600]);  % [x, y, width, height]

    % 1. Static Capcaity Time vs V/I
    subplot(2,1,1);
    yyaxis left
    plot(parsed.time_cap, parsed.voltage(mask), 'Color', color_voltage, 'LineWidth', 1.5);
    ylabel('Voltage [V]');
    ylim([2.5 4.5]);
    % ylim([min(parsed.voltage)*0.98, max(parsed.voltage)*1.02]);

    yyaxis right
    plot(parsed.time_cap, parsed.current(mask), 'Color', color_current, 'LineWidth', 1.5);
    ylabel('Current [A]');
    ylim([-50 50]);
    % ylim([min(parsed.current)*1.1, max(parsed.current)*1.1]);

    xlabel('Time [s]');
    % title(sprintf('%s - Time vs Voltage/Current (Static Capacity)', baseName));
    title(sprintf('%s - Time vs Voltage/Current [Static Capacity (RPT2)]', channelName),'FontSize', 12);
    % title(sprintf('%s - Capacity vs Time', channelName), 'FontSize', 14);


    grid on;

    % 2. Static Capcaity V vs Q
    subplot(2,1,2); hold on;

    legendList = {};

    stepCycleList = [1 1; 3 1; 1 2; 3 2];  % [StepIndex, CycleIndex] 조합
    for k = 1:size(stepCycleList,1)
        s = stepCycleList(k,1);  % Step index
        c = stepCycleList(k,2);  % Cycle index
        color = '#0073C2';  % 기본은 충전 색
        linestyle = '-';   % Cycle 1: 실선, Cycle 2: 점선

        if s == 3  % 방전이면 빨강
            color = '#CD534C';
        end
        if c == 2
            linestyle = '--';
        end

        % 해당 조합의 데이터 필터
        mask_k = (parsed.stepIndex == s) & (parsed.cycleIndex == c);
        if any(mask_k)
            % 마지막 용량값
            Qend = parsed.capacity(find(mask_k,1,'last'));

            % 그래프 그리기
            plot(parsed.capacity(mask_k), parsed.voltage(mask_k), ...
                'LineStyle', linestyle, 'Color', color, 'LineWidth', 1.5);

            % legend 문자열 생성
            if s == 1
                action = 'Charge';
            elseif s == 3
                action = 'Discharge';
            else
                action = sprintf('Step %d', s);  % 그 외 StepIndex 대응
            end
            legendList{end+1} = sprintf('%s(Cycle %d) \n Q = %.2f Ah', ...
                action, c, Qend);            
        end
    end

    xlabel('Capacity [Ah]');
    ylabel('Voltage [V]');
    xlim([0 70]);
    ylim([2.5 4.5]);
    
    title(sprintf('%s - Voltage vs Capacity [Static Capacity (RPT2)]', channelName), 'FontSize',12);    
    legend(legendList, 'Location', 'best');
    grid on;

    figSaveFolder = fullfile(saveFolder, 'figures');
    if ~exist(figSaveFolder, 'dir')
        mkdir(figSaveFolder);
    end

    figName_fig = fullfile(figSaveFolder, sprintf('%s.fig', baseName));
    savefig(gcf, figName_fig);

    % 충전 Step
    figure('Name', [channelName ' - OCV vs SOC (RPT2)'], 'Color', 'w');
    
    mask_chg = (stepIndex == 8) & (cycleIndex == 2);
    V_chg = voltage(mask_chg);
    cap_chg = capacity(mask_chg);
    Qchg_nom = cap_chg(end);  
    soc_chg = 100 * (cap_chg / Qchg_nom);

    % 방전 Step
    mask_dch = (stepIndex == 10) & (cycleIndex == 2);
    V_dch = voltage(mask_dch);
    cap_dch = capacity(mask_dch);
    Qdch_nom = cap_dch(end);
    soc_dch = 100*(1 - (cap_dch / Qdch_nom));

    % 방전 구간 내림차순으로
    [soc_dch_sorted, idx] = sort(soc_dch, 'descend');
    % [V_dch_sorted, idx]   = sort(V_dch, 'descend');
    V_dch_sorted = V_dch(idx);

    figure;
    hold on;
    plot(soc_chg, V_chg, 'b-', 'DisplayName', 'Charge ', 'LineWidth', 1.5);
    plot(soc_dch_sorted, V_dch_sorted, 'r-', 'DisplayName', 'Discharge', 'LineWidth', 1.5);
    xlabel('SOC [%]');
    ylabel('Voltage [V]');
    legend('Location', 'Best');
    title(sprintf('%s - OCV vs SOC [RPT2]', channelName),'FontSize', 12);        
    grid on;

    figName_fig = fullfile(figSaveFolder, sprintf('%s_OCV_SOC.fig', baseName));
    savefig(gcf, figName_fig);


%% ========== OCV Curve Average & Interpolation ==========
    soc_grid = 0:1:100;

    % 보간
    ocv_chg_interp = interp1(soc_chg, V_chg, soc_grid, 'linear');
    ocv_dch_interp = interp1(soc_dch_sorted, V_dch_sorted, soc_grid, 'linear');

    ocv_avg = (ocv_chg_interp + ocv_dch_interp) / 2;

    OCV_func = @(soc_query) interp1(soc_grid, ocv_avg, soc_query, 'linear');

    OCV_struct.soc_grid = soc_grid;
    OCV_struct.ocv_charge = ocv_chg_interp;
    OCV_struct.ocv_discharge = ocv_dch_interp;
    OCV_struct.ocv_avg = ocv_avg;
    OCV_struct.OCV_func = OCV_func;

    matName_OCV = fullfile(saveFolder, sprintf('%s_OCV_interp.mat', baseName));
    save(matName_OCV, 'OCV_struct');
    % 
    % figure('Name', [channelName ' - Averaged OCV'], 'Color', 'w');
    % hold on;
    % plot(soc_grid, ocv_chg_interp, 'b--', 'DisplayName', 'Charge OCV');
    % plot(soc_grid, ocv_dch_interp, 'r--', 'DisplayName', 'Discharge OCV');
    % plot(soc_grid, ocv_avg, 'k-', 'LineWidth', 2, 'DisplayName', 'Avg OCV');
    % xlabel('SOC [%]');
    % ylabel('OCV [V]');
    % title(sprintf('%s - Averaged OCV-SOC Curve', channelName), 'FontSize', 12);
    % legend('Location', 'best');
    % grid on;
    % 
    % figName_ocv_avg = fullfile(figSaveFolder, sprintf('%s_OCV_interp.fig', baseName));
    % savefig(gcf, figName_ocv_avg);



end