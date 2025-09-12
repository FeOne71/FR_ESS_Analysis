%% Static Capacity Summary

clc; close all; clear all;

% folderPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT';
% saveFolder = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary';
csvPath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\csv';
saveFolder = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary';


fileNames = {
    'Ch9_RPT_0cyc.csv';    'Ch9_RPT_200cyc.csv';
    'Ch10_RPT_0cyc.csv';   'Ch10_RPT_200cyc.csv';
    'Ch11_RPT_0cyc.csv';   'Ch11_RPT_200cyc.csv';
    'Ch12_RPT_0cyc.csv';   'Ch12_RPT_200cyc.csv';
    'Ch13_RPT_0cyc.csv';   'Ch13_RPT_200cyc.csv';
    'Ch14_RPT_0cyc.csv';   'Ch14_RPT_200cyc.csv';
    'Ch15_RPT_0cyc.csv';   'Ch15_RPT_200cyc.csv';
    'Ch16_RPT_0cyc.csv';   'Ch16_RPT_200cyc.csv';
};

channelList = unique(extractBetween(fileNames, 'Ch', '_'));
cycleList = {'0cyc', '200cyc'};

sc_stepIdx_chg = 1;
sc_stepIdx_dch = 3;

stepIdx_chg = 8;  % OCV-Charge
stepIdx_dch = 10; % OCV-Discharge
soc_grid = 0:1:100;  % 1% 간격

for i = 1:length(channelList)
    ch = channelList{i};

    % --- Figure 생성 또는 검색: SOC–OCV
    fig_soc = findobj('Type', 'figure', 'Name', ['SOC–OCV - Channel ' ch]);
    if isempty(fig_soc)
        fig_soc = figure('Name', ['SOC–OCV - Channel ' ch], 'Color', 'w');
    else
        figure(fig_soc);
    end
    hold on;

    % --- Figure 생성 또는 검색: SOC–OCV
    Q_voltage = findobj('Type', 'figure', 'Name', ['Q–OCV - Channel ' ch]);
    if isempty(Q_voltage)
        Q_voltage = figure('Name', ['Q–OCV - Channel ' ch], 'Color', 'w');
    else
        figure(Q_voltage);
    end
    hold on;

    % --- Figure 생성 또는 검색: Capacity
    fig_cap = findobj('Type', 'figure', 'Name', ['Capacity - Channel ' ch]);
    if isempty(fig_cap)
        fig_cap = figure('Name', ['Capacity - Channel ' ch], 'Color', 'w');
    else
        figure(fig_cap);
    end
    hold on;

    for j = 1:length(cycleList)
        cyc = cycleList{j};
        filename = sprintf('Ch%s_RPT_%s.csv', ch, cyc);
        filepath = fullfile(csvPath, filename);
        if ~isfile(filepath), fprintf('파일 없음: %s\n', filename); continue; end

        T = readtable(filepath);
        stepIndex   = T{:,2};
        cycleIndex  = T{:,4};
        current     = T{:,7};
        voltage     = T{:,8};
        capacity    = T{:,9};

%% SOC-COV
        % 충전 SOC-OCV
        mask_chg = (stepIndex == 8) & (cycleIndex == 2);
        V_chg = voltage(mask_chg);
        Q_chg = capacity(mask_chg);
        Q_nom_chg = Q_chg(end);
        SOC_chg = 100 * Q_chg / Q_nom_chg;

        % 방전 SOC-OCV
        mask_dch = (stepIndex == 10) & (cycleIndex == 2);
        V_dch = voltage(mask_dch);
        Q_dch = capacity(mask_dch);
        Q_nom_dch = Q_dch(end);
        SOC_dch = 100 * (1 - Q_dch / Q_nom_dch);
        [SOC_dch_sorted, idx_d] = sort(SOC_dch, 'ascend');
        V_dch_sorted = V_dch(idx_d);
        Q_dch_sorted = Q_dch(idx_d);  % Q_dch도 같은 순서로 정렬

        % 1% 간격 보간 후 평균
        V_chg_interp = interp1(SOC_chg, V_chg, soc_grid, 'linear');
        V_dch_interp = interp1(SOC_dch_sorted, V_dch_sorted, soc_grid, 'linear');
        V_avg = (V_chg_interp + V_dch_interp) / 2;

        % Interpolation function 저장 (기존 SOC 기준)
        OCV_func = @(soc_query) interp1(soc_grid, V_avg, soc_query, 'linear');
        fieldname = sprintf('ch%s_%s', ch, cyc);
        OCV_struct.(fieldname).OCV_func = OCV_func;
        OCV_struct.(fieldname).soc_grid = soc_grid;
        OCV_struct.(fieldname).V_OCV = V_avg;
        
        % ECM 피팅용 용량 정보와 Q-OCV 데이터 추가
        OCV_struct.(fieldname).capacity_chg = Q_nom_chg;  % 충전 용량 [Ah]
        OCV_struct.(fieldname).capacity_dch = Q_nom_dch;  % 방전 용량 [Ah]
        OCV_struct.(fieldname).Q_chg = Q_chg;             % 충전 Q 데이터 [Ah]
        OCV_struct.(fieldname).V_chg = V_chg;             % 충전 전압 데이터 [V]
        OCV_struct.(fieldname).Q_dch = Q_dch_sorted;      % 방전 Q 데이터 [Ah] (정렬됨)
        OCV_struct.(fieldname).V_dch = V_dch_sorted;      % 방전 전압 데이터 [V] (정렬됨)
        OCV_struct.(fieldname).SOC_chg = SOC_chg;         % 충전 SOC 데이터 [%]
        OCV_struct.(fieldname).SOC_dch = SOC_dch_sorted;  % 방전 SOC 데이터 [%] (정렬됨)
        
        % Q-OCV 평균 함수 생성 (방전 기준)
        Q_OCV_func = @(q_query) interp1(Q_dch_sorted, V_dch_sorted, q_query, 'linear');
        OCV_struct.(fieldname).Q_OCV_func = Q_OCV_func;
        
        fprintf('  - 처리 완료: %s (충전: %.2f Ah, 방전: %.2f Ah)\n', ...
                fieldname, Q_nom_chg, Q_nom_dch);

        % SOC–OCV 플롯
        figure(fig_soc); hold on;
        if strcmp(cyc, '0cyc')
            plot(soc_grid, V_avg, 'b-o', 'LineWidth', 1.5,'MarkerSize', 2,'DisplayName', 'BOL');
        else
            plot(soc_grid, V_avg, 'r-o', 'LineWidth', 1.5,'MarkerSize', 2,'DisplayName', 'RPT2');
        end
        xlabel('SOC [%]'); ylabel('Voltage [V]');
        title(sprintf('Channel %s - SOC vs OCV Comparison', ch));
        legend('Location', 'best'); grid on;
        savefig(fig_soc, fullfile(saveFolder, sprintf('Ch%s_SOC_OCV_.fig', ch)));

        % Q-OCV 플롯 (방전만)
        figure(Q_voltage); hold on;
        if strcmp(cyc, '0cyc')
            % plot(Q_dch, V_dch_sorted, 'b-o', 'LineWidth', 1.5,'MarkerSize', 2,'DisplayName', 'BOL');
            plot(Q_dch, V_dch_sorted, 'b-o', 'LineWidth', 1.5,'MarkerSize', 2,'DisplayName', sprintf('Discharge SOC (BOL) (Q=%.2fAh)', Q_nom_dch));

        else
            % plot(Q_dch, V_dch_sorted, 'r-o', 'LineWidth', 1.5,'MarkerSize', 2,'DisplayName', 'RPT2');
            plot(Q_dch, V_dch_sorted, 'r-o', 'LineWidth', 1.5,'MarkerSize', 2,'DisplayName', sprintf('Discharge SOC (RPT2) (Q=%.2fAh)', Q_nom_dch));
        end
        xlabel('Capacity [Ah]'); ylabel('Voltage [V]');
        title(sprintf('Channel %s - Q vs OCV Comparison', ch));
        legend('Location', 'best'); grid on;
        savefig(Q_voltage, fullfile(saveFolder, sprintf('Ch%s_Q_OCV_.fig', ch)));

%% CAPACITY
        % Static capacity용 step
        static_stepIdx_chg = 1;  
        static_stepIdx_dch = 3;  

        mask_cap_chg = (stepIndex == static_stepIdx_chg) & (cycleIndex == 2);
        mask_cap_dch = (stepIndex == static_stepIdx_dch) & (cycleIndex == 2);

        V_chg_cap = voltage(mask_cap_chg);
        Q_chg_cap = capacity(mask_cap_chg);
        V_dch_cap = voltage(mask_cap_dch);
        Q_dch_cap = capacity(mask_cap_dch);

        % Plot
        figure(fig_cap); hold on;
        Q_chg_end = Q_chg_cap(end);
        Q_dch_end = Q_dch_cap(end);

        if strcmp(cyc, '0cyc')
            plot(Q_chg_cap, V_chg_cap, 'b-', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Charge BOL (Q=%.2fAh)', Q_chg_end));
            plot(Q_dch_cap, V_dch_cap, 'r-', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Discharge BOL (Q=%.2fAh)', Q_dch_end));
        else
            plot(Q_chg_cap, V_chg_cap, 'b--', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Charge RPT2 (Q=%.2fAh)', Q_chg_end));
            plot(Q_dch_cap, V_dch_cap, 'r--', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Discharge RPT2 (Q=%.2fAh)', Q_dch_end));
        end
        xlabel('Capacity [Ah]'); ylabel('Voltage [V]');
        title(['Channel ' ch ' - Static Capacity (0.5C) vs Voltage']);
        legend('Location', 'best'); grid on;
        savefig(fig_cap, fullfile(saveFolder, sprintf('Ch%s_Capacity.fig', ch)));


    end
end

save(fullfile(saveFolder, 'OCV_interpolation_functions.mat'), 'OCV_struct');
