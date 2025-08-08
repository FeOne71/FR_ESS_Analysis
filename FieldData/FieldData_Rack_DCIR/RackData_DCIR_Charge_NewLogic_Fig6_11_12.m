%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rack Data DCIR Charge - NewLogic with Fig6_11_12 Logic (OPTIMIZED)
% Charging impedance analysis with NewLogic peak detection
% Based on fig_6_11_12.m logic for charging event extraction and analysis
% PERFORMANCE OPTIMIZED VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory (NewLogic4_5와 동일한 dataDir 사용)
dataDir = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Rack_raw2mat\New\2025';
yearList = {'2025'};
saveDir = fullfile('G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_Rack_DCIR\DCIR_Charge_NewLogic_Fig6_11_12');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters (fig_6_11_12와 동일한 파라미터)
C_nom_cell = 128;
thr = C_nom_cell * 0.01;   % Initial current threshold (A)
dI  = C_nom_cell * 0.2;    % Current change threshold (A)
ddI = 1;                   % Continuous current increase threshold (A)

%% Initialize optimized data structures
% 사전 할당으로 메모리 효율성 향상
max_events = 10000;  % 예상 최대 이벤트 수
max_points_per_event = 100000;  % 예상 최대 포인트 수

Time_n = cell(max_events, 1);
Current_n = cell(max_events, 1);
Voltage_n = cell(max_events, 1);
SoC_n = cell(max_events, 1);
Temp_n = cell(max_events, 1);
DatesVector = cell(max_events, 1);

%% Data Collection and Event Extraction (OPTIMIZED)
event_counter = 1;

% Process monthly folders
monthDirs = dir(fullfile(dataDir, '20*'));
fprintf('Found %d month directories\n', length(monthDirs));

for m = 1:length(monthDirs)
    if ~monthDirs(m).isdir, continue; end
    monthPath = fullfile(dataDir, monthDirs(m).name);
    matFiles = dir(fullfile(monthPath, '*.mat'));

    month_num = str2num(monthDirs(m).name(5:6));
    fprintf('Processing month: %d (%d files)\n', month_num, length(matFiles));
    
    for f = 1:length(matFiles)
        if mod(f, 10) == 0  % 진행상황 표시 최적화
            fprintf('  File %d/%d\n', f, length(matFiles));
        end
        
        matFilePath = fullfile(monthPath, matFiles(f).name);
        load(matFilePath);

        % Extract data from Raw structure
        t = Raw.Date_Time_seconds;
        I = Raw.DCCurrent;
        V = Raw.CVavg;
        T_batt = Raw.MTavg;
        soc = Raw.SOC_BMS;

        % Ensure column vectors
        if isrow(t), t = t'; end
        if isrow(I), I = I'; end
        if isrow(V), V = V'; end
        if isrow(T_batt), T_batt = T_batt'; end
        if isrow(soc), soc = soc'; end

        % OPTIMIZED: 벡터화된 시간 차이 계산
        time_diff = diff(t);
        gap_indices = find(time_diff >= 120);  % 2분 이상 간격
        
        % 이벤트 시작/끝 인덱스 계산
        if isempty(gap_indices)
            % 단일 이벤트
            event_starts = 1;
            event_ends = length(t);
        else
            % 여러 이벤트
            event_starts = [1; gap_indices + 1];
            event_ends = [gap_indices; length(t)];
        end
        
        % 각 이벤트 처리
        for evt_idx = 1:length(event_starts)
            start_idx = event_starts(evt_idx);
            end_idx = event_ends(evt_idx);
            
            if end_idx - start_idx > 10  % 최소 길이 체크
                % 이벤트 데이터 추출 (벡터화)
                Time_n{event_counter} = t(start_idx:end_idx);
                Current_n{event_counter} = I(start_idx:end_idx);
                Voltage_n{event_counter} = V(start_idx:end_idx);
                SoC_n{event_counter} = soc(start_idx:end_idx);
                Temp_n{event_counter} = T_batt(start_idx:end_idx);
                DatesVector{event_counter} = datetime(t(start_idx), 'convertfrom', 'posixtime');
                
                event_counter = event_counter + 1;
            end
        end
        
        clear Raw
    end
end

% 실제 사용된 셀만 유지
Time_n = Time_n(1:event_counter-1);
Current_n = Current_n(1:event_counter-1);
Voltage_n = Voltage_n(1:event_counter-1);
SoC_n = SoC_n(1:event_counter-1);
Temp_n = Temp_n(1:event_counter-1);
DatesVector = DatesVector(1:event_counter-1);

fprintf('Total events extracted: %d\n', length(Time_n));

%% Event Filtering (OPTIMIZED)
fprintf('Filtering charging events...\n');

% Parameters
initial = 180;  % 초기 18000개 포인트 제거
final = 180;    % 마지막 18000개 포인트 제거
interval = (initial + final) / 100;  % 360

% OPTIMIZED: 사전 할당
valid_events = false(length(Time_n), 1);

% 벡터화된 필터링 조건 체크
for i = 1:length(Time_n)
    if ~isempty(Time_n{i})
        event_time = Time_n{i}(end) - Time_n{i}(1);
        event_length = length(Time_n{i});
        mean_current = mean(Current_n{i});
        
        % 모든 조건을 한번에 체크
        valid_events(i) = (event_time > 2 * interval) && ...
                          (mean_current > 0) && ...
                          (event_length > initial + final);
    end
end

% 필터링된 이벤트 추출
Time_mod_clean = Time_n(valid_events);
Current_clean = Current_n(valid_events);
Voltage_clean = Voltage_n(valid_events);
SoC_clean = SoC_n(valid_events);
Temp_clean = Temp_n(valid_events);
DatesVector_clean = DatesVector(valid_events);

fprintf('Total filtered events: %d\n', length(Time_mod_clean));

%% Compute charging impedance (OPTIMIZED)
fprintf('Computing charging impedance...\n');

dt = [1, 10000];  % Time window: [0.1s, 100s]
R_pseudo = cell(length(dt), 1);

% OPTIMIZED: 벡터화된 DCIR 계산
for k_dt = 1:length(dt)
    R_pseudo{k_dt} = cell(length(Time_mod_clean), 1);
    
    for i = 1:length(Time_mod_clean)
        if length(Time_mod_clean{i}) > dt(k_dt)
            % 벡터화된 계산
            V1 = Voltage_clean{i}(1:end-dt(k_dt));
            V2 = Voltage_clean{i}(1+dt(k_dt):end);
            I1 = Current_clean{i}(1:end-dt(k_dt));
            
            % DCIR 계산 (벡터화)
            dV = V2 - V1;
            dI = I1;  % Current는 동일한 값 사용
            R_pseudo{k_dt}{i} = abs(dV ./ dI) * 1000;
        end
    end
end

%% Average current and charging levels (OPTIMIZED)
fprintf('Analyzing charging levels...\n');

% 벡터화된 평균 계산
I_ave = cellfun(@mean, Current_clean);

% Battery pack capacity [Ah]
Cap = 240;

% Sampling time [s]
Ts = 0.01;

%% Division depending on charging levels (OPTIMIZED)
% 벡터화된 분류
low_mask = I_ave < 5;
medium_mask = I_ave >= 5 & I_ave < 80;
high_mask = I_ave >= 80;

% 각 레벨별 데이터 추출
Time_mod_1 = Time_mod_clean(low_mask);
Curr_mod_1 = Current_clean(low_mask);
SoC_mod_1 = SoC_clean(low_mask);
Temp_mod_1 = Temp_clean(low_mask);
DatesVector_mod_1 = DatesVector_clean(low_mask);
R_pseudo_1 = {R_pseudo{1}(low_mask), R_pseudo{2}(low_mask)};
I_ave_1 = I_ave(low_mask);
lev1 = find(low_mask);

Time_mod_2 = Time_mod_clean(medium_mask);
Curr_mod_2 = Current_clean(medium_mask);
SoC_mod_2 = SoC_clean(medium_mask);
Temp_mod_2 = Temp_clean(medium_mask);
DatesVector_mod_2 = DatesVector_clean(medium_mask);
R_pseudo_2 = {R_pseudo{1}(medium_mask), R_pseudo{2}(medium_mask)};
I_ave_2 = I_ave(medium_mask);
lev2 = find(medium_mask);

Time_mod_3 = Time_mod_clean(high_mask);
Curr_mod_3 = Current_clean(high_mask);
SoC_mod_3 = SoC_clean(high_mask);
Temp_mod_3 = Temp_clean(high_mask);
DatesVector_mod_3 = DatesVector_clean(high_mask);
R_pseudo_3 = {R_pseudo{1}(high_mask), R_pseudo{2}(high_mask)};
I_ave_3 = I_ave(high_mask);
lev3 = find(high_mask);

fprintf('Low charging events (<5A): %d\n', length(lev1));
fprintf('Medium charging events (5-80A): %d\n', length(lev2));
fprintf('High charging events (>=80A): %d\n', length(lev3));

%% Figure 6a - C-rate vs Time (OPTIMIZED)
figure(1);
sz = 300; Line = 0.5;
coloredge = 'k';

% Low charging events
if ~isempty(I_ave_1)
    scatter([DatesVector_mod_1{:}], I_ave_1/Cap, sz, 'MarkerEdgeColor', coloredge, ...
        'MarkerFaceColor', [0.929 0.6940 0.1250], ...
        'LineWidth', Line, 'DisplayName', 'Low (<5A)');
    hold on;
end

% Medium charging events
if ~isempty(I_ave_2)
    scatter([DatesVector_mod_2{:}], I_ave_2/Cap, sz, 'MarkerEdgeColor', coloredge, ...
        'MarkerFaceColor', [0 0.4470 0.7410], ...
        'LineWidth', Line, 'DisplayName', 'Medium (5-80A)');
end

% High charging events
if ~isempty(I_ave_3)
    scatter([DatesVector_mod_3{:}], I_ave_3/Cap, sz, 'MarkerEdgeColor', coloredge, ...
        'MarkerFaceColor', [0.6350 0.0780 0.1840], ...
        'LineWidth', Line, 'DisplayName', 'High (>=80A)');
end

ylabel('C-rate [1/h]', 'interpreter', 'tex');
yticks([0 1/10 1/5 1/3.3 1/2.5 1/2]);
yticklabels({'' 'C/10' 'C/5' 'C/3.3' 'C/2.5' 'C/2'});
ylim([0 0.5]);
xtickangle(45);
set(gca, 'fontsize', 20);
title('Charging Events by C-rate (NewLogic - Optimized)');
box on;

%% Figure 6b - Event Count by Charging Level (OPTIMIZED)
figure(2);
x1 = 1; Bin1 = length(lev1);
x2 = 2; Bin2 = length(lev2);
x3 = 3; Bin3 = length(lev3);

bars = [length(lev1) length(lev2) length(lev3)];
bar(x1, Bin1, 'FaceColor', [0.9290 0.6940 0.1250], 'EdgeColor', [0 0 0], 'LineWidth', 0.5);
hold all;
bar(x2, Bin2, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', [0 0 0], 'LineWidth', 0.5);
bar(x3, Bin3, 'FaceColor', [0.6350 0.0780 0.1840], 'EdgeColor', [0 0 0], 'LineWidth', 0.5);

ylabel('Number of charging events [#]', 'interpreter', 'tex');
xticklabels({'Low (<5A)', 'Medium (5-80A)', 'High (>=80A)'});
xticks(1:3);
xlim([0 4]);
set(gca, 'ticklabelinterpreter', 'tex');
set(gca, 'fontsize', 20);
title('Charging Events Distribution (NewLogic - Optimized)');

%% Figure 11 - Low Charging SoC vs Impedance (OPTIMIZED)
if ~isempty(lev1)
    figure(3);
    step = 100;  % Plotting step
    remove_transient = 2000/Ts; sz = 20;
    
    for k = 1:length(dt)
        f = figure(3 + k);
        box on; hold all;
        
        for i = 1:length(R_pseudo_1{1,k})
            if ~isempty(R_pseudo_1{1,k}{i})
                scatter(SoC_mod_1{i}(remove_transient:step:end-dt(k)), ...
                    R_pseudo_1{1,k}{i}(remove_transient:step:end), ...
                    sz, Temp_mod_1{i}(remove_transient:step:end-dt(k)), 'filled', 'LineWidth', 2);
            end
        end
        
        xlabel('SoC [%]', 'interpreter', 'tex');
        ylabel('Z_{CHG} [m\Omega]', 'interpreter', 'tex');
        xlim([0 100]);
        ylim([0 0.40*dt(k)/100]);
        colormap(flipud(autumn));
        c = colorbar('eastoutside');
        caxis([10 32.5]);
        c.Label.String = 'Temperature [°C]';
        c.Label.Interpreter = 'tex';
        set(findall(gcf, '-property', 'interpreter'), 'interpreter', 'tex');
        set(findall(gcf, '-property', 'ticklabelinterpreter'), 'ticklabelinterpreter', 'tex');
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
        title(sprintf('Low Charging Impedance (dt=%d) - NewLogic Optimized', dt(k)));
    end
end

%% Figure 11 - Medium Charging SoC vs Impedance (OPTIMIZED)
if ~isempty(lev2)
    figure(6);
    remove_transient = 2000/Ts;
    
    for k = 1:length(dt)
        f = figure(6 + k);
        box on; hold all;
        
        for i = 1:length(R_pseudo_2{1,k})
            if ~isempty(R_pseudo_2{1,k}{i})
                scatter(SoC_mod_2{i}(remove_transient:step:end-dt(k)), ...
                    R_pseudo_2{1,k}{i}(remove_transient:step:end), ...
                    sz, Temp_mod_2{i}(remove_transient:step:end-dt(k)), 'filled', 'LineWidth', 2);
            end
        end
        
        xlabel('SoC [%]', 'interpreter', 'tex');
        ylabel('Z_{CHG} [m\Omega]', 'interpreter', 'tex');
        xlim([0 100]);
        ylim([0 0.25*dt(k)/100]);
        colormap(flipud(autumn));
        c = colorbar('eastoutside');
        caxis([10 32.5]);
        c.Label.String = 'Temperature [°C]';
        c.Label.Interpreter = 'tex';
        set(findall(gcf, '-property', 'interpreter'), 'interpreter', 'tex');
        set(findall(gcf, '-property', 'ticklabelinterpreter'), 'ticklabelinterpreter', 'tex');
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
        title(sprintf('Medium Charging Impedance (dt=%d) - NewLogic Optimized', dt(k)));
    end
end

%% Figure 11 - High Charging SoC vs Impedance (OPTIMIZED)
if ~isempty(lev3)
    figure(9);
    remove_transient = 50/Ts;
    
    for k = 1:length(dt)
        f = figure(9 + k);
        box on; hold all;
        
        for i = 1:length(R_pseudo_3{1,k})
            if ~isempty(R_pseudo_3{1,k}{i})
                scatter(SoC_mod_3{i}(remove_transient:step:end-dt(k)), ...
                    R_pseudo_3{1,k}{i}(remove_transient:step:end), ...
                    sz, Temp_mod_3{i}(remove_transient:step:end-dt(k)), 'filled', 'LineWidth', 2);
            end
        end
        
        xlabel('SoC [%]', 'interpreter', 'tex');
        ylabel('Z_{CHG} [m\Omega]', 'interpreter', 'tex');
        xlim([0 100]);
        ylim([0 0.25*dt(k)/100]);
        colormap(flipud(autumn));
        c = colorbar('eastoutside');
        caxis([10 32.5]);
        c.Label.String = 'Temperature [°C]';
        c.Label.Interpreter = 'tex';
        set(findall(gcf, '-property', 'interpreter'), 'interpreter', 'tex');
        set(findall(gcf, '-property', 'ticklabelinterpreter'), 'ticklabelinterpreter', 'tex');
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
        title(sprintf('High Charging Impedance (dt=%d) - NewLogic Optimized', dt(k)));
    end
end

%% Save results
save(fullfile(saveDir, 'Charging_Events_NewLogic_Fig6_11_12_Optimized.mat'), ...
    'Time_mod_clean', 'Current_clean', 'Voltage_clean', 'SoC_clean', 'Temp_clean', ...
    'R_pseudo', 'I_ave', 'lev1', 'lev2', 'lev3', ...
    'Time_mod_1', 'Curr_mod_1', 'SoC_mod_1', 'Temp_mod_1', 'R_pseudo_1', 'I_ave_1', ...
    'Time_mod_2', 'Curr_mod_2', 'SoC_mod_2', 'Temp_mod_2', 'R_pseudo_2', 'I_ave_2', ...
    'Time_mod_3', 'Curr_mod_3', 'SoC_mod_3', 'Temp_mod_3', 'R_pseudo_3', 'I_ave_3');

fprintf('Processing complete (OPTIMIZED VERSION)\n');
fprintf('Results saved to: %s\n', saveDir); 