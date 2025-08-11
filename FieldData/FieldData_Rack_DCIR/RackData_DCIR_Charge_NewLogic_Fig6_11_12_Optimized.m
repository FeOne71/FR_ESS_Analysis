%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rack Data DCIR Charge - NewLogic with Fig6_11_12 Logic (OPTIMIZED)
% Charging impedance analysis with NewLogic peak detection
% Based on fig_6_11_12.m logic for charging event extraction and analysis
% PERFORMANCE OPTIMIZED VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory (NewLogic4_5와 동일한 dataDir 사용)
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New\2025';
yearList = {'2025'};
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\DCIR_Charge_NewLogic_Fig6_11_12');
%% 
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters (CurrentClustering_Auto와 동일)
% Battery pack capacity [Ah]
C_nom_cell = 128;           % Ah
idle_thr = C_nom_cell * 0.02;   % Initial current threshold (A)
P_nom_rack = 112;           % 3.68*64/1000 [kW] | 0.23552 * 2P * 14S * 17S = 112
min_duration = 10;         % [s] - Charging duration (5 minutes or more)
max_P_std = P_nom_rack * 0.1;  % Max power standard deviation [kW] 1.12kW
max_I_std = C_nom_cell * 0.01;   % Max current standard deviation [A] 1.28A

% Sampling time [s]
Ts = 1.0;  % 1초 간격 데이터

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

% Process monthly folders (한 달치만 처리)
monthDirs = dir(fullfile(dataDir, '20*'));
fprintf('Found %d month directories\n', length(monthDirs));

%for m = 1:length(monthDirs)
for m = 1:1  % 1개 월만 처리
    if ~monthDirs(m).isdir, continue; end
    monthPath = fullfile(dataDir, monthDirs(m).name);
    matFiles = dir(fullfile(monthPath, '*.mat'));

    month_num = str2num(monthDirs(m).name(5:6));
    fprintf('Processing month: %d (%d files)\n', month_num, length(matFiles));
    
    for f = 1:length(matFiles)
        %if mod(f, 10) == 0  % 진행상황 표시 최적화
        %    fprintf('  File %d/%d\n', f, length(matFiles));
        %end
        
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

        % 충전 구간만 추출 (idle 제거)
        % 양수 전류 구간을 충전으로 정의
        charging_mask = I > 0;
        
        % 연속된 충전 구간 찾기
        charging_starts = find(diff([0; charging_mask; 0]) == 1);
        charging_ends = find(diff([0; charging_mask; 0]) == -1) - 1;
        
        % 각 충전 구간을 이벤트로 처리 (길이 제한 없음)
        for evt_idx = 1:length(charging_starts)
            start_idx = charging_starts(evt_idx);
            end_idx = charging_ends(evt_idx);
            
        % 모든 충전 구간을 이벤트로 저장 (길이 제한 없음)
                Time_n{event_counter} = t(start_idx:end_idx);
                Current_n{event_counter} = I(start_idx:end_idx);
                Voltage_n{event_counter} = V(start_idx:end_idx);
                SoC_n{event_counter} = soc(start_idx:end_idx);
                Temp_n{event_counter} = T_batt(start_idx:end_idx);
        DatesVector{event_counter} = t(start_idx);  % 이벤트 시작 시간 저장
                event_counter = event_counter + 1;
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

%% Event Filtering (Clustering_Auto 로직 적용)
fprintf('Filtering charging events using Clustering_Auto logic...\n');

% OPTIMIZED: 사전 할당
valid_events = false(length(Time_n), 1);

% 충전 구간 필터링 로직 (설정한 파라미터 사용)
filtered_by_duration = 0;
filtered_by_stability = 0;
filtered_by_current = 0;

for i = 1:length(Time_n)
    if ~isempty(Time_n{i})
        % 충전 duration 체크 (60초 이상)
        charge_duration = length(Time_n{i});
        if charge_duration < 60  % 60초 미만이면 제외
            filtered_by_duration = filtered_by_duration + 1;
            continue;
        end
        
        % 전류/전력 안정성 체크 (3번째 포인트부터)
        if charge_duration > 3
            current_std = std(Current_n{i}(3:end));
            % 전력 계산 (P = V * I)
            power = Voltage_n{i} .* Current_n{i};
            power_std = std(power(3:end));
            
            % 안정성 조건 체크 (설정한 파라미터 사용)
            if current_std >= max_I_std || power_std >= max_P_std
                filtered_by_stability = filtered_by_stability + 1;
                continue;
            end
        end
        

        
        % 평균 전류 체크 (충전 조건)
        mean_current = mean(Current_n{i});
        c_rate = mean_current / C_nom_cell;  % C-rate 계산
        
        if mean_current > 0 && c_rate >= 0.05  % 양수 전류 (충전) + 0.05C 이상
            valid_events(i) = true;
        else
            filtered_by_current = filtered_by_current + 1;
        end
    end
end

% 필터링 결과 출력
fprintf('Filtering results:\n');
fprintf('  - Filtered by duration (<60 s): %d\n', filtered_by_duration);
fprintf('  - Filtered by stability: %d\n', filtered_by_stability);
fprintf('  - Filtered by current (<=0 or <0.05C): %d\n', filtered_by_current);
fprintf('  - Valid events: %d\n', sum(valid_events));

% 필터링된 이벤트 추출
Time_clean = Time_n(valid_events);
Current_clean = Current_n(valid_events);
Voltage_clean = Voltage_n(valid_events);
SoC_clean = SoC_n(valid_events);
Temp_clean = Temp_n(valid_events);
DatesVector_clean = DatesVector(valid_events);

% Transient 제거 없이 바로 사용
fprintf('Events after filtering: %d\n', length(Time_clean));

fprintf('Total filtered events: %d\n', length(Time_clean));

%% Compute charging impedance (CurrentClustering_Auto와 동일)
fprintf('Computing charging impedance...\n');

% 레퍼런스와 동일한 시간 윈도우 사용
dt = [30, 120];  % [30초, 120초] (더 큰 변화 감지)
R_pseudo = cell(length(dt), length(Time_clean));

% 레퍼런스와 동일한 DCIR 계산
for k_dt = 1:length(dt)
    for i = 1:length(Time_clean)
        if length(Time_clean{i}) > dt(k_dt)
            R_pseudo{k_dt}{i} = [];
            
            % 레퍼런스와 동일한 계산 방식
            for j = 1:(length(Current_clean{i})-dt(k_dt))
                V1 = Voltage_clean{i}(j);
                V2 = Voltage_clean{i}(j+dt(k_dt));
                I1 = Current_clean{i}(j);
                
                % 레퍼런스 방식: (V2-V1) / I1 * (1000)
                if abs(V2 - V1) > 0 && abs(I1) > 0  % 최소 변화량 조건
                    R_pseudo{k_dt}{i}(j) = abs((V2 - V1)) / I1 * (1000);  % mΩ 단위
                else
                    R_pseudo{k_dt}{i}(j) = NaN;  % 조건 불만족시 NaN
                end
            end
            
            % 디버깅: 값 확인
            valid_values = R_pseudo{k_dt}{i}(~isnan(R_pseudo{k_dt}{i}));
            if ~isempty(valid_values)
                fprintf('Event %d, dt=%d DCIR: min=%.6f, max=%.6f, mean=%.6f\n', ...
                    i, dt(k_dt), min(valid_values), max(valid_values), mean(valid_values));
            end
        end
    end
end

%% Average current and charging levels (OPTIMIZED)
fprintf('Analyzing charging levels...\n');

% 벡터화된 평균 계산
I_ave = cellfun(@mean, Current_clean);

% Sampling time [s]
Ts = 1.0;  % 1초 간격 데이터

%% Division depending on charging levels (OPTIMIZED)
% 실제 데이터 분포에 맞는 전류 분류
fprintf('Current range: %.2f ~ %.2f A\n', min(I_ave), max(I_ave));

% 실제 데이터 분포에 따른 분류 (C-rate 기준)
low_mask = I_ave < C_nom_cell * 0.1;           % 0.1C 미만
medium_mask = I_ave >= C_nom_cell * 0.1 & I_ave < C_nom_cell * 0.5;  % 0.1C-0.5C
high_mask = I_ave >= C_nom_cell * 0.5;         % 0.5C 이상

fprintf('Low charging events (<%.1fA): %d\n', C_nom_cell * 0.1, sum(low_mask));
fprintf('Medium charging events (%.1f-%.1fA): %d\n', C_nom_cell * 0.1, C_nom_cell * 0.5, sum(medium_mask));
fprintf('High charging events (>=%.1fA): %d\n', C_nom_cell * 0.5, sum(high_mask));

% 가장 많은 이벤트가 있는 구간 선택
event_counts = [sum(low_mask), sum(medium_mask), sum(high_mask)];
[~, max_idx] = max(event_counts);

if max_idx == 1
    selected_mask = low_mask;
    selected_name = 'Low';
elseif max_idx == 2
    selected_mask = medium_mask;
    selected_name = 'Medium';
else
    selected_mask = high_mask;
    selected_name = 'High';
end

fprintf('Selected category: %s (%d events)\n', selected_name, sum(selected_mask));

%% Create Event Structure (Original Format)
fprintf('Creating event structure...\n');

% Initialize global event structure
global_eventStruct = struct();
eventStruct = struct();  % eventStruct 초기화 추가

% Create event structure for each charging event
eventCount = 0;
for i = 1:length(Time_clean)
    eventCount = eventCount + 1;
    evtName = sprintf('event%d', eventCount);
    
    % Extract date from Raw data
    if ~isempty(Time_clean{i})
        % Raw 데이터에서 직접 날짜 추출
        event_date = datetime(Time_clean{i}(1), 'convertfrom', 'posixtime');
        formatted_date = datestr(event_date, 'yyyy-mm-dd');
    else
        formatted_date = '2025-01-01';
    end
    
    % Calculate C-rate
    avg_current = I_ave(i);
    c_rate = avg_current / C_nom_cell;
    
    % Determine C-rate category
    if avg_current < C_nom_cell * 0.1
        c_rate_category = 'Low';
    elseif avg_current < C_nom_cell * 0.5 && avg_current >= C_nom_cell * 0.1
        c_rate_category = 'Medium';
    elseif avg_current >= C_nom_cell * 0.5
        c_rate_category = 'High';
    end
    
    % Create event structure
    eventStruct.(evtName).rack_name = 'Rack01';  % Default rack name
    eventStruct.(evtName).year = '2025';
    eventStruct.(evtName).date = formatted_date;
    eventStruct.(evtName).start_idx = 1;
    eventStruct.(evtName).end_idx = length(Time_clean{i});
    eventStruct.(evtName).charge_duration = length(Time_clean{i});
    eventStruct.(evtName).avg_current = avg_current;
    eventStruct.(evtName).c_rate = c_rate;
    eventStruct.(evtName).c_rate_category = c_rate_category;
    
    % Time series data
    eventStruct.(evtName).t_seq = Time_clean{i};
    eventStruct.(evtName).I_seq = Current_clean{i};
    eventStruct.(evtName).V_seq = Voltage_clean{i};
    eventStruct.(evtName).T_seq = Temp_clean{i};
    eventStruct.(evtName).soc_seq = SoC_clean{i};
    
    % DCIR values (single time window)
    if ~isempty(R_pseudo{i})
        eventStruct.(evtName).DCIR_dt10 = R_pseudo{i};
    end
    
    % Store in global structure by C-rate category
    if ~isfield(global_eventStruct, c_rate_category)
        global_eventStruct.(c_rate_category) = struct();
    end
    global_eventStruct.(c_rate_category).(evtName) = eventStruct.(evtName);
end

%% Generate Figures (Fig 6c, Fig 12만 표시)
fprintf('Generating figures...\n');

% 디버깅 정보 출력
fprintf('Low events: %d, Medium events: %d, High events: %d\n', ...
    sum(low_mask), sum(medium_mask), sum(high_mask));

% Parameters for plotting
sz = 5;  % 점 크기 축소

% 선택된 카테고리로 데이터 분류
% 모든 변수를 미리 초기화
Time_mod_selected = {}; Current_mod_selected = {}; SoC_mod_selected = {}; Temp_mod_selected = {};
R_pseudo_selected = {};

if sum(selected_mask) > 0
    Time_selected = Time_clean(selected_mask);
    Current_selected = Current_clean(selected_mask);
    Voltage_selected = Voltage_clean(selected_mask);
    SoC_selected = SoC_clean(selected_mask);
    Temp_selected = Temp_clean(selected_mask);
    R_pseudo_selected = R_pseudo{1}(selected_mask);  % dt(1) 사용
end
        
%% Figure 6c - Selected charging events (Z_CHG vs SoC)
% Z_CHG = R_pseudo (충전 임피던스)

if sum(selected_mask) > 0
    figure; box on; hold all
    fprintf('Fig 6c - %s events: %d\n', selected_name, length(R_pseudo_selected));
    fprintf('Debug: dt = [%d, %d], Ts = %.1f\n', dt(1), dt(2), Ts);
    
    plotted_events = 0;
    for i = 1:length(R_pseudo_selected)
        fprintf('Event %d: SoC length = %d, R_pseudo length = %d\n', ...
            i, length(SoC_selected{i}), length(R_pseudo_selected{i}));
        
        if length(SoC_selected{i}) > dt(1) && ...
           ~isempty(R_pseudo_selected{i}) && ...
           length(R_pseudo_selected{i}) > 0
            
            % R_pseudo와 해당 시점의 SoC 매칭 (dt(1) 사용)
            r_length = length(R_pseudo_selected{i});
            soc_length = length(SoC_selected{i});
            
            % 길이 확인 및 조정
            if soc_length > dt(1) && r_length > 0
                soc_indices = (1:r_length) + dt(1);  % dt(1)만큼 뒤로 이동
                r_indices = 1:r_length;
                temp_indices = (1:r_length) + dt(1);  % 온도도 동일하게 이동
                
                % 인덱스 범위 체크
                if max(soc_indices) <= soc_length
                    % 디버깅: SoC와 R_pseudo 값 범위 출력
                    soc_values = SoC_selected{i}(soc_indices);
                    r_values = R_pseudo_selected{i}(r_indices);
                    fprintf('  Event %d: SoC range [%.2f, %.2f], R_pseudo range [%.2f, %.2f]\n', ...
                        i, min(soc_values), max(soc_values), min(r_values), max(r_values));
                    
                    % 모든 데이터 포인트를 개별적으로 플롯
                    for k = 1:length(soc_indices)
                        scatter(SoC_selected{i}(soc_indices(k)), ...
                               R_pseudo_selected{i}(r_indices(k)), ...
                               sz, Temp_selected{i}(temp_indices(k)), ...
                               'filled', 'LineWidth', 2);
                    end
                    plotted_events = plotted_events + 1;
                else
                    fprintf('  Event %d: Index out of range (max_soc_idx=%d, soc_length=%d)\n', ...
                        i, max(soc_indices), soc_length);
                end
            else
                fprintf('  Event %d: Insufficient data (soc_length=%d, r_length=%d, dt=%d)\n', ...
                    i, soc_length, r_length, dt(1));
            end
        end
    end
        xlabel('SoC [%]', 'interpreter', 'tex');
        ylabel('Z_{CHG} [m\Omega]', 'interpreter', 'tex');
    xlim([40 80]);
        colormap(flipud(autumn));
    c = colorbar('eastoutside'); caxis([10 32.5]);
        c.Label.String = 'Temperature [°C]';
        c.Label.Interpreter = 'tex';
    set(findall(gcf,'-property','interpreter'),'interpreter','tex');
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex');
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    title(sprintf('Fig 6c - %s Charging Events (Z_CHG vs SoC)', selected_name));
    fprintf('Fig 6c: Plotted %d events\n', plotted_events);
end

%% Figure 6d - Selected charging events (Z_CHG vs Voltage)
% Z_CHG = R_pseudo (충전 임피던스)

if sum(selected_mask) > 0
    figure; box on; hold all
    fprintf('Fig 6d - %s events: %d\n', selected_name, length(R_pseudo_selected));
    fprintf('Debug: dt = [%d, %d], Ts = %.1f\n', dt(1), dt(2), Ts);
    
    plotted_events = 0;
    for i = 1:length(R_pseudo_selected)
        fprintf('Event %d: Voltage length = %d, R_pseudo length = %d\n', ...
            i, length(Voltage_selected{i}), length(R_pseudo_selected{i}));
        
        if length(Voltage_selected{i}) > dt(1) && ...
           ~isempty(R_pseudo_selected{i}) && ...
           length(R_pseudo_selected{i}) > 0
            
            % R_pseudo와 해당 시점의 Voltage 매칭 (전체 범위 사용)
            r_length = length(R_pseudo_selected{i});
            voltage_length = length(Voltage_selected{i});
            
            % 길이 맞추기
            min_length = min(voltage_length, r_length);
            voltage_indices = 1:min_length;
            r_indices = 1:min_length;
            temp_indices = 1:min_length;
            
            % 디버깅: Voltage와 R_pseudo 값 범위 출력
            voltage_values = Voltage_selected{i}(voltage_indices);
            r_values = R_pseudo_selected{i}(r_indices);
            fprintf('  Event %d: Voltage range [%.2f, %.2f], R_pseudo range [%.2f, %.2f]\n', ...
                i, min(voltage_values), max(voltage_values), min(r_values), max(r_values));
            
            % 모든 데이터 포인트를 개별적으로 플롯
            for k = 1:length(voltage_indices)
                scatter(Voltage_selected{i}(voltage_indices(k)), ...
                       R_pseudo_selected{i}(r_indices(k)), ...
                       sz, Temp_selected{i}(temp_indices(k)), ...
                       'filled', 'LineWidth', 2);
            end
            plotted_events = plotted_events + 1;
        end
    end
    xlabel('Voltage [V]', 'interpreter', 'tex');
        ylabel('Z_{CHG} [m\Omega]', 'interpreter', 'tex');
        colormap(flipud(autumn));
    c = colorbar('eastoutside'); caxis([10 32.5]);
        c.Label.String = 'Temperature [°C]';
        c.Label.Interpreter = 'tex';
    set(findall(gcf,'-property','interpreter'),'interpreter','tex');
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex');
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    title(sprintf('Fig 6d - %s Charging Events (Z_CHG vs Voltage)', selected_name));
    fprintf('Fig 6d: Plotted %d events\n', plotted_events);
end

%% Figure 12 - Selected charging events (DV vs SoC)
% DV = Z_CHG / (1000*(time_window*Ts)) * 3600 [V/Ah]

if sum(selected_mask) > 0
    figure; box on; hold all
    fprintf('Fig 12 - %s events: %d\n', selected_name, length(R_pseudo_selected));
    fprintf('Debug: dt = [%d, %d], Ts = %.1f\n', dt(1), dt(2), Ts);
    
    plotted_events = 0;
    for i = 1:length(R_pseudo_selected)
        fprintf('Event %d: SoC length = %d, R_pseudo length = %d\n', ...
            i, length(SoC_selected{i}), length(R_pseudo_selected{i}));
        
        if length(SoC_selected{i}) > dt(1) && ...
           ~isempty(R_pseudo_selected{i}) && ...
           length(R_pseudo_selected{i}) > 0
            
            % R_pseudo와 해당 시점의 SoC 매칭 (dt(1) 사용)
            r_length = length(R_pseudo_selected{i});
            soc_indices = (1:r_length) + dt(1);  % dt(1)만큼 뒤로 이동
            r_indices = 1:r_length;
            temp_indices = (1:r_length) + dt(1);  % 온도도 동일하게 이동
            
            % DV 계산: Z_CHG / (1000*(dt(1)*Ts)) * 3600 [V/Ah]
            DV_values = R_pseudo_selected{i}(r_indices) / (1000 * (dt(1) * Ts)) * 3600;
            
            scatter(SoC_selected{i}(soc_indices), ...
                   DV_values, ...
                   sz, Temp_selected{i}(temp_indices), ...
                   'filled', 'LineWidth', 2);
            plotted_events = plotted_events + 1;
        end
    end
        xlabel('SoC [%]', 'interpreter', 'tex');
    ylabel('DV [V/Ah]', 'interpreter', 'tex');
        colormap(flipud(autumn));
    c = colorbar('eastoutside'); caxis([10 32.5]);
        c.Label.String = 'Temperature [°C]';
        c.Label.Interpreter = 'tex';
    set(findall(gcf,'-property','interpreter'),'interpreter','tex');
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex');
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    title(sprintf('Fig 12 - %s Charging Events (DV vs SoC)', selected_name));
    fprintf('Fig 12: Plotted %d events\n', plotted_events);
end

%% Save results (Original Format)
save(fullfile(saveDir, 'Charging_Events_NewLogic_Fig6_11_12_Optimized.mat'), ...
    'global_eventStruct', 'eventStruct', 'eventCount', ...
    'Time_clean', 'Current_clean', 'Voltage_clean', 'SoC_clean', 'Temp_clean', ...
    'R_pseudo', 'R_pseudo_selected', 'I_ave', 'selected_mask', 'selected_name');

fprintf('Processing complete (OPTIMIZED VERSION)\n');
fprintf('Results saved to: %s\n', saveDir); 