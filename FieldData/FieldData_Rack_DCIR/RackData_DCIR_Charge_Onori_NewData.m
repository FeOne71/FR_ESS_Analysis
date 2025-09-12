%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rack Data DCIR Charge - New Logic (idle -> charging)
% 2023-2025 Data
% Peak detection: NewLogic (idle -> charging)
% 1초 샘플링 데이터, 새로운 피크 검출 로직 적용
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;

%% Directory
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New';
yearList = {'2023', '2024'}; %% Find -> Replace 할것
% Output roots for charge/discharge
yearStr = strjoin(yearList, '_');
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR','TotalPeak_NewData');

if ~exist(saveDir, 'dir'); mkdir(saveDir); end

%% Parameters
C_nom_cell = 128;
thr = [C_nom_cell*0.07, C_nom_cell*0.1]; % Idle entry thresholds (A): thr(1)=charge, thr(2)=discharge
dI_chg = C_nom_cell * 0.20; % Charge: ΔI > dI_chg (A)
dI_dch = C_nom_cell * 0.15; % Discharge: ΔI < -dI_dch (A)
ddI = 1;
te = 6;  % Moving average window: ±10 seconds (1초 샘플링 기준)
dt = 1;   % Peak window length in samples (1 Hz sampling)

%% Visualization variables
Fontsize = 10;
LineWidth = 3;

rackNames_all = {'Rack01'};
global_eventStruct = struct();
global_eventStruct.ChgPeak = struct();
global_eventStruct.DchPeak = struct();
for r = 1:length(rackNames_all)
    global_eventStruct.ChgPeak.(rackNames_all{r}) = struct();
    global_eventStruct.DchPeak.(rackNames_all{r}) = struct();
end



for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    fprintf('Processing year: %s\n', year);
    eventStruct = struct();
    eventCount = 0;
    yearPath = fullfile(dataDir, year);
    fprintf('Year path: %s\n', yearPath);
    fprintf('Year path exists: %d\n', exist(yearPath, 'dir'));
    monthDirs = dir(fullfile(yearPath, '20*'));
    fprintf('Found %d month directories\n', length(monthDirs));
    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        fprintf('Processing month: %s\n', monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, '*.mat'));
        fprintf('Found %d mat files in %s\n', length(matFiles), monthDirs(m).name);
        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            fprintf('Loading file: %s\n', matFiles(f).name);
            load(matFilePath);
            % fprintf('Raw structure fields: %s\n', strjoin(fieldnames(Raw), ', '));
            rackNames = rackNames_all;
            rackName = 'Rack01';  % NewData는 단일 랙 처리

            % Initialize daily tracking for this file
            daily_eventCount = 0; % 날짜별 이벤트 번호 초기화

            % 2023년도 데이터는 단일 Raw 구조체 (Rack 구조 없음)
            fprintf('Data length: %d points\n', length(Raw.Date_Time_seconds));
            t = Raw.Date_Time;
            I = Raw.DCCurrent;
            V = Raw.CVavg;
            soc = Raw.SOC_BMS;
            T_batt = Raw.MTavg;

            % Convert time (duration) to seconds-from-start (no conditionals)
            t = seconds(t - t(1));
            N = length(I);

            % Moving Average Filter on Current (fig_4_5.m 방식)
            I_filt = zeros(N,1);
            Pre = 0*ones(te,1);
            Post = 0*ones(te,1);
            I_calc = [Pre;I;Post];
            for i=1:N
                for m=0:2*te
                    I_filt(i) = I_filt(i)+I_calc(i+m);
                end
                I_filt(i) = I_filt(i)/(2*te+1);
            end

            % Derivative of Current (원본 전류의 미분값)
            for i = 1:length(I)-1
                dI_dt(i) = (I(i+1) - I(i)) / (t(i+1) - t(i));
            end

            % Moving Average on Current Derivative
            N_dI = length(dI_dt);
            filt_dI_dt = zeros(N_dI,1);
            Pre = 0*ones(te,1);
            Post = 0*ones(te,1);
            calc_dI_dt = [Pre;dI_dt';Post];
            for i=1:N_dI
                for m=0:2*te
                    filt_dI_dt(i) = filt_dI_dt(i)+calc_dI_dt(i+m);
                end
                filt_dI_dt(i) = filt_dI_dt(i)/(2*te+1);
            end

            % New Peak Detection Logic (레퍼런스 방식 적용)
            PeakTime = {};
            PeakCurrent = {};
            PeakVoltage = {};
            PeakSoC = {};
            PeakTemp = {};
            z = 1;
            prev_peak_end_idx = 0;  % 이전 피크 끝점 인덱스

%% Charge Peak Detection Logic (레퍼런스 방식 적용)

            for i = 1:(length(I)-dt)
                % 1단계: 휴지 상태 확인
                if I(i+dt) - I(i) > dI_chg
                    % if I(i) > -thr(1) && I(i) < thr(1)
                    % 2단계: dt초 후 전류 변화량 확인
                    % if I(i+dt) - I(i) > dI_chg
                    if I(i) > -thr(1) && I(i) < thr(1)
                        % 3단계: 전압 범위 필터링 추가
                        % if V(i) >= 3.5
                            % 피크 저장 (2-point DCIR: 시작점과 dt초 후)
                            PeakTime{z} = [t(i), t(i+dt)];
                            PeakCurrent{z} = [I(i), I(i+dt)];
                            PeakVoltage{z} = [V(i), V(i+dt)];
                            PeakSoC{z} = soc(i);
                            PeakTemp{z} = T_batt(i);
                            z = z + 1;
                        % end
                    end
                end
            end

            % Accumulate daily peaks (단일 데이터)
            % Charge/Discharge 합산 개수 기준으로 원본 저장(오버랩용)
            total_chg = length(PeakTime);
            % Discharge 피크는 아래에서 검출됨. 이후 합산하여 저장 판단

            % 원본 데이터 저장 (시간, 전압, 전류, 온도, SOC)
            daily_original_data = {struct('time', Raw.Date_Time_seconds, 'voltage', Raw.CVavg, ...
                'current', Raw.DCCurrent, 'temperature', Raw.MTavg, 'soc', Raw.SOC_BMS)};

            % Resistance Computation - CHARGE
            for i = 1:length(PeakTime)
                DV = PeakVoltage{i}(end) - PeakVoltage{i}(1);
                DI = PeakCurrent{i}(end) - PeakCurrent{i}(1);
                dcir_val = (DV / DI) * 1000;
                % Extract date from filename (for event naming before use)
                date_str = matFiles(f).name(5:12);  % Extract YYYYMMDD from Raw_YYYYMMDD.mat
                formatted_date = sprintf('%s-%s-%s', date_str(1:4), date_str(5:6), date_str(7:8));
                eventCount = eventCount + 1;
                daily_eventCount = daily_eventCount + 1;
                date_key = sprintf('E%s', strrep(formatted_date,'-',''));
                evtName = sprintf('Event%d', daily_eventCount);
                eventStruct.(evtName).rack_name = 'Rack01';
                eventStruct.(evtName).year = year;
                eventStruct.(evtName).date = formatted_date;
                eventStruct.(evtName).type = 'charge';
                eventStruct.(evtName).start_idx = find(t >= PeakTime{i}(1), 1);
                eventStruct.(evtName).end_idx = find(t >= PeakTime{i}(end), 1);
                eventStruct.(evtName).charge_duration = length(PeakTime{i});
                eventStruct.(evtName).avg_current = mean(PeakCurrent{i});
                eventStruct.(evtName).t_seq = Raw.Date_Time(eventStruct.(evtName).start_idx:eventStruct.(evtName).end_idx);
                eventStruct.(evtName).I_seq = PeakCurrent{i};
                eventStruct.(evtName).V_seq = PeakVoltage{i};
                eventStruct.(evtName).P_seq = Raw.DCPower(eventStruct.(evtName).start_idx:eventStruct.(evtName).end_idx);
                % Peak 저항 저장
                eventStruct.(evtName).PeakChgR = dcir_val;
            end

            % Debug peak detection
            if ~isempty(PeakTime)
                avg_voltage = mean(cellfun(@(x) x(1), PeakVoltage));
                fprintf('Rack %s: %d charge peaks detected (Voltage: %.2fV)\n', rackName, length(PeakTime), avg_voltage);
            else
                fprintf('Rack %s: No charge peaks detected\n', rackName);
            end

%% Discharge Peak Detection Logic (레퍼런스 방식 적용)
            PeakTimeDis = {}; PeakCurrentDis = {}; PeakVoltageDis = {}; PeakSoCDis = {}; PeakTempDis = {}; zd = 1; prev_end_dis = 0;

            for i = 1:(length(I)-dt)
                % 1단계: 휴지 상태 확인
                if I(i+dt) - I(i) < -dI_dch
                    % if I(i) > -thr(2) && I(i) < thr(2)
                    % 2단계: dt초 후 전류 변화량 확인
                    % if I(i+dt) - I(i) < -dI_dch
                    if I(i) > -thr(2) && I(i) < thr(2)
                        % 3단계: 전압 범위 필터링 추가
                        % if V(i) <= 3.9
                            % 피크 저장 (2-point DCIR: 시작점과 dt초 후)
                            PeakTimeDis{zd} = [t(i), t(i+dt)];
                            PeakCurrentDis{zd} = [I(i), I(i+dt)];
                            PeakVoltageDis{zd} = [V(i), V(i+dt)];
                            PeakSoCDis{zd} = soc(i);
                            PeakTempDis{zd} = T_batt(i);
                            zd = zd + 1;
                        % end
                    end
                end
            end

            % Resistance Computation - DISCHARGE
            for i = 1:length(PeakTimeDis)
                DV = PeakVoltageDis{i}(end) - PeakVoltageDis{i}(1);
                DI = PeakCurrentDis{i}(end) - PeakCurrentDis{i}(1);
                dis_val = (DV / DI) * 1000;

                % save discharge event (date for event naming before use)
                date_str = matFiles(f).name(5:12);
                formatted_date = sprintf('%s-%s-%s', date_str(1:4), date_str(5:6), date_str(7:8));
                eventCount = eventCount + 1;
                daily_eventCount = daily_eventCount + 1;
                date_key = sprintf('E%s', strrep(formatted_date,'-',''));
                evtName = sprintf('Event%d', daily_eventCount);
                eventStruct.(evtName).rack_name = 'Rack01';
                eventStruct.(evtName).year = year;
                eventStruct.(evtName).date = formatted_date;
                eventStruct.(evtName).type = 'discharge';
                eventStruct.(evtName).start_idx = find(t >= PeakTimeDis{i}(1), 1);
                eventStruct.(evtName).end_idx = find(t >= PeakTimeDis{i}(end), 1);
                eventStruct.(evtName).charge_duration = length(PeakTimeDis{i});
                eventStruct.(evtName).avg_current = mean(PeakCurrentDis{i});
                eventStruct.(evtName).t_seq = Raw.Date_Time(eventStruct.(evtName).start_idx:eventStruct.(evtName).end_idx);
                eventStruct.(evtName).I_seq = PeakCurrentDis{i};
                eventStruct.(evtName).V_seq = PeakVoltageDis{i};
                eventStruct.(evtName).P_seq = Raw.DCPower(eventStruct.(evtName).start_idx:eventStruct.(evtName).end_idx);
                eventStruct.(evtName).PeakDisR = dis_val;
            end

            % Debug discharge peak detection
            if ~isempty(PeakTimeDis)
                avg_voltage_dis = mean(cellfun(@(x) x(1), PeakVoltageDis));
                fprintf('Rack %s: %d discharge peaks detected (Voltage: %.2fV)\n', rackName, length(PeakTimeDis), avg_voltage_dis);
            else
                fprintf('Rack %s: No discharge peaks detected\n', rackName);
            end


        end

        % Save structure - CHARGE (simplified structure)
        eventNames = fieldnames(eventStruct);
        for i = 1:length(eventNames)
            evt = eventStruct.(eventNames{i});
            if isfield(evt,'type') && strcmp(evt.type,'charge')
                date_key = sprintf('E%s', strrep(evt.date, '-', ''));
                if ~isfield(global_eventStruct.ChgPeak, 'Rack01')
                    global_eventStruct.ChgPeak.Rack01 = struct();
                end
                if ~isfield(global_eventStruct.ChgPeak.Rack01, date_key)
                    global_eventStruct.ChgPeak.Rack01.(date_key) = struct();
                end
                global_eventStruct.ChgPeak.Rack01.(date_key).(eventNames{i}) = evt;
            end
        end

        % Save structure - DISCHARGE (simplified structure)
        for i = 1:length(eventNames)
            evt = eventStruct.(eventNames{i});
            if isfield(evt,'type') && strcmp(evt.type,'discharge')
                date_key = sprintf('E%s', strrep(evt.date, '-', ''));
                if ~isfield(global_eventStruct.DchPeak, 'Rack01')
                    global_eventStruct.DchPeak.Rack01 = struct();
                end
                if ~isfield(global_eventStruct.DchPeak.Rack01, date_key)
                    global_eventStruct.DchPeak.Rack01.(date_key) = struct();
                end
                global_eventStruct.DchPeak.Rack01.(date_key).(eventNames{i}) = evt;
            end
        end
    end
end % Close for year_idx loop

%% Yearly comparison - Histograms (overlay) and Boxplots for multiple years
% This section visualizes per-rack charge/discharge peak resistance across years
if length(yearList) >= 1
    % Colors for up to 3 years
    yearColors = {'#0073C2','#FF8C00','#77AC30'}; % extend if needed
    % Prepare date pattern for each year
    datePatterns = cell(1, length(yearList));
    for yi = 1:length(yearList)
        datePatterns{yi} = sprintf('E%s', yearList{yi}); % E2023, E2024 패턴
    end

    % ---------- Overlay Histograms: CHARGE ----------
    rackName = 'Rack01';
    perYearVals = cell(1, length(yearList));
    perYearStats = cell(1, length(yearList));
    % Collect values per year first (NaN removed, 3σ filter)
    for yi = 1:length(yearList)
        vals = [];
        if isfield(global_eventStruct.ChgPeak, rackName)
            date_names = fieldnames(global_eventStruct.ChgPeak.(rackName));
            year_dates = date_names(startsWith(date_names, datePatterns{yi}));
            fprintf('Rack %s, Year %s: Found %d dates\n', rackName, yearList{yi}, length(year_dates));
            for d = 1:length(year_dates)
                dayEvents = global_eventStruct.ChgPeak.(rackName).(year_dates{d});
                evtNamesDay = fieldnames(dayEvents);
                fprintf('  Date %s: Found %d events\n', year_dates{d}, length(evtNamesDay));
                for e = 1:length(evtNamesDay)
                    evt = dayEvents.(evtNamesDay{e});
                    if isfield(evt, 'PeakChgR')
                        vals = [vals, evt.PeakChgR]; %#ok<AGROW>
                    end
                end
            end
        end
        fprintf('Rack %s, Year %s: Collected %d charge values\n', rackName, yearList{yi}, length(vals));
        vals = vals(~isnan(vals) & vals ~= 0 & vals > 0);
        if ~isempty(vals)
            mu_v = mean(vals); sd_v = std(vals);
            inlier = (vals >= mu_v - 3*sd_v) & (vals <= mu_v + 3*sd_v);
            vals_f = vals(inlier);

            perYearVals{yi} = vals_f;
            perYearStats{yi} = struct('mu', mean(vals_f), 'sd', std(vals_f), 'n', numel(vals_f));
        else
            perYearVals{yi} = [];
            perYearStats{yi} = struct('mu', NaN, 'sd', NaN, 'n', 0);
        end
    end
    % Determine common bin edges from combined data
    allVals = [perYearVals{:}];
    if ~isempty(allVals)
        min_v = min(allVals); max_v = max(allVals);
        if min_v == max_v, min_v = 0; max_v = min_v + 1e-3; end
        edges = linspace(min_v, max_v, 31);
    else
        edges = linspace(0, 1.0, 31);
    end
    figName = sprintf('YearCompare_Hist_Chg_%s', rackName);
    figure('Name', figName, 'NumberTitle', 'off'); hold on; grid on; box on;
    for yi = 1:length(yearList)
        v = perYearVals{yi}; if isempty(v), continue; end
        st = perYearStats{yi}; col = yearColors{min(yi, numel(yearColors))};
        % filled histogram (opaque)
        histogram(v, edges, 'FaceColor', col, 'EdgeColor', col, 'LineWidth', 0.5);
        % normal curve scaled to counts
        bw = edges(2)-edges(1); xg = linspace(edges(1), edges(end), 200);
        yg = normpdf(xg, st.mu, max(st.sd, eps)) * st.n * bw;
        plot(xg, yg, 'Color', col, 'LineWidth', 2);
        % per-year stats text
        text(0.02, 0.95 - 0.06*yi, sprintf('%s: n=%d, mean= %.4f m\\Omega, std= %.4f m\\Omega', yearList{yi}, st.n, st.mu, st.sd), 'Units','normalized', 'Color', col, 'FontSize', Fontsize);
    end
    xlabel('Resistance [m\Omega]', 'interpreter','tex'); ylabel('Frequency','interpreter','tex');
    title(sprintf('Charge Peaks - %s (Years: %s)', rackName, strjoin(yearList, ', ')));
    xlim([0 1.0]); set(gca,'fontsize',Fontsize,'ticklabelinterpreter','tex');
    set(gcf,'Position',[100 100 800 600]);
    saveas(gcf, fullfile(saveDir, sprintf('YearCompare_Hist_Chg_%s_%s.fig', rackName, strjoin(yearList,'_'))));

    % ---------- Overlay Histograms: DISCHARGE ----------
    rackName = 'Rack01';
    perYearVals = cell(1, length(yearList));
    perYearStats = cell(1, length(yearList));
    for yi = 1:length(yearList)
        vals = [];
        if isfield(global_eventStruct.DchPeak, rackName)
            date_names = fieldnames(global_eventStruct.DchPeak.(rackName));
            year_dates = date_names(startsWith(date_names, datePatterns{yi}));
            fprintf('Rack %s, Year %s: Found %d dates (discharge)\n', rackName, yearList{yi}, length(year_dates));
            for d = 1:length(year_dates)
                dayEvents = global_eventStruct.DchPeak.(rackName).(year_dates{d});
                evtNamesDay = fieldnames(dayEvents);
                fprintf('  Date %s: Found %d events (discharge)\n', year_dates{d}, length(evtNamesDay));
                for e = 1:length(evtNamesDay)
                    evt = dayEvents.(evtNamesDay{e});
                    vals = [vals, evt.PeakDisR]; %#ok<AGROW>
                end
            end
        end
        fprintf('Rack %s, Year %s: Collected %d discharge values\n', rackName, yearList{yi}, length(vals));
        vals = vals(~isnan(vals) & vals ~= 0 & vals > 0);
        if ~isempty(vals)
            mu_v = mean(vals); sd_v = std(vals);
            inlier = (vals >= mu_v - 3*sd_v) & (vals <= mu_v + 3*sd_v);
            vals_f = vals(inlier);

            perYearVals{yi} = vals_f;
            perYearStats{yi} = struct('mu', mean(vals_f), 'sd', std(vals_f), 'n', numel(vals_f));
        else
            perYearVals{yi} = [];
            perYearStats{yi} = struct('mu', NaN, 'sd', NaN, 'n', 0);
        end
    end
    allVals = [perYearVals{:}];
    if ~isempty(allVals)
        min_v = min(allVals); max_v = max(allVals);
        if min_v == max_v, min_v = 0; max_v = min_v + 1e-3; end
        edges = linspace(min_v, max_v, 31);
    else
        edges = linspace(0, 1.0, 31);
    end
    figName = sprintf('YearCompare_Hist_Dch_%s', rackName);
    figure('Name', figName, 'NumberTitle', 'off'); hold on; grid on; box on;
    for yi = 1:length(yearList)
        v = perYearVals{yi}; if isempty(v), continue; end
        st = perYearStats{yi}; col = yearColors{min(yi, numel(yearColors))};
        histogram(v, edges, 'FaceColor', col, 'EdgeColor', col, 'LineWidth', 0.5);
        bw = edges(2)-edges(1); xg = linspace(edges(1), edges(end), 200);
        yg = normpdf(xg, st.mu, max(st.sd, eps)) * st.n * bw;
        plot(xg, yg, 'Color', col, 'LineWidth', 2);
        text(0.02, 0.95 - 0.06*yi, sprintf('%s: n=%d, mean= %.4f m\\Omega, std= %.4f m\\Omega', yearList{yi}, st.n, st.mu, st.sd), 'Units','normalized', 'Color', col, 'FontSize', Fontsize);
    end
    xlabel('Resistance [m\Omega]', 'interpreter','tex'); ylabel('Frequency','interpreter','tex');
    title(sprintf('Discharge Peaks - %s (Years: %s)', rackName, strjoin(yearList, ', ')));
    xlim([0 1.0]); set(gca,'fontsize',Fontsize,'ticklabelinterpreter','tex');
    set(gcf,'Position',[100 100 800 600]);
    saveas(gcf, fullfile(saveDir, sprintf('YearCompare_Hist_Dch_%s_%s.fig', rackName, strjoin(yearList,'_'))));

    % ---------- Boxplots: CHARGE (single plot for Rack01) ----------
    figure('Name','YearCompare_Boxplot_Charge','NumberTitle','off');
    hold on; grid on; box on;
    valsAll = []; grp = [];
    for yi = 1:length(yearList)
        vals = [];
        if isfield(global_eventStruct.ChgPeak, 'Rack01')
            date_names = fieldnames(global_eventStruct.ChgPeak.Rack01);
            year_dates = date_names(startsWith(date_names, datePatterns{yi}));
            for d = 1:length(year_dates)
                dayEvents = global_eventStruct.ChgPeak.Rack01.(year_dates{d});
                evtNamesDay = fieldnames(dayEvents);
                for e = 1:length(evtNamesDay)
                    evt = dayEvents.(evtNamesDay{e});
                    if isfield(evt, 'PeakChgR')
                        vals = [vals, evt.PeakChgR]; %#ok<AGROW>
                    end
                end
            end
        end
        vals = vals(~isnan(vals));
        valsAll = [valsAll, vals]; %#ok<AGROW>
        grp = [grp; repmat(string(yearList{yi}), numel(vals), 1)]; %#ok<AGROW>
    end
    if ~isempty(valsAll)
        % 연도별 데이터 분리 및 3σ 이상치 제거
        valsAll_f = []; grp_f = []; year_labels = [];
        start_idx = 1;
        for yi = 1:length(yearList)
            vals = [];
            if isfield(global_eventStruct.ChgPeak, 'Rack01')
                date_names = fieldnames(global_eventStruct.ChgPeak.Rack01);
                year_dates = date_names(startsWith(date_names, datePatterns{yi}));
                for d = 1:length(year_dates)
                    dayEvents = global_eventStruct.ChgPeak.Rack01.(year_dates{d});
                    evtNamesDay = fieldnames(dayEvents);
                    for e = 1:length(evtNamesDay)
                        evt = dayEvents.(evtNamesDay{e});
                        if isfield(evt, 'PeakChgR')
                            vals = [vals, evt.PeakChgR]; %#ok<AGROW>
                        end
                    end
                end
            end
            vals = vals(~isnan(vals) & vals ~= 0 & vals > 0);
            if ~isempty(vals)
                mu_v = mean(vals); sd_v = std(vals);
                inlier = (vals >= mu_v - 3*sd_v) & (vals <= mu_v + 3*sd_v);
                vals_filtered = vals(inlier);
                valsAll_f = [valsAll_f, vals_filtered];
                grp_f = [grp_f; repmat(string(yearList{yi}), numel(vals_filtered), 1)];
                year_labels = [year_labels, repmat(string(yearList{yi}), 1, numel(vals_filtered))];
            end
        end

        if ~isempty(valsAll_f)
            % Use daboxplot with empty boxes, whiskers, mean, and linkline
            daboxplot(valsAll_f, 'groups', grp_f, 'xtlabels', yearList, 'fill', 0, 'mean', 1, 'linkline', 1);
            % Add manual connection lines between years
            hold on;
            for yi = 1:length(yearList)
                if yi < length(yearList)
                    % Get mean values for consecutive years
                    vals1 = valsAll_f(grp_f == string(yearList{yi}));
                    vals2 = valsAll_f(grp_f == string(yearList{yi+1}));
                    if ~isempty(vals1) && ~isempty(vals2)
                        mean1 = mean(vals1);
                        mean2 = mean(vals2);
                        plot([yi, yi+1], [mean1, mean2], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    end
                end
            end
            hold off;
            % Add legend to clarify what lines represent
            legend('Median', 'Mean', 'Location', 'best');
        else
            daboxplot(valsAll, 'groups', grp, 'xtlabels', yearList, 'fill', 0, 'mean', 1, 'linkline', 1);
            % Add manual connection lines between years
            hold on;
            for yi = 1:length(yearList)
                if yi < length(yearList)
                    % Get mean values for consecutive years
                    vals1 = valsAll(grp == string(yearList{yi}));
                    vals2 = valsAll(grp == string(yearList{yi+1}));
                    if ~isempty(vals1) && ~isempty(vals2)
                        mean1 = mean(vals1);
                        mean2 = mean(vals2);
                        plot([yi, yi+1], [mean1, mean2], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    end
                end
            end
            hold off;
            % Add legend to clarify what lines represent
            legend('Median', 'Mean', 'Location', 'best');
        end

        ylabel('R [m\Omega]','interpreter','tex');
        ylim([0 Inf]);

        % Add statistical text
        grp_cat = categorical(year_labels, string(yearList));
        cats = categories(grp_cat);
        for ci = 1:numel(cats)
            vals_ci = valsAll_f(grp_cat == cats{ci});
            if isempty(vals_ci), continue; end
            mu_ci = mean(vals_ci); sd_ci = std(vals_ci); n_ci = numel(vals_ci);
            text(0.02, 0.95 - 0.08*ci, sprintf('%s: n=%d, mean=%.4f, std=%.4f', string(cats{ci}), n_ci, mu_ci, sd_ci), ...
                'Units','normalized', 'FontSize', Fontsize, 'BackgroundColor', 'white', 'Margin', 1);
        end
    else
        text(0.5,0.5,'No data','HorizontalAlignment','center');
    end
    % Add p-value to title: ANOVA for >=3 years, Welch t-test for 2 years
    pStr = '';
    if ~isempty(valsAll)
        % 원시 데이터로 p-value 계산 (이상치 포함)
        valsAll_raw = []; grp_raw = [];
        for yi = 1:length(yearList)
            vals = [];
            if isfield(global_eventStruct.ChgPeak, 'Rack01')
                date_names = fieldnames(global_eventStruct.ChgPeak.Rack01);
                year_dates = date_names(startsWith(date_names, datePatterns{yi}));
                for d = 1:length(year_dates)
                    dayEvents = global_eventStruct.ChgPeak.Rack01.(year_dates{d});
                    evtNamesDay = fieldnames(dayEvents);
                    for e = 1:length(evtNamesDay)
                        evt = dayEvents.(evtNamesDay{e});
                        if isfield(evt, 'PeakChgR')
                            vals = [vals, evt.PeakChgR]; %#ok<AGROW>
                        end
                    end
                end
            end
            vals = vals(~isnan(vals) & vals ~= 0 & vals > 0); % 0, NaN, 음수 제거 (이상치 포함)
            valsAll_raw = [valsAll_raw, vals];
            grp_raw = [grp_raw; repmat(string(yearList{yi}), numel(vals), 1)];
        end
        if ~isempty(valsAll_raw)
            grp_cat_raw = categorical(grp_raw, string(yearList));
            cats_raw = categories(grp_cat_raw);
            if numel(cats_raw) >= 3
                pOmni = anova1(valsAll_raw, grp_cat_raw, 'off');
                pStr = sprintf(' (ANOVA p=%.3f)', pOmni);
            else
                x1 = valsAll_raw(grp_cat_raw == cats_raw{1}); x2 = valsAll_raw(grp_cat_raw == cats_raw{2});
                [~, pPair] = ttest2(x1, x2, 'Vartype', 'unequal'); % Welch t-test
                pStr = sprintf(' (Welch p=%.3f)', pPair);
            end
        end
    end
    title(sprintf('Rack01 - Charge%s', pStr)); set(gca,'fontsize',Fontsize,'ticklabelinterpreter','tex');
    ylim([0 Inf]);
end
set(gcf,'Position',[100 100 800 600]);
% Add title using annotation at the top
annotation('textbox', [0.1 0.95 0.8 0.05], 'String', sprintf('Charge Peaks - Yearly Boxplots (%s)', strjoin(yearList, ', ')), ...
    'FontSize', 14, 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'BackgroundColor', 'none');
saveas(gcf, fullfile(saveDir, sprintf('YearCompare_Boxplot_Charge_%s.fig', strjoin(yearList,'_'))));

% ---------- Boxplots: DISCHARGE (single plot for Rack01) ----------
figure('Name','YearCompare_Boxplot_Discharge','NumberTitle','off');
hold on; grid on; box on;
valsAll = []; grp = [];
for yi = 1:length(yearList)
    vals = [];
    if isfield(global_eventStruct.DchPeak, 'Rack01')
        date_names = fieldnames(global_eventStruct.DchPeak.Rack01);
        year_dates = date_names(startsWith(date_names, datePatterns{yi}));
        for d = 1:length(year_dates)
            dayEvents = global_eventStruct.DchPeak.Rack01.(year_dates{d});
            evtNamesDay = fieldnames(dayEvents);
            for e = 1:length(evtNamesDay)
                evt = dayEvents.(evtNamesDay{e});
                vals = [vals, evt.PeakDisR]; %#ok<AGROW>
            end
        end
    end
    vals = vals(~isnan(vals));
    valsAll = [valsAll, vals]; %#ok<AGROW>
    grp = [grp; repmat(string(yearList{yi}), numel(vals), 1)]; %#ok<AGROW>
end

if ~isempty(valsAll)
    % 연도별 데이터 분리 및 3σ 이상치 제거
    valsAll_f = []; grp_f = []; year_labels = [];
    start_idx = 1;
    for yi = 1:length(yearList)
        vals = [];
        if isfield(global_eventStruct.DchPeak, 'Rack01')
            date_names = fieldnames(global_eventStruct.DchPeak.Rack01);
            year_dates = date_names(startsWith(date_names, datePatterns{yi}));
            for d = 1:length(year_dates)
                dayEvents = global_eventStruct.DchPeak.Rack01.(year_dates{d});
                evtNamesDay = fieldnames(dayEvents);
                for e = 1:length(evtNamesDay)
                    evt = dayEvents.(evtNamesDay{e});
                    vals = [vals, evt.PeakDisR]; %#ok<AGROW>
                end
            end
        end
        vals = vals(~isnan(vals) & vals ~= 0 & vals > 0);
        if ~isempty(vals)
            mu_v = mean(vals); sd_v = std(vals);
            inlier = (vals >= mu_v - 3*sd_v) & (vals <= mu_v + 3*sd_v);
            vals_filtered = vals(inlier);
            valsAll_f = [valsAll_f, vals_filtered];
            grp_f = [grp_f; repmat(string(yearList{yi}), numel(vals_filtered), 1)];
            year_labels = [year_labels, repmat(string(yearList{yi}), 1, numel(vals_filtered))];
        end
    end

    if ~isempty(valsAll_f)
        % Use daboxplot with empty boxes, whiskers, mean, and linkline
        daboxplot(valsAll_f, 'groups', grp_f, 'xtlabels', yearList, 'fill', 0, 'mean', 1, 'linkline', 1);
        % Add manual connection lines between years
        hold on;
        for yi = 1:length(yearList)
            if yi < length(yearList)
                % Get mean values for consecutive years
                vals1 = valsAll_f(grp_f == string(yearList{yi}));
                vals2 = valsAll_f(grp_f == string(yearList{yi+1}));
                if ~isempty(vals1) && ~isempty(vals2)
                    mean1 = mean(vals1);
                    mean2 = mean(vals2);
                    plot([yi, yi+1], [mean1, mean2], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
            end
        end
        hold off;
        % Add legend to clarify what lines represent
        legend('Median', 'Mean', 'Location', 'best');
    else
        daboxplot(valsAll, 'groups', grp, 'xtlabels', yearList, 'fill', 0, 'mean', 1, 'linkline', 1);
        % Add manual connection lines between years
        hold on;
        for yi = 1:length(yearList)
            if yi < length(yearList)
                % Get mean values for consecutive years
                vals1 = valsAll(grp == string(yearList{yi}));
                vals2 = valsAll(grp == string(yearList{yi+1}));
                if ~isempty(vals1) && ~isempty(vals2)
                    mean1 = mean(vals1);
                    mean2 = mean(vals2);
                    plot([yi, yi+1], [mean1, mean2], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
            end
        end
        hold off;
        % Add legend to clarify what lines represent
        legend('Median', 'Mean', 'Location', 'best');
    end

    ylabel('R [m\Omega]','interpreter','tex');
    ylim([0 Inf]);

    % Add statistical text
    grp_cat = categorical(year_labels, string(yearList));
    cats = categories(grp_cat);
    for ci = 1:numel(cats)
        vals_ci = valsAll_f(grp_cat == cats{ci});
        if isempty(vals_ci), continue; end
        mu_ci = mean(vals_ci); sd_ci = std(vals_ci); n_ci = numel(vals_ci);
        text(0.02, 0.95 - 0.08*ci, sprintf('%s: n=%d, mean=%.4f, std=%.4f', string(cats{ci}), n_ci, mu_ci, sd_ci), ...
            'Units','normalized', 'FontSize', Fontsize, 'BackgroundColor', 'white', 'Margin', 1);
    end
else
    text(0.5,0.5,'No data','HorizontalAlignment','center');
end
% Add p-value to title: ANOVA for >=3 years, Welch t-test for 2 years
pStr = '';
if ~isempty(valsAll)
    % 원시 데이터로 p-value 계산 (이상치 포함)
    valsAll_raw = []; grp_raw = [];
    for yi = 1:length(yearList)
        vals = [];
        if isfield(global_eventStruct.DchPeak, 'Rack01')
            date_names = fieldnames(global_eventStruct.DchPeak.Rack01);
            year_dates = date_names(startsWith(date_names, datePatterns{yi}));
            for d = 1:length(year_dates)
                dayEvents = global_eventStruct.DchPeak.Rack01.(year_dates{d});
                evtNamesDay = fieldnames(dayEvents);
                for e = 1:length(evtNamesDay)
                    evt = dayEvents.(evtNamesDay{e});
                    vals = [vals, evt.PeakDisR]; %#ok<AGROW>
                end
            end
        end
        vals = vals(~isnan(vals) & vals ~= 0 & vals > 0); % 0, NaN, 음수 제거 (이상치 포함)
        valsAll_raw = [valsAll_raw, vals];
        grp_raw = [grp_raw; repmat(string(yearList{yi}), numel(vals), 1)];
    end
    if ~isempty(valsAll_raw)
        grp_cat_raw = categorical(grp_raw, string(yearList));
        cats_raw = categories(grp_cat_raw);
        if numel(cats_raw) >= 3
            pOmni = anova1(valsAll_raw, grp_cat_raw, 'off');
            pStr = sprintf(' (ANOVA p=%.3f)', pOmni);
        else
            x1 = valsAll_raw(grp_cat_raw == cats_raw{1}); x2 = valsAll_raw(grp_cat_raw == cats_raw{2});
            [~, pPair] = ttest2(x1, x2, 'Vartype', 'unequal'); % Welch t-test
            pStr = sprintf(' (Welch p=%.3f)', pPair);
        end
    end
end
title(sprintf('Rack01 - Discharge%s', pStr)); set(gca,'fontsize',Fontsize,'ticklabelinterpreter','tex');
ylim([0 inf]);
set(gcf,'Position',[100 100 800 600]);
% Add title using annotation at the top
annotation('textbox', [0.1 0.95 0.8 0.05], 'String', sprintf('Discharge Peaks - Yearly Boxplots (%s)', strjoin(yearList, ', ')), ...
    'FontSize', 14, 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'BackgroundColor', 'none');
saveas(gcf, fullfile(saveDir, sprintf('YearCompare_Boxplot_Discharge_%s.fig', strjoin(yearList,'_'))));


% Close if length(yearList) >= 1

save(fullfile(saveDir, sprintf('Peak_NewData_%s_total.mat', yearStr)), 'global_eventStruct');
fprintf('Processing complete\n');
fprintf('Results saved to: %s\n', saveDir);

