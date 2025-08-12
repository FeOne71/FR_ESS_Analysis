%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rack Data Peak Charge - NewLogic with Fig4_5 Visualization
% Peak detection: NewLogic (idle -> charging)
% Visualization: fig_4_5 style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Rack_raw2mat\New\2025';
yearList = {'2025'};
saveDir = fullfile('G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_Rack_DCIR\DCIR_Charge_NewLogic_Fig4_5');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters (NewLogic 파라미터 사용)
C_nom_cell = 128;
thr = C_nom_cell * 0.01;   % Initial current threshold (A)
dI  = C_nom_cell * 0.2;    % Current change threshold (A)
ddI = 1;                   % Continuous current increase threshold (A)

%% Data Collection
Peaks_all.Time = [];
Peaks_all.Current = [];
Peaks_all.Voltage = [];
Peaks_all.Temperature = [];
Peaks_all.SoC = [];
Peaks_all.Date = [];
Peaks_all.R = [];
Peaks_all.Month = [];

Peaks_monthly = struct();

% Process monthly folders
monthDirs = dir(fullfile(dataDir, '20*'));
for m = 1:length(monthDirs)
    if ~monthDirs(m).isdir, continue; end
    monthPath = fullfile(dataDir, monthDirs(m).name);
    matFiles = dir(fullfile(monthPath, '*.mat'));

    month_num = str2num(monthDirs(m).name(5:6));
    fprintf('  Extracted month_num: %d\n', month_num);
    for f = 1:length(matFiles)
        fprintf('Processing daily file: %s (%d/%d)\n', matFiles(f).name, f, length(matFiles));
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

        N = length(I);

        % Derivative of Current (원본 전류의 미분값)
        for i = 1:length(I)-1
            dI_dt(i) = (I(i+1) - I(i)) / (t(i+1) - t(i));
        end

        % NewLogic Peak Detection (idle -> charging)
        PeakTime = {};
        PeakCurrent = {};
        PeakVoltage = {};
        PeakTemp = {};
        PeakSoC = {};
        z = 1;
        prev_peak_end_idx = 0;  % 이전 피크 끝점 인덱스

        for i = 1:(length(I)-1)
            % 1. idle 구간에서 charging으로 전환되는 지점 찾기
            if I(i) > -thr && I(i) < thr  % idle 구간
                if i+1 <= length(I) && I(i+1) >= thr  % 다음 시점이 charging
                    % 2. 피크 시작점 확인 (idle의 마지막 시점)
                    peak_start_idx = i;
                    
                    % 3. 피크 끝점 찾기 (전류가 더 이상 증가하지 않는 첫 시점)
                    peak_end_idx = peak_start_idx;
                    for j = peak_start_idx+1:length(I)-1
                        if I(j+1) <= I(j)  % 전류 증가가 멈춘 시점
                            peak_end_idx = j;
                            break;
                        end
                    end
                    
                    % 4. 피크 길이 확인 (최소 2개 포인트 이상)
                    if peak_end_idx > peak_start_idx
                        % 5. 중복 피크 방지 (겹침 확인)
                        if peak_start_idx > prev_peak_end_idx
                            % 6. 피크 검출 조건 확인 (NewLogic과 동일)
                            if (I(peak_end_idx) - I(peak_start_idx)) >= dI
                                if (I(peak_start_idx+1) - I(peak_start_idx)) > ddI
                                    flag = 1;
                                    % 구간 내 전류/미분 체크
                                    for zi = peak_start_idx:peak_end_idx-1
                                        if zi <= length(dI_dt) && (dI_dt(zi) < 0 || I(zi+1) < 0)
                                            flag = 0;
                                            break;
                                        end
                                    end
                                    
                                    if flag == 1
                                        % 피크 저장
                                        peak_length = peak_end_idx - peak_start_idx + 1;
                                        PeakTime{z} = t(peak_start_idx:peak_end_idx);
                                        PeakCurrent{z} = I(peak_start_idx:peak_end_idx);
                                        PeakVoltage{z} = V(peak_start_idx:peak_end_idx);
                                        PeakTemp{z} = T_batt(peak_start_idx:peak_end_idx);
                                        PeakSoC{z} = soc(peak_start_idx:peak_end_idx);
                                        z = z + 1;
                                        prev_peak_end_idx = peak_end_idx;  % 끝점 업데이트
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        % Resistance Computation
        for i = 1:length(PeakTime)
            if isempty(PeakTime{i}), continue; end

            DV = (PeakVoltage{i}(end) - PeakVoltage{i}(1));
            DI = (PeakCurrent{i}(end) - PeakCurrent{i}(1));

            if DI > 0 && PeakCurrent{i}(end) > 0
                PeakChg_R_val = (DV / DI) * 1000;
            else
                PeakChg_R_val = NaN;
            end

            if ~isnan(PeakChg_R_val)
                [~, filename, ~] = fileparts(matFiles(f).name);
                date_str = filename(5:end);
                file_date = datetime(date_str, 'InputFormat', 'yyyyMMdd');
                peak_datetime = file_date;
                month_str = sprintf('Month_%02d', month_num);

                Peaks_all.Time = [Peaks_all.Time PeakTime{i}(1)];
                Peaks_all.Current = [Peaks_all.Current PeakCurrent{i}(1)];
                Peaks_all.Voltage = [Peaks_all.Voltage PeakVoltage{i}(1)];
                Peaks_all.Temperature = [Peaks_all.Temperature mean(PeakTemp{i})];
                Peaks_all.SoC = [Peaks_all.SoC mean(PeakSoC{i})];
                Peaks_all.Date = [Peaks_all.Date peak_datetime];
                Peaks_all.R = [Peaks_all.R PeakChg_R_val];
                Peaks_all.Month = [Peaks_all.Month month_num];

                if ~isfield(Peaks_monthly, month_str)
                    Peaks_monthly.(month_str).Time = [];
                    Peaks_monthly.(month_str).Current = [];
                    Peaks_monthly.(month_str).Voltage = [];
                    Peaks_monthly.(month_str).Temperature = [];
                    Peaks_monthly.(month_str).SoC = [];
                    Peaks_monthly.(month_str).Date = [];
                    Peaks_monthly.(month_str).R = [];
                end

                Peaks_monthly.(month_str).Time = [Peaks_monthly.(month_str).Time PeakTime{i}(1)];
                Peaks_monthly.(month_str).Current = [Peaks_monthly.(month_str).Current PeakCurrent{i}(1)];
                Peaks_monthly.(month_str).Voltage = [Peaks_monthly.(month_str).Voltage PeakVoltage{i}(1)];
                Peaks_monthly.(month_str).Temperature = [Peaks_monthly.(month_str).Temperature mean(PeakTemp{i})];
                Peaks_monthly.(month_str).SoC = [Peaks_monthly.(month_str).SoC mean(PeakSoC{i})];
                Peaks_monthly.(month_str).Date = [Peaks_monthly.(month_str).Date peak_datetime];
                Peaks_monthly.(month_str).R = [Peaks_monthly.(month_str).R PeakChg_R_val];
            end
        end
        clear Raw
    end
end

%% Remove outliers
if ~isempty(Peaks_all.R)
    pd = fitdist(Peaks_all.R', 'Normal');
    out_max = pd.mu + 3*pd.sigma;
    out_min = pd.mu - 3*pd.sigma;

    r = 1;
    while r <= length(Peaks_all.R)
        if(Peaks_all.R(r) > out_max || Peaks_all.R(r) < out_min)
            f = fieldnames(Peaks_all);
            for j = 1:length(f)
                Peaks_all.(f{j})(r) = [];
            end
        else
            r = r + 1;
        end
    end
    pd = fitdist(Peaks_all.R', 'Normal');
end

%% Figure 4 - Monthly Histogram and Temperature vs Resistance
if ~isempty(Peaks_all.R)
    month_fields = fieldnames(Peaks_monthly);
    
    % 1. Monthly histogram overlap
    figure(1);
    hold on; box on;
    
    colors = lines(length(month_fields));
    
    for m_idx = 1:length(month_fields)
        month_field = month_fields{m_idx};
        month_R = Peaks_monthly.(month_field).R;
        
        if ~isempty(month_R)
            h = histogram(month_R, 20, 'EdgeColor', colors(m_idx,:), 'FaceColor', colors(m_idx,:), 'FaceAlpha', 0.7);
            hold on;
        end
    end
    
    xlabel('R_{Peak}_{CHG} [m\Omega]'); ylabel('Frequency [-]');
    legend(month_fields, 'Interpreter', 'none');
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    set(findall(gcf,'-property','interpreter'),'interpreter','tex')
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')
    title('Monthly Peak_{Chg} Distribution (NewLogic)');
    
    saveas(gcf, fullfile(saveDir, sprintf('DCIR_Monthly_Histogram_NewLogic_%s.fig', yearList{1})));
    
    % 2. Overall year histogram
    figure(3);
    hold on; box on;
    
    % 전체 연도의 DCIR 값 히스토그램
    h = histogram(Peaks_all.R, 30, 'FaceColor', '#0073C2', 'EdgeColor', 'black', 'LineWidth', 1);
    
    % 통계 정보 계산
    mean_dcir = mean(Peaks_all.R);
    std_dcir = std(Peaks_all.R);
    median_dcir = median(Peaks_all.R);
    
    % 정규분포 곡선 추가
    x_range = linspace(min(Peaks_all.R), max(Peaks_all.R), 100);
    normal_curve = normpdf(x_range, mean_dcir, std_dcir);
    
    % 히스토그램의 최대 높이에 맞춰 정규분포 곡선 스케일링
    max_hist = max(h.Values);
    max_normal = max(normal_curve);
    scaled_normal = normal_curve * (max_hist / max_normal);
    
    % 정규분포 곡선 플롯
    plot(x_range, scaled_normal, 'r-', 'LineWidth', 2, 'DisplayName', 'Normal Distribution');
    
    % 평균선 추가
    xline(mean_dcir, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Mean: %.2f mΩ', mean_dcir));
    
    % 중앙값선 추가
    xline(median_dcir, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Median: %.2f mΩ', median_dcir));
    
    xlabel('DCIR [mΩ]', 'interpreter', 'tex');
    ylabel('Frequency', 'interpreter', 'tex');
    title('Overall Year DCIR Distribution (NewLogic)', 'FontSize', 14);
    legend('Location', 'best');
    set(gca, 'fontsize', 12);
    set(gca, 'ticklabelinterpreter', 'tex');
    
    % 통계 정보 텍스트 박스 추가
    stats_text = sprintf(['Statistics:\n' ...
        'Total Events: %d\n' ...
        'Mean: %.2f mΩ\n' ...
        'Std: %.2f mΩ\n' ...
        'Median: %.2f mΩ\n' ...
        'Min: %.2f mΩ\n' ...
        'Max: %.2f mΩ'], ...
        length(Peaks_all.R), mean_dcir, std_dcir, median_dcir, ...
        min(Peaks_all.R), max(Peaks_all.R));
    
    text(0.02, 0.98, stats_text, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 10, ...
        'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % 히스토그램 저장
    saveas(gcf, fullfile(saveDir, sprintf('DCIR_Overall_Histogram_NewLogic_%s.fig', yearList{1})));
    
    % 2. Monthly temperature vs resistance
    figure(2);
    box on; hold on;
    
    month_labels = {};
    valid_months = 0;
    monthly_x_pos = [];
    monthly_R_values = [];
    monthly_T_values = [];
    
    for m_idx = 1:length(month_fields)
        month_field = month_fields{m_idx};
        month_R = Peaks_monthly.(month_field).R;
        month_T = Peaks_monthly.(month_field).Temperature;
        
        if ~isempty(month_R)
            valid_months = valid_months + 1;
            month_num = str2num(month_field(7:8));
            month_labels{valid_months} = sprintf('2025-%02d', month_num);
            
            x_pos = month_num;
            monthly_x_pos = [monthly_x_pos x_pos];
            monthly_R_values = [monthly_R_values month_R];
            monthly_T_values = [monthly_T_values month_T];
        end
    end
    
    % Plot monthly resistance values
    for m_idx = 1:length(month_fields)
        month_field = month_fields{m_idx};
        month_R = Peaks_monthly.(month_field).R;
        month_T = Peaks_monthly.(month_field).Temperature;
        
        if ~isempty(month_R)
            month_num = str2num(month_field(7:8));
            x_positions = month_num * ones(size(month_R));
            scatter(x_positions, month_R, 100, month_T, 'filled', 'LineWidth', 2);
            hold on;
        end
    end
    
    % Left y-axis: Resistance values
    ylabel('R_{Peak}_{CHG} [m\Omega]'); ylim([0 2]);
    set(gca, 'YTick',(0:0.5:2)); set(gca,'YColor','k');
    
    xlabel('Month');
    
    if ~isempty(monthly_x_pos)
        xticks(monthly_x_pos);
        xticklabels(month_labels);
        xtickangle(45);
        xlim([min(monthly_x_pos)-0.5 max(monthly_x_pos)+0.5]);
    end
    
    % Load monthly ambient temperature data
    ambTempPath = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\JeongEupSi_MonthlyAvg_AmbTemp.csv';
    if exist(ambTempPath, 'file')
        ambTempData = readtable(ambTempPath);
        target_year = 2025;
        
        if target_year == 2023
            target_column = 'Var2';
        elseif target_year == 2025
            target_column = 'Var3';
        elseif target_year == 2025
            target_column = 'Var4';
        end
        
        monthly_temp_data = ambTempData.(target_column);
        monthly_temp_data = monthly_temp_data(2:13);
        monthly_temp_data = double(monthly_temp_data);
        monthly_temp_data = monthly_temp_data(~isnan(monthly_temp_data));
        interp_avgT = monthly_temp_data;
    else
        interp_avgT = [];
    end
    
    % Right y-axis: Temperature
    yyaxis right;
    ylabel('Temperature [°C]');
    set(gca, 'YColor', 'b');
    ylim([0 30]);
    
    temp_x_pos = monthly_x_pos;
    temp_y_values = interp_avgT(monthly_x_pos);
    
    plot(temp_x_pos, temp_y_values, '-o', 'color', 'b', 'linewidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    
    % Calculate actual temperature range for colorbar
    all_temperatures = [];
    for m_idx = 1:length(month_fields)
        month_field = month_fields{m_idx};
        month_T = Peaks_monthly.(month_field).Temperature;
        all_temperatures = [all_temperatures month_T];
    end
    
    temp_min = min(all_temperatures);
    temp_max = max(all_temperatures);
    
    colormap(flipud(autumn));
    c = colorbar('eastoutside', 'Position', [0.95,0.165,0.016,0.75], 'Limits', [temp_min temp_max]);
    c.Label.String = 'Temperature [°C]';
    title('Peak_{Chg} vs Month with Temperature (NewLogic)');
    
    saveas(gcf, fullfile(saveDir, sprintf('Peak_Chg_Monthly_TimeSeries_NewLogic_%s.fig', yearList{1})));
else
    fprintf('No peaks detected.\n');
end

%% Save results
save(fullfile(saveDir, 'Peaks_all_NewLogic.mat'), 'Peaks_all', 'Peaks_monthly');
fprintf('Processing complete\n');
fprintf('Results saved to: %s\n', saveDir); 