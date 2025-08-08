%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rack Data Peak Charge - Onori Method with Fig4_5 Visualization
% Peak detection: Charge_Onori logic
% Visualization: fig_4_5 style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory
dataDir = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Rack_raw2mat\New\2025';
saveDir = fullfile('G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_Rack_DCIR\DCIR_Charge_Onori_Fig4_5');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters
C_nom_cell = 128;
dt = 5;                   % Current monotonic increase interval
thr = C_nom_cell * 0.05;  % Initial current threshold (A)
dI =  C_nom_cell * 0.2;    % Current change threshold after 5s (A)
ddI = 1;                  % Continuous current increase threshold (A)

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

        % Filter on Current
        te = 4;
        I_filt = zeros(N,1);
        Pre = 0*ones(te/2,1);
        Post = 0*ones(te/2,1);
        I_calc = [Pre; I; Post];
        for i = 1:N
            for m = 0:te
                I_filt(i) = I_filt(i) + I_calc(i+m);
            end
            I_filt(i) = I_filt(i)/(te+1);
        end

        % Derivative of Current
        for i = 1:length(I)-1
            dI_dt(i) = (I(i+1) - I(i)) / (t(i+1) - t(i));
        end

        % MA on Current Derivative
        N_dI = length(dI_dt);
        filt_dI_dt = zeros(N_dI,1);
        Pre = 0*ones(te/2,1);
        Post = 0*ones(te/2,1);
        calc_dI_dt = [Pre; dI_dt'; Post];
        for i = 1:N_dI
            for m = 0:te
                filt_dI_dt(i) = filt_dI_dt(i) + calc_dI_dt(i+m);
            end
            filt_dI_dt(i) = filt_dI_dt(i)/(te+1);
        end

        % Driving Peaks Identification
        PeakTime = {};
        PeakCurrent = {};
        PeakVoltage = {};
        PeakTemp = {};
        PeakSoC = {};
        z = 1;

        for i = 1:(length(I_filt)-dt)
            if (I_filt(i+dt) - I_filt(i)) > dI
                if I(i) > -thr && I(i) < thr
                    if (I_filt(i+1) - I_filt(i)) > ddI
                        flag = 1;
                        for zi = 1:dt
                            if filt_dI_dt(i+zi-1) < 0 || I(i+zi) < 0
                                flag = 0;
                            end
                        end
                        if flag == 1
                            if I(i+dt) > I(i) && I(i+dt) > 0
                                PeakTime{z} = t(i:i+dt-1);
                                PeakCurrent{z} = I(i:i+dt-1);
                                PeakVoltage{z} = V(i:i+dt-1);
                                PeakTemp{z} = T_batt(i:i+dt-1);
                                PeakSoC{z} = soc(i:i+dt-1);
                                z = z + 1;
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
    title('Monthly Peak_{Chg} Distribution');
    
    saveas(gcf, fullfile(saveDir, 'DCIR_Monthly_Histogram.fig'));
    
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
            month_labels{valid_months} = sprintf('2024-%02d', month_num);
            
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
    ylabel('R_{Peak}_{CHG} [m\Omega]','FontSize',15); ylim([0 2]);
    set(gca, 'YTick',(0:0.5:2)); set(gca,'YColor','k');    
    xlabel('Month','FontSize',15);
    
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
        target_year = 2023;
        
        if target_year == 2023
            target_column = 'Var2';
        elseif target_year == 2024
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
    ylabel('Monthly Temperature [°C]','FontSize',12);
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
    title('Peak_{Chg} vs Monthly Ambient Temperature');
    
    saveas(gcf, fullfile(saveDir, 'Peak_Chg_Monthly_TimeSeries.fig'));
else
    fprintf('No peaks detected.\n');
end

%% Save results
save(fullfile(saveDir, 'Peaks_all.mat'), 'Peaks_all', 'Peaks_monthly');
fprintf('Processing complete\n');
fprintf('Results saved to: %s\n', saveDir);