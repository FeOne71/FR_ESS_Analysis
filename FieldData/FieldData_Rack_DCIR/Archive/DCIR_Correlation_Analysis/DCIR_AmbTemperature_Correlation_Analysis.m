%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCIR Ambient Temperature Correlation Analysis
% Analysis of correlation between ambient temperature and DCIR values
% Based on 2023NewLogic results with Fig4_5 visualization style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Load the saved results from 2023NewLogic analysis
% Define the directory containing the results
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\Charge';
folder = 'DCIR_Charge_2023_NewLogic';
folder_path = fullfile(baseDir, folder);

fprintf('Loading data from: %s\n', folder);

if ~exist(folder_path, 'dir')
    error('Folder does not exist: %s', folder_path);
end

% Find all mat files in the folder
mat_files = dir(fullfile(folder_path, '*.mat'));
fprintf('Found %d mat files\n', length(mat_files));

% Initialize data collection
all_temperatures = [];
all_resistances = [];
all_dates = [];
all_years = [];
all_racks = [];
all_months = [];

for file_idx = 1:length(mat_files)
    file_path = fullfile(folder_path, mat_files(file_idx).name);
    fprintf('Loading file: %s\n', mat_files(file_idx).name);
    
    try
        load(file_path);
        
        % Debug: Check what variables are loaded
        fprintf('Variables loaded: %s\n', strjoin(who, ', '));
        
        % Extract data based on the structure
        if exist('global_eventStruct', 'var')
            fprintf('global_eventStruct found\n');
            
            % Process global_eventStruct
            racks = fieldnames(global_eventStruct);
            fprintf('Racks found: %s\n', strjoin(racks, ', '));
            
            for rack_idx = 1:length(racks)
                rack_name = racks{rack_idx};
                fprintf('Processing rack: %s\n', rack_name);
                
                if isstruct(global_eventStruct.(rack_name))
                    years = fieldnames(global_eventStruct.(rack_name));
                    fprintf('Years found: %s\n', strjoin(years, ', '));
                    
                    for year_idx = 1:length(years)
                        year_name = years{year_idx};
                        fprintf('Processing year: %s\n', year_name);
                        
                        if isstruct(global_eventStruct.(rack_name).(year_name))
                            dates = fieldnames(global_eventStruct.(rack_name).(year_name));
                            fprintf('Dates found: %s\n', strjoin(dates, ', '));
                            
                            for date_idx = 1:length(dates)
                                date_name = dates{date_idx};
                                fprintf('Processing date: %s\n', date_name);
                                
                                if isstruct(global_eventStruct.(rack_name).(year_name).(date_name))
                                    events = fieldnames(global_eventStruct.(rack_name).(year_name).(date_name));
                                    fprintf('Events found: %s\n', strjoin(events, ', '));
                                    
                                    for event_idx = 1:length(events)
                                        event_name = events{event_idx};
                                        event_data = global_eventStruct.(rack_name).(year_name).(date_name).(event_name);
                                        
                                        % Debug: Check event data structure
                                        if event_idx == 1
                                            fprintf('Event fields: %s\n', strjoin(fieldnames(event_data), ', '));
                                        end
                                        
                                        % Extract DCIR value
                                        if isfield(event_data, 'PeakChgR') && ~isnan(event_data.PeakChgR)
                                            % Extract temperature from T_seq field
                                            if isfield(event_data, 'T_seq') && ~isempty(event_data.T_seq)
                                                temp_value = mean(event_data.T_seq);
                                                fprintf('Event %s: T_seq length = %d, temp_value = %.2f\n', ...
                                                    event_name, length(event_data.T_seq), temp_value);
                                            else
                                                temp_value = NaN;
                                                fprintf('Event %s: T_seq field missing or empty\n', event_name);
                                            end
                                            
                                            % Extract month from date
                                            date_str = event_data.date;
                                            if length(date_str) >= 10
                                                month_str = date_str(6:7);
                                                month_num = str2double(month_str);
                                            else
                                                month_num = NaN;
                                            end
                                            
                                            % Store data
                                            all_temperatures = [all_temperatures; temp_value];
                                            all_resistances = [all_resistances; event_data.PeakChgR];
                                            all_dates = [all_dates; {event_data.date}];
                                            all_years = [all_years; {event_data.year}];
                                            all_racks = [all_racks; {event_data.rack_name}];
                                            all_months = [all_months; month_num];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        else
            fprintf('global_eventStruct not found in file\n');
        end
        
    catch ME
        fprintf('Error loading file %s: %s\n', mat_files(file_idx).name, ME.message);
    end
end

fprintf('Data loading complete. Total events found: %d\n', length(all_resistances));

% Debug: Check collected data
fprintf('Temperature data summary:\n');
fprintf('Total temperatures: %d\n', length(all_temperatures));
fprintf('NaN temperatures: %d\n', sum(isnan(all_temperatures)));
fprintf('Valid temperatures: %d\n', sum(~isnan(all_temperatures)));
if ~isempty(all_temperatures)
    fprintf('Temperature range: %.2f to %.2f °C\n', min(all_temperatures(~isnan(all_temperatures))), max(all_temperatures(~isnan(all_temperatures))));
end

%% Extract temperature and resistance data
% Remove any NaN values and zero resistance values
valid_indices = ~isnan(all_temperatures) & ~isnan(all_resistances) & (all_resistances > 0) & ~isnan(all_months);
temperatures_clean = all_temperatures(valid_indices);
resistances_clean = all_resistances(valid_indices);
dates_clean = all_dates(valid_indices);
years_clean = all_years(valid_indices);
racks_clean = all_racks(valid_indices);
months_clean = all_months(valid_indices);

fprintf('Total data points: %d\n', length(all_temperatures));
fprintf('Valid data points: %d\n', length(temperatures_clean));
if ~isempty(temperatures_clean)
    fprintf('Temperature range: %.2f to %.2f °C\n', min(temperatures_clean), max(temperatures_clean));
    fprintf('Resistance range: %.2f to %.2f mΩ\n', min(resistances_clean), max(resistances_clean));
end

% Display data summary by month
unique_months = unique(months_clean);
fprintf('\nData summary by month:\n');
for month_idx = 1:length(unique_months)
    month_data = months_clean == unique_months(month_idx);
    month_count = sum(month_data);
    month_temp_mean = mean(temperatures_clean(month_data));
    month_res_mean = mean(resistances_clean(month_data));
    fprintf('Month %d: %d events, Avg Temp = %.2f°C, Avg Res = %.3f mΩ\n', ...
        unique_months(month_idx), month_count, month_temp_mean, month_res_mean);
end

%% Load monthly ambient temperature data
ambTempPath = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\JeongEupSi_MonthlyAvg_AmbTemp.csv';
if exist(ambTempPath, 'file')
    ambTempData = readtable(ambTempPath);
    target_year = 2023;  % 2023NewLogic 데이터에 맞춤
    
    if target_year == 2023
        target_column = 'Var2';
    elseif target_year == 2024
        target_column = 'Var3';
    elseif target_year == 2025
        target_column = 'Var4';
    end
    
    monthly_temp_data = ambTempData.(target_column);
    monthly_temp_data = monthly_temp_data(2:13);  % 1월부터 12월
    monthly_temp_data = double(monthly_temp_data);
    monthly_temp_data = monthly_temp_data(~isnan(monthly_temp_data));
    
    fprintf('Loaded ambient temperature data for year %d\n', target_year);
    fprintf('Monthly ambient temperatures: ');
    for i = 1:length(monthly_temp_data)
        fprintf('%.1f°C ', monthly_temp_data(i));
    end
    fprintf('\n');
else
    error('Ambient temperature file not found: %s', ambTempPath);
end

%% Create monthly data structure (Fig4_5 style)
Peaks_monthly = struct();

for month_idx = 1:length(unique_months)
    month_num = unique_months(month_idx);
    month_data = months_clean == month_num;
    
    if sum(month_data) > 0
        month_str = sprintf('Month_%02d', month_num);
        Peaks_monthly.(month_str).R = resistances_clean(month_data);
        Peaks_monthly.(month_str).Temperature = temperatures_clean(month_data);
        Peaks_monthly.(month_str).AmbientTemp = monthly_temp_data(month_num) * ones(sum(month_data), 1);
    end
end

%% Statistical Analysis - Ambient Temperature vs Resistance
fprintf('\n=== Ambient Temperature vs Resistance Analysis ===\n');

% Get month fields
month_fields = fieldnames(Peaks_monthly);

% Calculate overall correlation using monthly averages
monthly_avg_dcir = [];
monthly_amb_temp = [];

for m_idx = 1:length(month_fields)
    month_field = month_fields{m_idx};
    month_R = Peaks_monthly.(month_field).R;
    month_amb_T = Peaks_monthly.(month_field).AmbientTemp;
    
    if ~isempty(month_R)
        monthly_avg_dcir = [monthly_avg_dcir; mean(month_R)];
        monthly_amb_temp = [monthly_amb_temp; mean(month_amb_T)];
    end
end

% Overall correlation using monthly averages
if length(monthly_avg_dcir) > 2
    [overall_corr, overall_p] = corrcoef(monthly_avg_dcir, monthly_amb_temp);
    overall_r = overall_corr(1,2);
    overall_p_value = overall_p(1,2);
    
    fprintf('\nOverall correlation (monthly averages):\n');
    fprintf('Correlation coefficient (r): %.4f\n', overall_r);
    fprintf('P-value: %.6f\n', overall_p_value);
    fprintf('Number of months: %d\n', length(monthly_avg_dcir));
    
    % Determine significance
    if overall_p_value < 0.001
        significance = ' (p < 0.001)';
    elseif overall_p_value < 0.01
        significance = ' (p < 0.01)';
    elseif overall_p_value < 0.05
        significance = ' (p < 0.05)';
    else
        significance = 'ns (p >= 0.05)';
    end
    fprintf('Significance: %s\n', significance);
else
    overall_r = NaN;
    overall_p_value = NaN;
    significance = 'Insufficient data';
    fprintf('\nInsufficient data for overall correlation analysis\n');
end

% Calculate monthly statistics for display
monthly_correlations = [];
monthly_p_values = [];
monthly_counts = [];
monthly_amb_temp = [];
monthly_res_mean = [];

for m_idx = 1:length(month_fields)
    month_field = month_fields{m_idx};
    month_R = Peaks_monthly.(month_field).R;
    month_amb_T = Peaks_monthly.(month_field).AmbientTemp;
    
    monthly_counts(m_idx) = length(month_R);
    monthly_amb_temp(m_idx) = mean(month_amb_T);
    monthly_res_mean(m_idx) = mean(month_R);
    
    fprintf('Month %s: n = %d, Amb Temp = %.1f°C, Avg Res = %.3f mΩ\n', ...
        month_field, monthly_counts(m_idx), monthly_amb_temp(m_idx), monthly_res_mean(m_idx));
end

%% Visualization - Fig4_5 style
if ~isempty(month_fields)
    % Monthly temperature vs resistance (Fig4_5 style)
    figure(1);
    box on; hold on;
    
    month_labels = {};
    valid_months = 0;
    monthly_x_pos = [];
    monthly_R_values = [];
    monthly_T_values = [];
    
    for m_idx = 1:length(month_fields)
        month_field = month_fields{m_idx};
        month_R = Peaks_monthly.(month_field).R;
        month_amb_T = Peaks_monthly.(month_field).AmbientTemp;
        
        if ~isempty(month_R)
            valid_months = valid_months + 1;
            month_num = str2num(month_field(7:8));
            month_labels{valid_months} = sprintf('2023-%02d', month_num);
            
            x_pos = month_num;
            monthly_x_pos = [monthly_x_pos x_pos];
            % Store as column vectors to avoid dimension mismatch
            monthly_R_values = [monthly_R_values; month_R(:)];
            monthly_T_values = [monthly_T_values; month_amb_T(:)];
        end
    end
    
    % Plot monthly resistance values
    for m_idx = 1:length(month_fields)
        month_field = month_fields{m_idx};
        month_R = Peaks_monthly.(month_field).R;
        month_amb_T = Peaks_monthly.(month_field).AmbientTemp;
        
        if ~isempty(month_R)
            month_num = str2num(month_field(7:8));
            x_positions = month_num * ones(size(month_R));
            scatter(x_positions, month_R, 100, month_amb_T, 'filled', 'LineWidth', 2);
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
    
    % Right y-axis: Ambient Temperature
    yyaxis right;
    ylabel('Monthly Ambient Temperature [°C]','FontSize',12);
    set(gca, 'YColor', 'b');
    ylim([0 30]);
    
    temp_x_pos = monthly_x_pos;
    temp_y_values = monthly_temp_data(monthly_x_pos);
    
    plot(temp_x_pos, temp_y_values, '-o', 'color', 'b', 'linewidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    
    % Calculate actual temperature range for colorbar
    all_amb_temperatures = [];
    for m_idx = 1:length(month_fields)
        month_field = month_fields{m_idx};
        month_amb_T = Peaks_monthly.(month_field).AmbientTemp;
        all_amb_temperatures = [all_amb_temperatures; month_amb_T];
    end
    
    temp_min = min(all_amb_temperatures);
    temp_max = max(all_amb_temperatures);
    
    colormap(flipud(autumn));
    c = colorbar('eastoutside', 'Position', [0.95,0.165,0.016,0.75], 'Limits', [temp_min temp_max]);
    c.Label.String = 'Battery Temperature [°C]';
    title(sprintf('Peak_{Chg} vs 2023 Monthly Ambient Temperature (r = %.3f, p = %.4f %s)', ...
        overall_r, overall_p_value, significance));
    
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    set(findall(gcf,'-property','interpreter'),'interpreter','tex')
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')
    
else
    fprintf('No peaks detected.\n');
end

%% Additional Analysis - Scatter plot of ambient temperature vs resistance
figure(2);
hold on; box on;

% Create scatter plot of all data points
all_amb_temp_scatter = [];
all_res_scatter = [];

for m_idx = 1:length(month_fields)
    month_field = month_fields{m_idx};
    month_R = Peaks_monthly.(month_field).R;
    month_amb_T = Peaks_monthly.(month_field).AmbientTemp;
    
    if ~isempty(month_R)
        all_amb_temp_scatter = [all_amb_temp_scatter; month_amb_T];
        all_res_scatter = [all_res_scatter; month_R];
    end
end

if ~isempty(all_amb_temp_scatter)
    scatter(all_amb_temp_scatter, all_res_scatter, 50, 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Add regression line
    valid_data = ~isnan(all_amb_temp_scatter) & ~isnan(all_res_scatter);
    if sum(valid_data) > 2
        p = polyfit(all_amb_temp_scatter(valid_data), all_res_scatter(valid_data), 1);
        x_range = linspace(min(all_amb_temp_scatter), max(all_amb_temp_scatter), 100);
        y_regression = polyval(p, x_range);
        plot(x_range, y_regression, 'r-', 'LineWidth', 2);
        
        % Calculate correlation
        [corr_coeff, corr_p] = corrcoef(all_amb_temp_scatter(valid_data), all_res_scatter(valid_data));
        r_value = corr_coeff(1,2);
        p_value = corr_p(1,2);
        
        title(sprintf('2023 Ambient Temperature vs DCIR\nr = %.3f, p = %.4f', r_value, p_value));
    end
end

xlabel('Ambient Temperature [°C]', 'FontSize', 12);
ylabel('DCIR [mΩ]', 'FontSize', 12);
grid on;

%% Save results
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\DCIR_Correlation_Analysis';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

save(fullfile(saveDir, 'Ambient_Temperature_Correlation_Analysis_2023NewLogic.mat'), ...
    'temperatures_clean', 'resistances_clean', 'dates_clean', 'years_clean', 'racks_clean', 'months_clean', ...
    'Peaks_monthly', 'monthly_correlations', 'monthly_p_values', 'monthly_counts', ...
    'monthly_amb_temp', 'monthly_res_mean', 'overall_r', 'overall_p_value', 'significance');

% Save figures
saveas(figure(1), fullfile(saveDir, 'Ambient_Temperature_Correlation_Fig4_5_Style.fig'));
saveas(figure(2), fullfile(saveDir, 'Ambient_Temperature_Scatter_Plot.fig'));

fprintf('\nAnalysis complete! Results saved to: %s\n', saveDir); 