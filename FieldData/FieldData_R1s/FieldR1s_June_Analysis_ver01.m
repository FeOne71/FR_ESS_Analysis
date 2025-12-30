%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldR1s_June_Analysis_ver01.m
% ESS Rack01 R1s June Data Analysis
% 
% Key Features:
% 1. Extract June data from each year
% 2. Visualize June data by year (5 subplots)
% 3. Boxplot visualization by year
% 4. Statistical analysis with p-value for trend validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

%% Configuration
% Data directory (from FieldR1s_Integrated_ver03.m)
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_R1s','R1s_Results_Integrated_ver03');

% Output directory
outputDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_R1s','June_Analysis_Results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Load Already Temperature-Corrected Data
fprintf('=== Loading Temperature-Corrected Data ===\n');

% Load the temperature corrected data saved by FieldR1s_Integrated_ver03.m
correctedDataFile = fullfile(saveDir, 'Temperature_Corrected_R1s_Data.mat');

if ~exist(correctedDataFile, 'file')
    fprintf('Error: Temperature_Corrected_R1s_Data.mat not found!\n');
    fprintf('Please run FieldR1s_Integrated_ver03.m first to generate corrected data.\n');
    return;
end

load(correctedDataFile, 'corrected_data');

fprintf('Loaded temperature corrected data:\n');
fprintf('  Total data points: %d\n', length(corrected_data.dates));
fprintf('  Date range: %s to %s\n', ...
    datestr(min(corrected_data.dates), 'yyyy-mm-dd'), ...
    datestr(max(corrected_data.dates), 'yyyy-mm-dd'));

%% Filter June Data Only
fprintf('\n=== Filtering June Data ===\n');

% Extract month from dates
june_idx = month(corrected_data.dates) == 6;

if sum(june_idx) == 0
    fprintf('No June data found!\n');
    return;
end

% Filter June data
june_dates = corrected_data.dates(june_idx);
june_R1s_corrected = corrected_data.R1s_corrected(june_idx);
june_R1s_uncorrected = corrected_data.R1s_uncorrected(june_idx);
june_T = corrected_data.T(june_idx);
june_years = corrected_data.years(june_idx);

fprintf('Found %d June data points\n', length(june_dates));

% Get unique years with June data
unique_years = unique(june_years);
years = cellstr(unique_years); % Keep as cell array of strings
years = sort(years);

fprintf('Years with June data: %s\n', strjoin(years, ', '));

%% Organize June Data by Year
fprintf('\n=== Organizing June Data by Year ===\n');

june_data = struct();

for i = 1:length(years)
    year = years{i};
    field_name = ['year_' year];
    
    % Find indices for this year
    year_str = year;
    year_idx = strcmp(june_years, year_str);
    
    if sum(year_idx) > 0
        june_data.(field_name).R1s_corrected = june_R1s_corrected(year_idx);
        june_data.(field_name).R1s_uncorrected = june_R1s_uncorrected(year_idx);
        june_data.(field_name).T = june_T(year_idx);
        june_data.(field_name).dates = june_dates(year_idx);
        
        fprintf('Year %s: %d June days, R1s corrected range: %.3f-%.3f mΩ\n', ...
            year, length(june_data.(field_name).R1s_corrected), ...
            min(june_data.(field_name).R1s_corrected)*1000, ...
            max(june_data.(field_name).R1s_corrected)*1000);
    end
end

%% Visualization 1: June Data by Year (5 Subplots)
fprintf('\n=== Creating June Data Visualization by Year ===\n');

figure('Position', [100, 100, 1600, 1000]);

nYears = length(years);
nCols = min(3, nYears);
nRows = ceil(nYears / nCols);

for i = 1:nYears
    year = years{i};
    field_name = ['year_' year];
    
    subplot(nRows, nCols, i);
    
    R1s_corr = june_data.(field_name).R1s_corrected;
    T_data = june_data.(field_name).T;
    dates = june_data.(field_name).dates;
    
    % Sort by date
    [dates_sorted, sortIdx] = sort(dates);
    R1s_corr_sorted = R1s_corr(sortIdx);
    T_data_sorted = T_data(sortIdx);
    
    % Create sequential x-axis (days within June)
    xAxis = (1:length(dates_sorted))';
    
    % Plot scatter points with temperature color
    scatter(xAxis, R1s_corr_sorted*1000, 50, T_data_sorted, 'filled', 'HandleVisibility', 'off');
    hold on;
    
    % Add trend line
    if length(R1s_corr_sorted) > 1
        mdl = fitlm(xAxis, R1s_corr_sorted*1000, 'Intercept', true);
        trendLine = mdl.Fitted;
        plot(xAxis, trendLine, 'g--', 'LineWidth', 3, 'DisplayName', sprintf('Trend (R²=%.4f)', mdl.Rsquared.Ordinary));
    end
    
    % Colorbar for battery temperature
    colormap(gca, flipud(autumn));
    c = colorbar(gca, 'eastoutside');
    c.Label.String = 'Temperature (°C)';
    c.Label.FontSize = 10;
    c.Label.FontWeight = 'bold';
    c.Limits = [min(T_data_sorted), max(T_data_sorted)];
    
    % Formatting
    xlabel('Day of June', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('R1s (mΩ) - Corrected to 25°C', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Year %s June Data', year), 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Set thick axes
    ax = gca;
    ax.LineWidth = 2;
    ax.FontWeight = 'bold';
    
    % Set y-axis limits
    ylim([0.7, 0.9]);
    
    % Show only trend legend
    legend('Location', 'best', 'FontSize', 10);
    
    % Format x-axis - show day numbers
    xticks(1:max(1, floor(length(xAxis)/5)):length(xAxis));
end

sgtitle('ESS Rack01 June R1s Data by Year (Temperature Corrected)', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, fullfile(outputDir, 'June_R1s_by_Year.fig'));
fprintf('June data visualization saved to: %s\n', outputDir);

%% Visualization 2: Boxplot of June Data by Year
fprintf('\n=== Creating June Boxplot by Year ===\n');

figure('Position', [100, 100, 1200, 600]);

% Prepare data for boxplot
boxplot_data = [];
boxplot_groups = [];
boxplot_labels = {};

for i = 1:length(years)
    year = years{i};
    field_name = ['year_' year];
    
    R1s_corr = june_data.(field_name).R1s_corrected * 1000; % Convert to mΩ
    if length(R1s_corr) > 0
        boxplot_data = [boxplot_data; R1s_corr];
        boxplot_groups = [boxplot_groups; repmat(i, length(R1s_corr), 1)];
        boxplot_labels{end+1} = year;
    end
end

% Create boxplot
if ~isempty(boxplot_data)
    boxplot(boxplot_data, boxplot_groups, 'Notch', 'on', 'Labels', boxplot_labels);
end

% Formatting
ylabel('R1s (mΩ) - Temperature Corrected to 25°C', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
title('ESS Rack01 June R1s Distribution by Year', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Set thick axes
ax = gca;
ax.LineWidth = 2;
ax.FontWeight = 'bold';

% Save figure
saveas(gcf, fullfile(outputDir, 'June_R1s_Boxplot_by_Year.fig'));
fprintf('June boxplot saved to: %s\n', outputDir);

%% Statistical Analysis with p-value
fprintf('\n=== Statistical Analysis with p-value ===\n');

% Prepare data for statistical analysis
all_years_numeric = [];
all_R1s_corrected = [];

for i = 1:length(years)
    year = years{i};
    field_name = ['year_' year];
    
    R1s_corr = june_data.(field_name).R1s_corrected * 1000; % Convert to mΩ
    year_num = str2double(year);
    
    all_years_numeric = [all_years_numeric; repmat(year_num, length(R1s_corr), 1)];
    all_R1s_corrected = [all_R1s_corrected; R1s_corr];
end

% Sort by year
[all_years_numeric, sortIdx] = sort(all_years_numeric);
all_R1s_corrected = all_R1s_corrected(sortIdx);

% 1. Linear Trend Test (Time series regression)
fprintf('\n--- 1. Linear Trend Test ---\n');
if length(all_R1s_corrected) > 2
    % Create time variable (years since first year)
    first_year = min(all_years_numeric);
    time_var = all_years_numeric - first_year;
    
    % Linear regression: R1s = a + b*time
    mdl_trend = fitlm(time_var, all_R1s_corrected);
    slope_trend = mdl_trend.Coefficients.Estimate(2);
    intercept_trend = mdl_trend.Coefficients.Estimate(1);
    pvalue_trend = mdl_trend.Coefficients.pValue(2);
    r_squared_trend = mdl_trend.Rsquared.Ordinary;
    
    fprintf('Linear Trend Model: R1s = %.4f + %.4f * (Year - %d)\n', ...
        intercept_trend, slope_trend, first_year);
    fprintf('Slope: %.4f mΩ/year\n', slope_trend);
    fprintf('p-value (slope): %.4e\n', pvalue_trend);
    fprintf('R²: %.4f\n', r_squared_trend);
    
    if pvalue_trend < 0.05
        fprintf('Result: SIGNIFICANT trend (p < 0.05) - R1s changes significantly over time\n');
        if slope_trend > 0
            fprintf('Direction: INCREASING trend (R1s increases with time)\n');
        else
            fprintf('Direction: DECREASING trend (R1s decreases with time)\n');
        end
    else
        fprintf('Result: NO SIGNIFICANT trend (p >= 0.05) - R1s does not change significantly over time\n');
    end
end

% 2. One-way ANOVA (Compare means across years) - Method A: Standard ANOVA
fprintf('\n--- 2. One-way ANOVA (Year Comparison) ---\n');
fprintf('Method A: ANOVA - Testing if mean R1s values differ significantly across years\n\n');

% Prepare data for ANOVA in the format requested:
% Y = june_data.R1s_corrected (vector of all corrected R1s values)
% Group = june_data.Year (group vector indicating which year each value belongs to)
Y = [];
Group = [];

for i = 1:length(years)
    year = years{i};
    field_name = ['year_' year];
    
    R1s_corr = june_data.(field_name).R1s_corrected * 1000; % Convert to mΩ
    if length(R1s_corr) > 0
        Y = [Y; R1s_corr];
        % Create group vector: use numeric year for grouping
        year_num = str2double(year);
        Group = [Group; repmat(year_num, length(R1s_corr), 1)];
    end
end

if length(unique(Group)) > 1
    % Perform ANOVA: "Do the means of different year groups differ significantly?"
    fprintf('Data Preparation:\n');
    fprintf('  Y (R1s_corrected): %d data points\n', length(Y));
    fprintf('  Group (Year): %d unique groups\n', length(unique(Group)));
    fprintf('  Groups: %s\n', num2str(unique(Group)'));
    fprintf('\n');
    
    % Perform ANOVA
    [p_anova, tbl_anova, stats_anova] = anova1(Y, Group, 'off');
    
    fprintf('ANOVA Results:\n');
    fprintf('  F-statistic: %.4f\n', tbl_anova{2, 5});
    fprintf('  p-value: %.4e\n', p_anova);
    fprintf('\n');
    
    % Interpretation
    if p_anova < 0.05
        fprintf('Interpretation:\n');
        fprintf('  p-value < 0.05: REJECT the null hypothesis\n');
        fprintf('  Conclusion: There is a STATISTICALLY SIGNIFICANT difference\n');
        fprintf('               between the mean R1s values across different years.\n');
        fprintf('  This is the desired result - it proves that R1s changes significantly over time.\n');
    else
        fprintf('Interpretation:\n');
        fprintf('  p-value >= 0.05: FAIL TO REJECT the null hypothesis\n');
        fprintf('  Conclusion: There is NO STATISTICALLY SIGNIFICANT difference\n');
        fprintf('               between the mean R1s values across different years.\n');
    end
    
    fprintf('\n');
    
    % Print summary statistics by year
    fprintf('Year-by-Year Statistics:\n');
    fprintf('Year\t\tMean (mΩ)\tStd (mΩ)\tN\n');
    fprintf('----\t\t---------\t--------\t-\n');
    unique_years = unique(Group);
    for i = 1:length(unique_years)
        year_num = unique_years(i);
        year_idx = (Group == year_num);
        year_data = Y(year_idx);
        fprintf('%d\t\t%.4f\t\t%.4f\t\t%d\n', ...
            year_num, mean(year_data), std(year_data), length(year_data));
    end
    fprintf('\n');
else
    fprintf('Insufficient groups for ANOVA (need at least 2 groups)\n');
end

% 3. Kruskal-Wallis Test (Non-parametric alternative)
fprintf('--- 3. Kruskal-Wallis Test (Non-parametric) ---\n');
if exist('Y', 'var') && length(unique(Group)) > 1
    % Use the same Y and Group data
    [p_kruskal, tbl_kruskal, stats_kruskal] = kruskalwallis(Y, Group, 'off');
    
    fprintf('Kruskal-Wallis Test Results:\n');
    fprintf('  Chi-square statistic: %.4f\n', tbl_kruskal{2, 5});
    fprintf('  p-value: %.4e\n', p_kruskal);
    
    if p_kruskal < 0.05
        fprintf('  Result: SIGNIFICANT difference between years (p < 0.05)\n');
    else
        fprintf('  Result: NO SIGNIFICANT difference between years (p >= 0.05)\n');
    end
    fprintf('\n');
end

% 4. Pairwise Comparisons (Tukey-Kramer)
fprintf('\n--- 4. Pairwise Comparisons (Tukey-Kramer) ---\n');
if exist('p_anova', 'var') && exist('stats_anova', 'var') && p_anova < 0.05 && length(unique(Group)) > 1
    % Perform multiple comparisons
    [c, m, h, gnames] = multcompare(stats_anova, 'CType', 'tukey-kramer', 'Display', 'off');
    
    fprintf('Pairwise Comparisons (Tukey-Kramer):\n');
    fprintf('Group 1\t\tGroup 2\t\tDifference\t95%% CI Lower\t95%% CI Upper\tp-value\n');
    fprintf('-------\t\t-------\t\t----------\t------------\t------------\t-------\n');
    
    unique_years_str = gnames; % multcompare에서 반환된 그룹 이름 사용
    
    for i = 1:size(c, 1)
        group1_idx = c(i, 1);
        group2_idx = c(i, 2);
        
        % CORRECT column assignments based on MATLAB documentation
        ci_lower   = c(i, 3);  % Column 3: CI Lower bound
        diff       = c(i, 4);  % Column 4: Estimated Difference
        ci_upper   = c(i, 5);  % Column 5: CI Upper bound
        p_pairwise = c(i, 6);  % Column 6: p-value
        
        % Get year strings from gnames
        group1_year_str = unique_years_str{group1_idx};
        group2_year_str = unique_years_str{group2_idx};
        
        % Format p-value: use scientific notation
        p_str = sprintf('%.4e', p_pairwise);
        
        fprintf('%s\t\t%s\t\t%+.4f\t\t%+.4f\t\t%+.4f\t\t%s\n', ...
            group1_year_str, group2_year_str, diff, ci_lower, ci_upper, p_str);
    end
    fprintf('\n');
    
elseif exist('p_anova', 'var') && p_anova >= 0.05
    fprintf('Skipped: ANOVA was not significant, so pairwise comparisons were not performed.\n');
    fprintf('\n');
else
    fprintf('Skipped: Less than 2 groups, pairwise comparison not applicable.\n');
    fprintf('\n');
end

% 5. Year-to-Year Change Rate
fprintf('\n--- 5. Year-to-Year Change Rate ---\n');
if length(years) >= 2
    fprintf('Year-to-Year Change:\n');
    fprintf('Year\t\tMean R1s (mΩ)\tChange from Previous (mΩ)\tChange Rate (%%)\n');
    fprintf('----\t\t-------------\t------------------------\t---------------\n');
    
    prev_mean = [];
    for i = 1:length(years)
        year = years{i};
        field_name = ['year_' year];
        
        R1s_corr = june_data.(field_name).R1s_corrected * 1000;
        curr_mean = mean(R1s_corr);
        
        if isempty(prev_mean)
            fprintf('%s\t\t%.4f\t\tN/A\t\t\tN/A\n', year, curr_mean);
        else
            change = curr_mean - prev_mean;
            change_rate = (change / prev_mean) * 100;
            fprintf('%s\t\t%.4f\t\t%.4f\t\t\t%.2f%%\n', year, curr_mean, change, change_rate);
        end
        
        prev_mean = curr_mean;
    end
end

%% Save Results
fprintf('\n=== Saving Results ===\n');

% Save June data structure
save(fullfile(outputDir, 'June_Data_Analysis.mat'), 'june_data', 'years');

% Save statistical results
if exist('pvalue_trend', 'var')
    stats_results = struct();
    stats_results.trend_slope = slope_trend;
    stats_results.trend_pvalue = pvalue_trend;
    stats_results.trend_r_squared = r_squared_trend;
    
    if exist('p_anova', 'var')
        stats_results.anova_pvalue = p_anova;
        stats_results.anova_fstat = tbl_anova{2, 5};
    end
    
    if exist('p_kruskal', 'var')
        stats_results.kruskal_pvalue = p_kruskal;
    end
    
    save(fullfile(outputDir, 'June_Statistical_Analysis.mat'), 'stats_results');
end

fprintf('Analysis completed!\n');
fprintf('Results saved to: %s\n', outputDir);

