%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCIR Temperature Correlation Analysis
% Analysis of correlation between ambient temperature and DCIR values
% Based on NewLogic Fig4_5 results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Load the saved results from NewLogic Fig4_5
% Load the saved peaks data
load('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\DCIR_Charge_Onori_Fig4_5\Peaks_all.mat');


%% Extract temperature and resistance data
% Extract temperature and resistance values
temperatures = Peaks_all.Temperature;
resistances = Peaks_all.R;

% Remove any NaN values and zero resistance values
valid_indices = ~isnan(temperatures) & ~isnan(resistances) & (resistances > 0);
temperatures_clean = temperatures(valid_indices);
resistances_clean = resistances(valid_indices);

fprintf('Total data points: %d\n', length(temperatures));
fprintf('Valid data points: %d\n', length(temperatures_clean));
fprintf('Temperature range: %.2f to %.2f °C\n', min(temperatures_clean), max(temperatures_clean));
fprintf('Resistance range: %.2f to %.2f mΩ\n', min(resistances_clean), max(resistances_clean));

%% Statistical Analysis
% Calculate correlation coefficient and p-value
[correlation_coeff, p_value] = corrcoef(temperatures_clean, resistances_clean);
r_value = correlation_coeff(1,2);
p_value = p_value(1,2);

% Calculate linear regression
p = polyfit(temperatures_clean, resistances_clean, 1);
slope = p(1);
intercept = p(2);

% Calculate R-squared
y_pred = polyval(p, temperatures_clean);
ss_res = sum((resistances_clean - y_pred).^2);
ss_tot = sum((resistances_clean - mean(resistances_clean)).^2);
r_squared = 1 - (ss_res / ss_tot);

fprintf('\n=== Statistical Analysis ===\n');
fprintf('Correlation coefficient (r): %.4f\n', r_value);
fprintf('P-value: %.6f\n', p_value);
fprintf('R-squared: %.4f\n', r_squared);
fprintf('Slope: %.4f mΩ/°C\n', slope);
fprintf('Intercept: %.4f mΩ\n', intercept);

% Determine significance
if p_value < 0.001
    significance = '*** (p < 0.001)';
elseif p_value < 0.01
    significance = '** (p < 0.01)';
elseif p_value < 0.05
    significance = '* (p < 0.05)';
else
    significance = 'ns (p >= 0.05)';
end

fprintf('Significance: %s\n', significance);

%% Visualization
figure('Position', [100, 100, 1200, 800]);

% Main scatter plot
subplot(2,2,1);
scatter(temperatures_clean, resistances_clean, 50, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;

% Add regression line
x_range = linspace(min(temperatures_clean), max(temperatures_clean), 100);
y_regression = polyval(p, x_range);
plot(x_range, y_regression, 'r-', 'LineWidth', 2);

% Add confidence intervals (95%)
% Calculate prediction intervals manually
n = length(temperatures_clean);
x_mean = mean(temperatures_clean);
Sxx = sum((temperatures_clean - x_mean).^2);

% Calculate residuals for confidence interval
residuals_for_ci = resistances_clean - y_pred;

% Standard error of prediction
se_pred = sqrt(1/n + (x_range - x_mean).^2 / Sxx) * sqrt(sum(residuals_for_ci.^2) / (n-2));

% 95% confidence interval (t-distribution)
t_critical = tinv(0.975, n-2);
ci_upper = y_regression + t_critical * se_pred;
ci_lower = y_regression - t_critical * se_pred;

plot(x_range, ci_upper, 'r--', 'LineWidth', 1);
plot(x_range, ci_lower, 'r--', 'LineWidth', 1);

xlabel('Ambient Temperature [°C]', 'FontSize', 12);
ylabel('DCIR [mΩ]', 'FontSize', 12);
title(sprintf('Temperature vs DCIR Correlation\nr = %.3f, p = %.4f %s', r_value, p_value, significance), 'FontSize', 14);
grid on;
legend('Data Points', 'Regression Line', '95% CI', 'Location', 'best');

% Residual plot
subplot(2,2,2);
residuals = resistances_clean - y_pred;
scatter(temperatures_clean, residuals, 50, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
yline(0, 'r-', 'LineWidth', 2);
xlabel('Ambient Temperature [°C]', 'FontSize', 12);
ylabel('Residuals [mΩ]', 'FontSize', 12);
title('Residual Plot', 'FontSize', 14);
grid on;

% Histogram of residuals
subplot(2,2,3);
histogram(residuals, 20, 'FaceColor', '#0073C2', 'EdgeColor', 'black');
xlabel('Residuals [mΩ]', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Residual Distribution', 'FontSize', 14);
grid on;

% Q-Q plot for normality test
subplot(2,2,4);
qqplot(residuals);
title('Q-Q Plot for Normality Test', 'FontSize', 14);
xlabel('Theoretical Quantiles', 'FontSize', 12);
ylabel('Sample Quantiles', 'FontSize', 12);

% Add statistics text box
stats_text = sprintf(['Statistics:\n' ...
    'Sample Size: %d\n' ...
    'Correlation (r): %.3f\n' ...
    'P-value: %.4f\n' ...
    'R-squared: %.3f\n' ...
    'Slope: %.3f mΩ/°C\n' ...
    'Intercept: %.2f mΩ\n' ...
    'Significance: %s'], ...
    length(temperatures_clean), r_value, p_value, r_squared, slope, intercept, significance);

annotation('textbox', [0.02, 0.02, 0.3, 0.15], 'String', stats_text, ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');

%% Monthly Analysis (Option 2: Monthly Averages)
% Group data by month and calculate monthly averages
if isfield(Peaks_all, 'Month')
    months = Peaks_all.Month;
    % Apply the same filtering to months as we did to temperatures and resistances
    months_clean = months(valid_indices);
    unique_months = unique(months_clean);
    
    fprintf('\n=== Monthly Analysis (Monthly Averages) ===\n');
    
    % Calculate monthly averages
    monthly_avg_temp = [];
    monthly_avg_res = [];
    monthly_std_temp = [];
    monthly_std_res = [];
    monthly_counts = [];
    
    for i = 1:length(unique_months)
        month_idx = months_clean == unique_months(i);
        month_temp = temperatures_clean(month_idx);
        month_res = resistances_clean(month_idx);
        
        if length(month_temp) > 0  % Only analyze months with data
            monthly_avg_temp(i) = mean(month_temp);
            monthly_avg_res(i) = mean(month_res);
            monthly_std_temp(i) = std(month_temp);
            monthly_std_res(i) = std(month_res);
            monthly_counts(i) = length(month_temp);
            
            fprintf('Month %d: Avg Temp = %.2f±%.2f°C, Avg Res = %.3f±%.3f mΩ, n = %d\n', ...
                unique_months(i), monthly_avg_temp(i), monthly_std_temp(i), ...
                monthly_avg_res(i), monthly_std_res(i), monthly_counts(i));
        else
            monthly_avg_temp(i) = NaN;
            monthly_avg_res(i) = NaN;
            monthly_std_temp(i) = NaN;
            monthly_std_res(i) = NaN;
            monthly_counts(i) = 0;
            fprintf('Month %d: No data\n', unique_months(i));
        end
    end
    
    % Calculate correlation using monthly averages
    valid_monthly_idx = ~isnan(monthly_avg_temp) & ~isnan(monthly_avg_res);
    if sum(valid_monthly_idx) > 2  % Need at least 3 months for correlation
        [monthly_corr, monthly_p] = corrcoef(monthly_avg_temp(valid_monthly_idx), monthly_avg_res(valid_monthly_idx));
        monthly_r_value = monthly_corr(1,2);
        monthly_p_value = monthly_p(1,2);
        
        fprintf('\nMonthly Average Correlation:\n');
        fprintf('Correlation coefficient (r): %.4f\n', monthly_r_value);
        fprintf('P-value: %.6f\n', monthly_p_value);
        
        % Calculate linear regression for monthly averages
        monthly_p_coeff = polyfit(monthly_avg_temp(valid_monthly_idx), monthly_avg_res(valid_monthly_idx), 1);
        monthly_slope = monthly_p_coeff(1);
        monthly_intercept = monthly_p_coeff(2);
        
        % Calculate R-squared for monthly averages
        monthly_y_pred = polyval(monthly_p_coeff, monthly_avg_temp(valid_monthly_idx));
        monthly_ss_res = sum((monthly_avg_res(valid_monthly_idx) - monthly_y_pred).^2);
        monthly_ss_tot = sum((monthly_avg_res(valid_monthly_idx) - mean(monthly_avg_res(valid_monthly_idx))).^2);
        monthly_r_squared = 1 - (monthly_ss_res / monthly_ss_tot);
        
        fprintf('R-squared: %.4f\n', monthly_r_squared);
        fprintf('Slope: %.4f mΩ/°C\n', monthly_slope);
        fprintf('Intercept: %.4f mΩ\n', monthly_intercept);
        
        % Determine significance for monthly averages
        if monthly_p_value < 0.001
            monthly_significance = '*** (p < 0.001)';
        elseif monthly_p_value < 0.01
            monthly_significance = '** (p < 0.01)';
        elseif monthly_p_value < 0.05
            monthly_significance = '* (p < 0.05)';
        else
            monthly_significance = 'ns (p >= 0.05)';
        end
        fprintf('Significance: %s\n', monthly_significance);
    else
        monthly_r_value = NaN;
        monthly_p_value = NaN;
        monthly_r_squared = NaN;
        monthly_slope = NaN;
        monthly_intercept = NaN;
        monthly_significance = 'Insufficient data';
        fprintf('\nInsufficient monthly data for correlation analysis\n');
    end
    
    % Visualization for monthly averages
    figure('Position', [200, 200, 1200, 800]);
    
    % Monthly averages scatter plot
    subplot(2,2,1);
    scatter(monthly_avg_temp(valid_monthly_idx), monthly_avg_res(valid_monthly_idx), 100, 'filled', 'MarkerFaceAlpha', 0.8);
    hold on;
    
    if sum(valid_monthly_idx) > 2
        % Add regression line for monthly averages
        x_range_monthly = linspace(min(monthly_avg_temp(valid_monthly_idx)), max(monthly_avg_temp(valid_monthly_idx)), 100);
        y_regression_monthly = polyval(monthly_p_coeff, x_range_monthly);
        plot(x_range_monthly, y_regression_monthly, 'r-', 'LineWidth', 2);
    end
    
    % Add error bars (standard deviation)
    errorbar(monthly_avg_temp(valid_monthly_idx), monthly_avg_res(valid_monthly_idx), ...
        monthly_std_res(valid_monthly_idx), monthly_std_res(valid_monthly_idx), ...
        monthly_std_temp(valid_monthly_idx), monthly_std_temp(valid_monthly_idx), 'o', 'LineWidth', 1);
    
    xlabel('Monthly Average Temperature [°C]', 'FontSize', 12);
    ylabel('Monthly Average DCIR [mΩ]', 'FontSize', 12);
    title(sprintf('Monthly Averages: Temperature vs DCIR\nr = %.3f, p = %.4f %s', ...
        monthly_r_value, monthly_p_value, monthly_significance), 'FontSize', 14);
    grid on;
    legend('Monthly Averages', 'Regression Line', 'Error Bars (±1σ)', 'Location', 'best');
    
    % Monthly trends
    subplot(2,2,2);
    plot(unique_months(valid_monthly_idx), monthly_avg_temp(valid_monthly_idx), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Month', 'FontSize', 12);
    ylabel('Average Temperature [°C]', 'FontSize', 12);
    title('Monthly Temperature Trends', 'FontSize', 14);
    grid on;
    
    subplot(2,2,3);
    plot(unique_months(valid_monthly_idx), monthly_avg_res(valid_monthly_idx), 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Month', 'FontSize', 12);
    ylabel('Average DCIR [mΩ]', 'FontSize', 12);
    title('Monthly DCIR Trends', 'FontSize', 14);
    grid on;
    
    % Monthly sample sizes
    subplot(2,2,4);
    bar(unique_months(valid_monthly_idx), monthly_counts(valid_monthly_idx));
    xlabel('Month', 'FontSize', 12);
    ylabel('Number of Data Points', 'FontSize', 12);
    title('Monthly Sample Sizes', 'FontSize', 14);
    grid on;
    
    % Add statistics text box for monthly analysis
    if sum(valid_monthly_idx) > 2
        monthly_stats_text = sprintf(['Monthly Analysis Statistics:\n' ...
            'Number of Months: %d\n' ...
            'Correlation (r): %.3f\n' ...
            'P-value: %.4f\n' ...
            'R-squared: %.3f\n' ...
            'Slope: %.3f mΩ/°C\n' ...
            'Intercept: %.2f mΩ\n' ...
            'Significance: %s'], ...
            sum(valid_monthly_idx), monthly_r_value, monthly_p_value, monthly_r_squared, ...
            monthly_slope, monthly_intercept, monthly_significance);
        
        annotation('textbox', [0.02, 0.02, 0.3, 0.15], 'String', monthly_stats_text, ...
            'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
end

%% Save results
saveDir = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\FieldData\FieldData_Rack_DCIR\DCIR_Charge_Onori_Fig4_5';
save(fullfile(saveDir, 'Temperature_Correlation_Analysis.mat'), ...
    'temperatures_clean', 'resistances_clean', 'r_value', 'p_value', ...
    'r_squared', 'slope', 'intercept', 'significance', 'residuals', ...
    'monthly_avg_temp', 'monthly_avg_res', 'monthly_std_temp', 'monthly_std_res', ...
    'monthly_counts', 'monthly_r_value', 'monthly_p_value', 'monthly_r_squared', ...
    'monthly_slope', 'monthly_intercept', 'monthly_significance');

% Save figures
saveas(gcf, fullfile(saveDir, 'Temperature_Correlation_Analysis.fig'));
saveas(gcf, fullfile(saveDir, 'Temperature_Correlation_Analysis.png'));

fprintf('\nAnalysis complete! Results saved to: %s\n', saveDir); 