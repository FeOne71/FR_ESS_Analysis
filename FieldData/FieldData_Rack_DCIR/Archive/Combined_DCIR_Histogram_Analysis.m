%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined DCIR Histogram Analysis
% Combines multiple years of DCIR data into a single histogram
% with normal distribution curves and statistics for each year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory and File Paths
% Load the saved data from the original analysis
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\Charge\Combined_DCIR_Analysis');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Load Data from Multiple Years
% Load the saved data files for each year
% You can modify this section to load data from different years
yearList = {'2023','2025'}; % Add more years as needed

all_dcir_data = struct();
all_years_data = [];

for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    
    % Create valid field name for struct (add 'year_' prefix)
    year_field_name = sprintf('year_%s', year);
    
    % Load the saved data file for this year from the specific year folder
    year_data_dir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_Rack_DCIR\Charge', sprintf('DCIR_Charge_%s_NewLogic', year));
    data_file = fullfile(year_data_dir, sprintf('all_chg_events_onori_newlogic_all_years_%s.mat', year));
    
    if exist(data_file, 'file')
        fprintf('Loading data for year: %s from %s\n', year, data_file);
        load(data_file);
        
        % Extract DCIR values from the loaded data
        if exist('global_eventStruct', 'var')
            % Extract DCIR values from the global_eventStruct
            dcir_values = [];
            
            % Navigate through the structure to find DCIR values
            if isfield(global_eventStruct, 'Rack01')
                rack_data = global_eventStruct.Rack01;
                year_fields = fieldnames(rack_data);
                
                for yf = 1:length(year_fields)
                    year_field = year_fields{yf};
                    if contains(year_field, year)
                        year_data = rack_data.(year_field);
                        date_fields = fieldnames(year_data);
                        
                        for df = 1:length(date_fields)
                            date_field = date_fields{df};
                            date_data = year_data.(date_field);
                            event_fields = fieldnames(date_data);
                            
                            for ef = 1:length(event_fields)
                                event_field = event_fields{ef};
                                event_data = date_data.(event_field);
                                
                                if isfield(event_data, 'PeakChgR')
                                    dcir_val = event_data.PeakChgR;
                                    if ~isnan(dcir_val) && dcir_val > 0
                                        dcir_values = [dcir_values; dcir_val];
                                    end
                                end
                            end
                        end
                    end
                end
                        end
            
            % 1σ 이상치 제거
            if ~isempty(dcir_values)
                mean_val = mean(dcir_values);
                std_val = std(dcir_values);
                lower_bound = mean_val - std_val;
                upper_bound = mean_val + std_val;
                
                % 이상치 제거 전 데이터 수
                original_count = length(dcir_values);
                
                % 1σ 범위 내 데이터만 유지
                dcir_values_filtered = dcir_values(dcir_values >= lower_bound & dcir_values <= upper_bound);
                
                % 이상치 제거 후 데이터 수
                filtered_count = length(dcir_values_filtered);
                removed_count = original_count - filtered_count;
                
                fprintf('Year %s: Original %d values, Removed %d outliers (%.1f%%), Final %d values\n', ...
                    year, original_count, removed_count, (removed_count/original_count)*100, filtered_count);
                
                % 필터링된 데이터로 통계 재계산
                dcir_values = dcir_values_filtered;
            end
            
            % Store the data for this year using valid field name
            all_dcir_data.(year_field_name).values = dcir_values;
            all_dcir_data.(year_field_name).mean = mean(dcir_values);
            all_dcir_data.(year_field_name).std = std(dcir_values);
            all_dcir_data.(year_field_name).count = length(dcir_values);
            all_dcir_data.(year_field_name).year = year; % Store original year string
            
            % Add to combined data
            all_years_data = [all_years_data; dcir_values];
            
            fprintf('Year %s: %d DCIR values, Mean=%.3f mΩ, Std=%.3f mΩ\n', ...
                year, length(dcir_values), mean(dcir_values), std(dcir_values));
        else
            fprintf('Warning: No valid data structure found for year %s\n', year);
        end
    else
        fprintf('Warning: Data file not found for year %s: %s\n', year, data_file);
    end
end

%% Create Combined Histogram
if ~isempty(all_years_data)
    figure('Position', [100, 100, 1200, 800]);
    
    % Set histogram parameters - Adaptive binning for normal-like distribution
    min_val = min(all_years_data);
    max_val = max(all_years_data);
    
    % Freedman-Diaconis rule for optimal bin width
    iqr_val = iqr(all_years_data);
    n_data = length(all_years_data);
    bin_width = 2 * iqr_val / (n_data^(1/3));
    
    % Ensure reasonable number of bins (not too few, not too many)
    n_bins = ceil((max_val - min_val) / bin_width);
    n_bins = max(10, min(50, n_bins)); % 10~50개 bin으로 제한
    
    % Recalculate bin width with adjusted number of bins
    bin_width = (max_val - min_val) / n_bins;
    bin_edges = min_val:bin_width:max_val;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    
    fprintf('Adaptive binning: Range [%.3f, %.3f] mΩ, Bin width: %.3f mΩ, Total bins: %d\n', ...
        min_val, max_val, bin_width, length(bin_edges)-1);
    fprintf('IQR: %.3f mΩ, Data points: %d\n', iqr_val, n_data);
    
    % Create subplot for histogram
    subplot(2, 1, 1);
    hold on;
    
    % Colors for different years (user specified) - ordered by year
    colors_by_year = [32/255, 133/255, 78/255;  % 2023: #20854E (green)
                      0/255, 115/255, 194/255;   % 2024: #0073C2 (blue)
                      146/255, 94/255, 159/255]; % 2025: #925E9F (purple)
    alpha_by_year = [0.4, 0.7, 1.0]; % Transparency: 2023(light) to 2025(solid)
    
    % Dynamic reordering based on number of years
    n_years = length(yearList);
    if n_years == 2
        % For 2 years: newer year first, older year last
        plot_indices = [n_years, n_years-1]; % e.g., [2, 1] for 2025, 2023
        colors = colors_by_year(plot_indices, :);
        alpha_values = alpha_by_year(plot_indices);
    else
        % For 3+ years: newest first, oldest last
        plot_indices = n_years:-1:1; % e.g., [3, 2, 1] for 2025, 2024, 2023
        colors = colors_by_year(plot_indices, :);
        alpha_values = alpha_by_year(plot_indices);
    end
    
    % Plot histogram for each year (newest first, oldest last - newest at bottom)
    if n_years == 2
        plot_order = [n_years, n_years-1]; % e.g., [2, 1] for 2025, 2023
    else
        plot_order = n_years:-1:1; % e.g., [3, 2, 1] for 2025, 2024, 2023
    end
    
    for plot_idx = 1:length(plot_order)
        year_idx = plot_order(plot_idx);
        year = yearList{year_idx};
        year_field_name = sprintf('year_%s', year);
        
        if isfield(all_dcir_data, year_field_name) && ~isempty(all_dcir_data.(year_field_name).values)
            dcir_values = all_dcir_data.(year_field_name).values;
            
            % Create histogram
            [counts, ~] = histcounts(dcir_values, bin_edges);
            
            % Plot histogram bars with year-specific styling
            bar(bin_centers, counts, 0.8, 'FaceColor', colors(plot_idx, :), ...
                'FaceAlpha', alpha_values(plot_idx), 'EdgeColor', colors(plot_idx, :), ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Year %s (n=%d)', year, length(dcir_values)));
            
            % Fit normal distribution to actual data
            try
                % Fit normal distribution to the data
                pd = fitdist(dcir_values, 'Normal');
                
                % Generate points for fitted normal distribution curve
                x_curve = linspace(min(bin_edges), max(bin_edges), 1000);
                y_curve = pdf(pd, x_curve);
                
                % Scale the curve to match histogram height
                scale_factor = max(counts) / max(y_curve);
                y_curve_scaled = y_curve * scale_factor;
                
                % Plot fitted normal distribution curve
                plot(x_curve, y_curve_scaled, 'Color', colors(plot_idx, :), ...
                    'LineWidth', 3, 'DisplayName', sprintf('Fitted Normal (μ=%.3f, σ=%.3f)', pd.mu, pd.sigma));
                
                % Store fitted parameters
                all_dcir_data.(year_field_name).fitted_mu = pd.mu;
                all_dcir_data.(year_field_name).fitted_sigma = pd.sigma;
                
            catch ME
                fprintf('Warning: Could not fit normal distribution for year %s: %s\n', year, ME.message);
                % Fallback to theoretical normal distribution
                mu = all_dcir_data.(year_field_name).mean;
                sigma = all_dcir_data.(year_field_name).std;
                x_curve = linspace(min(bin_edges), max(bin_edges), 1000);
                y_curve = normpdf(x_curve, mu, sigma);
                scale_factor = max(counts) / max(y_curve);
                y_curve_scaled = y_curve * scale_factor;
                plot(x_curve, y_curve_scaled, 'Color', colors(plot_idx, :), ...
                    'LineWidth', 3, 'DisplayName', sprintf('Theoretical Normal (μ=%.3f, σ=%.3f)', mu, sigma));
            end
        end
    end
    
    xlabel('DCIR [mΩ]', 'FontSize', 12);
    ylabel('Frequency', 'FontSize', 12);
    title('Combined DCIR Histogram with Fitted Normal Distribution Curves (1σ Outliers Removed)', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    xlim([0, 1]);
    
    % Add statistics text box
    stats_text = 'Statistics (1σ outliers removed):\n';
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_field_name = sprintf('year_%s', year);
        
        if isfield(all_dcir_data, year_field_name) && ~isempty(all_dcir_data.(year_field_name).values)
            if isfield(all_dcir_data.(year_field_name), 'fitted_mu')
                stats_text = [stats_text, sprintf('Year %s: n=%d, μ=%.3f, σ=%.3f (fitted)\n', ...
                    year, all_dcir_data.(year_field_name).count, all_dcir_data.(year_field_name).fitted_mu, all_dcir_data.(year_field_name).fitted_sigma)];
            else
                stats_text = [stats_text, sprintf('Year %s: n=%d, μ=%.3f, σ=%.3f\n', ...
                    year, all_dcir_data.(year_field_name).count, all_dcir_data.(year_field_name).mean, all_dcir_data.(year_field_name).std)];
            end
        end
    end
    
    text(0.02, 0.98, stats_text, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 10, ...
        'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    %% Create Box Plot for Comparison
    subplot(2, 1, 2);
    hold on;
    
    % Prepare data for box plot
    box_data = [];
    group_labels = {};
    
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_field_name = sprintf('year_%s', year);
        
        if isfield(all_dcir_data, year_field_name) && ~isempty(all_dcir_data.(year_field_name).values)
            box_data = [box_data; all_dcir_data.(year_field_name).values];
            group_labels = [group_labels; repmat({year}, length(all_dcir_data.(year_field_name).values), 1)];
        end
    end
    
    if ~isempty(box_data)
        % Create box plot
        boxplot(box_data, group_labels, 'Colors', colors(1:n_years, :), ...
            'Width', 0.7, 'Labels', yearList);
        
        % Add individual data points
        for year_idx = 1:n_years
            year = yearList{year_idx};
            year_field_name = sprintf('year_%s', year);
            
            if isfield(all_dcir_data, year_field_name) && ~isempty(all_dcir_data.(year_field_name).values)
                x_pos = year_idx + (rand(length(all_dcir_data.(year_field_name).values), 1) - 0.5) * 0.3;
                scatter(x_pos, all_dcir_data.(year_field_name).values, 20, colors(year_idx, :), ...
                    'filled', 'MarkerFaceAlpha', alpha_values(year_idx));
            end
        end
        
        ylabel('DCIR [mΩ]', 'FontSize', 12);
        title('DCIR Distribution Comparison by Year', 'FontSize', 14);
        grid on;
        ylim([min_val, max_val]);
    end
    
    % Save the combined figure
    saveas(gcf, fullfile(saveDir, 'Combined_DCIR_Histogram_Analysis.fig'));
    saveas(gcf, fullfile(saveDir, 'Combined_DCIR_Histogram_Analysis.png'));
    
    fprintf('Combined histogram saved to: %s\n', saveDir);
    
    %% Print Summary Statistics
    fprintf('\n=== Summary Statistics ===\n');
    fprintf('Total DCIR values across all years: %d\n', length(all_years_data));
    fprintf('Overall mean: %.3f mΩ\n', mean(all_years_data));
    fprintf('Overall std: %.3f mΩ\n', std(all_years_data));
    fprintf('Overall median: %.3f mΩ\n', median(all_years_data));
    
    % Save combined data
    save(fullfile(saveDir, 'combined_dcir_data.mat'), 'all_dcir_data', 'all_years_data', 'yearList');
    
else
    fprintf('No DCIR data found. Please check the data files.\n');
end

fprintf('\nAnalysis complete!\n');
