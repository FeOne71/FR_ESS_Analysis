%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined DCIR Histogram Analysis - Onori Method
% Loads existing .mat files and creates combined histogram with normal distribution
% Supports multiple years with different colors and transparency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory and file setup
currentDir = pwd;
fprintf('Current directory: %s\n', currentDir);

% ===== 연도 설정 (여기서 원하는 연도를 직접 입력하세요) =====
% 분석하고 싶은 연도들을 여기에 입력하세요
years_to_analyze = {'2021', '2022'}; % 예: {'2021', '2023'} 또는 {'2021', '2022', '2023'}
% ===============================================================

fprintf('Analyzing years: %s\n', strjoin(years_to_analyze, ', '));

% Find all .mat files in current directory
matFiles = dir(fullfile(currentDir, 'all_chg_events_onori_newlogic_all_years*.mat'));
fprintf('Found %d .mat files:\n', length(matFiles));
for i = 1:length(matFiles)
    fprintf('  %s\n', matFiles(i).name);
end

if isempty(matFiles)
    error('No .mat files found in current directory');
end

%% Load and process data from each .mat file
all_dcir_data = struct();
all_years_data = [];

for file_idx = 1:length(matFiles)
    filename = matFiles(file_idx).name;
    fprintf('\nLoading data from: %s\n', filename);
    
    % Load the .mat file
    loaded_data = load(fullfile(currentDir, filename));
    
    % Check if global_eventStruct exists
    if ~isfield(loaded_data, 'global_eventStruct')
        fprintf('Warning: global_eventStruct not found in %s, skipping...\n', filename);
        continue;
    end
    
    global_eventStruct = loaded_data.global_eventStruct;
    
    % Use the years specified by user
    years = years_to_analyze;
    
    fprintf('Processing years: %s\n', strjoin(years, ', '));
    
    % Process each year
    for year_idx = 1:length(years)
        year = years{year_idx};
        year_field_name = sprintf('year_%s', year);
        
        % Extract DCIR values for this year
        dcir_values = [];
        
        if isfield(global_eventStruct, 'Rack01') && isfield(global_eventStruct.Rack01, year_field_name)
            year_struct = global_eventStruct.Rack01.(year_field_name);
            date_names = fieldnames(year_struct);
            
            fprintf('  Year %s: Found %d dates\n', year, length(date_names));
            
            for d = 1:length(date_names)
                date_struct = year_struct.(date_names{d});
                event_names = fieldnames(date_struct);
                
                for e = 1:length(event_names)
                    event = date_struct.(event_names{e});
                    if isfield(event, 'PeakChgR') && ~isnan(event.PeakChgR)
                        dcir_values = [dcir_values; event.PeakChgR];
                    end
                end
            end
        end
        
        % 모든 데이터 사용 (이상치 제거하지 않음)
        if ~isempty(dcir_values)
            fprintf('  Year %s: Using all %d values (no outlier removal)\n', ...
                year, length(dcir_values));
        end
        
        % Store the data for this year using valid field name
        if ~isempty(dcir_values)
            all_dcir_data.(year_field_name).values = dcir_values;
            all_dcir_data.(year_field_name).count = length(dcir_values);
            all_dcir_data.(year_field_name).mean = mean(dcir_values);
            all_dcir_data.(year_field_name).std = std(dcir_values);
            all_dcir_data.(year_field_name).median = median(dcir_values);
            all_dcir_data.(year_field_name).min = min(dcir_values);
            all_dcir_data.(year_field_name).max = max(dcir_values);
            
            % Add to combined data
            all_years_data = [all_years_data; dcir_values];
        end
    end
end

%% Create combined histogram
if ~isempty(all_years_data)
    figure(1);
    hold on;
    grid on;
    box on;
    
    % Set histogram parameters - Fixed 30 bins
    min_val = min(all_years_data);
    max_val = max(all_years_data);
    
    % Fixed 30 bins
    n_bins = 30;
    bin_width = (max_val - min_val) / n_bins;
    
    % Create bin edges
    bin_edges = min_val:bin_width:max_val;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    
    fprintf('\nFixed binning: Range [%.3f, %.3f] mΩ, Bin width: %.3f mΩ, Total bins: %d\n', ...
        min_val, max_val, bin_width, n_bins);
    fprintf('Data points: %d\n', length(all_years_data));
    
    % Define colors and plotting order (newest year first, oldest last)
    colors_by_year = [32/255, 133/255, 78/255;  % 2021: #20854E (green)
                      0/255, 115/255, 194/255;   % 2022: #0073C2 (blue)
                      146/255, 94/255, 159/255]; % 2023: #925E9F (purple)
    alpha_by_year = [0.4, 0.7, 1.0]; % Transparency: 2021(light) to 2023(solid)
    
    % Get available years
    year_fields = fieldnames(all_dcir_data);
    yearList = {};
    for i = 1:length(year_fields)
        year_str = year_fields{i};
        if startsWith(year_str, 'year_')
            yearList{end+1} = year_str(6:end); % Remove 'year_' prefix
        end
    end
    
    n_years = length(yearList);
    if n_years == 2
        plot_order = [n_years, n_years-1]; % e.g., [2, 1] for 2023, 2021
    else
        plot_order = n_years:-1:1; % e.g., [3, 2, 1] for 2023, 2022, 2021
    end
    
    % Adjust colors and alpha arrays based on available years
    if n_years <= 3
        colors = colors_by_year(plot_order, :);
        alpha_values = alpha_by_year(plot_order);
    else
        % Generate colors for more years if needed
        colors = lines(n_years);
        alpha_values = linspace(0.4, 1.0, n_years);
    end
    
    % Plot histograms for each year
    for plot_idx = 1:length(plot_order)
        year_idx = plot_order(plot_idx);
        year = yearList{year_idx};
        year_field_name = sprintf('year_%s', year);
        
        if isfield(all_dcir_data, year_field_name) && ~isempty(all_dcir_data.(year_field_name).values)
            dcir_values = all_dcir_data.(year_field_name).values;
            
            % Create histogram for this year
            [counts, ~] = histcounts(dcir_values, bin_edges);
            
            % Plot histogram bars
            bar(bin_centers, counts, 1, 'FaceColor', colors(plot_idx, :), ...
                'FaceAlpha', alpha_values(plot_idx), 'EdgeColor', 'none', ...
                'LineWidth', 2, 'DisplayName', sprintf('Year %s (n=%d)', year, length(dcir_values)));
            
            % Fit normal distribution to actual data
            try
                pd = fitdist(dcir_values, 'Normal');
                x_curve = linspace(min(bin_edges), max(bin_edges), 1000);
                y_curve = pdf(pd, x_curve);
                scale_factor = max(counts) / max(y_curve);
                y_curve_scaled = y_curve * scale_factor;
                plot(x_curve, y_curve_scaled, 'Color', colors(plot_idx, :), 'LineWidth', 3, ...
                    'DisplayName', sprintf('Fitted Normal (μ=%.3f, σ=%.3f)', pd.mu, pd.sigma));
                all_dcir_data.(year_field_name).fitted_mu = pd.mu;
                all_dcir_data.(year_field_name).fitted_sigma = pd.sigma;
            catch ME
                % Fallback to theoretical normal distribution if fit fails
                mu = all_dcir_data.(year_field_name).mean;
                sigma = all_dcir_data.(year_field_name).std;
                x_curve = linspace(min(bin_edges), max(bin_edges), 1000);
                y_curve = normpdf(x_curve, mu, sigma);
                scale_factor = max(counts) / max(y_curve);
                y_curve_scaled = y_curve * scale_factor;
                plot(x_curve, y_curve_scaled, 'Color', colors(plot_idx, :), 'LineWidth', 3, ...
                    'DisplayName', sprintf('Theoretical Normal (μ=%.3f, σ=%.3f)', mu, sigma));
            end
        end
    end
    
    xlabel('DCIR [mΩ]', 'interpreter', 'tex');
    ylabel('Frequency', 'interpreter', 'tex');
    title('Combined DCIR Histogram with Fitted Normal Distribution Curves (All Data)', 'FontSize', 14);
    legend('Location', 'best');
    set(gca, 'fontsize', 12);
    set(gca, 'ticklabelinterpreter', 'tex');
    
    % Set axis limits based on data range
    xlim([min_val, max_val]);
    
    % Add statistics text box
    stats_text = {'Statistics (All data):'};
    for year_idx = 1:length(yearList)
        year = yearList{year_idx};
        year_field_name = sprintf('year_%s', year);
        if isfield(all_dcir_data, year_field_name) && ~isempty(all_dcir_data.(year_field_name).values)
            stats_text{end+1} = sprintf('Year %s: μ=%.3f, σ=%.3f, n=%d', ...
                year, all_dcir_data.(year_field_name).mean, ...
                all_dcir_data.(year_field_name).std, ...
                all_dcir_data.(year_field_name).count);
        end
    end
    
    text(0.02, 0.98, stats_text, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 10, ...
        'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Save the figure
    fig_filename = fullfile(currentDir, 'Combined_DCIR_Histogram_Onori.fig');
    saveas(gcf, fig_filename);
    fprintf('\nCombined histogram saved to: %s\n', fig_filename);
    
    fprintf('\nCombined DCIR Histogram complete: %d total valid DCIR values\n', length(all_years_data));
else
    fprintf('\nNo DCIR values found in any of the .mat files.\n');
end

fprintf('\nProcessing complete!\n');
