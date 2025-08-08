clear; clc; close all;

%% Configuration
base_path = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Experiment_DriveCycle_Fitting\Results';
channels = {'ch9'};
cycle_conditions = {'cyc0', 'cyc200'};
soc_conditions = {'SOC50'}; % , 'SOC70', 'SOC90'};
dc_name = 'DC1';

%% Standard PSD Analysis Function
function analyze_standard_psd(base_path, current_cycle, current_channel, current_soc, dc_name)
    % Construct file path for summary table
    summary_file = fullfile(base_path, current_cycle, current_channel, current_soc, ...
        sprintf('%s_%s_summary_table.mat', current_cycle, current_soc));
    
    if ~exist(summary_file, 'file')
        fprintf('Summary file not found: %s\n', summary_file);
        return;
    end
    
    fprintf('\n=== Standard PSD Analysis: %s %s %s %s ===\n', current_cycle, current_channel, current_soc, dc_name);
    
    % Load summary file
    load(summary_file);
    
    % Extract data from soc_results structure
    if exist('soc_results', 'var') && isfield(soc_results, dc_name)
        result = soc_results.(dc_name);
        measured_voltage = result.V_measured;
        predicted_voltage = result.V_model;
        dt = 0.1;  % Fixed time step
        
        % Extract RC parameters for reference
        tau1 = result.tau1;
        tau2 = result.tau2;
        R1 = result.R1;
        R2 = result.R2;
        C1 = result.C1;
        C2 = result.C2;
    else
        fprintf('Voltage data not found in summary file for %s\n', dc_name);
        return;
    end
    
    % Calculate model error
    voltage_error = measured_voltage - predicted_voltage;
    
    % Calculate time vector
    N = length(measured_voltage);
    time = (0:N-1) * dt;
    
    % Calculate PSD using standard method
    fs = 1/dt;  % Sampling frequency
    window_length = min(256, floor(N/4));
    window = hamming(window_length);
    [pxx, f] = pwelch(voltage_error, window, [], [], fs);
    
    % Standard frequency bands for battery analysis (simplified)
    lowest_mask = f < 0.1;           % < 0.1 Hz (SOC drift)
    low_mask = f >= 0.1 & f < 1;    % 0.1-1 Hz (OCV/RC network)
    mid_mask = f >= 1 & f < 5;      % 1-5 Hz (fast dynamics)
    high_mask = f >= 5;              % >= 5 Hz (noise)
    
    % Calculate power in each band
    
    if any(lowest_mask)
        f_lowest = f(lowest_mask);
        pxx_lowest = pxx(lowest_mask);
        if length(f_lowest) > 1
            lowest_power = trapz(f_lowest, pxx_lowest);
        else
            lowest_power = 0;
        end
    else
        lowest_power = 0;
    end
        
    if any(low_mask)
        f_low = f(low_mask);
        pxx_low = pxx(low_mask);
        if length(f_low) > 1
            low_power = trapz(f_low, pxx_low);
        else
            low_power = 0;
        end
    else
        low_power = 0;
    end
    
    if any(mid_mask)
        f_mid = f(mid_mask);
        pxx_mid = pxx(mid_mask);
        if length(f_mid) > 1
            mid_power = trapz(f_mid, pxx_mid);
        else
            mid_power = 0;
        end
    else
        mid_power = 0;
    end
    
    if any(high_mask)
        f_high = f(high_mask);
        pxx_high = pxx(high_mask);
        if length(f_high) > 1
            high_power = trapz(f_high, pxx_high);
        else
            high_power = 0;
        end
    else
        high_power = 0;
    end
    
    total_power = lowest_power + low_power + mid_power + high_power;
    
    % Display results
    fprintf('\nModel Error Statistics:\n');
    fprintf('  RMSE: %.6f V\n', sqrt(mean(voltage_error.^2)));
    fprintf('  MAE:  %.6f V\n', mean(abs(voltage_error)));
    fprintf('  Max Error: %.6f V\n', max(abs(voltage_error)));
    fprintf('  Std Error: %.6f V\n', std(voltage_error));
    
    fprintf('\nRC Network Parameters:\n');
    fprintf('  τ1 = %.2f s, τ2 = %.2f s\n', tau1, tau2);
    fprintf('  C1 = %.2f kF, C2 = %.2f kF\n', C1, C2);
    fprintf('  R1 = %.2f mΩ, R2 = %.2f mΩ\n', R1*1000, R2*1000);
    
    fprintf('\nStandard PSD Analysis:\n');
    fprintf('  Lowest frequency (<0.1 Hz): %.6f V² (%.1f%%)\n', lowest_power, lowest_power/total_power*100);
    fprintf('  Low frequency (0.1-1 Hz): %.6f V² (%.1f%%)\n', low_power, low_power/total_power*100);
    fprintf('  Mid frequency (1-5 Hz): %.6f V² (%.1f%%)\n', mid_power, mid_power/total_power*100);
    fprintf('  High frequency (>=5 Hz): %.6f V² (%.1f%%)\n', high_power, high_power/total_power*100);
    fprintf('  Total power: %.6f V²\n', total_power);
    
    % Standard diagnostics
    fprintf('\n=== Standard Diagnostics ===\n');
    
    % Determine dominant frequency component
    [~, max_idx] = max([lowest_power, low_power, mid_power, high_power]);
    issues = {'Lowest frequency issue (SOC drift)', 'Low frequency issue (OCV/RC)', 'Mid frequency issue (fast dynamics)', 'High frequency noise'};
    fprintf('Dominant component: %s\n', issues{max_idx});
    
    % Check if RC parameters are reasonable
    if tau1 > 100
        fprintf('⚠ τ1 = %.2f s is too large (typical: 10-50 s)\n', tau1);
    end
    
    if tau2 > 1000
        fprintf('⚠ τ2 = %.2f s is extremely large (typical: 100-500 s)\n', tau2);
    end
    
    if tau2/tau1 < 5
        fprintf('⚠ τ2/τ1 ratio = %.1f is too small (should be >5)\n', tau2/tau1);
    end
    
    % Create standard visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % Time domain error
    subplot(2,2,1);
    plot(time, voltage_error, 'b-', 'LineWidth', 1);
    title(sprintf('Model Error - %s %s %s %s', current_cycle, current_channel, current_soc, dc_name));
    xlabel('Time (s)'); ylabel('Error (V)');
    grid on;
    
    % PSD plot (standard)
    subplot(2,2,2);
    semilogy(f, pxx, 'r-', 'LineWidth', 2);
    title('Power Spectral Density');
    xlabel('Frequency (Hz)'); ylabel('PSD (V²/Hz)');
    grid on;
    
    % Frequency band comparison
    subplot(2,2,3);
    freq_bands = {'Lowest (<0.1 Hz)', 'Low (0.1-1 Hz)', 'Mid (1-5 Hz)', 'High (>=5 Hz)'};
    powers = [lowest_power, low_power, mid_power, high_power];
    bar(powers);
    title('Power in Frequency Bands');
    xlabel('Frequency Band'); ylabel('Power (V²)');
    set(gca, 'XTickLabel', freq_bands);
    grid on;
    
    % RC parameters
    subplot(2,2,4);
    rc_params = [tau1, tau2; C1, C2; R1*1000, R2*1000];
    bar(rc_params);
    title('RC Network Parameters');
    xlabel('RC Element'); ylabel('Value');
    set(gca, 'XTick', 1:2, 'XTickLabel', {'τ1', 'τ2'});
    legend('τ (s)', 'C (kF)', 'R (mΩ)', 'Location', 'best');
    grid on;
    
    % Save figure
    output_dir = fullfile(base_path, 'Standard_PSD_Analysis');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    fig_file = fullfile(output_dir, sprintf('Standard_PSD_%s_%s_%s_%s.fig', current_cycle, current_channel, current_soc, dc_name));
    savefig(gcf, fig_file);
    fprintf('Standard PSD analysis figure saved: %s\n', fig_file);
    
    % Save analysis results
    analysis_results = struct();
    analysis_results.condition = sprintf('%s_%s_%s_%s', current_cycle, current_channel, current_soc, dc_name);
    analysis_results.tau1 = tau1;
    analysis_results.tau2 = tau2;
    analysis_results.C1 = C1;
    analysis_results.C2 = C2;
    analysis_results.R1 = R1;
    analysis_results.R2 = R2;
    analysis_results.rmse = sqrt(mean(voltage_error.^2));
    analysis_results.mae = mean(abs(voltage_error));
    analysis_results.max_error = max(abs(voltage_error));
    analysis_results.std_error = std(voltage_error);
    analysis_results.lowest_power = lowest_power;    
    analysis_results.low_power = low_power;
    analysis_results.mid_power = mid_power;
    analysis_results.high_power = high_power;
    analysis_results.total_power = total_power;
    analysis_results.dominant_component = issues{max_idx};
    analysis_results.frequency = f;
    analysis_results.psd = pxx;
    analysis_results.voltage_error = voltage_error;
    analysis_results.time = time;
    
    mat_file = fullfile(output_dir, sprintf('Standard_PSD_%s_%s_%s_%s_analysis.mat', current_cycle, current_channel, current_soc, dc_name));
    save(mat_file, 'analysis_results');
    fprintf('Standard PSD analysis data saved: %s\n', mat_file);
    
    close(gcf);
end

%% Main Analysis Loop
fprintf('=== Standard PSD Analysis ===\n');
fprintf('Analyzing ch9 DC1 for all SOC conditions\n');

for current_cycle = cycle_conditions
    current_cycle = current_cycle{1};
    
    for current_soc = soc_conditions
        current_soc = current_soc{1};
        
        for current_channel = channels
            current_channel = current_channel{1};
            
            % Analyze standard PSD for this condition
            analyze_standard_psd(base_path, current_cycle, current_channel, current_soc, dc_name);
        end
    end
end

fprintf('\n=== Standard PSD Analysis Complete ===\n'); 