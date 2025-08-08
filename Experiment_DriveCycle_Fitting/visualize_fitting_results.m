function h_fig = visualize_fitting_results(DC_name, t_measured, I_measured, SOC, V_measured, V_model, V_OCV, terminal_error, terminal_rmse, terminal_mae, color_I, color_V, color_soc)
    % Visualize fitting results for ECM parameter optimization
    % Input:
    %   DC_name: Drive cycle name (e.g., 'DC1')
    %   t_measured: Time vector [s]
    %   I_measured: Current vector [A]
    %   SOC: State of charge vector [%]
    %   V_measured: Measured voltage [V]
    %   V_model: Model voltage [V]
    %   V_OCV: Open circuit voltage [V] (NEW)
    %   terminal_error: Voltage error [V]
    %   terminal_rmse: RMSE [V]
    %   terminal_mae: MAE [V]
    %   color_I, color_V, color_soc: Color codes for plots
    % Output:
    %   h_fig: Figure handle
    
    fprintf('Creating fitting results visualization with OCV for %s...\n', DC_name);
    
    h_fig = figure('Name', sprintf('%s - Terminal Voltage Fitting Results with OCV', DC_name), ...
                   'Position', [50, 50, 1200, 900]);
    
    % Combined Current and SOC profile
    subplot(3, 1, 1);
    yyaxis left
    plot(t_measured/60, I_measured, 'Color', color_I, 'LineWidth', 1.5);
    ylabel('Current (A)', 'Color', color_I);
    ylim([min(I_measured)*1.1, max(I_measured)*1.1]);
    
    yyaxis right
    plot(t_measured/60, SOC, 'Color', color_soc, 'LineWidth', 1.5);
    ylabel('SOC (%)', 'Color', color_soc);
    ylim([min(SOC)-2, max(SOC)+2]);
    
    xlabel('Time (min)');
    title(sprintf('%s: Current & SOC Profile', DC_name));
    grid on;
        
    % Terminal voltage comparison WITH OCV
    subplot(3, 1, 2);
    plot(t_measured/60, V_measured, 'Color', color_V, 'LineWidth', 1.5, 'DisplayName', 'Measured'); hold on;
    plot(t_measured/60, V_model, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Model');
    % plot(t_measured/60, V_OCV, 'Color', color_soc, 'LineWidth', 1.5, 'DisplayName', 'OCV');
    plot(t_measured/60, V_OCV, 'Color', color_soc, 'LineWidth', 1.5, 'LineStyle', '-', 'DisplayName', 'OCV');
    xlabel('Time (min)');
    ylabel('Voltage (V)');
    title(sprintf('%s: Terminal Voltage Fitting', DC_name));
    legend('Location', 'best');
    grid on;
    
    % Terminal voltage fitting errors 
    subplot(3, 1, 3);
    plot(t_measured/60, terminal_error*1000, 'Color', '#CD534C', 'LineWidth', 1.2);
    xlabel('Time (min)');
    ylabel('Error (mV)');
    title(sprintf('%s: Fitting Errors', DC_name));
    text(0.02, 0.98, sprintf('RMSE = %.3f mV\nMAE = %.3f mV', ...
        terminal_rmse*1000, terminal_mae*1000), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black');
    grid on;
    
    fprintf('Fitting results visualization with OCV completed for %s\n', DC_name);
end 