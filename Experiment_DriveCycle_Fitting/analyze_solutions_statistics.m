clear; clc; close all;

%% Configuration
target_channel = 'ch9';
target_soc = 'SOC90';
target_dc = 'DC5';
target_cycle = 'cyc200';

%% Load specific result file
base_path = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Experiment_DriveCycle_Fitting\Results';
result_file = fullfile(base_path, target_cycle, target_channel, target_soc, target_dc, ...
    sprintf('%s_%s_%s_result.mat', target_cycle, target_soc, target_dc));

if exist(result_file, 'file')
    fprintf('Loading results from: %s\n', result_file);
    load(result_file);
    
    %% Extract all solutions
    fprintf('\n=== MultiStart Solutions Analysis ===\n');
    fprintf('Condition: %s %s %s %s\n', target_cycle, target_channel, target_soc, target_dc);
    fprintf('Total solutions: %d\n', length(solutions));
    
    % Extract parameters from all solutions
    all_params = zeros(length(solutions), 5);
    all_costs = zeros(length(solutions), 1);
    all_exitflags = zeros(length(solutions), 1);
    
    for i = 1:length(solutions)
        all_params(i, :) = solutions(i).X;
        all_costs(i) = solutions(i).Fval;
        all_exitflags(i) = solutions(i).Exitflag;
    end
    
    %% Calculate statistics
    param_names = {'R0', 'R1', 'R2', 'tau1', 'tau2'};
    param_units = {'Ω', 'Ω', 'Ω', 's', 's'};
    param_units_display = {'mΩ', 'mΩ', 'mΩ', 's', 's'};
    
    fprintf('\n=== Parameter Statistics ===\n');
    for p = 1:5
        values = all_params(:, p);
        
        % Convert to appropriate units for display
        if p <= 3  % R0, R1, R2 -> mΩ
            values_display = values * 1000;
            unit = param_units_display{p};
        else  % tau1, tau2 -> s
            values_display = values;
            unit = param_units_display{p};
        end
        
        avg_val = mean(values_display);
        min_val = min(values_display);
        max_val = max(values_display);
        std_val = std(values_display);
        
        fprintf('%s: %.6f ± %.6f %s (%.6f - %.6f %s)\n', ...
            param_names{p}, avg_val, std_val, unit, min_val, max_val, unit);
    end
    
    %% Cost function statistics
    fprintf('\n=== Cost Function Statistics ===\n');
    cost_avg = mean(all_costs);
    cost_min = min(all_costs);
    cost_max = max(all_costs);
    cost_std = std(all_costs);
    
    fprintf('RMSE: %.6f ± %.6f V (%.6f - %.6f V)\n', ...
        cost_avg, cost_std, cost_min, cost_max);
    fprintf('RMSE: %.3f ± %.3f mV (%.3f - %.3f mV)\n', ...
        cost_avg*1000, cost_std*1000, cost_min*1000, cost_max*1000);
    
    %% Convergence statistics
    fprintf('\n=== Convergence Statistics ===\n');
    converged_solutions = sum(all_exitflags > 0);
    convergence_rate = converged_solutions / length(solutions) * 100;
    
    fprintf('Converged solutions: %d/%d (%.1f%%)\n', ...
        converged_solutions, length(solutions), convergence_rate);
    
    %% Best solution details
    [best_cost, best_idx] = min(all_costs);
    best_params = all_params(best_idx, :);
    
    fprintf('\n=== Best Solution ===\n');
    fprintf('Solution index: %d\n', best_idx);
    fprintf('RMSE: %.6f V (%.3f mV)\n', best_cost, best_cost*1000);
    fprintf('Exit flag: %d\n', all_exitflags(best_idx));
    fprintf('Parameters:\n');
    for p = 1:5
        if p <= 3
            fprintf('  %s: %.6f Ω (%.3f mΩ)\n', param_names{p}, best_params(p), best_params(p)*1000);
        else
            fprintf('  %s: %.6f s\n', param_names{p}, best_params(p));
        end
    end
    
    %% Capacitance calculation
    C1_best = (best_params(4) / best_params(2)) / 1000;  % kF
    C2_best = (best_params(5) / best_params(3)) / 1000;  % kF
    
    fprintf('Capacitances:\n');
    fprintf('  C1: %.2f kF\n', C1_best);
    fprintf('  C2: %.2f kF\n', C2_best);
    
    %% Save statistics
    statistics = struct();
    statistics.condition = sprintf('%s_%s_%s_%s', target_cycle, target_channel, target_soc, target_dc);
    statistics.total_solutions = length(solutions);
    statistics.converged_solutions = converged_solutions;
    statistics.convergence_rate = convergence_rate;
    statistics.best_solution_index = best_idx;
    statistics.best_cost = best_cost;
    statistics.best_params = best_params;
    statistics.all_params = all_params;
    statistics.all_costs = all_costs;
    statistics.all_exitflags = all_exitflags;
    
    % Save to file
    output_file = fullfile(base_path, target_cycle, target_channel, target_soc, target_dc, ...
        sprintf('%s_%s_%s_%s_solutions_statistics.mat', target_cycle, target_channel, target_soc, target_dc));
    save(output_file, 'statistics');
    fprintf('\nStatistics saved to: %s\n', output_file);
    
else
    fprintf('Error: Result file not found: %s\n', result_file);
end 