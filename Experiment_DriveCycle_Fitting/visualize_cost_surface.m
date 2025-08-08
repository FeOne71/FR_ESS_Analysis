function visualize_cost_surface(cost_func, options, tau1_vec, tau2_vec)
    % Visualize cost surface for parameter optimization (BSG_Reference 방식)
    % Input:
    %   cost_func: Cost function handle  
    %   optimal_params: Optimal parameter values [R0, R1, R2, tau1, tau2]
    %   param_names: Parameter names cell array
    %   options: fmincon optimization options (from main optimization)
    %   tau1_vec, tau2_vec: 그리드 벡터 (최적값 중심)
    % Features: tau1 < tau2 constraint applied
    
    fprintf('Creating cost surface visualization (tau1 vs tau2) with BSG_Reference constraints...\n');
    
    % Create figure for cost surface visualization
    figure('Position', [100, 100, 800, 600]);
    
    % 입력받은 tau1_vec, tau2_vec 사용
    fprintf('tau1 grid: %.2f ~ %.2f (%d points)\n', min(tau1_vec), max(tau1_vec), length(tau1_vec));
    fprintf('tau2 grid: %.2f ~ %.2f (%d points)\n', min(tau2_vec), max(tau2_vec), length(tau2_vec));
    fprintf('Grid size: %d x %d = %d evaluations\n', length(tau1_vec), length(tau2_vec), length(tau1_vec)*length(tau2_vec));
    
    options_surface = optimset(options, 'Display', 'off');
    
    [T1_grid, T2_grid] = meshgrid(tau1_vec, tau2_vec);
    tau1_flat = T1_grid(:);
    tau2_flat = T2_grid(:);
    cost_flat = zeros(size(tau1_flat));
    
    % 병렬 최적화 루프
    parfor k = 1:length(tau1_flat)
        tau1_k = tau1_flat(k);
        tau2_k = tau2_flat(k);
        p0_R = [0.0012, 0.0006, 0.0004, tau1_k, tau2_k];
        lb_R = [0, 0, 0, tau1_k, tau2_k];
        ub_R = [p0_R(1)*10, p0_R(2)*10, p0_R(3)*10, tau1_k, tau2_k];
        A_constraint = [0, 0, 0, 1, -1];
        b_constraint = 0;
        [~, cost_flat(k)] = fmincon(cost_func, p0_R, A_constraint, b_constraint, [], [], lb_R, ub_R, [], options_surface);
    end
    cost_surface = reshape(cost_flat, size(T1_grid));
    
    surf(T1_grid, T2_grid, cost_surface, 'EdgeColor','none','FaceColor','interp');
    view(3);
    set(gca, 'XScale', 'linear', 'YScale', 'linear');
    xlabel('\tau_1 [s]');
    ylabel('\tau_2 [s]');
    zlabel('RMSE [V]');
    colorbar;
    shading interp;
    hold on;
    
    [min_cost, linIdx] = min(cost_surface(:));
    [r, c] = ind2sub(size(cost_surface), linIdx);
    best_tau1 = tau1_vec(c);
    best_tau2 = tau2_vec(r);
    hStar = plot3(best_tau1, best_tau2, min_cost, 'r*', 'MarkerSize',12,'LineWidth',2);
    legend_str = sprintf('\\tau_1^* = %.3f s,  \\tau_2^* = %.3f s', best_tau1, best_tau2);
    legend(hStar, legend_str, 'Location', 'best');
    title('ECM 2RC Cost Surface: \\tau_1 vs \\tau_2', 'FontSize', 14, 'FontWeight', 'bold');
    hold off;
    fprintf('Cost surface visualization completed.\n');
    fprintf('Optimal tau1: %.3f s, tau2: %.3f s, RMSE: %.6f V\n', best_tau1, best_tau2, min_cost);
end 