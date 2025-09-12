function cost = cost_function_overpotential_improved(params, I_measured, V_overpotential_measured, t_measured, weight)
    % Improved cost function for overpotential with better initial condition matching
    % Input:
    %   params: [R0, R1, R2, tau1, tau2] parameter vector
    %   I_measured: Measured current [A]
    %   V_overpotential_measured: Measured overpotential [V]
    %   t_measured: Time vector [s]
    %   weight: Weight vector for data points
    % Output:
    %   cost: RMSE cost value
    
    try
        % Generate model overpotential
        V_overpotential_model = ECM_2RC_model_overpotential(params, I_measured, t_measured);
        
        % Calculate initial offset and apply correction
        initial_offset = V_overpotential_measured(1) - V_overpotential_model(1);
        V_overpotential_model_corrected = V_overpotential_model + initial_offset;
        
        % Calculate weighted RMSE
        residual = V_overpotential_measured - V_overpotential_model_corrected;
        cost = sqrt(sum((residual.^2) .* weight) / sum(weight));
        
        % Add penalty for unrealistic parameters
        R0 = params(1); R1 = params(2); R2 = params(3);
        tau1 = params(4); tau2 = params(5);
        
        % Penalty for parameter ordering violations
        penalty = 0;
        if tau1 >= tau2
            penalty = penalty + 0.1;  % tau1 should be < tau2
        end
        if R0 > max(R1, R2)
            penalty = penalty + 0.05;  % R0 should be small
        end
        
        % Penalty for extreme capacitance values
        C1 = tau1 / R1;
        C2 = tau2 / R2;
        if C1 > 100000 || C2 > 500000  % Very large capacitance
            penalty = penalty + 0.1;
        end
        if C1 < 1 || C2 < 1  % Very small capacitance
            penalty = penalty + 0.1;
        end
        
        cost = cost + penalty;
        
    catch ME
        % Return large penalty if calculation fails
        fprintf('Error in cost function: %s\n', ME.message);
        cost = 1e6;
    end
end 