function cost = cost_function_overpotential(params, I, V_overpotential_measured, t, weight)
    % Cost function for overpotential-based ECM parameter optimization with weight
    % Input:
    %   params: [R0, R1, R2, tau1, tau2] parameter vector
    %   I: Current vector [A] (positive: charge, negative: discharge)
    %   V_overpotential_measured: Measured overpotential [V] = V_terminal - V_OCV
    %   t: Time vector [s]
    %   weight: Weight vector for each data point (same size as V_overpotential_measured)
    % Output:
    %   cost: Weighted Root Mean Square Error (RMSE) [V]
    
        % Predict overpotential using 2RC ECM model
        V_overpotential_model = ECM_2RC_model_overpotential(params, I, t);
        
        % Calculate Weighted RMSE (same as func_cost in HNE_ECM_2RC_fitting.m)
        cost = sqrt(sum(((V_overpotential_measured - V_overpotential_model).^2) .* weight) / sum(weight));
end 