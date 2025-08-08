function cost = cost_function(params, I, SOC, V_measured, t, ocv_func)
    % Cost function for ECM parameter optimization (HM_2RC.m method)
    % Input:
    %   params: [R0, R1, R2, tau1, tau2] parameter vector
    %   I: Current vector [A] (positive: charge, negative: discharge)
    %   SOC: State of Charge vector [1-100%]
    %   V_measured: Measured terminal voltage [V]
    %   t: Time vector [s]
    %   ocv_func: OCV function handle
    % Output:
    %   cost: Root Mean Square Error (RMSE) [V]
    
    % Predict voltage using 2RC ECM model
    V_model = ECM_2RC_model(params, I, SOC, t, ocv_func);
    
    % Calculate RMSE (same as RMSE_2RC in HM_2RC.m)
    cost = sqrt(mean((V_measured - V_model).^2));
end 