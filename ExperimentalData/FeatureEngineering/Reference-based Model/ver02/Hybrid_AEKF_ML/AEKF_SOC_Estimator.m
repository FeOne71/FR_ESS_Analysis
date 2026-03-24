function [SOC_est, V_est, ECM_Params] = AEKF_SOC_Estimator(V_meas, I_meas, t_sec, SOC0, Q_cap_Ah, OCV_soc, OCV_v)
% AEKF_SOC_Estimator
% 1st-order Thevenin Equivalent Circuit Model (ECM) based Adaptive Extended Kalman Filter (AEKF)
% Inputs:
%   V_meas: Measured voltage array (V)
%   I_meas: Measured current array (A) (I > 0 for charge)
%   t_sec: Time array (seconds)
%   SOC0: Initial SOC (0~1)
%   Q_cap_Ah: Nominal Capacity (Ah)
%   OCV_soc: Reference SOC array for OCV lookup (0~1)
%   OCV_v: Reference OCV array (V)
% Outputs:
%   SOC_est: Estimated SOC array over time (0~1)
%   V_est: Estimated Voltage array over time (V)
%   ECM_Params: Estimated parameters (struct) if adapted, or fixed assumed.

    % Validate lengths
    N = length(V_meas);
    SOC_est = zeros(N, 1);
    V_est   = zeros(N, 1);
    
    % Ensure OCV lookup is strictly monotonic ascending for interpolation
    [OCV_soc_uq, uniq_idx] = unique(OCV_soc);
    OCV_v_uq = OCV_v(uniq_idx);
    
    % Initialize State: x = [SOC; V_c]
    x = [SOC0; 0];
    
    % Process & Measurement Noise Covariances
    P = diag([1e-4, 1e-4]);     % Initial error covariance
    Q = diag([1e-6, 1e-5]);     % Process noise covariance
    R = 1e-2;                   % Initial Measurement noise covariance
    alpha_R = 0.95;             % For fading memory adaptation of R
    
    % ECM Parameters (Assumed constants for typical LFP, could be parameterized)
    R0 = 0.015;  % Ohmic resistance (Ohms)
    R1 = 0.010;  % Polarization resistance (Ohms)
    C1 = 3000;   % Polarization capacitance (Farads)
    tau1 = R1 * C1;
    eta = 1.0;   % Coulombic efficiency
    
    % EKF Loop
    for k = 1:N
        % 1. Time Update (Prediction)
        if k == 1
            dt = 1; % Assume 1s if first step
            I_prev = I_meas(1);
        else
            dt = t_sec(k) - t_sec(k-1);
            if dt <= 0, dt = 1; end
            I_prev = I_meas(k-1);
        end
        
        % State Transition Matrix F
        F = [1, 0; 
             0, exp(-dt/tau1)];
        
        % Input Matrix B
        B = [eta * dt / (3600 * Q_cap_Ah); 
             R1 * (1 - exp(-dt/tau1))];
             
        % A priori state estimate
        x = F * x + B * I_prev;
        
        % Constrain SOC bounds
        x(1) = max(0, min(1, x(1)));
        
        % A priori covariance
        P = F * P * F' + Q;
        
        % 2. Measurement Update (Correction)
        % Predict measurement: V_t = OCV(SOC) + V_c + I_k * R0
        % Evaluate OCV and its derivative dOCV/dSOC
        ocv_k = interp1(OCV_soc_uq, OCV_v_uq, x(1), 'linear', 'extrap');
        
        % Numerical derivative for dOCV/dSOC
        delta_soc = 0.001;
        ocv_plus = interp1(OCV_soc_uq, OCV_v_uq, x(1) + delta_soc, 'linear', 'extrap');
        ocv_minus = interp1(OCV_soc_uq, OCV_v_uq, x(1) - delta_soc, 'linear', 'extrap');
        dOCV = (ocv_plus - ocv_minus) / (2 * delta_soc);
        
        V_pred = ocv_k + x(2) + I_meas(k) * R0;
        V_est(k) = V_pred;
        
        % Jacobian of measurement function H
        H = [dOCV, 1];
        
        % Innovation
        e = V_meas(k) - V_pred;
        
        % Adaptive R (AEKF component for robustness to field noise)
        R_hat = e^2 - H * P * H';
        R_hat = max(1e-4, R_hat); % bound lower to avoid singularity
        R = alpha_R * R + (1 - alpha_R) * R_hat;
        
        % Kalman Gain
        K = P * H' / (H * P * H' + R);
        
        % A posteriori state estimate
        x = x + K * e;
        x(1) = max(0, min(1, x(1))); % Constraint
        
        % A posteriori covariance
        P = (eye(2) - K * H) * P;
        
        % Save output
        SOC_est(k) = x(1);
    end
    
    ECM_Params.R0 = R0;
    ECM_Params.R1 = R1;
    ECM_Params.C1 = C1;
    ECM_Params.R_mat = R; % Final R
end
