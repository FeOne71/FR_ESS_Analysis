function [SOH_p, Delta_Q, Delta_SOC, SOC_est, V_est] = Physics_Prior_Calculator(V_meas, I_meas, t_sec, SOC0, Q_cap_Ah, mr_path)
% Physics_Prior_Calculator
% Estimates the Prior SOH (SOH_p) based on AEKF Delta SOC and Coulomb Counting Delta Q.
%
% Inputs:
%   V_meas: Measured voltage array (V)
%   I_meas: Measured current array (A)
%   t_sec: Time array (seconds)
%   SOC0: Initial SOC guess (0~1)
%   Q_cap_Ah: Nominal cell capacity (Ah)
%   mr_path: Path to MasterRulers_v3.mat for Fresh OCV extraction
% Outputs:
%   SOH_p: Prior SOH Estimate (0~1)
%   Delta_Q: Coulomb counted capacity change (Ah)
%   Delta_SOC: AEKF estimated SOC change (0~1)
%   SOC_est: Full AEKF SOC profile (array)
%   V_est: Full AEKF Voltage profile (array)

    % 1. Load Fresh OCV Reference
    try
        S = load(mr_path, 'MasterRulers');
        fns = fieldnames(S.MasterRulers);
        sample_ch = fns{1};
        Fresh_OCV = S.MasterRulers.(sample_ch).Fresh_OCV_Charge;
        
        % Normalize Q to get SOC (0 to 1)
        Q_max = max(Fresh_OCV.Q);
        OCV_soc = Fresh_OCV.Q / Q_max;
        OCV_v = Fresh_OCV.V_grid;
    catch
        error('Could not load Fresh OCV from %s', mr_path);
    end
    
    % 2. Run AEKF to get precise Delta SOC
    % Note: I_meas should be > 0 for charge (SOC increase).
    % If field data uses I_meas > 0 for charge, we are fine.
    % To be safe with the integration later, we'll use abs() for Delta_Q,
    % but AEKF needs signed I_meas matching convention.
    [SOC_est, V_est, ~] = AEKF_SOC_Estimator(V_meas, I_meas, t_sec, SOC0, Q_cap_Ah, OCV_soc, OCV_v);
    
    % Get the absolute SOC difference across the segment
    Delta_SOC = abs(SOC_est(end) - SOC_est(1));
    
    % 3. Calculate actual capacity change via Coulomb Counting
    % dt = diff(t_sec), trapz or cumsum
    % We handle variable time steps
    if length(t_sec) > 1
        dt_array = diff(t_sec);
        I_int = I_meas(1:end-1);
        % Integrated charge passed
        Delta_Q = sum(abs(I_int) .* dt_array) / 3600;
    else
        Delta_Q = 0;
    end
    
    % 4. Compute Prior SOH
    if Delta_SOC > 1e-4 % Avoid divide-by-zero
        SOH_p = (Delta_Q / Delta_SOC) / Q_cap_Ah;
    else
        SOH_p = NaN;
    end
    
    % Bound check purely for robustness
    if ~isnan(SOH_p)
        SOH_p = max(0.2, min(1.5, SOH_p));
    end
end
