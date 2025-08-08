function V_overpotential = ECM_2RC_model_overpotential(params, I, t)
    % Equivalent Circuit Model (ECM) with 2RC network for overpotential only
    % Input:
    %   params: [R0, R1, R2, tau1, tau2] parameter vector
    %   I: Current vector [A] (positive: charge, negative: discharge)
    %   t: Time vector [s]
    % Output:
    %   V_overpotential: Overpotential vector [V] = R0*I + V_C1 + V_C2
    
    % Extract parameters
    R0 = params(1);    % Ohmic resistance [Ω]
    R1 = params(2);    % Fast RC resistance [Ω]
    R2 = params(3);    % Slow RC resistance [Ω]
    tau1 = params(4);  % Fast time constant [s]
    tau2 = params(5);  % Slow time constant [s]
    
    % Initialize arrays
    N = length(t);
    V_overpotential = zeros(N, 1);
    
    % Calculate time step vector 
    dt_vec = diff(t);
    dt = [dt_vec(1); dt_vec];  % 첫 번째 dt = 실제 시간 간격 (0.1초)
    
    % Initialize RC voltages (scalar variables like BSG)
    Vrc1 = 0;
    Vrc2 = 0;
    
    % Time-domain simulation (BSG method)
    for k = 1:N
        % R0 voltage drop
        IR0 = R0 * I(k);
        
        % RC voltage update with time-varying dt
        alpha1 = exp(-dt(k)/tau1);
        alpha2 = exp(-dt(k)/tau2);
        
        if k == 1
            % Initial condition: RC voltages start from zero
            Vrc1 = 0;
            Vrc2 = 0;
        else
            % RC update using previous current
            Vrc1 = Vrc1 * alpha1 + R1 * (1 - alpha1) * I(k-1);
            Vrc2 = Vrc2 * alpha2 + R2 * (1 - alpha2) * I(k-1);
        end
        
        % Overpotential calculation (NO OCV component)
        V_overpotential(k) = IR0 + Vrc1 + Vrc2;
    end
end 