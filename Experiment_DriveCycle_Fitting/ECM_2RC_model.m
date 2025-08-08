function V_terminal = ECM_2RC_model(params, I, SOC, t, ocv_func)
    % 2RC model function (Reference 방식 + OCV 포함)
    
    % params: [R0, R1, R2, tau1, tau2]
    R0   = params(1);
    R1   = params(2);
    R2   = params(3);
    tau1 = params(4);
    tau2 = params(5);

    % dt 계산 (0.1초 간격 유지)
    dt = [0.1; diff(t)];

    N = length(t);
    V_terminal = zeros(N, 1);

    % OCV vector (07/14 17:27 Added)
    soc_rounded = round(SOC, 2);
    V_oc_vec = ocv_func(soc_rounded);

    for k = 1:N
        % (1) OCV 계산 (Reference와 다른 점: OCV 포함)
        % Handle SOC values that may not exist exactly in the grid
        % try
        %    V_oc = ocv_func(SOC(k));
        % catch
        %    % If function fails, try with nearest grid value
        %    % This is a fallback - ideally SOC should already be grid-aligned
        %    fprintf('Warning: SOC_OCV_func failed for SOC=%.3f%%, using nearest value\n', SOC(k));
        %    V_oc = ocv_func(round(SOC(k), 2));  % Round to 2 decimal places
        % end
        
        % V_oc = ocv_func(round(SOC(k), 2));  % Round to 2 decimal places
        V_oc = V_oc_vec(k);

        % (2) R0 전압강하
        IR0 = R0 * I(k);

        % (3) RC1, RC2 업데이트
        alpha1 = exp(-dt(k)/tau1);
        alpha2 = exp(-dt(k)/tau2);

        if k == 1
            % k-1번째 = 0번째 (RC에 대한 전압강하 없음)
            Vrc1 = 0;
            Vrc2 = 0;
        else
            % 이후는 기존 공식
            Vrc1 = Vrc1*alpha1 + R1*(1 - alpha1)*I(k-1);
            Vrc2 = Vrc2*alpha2 + R2*(1 - alpha2)*I(k-1);
        end

        % (4) 최종 전압
        V_terminal(k) = V_oc + IR0 + Vrc1 + Vrc2;
    end
end 