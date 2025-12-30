function daily_por = calculate_por(t_vec, P, V, I, fs, Cnom_cell)
% [개선] 충전 펄스와 방전 펄스를 명확하게 분리하여 POR을 계산합니다.

    % --- 파라미터 ---
    IDLE_THR_A = Cnom_cell * 0.05;
    MIN_DUR_SEC = 30;
    STABILITY_THR = 0.1;

    if isempty(I) || length(I) < floor(MIN_DUR_SEC * fs)
        daily_por = NaN;
        return;
    end
    
    % --- 충전 펄스와 방전 펄스를 별도로 처리 ---
    por_values_ch = find_por_in_pulses(P, V, I, fs, IDLE_THR_A, MIN_DUR_SEC, STABILITY_THR, 'charge');
    por_values_dch = find_por_in_pulses(P, V, I, fs, IDLE_THR_A, MIN_DUR_SEC, STABILITY_THR, 'discharge');
    
    % 계산된 모든 POR 값들을 합침
    all_por_values = [por_values_ch; por_values_dch];
    
    if isempty(all_por_values)
        daily_por = NaN;
    else
        daily_por = median(all_por_values, 'omitnan');
    end
end

% --- 헬퍼 함수: 특정 방향의 펄스를 찾아 POR 계산 ---
function por_values = find_por_in_pulses(P, V, I, fs, idle_thr, min_dur, stability_thr, direction)
    
    por_values = [];
    
    % 방향에 따라 펄스 후보군 정의
    if strcmp(direction, 'charge')
        is_event = I > idle_thr; % 충전 이벤트
    elseif strcmp(direction, 'discharge')
        is_event = I < -idle_thr; % 방전 이벤트
    else
        return;
    end
    
    diff_event = diff([0; is_event; 0]);
    pulse_starts = find(diff_event == 1);
    pulse_ends = find(diff_event == -1) - 1;
    
    for i = 1:length(pulse_starts)
        start_idx = pulse_starts(i);
        end_idx = pulse_ends(i);
        
        duration = (end_idx - start_idx) / fs;
        if duration < min_dur
            continue;
        end
        
        P_pulse = P(start_idx:end_idx);
        
        idx_3s = floor(3 * fs);
        idx_10s = floor(10 * fs);
        idx_11s = floor(11 * fs);
        idx_30s = floor(30 * fs);
        
        if length(P_pulse) < idx_30s || idx_3s < 1
            continue;
        end
        
        p_avg_3_10 = mean(P_pulse(idx_3s:idx_10s));
        p_avg_11_30 = mean(P_pulse(idx_11s:idx_30s));
        
        if abs(p_avg_3_10) < 1e-3
            continue;
        end
        
        if abs(p_avg_11_30 - p_avg_3_10) / abs(p_avg_3_10) < stability_thr
            V_pulse = V(start_idx:end_idx);
            
            v0_idx = 1;
            v1_idx = floor(1 * fs);
            v30_idx = floor(30 * fs);
            
            if v1_idx <= v0_idx, v1_idx = v0_idx + 1; end
            if length(V_pulse) < v30_idx || v30_idx <= v1_idx
                continue;
            end

            delta_V_ohmic = V_pulse(v1_idx) - V_pulse(v0_idx);
            delta_V_polarization = V_pulse(v30_idx) - V_pulse(v1_idx);
            
            if abs(delta_V_ohmic) < 1e-4
                continue;
            end
            
            current_por = abs(delta_V_polarization / delta_V_ohmic);
            
            if current_por > 0 && current_por < 100
                por_values = [por_values; current_por];
            end
        end
    end
end