function [SOC_initial, battery_capacity, validation_result, initial_rest_end, SOC_full, final_rest_end] = estimate_and_validate_soc(V_measured, I_measured, t_measured, OCV_data)
    % Estimate initial SOC and validate final SOC consistency
    % NEW: Calculate full SOC profile using linear interpolation + current integration ratio
    
    % Extract OCV data and create inverse function (V→SOC)
    soc_grid = OCV_data.SOC_grid;
    ocv_values = OCV_data.V_avg_SOC;
    
    % Create inverse OCV function (Voltage → SOC)
    % Ensure the OCV array is strictly monotonic for reliable inversion
    [ocv_values_sorted, uniqueIdx] = unique(ocv_values, 'stable');
    soc_grid_sorted = soc_grid(uniqueIdx);
    inverse_OCV_func = @(voltage) interp1(ocv_values_sorted, soc_grid_sorted, voltage, 'linear');
    
    % Helper function to find nearest SOC in grid
    find_nearest_soc = @(target_soc) soc_grid(min(abs(soc_grid - target_soc)) == abs(soc_grid - target_soc));
    
    % Find rest periods (|I| <= 2A) - 초기 휴지기 5991 포인트 포함
    rest_mask = abs(I_measured) <= 2.0;
    min_rest_time_points = 5990;  % 5991 포인트 포함하도록 기준 설정
    
    % Find continuous rest periods
    rest_periods = [];
    in_rest = false;
    rest_start = 0;
    
    for i = 1:length(rest_mask)
        if rest_mask(i) && ~in_rest
            rest_start = i;
            in_rest = true;
        elseif ~rest_mask(i) && in_rest
            rest_duration_points = (i-1) - rest_start + 1;
            if rest_duration_points >= min_rest_time_points
                rest_periods = [rest_periods; rest_start, i-1];
            end
            in_rest = false;
        end
    end
    
    % Check if last period is still ongoing
    if in_rest
        rest_duration_points = length(rest_mask) - rest_start + 1;
        if rest_duration_points >= min_rest_time_points
            rest_periods = [rest_periods; rest_start, length(rest_mask)];
        end
    end
    
    % Debug: Check rest periods found
    fprintf('Debug: Found %d rest periods\n', size(rest_periods, 1));
    for i = 1:size(rest_periods, 1)
        fprintf('Debug: Rest period %d: index %d to %d (duration: %d points)\n', ...
            i, rest_periods(i, 1), rest_periods(i, 2), rest_periods(i, 2) - rest_periods(i, 1) + 1);
    end
    
    % 초기 SOC (SOC1): 첫 10분 휴지기 마지막 시점
    initial_rest_end = rest_periods(1, 2);
    V_ocv_initial = V_measured(initial_rest_end);
    
    % Debug: Check voltage and OCV range
    fprintf('Debug: Initial voltage: %.3f V\n', V_ocv_initial);
    fprintf('Debug: OCV range: %.3f to %.3f V\n', min(ocv_values), max(ocv_values));
    
    % Check if voltage is within OCV range
    if V_ocv_initial < min(ocv_values) || V_ocv_initial > max(ocv_values)
        fprintf('WARNING: Initial voltage %.3f V is outside OCV range [%.3f, %.3f] V\n', ...
            V_ocv_initial, min(ocv_values), max(ocv_values));
    end
    
    SOC_initial_raw = inverse_OCV_func(V_ocv_initial);
    
    % Debug: Check SOC calculation
    fprintf('Debug: Initial SOC calculated: %.2f%%\n', SOC_initial_raw);
    
    % Use raw SOC value directly (no rounding to grid)
    SOC_initial = SOC_initial_raw;  % SOC1
    
    % Battery capacity: 구조체에 있는 mean_capacity 직접 사용
    battery_capacity = OCV_data.mean_capacity;  % 64.9767 Ah
    
    % 후기 SOC (SOC2): 뒤 10분 휴지기 마지막 시점 (데이터 끝)
    final_rest_end = rest_periods(end, 2);
    V_ocv_final = V_measured(final_rest_end);
    
    % Debug: Check final voltage
    fprintf('Debug: Final voltage: %.3f V\n', V_ocv_final);
    
    % Check if voltage is within OCV range
    if V_ocv_final < min(ocv_values) || V_ocv_final > max(ocv_values)
        fprintf('WARNING: Final voltage %.3f V is outside OCV range [%.3f, %.3f] V\n', ...
            V_ocv_final, min(ocv_values), max(ocv_values));
    end
    
    SOC_final_raw = inverse_OCV_func(V_ocv_final);
    
    % Debug: Check SOC calculation
    fprintf('Debug: Final SOC calculated: %.2f%%\n', SOC_final_raw);
    
    % Use raw SOC value directly (no rounding to grid)
    SOC_final = SOC_final_raw;  % SOC2
    
    % NEW: Calculate full SOC profile using new method
    % SOC(t) = SOC1 + (SOC2-SOC1) * (∫₀ᵗ I*dt) / (∫₀ᵉⁿᵈ I*dt)
    % where 0 = initial_rest_end, end = final_rest_end
    
    N = length(I_measured);
    dt = 0.1;  % Time step (seconds)
    SOC_full = zeros(N, 1);
    
    % Calculate cumulative current integration from initial_rest_end (SOC1 시점)
    cumulative_current = zeros(N, 1);
    for i = (initial_rest_end + 1):final_rest_end
        cumulative_current(i) = cumulative_current(i-1) + I_measured(i) * dt;
    end
    
    % Total current integration from initial_rest_end to final_rest_end
    total_current_integration = cumulative_current(final_rest_end);
    
    % Apply new SOC calculation method
    for i = 1:N
        if i <= initial_rest_end
            % Phase 1: Initial rest period - constant SOC1
            SOC_full(i) = SOC_initial;
        elseif i >= final_rest_end
            % Phase 3: Final rest period - constant SOC2
            SOC_full(i) = SOC_final;
        else
            % Phase 2: Active period - always use current integration ratio
            SOC_full(i) = SOC_initial + (SOC_final - SOC_initial) * (cumulative_current(i) / total_current_integration);
        end
    end
    
    % Validation result
    validation_result.has_validation = true;
    validation_result.initial_ocv = V_ocv_initial;
    validation_result.final_ocv = V_ocv_final;
    validation_result.initial_soc = SOC_initial;
    validation_result.final_soc = SOC_final;
    validation_result.total_current_integration = total_current_integration;
    validation_result.soc_method = 'linear_interpolation_current_ratio';
    
    fprintf('=== NEW SOC Estimation Method ===\n');
    fprintf('SOC1 (initial): %.2f%% at index %d (%.4f V)\n', SOC_initial, initial_rest_end, V_ocv_initial);
    fprintf('SOC2 (final): %.2f%% at index %d (%.4f V)\n', SOC_final, final_rest_end, V_ocv_final);
    fprintf('Total current integration: %.2f A·s\n', total_current_integration);
    fprintf('SOC calculation: SOC(t) = SOC1 + (SOC2-SOC1) * (∫₀ᵗ I*dt) / (∫₀ᵉⁿᵈ I*dt)\n');
    fprintf('SOC range: %.2f%% to %.2f%%\n', min(SOC_full), max(SOC_full));
end 