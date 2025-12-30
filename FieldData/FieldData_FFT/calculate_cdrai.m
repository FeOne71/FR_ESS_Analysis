function daily_cdrai = calculate_cdrai(I, V_cell, SOC, Cnom_cell, min_points)
% [최종 수정] 셀 전압을 사용하고, 이상치 제거를 강화하여 CDRAI 안정성 향상

    idle_thr = Cnom_cell * 0.05;
    if length(I) < 100, daily_cdrai = NaN; return; end

    dV = diff(V_cell);
    dI = diff(I);
    dI(abs(dI) < 0.1) = NaN; 
    R_dc = dV ./ dI;
    
    % 물리적 범위 필터링 (셀 저항은 보통 mOhm 단위)
    R_dc(R_dc <= 0 | R_dc > 0.01) = NaN; % 10mOhm 이상은 노이즈로 간주

    I_mid = (I(1:end-1) + I(2:end))/2;
    SOC_mid = (SOC(1:end-1) + SOC(2:end))/2;
    
    charge_idx = I_mid > idle_thr;
    dch_idx = I_mid < -idle_thr;
    
    cdrai_values = [];
    soc_bins = (40:5:95)';
    
    for i = 1:length(soc_bins)-1
        soc_target_idx = (SOC_mid >= soc_bins(i)) & (SOC_mid < soc_bins(i+1));
        
        R_ch_soc = R_dc(charge_idx & soc_target_idx);
        R_dch_soc = R_dc(dch_idx & soc_target_idx);
        
        if sum(~isnan(R_ch_soc)) < min_points || sum(~isnan(R_dch_soc)) < min_points
            continue;
        end

        % [이상치 제거 강화] IQR(사분위수 범위) 기반의 더 표준적인 방법 사용
        R_ch_soc_filtered = R_ch_soc(~isoutlier(R_ch_soc, 'quartiles'));
        R_dch_soc_filtered = R_dch_soc(~isoutlier(R_dch_soc, 'quartiles'));

        if isempty(R_ch_soc_filtered) || isempty(R_dch_soc_filtered)
            continue;
        end
        
        avg_R_ch = mean(R_ch_soc_filtered, 'omitnan');
        avg_R_dch = mean(R_dch_soc_filtered, 'omitnan');
        
        if avg_R_dch == 0, continue; end
        
        current_cdrai = avg_R_ch / avg_R_dch;

        if current_cdrai > 0 && current_cdrai < 50
            cdrai_values = [cdrai_values; current_cdrai];
        end
    end
    
    if isempty(cdrai_values)
        daily_cdrai = NaN;
    else
        daily_cdrai = median(cdrai_values, 'omitnan');
    end
end