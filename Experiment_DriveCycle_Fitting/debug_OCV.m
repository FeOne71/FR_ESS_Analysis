%% OCV 함수 테스트 스크립트
clear; clc;

% OCV 데이터 로드
ocv_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary\OCV_interpolation_functions.mat';

load(ocv_path);

fprintf('=== OCV 함수 테스트 ===\n');

% OCV 구조체 필드 확인
fprintf('OCV 구조체 필드들:\n');
fields = fieldnames(OCV_integrated_RPT0cyc);
for i = 1:length(fields)
    if strcmp(fields{i}, 'SOC_OCV_func')
        fprintf('  - %s: function handle\n', fields{i});
    elseif strcmp(fields{i}, 'Q_OCV_func')
        fprintf('  - %s: function handle\n', fields{i});
    else
        data = OCV_integrated_RPT0cyc.(fields{i});
        try
            if isnumeric(data)
                fprintf('  - %s: %s, range: %.4f ~ %.4f\n', fields{i}, mat2str(size(data)), min(data(:)), max(data(:)));
            else
                fprintf('  - %s: %s, type: %s\n', fields{i}, mat2str(size(data)), class(data));
            end
        catch
            fprintf('  - %s: %s, type: %s\n', fields{i}, mat2str(size(data)), class(data));
        end
    end
end

fprintf('\n=== 실제 데이터 확인 ===\n');
fprintf('SOC_grid 범위: %.2f ~ %.2f%%\n', min(OCV_integrated_RPT0cyc.SOC_grid), max(OCV_integrated_RPT0cyc.SOC_grid));
if isfield(OCV_integrated_RPT0cyc, 'V_avg_Q')
    fprintf('V_avg_Q 크기: %s\n', mat2str(size(OCV_integrated_RPT0cyc.V_avg_Q)));
    fprintf('V_avg_Q 범위: %.4f ~ %.4f V\n', min(OCV_integrated_RPT0cyc.V_avg_Q(:)), max(OCV_integrated_RPT0cyc.V_avg_Q(:)));
end

fprintf('\n=== 특정 SOC 테스트 (8.3%%) ===\n');
try
    test_soc = 8.3;  % 실제 initial SOC
    if isfield(OCV_integrated_RPT0cyc, 'SOC_OCV_func')
        ocv_val = OCV_integrated_RPT0cyc.SOC_OCV_func(test_soc);
        fprintf('SOC_OCV_func(%.1f%%) = %.4f V\n', test_soc, ocv_val);
    end
    
    if isfield(OCV_integrated_RPT0cyc, 'Q_OCV_func')
        capacity_val = OCV_integrated_RPT0cyc.Q_OCV_func(test_soc);
        fprintf('Q_OCV_func(%.1f%%) = %.2f Ah\n', test_soc, capacity_val);
    end
catch ME
    fprintf('에러: %s\n', ME.message);
end

fprintf('\n=== OCV 함수 호출 테스트 ===\n');

% 1. 0~100 범위로 테스트
fprintf('1. SOC 0~100 범위 테스트:\n');
try
    test_socs_01 = [10, 50, 90];  % 0-100 범위로 수정
    for soc = test_socs_01
        ocv_val = OCV_integrated_RPT0cyc.SOC_OCV_func(soc);
        fprintf('  SOC=%.1f%% -> OCV=%.3f V\n', soc, ocv_val);
    end
catch ME
    fprintf('  에러: %s\n', ME.message);
end

% 2. 1~100 범위로 테스트  
fprintf('\n2. SOC 1~100 범위 테스트:\n');
try
    test_socs_100 = [10, 50, 90];
    for soc = test_socs_100
        ocv_val = OCV_integrated_RPT0cyc.SOC_OCV_func(soc);
        fprintf('  SOC=%d%% -> OCV=%.3f V\n', soc, ocv_val);
    end
catch ME
    fprintf('  에러: %s\n', ME.message);
end

% 3. 실제 soc_grid 값들로 테스트
fprintf('\n3. 실제 SOC_grid 값들로 테스트:\n');
try
    test_indices = [1, 50, 100];
    for idx = test_indices
        if idx <= length(OCV_integrated_RPT0cyc.SOC_grid)
            soc_val = OCV_integrated_RPT0cyc.SOC_grid(idx);
            ocv_val = OCV_integrated_RPT0cyc.SOC_OCV_func(soc_val);
            if isfield(OCV_integrated_RPT0cyc, 'V_avg_Q')
                if size(OCV_integrated_RPT0cyc.V_avg_Q, 1) > 1
                    expected_ocv = OCV_integrated_RPT0cyc.V_avg_Q(1, idx);
                else
                    expected_ocv = OCV_integrated_RPT0cyc.V_avg_Q(idx);
                end
                fprintf('  SOC_grid(%d)=%.1f%% -> SOC_OCV_func=%.3f, V_avg_Q=%.3f\n', idx, soc_val, ocv_val, expected_ocv);
            end
        end
    end
catch ME
    fprintf('  에러: %s\n', ME.message);
end

fprintf('\n=== 배터리 용량 함수 테스트 ===\n');
try
    test_socs = [5, 10, 50, 90];
    for soc = test_socs
        if isfield(OCV_integrated_RPT0cyc, 'Q_OCV_func')
            capacity = OCV_integrated_RPT0cyc.Q_OCV_func(soc);
            fprintf('  SOC=%d%% -> Capacity=%.2f Ah\n', soc, capacity);
        end
    end
catch ME
    fprintf('  에러: %s\n', ME.message);
end

fprintf('\n=== 결론 ===\n');
fprintf('함수들이 올바른 값을 반환하는지 확인하세요!\n');

fprintf('\n=== SOC_grid 실제 값 확인 ===\n');
try
    % Find nearest SOC values to common test points
    test_socs = [0, 10, 50, 90, 100];
    soc_grid = OCV_integrated_RPT0cyc.SOC_grid;
    
    for target_soc = test_socs
        [~, nearest_idx] = min(abs(soc_grid - target_soc));
        nearest_soc = soc_grid(nearest_idx);
        fprintf('  Target %.0f%% -> Nearest: %.3f%% (index %d)\n', target_soc, nearest_soc, nearest_idx);
        
        % Test the function with the actual grid value
        if isfield(OCV_integrated_RPT0cyc, 'SOC_OCV_func')
            ocv_val = OCV_integrated_RPT0cyc.SOC_OCV_func(nearest_soc);
            fprintf('    SOC_OCV_func(%.3f%%) = %.4f V\n', nearest_soc, ocv_val);
        end
    end
    
    % Show some sample SOC_grid values around 90%
    fprintf('\n  SOC_grid values around 90%%:\n');
    idx_90 = find(soc_grid >= 89.5 & soc_grid <= 90.5);
    if ~isempty(idx_90)
        for i = 1:min(10, length(idx_90))
            fprintf('    [%d]: %.6f%%\n', idx_90(i), soc_grid(idx_90(i)));
        end
    end
    
catch ME
    fprintf('  에러: %s\n', ME.message);
end 