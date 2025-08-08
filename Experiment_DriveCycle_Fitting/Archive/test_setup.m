%% 파일 경로 및 데이터 로드 테스트 스크립트
clear; clc; close all;

fprintf('=== 2RC ECM Fitting Setup Test ===\n\n');

%% 1. 파일 경로 확인
fprintf('1. 파일 경로 확인...\n');

% 실부하 데이터 경로
drive_cycle_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\parsedDriveCycle\parsedDriveCycle_0cyc.mat';

% OCV 데이터 경로
ocv_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary\OCV_interpolation_functions.mat';

% 파일 존재 확인
if exist(drive_cycle_path, 'file')
    fprintf('✓ Drive cycle data file found\n');
else
    fprintf('✗ Drive cycle data file NOT found: %s\n', drive_cycle_path);
    return;
end

if exist(ocv_path, 'file')
    fprintf('✓ OCV data file found\n');
else
    fprintf('✗ OCV data file NOT found: %s\n', ocv_path);
    return;
end

%% 2. 데이터 로드 테스트
fprintf('\n2. 데이터 로드 테스트...\n');

try
    % 실부하 데이터 로드
    load(drive_cycle_path);
    fprintf('✓ Drive cycle data loaded successfully\n');
    
    % OCV 데이터 로드
    load(ocv_path);
    fprintf('✓ OCV data loaded successfully\n');
    
catch ME
    fprintf('✗ Error loading data: %s\n', ME.message);
    return;
end

%% 3. 데이터 구조 확인
fprintf('\n3. 데이터 구조 확인...\n');

% Drive cycle 데이터 구조 확인
if exist('parsedDriveCycle_0cyc', 'var')
    fprintf('✓ parsedDriveCycle_0cyc variable found\n');
    
    % ch9 데이터 확인
    if isfield(parsedDriveCycle_0cyc, 'ch9_Drive_0cyc')
        fprintf('✓ ch9_Drive_0cyc field found\n');
        
        ch9_data = parsedDriveCycle_0cyc.ch9_Drive_0cyc;
        
        % SOC90 구조체 확인
        if isfield(ch9_data, 'SOC90')
            fprintf('✓ SOC90 field found\n');
            
            % DC1~DC8 확인
            dc_count = 0;
            for i = 1:8
                dc_name = sprintf('DC%d', i);
                if isfield(ch9_data.SOC90, dc_name)
                    dc_count = dc_count + 1;
                    fprintf('✓ %s field found\n', dc_name);
                    
                    % 필요한 필드 확인
                    dc_data = ch9_data.SOC90.(dc_name);
                    required_fields = {'V', 'I', 't'};
                    for j = 1:length(required_fields)
                        if isfield(dc_data, required_fields{j})
                            if strcmp(required_fields{j}, 't')
                                % duration 타입인 경우 크기만 확인
                                fprintf('  - %s field found (size: %s, type: %s)\n', required_fields{j}, mat2str(size(dc_data.(required_fields{j}))), class(dc_data.(required_fields{j})));
                            else
                                fprintf('  - %s field found (size: %s)\n', required_fields{j}, mat2str(size(dc_data.(required_fields{j}))));
                            end
                        else
                            fprintf('  - %s field NOT found\n', required_fields{j});
                            return;
                        end
                    end
                else
                    fprintf('✗ %s field NOT found\n', dc_name);
                end
            end
            
            fprintf('Total DC datasets found: %d\n', dc_count);
            if dc_count == 0
                fprintf('✗ No DC data found in SOC90 structure\n');
                return;
            end
            
        else
            fprintf('✗ SOC90 field NOT found\n');
            return;
        end
        
    else
        fprintf('✗ ch9_Drive_0cyc field NOT found\n');
        return;
    end
else
    fprintf('✗ parsedDriveCycle_0cyc variable NOT found\n');
    return;
end

% OCV 데이터 구조 확인
if exist('OCV_integrated_RPT0cyc', 'var')
    fprintf('✓ OCV_integrated_RPT0cyc variable found\n');
    
    if isfield(OCV_integrated_RPT0cyc, 'SOC_OCV_func')
        fprintf('✓ SOC_OCV_func field found\n');
        
        % OCV 함수 테스트
        try
            test_soc = 50;  % 50% SOC (0-100 범위)
            test_ocv = OCV_integrated_RPT0cyc.SOC_OCV_func(test_soc);
            fprintf('✓ SOC_OCV function test: SOC=%.0f%% -> OCV=%.3f V\n', test_soc, test_ocv);
        catch ME
            fprintf('✗ SOC_OCV function test failed: %s\n', ME.message);
            return;
        end
        
    else
        fprintf('✗ SOC_OCV_func field NOT found\n');
        return;
    end
else
    fprintf('✗ OCV_integrated_RPT0cyc variable NOT found\n');
    return;
end

%% 4. 데이터 기본 정보 출력 (DC1 예시)
fprintf('\n4. 데이터 기본 정보 (DC1 예시)...\n');

if isfield(ch9_data.SOC90, 'DC1')
    dc1_data = ch9_data.SOC90.DC1;
    V_measured = dc1_data.V;
    I_measured = dc1_data.I;
    t_measured = seconds(dc1_data.t);  % duration을 seconds로 변환
    
    fprintf('  데이터 포인트 수: %d\n', length(V_measured));
    fprintf('  시간 범위: %.2f - %.2f s (%.2f hr)\n', min(t_measured), max(t_measured), (max(t_measured) - min(t_measured))/3600);
    fprintf('  전압 범위: %.3f - %.3f V\n', min(V_measured), max(V_measured));
    fprintf('  전류 범위: %.3f - %.3f A\n', min(I_measured), max(I_measured));
    fprintf('  평균 시간 간격: %.3f s\n', mean(diff(t_measured)));
else
    fprintf('  DC1 데이터 없음 - 기본 정보 확인 불가\n');
end

%% 5. 필요한 함수 파일 확인
fprintf('\n5. 필요한 함수 파일 확인...\n');

required_functions = {'ECM_2RC_model.m', 'cost_function.m', 'estimate_and_validate_soc.m', 'visualize_cost_surface.m'};

for i = 1:length(required_functions)
    if exist(required_functions{i}, 'file')
        fprintf('✓ %s found\n', required_functions{i});
    else
        fprintf('✗ %s NOT found\n', required_functions{i});
    end
end

%% 6. 간단한 시각화 (DC1 예시)
fprintf('\n6. 데이터 시각화 (DC1 예시)...\n');

if isfield(ch9_data.SOC90, 'DC1')
    figure('Position', [100, 100, 1200, 400]);
    
    subplot(1, 3, 1);
    plot(t_measured/3600, V_measured, 'b-', 'LineWidth', 1);
    xlabel('Time (hr)');
    ylabel('Voltage (V)');
    title('DC1 Measured Voltage');
    grid on;
    
    subplot(1, 3, 2);
    plot(t_measured/3600, I_measured, 'r-', 'LineWidth', 1);
    xlabel('Time (hr)');
    ylabel('Current (A)');
    title('DC1 Current Profile');
    grid on;
    
    subplot(1, 3, 3);
    % OCV 곡선 그리기 (원본 데이터 사용)
    plot(OCV_integrated_RPT0cyc.SOC_grid, OCV_integrated_RPT0cyc.V_avg_Q(1,:), 'g-', 'LineWidth', 2);
    xlabel('SOC (%)');
    ylabel('OCV (V)');
    title('OCV Curve');
    grid on;
else
    fprintf('  DC1 데이터 없음 - 시각화 불가\n');
end

fprintf('\n=== Setup test completed successfully! ===\n');
fprintf('모든 데이터와 함수가 준비되었습니다. ECM_2RC_Fitting.m을 실행할 수 있습니다.\n'); 