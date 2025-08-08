%% Cost Function 디버그 스크립트
clear; clc;

% 데이터 로드
drive_cycle_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\parsedDriveCycle\parsedDriveCycle_0cyc.mat';
load(drive_cycle_path);

ocv_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary\OCV_interpolation_functions.mat';

load(ocv_path);

ch9_data = parsedDriveCycle_0cyc.ch9_Drive_0cyc;
ocv_func = OCV_integrated_RPT0cyc.SOC_OCV_func;

fprintf('=== Cost Function 디버그 ===\n');

% DC1 데이터 추출
dc_data = ch9_data.SOC90.DC1;
V_measured = dc_data.V;
I_measured = dc_data.I;
t_measured = seconds(dc_data.t);

% 작은 샘플로 테스트 (처음 1000개 포인트만)
sample_size = 1000;
V_sample = V_measured(1:sample_size);
I_sample = I_measured(1:sample_size);
t_sample = t_measured(1:sample_size);

dt = 0.1;  % Fixed time step as specified
fprintf('Sample size: %d, dt: %.3f s\n', sample_size, dt);

% 초기 조건 추정
[SOC_initial, battery_capacity, validation_result] = estimate_and_validate_soc(V_sample, I_sample, t_sample, OCV_integrated_RPT0cyc);

% SOC 계산
SOC_sample = zeros(sample_size, 1);
SOC_sample(1) = SOC_initial;
for i = 2:sample_size
    SOC_sample(i) = SOC_sample(i-1) - (I_sample(i) * dt) / (battery_capacity * 3600) * 100;
end

fprintf('SOC range: %.2f ~ %.2f%%\n', min(SOC_sample), max(SOC_sample));

% 초기 파라미터 테스트
params_initial = [0.003, 0.005, 1500, 0.008, 8000];
fprintf('\n초기 파라미터: [%.6f, %.6f, %.1f, %.6f, %.1f]\n', params_initial);

% 1. OCV 테스트
fprintf('\n=== OCV 테스트 ===\n');
try
    for i = [1, sample_size/2, sample_size]
        i = round(i);
        soc_val = SOC_sample(i);
        if soc_val >= 0 && soc_val <= 100
            ocv_val = ocv_func(soc_val);
            fprintf('SOC(%.0f) = %.2f%% -> OCV = %.4f V\n', i, soc_val, ocv_val);
        else
            fprintf('SOC(%.0f) = %.2f%% -> 범위 벗어남!\n', i, soc_val);
        end
    end
catch ME
    fprintf('OCV 계산 에러: %s\n', ME.message);
end

% 2. ECM 모델 테스트
fprintf('\n=== ECM 모델 테스트 ===\n');
try
    V_model = ECM_2RC_model(params_initial, I_sample, SOC_sample, t_sample, ocv_func);
    fprintf('ECM 모델 계산 성공!\n');
    fprintf('V_model range: %.4f ~ %.4f V\n', min(V_model), max(V_model));
    fprintf('V_measured range: %.4f ~ %.4f V\n', min(V_sample), max(V_sample));
catch ME
    fprintf('ECM 모델 에러: %s\n', ME.message);
    return;
end

% 3. Cost function 테스트
fprintf('\n=== Cost Function 테스트 ===\n');
try
    cost_val = cost_function(params_initial, I_sample, SOC_sample, V_sample, t_sample, ocv_func);
    fprintf('Cost function 계산 성공!\n');
    fprintf('Initial cost (RMSE): %.6f V\n', cost_val);
catch ME
    fprintf('Cost function 에러: %s\n', ME.message);
    return;
end

% 4. 파라미터 범위 테스트
fprintf('\n=== 파라미터 범위 테스트 ===\n');
lb = [0.0005, 0.0005, 50, 0.0005, 100];
ub = [0.05, 0.05, 20000, 0.05, 50000];

for i = 1:length(params_initial)
    if params_initial(i) < lb(i) || params_initial(i) > ub(i)
        fprintf('파라미터 %d (%.6f)가 범위 [%.6f, %.6f] 벗어남!\n', i, params_initial(i), lb(i), ub(i));
    else
        fprintf('파라미터 %d (%.6f) 범위 OK\n', i, params_initial(i));
    end
end

fprintf('\n=== 결론 ===\n');
if exist('cost_val', 'var') && isfinite(cost_val)
    fprintf('모든 테스트 통과! Cost = %.6f\n', cost_val);
else
    fprintf('문제 발견됨!\n');
end 