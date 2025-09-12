clear; clc; close all;

% 시뮬레이션 파라미터
true_params = [0.002, 0.001, 0.0005, 10, 100]; % [R0, R1, R2, tau1, tau2]
dt = 1; % 1초 간격
N = 1000; % 1000초 시뮬레이션

% 전류 프로파일 생성 (step + pulse)
I = zeros(N, 1);
I(100:400) = -2; % 방전
I(500:600) = 1;  % 충전
I(700:800) = -3; % 강한 방전

% SOC 계산 (임의의 초기값)
SOC = zeros(N, 1);
SOC(1) = 0.5; % 50%에서 시작
capacity = 64; % Ah

for k = 2:N
    SOC(k) = SOC(k-1) + I(k-1) * dt / (capacity * 3600);
end

% OCV 함수 (선형 근사)
V_OCV = 3.0 + 1.2 * SOC;

% 실제 과전압 계산 (2RC 모델)
V_overpotential_true = calculate_overpotential(I, true_params, dt);

% 실제 터미널 전압
V_terminal_true = V_OCV + V_overpotential_true;

% 노이즈 추가
noise_level = 0.005; % 5mV
V_terminal_measured = V_terminal_true + noise_level * randn(N, 1);

% OCV 불확실성 (실제로는 OCV도 불확실함)
OCV_uncertainty = 0.002; % 2mV
V_OCV_estimated = V_OCV + OCV_uncertainty * randn(N, 1);

% 과전압 계산 (OCV 불확실성 포함)
V_overpotential_measured = V_terminal_measured - V_OCV_estimated;

%% 방법 1: 과전압 피팅
fprintf('=== 방법 1: 과전압 피팅 ===\n');

% 초기 추정값
p0 = [0.0025, 0.0012, 0.0006, 12, 120];
lb = [0, 0, 0, 1, 10];
ub = [0.01, 0.01, 0.01, 100, 1000];

% 최적화 옵션 (aligned with BSG_Reference)
options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 3000, 'StepTolerance', 1e-15);

% 과전압 피팅 - createOptimProblem 사용
cost_func_overpotential = @(p) cost_overpotential_func(p, V_overpotential_measured, I, dt);
problem_overpotential = createOptimProblem('fmincon', ...
    'objective', cost_func_overpotential, ...
    'x0', p0, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

tic;
[params_overpotential, cost_overpotential] = fmincon(problem_overpotential);
time_overpotential = toc;

fprintf('추정된 파라미터: [%.4f, %.4f, %.4f, %.1f, %.1f]\n', params_overpotential);
fprintf('실제 파라미터:   [%.4f, %.4f, %.4f, %.1f, %.1f]\n', true_params);
fprintf('RMSE: %.5f V\n', cost_overpotential);
fprintf('소요 시간: %.3f 초\n\n', time_overpotential);

%% 방법 2: 전체 터미널 전압 피팅
fprintf('=== 방법 2: 전체 터미널 전압 피팅 ===\n');

% 전체 터미널 전압 피팅 - createOptimProblem 사용
cost_func_terminal = @(p) cost_terminal_func(p, V_terminal_measured, I, V_OCV, dt);
problem_terminal = createOptimProblem('fmincon', ...
    'objective', cost_func_terminal, ...
    'x0', p0, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

tic;
[params_terminal, cost_terminal] = fmincon(problem_terminal);
time_terminal = toc;

fprintf('추정된 파라미터: [%.4f, %.4f, %.4f, %.1f, %.1f]\n', params_terminal);
fprintf('실제 파라미터:   [%.4f, %.4f, %.4f, %.1f, %.1f]\n', true_params);
fprintf('RMSE: %.5f V\n', cost_terminal);
fprintf('소요 시간: %.3f 초\n\n', time_terminal);

%% 결과 비교
fprintf('=== 파라미터 추정 오차 비교 ===\n');
error_overpotential = abs(params_overpotential - true_params) ./ true_params * 100;
error_terminal = abs(params_terminal - true_params) ./ true_params * 100;

param_names = {'R0', 'R1', 'R2', 'tau1', 'tau2'};
fprintf('파라미터\t과전압 방법\t터미널 전압 방법\n');
fprintf('-------\t-------\t-------\n');
for i = 1:5
    fprintf('%s\t\t%.1f%%\t\t%.1f%%\n', param_names{i}, error_overpotential(i), error_terminal(i));
end

%% 시각화
figure('Position', [100, 100, 1200, 800]);

% 전류 프로파일
subplot(3, 2, 1);
plot(1:N, I, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Current (A)');
title('Current Profile');
grid on;

% SOC 변화
subplot(3, 2, 2);
plot(1:N, SOC*100, 'g-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('SOC (%)');
title('SOC Profile');
grid on;

% 터미널 전압 비교
subplot(3, 2, 3);
plot(1:N, V_terminal_true, 'k-', 'LineWidth', 2); hold on;
plot(1:N, V_terminal_measured, 'r:', 'LineWidth', 1);
V_terminal_model_overpotential = V_OCV + calculate_overpotential(I, params_overpotential, dt);
V_terminal_model_terminal = V_OCV + calculate_overpotential(I, params_terminal, dt);
plot(1:N, V_terminal_model_overpotential, 'b--', 'LineWidth', 1.5);
plot(1:N, V_terminal_model_terminal, 'm--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Terminal Voltage (V)');
title('Terminal Voltage Comparison');
legend('True', 'Measured', 'Overpotential Method', 'Terminal Method', 'Location', 'best');
grid on;

% 과전압 비교
subplot(3, 2, 4);
plot(1:N, V_overpotential_true, 'k-', 'LineWidth', 2); hold on;
plot(1:N, V_overpotential_measured, 'r:', 'LineWidth', 1);
V_overpotential_model_overpotential = calculate_overpotential(I, params_overpotential, dt);
V_overpotential_model_terminal = calculate_overpotential(I, params_terminal, dt);
plot(1:N, V_overpotential_model_overpotential, 'b--', 'LineWidth', 1.5);
plot(1:N, V_overpotential_model_terminal, 'm--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Overpotential (V)');
title('Overpotential Comparison');
legend('True', 'Measured', 'Overpotential Method', 'Terminal Method', 'Location', 'best');
grid on;

% 에러 분석
subplot(3, 2, 5);
error_overpotential_fit = V_overpotential_measured - V_overpotential_model_overpotential;
error_terminal_fit = V_terminal_measured - V_terminal_model_terminal;
plot(1:N, error_overpotential_fit, 'b-', 'LineWidth', 1); hold on;
plot(1:N, error_terminal_fit, 'r-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Fitting Error (V)');
title('Fitting Error Comparison');
legend('Overpotential Method', 'Terminal Method', 'Location', 'best');
grid on;

% 파라미터 추정 오차
subplot(3, 2, 6);
x = 1:5;
bar(x-0.2, error_overpotential, 0.4, 'b'); hold on;
bar(x+0.2, error_terminal, 0.4, 'r');
set(gca, 'XTickLabel', param_names);
xlabel('Parameters');
ylabel('Relative Error (%)');
title('Parameter Estimation Error');
legend('Overpotential Method', 'Terminal Method', 'Location', 'best');
grid on;

sgtitle('Comparison: Overpotential vs Terminal Voltage Fitting');

%% 함수 정의
function V_overpotential = calculate_overpotential(I, params, dt)
    R0 = params(1);
    R1 = params(2);
    R2 = params(3);
    tau1 = params(4);
    tau2 = params(5);
    
    N = length(I);
    V_overpotential = zeros(N, 1);
    V_C1 = 0;
    V_C2 = 0;
    
    for k = 1:N
        % R0 전압강하
        V_R0 = R0 * I(k);
        
        % RC 전압 업데이트
        if k > 1
            alpha1 = exp(-dt/tau1);
            alpha2 = exp(-dt/tau2);
            V_C1 = V_C1 * alpha1 + R1 * (1 - alpha1) * I(k-1);
            V_C2 = V_C2 * alpha2 + R2 * (1 - alpha2) * I(k-1);
        end
        
        % 총 과전압
        V_overpotential(k) = V_R0 + V_C1 + V_C2;
    end
end

function cost = cost_overpotential_func(params, V_overpotential_measured, I, dt)
    V_overpotential_model = calculate_overpotential(I, params, dt);
    cost = sqrt(mean((V_overpotential_measured - V_overpotential_model).^2));
end

function cost = cost_terminal_func(params, V_terminal_measured, I, V_OCV, dt)
    V_overpotential_model = calculate_overpotential(I, params, dt);
    V_terminal_model = V_OCV + V_overpotential_model;
    cost = sqrt(mean((V_terminal_measured - V_terminal_model).^2));
end 