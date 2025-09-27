clear; clc; close all;

current_channel = 'ch9';  % Type Channel

%% Configuration
% Define processing conditions
cycle_conditions = {'cyc0', 'cyc200'};
soc_conditions = {'SOC70'}; % , 'SOC70', 'SOC90'};

% Base paths
base_path = 'G:\공유 드라이브\Battery Software Lab\Projects\KEPCO_ATB_Lab\ESS_Data_Preprocessing\Experiment_DriveCycle_Fitting';
results_path = fullfile(base_path, 'Results');

% Create results directory
if ~exist(results_path, 'dir')
    mkdir(results_path);
    fprintf('Created Results directory: %s\n', results_path);
end

% Drive cycle paths
drive_cycle_paths = struct();
drive_cycle_paths.('cyc0') = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\parsedDriveCycle\parsedDriveCycle_0cyc_filtered.mat';
drive_cycle_paths.('cyc200') = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\parsedDriveCycle\parsedDriveCycle_200cyc_filtered.mat';

% OCV path (contains both 0cyc and 200cyc OCV functions)
ocv_path = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\RPT\Summary\OCV_interpolation_functions.mat';

% Load OCV functions
load(ocv_path);
fprintf('OCV interpolation functions loaded successfully\n');

% Color settings
color_I   = '#CD534C';
color_V   = '#0073C2';
color_Vest= '#925E9F';
color_soc = '#20854E';
color_err = '#EE4C97';

%% Parameter Setting
% SOC별 초기 파라미터 설정
X_SOC50 = [0.001, 0.003, 0.003, 10, 100];     % SOC50
X_SOC70 = [0.0015, 0.004, 0.004, 15, 150];    % SOC70  
X_SOC90 = [0.002, 0.005, 0.005, 20, 200];     % SOC90

% 기본값 (나중에 SOC별로 업데이트됨)
para0 = X_SOC50;

lb = [0.0001, 0.0001, 0.0001, 1, 1];                 % lower
ub = [para0(1:3)*1000, para0(4)*50, para0(5)*100]; % upper 

% (tau1 < tau2)
% A * params <= b, where params = [R0, R1, R2, tau1, tau2]
% tau1 - tau2 <= 0  →  [0, 0, 0, 1, -1] * params <= 0
A_lin = [0, 0, 0, 1, -1];   % 1*tau1 + (-1)*tau2 <= 0
b_lin = 0;

fprintf('Parameter bounds:\n');
fprintf('  R0: %.1f - %.1f mΩ\n', lb(1)*1000, ub(1)*1000);
fprintf('  R1: %.1f - %.1f mΩ\n', lb(2)*1000, ub(2)*1000);
fprintf('  R2: %.1f - %.1f mΩ\n', lb(3)*1000, ub(3)*1000);
fprintf('  τ1: %.1f - %.1f s\n', lb(4), ub(4));
fprintf('  τ2: %.1f - %.1f s\n', lb(5), ub(5));

fprintf('Linear constraint: τ1 < τ2 \n');
fprintf('Initial value: R0=%.1f, R1=%.1f, R2=%.1f mΩ, τ1=%.0f, τ2=%.0f s\n', ...
    para0(1)*1000, para0(2)*1000, para0(3)*1000, para0(4), para0(5));

% Optimization options
options = optimset( ...
    'Display',     'off', ...
    'MaxIter',     3000, ...
    'MaxFunEvals', 1e6, ...
    'TolFun',      1e-16, ...
    'TolX',        1e-16);

%% Main Processing Loop
fprintf('\n=== Starting ECM 2RC Fitting for Multiple Conditions ===\n');

for current_cycle = cycle_conditions
    current_cycle = current_cycle{1};  % Extract string from cell array
    
    fprintf('\n=== Processing %s Condition ===\n', current_cycle);
    
    % Load drive cycle data for current cycle
    drive_cycle_path = drive_cycle_paths.(current_cycle);
    load(drive_cycle_path);
    
    % Select appropriate OCV function
    if strcmp(current_cycle, 'cyc0')
        ocv_func = OCV_integrated_RPT0cyc.OCV_SOC_func;
        channel_data = parsedDriveCycle_0cyc.(sprintf('%s_Drive_0cyc', current_channel));
    else  % cyc200
        ocv_func = OCV_integrated_RPT200cyc.OCV_SOC_func;
        channel_data = parsedDriveCycle_200cyc.(sprintf('%s_Drive_200cyc', current_channel));
    end
    
    fprintf('Loaded %s drive cycle data and %s OCV function\n', current_cycle, current_cycle);
    
    % Create cycle results directory
    cycle_results_path = fullfile(results_path, current_cycle);
    if ~exist(cycle_results_path, 'dir')
        mkdir(cycle_results_path);
    end
    
    % Create channel results directory
    channel_results_path = fullfile(cycle_results_path, current_channel);
    if ~exist(channel_results_path, 'dir')
        mkdir(channel_results_path);
    end
    
    for current_soc = soc_conditions
        current_soc = current_soc{1};  % Extract string from cell array
        
        fprintf('\n=== Processing %s ===\n', current_soc);
        
        % SOC별 초기 파라미터 설정
        if strcmp(current_soc, 'SOC50')
            para0 = X_SOC50;
            fprintf('Using SOC50 parameters: R0=%.1f, R1=%.1f, R2=%.1f mΩ, τ1=%.0f, τ2=%.0f s\n', ...
                para0(1)*1000, para0(2)*1000, para0(3)*1000, para0(4), para0(5));
        elseif strcmp(current_soc, 'SOC70')
            para0 = X_SOC70;
            fprintf('Using SOC70 parameters: R0=%.1f, R1=%.1f, R2=%.1f mΩ, τ1=%.0f, τ2=%.0f s\n', ...
                para0(1)*1000, para0(2)*1000, para0(3)*1000, para0(4), para0(5));
        elseif strcmp(current_soc, 'SOC90')
            para0 = X_SOC90;
            fprintf('Using SOC90 parameters: R0=%.1f, R1=%.1f, R2=%.1f mΩ, τ1=%.0f, τ2=%.0f s\n', ...
                para0(1)*1000, para0(2)*1000, para0(3)*1000, para0(4), para0(5));
        end
        
        % Create SOC results directory
        soc_results_path = fullfile(channel_results_path, current_soc);
        if ~exist(soc_results_path, 'dir')
            mkdir(soc_results_path);
        end
        
        % Initialize SOC results structure
        soc_results = struct();
        
        % Process DC1-DC8
        for dc_idx = 1:8
            DC_name = sprintf('DC%d', dc_idx);
            
            % Check if DC data exists
            if ~isfield(channel_data, current_soc) || ~isfield(channel_data.(current_soc), DC_name)
                fprintf('Skipping %s - not found in %s\n', DC_name, current_soc);
                continue;
            end
            
            fprintf('\n=== Processing %s ===\n', DC_name);
            
            % Create DC results directory
            dc_results_path = fullfile(soc_results_path, DC_name);
            if ~exist(dc_results_path, 'dir')
                mkdir(dc_results_path);
            end
            
            % Extract current DC data
            dc_data = channel_data.(current_soc).(DC_name);
            V_measured = dc_data.V;
            I_measured = dc_data.I;
            t_measured = seconds(dc_data.t);  % Convert duration to seconds
            
            % Calculate time interval and data size
            dt = 0.1;  % Fixed time step as specified
            N = length(t_measured);
            
            fprintf('Data points: %d, Time step: %.3f s\n', N, dt);
            fprintf('Time range: %.1f to %.1f s (%.1f minutes)\n', min(t_measured), max(t_measured), (max(t_measured) - min(t_measured))/60);
            fprintf('Current range: %.3f to %.3f A\n', min(I_measured), max(I_measured));
            fprintf('Voltage range: %.3f to %.3f V\n', min(V_measured), max(V_measured));
            
            % Calculate SOC using estimate_and_validate_soc function
            if strcmp(current_cycle, 'cyc0')
                [SOC_initial, battery_capacity, validation_result, initial_rest_end, SOC, final_rest_end] = estimate_and_validate_soc(V_measured, I_measured, t_measured, OCV_integrated_RPT0cyc);
            else  % cyc200
                [SOC_initial, battery_capacity, validation_result, initial_rest_end, SOC, final_rest_end] = estimate_and_validate_soc(V_measured, I_measured, t_measured, OCV_integrated_RPT200cyc);
            end
            
            % Validation result available
            fprintf(' SOC validation: Initial %.2f%%, Final %.2f%%\n', validation_result.initial_soc, validation_result.final_soc);
            
            % Debug: Check SOC range
            fprintf(' SOC range: %.2f%% to %.2f%%\n', min(SOC), max(SOC));
            fprintf(' SOC std: %.3f%%\n', std(SOC));
            
            %% KEY FEATURE: Calculate OCV and Overpotential for Analysis
            fprintf('\n=== Calculating OCV and Overpotential ===\n');
            
            % Calculate OCV at each SOC point
            V_OCV = ocv_func(SOC);
            
            % Calculate overpotential = V_measured - V_OCV (for analysis)
            V_overpotential = V_measured - V_OCV;
            
            fprintf('OCV range: %.3f to %.3f V\n', min(V_OCV), max(V_OCV));
            fprintf('Overpotential range: %.3f to %.3f V\n', min(V_overpotential), max(V_overpotential));
            
            % Show some examples (commented out)
            % fprintf('\nExample data points:\n');
            % for i = 1:min(5, length(t_measured))
            %    fprintf('  t=%.1fs: I=%.3fA, V=%.3fV, OCV=%.3fV, Overpotential=%.3fV, SOC=%.1f%%\n', ...
            %        t_measured(i), I_measured(i), V_measured(i), V_OCV(i), V_overpotential(i), SOC(i));
            % end
            
            %% Parameter Optimization for Terminal Voltage (OCV + Overpotential)
            fprintf('\n=== Parameter Optimization for Terminal Voltage (%s) ===\n', DC_name);
            
            % Initial condition matching
            fprintf('Initial overpotential: measured = %.4f V\n', V_overpotential(1));
            
            % Create cost function handle for terminal voltage (OCV + overpotential)
            cost_handle = @(params) cost_function(params, I_measured, SOC, V_measured, t_measured, ocv_func);
            
            % Initial state output
            fprintf('Initial parameters: R0=%.6f, R1=%.6f, R2=%.6f, tau1=%.2f, tau2=%.2f\n', ...
                para0(1), para0(2), para0(3), para0(4), para0(5));
            
            initial_cost = cost_handle(para0);
            fprintf('Initial cost (RMSE): %.6f V\n', initial_cost);
            
            % MultiStart 사용
            fprintf('\n=== Starting optimization with MultiStart ===\n');
            
            total_tic = tic;
            
            % MultiStart 설정
            ms = MultiStart( ...
                "UseParallel", true, ...
                "Display", "iter");
            nStartPts = 100;  
            
            % 초기값 생성
            startPoints = zeros(nStartPts, 5);
            
            % 1. 기본 초기값 (para0)
            startPoints(1, :) = para0;
            
            % 2. 시나리오 생성
            for i = 2:nStartPts
                scenario = mod(i-2, 8) + 1;  % 8가지 시나리오 반복
                
                switch scenario
                    case 1  % 저저항 시나리오
                        r0_rand = para0(1) * (0.5 + 0.3*rand());
                        r1_rand = para0(2) * (0.3 + 0.4*rand());
                        r2_rand = para0(3) * (0.3 + 0.4*rand());
                        tau1_rand = para0(4) * (0.8 + 0.4*rand());
                        tau2_rand = para0(5) * (0.8 + 0.6*rand());
                        
                    case 2  % 고저항 시나리오
                        r0_rand = para0(1) * (1.5 + 0.5*rand());
                        r1_rand = para0(2) * (1.2 + 0.8*rand());
                        r2_rand = para0(3) * (1.2 + 0.8*rand());
                        tau1_rand = para0(4) * (0.6 + 0.4*rand());
                        tau2_rand = para0(5) * (0.6 + 0.6*rand());
                        
                    case 3  % 빠른 응답 시나리오 (작은 tau)
                        r0_rand = para0(1) * (0.8 + 0.4*rand());
                        r1_rand = para0(2) * (0.6 + 0.4*rand());
                        r2_rand = para0(3) * (0.6 + 0.4*rand());
                        tau1_rand = para0(4) * (0.3 + 0.3*rand());
                        tau2_rand = para0(5) * (0.4 + 0.4*rand());
                        
                    case 4  % 느린 응답 시나리오 (큰 tau)
                        r0_rand = para0(1) * (0.8 + 0.4*rand());
                        r1_rand = para0(2) * (0.8 + 0.4*rand());
                        r2_rand = para0(3) * (0.8 + 0.4*rand());
                        tau1_rand = para0(4) * (1.5 + 0.5*rand());
                        tau2_rand = para0(5) * (2.0 + 1.0*rand());
                        
                    case 5  % 균형 시나리오
                        r0_rand = para0(1) * (0.8 + 0.4*rand());
                        r1_rand = para0(2) * (0.8 + 0.4*rand());
                        r2_rand = para0(3) * (0.8 + 0.4*rand());
                        tau1_rand = para0(4) * (0.8 + 0.4*rand());
                        tau2_rand = para0(5) * (1.0 + 0.6*rand());
                        
                    case 6  % R1 우세 시나리오
                        r0_rand = para0(1) * (0.6 + 0.3*rand());
                        r1_rand = para0(2) * (1.5 + 0.5*rand());
                        r2_rand = para0(3) * (0.5 + 0.3*rand());
                        tau1_rand = para0(4) * (1.0 + 0.5*rand());
                        tau2_rand = para0(5) * (0.8 + 0.4*rand());
                        
                    case 7  % R2 우세 시나리오
                        r0_rand = para0(1) * (0.6 + 0.3*rand());
                        r1_rand = para0(2) * (0.5 + 0.3*rand());
                        r2_rand = para0(3) * (1.5 + 0.5*rand());
                        tau1_rand = para0(4) * (0.8 + 0.4*rand());
                        tau2_rand = para0(5) * (1.2 + 0.6*rand());
                        
                    case 8  % 극단적 시나리오 (경계값 근처)
                        r0_rand = para0(1) * (0.1 + 0.2*rand());
                        r1_rand = para0(2) * (0.1 + 0.2*rand());
                        r2_rand = para0(3) * (0.1 + 0.2*rand());
                        tau1_rand = para0(4) * (0.2 + 0.3*rand());
                        tau2_rand = para0(5) * (0.3 + 0.4*rand());
                end
                
                % tau1 < tau2 제약 확인 및 조정
                if tau1_rand >= tau2_rand
                    tau2_rand = tau1_rand * (1.2 + 0.3*rand());
                end
                
                % 물리적 제약 확인
                if r0_rand < lb(1), r0_rand = lb(1) * (1 + 0.1*rand()); end
                if r1_rand < lb(2), r1_rand = lb(2) * (1 + 0.1*rand()); end
                if r2_rand < lb(3), r2_rand = lb(3) * (1 + 0.1*rand()); end
                if tau1_rand < lb(4), tau1_rand = lb(4) * (1 + 0.1*rand()); end
                if tau2_rand < lb(5), tau2_rand = lb(5) * (1 + 0.1*rand()); end
                
                startPoints(i, :) = [r0_rand, r1_rand, r2_rand, tau1_rand, tau2_rand];
            end
            
            startPts = CustomStartPointSet(startPoints);
            
            fprintf('Using MultiStart with %d smart starting points\n', nStartPts);
            
            % Create optimization problem
            problem = createOptimProblem('fmincon', ...
                'objective', cost_handle, ...
                'x0', para0, ...
                'lb', lb, ...
                'ub', ub, ...
                'Aineq', A_lin, ...
                'bineq', b_lin, ...
                'options', options);
            
            % Run MultiStart optimization
            [params_optimal, cost_optimal, exitflag, output, solutions] = run(ms, problem, startPts);
            
            % 결과 변환
            all_results = [];
            for i = 1:length(solutions)
                all_results(i).params = solutions(i).X;
                all_results(i).cost = solutions(i).Fval;
                all_results(i).converged = (solutions(i).Exitflag > 0);
                all_results(i).initial_params = [];  
            end
            
            %  출력
            fprintf('\n=== MultiStart Optimization Results ===\n');
            fprintf('Best solution found:\n');
            fprintf('  RMSE: %.6f V\n', cost_optimal);
            fprintf('  Exit flag: %d\n', exitflag);
            fprintf('  Parameters: [%.6g %.6g %.6g %.6g %.6g]\n', params_optimal(1), params_optimal(2), params_optimal(3), params_optimal(4), params_optimal(5));
            
            % Show top 5 solutions
            fprintf('\nTop 5 solutions:\n');
            for i = 1:min(5, length(solutions))
                fprintf('Solution %d: RMSE=%.6f V, exitflag=%d, params=[%.6g %.6g %.6g %.6g %.6g]\n', ...
                    i, solutions(i).Fval, solutions(i).Exitflag, solutions(i).X(1), solutions(i).X(2), solutions(i).X(3), solutions(i).X(4), solutions(i).X(5));
            end
            
            total_time = toc(total_tic);
            
            %% Process Multiple Solutions for Statistics
            fprintf('\n=== Processing Multiple Solutions ===\n');
            
            % MultiStart 결과 사용
            num_runs = length(solutions);
            
            fprintf('\n=== Optimization completed for %s ===\n', DC_name);
            fprintf('Total runs: %d\n', num_runs);
            fprintf('Optimization with MultiStart and tau1 < tau2 constraint\n');
            fprintf('Storing %d solutions for statistical analysis\n', num_runs);
            
            % Store individual runs
            runs_data = all_results;
            
            % Calculate statistics
            all_params = vertcat(runs_data.params);  % Nx5 matrix
            all_costs = [runs_data.cost];
            converged_mask = [runs_data.converged];
            
            statistics = struct( ...
                'mean_params', mean(all_params, 1), ...
                'std_params', std(all_params, 0, 1), ...
                'min_params', min(all_params, [], 1), ...
                'max_params', max(all_params, [], 1), ...
                'mean_cost', mean(all_costs), ...
                'std_cost', std(all_costs), ...
                'min_cost', min(all_costs), ...
                'max_cost', max(all_costs), ...
                'num_runs', num_runs, ...
                'num_converged', sum(converged_mask), ...
                'convergence_rate', sum(converged_mask) / num_runs ...
                );
            
            % Display statistics
            fprintf('\n--- Statistics Summary ---\n');
            fprintf('Convergence rate: %.1f%% (%d/%d)\n', statistics.convergence_rate*100, statistics.num_converged, num_runs);
            fprintf('RMSE range: %.6f - %.6f V (std: %.6f V)\n', statistics.min_cost, statistics.max_cost, statistics.std_cost);
            
            param_names = {'R0', 'R1', 'R2', 'τ1', 'τ2'};
            param_units = {'Ω', 'Ω', 'Ω', 's', 's'};
            for p = 1:5
                fprintf('%s: %.6f ± %.6f %s (range: %.6f - %.6f %s)\n', ...
                    param_names{p}, statistics.mean_params(p), statistics.std_params(p), param_units{p}, ...
                    statistics.min_params(p), statistics.max_params(p), param_units{p});
            end
            
            % Result validation
            if cost_optimal < initial_cost
                fprintf('✓ Optimization successful: cost reduced from %.6f to %.6f V\n', initial_cost, cost_optimal);
            else
                fprintf('⚠ Optimization may have failed: cost increased from %.6f to %.6f V\n', initial_cost, cost_optimal);
            end
            
            % Extract optimal parameters
            R0 = params_optimal(1);
            R1 = params_optimal(2);
            R2 = params_optimal(3);
            tau1 = params_optimal(4);
            tau2 = params_optimal(5);
            
            % Calculate capacitors
            C1 = (tau1 / R1) / 1000;
            C2 = (tau2 / R2) / 1000;
            
            fprintf('Optimization completed in %.2f seconds\n', total_time);
            fprintf('Best RMSE achieved: %.6f V\n', cost_optimal);
            
            % Detailed parameter output
            fprintf('\n=== Optimized Parameters for %s ===\n', DC_name);
            fprintf('R0 = %.6f Ω (%.1f mΩ)\n', R0, R0*1000);
            fprintf('R1 = %.6f Ω (%.1f mΩ)\n', R1, R1*1000);
            fprintf('R2 = %.6f Ω (%.1f mΩ)\n', R2, R2*1000);
            fprintf('τ1 = %.2f s\n', tau1);
            fprintf('τ2 = %.2f s\n', tau2);
            fprintf('C1 = %.2f kF\n', C1);
            fprintf('C2 = %.2f kF\n', C2);
            
            % Parameter sanity check
            fprintf('\n--- Parameter Sanity Check ---\n');
            
            % Check bounds
            if R0 >= ub(1)*0.95, fprintf('   WARNING: R0 near upper bound (%.1f mΩ)\n', ub(1)*1000); end
            if R1 >= ub(2)*0.95, fprintf('   WARNING: R1 near upper bound (%.1f mΩ)\n', ub(2)*1000); end
            if R2 >= ub(3)*0.95, fprintf('   WARNING: R2 near upper bound (%.1f mΩ)\n', ub(3)*1000); end
            if tau1 >= ub(4)*0.95, fprintf('   WARNING: τ1 near upper bound (%.0f s)\n', ub(4)); end
            if tau2 >= ub(5)*0.95, fprintf('   WARNING: τ2 near upper bound (%.0f s = %.1f min)\n', ub(5), ub(5)/60); end
            
            if R0 <= lb(1)*1.05, fprintf('   WARNING: R0 near lower bound (%.1f mΩ)\n', lb(1)*1000); end
            if R1 <= lb(2)*1.05, fprintf('   WARNING: R1 near lower bound (%.1f mΩ)\n', lb(2)*1000); end
            if R2 <= lb(3)*1.05, fprintf('   WARNING: R2 near lower bound (%.1f mΩ)\n', lb(3)*1000); end
            if tau1 <= lb(4)*1.05, fprintf('   WARNING: τ1 near lower bound (%.0f s)\n', lb(4)); end
            if tau2 <= lb(5)*1.05, fprintf('   WARNING: τ2 near lower bound (%.0f s)\n', lb(5)); end
            
            % Check capacitance values (STRICTER LIMITS)
            if C1 > 50000, fprintf('   WARNING: C1 = %.0f F is unusually large (>50kF)\n', C1); end
            if C2 > 200000, fprintf('   CRITICAL: C2 = %.0f F is extremely large (>200kF) - Check τ2 and R2!\n', C2); end
            if C1 < 10, fprintf('   WARNING: C1 = %.2f F is unusually small (<10F)\n', C1); end
            if C2 < 100, fprintf('   WARNING: C2 = %.2f F is unusually small (<100F)\n', C2); end
            
            % Additional physical constraint checks
            if tau2 > 1800, fprintf('   CRITICAL: τ2 = %.0f s (%.1f min) exceeds 30 minutes!\n', tau2, tau2/60); end
            if C2 > 10*C1, fprintf('   WARNING: C2/C1 ratio = %.1f is very high (>10)\n', C2/C1); end
            
            % ECM parameter relationships
            tau_ratio = tau2/tau1;
            if tau_ratio > 100, fprintf('   WARNING: τ2/τ1 ratio = %.1f is very high (>100)\n', tau_ratio); end
            if tau_ratio < 2, fprintf('   WARNING: τ2/τ1 ratio = %.1f is too small (<2)\n', tau_ratio); end
            
            % Check time constant ordering
            if tau1 >= tau2
                fprintf('   ERROR: τ1 ≥ τ2 (%.2f ≥ %.2f) - VIOLATES CONSTRAINT!\n', tau1, tau2);
            else
                fprintf('   ✓ CONSTRAINT SATISFIED: τ1 < τ2 (%.2f < %.2f)\n', tau1, tau2);
            end
            
            % Check resistance ordering (typical: R0 < R1, R2)
            if R0 > R1, fprintf('   WARNING: R0 > R1 (unusual for ECM)\n'); end
            if R0 > R2, fprintf('   WARNING: R0 > R2 (unusual for ECM)\n'); end
            
            fprintf('Parameter check completed.\n');
            
            %% Model Validation
            fprintf('\n=== Model Validation ===\n');
            
            % Calculate model terminal voltage (OCV + overpotential)
            V_terminal_model = ECM_2RC_model(params_optimal, I_measured, SOC, t_measured, ocv_func);
            
            % Calculate model overpotential for comparison
            V_overpotential_model = V_terminal_model - V_OCV;
            
            % Calculate performance metrics
            overpotential_error = V_overpotential - V_overpotential_model;
            terminal_error = V_measured - V_terminal_model;
            
            overpotential_rmse = sqrt(mean(overpotential_error.^2));
            terminal_rmse = sqrt(mean(terminal_error.^2));
            overpotential_mae = mean(abs(overpotential_error));
            terminal_mae = mean(abs(terminal_error));
            
            fprintf('Terminal voltage fitting performance (primary metric):\n');
            fprintf('  RMSE: %.6f V\n', terminal_rmse);
            fprintf('  MAE:  %.6f V\n', terminal_mae);
            fprintf('  Max error: %.6f V\n', max(abs(terminal_error)));
            
            fprintf('Overpotential prediction performance (secondary):\n');
            fprintf('  RMSE: %.6f V\n', overpotential_rmse);
            fprintf('  MAE:  %.6f V\n', overpotential_mae);
            fprintf('  Max error: %.6f V\n', max(abs(overpotential_error)));
            
            %% 피팅 결과 시각화 (Cost Surface 생성 전에 바로 표시)
            fprintf('\n=== Creating Fitting Results Visualization ===\n');
            
            % Use dedicated visualization function with OCV
            h_fitting = visualize_fitting_results(DC_name, t_measured, I_measured, SOC, V_measured, V_terminal_model, V_OCV, terminal_error, terminal_rmse, terminal_mae, color_I, color_V, color_soc);
            drawnow; % 그래프 완성 보장
            
            % Save figure to new Results structure
            fig_file = fullfile(dc_results_path, sprintf('%s_%s_%s_fitting_result.fig', current_cycle, current_soc, DC_name));
            savefig(h_fitting, fig_file);
            fprintf('Fitting result figure saved: %s\n', fig_file);
            
            % Save results to new Results structure
            dc_result_file = fullfile(dc_results_path, sprintf('%s_%s_%s_result.mat', current_cycle, current_soc, DC_name));
            save(dc_result_file, 'params_optimal', 'cost_optimal', 'exitflag', 'output', 'solutions');
            fprintf('Fitting result .mat saved: %s\n', dc_result_file);
            
            % Close figure to free memory
            close(h_fitting);
            
            %% Cost Surface Visualization
            tau1_center = params_optimal(4);
            tau2_center = params_optimal(5);
            tau1_vec = linspace(tau1_center * 0.5, tau1_center * 1.5, 11);
            tau2_vec = linspace(tau2_center * 0.5, tau2_center * 1.5, 11);
            fprintf('\n=== Creating Cost Surface for %s (centered at optimal) ===\n', DC_name);
            figure;
            param_names = {'R0', 'R1', 'R2', 'τ1', 'τ2'};
            visualize_cost_surface(cost_handle, options, tau1_vec, tau2_vec);
            cost_surface_file = fullfile(dc_results_path, sprintf('%s_%s_%s_cost_surface.fig', current_cycle, current_soc, DC_name));
            saveas(gcf, cost_surface_file);
            fprintf('Cost surface saved to: %s\n', cost_surface_file);
            close(gcf);
            
            %% Store results in SOC results structure
            soc_results.(DC_name) = struct( ...
                'params', params_optimal, ...
                'R0', R0, ...
                'R1', R1, ...
                'R2', R2, ...
                'tau1', tau1, ...
                'tau2', tau2, ...
                'C1', C1, ...
                'C2', C2, ...
                'rmse', terminal_rmse, ...
                'mae', terminal_mae, ...
                'max_error', max(abs(terminal_error)), ...
                'V_model', V_terminal_model, ...
                'V_measured', V_measured, ...
                'I_measured', I_measured, ...
                'SOC', SOC, ...
                't_measured', t_measured, ...
                'V_OCV', V_OCV, ...
                'V_overpotential', V_overpotential, ...
                'initial_cost', initial_cost, ...
                'final_cost', cost_optimal, ...
                'exitflag', exitflag, ...
                'solutions', solutions, ...
                'optimization_time', total_time ...
                );
            
            fprintf('\n=== %s Processing Complete ===\n', DC_name);
        end  % DC loop end
        
        %% Create summary table for current SOC
        fprintf('\n=== Creating Summary Table for %s ===\n', current_soc);
        
        % Create summary data table
        dc_names = fieldnames(soc_results);
        summary_data = table();
        
        for i = 1:length(dc_names)
            dc_name = dc_names{i};
            result = soc_results.(dc_name);
            
            % Add row to summary table
            new_row = table({dc_name}, result.R0*1000, result.R1*1000, result.R2*1000, ...
                result.tau1, result.tau2, result.C1, result.C2, ...
                result.rmse*1000, result.mae*1000, result.max_error*1000, ...
                result.optimization_time, ...
                'VariableNames', {'DC', 'R0_mOhm', 'R1_mOhm', 'R2_mOhm', ...
                'tau1_s', 'tau2_s', 'C1_kF', 'C2_kF', ...
                'RMSE_mV', 'MAE_mV', 'MaxError_mV', 'Time_s'});
            
            summary_data = [summary_data; new_row];
        end
        
        % Save summary table
        summary_file = fullfile(soc_results_path, sprintf('%s_%s_summary_table.mat', current_cycle, current_soc));
        save(summary_file, 'summary_data', 'soc_results');
        fprintf('Summary table saved: %s\n', summary_file);
        
        fprintf('\n=== %s Processing Complete ===\n', current_soc);
    end  % SOC loop end
    
    fprintf('\n=== %s Processing Complete ===\n', current_cycle);
end  % Cycle loop end

fprintf('\n=== All Processing Complete ===\n');