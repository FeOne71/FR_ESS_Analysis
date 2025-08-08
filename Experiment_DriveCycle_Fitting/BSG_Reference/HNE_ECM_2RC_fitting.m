clear; clc; close all
c_mat = lines(9);

% data
load("H:\공유 드라이브\BSL_Data4 (1)\HNE_agedcell_2025_processed\RPT_DCIR_processed\HNE_RPT_fresh_4_1_postprocessing_HPPC.mat")
SOC_array = table2array(NE_OCV_linear(:,"SOC"));
V_array = table2array(NE_OCV_linear(:,"V"));

for i = 1:size(n1C_pulse,1)
SOC_val = cell2mat(n1C_pulse.SOC(i)); 
OCV_vec = interp1(SOC_array,V_array,SOC_val,'linear','extrap');
n1C_pulse.OCV{i} = OCV_vec;
end

para_hats = zeros(size(n1C_pulse,1), 5);
for i_pulse = 1:size(n1C_pulse,1)

    x = n1C_pulse.t{i_pulse,1}-n1C_pulse.t{i_pulse,1}(1);
    y1 = n1C_pulse.V{i_pulse,1}-n1C_pulse.OCV{i_pulse,1}; % dV from OCV
    y2 = n1C_pulse.I{i_pulse,1};
   
    time_diff = n1C_pulse.t{i_pulse,1}(end) - n1C_pulse.t{i_pulse,1}(1);
if 150 <= time_diff && time_diff<190
   tau1_0 = 5; % [sec]
   tau2_0 = 60;
   
elseif time_diff <150
   scale = 180/time_diff;
   tau1_0 = 5/scale; % [sec]
   tau2_0 = 60/scale;
end


figure(1)
subplot(5,2,i_pulse)
%yyaxis left
scatter(x, y1, 15, 'o', 'MarkerEdgeColor', c_mat(1,:), 'MarkerFaceColor', c_mat(1,:), ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);

ylim([1.1*min(y1) 0])



% model visualization
para0 = [abs(y1(1))/abs(y2(1)) abs(y1(end) - y1(1))/abs(y2(1))/2 tau1_0 abs(y1(end) - y1(1))/abs(y2(1))/2 tau2_0]; % initial guess
y_model = func_2RC(x,y2,para0);

hold on
%yyaxis left
plot(x,y_model,'-','Color',c_mat(3,:))

% legend({})


%% fitting CASE 1
% initial guess
    %para0
% bound 
    lb= [0 0 1 0 1];
    ub = para0*3;
% weight
    weight = ones(size(y1)); 

% fitting
        options = optimset('display','iter', 'MaxIter',400, 'MaxFunEvals',1e5, ...
                       'TolFun',1e-10, 'TolX',1e-8, 'FinDiffType','central');
        ms = MultiStart('Display', 'iter', 'UseParallel', true); 
        problem = createOptimProblem('fmincon', ...
        'objective', @(para) func_cost(y1, para, x, y2, ones(size(y1))), ...
        'x0', para0, 'lb', lb, 'ub', ub, 'options', options);
        num_start_points = 100;
        [para_hat, fval, exitflag, output, solutions] = run(ms, problem, num_start_points);
        para_hats(i_pulse, :) = para_hat;
    
% visualize
    y_model_hat = func_2RC(x,y2,para_hat);

%    yyaxis left
    plot(x,y_model_hat,'-','Color',c_mat(4,:))

end

lgd = legend({'Experimental Data', 'Initial Guess', 'Fitted Model'}, ...
    'Orientation', 'horizontal', 'FontSize', 10, 'Box', 'on');

lgd.Position = [0.4, 0.95, 0.2, 0.05];




% model
function y = func_2RC(t,I,para)

R0 = para(1);
R1 = para(2);
tau1 = para(3);
R2 = para(4);
tau2 = para(5);
y = I*R0 + I*R1.*(1-exp(-t/tau1))+I*R2.*(1-exp(-t/tau2));

end


% cost (weight)
function cost = func_cost(y_data,para,t,I,weight)
% this is a cost function to be minimized
y_model = func_2RC(t,I,para);
cost = sqrt(sum(((y_data - y_model).^2) .* weight) / sum(weight)); % RMSE error

end