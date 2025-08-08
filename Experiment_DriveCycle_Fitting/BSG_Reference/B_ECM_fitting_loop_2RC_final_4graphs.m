clear; clc; close all

% hyper parameters
% filename = 'postprocessing_HPPC.mat';
c_mat = lines(9);

% data
load("G:\공유 드라이브\GSP_Data\postprocessing_HPPC.mat")
load("G:\공유 드라이브\GSP_Data\ecm_code\2RC_cost_initial.mat")
SOC_array = table2array(NE_OCV_linear(:,"SOC"));
V_array = table2array(NE_OCV_linear(:,"V"));

for i = 1:size(n1C_pulse,1)
SOC_val = cell2mat(n1C_pulse.SOC(i)); 
OCV_vec = interp1(SOC_array,V_array,SOC_val,'linear','extrap');
n1C_pulse.OCV{i} = OCV_vec;
end

selected_pulses = [1,2,7,8]; 
for i_pulse = 1:size(n1C_pulse,1)

    set(gcf, 'Units', 'centimeters', 'Position', [3, 3, 26, 20]);
    x = n1C_pulse.t{i_pulse,1}-n1C_pulse.t{i_pulse,1}(1);
    % y1 = n1C_pulse.V{i_pulse,1}-n1C_pulse.V_final(i_pulse); % dV from OCV
    y1 = n1C_pulse.V{i_pulse,1}-n1C_pulse.OCV{i_pulse,1}; % dV from OCV
    y2 = n1C_pulse.I{i_pulse,1};  


if ismember(i_pulse, selected_pulses)
figure(1)
grid on;
box on;
subplot_index = find(selected_pulses == i_pulse); 
subplot(2,2,subplot_index);

% subplot(2,2,i_pulse)
%yyaxis left
scatter(x, y1, 10, 'o', 'MarkerEdgeColor', c_mat(1,:), 'MarkerFaceColor', c_mat(1,:), ...
    'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2);
% plot(x, y1,'Color',c_mat(1,:));
ylim([1.1*min(y1) 0])
%yyaxis right
%plot(x, y2,'o')
%ylim([1.1*min(y2) 0])
end

% model visualization
%para0_mat =[para0; 0.001 0.001 65; 0.0013 0.0018*2 60; 0.0013 0.0018 60*1.5; 0.0013*0.6 0.0018*0.5 60];
para0 = [abs(y1(1))/abs(y2(1)) abs(y1(end) - y1(1))/abs(y2(1))/2 opt_tau1(i_pulse) abs(y1(end) - y1(1))/abs(y2(1))/2 opt_tau2(i_pulse)]; % initial guess
%para0 = [0.000823014188098458 0.000877899605928098 5.00196884760841 0.00187505094799395 65.0006027524277];
%para0 = [0.001 0.001 5 0.002 65];
y_model = func_2RC(x,y2,para0);

hold on
%yyaxis left
if ismember(i_pulse, selected_pulses)
figure(1)
subplot_index = find(selected_pulses == i_pulse); 
subplot(2,2,subplot_index);
plot(x,y_model,'-','Color',c_mat(2,:),'LineWidth',1.5)
end
% legend({})


%% fitting CASE 1
% initial guess
    %para0
% bound 
    lb= [1e-6 1e-6 1e-6 1e-6 1e-6];
    ub = para0*2;
% weight
 weight = ones(size(y1)); % uniform weighting
%weight = ones(size(y1)).*(exp(-x./para_hats(i_pulse,5)));


% fitting
        % options = optimset('display','iter','MaxIter',400,'MaxFunEvals',1e5,...
        % 'TolFun',1e-10,'TolX',1e-8,'FinDiffType','central');    
        % para_hat = fmincon(@(para)func_cost(y1,para,x,y2,weight),para0,[],[],[],[],lb,ub,[],options);
        options = optimset('display','iter', 'MaxIter',400, 'MaxFunEvals',1e5, ...
                       'TolFun',1e-10, 'TolX',1e-8, 'FinDiffType','central');
        ms = MultiStart('Display', 'iter', 'UseParallel', true); 
        problem = createOptimProblem('fmincon', ...
        'objective', @(para) func_cost(y1, para, x, y2, ones(size(y1))), ...
        'x0', para0, 'lb', lb, 'ub', ub, 'options', options);
        num_start_points = 1000;
        [para_hat, fval, exitflag, output, solutions] = run(ms, problem, num_start_points);
        para_hats(i_pulse, :) = para_hat;
    
% visualize
    y_model_hat = func_2RC(x,y2,para_hat);

%    yyaxis left
if ismember(i_pulse, selected_pulses)

subplot_index = find(selected_pulses == i_pulse); 
subplot(2,2,subplot_index); 
    plot(x,y_model_hat,'-','Color',c_mat(3,:),'LineWidth',1.5); 


legend({'Experimental Data', 'Initial Guess', 'Fitted Model'}, ...
    'Orientation', 'horizontal', 'FontSize', 6, 'Box', 'on');
grid on;
box on;


xlabel('Time (sec)');
ylabel('Voltage [V]');
end
end


% lgd = legend({'Experimental Data', 'Initial Guess', 'Fitted Model'}, ...
%     'Orientation', 'horizontal', 'FontSize', 10, 'Box', 'on');

% cd('G:\공유 드라이브\GSP_Data\driving_sample')
% save('2RC_para_2_scaled_ocv_vec_multi_1000_tot_last','para_hats')
% savefig('2RC_para_2_scaled_ocv_vec_multi_1000_tot_last')
% print('2RC_para_2_scaled_ocv_vec_multi_1000_tot_last','-dtiff','-r1200')
% %% initial guess
% 
% para0_mat =[para0; 0.001 0.001 65; 0.0013 0.0018*2 60; 0.0013 0.0018 60*1.5; 0.0013*0.6 0.0018*0.5 60];
% 
% 
% for i = 1: size(para0_mat,1)
% 
% 
%     para0 = para0_mat(i,:);
%     lb= [0 0 1];
%     ub = para0*10;
% 
%     % initial model 
%     y_model = func_1RC(x,y2,para0);
%     figure(4)
%     hold on
%     %yyaxis left
%     plot(x,y_model,'-','Color',c_mat(2,:));
% 
%     % fitted model
%     para_hat = fmincon(@(para)func_cost(y1,para,x,y2,weight),para0,[],[],[],[],lb,ub,[],options);
% 
%     % visualize
%     y_model_hat = func_1RC(x,y2,para_hat);
%     figure(4)
%     hold on
%     %    yyaxis left
%     plot(x,y_model_hat,'-','Color',c_mat(3,:))
% 
% end
% 
% 


% model
function y = func_2RC(t,I,para)
% x; time in sec
% para(1) = R0 [ohm]
% para (2) = R1 [ohm]
% para (3) = tau1 [sec]
% para (4) = R2 [ohm]
% para (5) = tau2 [sec]
% y = overpotential (V - OCV) [V]

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
cost = sqrt(mean((y_data - y_model).*weight).^2); % RMSE error

end

