
clear; clc; close all

% config
filename = 'ex_data_rOCV.csv';
filename_out = 'ex_data_rOCV.mat';
I_1C = 55.6; %[A]


%% Load

data = readtable(filename);


% Gather

intrim.t = data.TotTime_sec_; %[sec]
intrim.V = data.Voltage_V_;
intrim.I = data.Current_A_; %[A]
% intrim.T = data.

intrim.cycle = data.TotCycle;
intrim.step_p = data.StepNo;
intrim.type_p = data.Type;




% Step assign
    % charging/rest/discharging
intrim.type = char(zeros(size(intrim.t))); % assign space
intrim.type(intrim.I>0) = 'C';
intrim.type(intrim.I==0) = 'R';
intrim.type(intrim.I<0) = 'D';

    %step
n =1; % initial step
intrim.step(1) = n;

for i = 2:size(intrim.t,1)

    if intrim.type(i) == intrim.type(i-1)
        intrim.step(i) = n;
    else
        n=n+1;
        intrim.step(i) = n;
    end

end
    %check
% plot(intrim.step)

%% parsing
step_vec = unique(intrim.step);

for j = 1:length(step_vec)
    pdata(j).V = intrim.V(intrim.step==step_vec(j));
    pdata(j).t = intrim.t(intrim.step==step_vec(j));
    pdata(j).I = intrim.I(intrim.step==step_vec(j));
    pdata(j).Crate = pdata(j).I/I_1C;
    pdata(j).step = intrim.step(intrim.step==step_vec(j));
    pdata(j).type = intrim.type(intrim.step==step_vec(j));

    pdata(j).step = pdata(j).step(1);
    pdata(j).type = pdata(j).type(1);
end


% Save
save(filename_out,"pdata")