
clear;clc;close all

% Config
filename = "ex_data_RPT.csv";
filename_out = 'ex_data_RPT.mat';

I_1C = 2.0/1000; %[A]



%% load

data = readtable(filename);


% Gather

intrim.t = seconds(data.TotalTime); %[sec]
intrim.V = data.Voltage_V_;
intrim.I = data.Current_mA_/1000; %[A]
% intrim.T = data.

intrim.cycle = data.CycleIndex;
intrim.step_p = data.StepIndex;
intrim.type_p = data.StepType;


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