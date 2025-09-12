

% config
filename_in = 'ex_data_RPT.mat';
filename_out = 'pOCV.mat';
Crate_detect = 0.05;
delta = 0.05;
min_duration = 10;%[hr]
num_point = 201;


% load pdata
load(filename_in) % pdata;


% identify OCV step


for i = 1:size(pdata,2)
    pdata(i).Crate_avg = mean(pdata(i).Crate);
    pdata(i).step_duration = pdata(i).t(end)-pdata(i).t(1);
end

detect = abs([pdata.Crate_avg]') < Crate_detect*(1+delta)...
    & abs([pdata.Crate_avg]') > Crate_detect*(1-delta)...
    & [pdata.step_duration]' > min_duration*3600;

ocv_step = find(detect);

rdata = pdata(detect); %reduced data


%% SOC define

for j = 1:size(rdata,2)
    
    if rdata(j).Crate_avg > 0 %charging
    rdata(j).SOC = cumtrapz(rdata(j).t,rdata(j).I)/trapz(rdata(j).t,rdata(j).I);
    else % discharging
    rdata(j).SOC = 1-cumtrapz(rdata(j).t,rdata(j).I)/trapz(rdata(j).t,rdata(j).I);
    end

end




%% Interpolation (pseudo ocv )
SOC_vec = linspace(0,1,num_point);
OCVc_vec =  interp1(rdata(1).SOC,rdata(1).V,SOC_vec); %either look-up or interp1
OCVd_vec =  interp1(rdata(2).SOC,rdata(2).V,SOC_vec); %either look-up or interp1


plot(SOC_vec,OCVc_vec,'o'); hold on
plot(SOC_vec,OCVd_vec,'o')


%% output
OCV.OCVc = OCVc_vec;
OCV.OCVd = OCVd_vec;
OCV.SOC = SOC_vec;

save(filename_out,'OCV')