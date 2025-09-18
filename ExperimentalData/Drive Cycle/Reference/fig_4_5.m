%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright @ 2023 Stanford Energy Control Laboratory 
% (PI: Simona Onori, sonori@stanford.edu), 
% Stanford University. All Rights Reserved. 
% Developed by Luca Pulvirenti and Dr. Gabriele Pozzato 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%% Resistance during braking and acceleration
% Folders
folders_driving = {'2' '4' '6' '8' '10' '12' '14' '16'};

for i_folder = 1:length(folders_driving)
    load(['D:\JCW\Projects\Reference\Real-world electric vehicle data driving and charging\Drive\Folder' folders_driving{i_folder} '/Raw.mat'])

    % Temperature, SOC, and Epoch must be interpolated to have the same length of Current
    Temp_int = interp1(Raw.TimeTemp,Raw.Temp,Raw.TimeCurr, 'linear','extrap');
    SoC_int = interp1(Raw.TimeSoC,Raw.SoC,Raw.TimeCurr, 'linear','extrap');
    SignalEpoch = interp1(Raw.TimeEpoch,Raw.Epoch,Raw.TimeCurr, 'linear','extrap');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Acceleration %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Algorithm settings (dt is the number of samples in 1s, the sampling time is 0.01s)
    dt = 100;                       
    thr = 5;
    dI = 100;

    % Auxiliary variable 
    ddI = 1;
    te = 100;

    % Filter on Current
    % Moving average of 1s --> 100 sampling points    
    N=length(Raw.Curr);
    I=Raw.Curr;
    I_filt = zeros(N,1);
    Pre = 0*ones(te/2,1);
    Post = 0*ones(te/2,1);
    I_calc = [Pre;I;Post];
    for i=1:N
        for m=0:te
            I_filt(i) = I_filt(i)+I_calc(i+m);
        end
        I_filt(i) = I_filt(i)/(te+1);
    end    
    
    % Derivative of Current
    for i=1:length(Raw.Curr)-1
        dI_dt(i) = (Raw.Curr(i+1)-Raw.Curr(i))/(Raw.TimeCurr(i+1)-Raw.TimeCurr(i));
    end
    
    % MA on Current Derivative
    N=length(dI_dt);
    filt_dI_dt = zeros(N,1);
    Pre = 0*ones(te/2,1);
    Post = 0*ones(te/2,1);
    calc_dI_dt  = [Pre;dI_dt';Post];
    for i=1:N
        for m=0:te
            filt_dI_dt(i) = filt_dI_dt(i)+calc_dI_dt(i+m);
        end
        filt_dI_dt(i) = filt_dI_dt(i)/(te+1);
    end
    
    % Driving Peaks Identification
    clear PeakTime PeakVoltage PeakSoC PeakCurrentMax PeakCurrent PeakT PeakEpoch ...
        R TimeR TimeEpoch PeakVoltageIn PeakCurrentIn PeakCurrentMax Capacity;
    z = 1;
    for i=1:(length(I)-dt)
        if (I_filt(i+dt)-I_filt(i))>dI
            if Raw.Curr(i)>-thr && Raw.Curr(i)<thr
                % The variariable ddI is introduced to ensure that a peak has an initial increase 
                % of at least 1A (peaks with initial flat plateaus are discarded)
                if (I_filt(i+1)-I_filt(i))>ddI
                    % The variable flag is introduced to check if the peak is always increasing and the current never changes sign 
                    flag=1;
                    for zi=1:dt
                        if filt_dI_dt(i+zi-1)<0 || Raw.Curr(i+zi)<0
                            flag=0;
                        end
                    end
                    if flag==1
                        if z==1
                            % The current, voltage, time, and epoch signals are saved for 100 time steps once a peak has been detected
                            for j=1:dt
                                PeakTime{z}(j)=Raw.TimeCurr(i+j-1);
                                PeakCurrent{z}(j)=Raw.Curr(i+j-1);
                                PeakVoltage{z}(j)=Raw.Volt(i+j-1);
                                PeakEpoch{z}(j)=SignalEpoch(i+j-1);
                            end
                            PeakSoC(z)=SoC_int(i);
                            PeakT(z)=Temp_int(i);
                            z=z+1;
                        else
                            % Additional condition that discards a new peak if a sufficient amount of time (dt=1s) is not passed from
                            % the beginning of the last peak (in this way the same peak cannot be saved twice in case of numerical errors) 
                            if(Raw.TimeCurr(i)-PeakTime{z-1}(1))>dt
                                for j=1:dt
                                    PeakTime{z}(j)=Raw.TimeCurr(i+j-1);
                                    PeakCurrent{z}(j)=Raw.Curr(i+j-1);
                                    PeakVoltage{z}(j)=Raw.Volt(i+j-1);
                                    PeakEpoch{z}(j)=SignalEpoch(i+j-1);
                                end
                                PeakSoC(z)=SoC_int(i);
                                PeakT(z)=Temp_int(i);
                                z=z+1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Resistance Computation
    for i = 1:length(PeakTime)
        DV = -(PeakVoltage{i}(end)-PeakVoltage{i}(1));
        DI = PeakCurrent{i}(end)-PeakCurrent{i}(1);
        R(i) = DV/DI*1000;
    end
    
    for i=1:length(PeakTime)
        TimeEpoch(i)=PeakEpoch{i}(1);
        TimeR(i)=PeakTime{i}(1);
        PeakVoltageIn(i)=PeakVoltage{i}(1);
        PeakCurrentIn(i)=PeakCurrent{i}(1);
        PeakCurrentMax(i)=PeakCurrent{i}(end);
    end
    DatesVector = datetime(TimeEpoch, 'convertfrom','posixtime');
    
    % Saving auxiliary .mat files
    PeaksACC.Time=TimeR;
    PeaksACC.Current=PeakCurrentIn;
    PeaksACC.Voltage=PeakVoltageIn;
    PeaksACC.Temperature=PeakT;
    PeaksACC.SoC=PeakSoC;
    PeaksACC.Epoch=TimeEpoch;
    PeaksACC.Date=DatesVector;
    PeaksACC.R=R;
    save(['mat/PeaksACC_' folders_driving{i_folder}],'PeaksACC');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Braking %%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    clear PeakTime PeakVoltage PeakSoC PeakCurrentMax PeakCurrent PeakT ...
        PeakEpoch R TimeR TimeEpoch PeakVoltageIn PeakCurrentIn PeakCurrentMax Capacity
    
    % Algorithm settings (dt is the number of samples in 1s, the sampling time is 0.01s)
    dt = 100;
    thr = 2;
    dI= -100;

    % Auxiliary variables
    ddI = -1;
    z = 1;
    
    for i=1:(length(I_filt)-dt)
        if (I_filt(i+dt)-I_filt(i))<dI
            if Raw.Curr(i)>-thr && Raw.Curr(i)<thr
                % The variariable ddI is introduced to ensure that a peak has an initial decrease 
                % of at least -1A (peaks with initial flat plateaus are discarded)                
                if (I_filt(i+1)-I_filt(i))<ddI
                    % The variable flag is introduced to check if the peak is always decreasing and the current never changes sign 
                    flag=1;
                    for zi=1:dt
                        if filt_dI_dt(i+zi-1)>0 || Raw.Curr(i+zi)>0
                            flag=0;
                        end
                    end
                    if flag==1
                        if z==1
                            for j=1:dt
                                PeakTime{z}(j)=Raw.TimeCurr(i+j-1);
                                PeakCurrent{z}(j)=Raw.Curr(i+j-1);
                                PeakVoltage{z}(j)=Raw.Volt(i+j-1);
                                PeakEpoch{z}(j)=SignalEpoch(i+j-1);
                            end
                            PeakSoC(z)=SoC_int(i);
                            PeakT(z)=Temp_int(i);
                            z=z+1;
                        else
                            if(Raw.TimeCurr(i)-PeakTime{z-1}(1))>dt
                                for j=1:dt
                                    PeakTime{z}(j)=Raw.TimeCurr(i+j-1);
                                    PeakCurrent{z}(j)=Raw.Curr(i+j-1);
                                    PeakVoltage{z}(j)=Raw.Volt(i+j-1);
                                    PeakEpoch{z}(j)=SignalEpoch(i+j-1);
                                end
                                PeakSoC(z)=SoC_int(i);
                                PeakT(z)=Temp_int(i);
                                z=z+1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Resistance Computation
    for i=1:length(PeakTime)
        DV = -(PeakVoltage{i}(end)-PeakVoltage{i}(1));
        DI = PeakCurrent{i}(end)-PeakCurrent{i}(1);
        R(i)=DV/DI*1000;
    end
    
    for i=1:length(PeakTime)
        TimeEpoch(i)=PeakEpoch{i}(1);
        TimeR(i)=PeakTime{i}(1);
        PeakVoltageIn(i)=PeakVoltage{i}(1);
        PeakCurrentIn(i)=PeakCurrent{i}(1);
        PeakCurrentMax(i)=PeakCurrent{i}(end);
    end
    DatesVector=datetime(TimeEpoch, 'convertfrom','posixtime');
    
    % Saving auxiliary .mat files
    PeaksBR.Time=TimeR;
    PeaksBR.Current=PeakCurrentIn;
    PeaksBR.Voltage=PeakVoltageIn;
    PeaksBR.Temperature=PeakT;
    PeaksBR.SoC=PeakSoC;
    PeaksBR.Epoch=TimeEpoch;
    PeaksBR.Date=DatesVector;
    PeaksBR.R=R;
    save(['mat/PeaksBR_' folders_driving{i_folder}],'PeaksBR');

    clearvars -except i_folder folders_driving 
end

%% Create one mat file
% Braking
PeaksBR_all.Time = [];
PeaksBR_all.Current = [];
PeaksBR_all.Voltage = [];
PeaksBR_all.Temperature = [];
PeaksBR_all.SoC = [];
PeaksBR_all.Epoch = [];
PeaksBR_all.Date = [];
PeaksBR_all.R = [];
PeaksBR_all.Folder = [];

% Acceleration
PeaksACC_all.Time = [];
PeaksACC_all.Current = [];
PeaksACC_all.Voltage = [];
PeaksACC_all.Temperature = [];
PeaksACC_all.SoC = [];
PeaksACC_all.Epoch = [];
PeaksACC_all.Date = [];
PeaksACC_all.R = [];
PeaksACC_all.Folder = [];

for i_folder = 1:length(folders_driving)
    % Load partial .mat
    load(['mat/PeaksBR_' folders_driving{i_folder} '.mat'])
    load(['mat/PeaksACC_' folders_driving{i_folder} '.mat'])

    % Braking
    PeaksBR_all.Time = [PeaksBR_all.Time PeaksBR.Time];
    PeaksBR_all.Current = [PeaksBR_all.Current PeaksBR.Current];
    PeaksBR_all.Voltage = [PeaksBR_all.Voltage PeaksBR.Voltage];
    PeaksBR_all.Temperature = [PeaksBR_all.Temperature PeaksBR.Temperature];
    PeaksBR_all.SoC = [PeaksBR_all.SoC PeaksBR.SoC];
    PeaksBR_all.Epoch = [PeaksBR_all.Epoch PeaksBR.Epoch];
    PeaksBR_all.Date = [ PeaksBR_all.Date PeaksBR.Date];
    PeaksBR_all.R = [PeaksBR_all.R PeaksBR.R];
    PeaksBR_all.Folder = [PeaksBR_all.Folder str2num(folders_driving{i_folder})*ones(1,length(PeaksBR.Time))];

    % Acceleration
    PeaksACC_all.Time = [PeaksACC_all.Time PeaksACC.Time];
    PeaksACC_all.Current = [PeaksACC_all.Current PeaksACC.Current];
    PeaksACC_all.Voltage = [PeaksACC_all.Voltage PeaksACC.Voltage];
    PeaksACC_all.Temperature = [PeaksACC_all.Temperature PeaksACC.Temperature];
    PeaksACC_all.SoC = [PeaksACC_all.SoC PeaksACC.SoC];
    PeaksACC_all.Epoch = [PeaksACC_all.Epoch PeaksACC.Epoch];
    PeaksACC_all.Date = [ PeaksACC_all.Date PeaksACC.Date];
    PeaksACC_all.R = [PeaksACC_all.R PeaksACC.R];
    PeaksACC_all.Folder = [PeaksACC_all.Folder str2num(folders_driving{i_folder})*ones(1,length(PeaksACC.Time))];

    clearvars -except folders_driving i_folder PeaksBR_all PeaksACC_all
end
save('mat/Peaks_all.mat','PeaksACC_all','PeaksBR_all')

%% Load all data and remove outliers
% Braking
pd_BR = fitdist(PeaksBR_all.R', 'Normal');
out_max_BR = pd_BR.mu + 3*pd_BR.sigma;
out_min_BR = pd_BR.mu - 3*pd_BR.sigma;
r = 1;
while r <= length(PeaksBR_all.R) % Removes outliers
    if(PeaksBR_all.R(r) > out_max_BR || PeaksBR_all.R(r) < out_min_BR) 
        f = fieldnames(PeaksBR_all);
        for j = 1:length(f)
            PeaksBR_all.(f{j})(r) = [];
            
        end
    else
        r = r + 1;
    end
end
pd_BR = fitdist(PeaksBR_all.R.', 'Normal');

% Acceleration
pd_ACC = fitdist(PeaksACC_all.R', 'Normal');
out_max_ACC = pd_ACC.mu + 3*pd_ACC.sigma;
out_min_ACC = pd_ACC.mu - 3*pd_ACC.sigma;
r = 1;
while r <= length(PeaksACC_all.R) % Removes outliers
    if(PeaksACC_all.R(r) > out_max_ACC || PeaksACC_all.R(r) < out_min_ACC) 
        f = fieldnames(PeaksACC_all);
        for j = 1:length(f)
            PeaksACC_all.(f{j})(r) = [];
            
        end
    else
        r = r + 1;
    end
end
pd_ACC = fitdist(PeaksACC_all.R.', 'Normal');
save('mat/Peaks_all','PeaksACC_all','PeaksBR_all');

%% Load temperature data from https://weatherspark.com/y/545/Average-Weather-in-Palo-Alto-California-United-States-Year-Round
load('mat/temperature_data.mat')
if i == 1
    low.month = round(Data001(1:41,1),1);
    low.temp = Data001(1:41,2);
    high.month = round(Data001(42:end,1),1);
    high.temp = Data001(42:end,2);
else
    low.month = round(Data002(1:47,1),1);
    low.temp = Data002(1:47,2);
    high.month = round(Data002(48:end,1),1);
    high.temp = Data002(48:end,2);    
end

% Avg
interp_lowT = interp1(my_quasi_unique([low.month; 12+low.month]),[fa_to_ce(low.temp); fa_to_ce(low.temp)],1:0.01:25);
interp_highT = interp1(my_quasi_unique([high.month; 12+high.month]),[fa_to_ce(high.temp); fa_to_ce(high.temp)],1:0.01:25);

interp_avgT = interp1(1:0.01:25,mean([interp_lowT' interp_highT'],2),1.5:1:23.5);
for i = 1:length(1.5:1:23.5)
    if i <= 12
        x_date(i) = datetime(2019,i,1,00,0,0);
    else
        x_date(i) = datetime(2020,i-12,1,0,0,0);
    end
end

%% Figure 4
% a
figure;
h = histfit(PeaksBR_all.R);
h(1).EdgeColor = [0.3010 0.7450 0.9330];
h(2).Color = [0.9290 0.6940 0.1250];
h(2).LineWidth = 4;
xlabel('R_{BR} [m\Omega]'); ylabel('Frequency [-]'); xlim([20 45]); 
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','interpreter'),'interpreter','tex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')

% b
figure; box on; 
h = histfit(PeaksACC_all.R);
h(1).EdgeColor = [0.3010 0.7450 0.9330];
h(2).Color = [0.9290 0.6940 0.1250];
h(2).LineWidth = 4;
xlabel('R_{ACC} [m\Omega]'); ylabel('Frequency [-]'); xlim([20 45]);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','interpreter'),'interpreter','tex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')

% c
figure; subplot(2,1,1); box on; hold on; yyaxis left; 
for i=1:length(folders_driving)
    load(['mat/PeaksBR_' folders_driving{i} '.mat'])
    scatter(PeaksBR.Date,PeaksBR.R,100,PeaksBR.Temperature,'filled','LineWidth',2); 
    hold on;
end
ylabel('R_{BR} [m\Omega]'); ylim([10 50]); xticks([])
set(gca,'YColor','k'); set(gca, 'YTick',(10:10:50)); 
xticks([datetime(2019,11,1) datetime(2019,12,1) datetime(2020,1,1) datetime(2020,2,1) datetime(2020,3,1) datetime(2020,4,1) datetime(2020,5,1) ...
    datetime(2020,6,1) datetime(2020,7,1) datetime(2020,8,1) datetime(2020,9,1) datetime(2020,10,1) datetime(2020,11,1)]); xticklabels([])
yyaxis right; 
plot(x_date,interp_avgT,'-o','color','b','linewidth',2,'MarkerSize',8,'MarkerFaceColor','b')
ylabel('Monthly temperature [°C]','Position',[383,8.5,-1]); xlim([x_date(11) x_date(23)]);
set(gca,'YColor','b'); colormap(flipud(autumn)); 
c = colorbar('eastoutside','Position',[0.9,0.165,0.016,0.75],'Limits',[10 40]); 
c.Label.String = 'Temperature [°C]';
set(gca,'InnerPosition',[0.13,0.60,0.68,0.32])
subplot(2,1,2); box on; hold on; yyaxis left
for i=1:length(folders_driving)
    load(['mat/PeaksACC_' folders_driving{i} '.mat'])
    scatter(PeaksACC.Date,PeaksACC.R,100,PeaksACC.Temperature,'filled','LineWidth',2); 
    hold on;
end
ylabel('R_{ACC} [m\Omega]'); ylim([10 50]); xtickangle(45);
set(gca, 'YTick',(10:10:50)); set(gca,'YColor','k');
set(gca,'InnerPosition',[0.13,0.15,0.68,0.32])
xticks([datetime(2019,11,1) datetime(2019,12,1) datetime(2020,1,1) datetime(2020,2,1) datetime(2020,3,1) datetime(2020,4,1) datetime(2020,5,1) ...
    datetime(2020,6,1) datetime(2020,7,1) datetime(2020,8,1) datetime(2020,9,1) datetime(2020,10,1) datetime(2020,11,1)])
xticklabels({'Nov 2019','Dec 2019','Jan 2020','Feb 2020','Mar 2020','Apr 2020','May 2020','Jun 2020','Jul 2020','Aug 2020','Sep 2020','Oct 2020'});
yyaxis right; 
plot(x_date,interp_avgT,'-o','color','b','linewidth',2,'MarkerSize',8,'MarkerFaceColor','b'); xlim([x_date(11) x_date(23)]);
set(gca,'YColor','b'); colormap(flipud(autumn)); 
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','interpreter'),'interpreter','tex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')
set(gcf,'Position',[100 100 1000 600])

%% Figure 5
months = 1:1:12;
temp_step = 1;

% Braking
temp_low = floor(min(PeaksBR_all.Temperature));
temp_high = ceil(max((PeaksBR_all.Temperature)));
temps = temp_low:temp_step:temp_high;
data_BR = zeros(length(months), length(temps));
entries_BR = zeros(length(months), length(temps));
for i = 1:length(PeaksBR_all.R)
    m = month(PeaksBR_all.Date(i));
    if (m > 10) % shift months index
        M_index = mod(m, 10);
    else
        M_index = m + 2;
    end
    T_index = floor(PeaksBR_all.Temperature(i)) - temp_low + 1;
    data_BR(M_index, T_index) = data_BR(M_index, T_index) + PeaksBR_all.R(i);
    entries_BR(M_index, T_index) = entries_BR(M_index, T_index) + 1;
end
AvgR_BR = data_BR./entries_BR;


fig = figure; grid on; hold on
my_map = colormap("winter");
for i = 1:size(AvgR_BR,1)
    stem3((temps),months(i)*ones(1,length(temps)),(AvgR_BR(i,:)),':d','markersize',9,'color',[1 1 1]*0.3,'linewidth',1,'markerfacecolor',my_map((i-1+1)*21,:))
end
xlabel('Temperature [°C]'); zlabel('R_{BR} [m\Omega]')
yticklabels({'','Nov 2019', '', 'Jan 2020', '', 'Mar 2020', '', 'May 2020', '', 'Jul 2020', '', 'Sep 2020', '',''}); yticks([0:13]);
xlim([10 35]); ylim([0 13]); zlim([20 45]); view(52,51)
set(findall(gcf,'-property','FontSize'),'FontSize',20); 
set(findall(gcf,'-property','interpreter'),'interpreter','tex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')

fig = figure; box on; hold on
for i = 1:size(AvgR_BR,1)
    plot(temps,AvgR_BR(i,:),'d','markersize',9,'color',[1 1 1]*0.3,'linewidth',1.5,'markerfacecolor',my_map((i-1+1)*21,:))
end
xlim([10 35]); ylim([20 45])
xlabel('Temperature [°C]'); ylabel('R_{BR} [m\Omega]')
set(findall(gcf,'-property','FontSize'),'FontSize',20); 
set(findall(gcf,'-property','interpreter'),'interpreter','tex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')

% Acceleration
temp_low = floor(min(PeaksACC_all.Temperature));
temp_high = ceil(max((PeaksACC_all.Temperature)));
temps = temp_low:temp_step:temp_high;
data_ACC = zeros(length(months), length(temps));
entries_ACC = zeros(length(months), length(temps));
for i = 1:length(PeaksACC_all.R)
    m = month(PeaksACC_all.Date(i));
    if (m > 10) % shift months index
        M_index = mod(m, 10);
    else
        M_index = m + 2;
    end
    T_index = floor(PeaksACC_all.Temperature(i)) - temp_low + 1;
    data_ACC(M_index, T_index) = data_ACC(M_index, T_index) + PeaksACC_all.R(i);
    entries_ACC(M_index, T_index) = entries_ACC(M_index, T_index) + 1;
end
AvgR_ACC = data_ACC./entries_ACC;

fig = figure; grid on; hold on
my_map = colormap("autumn");
for i = 1:size(AvgR_ACC,1)
    stem3((temps),months(i)*ones(1,length(temps)),(AvgR_ACC(i,:)),':d','markersize',9,'color',[1 1 1]*0.3,'linewidth',1,'markerfacecolor',my_map((i-1+1)*21,:))
end
xlabel('Temperature [°C]'); zlabel('R_{ACC} [m\Omega]'); colormap(winter)
yticklabels({'','Nov 2019', '', 'Jan 2020', '', 'Mar 2020', '', 'May 2020', '', 'Jul 2020', '', 'Sep 2020', '',''}); yticks([0:13]);
xlim([10 35]); ylim([0 13]); zlim([20 45]); view(52,51)
set(findall(gcf,'-property','FontSize'),'FontSize',20); 
set(findall(gcf,'-property','interpreter'),'interpreter','tex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')

fig = figure; box on; hold on
for i = 1:size(AvgR_ACC,1)
    plot(temps,AvgR_ACC(i,:),'d','markersize',9,'color',[1 1 1]*0.3,'linewidth',1.5,'markerfacecolor',my_map((i-1+1)*21,:))
end
xlim([10 35]); ylim([20 45])
xlabel('Temperature [°C]'); ylabel('R_{ACC} [m\Omega]')
set(findall(gcf,'-property','FontSize'),'FontSize',20); 
set(findall(gcf,'-property','interpreter'),'interpreter','tex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')

%% Functions
% Fahrenheit to Celsius
function deg_ce = fa_to_ce(deg_f)
    deg_ce = (5/9)*(deg_f-32);
end