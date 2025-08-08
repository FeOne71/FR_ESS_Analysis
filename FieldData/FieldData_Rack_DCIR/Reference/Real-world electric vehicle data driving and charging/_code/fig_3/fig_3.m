%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright @ 2023 Stanford Energy Control Laboratory 
% (PI: Simona Onori, sonori@stanford.edu), 
% Stanford University. All Rights Reserved. 
% Developed by Luca Pulvirenti and Dr. Gabriele Pozzato 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

% Variables for plotting
Markersize = 18;
Fontsize = 20;
LineWidth = 5;
BarLineWidth = 1.5;

%% Load Variables
load('../../Drive/Folder2/Raw.mat')

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
I_calc  = [Pre;I;Post];
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
clear PeakTime PeakVoltage PeakCurrent R TimeR ;
z = 1;
for i=1:(length(I)-dt)
    if (I_filt(i+dt)-I_filt(i))>dI
        if Raw.Curr(i)>-thr && Raw.Curr(i)<thr
            if (I_filt(i+1)-I_filt(i))>ddI
                flag=1;
                for zi=1:dt
                    if filt_dI_dt(i+zi-1)<0 || Raw.Curr(i+zi)<0
                        flag=0;
                    end
                end
                if flag==1
                    if z==1
                        for j=1:dt
                            PeakTime{z}(j)=Raw.TimeCurr(i+j-1);
                            PeakCurrent{z}(j)=Raw.Curr(i+j-1);
                            PeakVoltage{z}(j)=Raw.Volt(i+j-1);
                        end
                        z=z+1;
                    else
                        if(Raw.TimeCurr(i)-PeakTime{z-1}(1))>dt
                            for j=1:dt
                                PeakTime{z}(j)=Raw.TimeCurr(i+j-1);
                                PeakCurrent{z}(j)=Raw.Curr(i+j-1);
                                PeakVoltage{z}(j)=Raw.Volt(i+j-1);
                            end
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

%% Figure 3 - Acceleration
% Auxiliary variables
t = Raw.TimeCurr;
V = Raw.Volt;
time_offset = 4335;

% Peaks on Voltage
color = [0 0 0]+0.03*15;
figure(1); hold on; box on
plot(t/60-time_offset,V,'k-','Linewidth',LineWidth-3,'Color',color,'HandleVisibility','off');
for i=1:length(PeakTime)
    if i == 1
        plot(PeakTime{i}/60-time_offset,PeakVoltage{i},'r-','Linewidth',LineWidth+3,'Color',[0 0 1]);
    else
        plot(PeakTime{i}/60-time_offset,PeakVoltage{i},'r-','Linewidth',LineWidth+3,'Color',[0 0 1],'HandleVisibility','off');
    end
end
xlabel('Time [min]','interpreter','tex');
ylabel('Voltage [V]','interpreter','tex');
xlim([0 45]); ylim([410 460])
set(gca,'fontsize',Fontsize);
set(gca,'ticklabelinterpreter','tex');
set(gca, 'XTick',(0:10:45));

% Peaks on Current
figure(2); hold on; box on
plot(t/60-time_offset,I,'k-','Linewidth',LineWidth-3,'Color',color)
for i=1:length(PeakTime)
    plot(PeakTime{i}/60-time_offset,PeakCurrent{i},'r-','Linewidth',LineWidth+3,'Color',[0 0 1]);
end
xlabel('Time [min]','interpreter','tex');
ylabel('Current [A]','interpreter','tex');
xlim([0 45]); ylim([-600 600])
set(gca,'fontsize',Fontsize);
set(gca,'ticklabelinterpreter','tex');
set(gca, 'XTick',(0:10:45));

% Peak "zoomed"
i = 16;
figure;
subplot(2,1,2); hold on; box on; 
plot(t/60-time_offset,I,'k-','Linewidth',LineWidth-3,'Color',color);
plot(PeakTime{i}/60-time_offset,PeakCurrent{i},'linewidth',7,'color','b');
ylabel('Current [A]','interpreter','tex');
xlabel('Time [min]','interpreter','tex');
set(gca,'fontsize',Fontsize);
xlim([7.14 7.185]); ylim([-5 120])
subplot(2,1,1); hold on; box on;
plot(t/60-time_offset,V,'k-','Linewidth',LineWidth-3,'Color',color);
plot(PeakTime{i}/60-time_offset,PeakVoltage{i},'linewidth',7,'color','b');
ylabel('Voltage [V]','interpreter','tex');
xlim([7.14 7.185]); ylim([446 450]);
set(gca,'ticklabelinterpreter','tex');
set(gca,'fontsize',Fontsize);

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Braking %%%%%%%%%%%%%%%%%%%%%%%%%%%%       
clear PeakTime PeakVoltage  PeakCurrent  ...
     R TimeR TimeEpoch

% Algorithm settings (dt is the number of samples in 1s, the sampling time is 0.01s)
dt = 100;
thr = 2;
dI = -100;

% Auxiliary variables
ddI = -1;
z = 1;

for i=1:(length(I_filt)-dt)
    if (I_filt(i+dt)-I_filt(i))<dI
        if Raw.Curr(i)>-thr && Raw.Curr(i)<thr
            if (I_filt(i+1)-I_filt(i))<ddI
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
                        end
                        z=z+1;
                    else
                        if(Raw.TimeCurr(i)-PeakTime{z-1}(1))>dt
                            for j=1:dt
                                PeakTime{z}(j)=Raw.TimeCurr(i+j-1);
                                PeakCurrent{z}(j)=Raw.Curr(i+j-1);
                                PeakVoltage{z}(j)=Raw.Volt(i+j-1);
                            end
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
    R(i) = DV/DI*1000;
end
    
%% Figure 3 - Braking
% Peaks on Voltage
figure(1); hold on
for i = 1:length(PeakTime)
    if i == 1
        plot(PeakTime{i}/60-time_offset,PeakVoltage{i},'r-','Linewidth',LineWidth+3,'Color',[1 0 0]);
    else
        plot(PeakTime{i}/60-time_offset,PeakVoltage{i},'r-','Linewidth',LineWidth+3,'Color',[1 0 0],'HandleVisibility','off');
    end
end
labels = {'Acceleration','Braking'};
legend(labels,'interpreter','tex');

% Peaks on Current
figure(2); hold on
for i=1:length(PeakTime)
    plot(PeakTime{i}/60-time_offset,PeakCurrent{i},'r-','Linewidth',LineWidth+3,'Color',[1 0 0]);
end

% Peak "zoomed"
i = 11;
figure;
subplot(2,1,2); hold on; box on; 
plot(t/60-time_offset,I,'k-','Linewidth',LineWidth-3,'Color',color);
plot(PeakTime{i}/60-time_offset,PeakCurrent{i},'linewidth',7,'color','r');
ylabel('Current [A]','interpreter','tex');
xlabel('Time [min]','interpreter','tex');
set(gca,'fontsize',Fontsize);
xlim([9.495 9.53]); ylim([-200 5])
subplot(2,1,1); hold on; box on;
plot(t/60-time_offset,V,'k-','Linewidth',LineWidth-3,'Color',color);
plot(PeakTime{i}/60-time_offset,PeakVoltage{i},'linewidth',7,'color','r');
ylabel('Voltage [V]','interpreter','tex');
xlim([9.495 9.53]); ylim([440 455]);
set(gca,'ticklabelinterpreter','tex');
set(gca,'fontsize',Fontsize);