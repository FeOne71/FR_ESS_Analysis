%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright @ 2023 Stanford Energy Control Laboratory 
% (PI: Simona Onori, sonori@stanford.edu), 
% Stanford University. All Rights Reserved. 
% Developed by Luca Pulvirenti and Dr. Gabriele Pozzato 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%% Charging impedance 
folders_ch = {'1' '3' '5' '7' '9' '11' '13' '15'};
 
for i_folder = 1:length(folders_ch)
    load(['../../Charge/Folder' folders_ch{i_folder} '/Raw.mat']) 
    
    % Signals pre-processing/interpolation 
    Temp_int = interp1(Raw.TimeTemp,Raw.Temp,Raw.TimeCurr, 'linear','extrap');
    SoC_int = interp1(Raw.TimeSoC,Raw.SoC,Raw.TimeCurr, 'linear','extrap');
    SignalEpoch = interp1(Raw.TimeEpoch,Raw.Epoch,Raw.TimeCurr, 'linear','extrap');

    j=1; k=1;
    for i=1:length(Raw.TimeCurr)-1
    % A profile is considered as a single charging event if it is separated 
    % from previous and next charging events by at least 2 minutes    
        if (Raw.TimeCurr(i+1)-Raw.TimeCurr(i)) < 120
            Time_n_all{i_folder}{j}(k) = Raw.TimeCurr(i);
            Current_n_all{i_folder}{j}(k) = Raw.Curr(i);
            Voltage_n_all{i_folder}{j}(k) = Raw.Volt(i);
            SoC_n_all{i_folder}{j}(k) = SoC_int(i);
            Temp_n_all{i_folder}{j}(k) = Temp_int(i);
            % We rebuild the epoch signal since it is not well interpolated
            if k==1
                Epoch_n_all{i_folder}{j}(k) = SignalEpoch(i);
            else
                Epoch_n_all{i_folder}{j}(k) = Epoch_n_all{i_folder}{j}(k-1)+0.01;
            end
            k = k+1;
        else
            j = j+1;
            k = 1;
        end
    end
    for i=1:length(Epoch_n_all{i_folder})
        DatesVector_all{i_folder}{i} = datetime(Epoch_n_all{i_folder}{i}, 'convertfrom','posixtime');
    end

    clear Raw.TimeCurr Raw.Curr Raw.TimeVolt Raw.Volt Raw.TimeSoC Raw.SoC Raw.TimeTemp Raw.Temp Raw.TimeEpoch Raw.Epoch
    
end

%% Convolute all the charging events in one structure
k=1;
for i_folder=1:length(folders_ch)
    for i=1:length(Time_n_all{i_folder})
        Time_n{k}=Time_n_all{i_folder}{i};
        Current_n{k}=Current_n_all{i_folder}{i};
        Voltage_n{k}=Voltage_n_all{i_folder}{i};
        SoC_n{k}=SoC_n_all{i_folder}{i};
        Temp_n{k}=Temp_n_all{i_folder}{i};        
        Epoch_n{k}=Epoch_n_all{i_folder}{i};
        DatesVector{k}=DatesVector_all{i_folder}{i};        
        k=k+1;
    end
end
    
%% Modify the data       
% Discard transient phases
zi=1;
initial = 18000;
final = 18000;

interval = (initial+final)/(100);

for i=1:length(Time_n)
    if Time_n{i}(end)-Time_n{i}(1) > 2*interval %Time sufficient to delete the initial and final sections
        if mean(Current_n{i}) < 0
            Time_mod{zi} = Time_n{i}(initial:end-final);
            Current_mod{zi} = Current_n{i}(initial:end-final);
            Voltage_mod{zi} = Voltage_n{i}(initial:end-final);
            SoC_mod{zi} = SoC_n{i}(initial:end-final);
            Temp_mod{zi} = Temp_n{i}(initial:end-final);
            DatesVector_mod{zi} = DatesVector{i}(initial:end-final);
            Time{zi} = (Time_mod{zi}-Time_mod{zi}(1))/3600;

            clear x_mark;
            clear y_mark;

            % Linear Fitting of voltage signal
            k=1;
            for j=1:length(Time_mod{zi})
                if j==1
                    x_mark(k)=Time_mod{zi}(j);
                    y_mark(k)=Voltage_mod{zi}(j);
                    k=k+1;
                else
                    flag=1;
                    for xi=1:k-1
                        if Voltage_mod{zi}(j)<=y_mark(xi)
                            flag=0;
                        end
                    end

                    if flag==1
                        x_mark(k)=Time_mod{zi}(j);
                        y_mark(k)=Voltage_mod{zi}(j);
                        k=k+1;
                    end
                end
            end
            yfit{zi}= interp1(x_mark,y_mark,Time_mod{zi}, 'linear','extrap');
            zi=zi+1;
        end
    end
end


%% Discard not relevant charging events
% One fast-charging event (18) at 1.5C is present and not considered in the analysis 
% because of its low statistical significance. Additionally, four charging events 
% showing discontinuous current patterns are discarded from the analysis (8,20,36,53)
k=1; z=1;
for i=1:length(Time_mod)
    if i == 8 || i == 18 || i == 20 || i == 36 || i == 53
        Time_mod_disc{z} = Time_mod{i};
        Current_disc{z} = Current_mod{i};
        Volt_disc{z} = yfit{i};
        SoC_disc{z} = SoC_mod{i};
        Temp_disc{z} = Temp_mod{i};
        DatesVector_disc{z} = DatesVector_mod{i};
        Time_disc{z} = Time{i};
        z = z+1;         
    else
        Time_mod_clean{k} = Time_mod{i};
        Current_clean{k} = Current_mod{i};
        Volt_clean{k} = yfit{i};
        SoC_clean{k} = SoC_mod{i};
        Temp_clean{k} = Temp_mod{i};
        DatesVector_clean{k} = DatesVector_mod{i};
        Time_clean{k} = Time{i};
        k = k+1;
    end
end

%% Compute charge impedance
dt = [1 10000];                                                             % Time window: [0.1s 100s]                   
for k_dt = 1:length(dt)
    for i=1:length(Time_mod_clean)
        
        if length(Time_mod_clean{i}) > dt                      
            % The charging impedance curve is computed by dividing the voltage drop over a time window Δt in the future
            % The last charging impedance value is computed a time window Δt before the end of the vector 
            for j=1:(length(Current_clean{i})-dt(k_dt))
                R_pseudo{k_dt}{i}(j) = abs((Volt_clean{i}(j+dt(k_dt))-Volt_clean{i}(j)))/Current_clean{i}(j)*(-1000);
            end
        end
    end
end

%% Average current, charging levels, and battery pack capacity
% Average current
for i=1:length(Time_mod_clean)
    I_ave(i) = mean(Current_clean{i});
    DatesVector_first(i) = DatesVector_clean{i}(1);
end

% Battery pack capacity [Ah]
Cap = 240;

% Sampling time [s]
Ts = 0.01;

%% Division depending on charging levels
k=1; h=1; j=1;
for i=1:length(I_ave)
    if I_ave(i)>-5
        Time_mod_1{k} = Time_mod_clean{i};
        Curr_mod_1{k} = Current_clean{i};
        SoC_mod_1{k} = SoC_clean{i};
        Temp_mod_1{k} = Temp_clean{i};
        DatesVector_mod_1{k} = DatesVector_clean{i};
        R_pseudo_1{1}{k} = R_pseudo{1}{i};
        R_pseudo_1{2}{k} = R_pseudo{2}{i};
        I_ave_1(k)=I_ave(i);
        lev1(k)=i;
        k=k+1;
    elseif I_ave(i)>-80 && I_ave(i)<=-5
        Time_mod_2{h} = Time_clean{i};
        Curr_mod_2{h} = Current_clean{i};
        SoC_mod_2{h} = SoC_clean{i};
        Temp_mod_2{h} = Temp_clean{i};
        DatesVector_mod_2{h} = DatesVector_clean{i};
        R_pseudo_2{1}{h} = R_pseudo{1}{i};
        R_pseudo_2{2}{h} = R_pseudo{2}{i}; 
        I_ave_2(h)=I_ave(i);
        lev2(h)=i;
        h=h+1;
    elseif I_ave(i)<=-80
        Time_mod_3{j} = Time_clean{i};
        Curr_mod_3{j} = Current_clean{i};
        SoC_mod_3{j} = SoC_clean{i};
        Temp_mod_3{j} = Temp_clean{i};
        DatesVector_mod_3{j} = DatesVector_clean{i};
        R_pseudo_3{1}{j} = R_pseudo{1}{i};
        R_pseudo_3{2}{j} = R_pseudo{2}{i}; 
        I_ave_3(j)=I_ave(i);
        lev3(j)=i;
        j=j+1;
    end
end

%% Figure 6a
figure; 
sz = 300; Line = 0.5;
coloredge = 'k';
for i = 1:length(I_ave_1)
    if i == 1 
        scatter(DatesVector_mod_1{i}(1),-I_ave_1(i)/Cap,sz,'MarkerEdgeColor',coloredge,...
                      'MarkerFaceColor',[0.929 0.6940 0.1250],...
                      'LineWidth',Line,'DisplayName','C/240')
                  hold on;
    else
        scatter(DatesVector_mod_1{i}(1),-I_ave_1(i)/Cap,sz,'MarkerEdgeColor',coloredge,...
              'MarkerFaceColor',[0.9290 0.6940 0.1250],...
              'LineWidth',Line,'HandleVisibility','off')
    end
end
for i = 1:length(I_ave_2)
    if i == 1
        scatter(DatesVector_mod_2{i}(1),-I_ave_2(i)/Cap,sz,'MarkerEdgeColor',coloredge,...
                      'MarkerFaceColor',[0 0.4470 0.7410],...
                      'LineWidth',Line,'DisplayName','C/20')
    else
        scatter(DatesVector_mod_2{i}(1),-I_ave_2(i)/Cap,sz,'MarkerEdgeColor',coloredge,...
                      'MarkerFaceColor',[0 0.4470 0.7410],...
                      'LineWidth',Line,'HandleVisibility','off')
    end
end
hold on;
for i = 1:length(I_ave_3)
    if i == 1
        scatter(DatesVector_mod_3{i}(1),-I_ave_3(i)/Cap,sz,'MarkerEdgeColor',coloredge,...
                      'MarkerFaceColor',[0.6350 0.0780 0.1840],...
                      'LineWidth',Line,'DisplayName','C/2')
    else
        scatter(DatesVector_mod_3{i}(1),-I_ave_3(i)/Cap,sz,'MarkerEdgeColor',coloredge,...
                  'MarkerFaceColor',[0.6350 0.0780 0.1840],...
                  'LineWidth',Line,'HandleVisibility','off')
    end
end
ylabel('C-rate [1/h]','interpreter','tex'); 
yticks([0 1/10 1/5 1/3.3 1/2.5 1/2]); yticklabels({'' 'C/10' 'C/5' 'C/3.3' 'C/2.5' 'C/2'})
xticks([datetime(2019,11,1) datetime(2020,1,1) datetime(2020,3,1) datetime(2020,5,1) datetime(2020,7,1) datetime(2020,9,1)])
ylabel('C-rate [1/h]','interpreter','tex'); box on
ylim([0 0.5]); xtickangle(45); set(gca,'fontsize',20);

%% Figure 6b
x1 = 1;
Bin1 = length(lev1);
x2 = 2;
Bin2 = length(lev2);
x3 = 3;
Bin3 = length(lev3);
Bin_axis = [1;2;3];

bars = [length(lev1) length(lev2) length(lev3)];
% Histograms 
figure;
bar(x1,Bin1,'FaceColor','[0.9290 0.6940 0.1250]','EdgeColor',[0 0 0],'LineWidth',0.5);
hold all;
bar(x2,Bin2,'FaceColor','[0 0.4470 0.7410]','EdgeColor',[0 0 0],'LineWidth',0.5);
bar(x3,Bin3,'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0 0 0],'LineWidth',0.5);
ylabel('Number of charging events [#]','interpreter','tex');
xticklabels({'C/240','C/20','C/2'}); xticks(1:3); xlim([0 4]);
set(gca,'ticklabelinterpreter','tex');
set(gca,'fontsize',20);

%% Figure 11 - C/240
% The step (tuned by the user) is used to reduce the plotting time and figure size
step = 100; 

% Plot
remove_transient = 2000/Ts; sz = 20;
for k = 1:length(dt)
    f = figure; box on; hold all
    for i=1:length(R_pseudo_1{1,k})
        scatter(SoC_mod_1{i}(remove_transient:step:end-dt(k)),R_pseudo_1{1,k}{1,i}(remove_transient:step:end),...
            sz,Temp_mod_1{i}(remove_transient:step:end-dt(k)),'filled','LineWidth',2); 
    end
    xlabel('SoC [%]','interpreter','tex');
    ylabel('Z_{CHG} [m\Omega]','interpreter','tex');
    xlim([0 100]); ylim([0 0.40*dt(k)/100])
    colormap(flipud(autumn));
    c=colorbar('eastoutside'); caxis([10 32.5]);
    c.Label.String = 'Temperature [°C]';
    c.Label.Interpreter = 'tex';
    set(findall(gcf,'-property','interpreter'),'interpreter','tex')
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
end

%% Figure 6c, 11, 12 - C/20
% Figure 6c
remove_transient = 3600/Ts; 
for ce_indx = 1:size(R_pseudo_2{2},2)
    R_pseudo_2_filt{ce_indx} = movmean(R_pseudo_2{2}{ce_indx}(1:end),50000);
end

figure; box on; hold all
for i=1:length(R_pseudo_2{1,k})
    scatter(SoC_mod_2{i}(remove_transient:step:end-dt(2)),R_pseudo_2_filt{1,i}(remove_transient:step:end),...
        sz,Temp_mod_2{i}(remove_transient:step:end-dt(2)),'filled','LineWidth',2); 
end
xlabel('SoC [%]','interpreter','tex');
ylabel('Z_{CHG} [m\Omega]','interpreter','tex');
xlim([0 100]); ylim([0 0.25*dt(2)/100])
colormap(flipud(autumn));
c=colorbar('eastoutside'); caxis([10 32.5]);
c.Label.String = 'Temperature [°C]';
c.Label.Interpreter = 'tex';
set(findall(gcf,'-property','interpreter'),'interpreter','tex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')
set(findall(gcf,'-property','FontSize'),'FontSize',20);

% Figure 11
remove_transient = 2000/Ts; 
for k = 1:length(dt)
    f=figure; box on; hold all
    for i=1:length(R_pseudo_2{1,k})
        scatter(SoC_mod_2{i}(remove_transient:step:end-dt(k)),R_pseudo_2{1,k}{1,i}(remove_transient:step:end),...
            sz,Temp_mod_2{i}(remove_transient:step:end-dt(k)),'filled','LineWidth',2); 
    end
    xlabel('SoC [%]','interpreter','tex');
    ylabel('Z_{CHG} [m\Omega]','interpreter','tex');
    xlim([0 100]); ylim([0 0.25*dt(k)/100])
    colormap(flipud(autumn));
    c=colorbar('eastoutside'); caxis([10 32.5]);
    c.Label.String = 'Temperature [°C]';
    c.Label.Interpreter = 'tex';
    set(findall(gcf,'-property','interpreter'),'interpreter','tex')
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
end

% Figure 12
remove_transient = 3600/Ts; 
figure; box on; hold all
for i=1:length(R_pseudo_2{1,k})
    scatter(SoC_mod_2{i}(remove_transient:step:end-dt(2)),R_pseudo_2_filt{1,i}(remove_transient:step:end)/(1000*(dt(2)*Ts))*3600,...
        sz,Temp_mod_2{i}(remove_transient:step:end-dt(2)),'filled','LineWidth',2); 
end
xlabel('SoC [%]','interpreter','tex');
ylabel('DV [V/Ah]','interpreter','tex');
xlim([0 100]); ylim([0 1])
colormap(flipud(autumn));
c=colorbar('eastoutside'); caxis([10 32.5]);
c.Label.String = 'Temperature [°C]';
c.Label.Interpreter = 'tex';
set(findall(gcf,'-property','interpreter'),'interpreter','tex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')
set(findall(gcf,'-property','FontSize'),'FontSize',20);

%% Figure 11 - C/2
remove_transient = 50/Ts;
for k = 1:length(dt)
    f=figure; box on; hold all
    for i=1:length(R_pseudo_3{1,k})
        scatter(SoC_mod_3{i}(remove_transient:step:end-dt(k)),R_pseudo_3{1,k}{1,i}(remove_transient:step:end),...
            sz,Temp_mod_3{i}(remove_transient:step:end-dt(k)),'filled','LineWidth',2); 
    end
    xlabel('SoC [%]','interpreter','tex');
    ylabel('Z_{CHG} [m\Omega]','interpreter','tex');
    xlim([0 100]); ylim([0 0.25*dt(k)/100])
    colormap(flipud(autumn));
    c=colorbar('eastoutside'); caxis([10 32.5]);
    c.Label.String = 'Temperature [°C]';
    c.Label.Interpreter = 'tex';
    set(findall(gcf,'-property','interpreter'),'interpreter','tex')
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','tex')
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
end