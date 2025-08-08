%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright @ 2023 Stanford Energy Control Laboratory 
% (PI: Simona Onori, sonori@stanford.edu), 
% Stanford University. All Rights Reserved. 
% Developed by Pietro Bosoni and Dr. Gabriele Pozzato 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

folders_ch  = {'1' '3' '5' '7' '9' '11' '13' '15'};
folders_dis = {'2' '4' '6' '8' '10' '12' '14' '16'};

Total_Driving_Time  = 0;
Total_Charging_Time = 0;

% Create one structure for each charging and driving folder
Charging_Times      = {};
Driving_Times       = {};

%% Compute Charging Time for each folder 
for h=1:2:15
    selected_folder = h;
    
    ind_ch  = find(ismember(folders_ch, num2str(selected_folder)));
    ind_dis = find(ismember(folders_dis,num2str(selected_folder)));
   
    if ~isempty(ind_ch)
        load(['../../Charge/Folder' num2str(folders_ch{ind_ch}) '/Raw.mat']);
        time_0      = (Raw.TimeCurr(:)-Raw.TimeCurr(1));
        epoch_vec   = (Raw.Epoch(1)+(Raw.TimeCurr(1)-Raw.TimeEpoch(1))) + time_0;
        time_vector = epoch_vec - epoch_vec(1);
        
        savedates{str2num(folders_ch{ind_ch})}.start = datetime(epoch_vec(1),'ConvertFrom','epochtime');
        savedates{str2num(folders_ch{ind_ch})}.end   = datetime(epoch_vec(end),'ConvertFrom','epochtime');

        Charging_Times{h}.epoch(1) = epoch_vec(1);
        Charging_Times{h}.epoch(2) = epoch_vec(end);
        clear epoch_vec
    end
   
    % Parameters
    Fs               = 0.01;                                                   % Signal sampling frequency
    time_thr         = Fs*2;                                                   % [sec] time treshold to recognize events splitted by a some NaN
    minimum_counter  = 1/Fs*60;                                                % [Fs*sec] discard all null current segments longer than 'minimum_counter', now 60 sec.
    
    Joint_Idles = Detect_Joint_Idles(time_vector, time_thr);
    Curr_Idles  = Detect_Curr_Idles(time_vector, Raw.Curr, minimum_counter, time_thr);
    
    % Compute Charging Time for each event 
    Idle_Time  = 0;
    Curr_Idle  = 0;
    Joint_Idle = 0;

    for i=1:length(Joint_Idles)
        Idle_Time       = Idle_Time + Joint_Idles{i}(2)-Joint_Idles{i}(1);
        Joint_Idle      = Joint_Idle + Joint_Idles{i}(2)-Joint_Idles{i}(1);

    end
    for i=1:length(Curr_Idles)
        Idle_Time       = Idle_Time + Curr_Idles{i}(2)-Curr_Idles{i}(1);
        Curr_Idle       = Curr_Idle + Curr_Idles{i}(2)-Curr_Idles{i}(1);
    end

    Charging_Times{h}.Charge = (time_vector(end)-time_vector(1))-Idle_Time;
    Charging_Times{h}.Joint_Idle   = Joint_Idle;
    Charging_Times{h}.Curr_Idle    = Curr_Idle;


    Total_Charging_Time = Total_Charging_Time + Charging_Times{h}.Charge;
end

%% Compute Driving Time for each folder
for h=2:2:16
    selected_folder = h;
    
    ind_ch  = find(ismember(folders_ch, num2str(selected_folder)));
    ind_dis = find(ismember(folders_dis,num2str(selected_folder)));
   
    if  ~isempty(ind_dis)
        load(['../../Drive/Folder' num2str(folders_dis{ind_dis}) '/Raw.mat']);
        time_0      = (Raw.TimeCurr(:)-Raw.TimeCurr(1));
        epoch_vec   = (Raw.Epoch(1)+(Raw.TimeCurr(1)-Raw.TimeEpoch(1))) + time_0;
        time_vector = epoch_vec - epoch_vec(1);
    
        savedates{str2num(folders_ch{ind_dis})}.start = datetime(epoch_vec(1),'ConvertFrom','epochtime');
        savedates{str2num(folders_ch{ind_dis})}.end   = datetime(epoch_vec(end),'ConvertFrom','epochtime');

        Driving_Times{h}.epoch(1) = epoch_vec(1);
        Driving_Times{h}.epoch(2) = epoch_vec(end);

        clear epoch_vec
    end
    
    % Parameters
    Fs               = 0.01;                                                   % Signal sampling frequency
    time_thr         = Fs*2;
    minimum_counter  = 1/Fs*60;                                                % Fs*sec.
    
    Joint_Idles = Detect_Joint_Idles(time_vector, time_thr);
    Curr_Idles  = Detect_Curr_Idles(time_vector, Raw.Curr, minimum_counter, time_thr);
    
    % Compute Driving Time for each event    
    Idle_Time  = 0;
    Curr_Idle  = 0;
    Joint_Idle = 0;

    for i=1:length(Joint_Idles)
        Idle_Time       = Idle_Time + Joint_Idles{i}(2)-Joint_Idles{i}(1);
        Joint_Idle      = Joint_Idle + Joint_Idles{i}(2)-Joint_Idles{i}(1);

    end
    for i=1:length(Curr_Idles)
        Idle_Time       = Idle_Time + Curr_Idles{i}(2)-Curr_Idles{i}(1);
        Curr_Idle       = Curr_Idle + Curr_Idles{i}(2)-Curr_Idles{i}(1);

    end

    Driving_Times{h}.Drive = (time_vector(end)-time_vector(1))-Idle_Time;
    Driving_Times{h}.Joint_Idle   = Joint_Idle;
    Driving_Times{h}.Curr_Idle    = Curr_Idle;

    Total_Driving_Time = Total_Driving_Time + Driving_Times{h}.Drive;
end

%% Compute the total logged time NOT considering the NaN as logged events
Total_Time     = 0;
Total_Driving  = 0;
Total_Charging = 0;
Total_Idle     = 0;

for i=1:2:15
    Total_Time     = Total_Time + Charging_Times{i}.Charge + Driving_Times{i+1}.Drive + Charging_Times{i}.Curr_Idle + Driving_Times{i+1}.Curr_Idle;
    Total_Driving  = Total_Driving + Driving_Times{i+1}.Drive;
    Total_Charging = Total_Charging + Charging_Times{i}.Charge;
    Total_Idle     = Total_Idle+Charging_Times{i}.Curr_Idle + Driving_Times{i+1}.Curr_Idle;
end

%% Total time
initial_epoch = [];
final_epoch   = [];
logged_time   = [];

for i=1:2:15
    initial_epoch(i) = min(Charging_Times{i}.epoch(1),Driving_Times{i+1}.epoch(1)); 
    final_epoch(i)   = max(Charging_Times{i}.epoch(2),Driving_Times{i+1}.epoch(2));
    logged_time(i)   = final_epoch(i) - initial_epoch(i);
end

Total_Time      = sum(logged_time);
Total_Idle_Time = Total_Time- Total_Driving_Time-Total_Charging_Time;

%% Pie Plot
explode = [1 1 0];
figure
p      = pie([Total_Charging_Time, Total_Driving_Time, Total_Idle_Time], explode);
labels = {'Charging Time', 'Driving Time', 'Idle Time'};
lgd    = legend(labels, "Location","east", 'FontSize', 16);
Ptext  = findobj(p, 'Type', 'text');
set(findall(gcf,'-property','FontSize'),'FontSize',20);