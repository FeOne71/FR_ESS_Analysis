% parsing code for practice

% structure
% 1. config
% 2. load
% 3. change variable type
% 4. parse and reconstruct
% 5. save

% may add
% looping over a folder (dir and regexp)


%% config
I_1C = 55.6; %[A]
fullpath_load = 'ex_data_rOCV';

%% load
data = readtable(fullpath_load);


%% change variable type
intrim.t = data{:,7};
intrim.V = data{:,8};
intrim.I = data{:,9};

plot(data{:,9})

%% parse and reconstruct

% check stepindex
plot(data.StepIndex) % (1) repeating (2) single pulse has multple steps

% assign CRDR types



