%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright @ 2023 Stanford Energy Control Laboratory 
% (PI: Simona Onori, sonori@stanford.edu), 
% Stanford University. All Rights Reserved. 
% Developed by Dr. Gabriele Pozzato 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%% Load folder #9 for current quantization
load('../../Charge/Folder9/Raw.mat');
time_0 = (Raw.TimeCurr(:)-Raw.TimeCurr(1));
epoch_vec = (Raw.Epoch(1)+(Raw.TimeCurr(1)-Raw.TimeEpoch(1))) + time_0;
time_vec_ch = epoch_vec - epoch_vec(1);

%% Figure 10
% Time offset in seconds [s]
time_offset = 179826;

% Plot
figure; hold on; box on; 
plot(time_vec_ch-time_offset,Raw.Curr,'linewidth',2)
xlim([0 24]); ylim([-2.1 0.1]); xlabel('Time [s]'); ylabel('Current [A]')
set(findall(gcf,'-property','FontSize'),'FontSize',20);