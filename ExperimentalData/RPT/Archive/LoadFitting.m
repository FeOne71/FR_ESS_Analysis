%% Operation Load Fitting

clc; close all; clear;

folderPath = 'C:\Users\Chulwon Jung\내 드라이브\JCW_GDrive\Experiment Data\RPT';
saveFolder = 'C:\Users\Chulwon Jung\내 드라이브\JCW_GDrive\Experiment Data\RPT\FitResult';

load(fullfile(folderPath, 'selected_profiles.mat'));       
load(fullfile(folderPath, 'OCV_interpolation_functions.mat'));

% 