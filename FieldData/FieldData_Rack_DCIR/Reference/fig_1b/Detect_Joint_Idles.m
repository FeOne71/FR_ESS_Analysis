%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright @ 2023 Stanford Energy Control Laboratory 
% (PI: Simona Onori, sonori@stanford.edu), 
% Stanford University. All Rights Reserved. 
% Developed by Pietro Bosoni (pbosoni@stanford.edu) and
% Dr. Gabriele Pozzato (gpozzato@stanford.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Joint_Idles = Detect_Joint_Idles(time, time_thr)

Joint_Idles = {};

sampling_periods = diff(time);

for  i =1:length(sampling_periods)
    if sampling_periods(i)>time_thr
        Joint_Idles{end+1} = [time(i);time(i+1)];
    end
end
end