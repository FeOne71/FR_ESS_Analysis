%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright @ 2023 Stanford Energy Control Laboratory 
% (PI: Simona Onori, sonori@stanford.edu), 
% Stanford University. All Rights Reserved. 
% Developed by Pietro Bosoni (pbosoni@stanford.edu) and 
% Dr. Gabriele Pozzato (gpozzato@stanford.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Curr_Idles = Detect_Curr_Idles(time, curr, minimum_counter, time_thr)

Curr_Idles      = {};
null_Curr       = curr==0;

i=1;
while i<=length(null_Curr)
    counter=0;
    while (null_Curr(i+counter)==1) && (i+counter<length(null_Curr)) && (time(i+counter+1)-time(i+counter)< time_thr)
        counter = counter+1;
    end
    if counter>= minimum_counter
        Curr_Idles{end+1} = [time(i);time(i+counter)];
    end
    i = i+counter+1;
end

if isempty(Curr_Idles)
   Curr_Idles{end+1}= [0;0];
end
end