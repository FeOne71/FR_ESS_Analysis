function t_sec = feat_prepare_time(t)
%FEAT_PREPARE_TIME Convert time vector to seconds.
if isa(t, 'duration')
    t_sec = seconds(t);
else
    t_sec = t;
end
end
