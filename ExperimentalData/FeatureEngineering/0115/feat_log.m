function feat_log(log_path, msg)
%FEAT_LOG Append message to log file with timestamp.
if nargin < 2
    return;
end
ts = datestr(now, 'yyyy-mm-dd HH:MM:SS');
line = sprintf('[%s] %s\n', ts, msg);
fid = fopen(log_path, 'a');
if fid > 0
    fprintf(fid, '%s', line);
    fclose(fid);
end
end
