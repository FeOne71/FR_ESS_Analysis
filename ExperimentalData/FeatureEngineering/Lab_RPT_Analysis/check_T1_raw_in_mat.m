%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_T1_raw_in_mat.m
% RPT_VQ_grid.mat 내 C-rate별 충/방전 세그먼트의 T1_raw, t_raw 끝부분 확인
% (그래프 끝에서 온도 급상승이 나오는지 MAT 원본에는 없는지 검사)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_mat = fullfile(baseDir, 'RPT', 'Postprocessing', 'OCV_integrated', 'RPT_VQ_grid.mat');
if ~exist(path_mat, 'file')
    fprintf('File not found: %s\n', path_mat);
    return;
end
S = load(path_mat, 'RPT_VQ_grid');
RPT = S.RPT_VQ_grid;
cyc_fields = fieldnames(RPT);
channels = fieldnames(RPT.(cyc_fields{1}));
channels = channels(strncmp(channels, 'Ch', 2));
crate_segments = {'c01_charge','c01_discharge','c05_charge','c05_discharge','c1_charge','c1_discharge','c2_charge','c2_discharge','c3_charge','c3_discharge'};
n_tail = 30;  % 끝에서 몇 점 확인할지

fprintf('=== Checking RPT_VQ_grid.mat T1_raw / t_raw ends (Last %d points) ===\n', n_tail);
fprintf('Path: %s\n\n', path_mat);

for cy = 1:min(2, numel(cyc_fields))   % 처음 2 사이클만 샘플
    cyc_name = cyc_fields{cy};
    for seg = 1:numel(crate_segments)
        seg_name = crate_segments{seg};
        fprintf('--- %s / %s ---\n', cyc_name, seg_name);
        ch_ok = 0;
        for ch = 1:numel(channels)
            ch_name = channels{ch};
            if ~isfield(RPT.(cyc_name), ch_name), continue; end
            C = RPT.(cyc_name).(ch_name);
            if ~isfield(C, seg_name), continue; end
            D = C.(seg_name);
            if ~isfield(D, 'T1_raw')
                if ch == 1, fprintf('  (No T1_raw)\n'); end
                continue;
            end
            T1 = D.T1_raw(:);
            if isfield(D, 't_raw'), tr = D.t_raw(:); else, tr = []; end
            n = numel(T1);
            if n < n_tail, idx = 1:n; else, idx = (n - n_tail + 1) : n; end
            T_tail = T1(idx);
            if ch == 1 || ch_ok == 0
                fprintf('  Ch%s: n=%d, Last %d pts T1_raw min=%.2f max=%.2f, last=%.2f\n', ...
            ch_name, n, numel(T_tail), ...
            min(T_tail), max(T_tail), T_tail(end));
            end
            ch_ok = ch_ok + 1;
        end
        if ch_ok == 0, fprintf('  (No Data)\n'); end
        fprintf('\n');
    end
end
fprintf('=== END ===\n');
