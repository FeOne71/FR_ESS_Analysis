%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_DCIR_v2.m
% - structure.m + byd.m 로직을 적용한 DCIR 처리 (Ch9~Ch14용)
% - n1C/0.2C step 감지 및 윈도우 기반 SOC 계산
% - n1C step만 사용하여 저항 계산
% - Ch9, Ch10, Ch11, Ch12, Ch13, Ch14 채널 처리
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Paths and settings
parsedDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_v2';

if ~exist(saveDir,'dir'); mkdir(saveDir); end

channels = {'Ch9','Ch10','Ch11','Ch12','Ch13','Ch14'};
rpt_cycles = {'0cyc','200cyc','400cyc'};

%% 파라미터 설정 (structure.m과 동일)
I_1C = 64;              % 1C 전류값 [A]
tolC_main = 0.05;       % 전류 판정 허용오차 (5%)
dt_thresh_n1C = 0.2;    % n1C 샘플 평균 간격 조건 [초] (0.1초 간격)
dt_thresh_02C = 0.2;    % n0.2C 샘플 평균 간격 조건 [초] (0.1초 간격)
dI_min_A = 0.1;         % |dI|<이 값이면 스킵
markerList = {'o','s','^','d','v','>','<','p','h','x','+'}; % RPT별 마커

%% 출력 데이터 구조
dcir_soc_data = struct();

%% Color mapping for cycles
function color_val = getCycleColor(rpt_cycle)
    if strcmp(rpt_cycle,'0cyc')
        color_val = [32/255, 133/255, 78/255];  % green #20854E
    elseif strcmp(rpt_cycle,'200cyc')
        color_val = [0/255, 115/255, 194/255];  % blue #0073C2
    elseif strcmp(rpt_cycle,'400cyc')
        color_val = [205/255, 83/255, 76/255];  % red #CD534C   
    end
end

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    
    % 피겨 준비 - subplot 형태
    fig = figure('Name', sprintf('DCIR vs SOC - %s', channel), 'Position', [100 100 1200 800]);
    tlo = tiledlayout(fig, 2, 1, 'TileSpacing','compact','Padding','compact');
    
    % 충전 subplot
    nexttile; hold on; grid on; title(sprintf('%s - Charge (DCIR)', channel));
    xlabel('SOC (%)'); ylabel('R (m\Omega)');
    xlim([0 100]); xticks(0:10:100);
    ylim([0 4]); yticks(0:0.5:4);
    
    % 방전 subplot
    nexttile; hold on; grid on; title(sprintf('%s - Discharge (DCIR)', channel));
    xlabel('SOC (%)'); ylabel('R (m\Omega)');
    xlim([0 100]); xticks(0:10:100);
    ylim([0 4]); yticks(0:0.5:4);
    
    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        
        % Map to valid struct field name
        if strcmp(rpt_cycle, '0cyc')
            rpt_key = 'cyc0';
        elseif strcmp(rpt_cycle, '200cyc')
            rpt_key = 'cyc200';
        elseif strcmp(rpt_cycle, '400cyc')
            rpt_key = 'cyc400';
        else
            rpt_key = matlab.lang.makeValidName(rpt_cycle);
        end
        
        % Load parsed MAT file
        cyc_num = strrep(rpt_cycle, 'cyc', '');
        filename = sprintf('RPT%s_ch%s_parsed.mat', cyc_num, lower(channel(3:end)));
        filepath = fullfile(parsedDataPath, filename);
        if ~exist(filepath,'file')
            fprintf('Skip missing file: %s\n', filepath);
            continue;
        end
        
        load(filepath);  % Loads: data, intrim, pdata, channel_var, cyc_var
        fprintf('Loaded file: %s\n', filename);
        
        % 디버깅: 데이터 기본 정보
        fprintf('  Total steps: %d\n', length(pdata));
        
        % Convert pdata to structure.m format
        data = [];
        for i = 1:length(pdata)
            data(i).t = pdata(i).t;
            data(i).V = pdata(i).V;
            data(i).I = pdata(i).I;
            data(i).type = char(pdata(i).type);
        end
        
        %% 1) step 통계(dQ_Ah 포함) - structure.m 로직
        for k = 1:numel(data)
            Ik = data(k).I(:);
            tk = data(k).t(:);
            data(k).avg_I = mean(Ik);
            n = min(numel(tk), numel(Ik));
            if n >= 2
                data(k).dt_mean = mean(diff(tk));
                data(k).dur = tk(end) - tk(1);
                data(k).dQ_Ah = trapz(tk, Ik) / 3600;   % Ah
            else
                data(k).dt_mean = Inf;
                data(k).dur = 0;
                data(k).dQ_Ah = 0;
            end
        end
        
        %% 2) 방전 n1C 플래그 - structure.m 로직 (0.1초 간격)
        n1C_low = -(1+tolC_main) * I_1C;
        n1C_high = -(1-tolC_main) * I_1C;
        is_n1C = arrayfun(@(x) ...
            (x.avg_I >= n1C_low && x.avg_I <= n1C_high) && ...
            (x.dt_mean <= dt_thresh_n1C) && numel(x.t) >= 2, data);
        for k = 1:numel(data); data(k).n1C_flag = double(is_n1C(k)); end
        n1C_idx = find(is_n1C);
        
        % 디버깅: n1C 검출 결과
        fprintf('  Discharge n1C steps found: %d (indices: %s)\n', length(n1C_idx), mat2str(n1C_idx));
        
        %% 3) 충전 n1C 플래그 (0.1초 간격)
        charge_n1C_low = (1-tolC_main) * I_1C;
        charge_n1C_high = (1+tolC_main) * I_1C;
        is_charge_n1C = arrayfun(@(x) ...
            (x.avg_I >= charge_n1C_low && x.avg_I <= charge_n1C_high) && ...
            (x.dt_mean <= dt_thresh_n1C) && numel(x.t) >= 2, data);
        for k = 1:numel(data); data(k).charge_n1C_flag = double(is_charge_n1C(k)); end
        charge_n1C_idx = find(is_charge_n1C);
        
        % 디버깅: 충전 n1C 검출 결과
        fprintf('  Charge n1C steps found: %d (indices: %s)\n', length(charge_n1C_idx), mat2str(charge_n1C_idx));
        
        %% 4) 방전 0.2C 플래그 - structure.m 로직 (1.1초 간격)
        c02_low = -(0.2 + 0.2*tolC_main) * I_1C;   % 0.2C ±5%
        c02_high = -(0.2 - 0.2*tolC_main) * I_1C;
        is_c02_neg = arrayfun(@(x) ...
            (x.avg_I >= c02_low && x.avg_I <= c02_high) && ...
            (x.avg_I < 0) && (x.dt_mean <= dt_thresh_02C) && numel(x.t) >= 2, data);
        c02_neg_idx = find(is_c02_neg);
        
        % 디버깅: 방전 0.2C 검출 결과
        fprintf('  Discharge 0.2C steps found: %d (indices: %s)\n', length(c02_neg_idx), mat2str(c02_neg_idx));
        
        %% 5) 충전 0.2C 플래그 (1.1초 간격)
        charge_c02_low = (0.2 - 0.2*tolC_main) * I_1C;
        charge_c02_high = (0.2 + 0.2*tolC_main) * I_1C;
        is_charge_c02 = arrayfun(@(x) ...
            (x.avg_I >= charge_c02_low && x.avg_I <= charge_c02_high) && ...
            (x.avg_I > 0) && (x.dt_mean <= dt_thresh_02C) && numel(x.t) >= 2, data);
        charge_c02_idx = find(is_charge_c02);
        
        % 디버깅: 충전 0.2C 검출 결과
        fprintf('  Charge 0.2C steps found: %d (indices: %s)\n', length(charge_c02_idx), mat2str(charge_c02_idx));
        
        %% 6) 방전 윈도우 설정
        discharge_integ_set = [];
        if ~isempty(n1C_idx) && ~isempty(c02_neg_idx)
            first_n1C = n1C_idx(1);
            last_c02_neg = c02_neg_idx(end);
            if last_c02_neg < first_n1C
                last_c02_neg = n1C_idx(end);
            end
            N = numel(data);
            first_n1C = max(1, min(N, first_n1C));
            last_c02_neg = max(1, min(N, last_c02_neg));
            if last_c02_neg < first_n1C
                tmp = first_n1C; first_n1C = last_c02_neg; last_c02_neg = tmp;
            end
            discharge_integ_set = first_n1C:last_c02_neg;
        end
        
        % 디버깅: 방전 윈도우 설정 결과
        fprintf('  Discharge window: %s\n', mat2str(discharge_integ_set));
        
        %% 7) 충전 윈도우 설정
        charge_integ_set = [];
        if ~isempty(charge_n1C_idx) && ~isempty(charge_c02_idx)
            first_charge_n1C = charge_n1C_idx(1);
            last_charge_c02 = charge_c02_idx(end);
            if last_charge_c02 < first_charge_n1C
                last_charge_c02 = charge_n1C_idx(end);
            end
            N = numel(data);
            first_charge_n1C = max(1, min(N, first_charge_n1C));
            last_charge_c02 = max(1, min(N, last_charge_c02));
            if last_charge_c02 < first_charge_n1C
                tmp = first_charge_n1C; first_charge_n1C = last_charge_c02; last_charge_c02 = tmp;
            end
            charge_integ_set = first_charge_n1C:last_charge_c02;
        end
        
        % 디버깅: 충전 윈도우 설정 결과
        fprintf('  Charge window: %s\n', mat2str(charge_integ_set));
        
        %% 8) 방전 용량 계산
        discharge_Cn_Ah = 0;
        if ~isempty(discharge_integ_set)
            N = numel(data);
            discharge_in_win_mask = false(1, N);
            discharge_in_win_mask(discharge_integ_set) = true;
            discharge_mask = reshape([data.avg_I] < 0, 1, []);
            discharge_active_mask = discharge_in_win_mask & discharge_mask;
            
            if any(discharge_active_mask)
                discharge_Cn_Ah = abs(sum([data(discharge_active_mask).dQ_Ah]));
            end
        end
        
        %% 9) 충전 용량 계산
        charge_Cn_Ah = 0;
        if ~isempty(charge_integ_set)
            N = numel(data);
            charge_in_win_mask = false(1, N);
            charge_in_win_mask(charge_integ_set) = true;
            charge_mask = reshape([data.avg_I] > 0, 1, []);
            charge_active_mask = charge_in_win_mask & charge_mask;
            
            if any(charge_active_mask)
                charge_Cn_Ah = abs(sum([data(charge_active_mask).dQ_Ah]));
            end
        end
        
        %% 10) SOC 계산 - 충전/방전 구분
        N = numel(data);
        
        % 방전 SOC 계산 (100%에서 시작)
        if discharge_Cn_Ah > 0
            SOC_now = 100; started = false;
            for k = 1:N
                if ~discharge_in_win_mask(k)
                    data(k).discharge_SOC_start = NaN;
                    data(k).discharge_SOC_end = NaN;
                    continue
                end
                
                if ~started
                    SOC_now = 100; started = true;  % 100% SOC
                end
                data(k).discharge_SOC_start = SOC_now;
                
                if discharge_active_mask(k)   % 방전 step만 SOC 변화
                    dSOC = (data(k).dQ_Ah / discharge_Cn_Ah) * 100;   % 방전(dQ<0) → dSOC<0
                else
                    dSOC = 0;                                % 휴지/충전은 유지
                end
                
                SOC_now = SOC_now + dSOC;
                data(k).discharge_SOC_end = SOC_now;
            end
        end
        
        % 충전 SOC 계산 (0%에서 시작)
        if charge_Cn_Ah > 0
            SOC_now = 0; started = false;
            for k = 1:N
                if ~charge_in_win_mask(k)
                    data(k).charge_SOC_start = NaN;
                    data(k).charge_SOC_end = NaN;
                    continue
                end
                
                if ~started
                    SOC_now = 0; started = true;  % 0% SOC
                end
                data(k).charge_SOC_start = SOC_now;
                
                if charge_active_mask(k)   % 충전 step만 SOC 변화
                    dSOC = abs(data(k).dQ_Ah / charge_Cn_Ah) * 100;   % 충전 → dSOC>0
                else
                    dSOC = 0;                                % 휴지/방전은 유지
                end
                
                SOC_now = SOC_now + dSOC;
                data(k).charge_SOC_end = SOC_now;
            end
        end
        
        %% 11) 충전 DCIR step 사용하여 저항 계산
        charge_soc_vals = []; charge_R1_vals = []; charge_R5_vals = []; charge_R30_vals = []; charge_R60_vals = [];
        charge_dQ_vals = []; charge_abs_dQ_vals = [];
        
        % 디버깅: 충전 DCIR step 처리 시작
        fprintf('  Processing charge DCIR steps...\n');
        fprintf('  Charge n1C steps in window: %d out of %d\n', sum(ismember(charge_n1C_idx, charge_integ_set)), length(charge_n1C_idx));
        
        for k = 1:numel(data)
            % 충전 DCIR step만
            if ~isfield(data(k),'charge_n1C_flag') || data(k).charge_n1C_flag ~= 1, continue; end
            if ~isfield(data(k),'I') || ~isfield(data(k),'V') || numel(data(k).I)<2, continue; end
            
            % 디버깅: 첫 번째 충전 n1C step에서만 상세 정보 출력
            if k == charge_n1C_idx(1)
                fprintf('    Debug - First charge n1C step %d:\n', k);
                fprintf('      I range: %.2f to %.2f A\n', min(data(k).I), max(data(k).I));
                fprintf('      V range: %.3f to %.3f V\n', min(data(k).V), max(data(k).V));
                fprintf('      t range: %.2f to %.2f sec\n', min(data(k).t), max(data(k).t));
                fprintf('      dt_mean: %.3f sec\n', data(k).dt_mean);
            end
            
            % 시간/전압/전류 벡터
            t_vec = data(k).t(:);
            V_vec = data(k).V(:);
            I_vec = data(k).I(:);
            
            % 시간 벡터가 0부터 시작하도록 조정
            if ~isempty(t_vec) && t_vec(1) ~= 0
                t_vec = t_vec - t_vec(1);
            end
            
            % 1초, 5초, 30초, 60초 시점 인덱스 찾기
            idx_1s = find(t_vec >= 1, 1, 'first');
            idx_5s = find(t_vec >= 5, 1, 'first');
            idx_30s = find(t_vec >= 30, 1, 'first');
            idx_60s = find(t_vec > 59.8, 1, 'last');  % 60초 대신 59.8초 이후 마지막 값
            
            % 디버깅: 첫 번째 n1C step에서만 출력
            if k == charge_n1C_idx(1)
                fprintf('    Debug - Time indices for step %d:\n', k);
                fprintf('      t_vec range: %.2f to %.2f sec\n', min(t_vec), max(t_vec));
                fprintf('      idx_1s: %d, idx_5s: %d, idx_30s: %d, idx_60s: %d\n', ...
                    idx_1s, idx_5s, idx_30s, idx_60s);
                fprintf('      V_vec range: %.3f to %.3f V\n', min(V_vec), max(V_vec));
                fprintf('      I_vec range: %.1f to %.1f A\n', min(I_vec), max(I_vec));
            end
            
            % 인덱스가 존재하는지 확인
            if isempty(idx_1s) || isempty(idx_5s) || isempty(idx_30s) || isempty(idx_60s), continue; end
            
            % 각 시점의 저항 계산
            V_start = V_vec(1);
            R1 = abs((V_vec(idx_1s) - V_start) / I_vec(idx_1s));
            R5 = abs((V_vec(idx_5s) - V_start) / I_vec(idx_5s));
            R30 = abs((V_vec(idx_30s) - V_start) / I_vec(idx_30s));
            R60 = abs((V_vec(idx_60s) - V_start) / I_vec(idx_60s));
            
            % 유효성 검사
            if ~isfinite(R1) || ~isfinite(R5) || ~isfinite(R30) || ~isfinite(R60)
                if k == charge_n1C_idx(1)
                    fprintf('    Debug - Invalid resistance values for step %d: R1=%.4f, R5=%.4f, R30=%.4f, R60=%.4f\n', k, R1, R5, R30, R60);
                end
                continue; 
            end
            if abs(I_vec(idx_1s)) < dI_min_A || abs(I_vec(idx_5s)) < dI_min_A || abs(I_vec(idx_30s)) < dI_min_A || abs(I_vec(idx_60s)) < dI_min_A
                if k == charge_n1C_idx(1)
                    fprintf('    Debug - Current too small for step %d: I1s=%.4f, I5s=%.4f, I30s=%.4f, I60s=%.4f (min=%.4f)\n', k, abs(I_vec(idx_1s)), abs(I_vec(idx_5s)), abs(I_vec(idx_30s)), abs(I_vec(idx_60s)), dI_min_A);
                end
                continue; 
            end
            
            % 충전 SOC 값
            if isfield(data(k),'charge_SOC_start') && ~isnan(data(k).charge_SOC_start)
                charge_soc_vals(end+1,1) = data(k).charge_SOC_start;
                charge_R1_vals(end+1,1) = R1;
                charge_R5_vals(end+1,1) = R5;
                charge_R30_vals(end+1,1) = R30;
                charge_R60_vals(end+1,1) = R60;
                charge_dQ_vals(end+1,1) = data(k).dQ_Ah;
                charge_abs_dQ_vals(end+1,1) = abs(data(k).dQ_Ah);
            else
                % 디버깅: SOC가 NaN인 경우
                if k == charge_n1C_idx(1)
                    fprintf('    Debug - SOC is NaN for step %d: charge_SOC_start = %s\n', k, mat2str(data(k).charge_SOC_start));
                end
            end
        end
        
        % 디버깅: 충전 DCIR 결과
        fprintf('  Charge DCIR results: %d valid steps\n', length(charge_soc_vals));
        
        %% 12) 방전 DCIR step 사용하여 저항 계산
        discharge_soc_vals = []; discharge_R1_vals = []; discharge_R5_vals = []; discharge_R30_vals = []; discharge_R60_vals = [];
        discharge_dQ_vals = []; discharge_abs_dQ_vals = [];
        
        % 디버깅: 방전 DCIR step 처리 시작
        fprintf('  Processing discharge DCIR steps...\n');
        fprintf('  Discharge n1C steps in window: %d out of %d\n', sum(ismember(n1C_idx, discharge_integ_set)), length(n1C_idx));
        
        for k = 1:numel(data)
            % 방전 DCIR step만
            if ~isfield(data(k),'n1C_flag') || data(k).n1C_flag ~= 1, continue; end
            if ~isfield(data(k),'I') || ~isfield(data(k),'V') || numel(data(k).I)<2, continue; end
            
            % 디버깅: 첫 번째 방전 n1C step에서만 상세 정보 출력
            if k == n1C_idx(1)
                fprintf('    Debug - First discharge n1C step %d:\n', k);
                fprintf('      I range: %.2f to %.2f A\n', min(data(k).I), max(data(k).I));
                fprintf('      V range: %.3f to %.3f V\n', min(data(k).V), max(data(k).V));
                fprintf('      t range: %.2f to %.2f sec\n', min(data(k).t), max(data(k).t));
                fprintf('      dt_mean: %.3f sec\n', data(k).dt_mean);
            end
            
            % 시간/전압/전류 벡터
            t_vec = data(k).t(:);
            V_vec = data(k).V(:);
            I_vec = data(k).I(:);
            
            % 시간 벡터가 0부터 시작하도록 조정
            if ~isempty(t_vec) && t_vec(1) ~= 0
                t_vec = t_vec - t_vec(1);
            end
            
            % 1초, 5초, 30초, 60초 시점 인덱스 찾기
            idx_1s = find(t_vec >= 1, 1, 'first');
            idx_5s = find(t_vec >= 5, 1, 'first');
            idx_30s = find(t_vec >= 30, 1, 'first');
            idx_60s = find(t_vec > 59.8, 1, 'last');  % 60초 대신 59.8초 이후 마지막 값
            
            % 디버깅: 첫 번째 n1C step에서만 출력
            if k == n1C_idx(1)
                fprintf('    Debug - Time indices for step %d:\n', k);
                fprintf('      t_vec range: %.2f to %.2f sec\n', min(t_vec), max(t_vec));
                fprintf('      idx_1s: %d, idx_5s: %d, idx_30s: %d, idx_60s: %d\n', ...
                    idx_1s, idx_5s, idx_30s, idx_60s);
                fprintf('      V_vec range: %.3f to %.3f V\n', min(V_vec), max(V_vec));
                fprintf('      I_vec range: %.1f to %.1f A\n', min(I_vec), max(I_vec));
            end
            
            % 인덱스가 존재하는지 확인
            if isempty(idx_1s) || isempty(idx_5s) || isempty(idx_30s) || isempty(idx_60s), continue; end
            
            % 각 시점의 저항 계산
            V_start = V_vec(1);
            R1 = abs((V_vec(idx_1s) - V_start) / I_vec(idx_1s));
            R5 = abs((V_vec(idx_5s) - V_start) / I_vec(idx_5s));
            R30 = abs((V_vec(idx_30s) - V_start) / I_vec(idx_30s));
            R60 = abs((V_vec(idx_60s) - V_start) / I_vec(idx_60s));
            
            % 유효성 검사
            if ~isfinite(R1) || ~isfinite(R5) || ~isfinite(R30) || ~isfinite(R60)
                if k == n1C_idx(1)
                    fprintf('    Debug - Invalid resistance values for step %d: R1=%.4f, R5=%.4f, R30=%.4f, R60=%.4f\n', k, R1, R5, R30, R60);
                end
                continue; 
            end
            if abs(I_vec(idx_1s)) < dI_min_A || abs(I_vec(idx_5s)) < dI_min_A || abs(I_vec(idx_30s)) < dI_min_A || abs(I_vec(idx_60s)) < dI_min_A
                if k == n1C_idx(1)
                    fprintf('    Debug - Current too small for step %d: I1s=%.4f, I5s=%.4f, I30s=%.4f, I60s=%.4f (min=%.4f)\n', k, abs(I_vec(idx_1s)), abs(I_vec(idx_5s)), abs(I_vec(idx_30s)), abs(I_vec(idx_60s)), dI_min_A);
                end
                continue; 
            end
            
            % 방전 SOC 값
            if isfield(data(k),'discharge_SOC_start') && ~isnan(data(k).discharge_SOC_start)
                discharge_soc_vals(end+1,1) = data(k).discharge_SOC_start;
                discharge_R1_vals(end+1,1) = R1;
                discharge_R5_vals(end+1,1) = R5;
                discharge_R30_vals(end+1,1) = R30;
                discharge_R60_vals(end+1,1) = R60;
                discharge_dQ_vals(end+1,1) = data(k).dQ_Ah;
                discharge_abs_dQ_vals(end+1,1) = abs(data(k).dQ_Ah);
            else
                % 디버깅: SOC가 NaN인 경우
                if k == n1C_idx(1)
                    fprintf('    Debug - SOC is NaN for step %d: discharge_SOC_start = %s\n', k, mat2str(data(k).discharge_SOC_start));
                end
            end
        end
        
        % 디버깅: 방전 DCIR 결과
        fprintf('  Discharge DCIR results: %d valid steps\n', length(discharge_soc_vals));
        
        %% 13) 충전 DCIR 데이터 저장 및 시각화
        if ~isempty(charge_soc_vals)
            [charge_soc_sorted, charge_idx] = sort(charge_soc_vals, 'ascend');  % 0%에서 시작
            charge_R1_sorted = charge_R1_vals(charge_idx);
            charge_R5_sorted = charge_R5_vals(charge_idx);
            charge_R30_sorted = charge_R30_vals(charge_idx);
            charge_R60_sorted = charge_R60_vals(charge_idx);
            charge_dQ_sorted = charge_dQ_vals(charge_idx);
            charge_abs_dQ_sorted = charge_abs_dQ_vals(charge_idx);
            
            % Table 생성 (mOhm 단위로 저장)
            charge_Cn_Ah_vec = repmat(charge_Cn_Ah, length(charge_soc_sorted), 1);
            charge_table = table(charge_soc_sorted, charge_dQ_sorted, charge_Cn_Ah_vec, ...
                                1e3*charge_R1_sorted, 1e3*charge_R5_sorted, 1e3*charge_R30_sorted, 1e3*charge_R60_sorted, ...
                                'VariableNames', {'SOC', 'dQ_Ah', 'Cn_Ah', 'R1_mOhm', 'R5_mOhm', 'R30_mOhm', 'R60_mOhm'});
            
            % 데이터 저장
            dcir_soc_data.(channel).(rpt_key).charge_table = charge_table;
            
            % 충전 DCIR 시각화 (subplot 1) - R30 사용 (mOhm)
            nexttile(1);
            color_val = getCycleColor(rpt_cycle);
            mk = markerList{mod(rpt_idx-1,numel(markerList))+1};
            plot(charge_soc_sorted, 1e3*charge_R30_sorted, ['-' mk], ...
                'LineWidth',1.5,'MarkerSize',8, 'Color', color_val, ...
                'DisplayName', sprintf('%s', rpt_cycle));
            
            % 충전 DCIR 디버깅 표 출력
            fprintf('\n=== %s - Charge DCIR ===\n', filename);
            fprintf('SOC | dQ_Ah | Cn_Ah | R1 | R5 | R30 | R60\n');
            for i = 1:height(charge_table)
                fprintf('%.2f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f\n', ...
                    charge_table.SOC(i), charge_table.dQ_Ah(i), charge_table.Cn_Ah(i), ...
                    charge_table.R1_mOhm(i), charge_table.R5_mOhm(i), charge_table.R30_mOhm(i), charge_table.R60_mOhm(i));
            end
        end
        
        %% 14) 방전 DCIR 데이터 저장 및 시각화
        if ~isempty(discharge_soc_vals)
            [discharge_soc_sorted, discharge_idx] = sort(discharge_soc_vals, 'descend');  % 100%에서 시작
            discharge_R1_sorted = discharge_R1_vals(discharge_idx);
            discharge_R5_sorted = discharge_R5_vals(discharge_idx);
            discharge_R30_sorted = discharge_R30_vals(discharge_idx);
            discharge_R60_sorted = discharge_R60_vals(discharge_idx);
            discharge_dQ_sorted = discharge_dQ_vals(discharge_idx);
            discharge_abs_dQ_sorted = discharge_abs_dQ_vals(discharge_idx);
            
            % Table 생성 (mOhm 단위로 저장)
            discharge_Cn_Ah_vec = repmat(discharge_Cn_Ah, length(discharge_soc_sorted), 1);
            discharge_table = table(discharge_soc_sorted, discharge_dQ_sorted, discharge_Cn_Ah_vec, ...
                                   1e3*discharge_R1_sorted, 1e3*discharge_R5_sorted, 1e3*discharge_R30_sorted, 1e3*discharge_R60_sorted, ...
                                   'VariableNames', {'SOC', 'dQ_Ah', 'Cn_Ah', 'R1_mOhm', 'R5_mOhm', 'R30_mOhm', 'R60_mOhm'});
            
            % 데이터 저장
            dcir_soc_data.(channel).(rpt_key).discharge_table = discharge_table;
            
            % 방전 DCIR 시각화 (subplot 2) - R30 사용 (mOhm)
            nexttile(2);
            color_val = getCycleColor(rpt_cycle);
            mk = markerList{mod(rpt_idx-1,numel(markerList))+1};
            plot(discharge_soc_sorted, 1e3*discharge_R30_sorted, ['-' mk], ...
                'LineWidth',1.5,'MarkerSize',8, 'Color', color_val, ...
                'DisplayName', sprintf('%s', rpt_cycle));
            
            % 방전 DCIR 디버깅 표 출력
            fprintf('\n=== %s - Discharge DCIR ===\n', filename);
            fprintf('SOC | dQ_Ah | Cn_Ah | R1 | R5 | R30 | R60\n');
            for i = 1:height(discharge_table)
                fprintf('%.2f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f\n', ...
                    discharge_table.SOC(i), discharge_table.dQ_Ah(i), discharge_table.Cn_Ah(i), ...
                    discharge_table.R1_mOhm(i), discharge_table.R5_mOhm(i), discharge_table.R30_mOhm(i), discharge_table.R60_mOhm(i));
            end
        end
    end
    
    % 범례 설정
    nexttile(1); legend('Location','northeast','Interpreter','none');
    nexttile(2); legend('Location','northeast','Interpreter','none');
    
    % Save figure
    savefig(fig, fullfile(saveDir, sprintf('DCIR_%s_integrated_v2.fig', channel)));
    close(fig);
end

% Save MAT file
matFile = fullfile(saveDir, 'DCIR_SOC_data_9to14_v2.mat');
save(matFile, 'dcir_soc_data');
fprintf('\nSaved: %s\n', matFile);

fprintf('\n=== RPT_DCIR_v2.m 완료 ===\n');

% Run RPT_DCIR_1s.m for 1s-sampled channels
RPT_DCIR_1s_v2