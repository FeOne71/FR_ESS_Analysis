clear; clc; close all

data_folder = "G:\공유 드라이브\BSL_Data4 (1)\HNE_reduced_RPT_0812_processed\processed";
splitPath = split(data_folder, filesep);
index = find(strcmp('processed', splitPath), 1);
splitPath{index} = 'processed_DCIR';

save_path = strjoin(splitPath, filesep);
if ~exist(save_path, 'dir'); mkdir(save_path); end

files = dir(fullfile(data_folder, '*R*.mat'));
I_1C = 16;

%% 사용자 설정
SOC0_default   = 100;   % 정보용(계산에는 사용하지 않음)
tolC_main      = 0.05;  % 1C/0.2C 전류 판정 허용오차(±5%)
dt_thresh      = 2;     % 샘플 평균 간격 조건(2초 미만)
Vmin           = 2.0;   
time_overrun_s = 20;   

for fi = 1:length(files)
    fullpath_now = fullfile(data_folder, files(fi).name);
    [~, filename, ~] = fileparts(files(fi).name);
    load(fullpath_now, 'data');   % data struct array

    if isfield(data(1), 'SOC0') && ~isempty(data(1).SOC0)
        SOC0_meta = data(1).SOC0;
    else
        SOC0_meta = SOC0_default;
    end

    %% 1) step 통계(dQ_Ah 포함)
    for k = 1:numel(data)
        Ik = data(k).I(:);
        tk = data(k).t(:);
        data(k).avg_I = mean(Ik);
        n = min(numel(tk), numel(Ik));
        if n >= 2
            data(k).dt_mean = mean(diff(tk));
            data(k).dur     = tk(end) - tk(1);
            data(k).dQ_Ah   = trapz(tk, Ik) / 3600;   
            data(k).dt_mean = Inf;
            data(k).dur     = 0;
            data(k).dQ_Ah   = 0;
        end
    end

    %% 2) n1C 플래그
    n1C_low  = -(1+tolC_main) * I_1C;
    n1C_high = -(1-tolC_main) * I_1C;
    is_n1C = arrayfun(@(x) ...
        (x.avg_I >= n1C_low && x.avg_I <= n1C_high) && ...
        (x.dt_mean < dt_thresh) && numel(x.t) >= 2, data);
    for k = 1:numel(data); data(k).n1C_flag = double(is_n1C(k)); end
    n1C_idx = find(is_n1C);

    if isempty(n1C_idx)
        warning('(%s) n1C 구간 없음 → 스킵', filename);
        save(fullfile(save_path, [filename '_SOC.mat']), 'data', 'I_1C');
        continue
    end

    %% 3) 0.2C 방전 플래그
    c02_low  = -(0.2 + 0.2*tolC_main) * I_1C;   % 0.2C ±5%
    c02_high = -(0.2 - 0.2*tolC_main) * I_1C;
    is_c02_neg = arrayfun(@(x) ...
        (x.avg_I >= c02_low && x.avg_I <= c02_high) && ...
        (x.avg_I < 0) && (x.dt_mean < dt_thresh) && numel(x.t) >= 2, data);
    c02_neg_idx = find(is_c02_neg);

    %% 4) 적산 시작/종료 step  (첫 n1C ~ 마지막 -0.2C, 없으면 마지막 n1C)
    first_n1C = n1C_idx(1);
    if ~isempty(c02_neg_idx)
        last_c02_neg = c02_neg_idx(end);
        if last_c02_neg < first_n1C
            last_c02_neg = n1C_idx(end);
        end
    else
        last_c02_neg = n1C_idx(end);
        warning('(%s) 0.2C- 없음 → 마지막 n1C까지', filename);
    end


    N = numel(data);
    first_n1C    = max(1, min(N, first_n1C));
    last_c02_neg = max(1, min(N, last_c02_neg));
    if last_c02_neg < first_n1C
        warning('(%s) last_c02_neg < first_n1C → 교환', filename);
        tmp = first_n1C; first_n1C = last_c02_neg; last_c02_neg = tmp;
    end

    integ_set = first_n1C:last_c02_neg;        
    if isempty(integ_set)
        warning('(%s) integ_set 비어 있음 → 스킵', filename);
        save(fullfile(save_path, [filename '_SOC.mat']), 'data', 'I_1C');
        continue
    end

    %% 5) 용량 Cn_Ah = 윈도우 내 "방전 step"의 dQ_Ah 누적 절대값
    in_win_mask     = false(1, N);    in_win_mask(integ_set) = true;        
    discharge_mask  = reshape([data.avg_I] < 0, 1, []);                      
    active_mask     = in_win_mask & discharge_mask;                          

    if ~any(active_mask)
        warning('(%s) 윈도우 내 방전 step이 없습니다 → 스킵', filename);
        save(fullfile(save_path, [filename '_SOC.mat']), 'data', 'I_1C');
        continue
    end

    Cn_Ah = abs(sum([data(active_mask).dQ_Ah]));  % 누적 방전 용량(Ah)
    if ~isfinite(Cn_Ah) || Cn_Ah <= 0
        warning('(%s) Cn_Ah 비정상: %.6f Ah → 스킵', filename, Cn_Ah);
        save(fullfile(save_path, [filename '_SOC.mat']), 'data', 'I_1C');
        continue
    end

    %% 6)  SOC 계산 (윈도우 입장 시점 100%, active만 변화)
    SOC_now = NaN; started = false;
    for k = 1:N
        if ~in_win_mask(k)
            data(k).dSOC_pct  = NaN;
            data(k).SOC_start = NaN;
            data(k).SOC_end   = NaN;
            continue
        end

        if ~started
            SOC_now = 1; started = true;  
        end
        data(k).SOC_start = SOC_now;

        if active_mask(k)   % 방전 step만 SOC 변화
            dSOC = (data(k).dQ_Ah / Cn_Ah) * 100;   % 방전(dQ<0) → dSOC<0
        else
            dSOC = 0;                                % 휴지/충전은 유지
        end

        data(k).dSOC_pct = dSOC;
        SOC_now = SOC_now + dSOC;
        data(k).SOC_end = SOC_now;
    end


    total_dQ_Ah_discharge = sum([data(active_mask).dQ_Ah]);   
    total_dQ_Ah_net       = sum([data(integ_set).dQ_Ah]);     
    total_dSOC            = (total_dQ_Ah_discharge / Cn_Ah) * 100;
    SOC_after             = 100 + total_dSOC;

    summary = struct();
    summary.filename        = filename;
    summary.n1C_idx         = n1C_idx;
    summary.first_n1C       = first_n1C;
    summary.c02_neg_idx     = c02_neg_idx;
    summary.last_c02_neg    = last_c02_neg;
    summary.integ_set       = integ_set;
    summary.Cn_Ah           = Cn_Ah;                       
    summary.total_dQ_Ah_neg = total_dQ_Ah_discharge;      
    summary.total_dQ_Ah_net = total_dQ_Ah_net;           
    summary.total_dSOC      = total_dSOC;
    summary.SOC0_meta       = SOC0_meta;
    summary.SOC0_effective  = 100;                      
    summary.SOC_after       = SOC_after;

    
end
