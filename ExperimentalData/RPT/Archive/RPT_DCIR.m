%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT DCIR Processing (exact DCIR.m logic)
% - Load parsed MAT files from RPT_Parsing.m
% - Detect DCIR by pattern matching ('C','R' and 'D','R' sequences)
% - Compute SOC via current integration
% - Calculate R at 1/5/30s using exact time matching
% - Visualize per-channel Charge/Discharge DCIR vs SOC
% - mOhm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Paths and settings
parsedDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR';

if ~exist(saveDir,'dir'); mkdir(saveDir); end

channels = {'Ch9','Ch10','Ch11','Ch12','Ch13','Ch14'};%,'Ch15','Ch16'};
rpt_cycles = {'0cyc','200cyc','400cyc'};

% Minimum repeat count for DCIR detection (same as DCIR.m)
min_repeat_CR = 6;  % 'C','R' pattern
min_repeat_DR = 6;  % 'D','R' pattern

% Output struct
dcir_soc_data = struct();

% Debug flag
isDebug = true;

% Color mapping for cycles
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

    fig = figure('Name', sprintf('DCIR vs SOC - %s', channel), 'Position', [100 100 1200 800]);
    tlo = tiledlayout(fig, 2, 1, 'TileSpacing','compact','Padding','compact');

    % Charge subplot
    nexttile; hold on; grid on; title(sprintf('%s - Charge (DCIR)', channel));
    xlabel('SOC [%]'); ylabel('DCIR [m\Omega]');
    ylim([0 5]);

    % Discharge subplot
    nexttile; hold on; grid on; title(sprintf('%s - Discharge (DCIR)', channel));
    xlabel('SOC [%]'); ylabel('DCIR [m\Omega]');
    ylim([0 10]);

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
        fprintf('pdata length: %d\n', length(pdata));

        % Convert pdata to DCIR.m format (struct array with fields: t, V, I, type)
        data = [];
        for i = 1:length(pdata)
            data(i).t = pdata(i).t;
            data(i).V = pdata(i).V;
            data(i).I = pdata(i).I;
            data(i).type = char(pdata(i).type);  % Convert to char array
            data(i).steptime = pdata(i).steptime;
            data(i).steptime_double = pdata(i).steptime_double;
        end
        fprintf('Converted data length: %d\n', length(data));

        %% Charge DCIR (exact DCIR.m logic)
        % 'C','R'이 6번 이상 반복되는 구간을 Charge DCIR 이라 지정
        i = 1;
        data_len = length(data);
        all_types = {data.type};

        charge_start_index = NaN;
        charge_end_index = NaN;

        while i <= data_len - 1
            repeat_CR = 0;
            start_idx = i;

            % 'C','R'로 반복되는 구간을 탐색
            while i + 1 <= data_len
                if strcmp(all_types{i},'C') && strcmp(all_types{i+1},'R')
                    repeat_CR = repeat_CR + 1;
                    i = i + 2;
                else 
                    break
                end
            end
            % 'C','R'이 반복되는 구간의 시작 index 설정
            if repeat_CR >= min_repeat_CR
                charge_start_index = start_idx;
                charge_end_index = i - 1;  % 마지막 R까지 포함 
                break;  % 조건 만족 시 종료
            else
                i = start_idx + 1;  % 불충족 시 한 칸 뒤로 이동해서 재탐색
            end
        end

        % 구조체 추출
        if ~isnan(charge_start_index)
            data_DCIR_chg = data(charge_start_index : charge_end_index);
            fprintf('Charge DCIR found: %d steps\n', length(data_DCIR_chg));
        else
            fprintf('No Charge DCIR pattern found\n');
            data_DCIR_chg = [];
        end

        %% Discharge DCIR (exact DCIR.m logic)
        % 'D','R'이 6번 이상 반복되는 구간을 Discharge DCIR 이라 지정
        i = 1;
        data_len = length(data);
        all_types = {data.type};

        Discharge_start_index = NaN;
        Discharge_end_index = NaN;

        while i <= data_len - 1
            repeat_DR = 0;
            start_idx = i;
            % 'D','R'로 반복되는 구간을 탐색
            while i + 1 <= data_len
                if strcmp(all_types{i},'D') && strcmp(all_types{i+1},'R')
                    repeat_DR = repeat_DR + 1;
                    i = i + 2;
                else 
                    break
                end
            end
            % 'D','R'이 반복되는 구간의 시작 index 설정
            if repeat_DR >= min_repeat_DR
                Discharge_start_index = start_idx;
                Discharge_end_index = i - 1;  % 마지막 R까지 포함
                break;  % 조건 만족 시 종료
            else
                i = start_idx + 1;  % 불충족 시 한 칸 뒤로 이동해서 재탐색
            end
        end

        % 구조체 추출
        if ~isnan(Discharge_start_index)
            data_DCIR_dis = data(Discharge_start_index : Discharge_end_index);
            fprintf('Discharge DCIR found: %d steps\n', length(data_DCIR_dis));
        else
            fprintf('No Discharge DCIR pattern found\n');
            data_DCIR_dis = [];
        end

        %% Q 계산하기 (exact DCIR.m logic)
        % charge 
        % charge 상태에서 t, I 계산하기 위해 반출할 공간 생성
        t_chg_all = [];
        I_chg_all = [];

        for j = 1:length(data_DCIR_chg)
            t_tmp = data_DCIR_chg(j).t(:);
            I_tmp = data_DCIR_chg(j).I(:)*1000; %[mAh]

            t_chg_all = [t_chg_all; t_tmp];
            I_chg_all = [I_chg_all; I_tmp];
        end
        % cumQ, Q 계산
        cumQ_chg = cumtrapz(t_chg_all, I_chg_all) / 3600;
        Q_chg = trapz(t_chg_all, I_chg_all) / 3600;

        % discharge 상태
        t_dis_all = [];
        I_dis_all = [];

        % Discharge 상태에서 t, I 계산하기 위해 반출할 공간 생성
        for j = 1:length(data_DCIR_dis)
            t_dis_tmp = data_DCIR_dis(j).t(:);
            I_dis_tmp = abs(data_DCIR_dis(j).I(:))*1000; %[mAh]

            t_dis_all = [t_dis_all; t_dis_tmp];
            I_dis_all = [I_dis_all; I_dis_tmp];
        end

        % cumQ, Q 계산
        cumQ_dis = cumtrapz(t_dis_all, I_dis_all) / 3600;
        Q_dis = trapz(t_dis_all, I_dis_all) / 3600;

        % SOC 계산 (용량이 큰 것을 기준으로)
        if Q_dis > Q_chg
            Qref = Q_dis;
        else
            Qref = Q_chg;
        end

        % Charge SOC
        data_DCIR_chg_soc = cumQ_chg/Qref;

        % Discharge SOC
        data_DCIR_dis_soc = cumQ_dis/Qref;

        %% step 시간별로 SOC 잘라 넣기 (exact DCIR.m logic)
        % SOC 삽입
        % charge 
        % 계산한 SOC 를 step 별로 잘라 넣음
        start_idx = 1;
        for j = 1:length(data_DCIR_chg)
            n = length(data_DCIR_chg(j).t);
            idx_range = start_idx : start_idx + n - 1;

            data_DCIR_chg(j).soc = data_DCIR_chg_soc(idx_range);  % row로 맞춤
            start_idx = start_idx + n;
        end

        % discharge
        % 계산한 SOC 를 step 별로 잘라 넣음
        start_idx = 1;
           
        for j = 1:length(data_DCIR_dis)
            n = length(data_DCIR_dis(j).t);
            idx_range = start_idx : start_idx + n - 1;

            data_DCIR_dis(j).soc = data_DCIR_dis_soc(idx_range);  % row로 맞춤
            start_idx = start_idx + n;
        end

        %% DCIR 계산하기(charge) (exact DCIR.m logic)
        % DCIR 경과 시간 계산
        chg_data_DCIR_selected = data_DCIR_chg([]);  
        count = 0;

        % DCIR test C 만 추출 DCIR_time이 60초 지속인 것
        for i = 1:length(data_DCIR_chg)
            if isfield(data_DCIR_chg(i), 'steptime_double') && ~isempty(data_DCIR_chg(i).steptime_double)
                st = data_DCIR_chg(i).steptime_double;
                if isduration(st), st = seconds(st); end
                t_max = max(st);
                if t_max >= 59.9 && t_max <= 60
                    count = count + 1;
                    chg_data_DCIR_selected(count) = data_DCIR_chg(i);
                end
            end
        end
        % Debug: check selected pulse lengths
        fprintf('Charge DCIR selected pulses: %d\n', length(chg_data_DCIR_selected));
        for i = 1:length(chg_data_DCIR_selected)
            st = chg_data_DCIR_selected(i).steptime_double;
            if isduration(st), st = seconds(st); end
            t_max = max(st);
            fprintf('  Pulse %d: length=%d, t_max=%.1f\n', i, length(st), t_max);
        end
        
        % steptime already provided; ensure numeric type
        for i = 1:length(chg_data_DCIR_selected)
            st = chg_data_DCIR_selected(i).steptime_double;
            if isduration(st)
                chg_data_DCIR_selected(i).steptime_double = seconds(st);
            end
        end

        
        target_times = [1, 5, 30, 60];  % 사용자 지정 시점
        s = chg_data_DCIR_selected;

        % target time에 따라 추출하기
        for i = 1:length(s)
            t_vec = s(i).steptime_double;  % 이미 double
            V_vec = s(i).V;
            soc_vec = s(i).soc;
            I_1C_chg = mean(s(i).I);


            idx_list = zeros(size(target_times));

            % 시계열 해상도로 안전하게 인덱싱(스텝 상대시간이 0부터 시작한다고 가정)
            t_num = t_vec;
            if isduration(t_num)
                t_num = seconds(t_num);
            end
            for k = 1:length(target_times)
                idx_k = find(t_num == target_times(k), 1, 'first');
                idx_list(k) = idx_k; % will error if empty (keeps original logic)
            end

            % 인덱스로 추출
            chg_data_DCIR_selected(i).Vstart = V_vec(1);  % 첫 번째 시점 전압 (펄스 시작 전 전압)
            chg_data_DCIR_selected(i).V1s = V_vec(idx_list(1));   % 1초 시점 전압
            chg_data_DCIR_selected(i).V5s = V_vec(idx_list(2));   % 5초 시점 전압
            chg_data_DCIR_selected(i).V30s = V_vec(idx_list(3));  % 30초 시점 전압
            chg_data_DCIR_selected(i).V60s = V_vec(idx_list(4));  % 60초 시점 전압
            chg_data_DCIR_selected(i).SOC = soc_vec(1)*100;

            % (요청된 로그만 남김)
        end

        % 계산하기
        for i = 1:length(chg_data_DCIR_selected)
            Vstart = chg_data_DCIR_selected(i).Vstart;
            V1s = chg_data_DCIR_selected(i).V1s;
            V5s = chg_data_DCIR_selected(i).V5s;
            V30s = chg_data_DCIR_selected(i).V30s;
            V60s = chg_data_DCIR_selected(i).V60s;
              
            % ΔV (구간별 전압 변화) → 요구사항: 1s-0s, 5s-1s, 60s-5s
            dV_1s_0s  = V1s - Vstart;   % 1s - 0s
            dV_5s_1s  = V5s - V1s;      % 5s - 1s
            dV_60s_5s = V60s - V5s;     % 60s - 5s
            
            % 구간별 저항 (증분 저항) [mΩ]
            chg_data_DCIR_selected(i).dv = [dV_1s_0s, dV_5s_1s, dV_60s_5s];
            chg_data_DCIR_selected(i).dR = 1000*[dV_1s_0s/I_1C_chg, dV_5s_1s/I_1C_chg, dV_60s_5s/I_1C_chg];
            
            % 시점별 저항 (해당 시점에서의 저항값) [mΩ]
            chg_data_DCIR_selected(i).R = 1000*[(V1s-Vstart)/I_1C_chg, (V5s-Vstart)/I_1C_chg, (V30s-Vstart)/I_1C_chg, (V60s-Vstart)/I_1C_chg];
        end

        %% DCIR 계산하기(Discharge) (exact DCIR.m logic)
        % DCIR 경과 시간 계산
        dis_data_DCIR_selected = data_DCIR_dis([]);  
        count = 0;

        % DCIR test D만 추출 DCIR_time이 60초 지속인 것
        for i = 1:length(data_DCIR_dis)
            if isfield(data_DCIR_dis(i), 'steptime_double') && ~isempty(data_DCIR_dis(i).steptime_double)
                st = data_DCIR_dis(i).steptime_double;
                if isduration(st), st = seconds(st); end
                t_max = max(st);
                if t_max >= 59.9 && t_max <= 60
                    count = count + 1;
                    dis_data_DCIR_selected(count) = data_DCIR_dis(i);
                end
            end
        end
        % Debug: check selected pulse lengths
        fprintf('Discharge DCIR selected pulses: %d\n', length(dis_data_DCIR_selected));
        for i = 1:length(dis_data_DCIR_selected)
            st = dis_data_DCIR_selected(i).steptime_double;
            if isduration(st), st = seconds(st); end
            t_max = max(st);
            fprintf('  Pulse %d: length=%d, t_max=%.1f\n', i, length(st), t_max);
        end
        
        for i = 1:length(dis_data_DCIR_selected)
            st = dis_data_DCIR_selected(i).steptime_double;
            if isduration(st)
                dis_data_DCIR_selected(i).steptime_double = seconds(st);
            end
        end

        %% R 계산하기 (Discharge) (exact DCIR.m logic)
        target_times = [1, 5, 30, 60];  % 사용자 지정 시점
        s = dis_data_DCIR_selected;

        % target time에 따라 추출하기
        for i = 1:length(s)
            t_vec = s(i).steptime_double;  % 이미 double
            V_vec = s(i).V;
            soc_vec = s(i).soc;
            I_1C_dis = mean(s(i).I);

            idx_list = zeros(size(target_times));

            t_num = t_vec;
            if isduration(t_num)
                t_num = seconds(t_num);
            end
            for k = 1:length(target_times)
                idx_k = find(t_num == target_times(k), 1, 'first');
                idx_list(k) = idx_k; % will error if empty (keeps original logic)
            end

            % 인덱스로 추출
            dis_data_DCIR_selected(i).Vstart = V_vec(1);  % 첫 번째 시점 전압 (펄스 시작 전 전압)
            dis_data_DCIR_selected(i).V1s = V_vec(idx_list(1));   % 1초 시점 전압
            dis_data_DCIR_selected(i).V5s = V_vec(idx_list(2));   % 5초 시점 전압
            dis_data_DCIR_selected(i).V30s = V_vec(idx_list(3));  % 30초 시점 전압
            dis_data_DCIR_selected(i).V60s = V_vec(idx_list(4));  % 60초 시점 전압
            dis_data_DCIR_selected(i).SOC = soc_vec(1)*100;

            % (요청된 로그만 남김)
        end

        % 계산하기
        for i = 1:length(dis_data_DCIR_selected)
            Vstart = dis_data_DCIR_selected(i).Vstart;
            V1s = dis_data_DCIR_selected(i).V1s;
            V5s = dis_data_DCIR_selected(i).V5s;
            V30s = dis_data_DCIR_selected(i).V30s;
            V60s = dis_data_DCIR_selected(i).V60s;

            % ΔV (구간별 전압 변화) → 1s-0s, 5s-1s, 60s-5s
            dV1 = V1s - Vstart; % 1s - 0s
            dV2 = V5s - V1s;    % 5s - 1s
            dV3 = V60s - V5s;   % 60s - 5s

            % 구간별 저항 (증분 저항) [mΩ]
            dis_data_DCIR_selected(i).dv = [dV1, dV2, dV3];
            dis_data_DCIR_selected(i).dR = 1000*[dV1/I_1C_dis, dV2/I_1C_dis, dV3/I_1C_dis];
            
            % 시점별 저항 (해당 시점에서의 저항값) [mΩ]
            dis_data_DCIR_selected(i).R = 1000*[(V1s-Vstart)/I_1C_dis, (V5s-Vstart)/I_1C_dis, (V30s-Vstart)/I_1C_dis, (V60s-Vstart)/I_1C_dis];
        end

        %% Store and plot data
        % Charge data
        if ~isempty(chg_data_DCIR_selected)
            n_chg = length(chg_data_DCIR_selected);
            SOC_chg = zeros(n_chg, 1);
            R_chg = zeros(n_chg, 4);      % 시점별 저항 (1s,5s,30s,60s)
            dR_chg = zeros(n_chg, 3);     % 구간별 저항 (1→5, 5→30, 30→60)

            for i = 1:n_chg
                SOC_chg(i) = chg_data_DCIR_selected(i).SOC;  % 실제 SOC 값
                R_chg(i,:) = chg_data_DCIR_selected(i).R;    % 시점별 저항
                dR_chg(i,:) = chg_data_DCIR_selected(i).dR;  % 구간별 저항
            end
            
            % cumQ, Q 데이터 추출 (레퍼런스와 동일 로직: 전체 DCIR 구간)
            cumQ_chg_array = zeros(n_chg, 1);
            Q_chg_array = zeros(n_chg, 1);
            
            % steptime_double 사용 (상대시간, 0부터 시작)
            t_chg_all = [];
            I_chg_all = [];
            for i = 1:n_chg
                t_step = chg_data_DCIR_selected(i).steptime_double(:);
                I_step = chg_data_DCIR_selected(i).I(:);
                
                % 60초까지만 사용 (max 값의 인덱스)
                max_idx = find(t_step == max(t_step), 1, 'last');
                t_step = t_step(1:max_idx);
                I_step = I_step(1:max_idx);
                
                % 시간을 연속적으로 만들기 위해 이전 펄스의 마지막 시간에 더함
                if ~isempty(t_chg_all)
                    t_step = t_step + t_chg_all(end);
                end
                
                t_chg_all = [t_chg_all; t_step]; %#ok<AGROW>
                I_chg_all = [I_chg_all; I_step]; %#ok<AGROW>
            end
            
            % 레퍼런스와 동일: 전체 구간에서 cumQ, Q 계산
            if ~isempty(t_chg_all)
                cumQ_curve = cumtrapz(t_chg_all, I_chg_all) / 3600; % Ah
                Q_chg_total = trapz(t_chg_all, I_chg_all) / 3600; % Ah
                
                % 각 펄스의 끝 시점 cumQ 값 사용 (누적 용량)
                start_idx = 1;
                for i = 1:n_chg
                    t_step = chg_data_DCIR_selected(i).steptime_double(:);
                    max_idx = find(t_step == max(t_step), 1, 'last');
                    n = max_idx;
                    end_idx = start_idx + n - 1;
                    cumQ_chg_array(i) = cumQ_curve(end_idx);
                    Q_chg_array(i) = Q_chg_total; % 모든 펄스에서 동일한 총 용량
                    start_idx = start_idx + n;
                end
            else
                Q_chg_total = 0;
            end
            % (요청된 로그만 남김)

            % Store charge data - 시점별 저항
            dcir_soc_data.(channel).(rpt_key).charge_soc = SOC_chg;
            dcir_soc_data.(channel).(rpt_key).charge_cumQ = cumQ_chg_array;
            dcir_soc_data.(channel).(rpt_key).charge_Q = Q_chg_array;
            dcir_soc_data.(channel).(rpt_key).charge_r_1s = R_chg(:,1);
            dcir_soc_data.(channel).(rpt_key).charge_r_5s = R_chg(:,2);
            dcir_soc_data.(channel).(rpt_key).charge_r_30s = R_chg(:,3);
            dcir_soc_data.(channel).(rpt_key).charge_r_60s = R_chg(:,4);
            
            % Store charge data - 구간별 저항 (1s-0s, 5s-1s, 60s-5s)
            dcir_soc_data.(channel).(rpt_key).charge_dr_1s0s = dR_chg(:,1);
            dcir_soc_data.(channel).(rpt_key).charge_dr_5s1s = dR_chg(:,2);
            dcir_soc_data.(channel).(rpt_key).charge_dr_60s5s = dR_chg(:,3);
            
            % Store charge summary table
            dcir_soc_data.(channel).(rpt_key).charge_summary_table = table(...
                SOC_chg, cumQ_chg_array, Q_chg_array, R_chg(:,1), R_chg(:,2), R_chg(:,3), R_chg(:,4), ...
                'VariableNames', {'SOC_pct', 'cumQ_Ah', 'Q_Ah', 'R_1s_mOhm', 'R_5s_mOhm', 'R_30s_mOhm', 'R_60s_mOhm'});

            % Plot charge data (60s 시점 저항)
            nexttile(1);
            color_val = getCycleColor(rpt_cycle);
            scatter(SOC_chg, R_chg(:,4), 50, color_val, 'filled', 'DisplayName', sprintf('%s 60s', rpt_cycle));
        end

        % Discharge data
        if ~isempty(dis_data_DCIR_selected)
            n_dis = length(dis_data_DCIR_selected);
            SOC_dis = zeros(n_dis, 1);
            R_dis = zeros(n_dis, 4);      % 시점별 저항 (1s,5s,30s,60s)
            dR_dis = zeros(n_dis, 3);     % 구간별 저항 (1→5, 5→30, 30→60)

            for i = 1:n_dis
                SOC_dis(i) = dis_data_DCIR_selected(i).SOC;  % 실제 SOC 값
                R_dis(i,:) = dis_data_DCIR_selected(i).R;    % 시점별 저항
                dR_dis(i,:) = dis_data_DCIR_selected(i).dR;  % 구간별 저항
            end
            
            % cumQ, Q 데이터 추출 (레퍼런스와 동일 로직: 전체 DCIR 구간)
            cumQ_dis_array = zeros(n_dis, 1);
            Q_dis_array = zeros(n_dis, 1);
            
            % steptime_double 사용 (상대시간, 0부터 시작)
            t_dis_all = [];
            I_dis_all = [];
            for i = 1:n_dis
                t_step = dis_data_DCIR_selected(i).steptime_double(:);
                I_step = abs(dis_data_DCIR_selected(i).I(:));
                
                % 60초까지만 사용 (max 값의 인덱스)
                max_idx = find(t_step == max(t_step), 1, 'last');
                t_step = t_step(1:max_idx);
                I_step = I_step(1:max_idx);
                
                % 시간을 연속적으로 만들기 위해 이전 펄스의 마지막 시간에 더함
                if ~isempty(t_dis_all)
                    t_step = t_step + t_dis_all(end);
                end
                
                t_dis_all = [t_dis_all; t_step]; %#ok<AGROW>
                I_dis_all = [I_dis_all; I_step]; %#ok<AGROW>
            end
            
            % 레퍼런스와 동일: 전체 구간에서 cumQ, Q 계산
            if ~isempty(t_dis_all)
                cumQ_curve = cumtrapz(t_dis_all, I_dis_all) / 3600; % Ah
                Q_dis_total = trapz(t_dis_all, I_dis_all) / 3600; % Ah
                
                % 각 펄스의 끝 시점 cumQ 값 사용 (누적 용량)
                start_idx = 1;
                for i = 1:n_dis
                    t_step = dis_data_DCIR_selected(i).steptime_double(:);
                    max_idx = find(t_step == max(t_step), 1, 'last');
                    n = max_idx;
                    end_idx = start_idx + n - 1;
                    cumQ_dis_array(i) = cumQ_curve(end_idx);
                    Q_dis_array(i) = Q_dis_total; % 모든 펄스에서 동일한 총 용량
                    start_idx = start_idx + n;
                end
            else
                Q_dis_total = 0;
            end
            
            % 레퍼런스와 동일: Qref = max(Q_chg, Q_dis)로 SOC 계산
            if Q_dis_total > Q_chg_total
                Qref = Q_dis_total;
            else
                Qref = Q_chg_total;
            end
            
            % SOC 재계산 (공통 Qref 사용)
            if Qref > 0
                SOC_chg = (cumQ_chg_array / Qref) * 100; % %
                SOC_dis = (cumQ_dis_array / Qref) * 100; % %
            else
                SOC_chg = zeros(size(cumQ_chg_array));
                SOC_dis = zeros(size(cumQ_dis_array));
            end
            % (요청된 로그만 남김)

            % Store discharge data - 시점별 저항
            dcir_soc_data.(channel).(rpt_key).discharge_soc = SOC_dis;
            dcir_soc_data.(channel).(rpt_key).discharge_cumQ = cumQ_dis_array;
            dcir_soc_data.(channel).(rpt_key).discharge_Q = Q_dis_array;
            dcir_soc_data.(channel).(rpt_key).discharge_r_1s = R_dis(:,1);
            dcir_soc_data.(channel).(rpt_key).discharge_r_5s = R_dis(:,2);
            dcir_soc_data.(channel).(rpt_key).discharge_r_30s = R_dis(:,3);
            dcir_soc_data.(channel).(rpt_key).discharge_r_60s = R_dis(:,4);
            
            % Store discharge data - 구간별 저항 (1s-0s, 5s-1s, 60s-5s)
            dcir_soc_data.(channel).(rpt_key).discharge_dr_1s0s = dR_dis(:,1);
            dcir_soc_data.(channel).(rpt_key).discharge_dr_5s1s = dR_dis(:,2);
            dcir_soc_data.(channel).(rpt_key).discharge_dr_60s5s = dR_dis(:,3);
            
            % Store discharge summary table
            dcir_soc_data.(channel).(rpt_key).discharge_summary_table = table(...
                SOC_dis, cumQ_dis_array, Q_dis_array, R_dis(:,1), R_dis(:,2), R_dis(:,3), R_dis(:,4), ...
                'VariableNames', {'SOC_pct', 'cumQ_Ah', 'Q_Ah', 'R_1s_mOhm', 'R_5s_mOhm', 'R_30s_mOhm', 'R_60s_mOhm'});

            % Plot discharge data (60s 시점 저항)
            nexttile(2);
            color_val = getCycleColor(rpt_cycle);
            scatter(SOC_dis, R_dis(:,4), 50, color_val, 'filled', 'DisplayName', sprintf('%s 60s', rpt_cycle));
        end
    end

    % legends
    nexttile(1); legend('Location','north');
    nexttile(2); legend('Location','north');

    % Save integrated figure (all cycles)
    savefig(fig, fullfile(saveDir, sprintf('DCIR_%s_integrated.fig', channel)));
    
    % Save individual cycle figures
    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        cyc_num = strrep(rpt_cycle, 'cyc', '');
        
        % Create individual figure for each cycle (histogram-style of point-in-time R)
        fig_individual = figure('Name', sprintf('DCIR vs SOC - %s %s', channel, rpt_cycle), 'Position', [100 100 1200 800]);
        tlo_individual = tiledlayout(fig_individual, 2, 1, 'TileSpacing','compact','Padding','compact');
        
        % Charge subplot
        nexttile; hold on; grid on; title(sprintf('%s - Charge (DCIR) - %s', channel, rpt_cycle));
        xlabel('SOC [%]'); ylabel('DCIR [m\Omega]');
        ylim([0 5]);
        
        % Discharge subplot
        nexttile; hold on; grid on; title(sprintf('%s - Discharge (DCIR) - %s', channel, rpt_cycle));
        xlabel('SOC [%]'); ylabel('DCIR [m\Omega]');
        ylim([0 10]);
        
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
        
        % Plot charge data for this cycle as stacked bar of interval dR (mΩ)
        if isfield(dcir_soc_data, channel) && isfield(dcir_soc_data.(channel), rpt_key)
            if isfield(dcir_soc_data.(channel).(rpt_key), 'charge_soc')
                SOC_chg = dcir_soc_data.(channel).(rpt_key).charge_soc;
                DR1 = dcir_soc_data.(channel).(rpt_key).charge_dr_1s0s(:);   % mΩ
                DR2 = dcir_soc_data.(channel).(rpt_key).charge_dr_5s1s(:);   % mΩ
                DR3 = dcir_soc_data.(channel).(rpt_key).charge_dr_60s5s(:);  % mΩ
                nexttile(1);
                x = (1:numel(SOC_chg)).';
                Rstack = [DR1, DR2, DR3];
                bar(x, Rstack, 'stacked');
                xticks(x);
                xticklabels(arrayfun(@(v) sprintf('%.2f', v), SOC_chg, 'UniformOutput', false));
                xtickangle(45);
                legend({'1s-0s','5s-1s','60s-5s'}, 'Location','best');
            end
        end
        
        % Plot discharge data for this cycle as stacked bar of interval dR (mΩ)
        if isfield(dcir_soc_data, channel) && isfield(dcir_soc_data.(channel), rpt_key)
            if isfield(dcir_soc_data.(channel).(rpt_key), 'discharge_soc')
                SOC_dis = dcir_soc_data.(channel).(rpt_key).discharge_soc;
                DR1 = dcir_soc_data.(channel).(rpt_key).discharge_dr_1s0s(:);   % mΩ
                DR2 = dcir_soc_data.(channel).(rpt_key).discharge_dr_5s1s(:);   % mΩ
                DR3 = dcir_soc_data.(channel).(rpt_key).discharge_dr_60s5s(:);  % mΩ
                nexttile(2);
                x = (1:numel(SOC_dis)).';
                Rstack = [DR1, DR2, DR3];
                bar(x, Rstack, 'stacked');
                xticks(x);
                xticklabels(arrayfun(@(v) sprintf('%.2f', v), SOC_dis, 'UniformOutput', false));
                xtickangle(45);
                legend({'1s-0s','5s-1s','60s-5s'}, 'Location','best');
            end
        end
        
        % Add legends
        nexttile(1); legend('Location','best');
        nexttile(2); legend('Location','best');
        
        % Save individual cycle figure
        savefig(fig_individual, fullfile(saveDir, sprintf('DCIR_%s_%scyc.fig', channel, cyc_num)));
        close(fig_individual);
    end
end

% Create summary table for each cycle
fprintf('\n=== DCIR Summary Table ===\n');
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    fprintf('\n--- %s ---\n', channel);
    
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
        
        if isfield(dcir_soc_data, channel) && isfield(dcir_soc_data.(channel), rpt_key)
            fprintf('\n%s:\n', rpt_cycle);
            
            % Charge data
            if isfield(dcir_soc_data.(channel).(rpt_key), 'charge_soc')
                SOC_chg = dcir_soc_data.(channel).(rpt_key).charge_soc;
                R_chg_1s = dcir_soc_data.(channel).(rpt_key).charge_r_1s;
                R_chg_5s = dcir_soc_data.(channel).(rpt_key).charge_r_5s;
                R_chg_30s = dcir_soc_data.(channel).(rpt_key).charge_r_30s;
                R_chg_60s = dcir_soc_data.(channel).(rpt_key).charge_r_60s;
                
                cumQ_chg = dcir_soc_data.(channel).(rpt_key).charge_cumQ;
                Q_chg = dcir_soc_data.(channel).(rpt_key).charge_Q;
                
                fprintf('  Charge DCIR:\n');
                fprintf('    SOC[%%]    cumQ[Ah]    Q[Ah]    R_1s[mΩ]    R_5s[mΩ]    R_30s[mΩ]    R_60s[mΩ]\n');
                for i = 1:length(SOC_chg)
                    fprintf('    %6.2f    %8.3f    %6.3f    %8.4f    %8.4f    %8.4f    %8.4f\n', ...
                        SOC_chg(i), cumQ_chg(i), Q_chg(i), R_chg_1s(i), R_chg_5s(i), R_chg_30s(i), R_chg_60s(i));
                end
            end
            
            % Discharge data
            if isfield(dcir_soc_data.(channel).(rpt_key), 'discharge_soc')
                SOC_dis = dcir_soc_data.(channel).(rpt_key).discharge_soc;
                R_dis_1s = dcir_soc_data.(channel).(rpt_key).discharge_r_1s;
                R_dis_5s = dcir_soc_data.(channel).(rpt_key).discharge_r_5s;
                R_dis_30s = dcir_soc_data.(channel).(rpt_key).discharge_r_30s;
                R_dis_60s = dcir_soc_data.(channel).(rpt_key).discharge_r_60s;
                
                cumQ_dis = dcir_soc_data.(channel).(rpt_key).discharge_cumQ;
                Q_dis = dcir_soc_data.(channel).(rpt_key).discharge_Q;
                
                fprintf('  Discharge DCIR:\n');
                fprintf('    SOC[%%]    cumQ[Ah]    Q[Ah]    R_1s[mΩ]    R_5s[mΩ]    R_30s[mΩ]    R_60s[mΩ]\n');
                for i = 1:length(SOC_dis)
                    fprintf('    %6.2f    %8.3f    %6.3f    %8.4f    %8.4f    %8.4f    %8.4f\n', ...
                        SOC_dis(i), cumQ_dis(i), Q_dis(i), R_dis_1s(i), R_dis_5s(i), R_dis_30s(i), R_dis_60s(i));
                end
            end
        end
    end
end

% Save MAT with structured data
matFile = fullfile(saveDir, 'DCIR_SOC_data_9to14.mat');
save(matFile, 'dcir_soc_data');
fprintf('\nSaved: %s\n', matFile);

% Run RPT_DCIR_1s.m for 1s-sampled channels
RPT_DCIR_1s

