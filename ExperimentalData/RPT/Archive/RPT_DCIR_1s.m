%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT DCIR Processing for 1s-sampled channels (Ch15, Ch16 only)
% - Load parsed MAT files from RPT_Parsing.m
% - Detect DCIR by pattern matching ('C','R' and 'D','R' sequences)
% - Compute SOC via current integration (exact DCIR.m method)
% - Calculate point-in-time R at 1/5/30/60s and interval dR (1-0,5-1,60-5)
% - Robust to 1 s sampling via nearest-time indexing on step-relative time
% - Save integrated/per-cycle figures and MAT per channel
% - All resistance in mOhm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Paths and settings
parsedDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR';

if ~exist(saveDir,'dir'); mkdir(saveDir); end

channels = {'Ch15','Ch16'};  % Only 1 s-sampled channels
rpt_cycles = {'0cyc','200cyc','400cyc'};

% Minimum repeat count for DCIR detection (same as DCIR.m)
min_repeat_CR = 6;  % 'C','R' pattern
min_repeat_DR = 6;  % 'D','R' pattern

% Output struct
dcir_soc_data = struct();

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

        %% Charge DCIR (exact DCIR.m logic)
        i = 1;
        data_len = length(data);
        all_types = {data.type};

        charge_start_index = NaN;
        charge_end_index = NaN;

        while i <= data_len - 1
            repeat_CR = 0;
            start_idx = i;
            while i + 1 <= data_len
                if strcmp(all_types{i},'C') && strcmp(all_types{i+1},'R')
                    repeat_CR = repeat_CR + 1;
                    i = i + 2;
                else 
                    break
                end
            end
            if repeat_CR >= min_repeat_CR
                charge_start_index = start_idx;
                charge_end_index = i - 1;  % include last R
                break;
            else
                i = start_idx + 1;
            end
        end

        if ~isnan(charge_start_index)
            data_DCIR_chg = data(charge_start_index : charge_end_index);
        else
            data_DCIR_chg = [];
        end

        %% Discharge DCIR (exact DCIR.m logic)
        i = 1;
        data_len = length(data);
        all_types = {data.type};

        Discharge_start_index = NaN;
        Discharge_end_index = NaN;

        while i <= data_len - 1
            repeat_DR = 0;
            start_idx = i;
            while i + 1 <= data_len
                if strcmp(all_types{i},'D') && strcmp(all_types{i+1},'R')
                    repeat_DR = repeat_DR + 1;
                    i = i + 2;
                else 
                    break
                end
            end
            if repeat_DR >= min_repeat_DR
                Discharge_start_index = start_idx;
                Discharge_end_index = i - 1;  % include last R
                break;
            else
                i = start_idx + 1;
            end
        end

        if ~isnan(Discharge_start_index)
            data_DCIR_dis = data(Discharge_start_index : Discharge_end_index);
        else
            data_DCIR_dis = [];
        end

        %% Q 계산하기 (exact DCIR.m logic)
        % charge
        t_chg_all = [];
        I_chg_all = [];
        for j = 1:length(data_DCIR_chg)
            t_tmp = data_DCIR_chg(j).t(:);
            I_tmp = data_DCIR_chg(j).I(:)*1000; %[mAh]
            t_chg_all = [t_chg_all; t_tmp];
            I_chg_all = [I_chg_all; I_tmp];
        end
        cumQ_chg = cumtrapz(t_chg_all, I_chg_all) / 3600;
        Q_chg = trapz(t_chg_all, I_chg_all) / 3600;

        % discharge
        t_dis_all = [];
        I_dis_all = [];
        for j = 1:length(data_DCIR_dis)
            t_dis_tmp = data_DCIR_dis(j).t(:);
            I_dis_tmp = abs(data_DCIR_dis(j).I(:))*1000; %[mAh]
            t_dis_all = [t_dis_all; t_dis_tmp];
            I_dis_all = [I_dis_all; I_dis_tmp];
        end
        cumQ_dis = cumtrapz(t_dis_all, I_dis_all) / 3600;
        Q_dis = trapz(t_dis_all, I_dis_all) / 3600;

        % SOC 기준 (큰 용량)
        if Q_dis > Q_chg
            Qref = Q_dis;
        else
            Qref = Q_chg;
        end
        data_DCIR_chg_soc = cumQ_chg/Qref;
        data_DCIR_dis_soc = cumQ_dis/Qref;

        %% step 시간별로 SOC 잘라 넣기 (exact DCIR.m logic)
        start_idx = 1;
        for j = 1:length(data_DCIR_chg)
            n = length(data_DCIR_chg(j).t);
            idx_range = start_idx : start_idx + n - 1;
            data_DCIR_chg(j).soc = data_DCIR_chg_soc(idx_range);
            start_idx = start_idx + n;
        end

        start_idx = 1;
        for j = 1:length(data_DCIR_dis)
            n = length(data_DCIR_dis(j).t);
            idx_range = start_idx : start_idx + n - 1;
            data_DCIR_dis(j).soc = data_DCIR_dis_soc(idx_range);
            start_idx = start_idx + n;
        end

        %% DCIR 계산하기 (Charge)
        chg_data_DCIR_selected = data_DCIR_chg([]);
        count = 0;
        for i = 1:length(data_DCIR_chg)
            if isfield(data_DCIR_chg(i), 't') && ~isempty(data_DCIR_chg(i).t)
                t_vec = data_DCIR_chg(i).t;
                t_elapsed = t_vec(end) - t_vec(1);
                if t_elapsed >= 59.9 && t_elapsed <= 60
                    count = count + 1;
                    chg_data_DCIR_selected(count) = data_DCIR_chg(i);
                end
            end
        end
        % steptime already provided; ensure numeric type
        for i = 1:length(chg_data_DCIR_selected)
            st = chg_data_DCIR_selected(i).steptime_double;
            if isduration(st)
                chg_data_DCIR_selected(i).steptime_double = seconds(st);
            end
        end

        target_times = [2, 5, 30, 60];
        s = chg_data_DCIR_selected;
        for i = 1:length(s)
            t_vec = s(i).steptime_double;
            V_vec = s(i).V;
            soc_vec = s(i).soc;
            I_mean = mean(s(i).I);
            idx_list = zeros(size(target_times));
            t_num = t_vec;
            if isduration(t_num)
                t_num = seconds(t_num);
            end
            for k = 1:length(target_times)
                idx_k = find(t_num == target_times(k), 1, 'first');
                idx_list(k) = idx_k; % will error if empty (keeps original logic)
            end
            chg_data_DCIR_selected(i).Vstart = V_vec(1);
            chg_data_DCIR_selected(i).V1s = V_vec(idx_list(1));
            chg_data_DCIR_selected(i).V5s = V_vec(idx_list(2));
            chg_data_DCIR_selected(i).V30s = V_vec(idx_list(3));
            chg_data_DCIR_selected(i).V60s = V_vec(idx_list(4));
            chg_data_DCIR_selected(i).SOC = soc_vec(1)*100;
            % dV and R/dR in mOhm
            dV_1s_0s  = chg_data_DCIR_selected(i).V1s - chg_data_DCIR_selected(i).Vstart;
            dV_5s_1s  = chg_data_DCIR_selected(i).V5s - chg_data_DCIR_selected(i).V1s;
            dV_60s_5s = chg_data_DCIR_selected(i).V60s - chg_data_DCIR_selected(i).V5s;
            chg_data_DCIR_selected(i).dv = [dV_1s_0s, dV_5s_1s, dV_60s_5s];
            chg_data_DCIR_selected(i).dR = 1000*[dV_1s_0s/I_mean, dV_5s_1s/I_mean, dV_60s_5s/I_mean];
            chg_data_DCIR_selected(i).R  = 1000*[(chg_data_DCIR_selected(i).V1s-chg_data_DCIR_selected(i).Vstart)/I_mean,
                                                 (chg_data_DCIR_selected(i).V5s-chg_data_DCIR_selected(i).Vstart)/I_mean,
                                                 (chg_data_DCIR_selected(i).V30s-chg_data_DCIR_selected(i).Vstart)/I_mean,
                                                 (chg_data_DCIR_selected(i).V60s-chg_data_DCIR_selected(i).Vstart)/I_mean];
        end

        %% DCIR 계산하기 (Discharge)
        dis_data_DCIR_selected = data_DCIR_dis([]);
        count = 0;
        for i = 1:length(data_DCIR_dis)
            if isfield(data_DCIR_dis(i), 't') && ~isempty(data_DCIR_dis(i).t)
                t_vec = data_DCIR_dis(i).t;
                t_elapsed = t_vec(end) - t_vec(1);
                if t_elapsed >= 59.9 && t_elapsed <= 60
                    count = count + 1;
                    dis_data_DCIR_selected(count) = data_DCIR_dis(i);
                end
            end
        end
        for i = 1:length(dis_data_DCIR_selected)
            st = dis_data_DCIR_selected(i).steptime_double;
            if isduration(st)
                dis_data_DCIR_selected(i).steptime_double = seconds(st);
            end
        end

        target_times = [2, 5, 30, 60];
        s = dis_data_DCIR_selected;
        for i = 1:length(s)
            t_vec = s(i).steptime_double;
            V_vec = s(i).V;
            soc_vec = s(i).soc;
            I_mean = mean(s(i).I);
            idx_list = zeros(size(target_times));
            t_num = t_vec;
            if isduration(t_num)
                t_num = seconds(t_num);
            end
            for k = 1:length(target_times)
                idx_k = find(t_num == target_times(k), 1, 'first');
                idx_list(k) = idx_k; % will error if empty (keeps original logic)
            end
            dis_data_DCIR_selected(i).Vstart = V_vec(1);
            dis_data_DCIR_selected(i).V1s = V_vec(idx_list(1));
            dis_data_DCIR_selected(i).V5s = V_vec(idx_list(2));
            dis_data_DCIR_selected(i).V30s = V_vec(idx_list(3));
            dis_data_DCIR_selected(i).V60s = V_vec(idx_list(4));
            dis_data_DCIR_selected(i).SOC = soc_vec(1)*100;
            dV1 = dis_data_DCIR_selected(i).V1s - dis_data_DCIR_selected(i).Vstart;
            dV2 = dis_data_DCIR_selected(i).V5s - dis_data_DCIR_selected(i).V1s;
            dV3 = dis_data_DCIR_selected(i).V60s - dis_data_DCIR_selected(i).V5s;
            dis_data_DCIR_selected(i).dv = [dV1, dV2, dV3];
            dis_data_DCIR_selected(i).dR = 1000*[dV1/I_mean, dV2/I_mean, dV3/I_mean];
            dis_data_DCIR_selected(i).R  = 1000*[(dis_data_DCIR_selected(i).V1s-dis_data_DCIR_selected(i).Vstart)/I_mean,
                                                 (dis_data_DCIR_selected(i).V5s-dis_data_DCIR_selected(i).Vstart)/I_mean,
                                                 (dis_data_DCIR_selected(i).V30s-dis_data_DCIR_selected(i).Vstart)/I_mean,
                                                 (dis_data_DCIR_selected(i).V60s-dis_data_DCIR_selected(i).Vstart)/I_mean];
        end

        %% Store and plot data
        % Charge data
        if ~isempty(chg_data_DCIR_selected)
            n_chg = length(chg_data_DCIR_selected);
            SOC_chg = zeros(n_chg, 1);
            R_chg = zeros(n_chg, 4);
            dR_chg = zeros(n_chg, 3);
            for i = 1:n_chg
                SOC_chg(i) = chg_data_DCIR_selected(i).SOC;
                R_chg(i,:) = chg_data_DCIR_selected(i).R;
                dR_chg(i,:) = chg_data_DCIR_selected(i).dR;
            end
            % cumQ, Q (concatenate steps)
            cumQ_chg_array = zeros(n_chg, 1);
            Q_chg_array = zeros(n_chg, 1);
            t_all = [];
            I_all_A = [];
            step_lengths = zeros(n_chg,1);
            for i = 1:n_chg
                t_rel = chg_data_DCIR_selected(i).steptime_double(:);
                % 60초까지만 사용 (max 값의 인덱스)
                max_idx = find(t_rel == max(t_rel), 1, 'last');
                step_lengths(i) = max_idx;
            end
            abs_offset = 0;
            for i = 1:n_chg
                t_rel = chg_data_DCIR_selected(i).steptime_double(:);
                if isempty(t_rel), t_rel = 0; end
                
                % 60초까지만 사용 (max 값의 인덱스)
                max_idx = find(t_rel == max(t_rel), 1, 'last');
                t_rel = t_rel(1:max_idx);
                I_step_A = chg_data_DCIR_selected(i).I(1:max_idx);
                
                t_abs = abs_offset + t_rel;
                t_all = [t_all; t_abs]; %#ok<AGROW>
                I_all_A = [I_all_A; I_step_A]; %#ok<AGROW>
                abs_offset = t_abs(end);
            end
            if isduration(t_all), t_all = seconds(t_all); end
            cumQ_curve = cumtrapz(t_all, I_all_A) / 3600;
            Q_chg_total = trapz(t_all, I_all_A) / 3600; % Ah
            ends = cumsum(step_lengths);
            prev_end = 0;
            for i = 1:n_chg
                cumQ_chg_array(i) = cumQ_curve(ends(i));
                Q_chg_array(i) = Q_chg_total; % 모든 펄스에서 동일한 총 용량
                prev_end = ends(i);
            end
            % store
            dcir_soc_data.(channel).(rpt_key).charge_soc = SOC_chg;
            dcir_soc_data.(channel).(rpt_key).charge_cumQ = cumQ_chg_array;
            dcir_soc_data.(channel).(rpt_key).charge_Q = Q_chg_array;
            dcir_soc_data.(channel).(rpt_key).charge_r_1s = R_chg(:,1);
            dcir_soc_data.(channel).(rpt_key).charge_r_5s = R_chg(:,2);
            dcir_soc_data.(channel).(rpt_key).charge_r_30s = R_chg(:,3);
            dcir_soc_data.(channel).(rpt_key).charge_r_60s = R_chg(:,4);
            dcir_soc_data.(channel).(rpt_key).charge_dr_1s0s = dR_chg(:,1);
            dcir_soc_data.(channel).(rpt_key).charge_dr_5s1s = dR_chg(:,2);
            dcir_soc_data.(channel).(rpt_key).charge_dr_60s5s = dR_chg(:,3);
            
            % Store charge summary table
            dcir_soc_data.(channel).(rpt_key).charge_summary_table = table(...
                SOC_chg, cumQ_chg_array, Q_chg_array, R_chg(:,1), R_chg(:,2), R_chg(:,3), R_chg(:,4), ...
                'VariableNames', {'SOC_pct', 'cumQ_Ah', 'Q_Ah', 'R_1s_mOhm', 'R_5s_mOhm', 'R_30s_mOhm', 'R_60s_mOhm'});
            
            % plot 60s point
            nexttile(1);
            color_val = getCycleColor(rpt_cycle);
            scatter(SOC_chg, R_chg(:,4), 50, color_val, 'filled', 'DisplayName', sprintf('%s 60s', rpt_cycle));
        end

        % Discharge data
        if ~isempty(dis_data_DCIR_selected)
            n_dis = length(dis_data_DCIR_selected);
            SOC_dis = zeros(n_dis, 1);
            R_dis = zeros(n_dis, 4);
            dR_dis = zeros(n_dis, 3);
            for i = 1:n_dis
                SOC_dis(i) = dis_data_DCIR_selected(i).SOC;
                R_dis(i,:) = dis_data_DCIR_selected(i).R;
                dR_dis(i,:) = dis_data_DCIR_selected(i).dR;
            end
            % cumQ, Q (concatenate steps)
            cumQ_dis_array = zeros(n_dis, 1);
            Q_dis_array = zeros(n_dis, 1);
            t_all = [];
            I_all_A = [];
            step_lengths = zeros(n_dis,1);
            for i = 1:n_dis
                t_rel = dis_data_DCIR_selected(i).steptime_double(:);
                % 60초까지만 사용 (max 값의 인덱스)
                max_idx = find(t_rel == max(t_rel), 1, 'last');
                step_lengths(i) = max_idx;
            end
            abs_offset = 0;
            for i = 1:n_dis
                t_rel = dis_data_DCIR_selected(i).steptime_double(:);
                if isempty(t_rel), t_rel = 0; end
                
                % 60초까지만 사용 (max 값의 인덱스)
                max_idx = find(t_rel == max(t_rel), 1, 'last');
                t_rel = t_rel(1:max_idx);
                I_step_A = abs(dis_data_DCIR_selected(i).I(1:max_idx));
                
                t_abs = abs_offset + t_rel;
                t_all = [t_all; t_abs]; %#ok<AGROW>
                I_all_A = [I_all_A; I_step_A]; %#ok<AGROW>
                abs_offset = t_abs(end);
            end
            if isduration(t_all), t_all = seconds(t_all); end
            cumQ_curve = cumtrapz(t_all, I_all_A) / 3600;
            Q_dis_total = trapz(t_all, I_all_A) / 3600; % Ah
            ends = cumsum(step_lengths);
            prev_end = 0;
            for i = 1:n_dis
                cumQ_dis_array(i) = cumQ_curve(ends(i));
                Q_dis_array(i) = Q_dis_total; % 모든 펄스에서 동일한 총 용량
                prev_end = ends(i);
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
            % store
            dcir_soc_data.(channel).(rpt_key).discharge_soc = SOC_dis;
            dcir_soc_data.(channel).(rpt_key).discharge_cumQ = cumQ_dis_array;
            dcir_soc_data.(channel).(rpt_key).discharge_Q = Q_dis_array;
            dcir_soc_data.(channel).(rpt_key).discharge_r_1s = R_dis(:,1);
            dcir_soc_data.(channel).(rpt_key).discharge_r_5s = R_dis(:,2);
            dcir_soc_data.(channel).(rpt_key).discharge_r_30s = R_dis(:,3);
            dcir_soc_data.(channel).(rpt_key).discharge_r_60s = R_dis(:,4);
            dcir_soc_data.(channel).(rpt_key).discharge_dr_1s0s = dR_dis(:,1);
            dcir_soc_data.(channel).(rpt_key).discharge_dr_5s1s = dR_dis(:,2);
            dcir_soc_data.(channel).(rpt_key).discharge_dr_60s5s = dR_dis(:,3);
            
            % Store discharge summary table
            dcir_soc_data.(channel).(rpt_key).discharge_summary_table = table(...
                SOC_dis, cumQ_dis_array, Q_dis_array, R_dis(:,1), R_dis(:,2), R_dis(:,3), R_dis(:,4), ...
                'VariableNames', {'SOC_pct', 'cumQ_Ah', 'Q_Ah', 'R_1s_mOhm', 'R_5s_mOhm', 'R_30s_mOhm', 'R_60s_mOhm'});
            
            % plot 60s point
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
    % (RPT_DCIR과 동일하게 채널 단위 MAT 저장은 생략)

    % Save individual cycle figures
    for rpt_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{rpt_idx};
        cyc_num = strrep(rpt_cycle, 'cyc', '');
        fig_individual = figure('Name', sprintf('DCIR vs SOC - %s %s', channel, rpt_cycle), 'Position', [100 100 1200 800]);
        tlo_individual = tiledlayout(fig_individual, 2, 1, 'TileSpacing','compact','Padding','compact');
        nexttile; hold on; grid on; title(sprintf('%s - Charge (DCIR) - %s', channel, rpt_cycle));
        xlabel('SOC [%]'); ylabel('DCIR [m\Omega]');
        ylim([0 5]);

        nexttile; hold on; grid on; title(sprintf('%s - Discharge (DCIR) - %s', channel, rpt_cycle));
        xlabel('SOC [%]'); ylabel('DCIR [m\Omega]');
        ylim([0 10]);

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
            if isfield(dcir_soc_data.(channel).(rpt_key), 'charge_soc')
                SOC_chg = dcir_soc_data.(channel).(rpt_key).charge_soc;
                DR1 = dcir_soc_data.(channel).(rpt_key).charge_dr_1s0s(:);
                DR2 = dcir_soc_data.(channel).(rpt_key).charge_dr_5s1s(:);
                DR3 = dcir_soc_data.(channel).(rpt_key).charge_dr_60s5s(:);
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

        if isfield(dcir_soc_data, channel) && isfield(dcir_soc_data.(channel), rpt_key)
            if isfield(dcir_soc_data.(channel).(rpt_key), 'discharge_soc')
                SOC_dis = dcir_soc_data.(channel).(rpt_key).discharge_soc;
                DR1 = dcir_soc_data.(channel).(rpt_key).discharge_dr_1s0s(:);
                DR2 = dcir_soc_data.(channel).(rpt_key).discharge_dr_5s1s(:);
                DR3 = dcir_soc_data.(channel).(rpt_key).discharge_dr_60s5s(:);
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

% Save combined MAT (RPT_DCIR과 동일한 이름)
matFile = fullfile(saveDir, 'DCIR_SOC_data_15to16.mat');
save(matFile, 'dcir_soc_data');
fprintf('\nSaved: %s\n', matFile);


