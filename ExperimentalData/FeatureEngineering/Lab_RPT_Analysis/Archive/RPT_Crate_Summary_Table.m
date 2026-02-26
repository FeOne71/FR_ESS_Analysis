%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_Crate_Summary_Table.m
% C-rate_integrated.mat 파일에서 데이터를 읽어서
% 채널별/사이클별/C-rate별로 표를 생성
% - Charge/Discharge 각각 표 생성
% - 용량 범위, 전압 범위, 데이터 포인트 수 등 표시
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% ========================================================================
%  Configuration
% =========================================================================

% C-rate 데이터 파일 경로
crateDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Crate_integrated';
crateMatFile = fullfile(crateDataPath, 'Crate_integrated.mat');

% 저장 경로
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\Lab_RPT_Analysis\RPT_Crate_Postprocessing';

if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% C-rate 설정
crate_labels = {'c01', 'c05', 'c1', 'c2', 'c3'};
crate_names = {'0.1C', '0.5C', '1C', '2C', '3C'};
channels = {'Ch09', 'Ch10', 'Ch11', 'Ch12', 'Ch13', 'Ch14', 'Ch15', 'Ch16'};
rpt_cycles = {'0cyc', '200cyc', '400cyc', '600cyc', '800cyc', '1000cyc'};

%% ========================================================================
%  데이터 로드
% =========================================================================

fprintf('\n=== Loading C-rate Data ===\n');

if ~isfile(crateMatFile)
    error('C-rate data file not found: %s\nRun RPT_Postprocessing.m first.', crateMatFile);
end

load(crateMatFile, 'Crate_data');
fprintf('Loaded: %s\n', crateMatFile);

%% ========================================================================
%  원시 데이터 로드 검증
% =========================================================================

fprintf('\n=== Raw Data Load Verification ===\n');

% 상위 구조 확인
top_fields = fieldnames(Crate_data);
fprintf('  Top-level fields: %s\n', strjoin(top_fields, ', '));

% 채널별로 존재 여부
for ch_idx = 1:length(channels)
    channel_key = sprintf('Ch%s', channels{ch_idx}(3:end));
    if isfield(Crate_data, channel_key)
        cyc_fields = fieldnames(Crate_data.(channel_key));
        fprintf('  %s: %d cycle keys (%s ...)\n', channel_key, length(cyc_fields), strjoin(cyc_fields(1:min(2,end)), ', '));
    else
        fprintf('  %s: MISSING\n', channel_key);
    end
end

% 샘플 1개: Ch09 / 0cyc / c01 의 charge, discharge 원시 필드 확인
sample_ch = 'Ch09';
sample_cyc = 'cyc0';
sample_label = 'c01';
fprintf('\n  Sample path: %s / %s / %s\n', sample_ch, sample_cyc, sample_label);

if isfield(Crate_data, sample_ch) && isfield(Crate_data.(sample_ch), sample_cyc) && ...
   isfield(Crate_data.(sample_ch).(sample_cyc), sample_label)
    node = Crate_data.(sample_ch).(sample_cyc).(sample_label);
    for side = {'charge', 'discharge'}
        s = side{1};
        if isfield(node, s)
            d = node.(s);
            fns = fieldnames(d);
            fprintf('    %s fields: %s\n', s, strjoin(fns, ', '));
            if isfield(d, 'Q') && isfield(d, 'V')
                fprintf('      Q: n=%d, range [%.4f, %.4f] Ah\n', numel(d.Q), min(d.Q(:)), max(d.Q(:)));
                fprintf('      V: n=%d, range [%.4f, %.4f] V\n', numel(d.V), min(d.V(:)), max(d.V(:)));
            end
        else
            fprintf('    %s: MISSING\n', s);
        end
    end
else
    fprintf('    Sample path not found in Crate_data.\n');
end

% 전체 조합 중 유효 데이터 개수 요약
n_charge_ok = 0;
n_discharge_ok = 0;
for ch_idx = 1:length(channels)
    channel_key = sprintf('Ch%s', channels{ch_idx}(3:end));
    for cyc_idx = 1:length(rpt_cycles)
        cycle_key = sprintf('cyc%s', rpt_cycles{cyc_idx}(1:end-3));
        for rate_idx = 1:length(crate_labels)
            label = crate_labels{rate_idx};
            if isfield(Crate_data, channel_key) && isfield(Crate_data.(channel_key), cycle_key) && ...
               isfield(Crate_data.(channel_key).(cycle_key), label)
                chg = Crate_data.(channel_key).(cycle_key).(label);
                if isfield(chg, 'charge') && isfield(chg.charge, 'Q') && ~isempty(chg.charge.Q)
                    n_charge_ok = n_charge_ok + 1;
                end
                if isfield(chg, 'discharge') && isfield(chg.discharge, 'Q') && ~isempty(chg.discharge.Q)
                    n_discharge_ok = n_discharge_ok + 1;
                end
            end
        end
    end
end
n_total = length(channels) * length(rpt_cycles) * length(crate_labels);
fprintf('\n  Valid raw data count: charge %d / %d, discharge %d / %d\n', n_charge_ok, n_total, n_discharge_ok, n_total);

%% ========================================================================
%  Charge 데이터 표 생성
% =========================================================================

fprintf('\n=== Creating Charge Summary Table ===\n');

% 표 데이터 저장용 셀 배열
charge_table_data = {};
charge_table_data{1,1} = 'Channel';
charge_table_data{1,2} = 'Cycle';
charge_table_data{1,3} = 'C-rate';
charge_table_data{1,4} = 'Capacity Min [Ah]';
charge_table_data{1,5} = 'Capacity Max [Ah]';
charge_table_data{1,6} = 'Voltage Min [V]';
charge_table_data{1,7} = 'Voltage Max [V]';
charge_table_data{1,8} = 'Data Points';

row_idx = 2;

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    channel_key = sprintf('Ch%s', channel(3:end));  % 'Ch09' -> 'Ch09' (동일)
    
    for cyc_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{cyc_idx};
        cycle_key = sprintf('cyc%s', rpt_cycle(1:end-3));  % '0cyc' -> 'cyc0'
        
        for rate_idx = 1:length(crate_labels)
            label = crate_labels{rate_idx};
            rate_name = crate_names{rate_idx};
            
            % 데이터 확인 (원시 Q, V만 사용)
            chg = [];
            if isfield(Crate_data, channel_key) && ...
               isfield(Crate_data.(channel_key), cycle_key) && ...
               isfield(Crate_data.(channel_key).(cycle_key), label) && ...
               isfield(Crate_data.(channel_key).(cycle_key).(label), 'charge')
                chg = Crate_data.(channel_key).(cycle_key).(label).charge;
            end
            if ~isempty(chg) && isfield(chg, 'Q') && ~isempty(chg.Q)
                cap_vals = chg.Q;
                volt_vals = chg.V;
            else
                cap_vals = [];
                volt_vals = [];
            end
            if ~isempty(cap_vals) && ~isempty(volt_vals)
                % 유효한 데이터만
                valid_mask = ~isnan(cap_vals) & ~isnan(volt_vals);
                cap_valid = cap_vals(valid_mask);
                volt_valid = volt_vals(valid_mask);
                
                if ~isempty(cap_valid)
                    charge_table_data{row_idx,1} = channel;
                    charge_table_data{row_idx,2} = rpt_cycle;
                    charge_table_data{row_idx,3} = rate_name;
                    charge_table_data{row_idx,4} = min(cap_valid);
                    charge_table_data{row_idx,5} = max(cap_valid);
                    charge_table_data{row_idx,6} = min(volt_valid);
                    charge_table_data{row_idx,7} = max(volt_valid);
                    charge_table_data{row_idx,8} = length(cap_valid);
                    row_idx = row_idx + 1;
                else
                    charge_table_data{row_idx,1} = channel;
                    charge_table_data{row_idx,2} = rpt_cycle;
                    charge_table_data{row_idx,3} = rate_name;
                    charge_table_data{row_idx,4} = NaN;
                    charge_table_data{row_idx,5} = NaN;
                    charge_table_data{row_idx,6} = NaN;
                    charge_table_data{row_idx,7} = NaN;
                    charge_table_data{row_idx,8} = 0;
                    row_idx = row_idx + 1;
                end
            else
                % 데이터 없음
                charge_table_data{row_idx,1} = channel;
                charge_table_data{row_idx,2} = rpt_cycle;
                charge_table_data{row_idx,3} = rate_name;
                charge_table_data{row_idx,4} = NaN;
                charge_table_data{row_idx,5} = NaN;
                charge_table_data{row_idx,6} = NaN;
                charge_table_data{row_idx,7} = NaN;
                charge_table_data{row_idx,8} = 0;
                row_idx = row_idx + 1;
            end
        end
    end
end

% 표 생성 및 저장
charge_table = cell2table(charge_table_data(2:end,:), 'VariableNames', charge_table_data(1,:));
writetable(charge_table, fullfile(saveDir, 'Crate_Charge_Summary_Table.csv'));
fprintf('Saved: Crate_Charge_Summary_Table.csv\n');

% 표 출력
fprintf('\n--- Charge Summary Table (first 20 rows) ---\n');
disp(charge_table(1:min(20, height(charge_table)), :));

%% ========================================================================
%  Discharge 데이터 표 생성
% =========================================================================

fprintf('\n=== Creating Discharge Summary Table ===\n');

% 표 데이터 저장용 셀 배열
discharge_table_data = {};
discharge_table_data{1,1} = 'Channel';
discharge_table_data{1,2} = 'Cycle';
discharge_table_data{1,3} = 'C-rate';
discharge_table_data{1,4} = 'Capacity Min [Ah]';
discharge_table_data{1,5} = 'Capacity Max [Ah]';
discharge_table_data{1,6} = 'Voltage Min [V]';
discharge_table_data{1,7} = 'Voltage Max [V]';
discharge_table_data{1,8} = 'Data Points';

row_idx = 2;

for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    channel_key = sprintf('Ch%s', channel(3:end));
    
    for cyc_idx = 1:length(rpt_cycles)
        rpt_cycle = rpt_cycles{cyc_idx};
        cycle_key = sprintf('cyc%s', rpt_cycle(1:end-3));
        
        for rate_idx = 1:length(crate_labels)
            label = crate_labels{rate_idx};
            rate_name = crate_names{rate_idx};
            
            % 데이터 확인 (원시 Q, V만 사용)
            dch = [];
            if isfield(Crate_data, channel_key) && ...
               isfield(Crate_data.(channel_key), cycle_key) && ...
               isfield(Crate_data.(channel_key).(cycle_key), label) && ...
               isfield(Crate_data.(channel_key).(cycle_key).(label), 'discharge')
                dch = Crate_data.(channel_key).(cycle_key).(label).discharge;
            end
            if ~isempty(dch) && isfield(dch, 'Q') && ~isempty(dch.Q)
                cap_vals = dch.Q;
                volt_vals = dch.V;
            else
                cap_vals = [];
                volt_vals = [];
            end
            if ~isempty(cap_vals) && ~isempty(volt_vals)
                % 유효한 데이터만
                valid_mask = ~isnan(cap_vals) & ~isnan(volt_vals);
                cap_valid = cap_vals(valid_mask);
                volt_valid = volt_vals(valid_mask);
                
                if ~isempty(cap_valid)
                    discharge_table_data{row_idx,1} = channel;
                    discharge_table_data{row_idx,2} = rpt_cycle;
                    discharge_table_data{row_idx,3} = rate_name;
                    discharge_table_data{row_idx,4} = min(cap_valid);
                    discharge_table_data{row_idx,5} = max(cap_valid);
                    discharge_table_data{row_idx,6} = min(volt_valid);
                    discharge_table_data{row_idx,7} = max(volt_valid);
                    discharge_table_data{row_idx,8} = length(cap_valid);
                    row_idx = row_idx + 1;
                else
                    discharge_table_data{row_idx,1} = channel;
                    discharge_table_data{row_idx,2} = rpt_cycle;
                    discharge_table_data{row_idx,3} = rate_name;
                    discharge_table_data{row_idx,4} = NaN;
                    discharge_table_data{row_idx,5} = NaN;
                    discharge_table_data{row_idx,6} = NaN;
                    discharge_table_data{row_idx,7} = NaN;
                    discharge_table_data{row_idx,8} = 0;
                    row_idx = row_idx + 1;
                end
            else
                % 데이터 없음
                discharge_table_data{row_idx,1} = channel;
                discharge_table_data{row_idx,2} = rpt_cycle;
                discharge_table_data{row_idx,3} = rate_name;
                discharge_table_data{row_idx,4} = NaN;
                discharge_table_data{row_idx,5} = NaN;
                discharge_table_data{row_idx,6} = NaN;
                discharge_table_data{row_idx,7} = NaN;
                discharge_table_data{row_idx,8} = 0;
                row_idx = row_idx + 1;
            end
        end
    end
end

% 표 생성 및 저장
discharge_table = cell2table(discharge_table_data(2:end,:), 'VariableNames', discharge_table_data(1,:));
writetable(discharge_table, fullfile(saveDir, 'Crate_Discharge_Summary_Table.csv'));
fprintf('Saved: Crate_Discharge_Summary_Table.csv\n');

% 표 출력
fprintf('\n--- Discharge Summary Table (first 20 rows) ---\n');
disp(discharge_table(1:min(20, height(discharge_table)), :));

%% ========================================================================
%  Total8 (8채널 평균) 데이터 표 생성
% =========================================================================

fprintf('\n=== Creating Total8 (8-Channel Average) Summary Table ===\n');

% Charge Total8
total8_charge_table_data = {};
total8_charge_table_data{1,1} = 'Cycle';
total8_charge_table_data{1,2} = 'C-rate';
total8_charge_table_data{1,3} = 'Capacity Min [Ah]';
total8_charge_table_data{1,4} = 'Capacity Max [Ah]';
total8_charge_table_data{1,5} = 'Voltage Min [V]';
total8_charge_table_data{1,6} = 'Voltage Max [V]';
total8_charge_table_data{1,7} = 'Data Points';

row_idx = 2;

for cyc_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{cyc_idx};
    cycle_key = sprintf('cyc%s', rpt_cycle(1:end-3));
    
    for rate_idx = 1:length(crate_labels)
        label = crate_labels{rate_idx};
        rate_name = crate_names{rate_idx};
        
        chg = [];
        if isfield(Crate_data, 'total8') && ...
           isfield(Crate_data.total8, cycle_key) && ...
           isfield(Crate_data.total8.(cycle_key), label) && ...
           isfield(Crate_data.total8.(cycle_key).(label), 'charge')
            chg = Crate_data.total8.(cycle_key).(label).charge;
        end
        if ~isempty(chg) && isfield(chg, 'Q') && ~isempty(chg.Q)
            cap_vals = chg.Q;
            volt_vals = chg.V;
        else
            cap_vals = [];
            volt_vals = [];
        end
        if ~isempty(cap_vals) && ~isempty(volt_vals)
            valid_mask = ~isnan(cap_vals) & ~isnan(volt_vals);
            cap_valid = cap_vals(valid_mask);
            volt_valid = volt_vals(valid_mask);
            if ~isempty(cap_valid)
                total8_charge_table_data{row_idx,1} = rpt_cycle;
                total8_charge_table_data{row_idx,2} = rate_name;
                total8_charge_table_data{row_idx,3} = min(cap_valid);
                total8_charge_table_data{row_idx,4} = max(cap_valid);
                total8_charge_table_data{row_idx,5} = min(volt_valid);
                total8_charge_table_data{row_idx,6} = max(volt_valid);
                total8_charge_table_data{row_idx,7} = length(cap_valid);
                row_idx = row_idx + 1;
            end
        end
    end
end

if row_idx > 2
    total8_charge_table = cell2table(total8_charge_table_data(2:end,:), 'VariableNames', total8_charge_table_data(1,:));
    writetable(total8_charge_table, fullfile(saveDir, 'Crate_Total8_Charge_Summary_Table.csv'));
    fprintf('Saved: Crate_Total8_Charge_Summary_Table.csv\n');
    disp(total8_charge_table);
end

% Discharge Total8
total8_discharge_table_data = {};
total8_discharge_table_data{1,1} = 'Cycle';
total8_discharge_table_data{1,2} = 'C-rate';
total8_discharge_table_data{1,3} = 'Capacity Min [Ah]';
total8_discharge_table_data{1,4} = 'Capacity Max [Ah]';
total8_discharge_table_data{1,5} = 'Voltage Min [V]';
total8_discharge_table_data{1,6} = 'Voltage Max [V]';
total8_discharge_table_data{1,7} = 'Data Points';

row_idx = 2;

for cyc_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{cyc_idx};
    cycle_key = sprintf('cyc%s', rpt_cycle(1:end-3));
    
    for rate_idx = 1:length(crate_labels)
        label = crate_labels{rate_idx};
        rate_name = crate_names{rate_idx};
        
        dch = [];
        if isfield(Crate_data, 'total8') && ...
           isfield(Crate_data.total8, cycle_key) && ...
           isfield(Crate_data.total8.(cycle_key), label) && ...
           isfield(Crate_data.total8.(cycle_key).(label), 'discharge')
            dch = Crate_data.total8.(cycle_key).(label).discharge;
        end
        if ~isempty(dch) && isfield(dch, 'Q') && ~isempty(dch.Q)
            cap_vals = dch.Q;
            volt_vals = dch.V;
        else
            cap_vals = [];
            volt_vals = [];
        end
        if ~isempty(cap_vals) && ~isempty(volt_vals)
            valid_mask = ~isnan(cap_vals) & ~isnan(volt_vals);
            cap_valid = cap_vals(valid_mask);
            volt_valid = volt_vals(valid_mask);
            if ~isempty(cap_valid)
                total8_discharge_table_data{row_idx,1} = rpt_cycle;
                total8_discharge_table_data{row_idx,2} = rate_name;
                total8_discharge_table_data{row_idx,3} = min(cap_valid);
                total8_discharge_table_data{row_idx,4} = max(cap_valid);
                total8_discharge_table_data{row_idx,5} = min(volt_valid);
                total8_discharge_table_data{row_idx,6} = max(volt_valid);
                total8_discharge_table_data{row_idx,7} = length(cap_valid);
                row_idx = row_idx + 1;
            end
        end
    end
end

if row_idx > 2
    total8_discharge_table = cell2table(total8_discharge_table_data(2:end,:), 'VariableNames', total8_discharge_table_data(1,:));
    writetable(total8_discharge_table, fullfile(saveDir, 'Crate_Total8_Discharge_Summary_Table.csv'));
    fprintf('Saved: Crate_Total8_Discharge_Summary_Table.csv\n');
    disp(total8_discharge_table);
end

fprintf('\n=== C-rate Summary Table Generation Complete ===\n');
