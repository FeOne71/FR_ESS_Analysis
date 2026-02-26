%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_DCIR_Analysis.m
% - DCIR 데이터 시각화 스크립트
% - 채널별 저항값 추이 시각화 (0cyc, 200cyc, 400cyc, 600cyc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Paths and settings
dataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_Final\DCIR_SOC_data_all_channels_final.mat';
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_Analysis';

if ~exist(saveDir,'dir'); mkdir(saveDir); end

%% Load DCIR data
if ~exist(dataPath,'file')
    error('DCIR data file not found: %s', dataPath);
end

load(dataPath);  % Loads: dcir_soc_data
fprintf('Loaded DCIR data from: %s\n', dataPath);

%% Configuration
channels = {'Ch09','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16'};
rpt_cycles = {'0cyc','200cyc','400cyc','600cyc','800cyc'};
resistance_types = {'R1_mOhm','R5_mOhm','R30_mOhm','R60_mOhm'};
resistance_names = {'R1 (1s)','R5 (5s)','R30 (30s)','R60 (60s)'};

%% Helper function: SOC 매칭
function [matched_baseline, matched_current, soc_matched] = matchSOCData(baseline_data, current_data)
    soc_tolerance = 2.0;  % SOC 허용 오차 (±2%)
    
    matched_baseline = [];
    matched_current = [];
    soc_matched = [];
    
    for i = 1:height(baseline_data)
        baseline_soc = baseline_data.SOC(i);
        
        % 가장 가까운 SOC를 가진 current 데이터 찾기
        soc_diffs = abs(current_data.SOC - baseline_soc);
        [min_diff, min_idx] = min(soc_diffs);
        
        if min_diff <= soc_tolerance
            matched_baseline = [matched_baseline; baseline_data(i, :)];
            matched_current = [matched_current; current_data(min_idx, :)];
            soc_matched = [soc_matched; baseline_soc];
        end
    end
end

%% 시각화 생성
fprintf('\n=== 시각화 생성 ===\n');

% 저항 타입별로 피겨 생성
for res_idx = 1:length(resistance_types)
    res_type = resistance_types{res_idx};
    res_name = resistance_names{res_idx};
    
    fig = figure('Name', sprintf('DCIR Analysis - %s', res_name), 'Position', [100 100 2000 1200]);
    tlo = tiledlayout(fig, 2, 8, 'TileSpacing','compact','Padding','compact');
    
    for ch_idx = 1:length(channels)
        channel = channels{ch_idx};
        
        % 충전 subplot (1행)
        nexttile(ch_idx);
        hold on; grid on;
        title(sprintf('%s - Charge', channel));
        xlabel('Cycle'); ylabel(sprintf('%s (mΩ)', res_name));
        xticks([1 2 3 4]); xticklabels({'0cyc','200cyc','400cyc','600cyc'});
        
        % 방전 subplot (2행)
        discharge_tile_idx = ch_idx + 8;
        nexttile(discharge_tile_idx);
        hold on; grid on;
        title(sprintf('%s - Discharge', channel));
        xlabel('Cycle'); ylabel(sprintf('%s (mΩ)', res_name));
        xticks([1 2 3 4]); xticklabels({'0cyc','200cyc','400cyc','600cyc'});
        
        % 충전 데이터 처리
        if isfield(dcir_soc_data, channel) && isfield(dcir_soc_data.(channel), 'cyc0') && ...
           isfield(dcir_soc_data.(channel).cyc0, 'charge_table')
            
            baseline_data = dcir_soc_data.(channel).cyc0.charge_table;
            cycle_means = [];
            cycle_mins = [];
            cycle_maxs = [];
            
            % 각 사이클별 평균 저항값 계산
            for rpt_idx = 1:length(rpt_cycles)
                rpt_cycle = rpt_cycles{rpt_idx};
                
                if strcmp(rpt_cycle, '0cyc')
                    rpt_key = 'cyc0';
                elseif strcmp(rpt_cycle, '200cyc')
                    rpt_key = 'cyc200';
                elseif strcmp(rpt_cycle, '400cyc')
                    rpt_key = 'cyc400';
                elseif strcmp(rpt_cycle, '600cyc')
                    rpt_key = 'cyc600';
                end
                
                if isfield(dcir_soc_data, channel) && isfield(dcir_soc_data.(channel), rpt_key) && ...
                   isfield(dcir_soc_data.(channel).(rpt_key), 'charge_table')
                    
                    current_data = dcir_soc_data.(channel).(rpt_key).charge_table;
                    
                    if strcmp(rpt_cycle, '0cyc')
                        % 0cyc는 직접 사용
                        vals = current_data.(res_type);
                    else
                        % 다른 사이클은 SOC 매칭
                        [matched_baseline, matched_current, ~] = matchSOCData(baseline_data, current_data);
                        if ~isempty(matched_current)
                            vals = matched_current.(res_type);
                        else
                            vals = [];
                        end
                    end
                    
                    % 디버깅: 값 확인
                    % fprintf('DEBUG: %s %s %s (Charge) - vals length: %d\n', channel, rpt_cycle, res_type, length(vals));
                    % if ~isempty(vals)
                    %     fprintf('  vals range: [%.4f, %.4f]\n', min(vals), max(vals));
                    %     if any(vals < 0)
                    %         fprintf('WARNING: Negative resistance values found!\n');
                    %         fprintf('  Negative values: %s\n', mat2str(vals(vals < 0)));
                    %     end
                    % end
                    
                    if ~isempty(vals)
                        cycle_means(end+1) = mean(vals);
                        cycle_mins(end+1) = min(vals);
                        cycle_maxs(end+1) = max(vals);
                    else
                        cycle_means(end+1) = NaN;
                        cycle_mins(end+1) = NaN;
                        cycle_maxs(end+1) = NaN;
                    end
                else
                    cycle_means(end+1) = NaN;
                    cycle_mins(end+1) = NaN;
                    cycle_maxs(end+1) = NaN;
                end
            end
            
            % 충전 플롯
            if ~isempty(cycle_means)
                nexttile(ch_idx);
                x_vals = 1:length(cycle_means);
                color_val = [0.2 0.2 0.2];
                
                plot(x_vals, cycle_means, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
                    'Color', color_val, 'MarkerFaceColor', color_val);
                
                % 실제 최솟값/최댓값으로 에러바 표시
                lower_errors = cycle_means - cycle_mins;  % 평균 - 최솟값
                upper_errors = cycle_maxs - cycle_means;  % 최댓값 - 평균
                errorbar(x_vals, cycle_means, lower_errors, upper_errors, 'Color', color_val, 'LineWidth', 1);
                
                % Y축 범위 설정 (실제 데이터 범위)
                valid_mins = cycle_mins(~isnan(cycle_mins));
                valid_maxs = cycle_maxs(~isnan(cycle_maxs));
                if ~isempty(valid_mins)
                    ylim_min = min(valid_mins);
                    % ylim_max = max(valid_maxs);
                    ylim([0, inf]);
                end
            end
        end
        
        % 방전 데이터 처리
        if isfield(dcir_soc_data, channel) && isfield(dcir_soc_data.(channel), 'cyc0') && ...
           isfield(dcir_soc_data.(channel).cyc0, 'discharge_table')
            
            baseline_data = dcir_soc_data.(channel).cyc0.discharge_table;
            cycle_means = [];
            cycle_mins = [];
            cycle_maxs = [];
            
            % 각 사이클별 평균 저항값 계산
            for rpt_idx = 1:length(rpt_cycles)
                rpt_cycle = rpt_cycles{rpt_idx};
                
                if strcmp(rpt_cycle, '0cyc')
                    rpt_key = 'cyc0';
                elseif strcmp(rpt_cycle, '200cyc')
                    rpt_key = 'cyc200';
                elseif strcmp(rpt_cycle, '400cyc')
                    rpt_key = 'cyc400';
                elseif strcmp(rpt_cycle, '600cyc')
                    rpt_key = 'cyc600';
                end
                
                if isfield(dcir_soc_data, channel) && isfield(dcir_soc_data.(channel), rpt_key) && ...
                   isfield(dcir_soc_data.(channel).(rpt_key), 'discharge_table')
                    
                    current_data = dcir_soc_data.(channel).(rpt_key).discharge_table;
                    
                    if strcmp(rpt_cycle, '0cyc')
                        % 0cyc는 직접 사용
                        vals = current_data.(res_type);
                    else
                        % 다른 사이클은 SOC 매칭
                        [matched_baseline, matched_current, ~] = matchSOCData(baseline_data, current_data);
                        if ~isempty(matched_current)
                            vals = matched_current.(res_type);
                        else
                            vals = [];
                        end
                    end
                    
                    % % 디버깅: 값 확인
                    % fprintf('DEBUG: %s %s %s (Discharge) - vals length: %d\n', channel, rpt_cycle, res_type, length(vals));
                    % if ~isempty(vals)
                    %     fprintf('  vals range: [%.4f, %.4f]\n', min(vals), max(vals));
                    %     if any(vals < 0)
                    %         fprintf('WARNING: Negative resistance values found!\n');
                    %         fprintf('  Negative values: %s\n', mat2str(vals(vals < 0)));
                    %     end
                    % end
                    
                    if ~isempty(vals)
                        cycle_means(end+1) = mean(vals);
                        cycle_mins(end+1) = min(vals);
                        cycle_maxs(end+1) = max(vals);
                    else
                        cycle_means(end+1) = NaN;
                        cycle_mins(end+1) = NaN;
                        cycle_maxs(end+1) = NaN;
                    end
                else
                    cycle_means(end+1) = NaN;
                    cycle_mins(end+1) = NaN;
                    cycle_maxs(end+1) = NaN;
                end
            end
            
            % 방전 플롯
            if ~isempty(cycle_means)
                nexttile(discharge_tile_idx);
                x_vals = 1:length(cycle_means);
                color_val = [0.2 0.2 0.2];
                
                plot(x_vals, cycle_means, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
                    'Color', color_val, 'MarkerFaceColor', color_val);
                
                % 실제 최솟값/최댓값으로 에러바 표시
                lower_errors = cycle_means - cycle_mins;  % 평균 - 최솟값
                upper_errors = cycle_maxs - cycle_means;  % 최댓값 - 평균
                errorbar(x_vals, cycle_means, lower_errors, upper_errors, 'Color', color_val, 'LineWidth', 1);
                
                % Y축 범위 설정 (실제 데이터 범위)
                valid_mins = cycle_mins(~isnan(cycle_mins));
                valid_maxs = cycle_maxs(~isnan(cycle_maxs));
                if ~isempty(valid_mins)
                    ylim_min = min(valid_mins);
                    ylim_max = max(valid_maxs);
                    ylim([0, inf]);
                end
            end
        end
    end
    
    % 피겨 저장
    savefig(fig, fullfile(saveDir, sprintf('DCIR_Analysis_%s.fig', res_type)));
    close(fig);
    fprintf('Saved: DCIR_Analysis_%s.fig\n', res_type);
end

fprintf('\n=== RPT_DCIR_Analysis.m 완료 ===\n');