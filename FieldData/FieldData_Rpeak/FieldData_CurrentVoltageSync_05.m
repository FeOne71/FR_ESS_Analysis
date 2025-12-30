%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FieldData_CurrentVoltageSync_05.m
% 전류 증가 시 전압도 동시에 증가하는지 확인하는 분석 스크립트
% 
% 목적: 
% - 검출된 이벤트에서 전류가 증가할 때, 전압도 동일한 시간(인덱스)에 증가하는지 확인
% - 예) 전류: [0 10 12 10 10 10] -> 전압: [3.65 3.66 3.67 3.68 3.68 3.68]
%
% 출력:
% - 분석 결과 요약 (콘솔 출력)
% - 상세 분석 결과 (MAT 파일)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Configuration
% =========================================================================
scriptDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(scriptDir, 'EventsResults');
outputDir = fullfile(resultsDir, 'FieldData_Analysis');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

targetRack = 'Rack01';
targetEventType = 'charge';  % 'charge' or 'discharge'

% 정리된 이벤트 파일명
resultStructName = sprintf('FieldData_%s_%s', targetRack, targetEventType);
eventsFile = sprintf('%s_Events.mat', resultStructName);
% =========================================================================

fprintf('=== Current-Voltage Synchronization Analysis ===\n');
fprintf('Target Rack: %s\n', targetRack);
fprintf('Target Event Type: %s\n', targetEventType);
fprintf('\n');

%% Load Processed Events Data
fprintf('=== Loading Processed Events Data ===\n');
eventsPath = fullfile(outputDir, eventsFile);
if ~exist(eventsPath, 'file')
    error('Processed events file not found: %s', eventsPath);
end

load(eventsPath, resultStructName);
fprintf('Data loaded successfully.\n');

% 동적으로 구조체 접근
eventsData = eval(resultStructName);
if ~isstruct(eventsData)
    error('Invalid data structure in file: %s', eventsPath);
end

%% Collect Events
fprintf('=== Collecting Events ===\n');
all_events_list = {};

% 그룹별로 순회
groupNames = fieldnames(eventsData);
for g_idx = 1:length(groupNames)
    groupName = groupNames{g_idx};
    groupData = eventsData.(groupName);
    
    % 연도별로 순회
    yearNames = fieldnames(groupData);
    for y_idx = 1:length(yearNames)
        yearName = yearNames{y_idx};
        yearData = groupData.(yearName);
        
        % 이벤트별로 순회
        eventNames = fieldnames(yearData);
        for e_idx = 1:length(eventNames)
            eventName = eventNames{e_idx};
            eventData = yearData.(eventName);
            
            % voltage_seq와 current_seq가 있는지 확인
            if isfield(eventData, 'voltage_seq') && isfield(eventData, 'current_seq')
                % 이벤트 정보를 구조체로 저장
                evt = struct();
                evt.group = groupName;
                evt.year = yearName;
                evt.eventName = eventName;
                evt.voltage_seq = eventData.voltage_seq;
                evt.current_seq = eventData.current_seq;
                if isfield(eventData, 'timestamp')
                    evt.timestamp = eventData.timestamp;
                end
                if isfield(eventData, 'time_seq')
                    evt.time_seq = eventData.time_seq;
                end
                
                all_events_list{end+1} = evt;
            end
        end
    end
end

fprintf('Total events collected: %d\n', length(all_events_list));
fprintf('\n');

%% Analyze Current-Voltage Synchronization
fprintf('=== Analyzing Current-Voltage Synchronization ===\n');

% 결과 저장용 구조체
analysis_results = struct();
analysis_results.total_events = length(all_events_list);
analysis_results.events_with_current_increase = 0;
analysis_results.events_with_sync_increase = 0;  % 전류 증가 시 전압도 증가
analysis_results.events_with_async_increase = 0;  % 전류 증가 시 전압은 증가하지 않음
analysis_results.event_details = [];

% 각 이벤트 분석
for e_idx = 1:length(all_events_list)
    evt = all_events_list{e_idx};
    
    % 이벤트 데이터 추출
    if isfield(evt, 'current_seq') && isfield(evt, 'voltage_seq')
        I = evt.current_seq;
        V = evt.voltage_seq;
    else
        continue;  % 데이터가 없으면 스킵
    end
    
    if length(I) ~= length(V) || length(I) < 2
        continue;  % 길이가 다르거나 너무 짧으면 스킵
    end
    
    % 각 인덱스에서 전류가 증가하면 전압도 증가하는지 확인
    sync_count = 0;  % 동기화된 인덱스 수 (전류 증가 → 전압 증가)
    async_count = 0;  % 비동기화된 인덱스 수 (전류 증가 → 전압 증가 안 함)
    total_increase_count = 0;  % 전류가 증가한 총 인덱스 수
    
    % 각 인덱스별 상세 정보 저장
    index_details = [];
    
    for i = 2:length(I)
        % 전류가 증가하는 경우만 확인
        if I(i) > I(i-1)
            total_increase_count = total_increase_count + 1;
            
            % 해당 인덱스에서 전압도 증가하는지 확인
            if V(i) > V(i-1)
                sync_count = sync_count + 1;
                is_sync_idx = true;
            else
                async_count = async_count + 1;
                is_sync_idx = false;
            end
            
            % 상세 정보 저장
            idx_detail = struct();
            idx_detail.index = i;
            idx_detail.I_prev = I(i-1);
            idx_detail.I_curr = I(i);
            idx_detail.I_delta = I(i) - I(i-1);
            idx_detail.V_prev = V(i-1);
            idx_detail.V_curr = V(i);
            idx_detail.V_delta = V(i) - V(i-1);
            idx_detail.is_sync = is_sync_idx;
            index_details = [index_details; idx_detail];
        end
    end
    
    % 전류 증가가 없는 이벤트는 스킵
    if total_increase_count == 0
        continue;
    end
    
    analysis_results.events_with_current_increase = analysis_results.events_with_current_increase + 1;
    
    % 이벤트 전체적으로 동기화 여부 판단
    % 전류 증가 인덱스 중 절반 이상에서 전압도 증가하면 동기화된 것으로 판단
    is_sync = sync_count > async_count;
    
    if is_sync
        analysis_results.events_with_sync_increase = analysis_results.events_with_sync_increase + 1;
    else
        analysis_results.events_with_async_increase = analysis_results.events_with_async_increase + 1;
    end
    
    % 상세 정보 저장
    event_detail = struct();
    event_detail.event_idx = e_idx;
    if isfield(evt, 'group')
        event_detail.group = evt.group;
    else
        event_detail.group = 'unknown';
    end
    if isfield(evt, 'year')
        event_detail.year = evt.year;
    else
        event_detail.year = 'unknown';
    end
    if isfield(evt, 'eventName')
        event_detail.eventName = evt.eventName;
    else
        event_detail.eventName = 'unknown';
    end
    if isfield(evt, 'timestamp')
        event_detail.timestamp = evt.timestamp;
    end
    event_detail.total_increase_indices = total_increase_count;
    event_detail.sync_indices = sync_count;
    event_detail.async_indices = async_count;
    event_detail.is_sync = is_sync;
    event_detail.index_details = index_details;
    
    % 시각화를 위한 원본 데이터 저장
    event_detail.I_data = I;
    event_detail.V_data = V;
    if isfield(evt, 'time_seq')
        event_detail.time_seq = evt.time_seq;
    end
    
    analysis_results.event_details = [analysis_results.event_details; event_detail];
    
    % 진행 상황 출력 (100개마다)
    if mod(e_idx, 100) == 0
        fprintf('  Processed %d / %d events...\n', e_idx, length(all_events_list));
    end
end

fprintf('Analysis complete!\n');
fprintf('\n');

%% Display Results
fprintf('=== Analysis Results ===\n');
fprintf('Total events analyzed: %d\n', analysis_results.total_events);
fprintf('Events with current increase: %d\n', analysis_results.events_with_current_increase);
fprintf('Events with synchronized increase (current ↑ → voltage ↑): %d (%.1f%%)\n', ...
    analysis_results.events_with_sync_increase, ...
    100 * analysis_results.events_with_sync_increase / max(analysis_results.events_with_current_increase, 1));
fprintf('Events with async increase (current ↑ → voltage not ↑): %d (%.1f%%)\n', ...
    analysis_results.events_with_async_increase, ...
    100 * analysis_results.events_with_async_increase / max(analysis_results.events_with_current_increase, 1));
fprintf('\n');

% 상세 통계
if ~isempty(analysis_results.event_details)
    total_indices_all = sum([analysis_results.event_details.total_increase_indices]);
    sync_indices_all = sum([analysis_results.event_details.sync_indices]);
    async_indices_all = sum([analysis_results.event_details.async_indices]);
    
    fprintf('=== Index-level Statistics ===\n');
    fprintf('Total current increase indices: %d\n', total_indices_all);
    fprintf('Synchronized indices (current ↑ → voltage ↑): %d (%.1f%%)\n', ...
        sync_indices_all, 100 * sync_indices_all / max(total_indices_all, 1));
    fprintf('Async indices (current ↑ → voltage not ↑): %d (%.1f%%)\n', ...
        async_indices_all, 100 * async_indices_all / max(total_indices_all, 1));
    fprintf('\n');
    
    % 비동기화된 이벤트 예시 출력 (최대 5개)
    async_events = analysis_results.event_details([analysis_results.event_details.is_sync] == false);
    if ~isempty(async_events)
        fprintf('=== Example Async Events (first 5) ===\n');
        num_examples = min(5, length(async_events));
        for i = 1:num_examples
            evt_detail = async_events(i);
            fprintf('\nEvent %d:\n', evt_detail.event_idx);
            fprintf('  Group: %s\n', evt_detail.group);
            fprintf('  Year: %s, Event: %s\n', evt_detail.year, evt_detail.eventName);
            fprintf('  Total increase indices: %d (sync: %d, async: %d)\n', ...
                evt_detail.total_increase_indices, evt_detail.sync_indices, evt_detail.async_indices);
            
            % 첫 번째 비동기화 인덱스 상세 정보
            if ~isempty(evt_detail.index_details)
                for idx = 1:min(3, length(evt_detail.index_details))
                    idx_detail = evt_detail.index_details(idx);
                    if ~idx_detail.is_sync
                        fprintf('  Index %d (async): I [%.2f → %.2f] (Δ=%.2f A), V [%.3f → %.3f] (Δ=%.3f V)\n', ...
                            idx_detail.index, idx_detail.I_prev, idx_detail.I_curr, idx_detail.I_delta, ...
                            idx_detail.V_prev, idx_detail.V_curr, idx_detail.V_delta);
                        break;
                    end
                end
            end
        end
        fprintf('\n');
    end
end

%% Visualization
fprintf('=== Creating Visualizations ===\n');

% 비동기화된 이벤트 찾기
async_events = analysis_results.event_details([analysis_results.event_details.is_sync] == false);
if ~isempty(async_events)
    fprintf('Found %d async events. Creating visualizations...\n', length(async_events));
    
    % 시각화할 이벤트 수 (최대 20개)
    num_events_to_plot = min(20, length(async_events));
    
    % Figure 저장 경로
    figDir = fullfile(outputDir, 'figures', 'CurrentVoltageSync');
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
    
    % 각 비동기화 이벤트 시각화
    for plot_idx = 1:num_events_to_plot
        evt_detail = async_events(plot_idx);
        
        % 저장된 데이터 사용
        I_full = evt_detail.I_data;
        V_full = evt_detail.V_data;
        
        % 초반 60초(또는 60개 인덱스)만 사용
        max_duration_sec = 60;
        max_indices = 60;
        
        % 시간 시퀀스 확인 및 데이터 자르기
        if isfield(evt_detail, 'time_seq') && ~isempty(evt_detail.time_seq)
            if isdatetime(evt_detail.time_seq(1))
                time_start = evt_detail.time_seq(1);
                t_rel_full = seconds(evt_detail.time_seq - time_start);
                
                % 60초 이내 인덱스 찾기
                valid_idx = t_rel_full <= max_duration_sec;
                if ~any(valid_idx)
                    valid_idx = 1:min(max_indices, length(I_full));
                end
                
                I = I_full(valid_idx);
                V = V_full(valid_idx);
                t_rel = t_rel_full(valid_idx);
                xlabel_str = 'Time (seconds)';
            else
                % 인덱스 기반: 처음 60개만
                valid_idx = 1:min(max_indices, length(I_full));
                I = I_full(valid_idx);
                V = V_full(valid_idx);
                t_rel = 1:length(I);
                xlabel_str = 'Index';
            end
        else
            % 인덱스 기반: 처음 60개만
            valid_idx = 1:min(max_indices, length(I_full));
            I = I_full(valid_idx);
            V = V_full(valid_idx);
            t_rel = 1:length(I);
            xlabel_str = 'Index';
        end
        
        % Figure 생성
        fig = figure('Position', [100, 100, 1200, 800]);
        
        % 서브플롯 1: 전류와 전압 함께 플롯
        ax1 = subplot(2, 1, 1);
        yyaxis left
        plot(t_rel, I, 'b-', 'LineWidth', 1.5);
        ylabel('Current (A)', 'FontSize', 12, 'FontWeight', 'bold');
        yyaxis right
        plot(t_rel, V, 'r-', 'LineWidth', 1.5);
        ylabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel(xlabel_str, 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Event %d: %s / %s / %s (First 60s)', evt_detail.event_idx, ...
            evt_detail.group, evt_detail.year, evt_detail.eventName), ...
            'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        legend('Current', 'Voltage', 'Location', 'best');
        
        % 동기화/비동기화 인덱스 표시 (초반 60초/60개 인덱스 내에서만)
        sync_indices = [];
        async_indices = [];
        if ~isempty(evt_detail.index_details)
            for idx = 1:length(evt_detail.index_details)
                idx_detail = evt_detail.index_details(idx);
                % 초반 60초/60개 인덱스 내에 있는지 확인
                if idx_detail.index <= length(I)
                    if idx_detail.is_sync
                        sync_indices = [sync_indices, idx_detail.index];
                    else
                        async_indices = [async_indices, idx_detail.index];
                    end
                end
            end
        end
        
        % 비동기화 인덱스에 빨간 세로선 표시
        if ~isempty(async_indices)
            hold on;
            for idx = 1:length(async_indices)
                idx_val = async_indices(idx);
                if idx_val <= length(t_rel)
                    xline(t_rel(idx_val), 'r--', 'LineWidth', 1.5);
                end
            end
            hold off;
        end
        
        % 서브플롯 2: 전류 증가와 전압 변화 비교
        ax2 = subplot(2, 1, 2);
        
        % 전류 변화량
        I_delta = diff(I);
        I_delta = [I_delta(1); I_delta];  % 첫 번째 값 유지
        
        % 전압 변화량
        V_delta = diff(V);
        V_delta = [V_delta(1); V_delta];  % 첫 번째 값 유지
        
        % 전류 증가 인덱스만 표시
        I_increase_mask = I_delta > 0;
        t_increase = t_rel(I_increase_mask);
        I_delta_increase = I_delta(I_increase_mask);
        V_delta_increase = V_delta(I_increase_mask);
        
        % 동기화/비동기화 구분
        sync_mask = V_delta_increase > 0;
        async_mask = ~sync_mask;
        
        hold on;
        % 동기화된 인덱스: 파란색
        if any(sync_mask)
            scatter(t_increase(sync_mask), I_delta_increase(sync_mask), 50, 'b', 'filled', 'DisplayName', 'Sync (I↑, V↑)');
        end
        % 비동기화된 인덱스: 빨간색
        if any(async_mask)
            scatter(t_increase(async_mask), I_delta_increase(async_mask), 50, 'r', 'filled', 'DisplayName', 'Async (I↑, V not ↑)');
        end
        hold off;
        
        xlabel(xlabel_str, 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Current Change (A)', 'FontSize', 12, 'FontWeight', 'bold');
        title('Current Increase Points: Sync vs Async', 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'best');
        grid on;
        
        % 전체 제목
        sgtitle(sprintf('Current-Voltage Sync Analysis (Event %d: %d sync, %d async)', ...
            evt_detail.event_idx, evt_detail.sync_indices, evt_detail.async_indices), ...
            'FontSize', 16, 'FontWeight', 'bold');
        
        % Figure 저장
        figName = sprintf('CurrentVoltageSync_Event%d_%s_%s_%s.fig', ...
            evt_detail.event_idx, evt_detail.group, evt_detail.year, evt_detail.eventName);
        figPath = fullfile(figDir, figName);
        savefig(fig, figPath);
        fprintf('  Saved: %s\n', figName);
        
        close(fig);
    end
    
    fprintf('Visualizations saved to: %s\n', figDir);
    fprintf('\n');
else
    fprintf('No async events found. Skipping visualization.\n');
    fprintf('\n');
end

%% Save Results
fprintf('=== Saving Results ===\n');
savePath = fullfile(outputDir, sprintf('CurrentVoltageSync_Analysis_%s.mat', targetEventType));
save(savePath, 'analysis_results', '-v7.3');
fprintf('Results saved to: %s\n', savePath);
fprintf('\n');

fprintf('=== Analysis Complete ===\n');

