%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT_Discharge_dQdV_Visualization.m
% RPT 데이터에서 방전(Discharge) 부분만 추출하여 dQ/dV 곡선 시각화
% - RPT 데이터만 사용 (Aging 데이터 제외)
% - Step 10 (Discharge OCV), CycleIdx 2 사용
% - 여러 사이클(0cyc, 200cyc, 400cyc, 600cyc) 비교
% - 8채널 평균 및 개별 채널 시각화
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% ========================================================================
%  Configuration
% =========================================================================

% 데이터 경로
dataDir = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data\RPT';
saveDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\dQdV_Discharge';
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% 채널 및 사이클 설정
channels = {'09', '10', '11', '12', '13', '14', '15', '16'};
rptCycles = [0, 200, 400, 600,800];  % 분석할 RPT 사이클

% RPT 데이터 설정
discharge_step = 3;    % Step 10: Discharge OCV || Step 3: Discharge Static Cap
cycle_idx = 2;          % Cycle Index 2          || Cycle Index 2

%% ========================================================================
%  데이터 로드 및 dQ/dV 계산
% =========================================================================

fprintf('\n=== Loading RPT Discharge Data ===\n');

% 각 채널별, 사이클별 데이터 저장 구조체
all_dQdV_data = struct();

for chIdx = 1:length(channels)
    ch = channels{chIdx};
    channel_key = sprintf('Ch%s', ch);
    all_dQdV_data.(channel_key) = struct();
    
    for cycIdx = 1:length(rptCycles)
        cycle = rptCycles(cycIdx);
        cycle_key = sprintf('cyc%d', cycle);
        
        % CSV 파일 경로
        csvFile = fullfile(dataDir, sprintf('Ch%s_RPT_%dcyc.csv', ch, cycle));
        
        % CSV 파일 읽기 (파일이 없거나 읽을 수 없으면 건너뛰기)
        try
            if ~isfile(csvFile)
                fprintf('Warning: File not found - %s\n', csvFile);
                all_dQdV_data.(channel_key).(cycle_key).V = [];
                all_dQdV_data.(channel_key).(cycle_key).Q = [];
                all_dQdV_data.(channel_key).(cycle_key).V_dQdV = [];
                all_dQdV_data.(channel_key).(cycle_key).dQdV = [];
                continue;
            end
            
            T = readmatrix(csvFile);
        catch ME
            fprintf('Warning: Unable to read file - %s\n', csvFile);
            fprintf('  Error: %s\n', ME.message);
            all_dQdV_data.(channel_key).(cycle_key).V = [];
            all_dQdV_data.(channel_key).(cycle_key).Q = [];
            all_dQdV_data.(channel_key).(cycle_key).V_dQdV = [];
            all_dQdV_data.(channel_key).(cycle_key).dQdV = [];
            continue;
        end
        
        % Step Index (컬럼 2), Cycle Index (컬럼 4)
        % Discharge OCV 데이터 추출 (Step 10, CycleIdx 2)
        mask = (T(:,2) == discharge_step) & (T(:,4) == cycle_idx);
        
        if sum(mask) == 0
            fprintf('Warning: No discharge data found in Ch%s RPT%dcyc\n', ch, cycle);
            all_dQdV_data.(channel_key).(cycle_key).V = [];
            all_dQdV_data.(channel_key).(cycle_key).Q = [];
            all_dQdV_data.(channel_key).(cycle_key).V_dQdV = [];
            all_dQdV_data.(channel_key).(cycle_key).dQdV = [];
            continue;
        end
        
        % Voltage (컬럼 8), Capacity (컬럼 9)
        V_raw = T(mask, 8);  % Voltage [V]
        Q_raw = T(mask, 9);  % Capacity [Ah]
        
        % Sort by voltage (descending for discharge)
        % 전압 순서로 내림차순 정렬 (방전은 전압이 감소하므로)
        [V_sorted, sort_idx] = sort(V_raw, 'descend');
        Q_sorted = Q_raw(sort_idx);
        
        % Remove duplicate voltage values (보간을 위해)
        [V_unique, unique_idx, ~] = unique(V_sorted, 'stable');
        Q_unique = Q_sorted(unique_idx);
        
        % dQ/dV 계산 (FieldQmax_dQdV.m과 동일한 방식)
        if length(V_unique) > 1
            % Grid 보간법: 고정된 전압 간격(0.005V)으로 Grid 생성
            V_max = max(V_unique);
            V_min = min(V_unique);
            dV_grid = 0.005;   % 5mV (필드 데이터 노이즈 감소)
            V_grid = V_max:-dV_grid:V_min;  % Grid 전압 벡터 (높은 전압부터 낮은 전압으로)
            
            % 보간: 불규칙한 원본 데이터를 고른 Grid 위로 보간
            Q_grid = interp1(V_unique, Q_unique, V_grid, 'linear');
            
            % [핵심] Q 데이터를 부드럽게 다듬기 (전압 도메인 스무딩)
            Q_grid = smoothdata(Q_grid, 'gaussian', 15);
            
            % dQ/dV 계산: diff(Q_grid) / dV_grid (분모가 상수!)
            if length(Q_grid) > 1
                dQ = diff(Q_grid);  % ΔQ
                dV = dV_grid;       % 고정된 간격 0.005V
                
                % 중간 전압 값 사용 (V_mid = V_grid의 중간점)
                V_mid = (V_grid(1:end-1) + V_grid(2:end)) / 2;
                
                % dQ/dV 계산 (분모가 상수이므로 안정적)
                dQdV = dQ ./ dV;  % Ah/V (positive for discharge)
                
                all_dQdV_data.(channel_key).(cycle_key).V = V_raw;
                all_dQdV_data.(channel_key).(cycle_key).Q = Q_raw;
                all_dQdV_data.(channel_key).(cycle_key).V_dQdV = V_mid;
                all_dQdV_data.(channel_key).(cycle_key).dQdV = dQdV;
            else
                all_dQdV_data.(channel_key).(cycle_key).V = [];
                all_dQdV_data.(channel_key).(cycle_key).Q = [];
                all_dQdV_data.(channel_key).(cycle_key).V_dQdV = [];
                all_dQdV_data.(channel_key).(cycle_key).dQdV = [];
            end
        else
            all_dQdV_data.(channel_key).(cycle_key).V = [];
            all_dQdV_data.(channel_key).(cycle_key).Q = [];
            all_dQdV_data.(channel_key).(cycle_key).V_dQdV = [];
            all_dQdV_data.(channel_key).(cycle_key).dQdV = [];
        end
        
        fprintf('Processed: Ch%s RPT%dcyc - %d data points\n', ch, cycle, length(V_unique));
    end
end

%% ========================================================================
%  시각화: 개별 채널별 dQ/dV 곡선
% =========================================================================

fprintf('\n=== Creating Visualization ===\n');

% 개별 채널별 dQ/dV 곡선 (각 채널별 subplot)
fig = figure('Name', 'RPT Discharge dQ/dV - Individual Channels', ...
    'Position', [100 100 1600 1000], 'Color', 'w');

colors = lines(length(rptCycles));

% 2x4 서브플롯 (8채널)
for chIdx = 1:length(channels)
    ch = channels{chIdx};
    channel_key = sprintf('Ch%s', ch);
    
    subplot(2, 4, chIdx);
    hold on; grid on;
    
    for cycIdx = 1:length(rptCycles)
        cycle = rptCycles(cycIdx);
        cycle_key = sprintf('cyc%d', cycle);
        
        if isfield(all_dQdV_data, channel_key) && ...
           isfield(all_dQdV_data.(channel_key), cycle_key) && ...
           ~isempty(all_dQdV_data.(channel_key).(cycle_key).V_dQdV)
            
            V_plot = all_dQdV_data.(channel_key).(cycle_key).V_dQdV;
            dQdV_plot = all_dQdV_data.(channel_key).(cycle_key).dQdV * 1000;  % Ah/V -> mAh/V
            
            % NaN 제거
            valid_mask = ~isnan(V_plot) & ~isnan(dQdV_plot);
            V_plot = V_plot(valid_mask);
            dQdV_plot = dQdV_plot(valid_mask);
            
            if ~isempty(V_plot)
                plot(V_plot, dQdV_plot, '-', 'LineWidth', 1.5, ...
                    'Color', colors(cycIdx,:), 'DisplayName', sprintf('RPT%dcyc', cycle));
            end
        end
    end
    
    xlabel('Voltage [V]', 'FontSize', 10);
    ylabel('dQ/dV [mAh/V]', 'FontSize', 10);
    title(sprintf('Ch%s', ch), 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 8);
    xlim([2.5 4.2]);
    hold off;
end

% 그래프 저장
savefig(fig, fullfile(saveDir, 'RPT_Discharge_dQdV_Individual_Channels.fig'));
fprintf('Saved: RPT_Discharge_dQdV_Individual_Channels.fig\n');

%% ========================================================================
%  데이터 저장
% =========================================================================

matFile = fullfile(saveDir, 'RPT_Discharge_dQdV_data.mat');
save(matFile, 'all_dQdV_data', 'rptCycles', 'channels');
fprintf('\nSaved data to: %s\n', matFile);

fprintf('\n=== RPT Discharge dQ/dV Visualization Complete ===\n');

