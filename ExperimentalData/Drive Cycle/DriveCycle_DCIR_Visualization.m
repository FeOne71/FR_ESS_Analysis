%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drive Cycle DCIR Visualization - All Channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Lab_DC_DCIR_200cyc_Events.mat');
dt_list = [1, 3, 5, 10, 30, 50];
dt_box = [1, 5, 10, 30]; % 사용할 시간 구간 (초)
soc_levels = {'SOC90', 'SOC70', 'SOC50'};
dc_profiles = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

% 저장 폴더 생성
if ~exist('figures/DC_Events', 'dir')
    mkdir('figures/DC_Events');
end

% 채널 목록 추출
channels = {};
fields = fieldnames(Lab_DC_DCIR_200cyc);
for i = 1:length(fields)
    if contains(fields{i}, '_ChgEvent')
        channel_name = strrep(fields{i}, '_Drive_200cyc_ChgEvent', '');
        channels{end+1} = channel_name;
    end
end
channels = unique(channels);

fprintf('Processing %d channels: ', length(channels));
for i = 1:length(channels)
    fprintf('%s ', channels{i});
end
fprintf('\n');

%% 모든 채널에 대한 시간-전류, 시간-전압 플롯 (SOC50만)
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    
    % 충전 이벤트: 시간-전류 플롯
    figure('Name', sprintf('%s Charging Events: Time vs Current (SOC50)', channel), ...
           'Position', [50 + (ch_idx-1)*50, 100 + (ch_idx-1)*50, 800, 600]);
    chg_struct_name = sprintf('%s_Drive_200cyc_ChgEvent', channel);
    soc_level = 'SOC50';
    
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        events = {};
        if isfield(Lab_DC_DCIR_200cyc, chg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(chg_struct_name), soc_level) && isfield(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level), dc_profile)
            events = fieldnames(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level).(dc_profile));
        end
        subplot(4, 2, dc_idx);
        hold on;
        for evt_idx = 1:length(events)
            evt_name = events{evt_idx};
            evt_data = Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level).(dc_profile).(evt_name);
            t_normalized = evt_data.t - evt_data.t(1);
            plot(t_normalized, evt_data.I, 'LineWidth', 1.5, 'DisplayName', sprintf('Event%d', evt_data.event_number));
        end
        xlabel('Time (s)'); ylabel('Current (A)'); 
        ylim([-5 20]);
        title(sprintf('%s-%s: Time vs Current (SOC50, Charging)', channel, dc_profile));
        grid on; legend('Location', 'best');
    end
    sgtitle(sprintf('%s Charging Events: Time vs Current by DC Profile (SOC50 only)', channel));
    
    % 파일로 저장
    saveas(gcf, sprintf('figures/DC_Events/%s_Charging_Time_Current_SOC50.fig', channel));
    close(gcf);
    
    % 충전 이벤트: 시간-전압 플롯
    figure('Name', sprintf('%s Charging Events: Time vs Voltage (SOC50)', channel), ...
           'Position', [1000 + (ch_idx-1)*50, 100 + (ch_idx-1)*50, 800, 600]);
    
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        events = {};
        if isfield(Lab_DC_DCIR_200cyc, chg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(chg_struct_name), soc_level) && isfield(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level), dc_profile)
            events = fieldnames(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level).(dc_profile));
        end
        subplot(4, 2, dc_idx);
        hold on;
        for evt_idx = 1:length(events)
            evt_name = events{evt_idx};
            evt_data = Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level).(dc_profile).(evt_name);
            t_normalized = evt_data.t - evt_data.t(1);
            plot(t_normalized, evt_data.V, 'LineWidth', 1.5, 'DisplayName', sprintf('Event%d', evt_data.event_number));
        end
        xlabel('Time (s)'); ylabel('Voltage (V)'); 
        ylim([3.0 4.2]);
        title(sprintf('%s-%s: Time vs Voltage (SOC50, Charging)', channel, dc_profile));
        grid on; legend('Location', 'best');
    end
    sgtitle(sprintf('%s Charging Events: Time vs Voltage by DC Profile (SOC50 only)', channel));
    
    % 파일로 저장
    saveas(gcf, sprintf('figures/DC_Events/%s_Charging_Time_Voltage_SOC50.fig', channel));
    close(gcf);
    
    % 방전 이벤트: 시간-전류 플롯
    figure('Name', sprintf('%s Discharging Events: Time vs Current (SOC50)', channel), ...
           'Position', [50 + (ch_idx-1)*50, 700 + (ch_idx-1)*50, 800, 600]);
    dchg_struct_name = sprintf('%s_Drive_200cyc_DchEvent', channel);
    
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        events = {};
        if isfield(Lab_DC_DCIR_200cyc, dchg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name), soc_level) && isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level), dc_profile)
            events = fieldnames(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level).(dc_profile));
        end
        subplot(4, 2, dc_idx);
        hold on;
        for evt_idx = 1:length(events)
            evt_name = events{evt_idx};
            evt_data = Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level).(dc_profile).(evt_name);
            t_normalized = evt_data.t - evt_data.t(1);
            plot(t_normalized, evt_data.I, 'LineWidth', 1.5, 'DisplayName', sprintf('Event%d', evt_data.event_number));
        end
        xlabel('Time (s)'); ylabel('Current (A)'); 
        ylim([-20 5]);
        title(sprintf('%s-%s: Time vs Current (SOC50, Discharging)', channel, dc_profile));
        grid on; legend('Location', 'best');
    end
    sgtitle(sprintf('%s Discharging Events: Time vs Current by DC Profile (SOC50 only)', channel));
    
    % 파일로 저장
    saveas(gcf, sprintf('figures/DC_Events/%s_Discharging_Time_Current_SOC50.fig', channel));
    close(gcf);
    
    % 방전 이벤트: 시간-전압 플롯
    figure('Name', sprintf('%s Discharging Events: Time vs Voltage (SOC50)', channel), ...
           'Position', [1000 + (ch_idx-1)*50, 700 + (ch_idx-1)*50, 800, 600]);
    
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        events = {};
        if isfield(Lab_DC_DCIR_200cyc, dchg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name), soc_level) && isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level), dc_profile)
            events = fieldnames(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level).(dc_profile));
        end
        subplot(4, 2, dc_idx);
        hold on;
        for evt_idx = 1:length(events)
            evt_name = events{evt_idx};
            evt_data = Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level).(dc_profile).(evt_name);
            t_normalized = evt_data.t - evt_data.t(1);
            plot(t_normalized, evt_data.V, 'LineWidth', 1.5, 'DisplayName', sprintf('Event%d', evt_data.event_number));
        end
        xlabel('Time (s)'); ylabel('Voltage (V)'); 
        ylim([3.0 4.2]);
        title(sprintf('%s-%s: Time vs Voltage (SOC50, Discharging)', channel, dc_profile));
        grid on; legend('Location', 'best');
    end
    sgtitle(sprintf('%s Discharging Events: Time vs Voltage by DC Profile (SOC50 only)', channel));
    
    % 파일로 저장
    saveas(gcf, sprintf('figures/DC_Events/%s_Discharging_Time_Voltage_SOC50.fig', channel));
    close(gcf);
    
    fprintf('Created and saved time-current and time-voltage plots for %s (SOC50)\n', channel);
end

%% 채널별 DCIR boxplot (충전 이벤트) - Visible off
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    
    figure('Name', sprintf('DCIR Boxplot by DC/SOC90/time - %s (Charging)', channel), ...
           'Position', [100 + (ch_idx-1)*50, 100 + (ch_idx-1)*50, 1200, 800], 'Visible', 'off');
    chg_struct_name = sprintf('%s_Drive_200cyc_ChgEvent', channel);
    soc_level = 'SOC90';
    
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        for dt_col = 1:length(dt_box) % 1s, 5s, 10s, 30s
            dt_sec = dt_box(dt_col);
            subplot_idx = (dc_idx-1)*4 + dt_col;
            subplot(8,4,subplot_idx);
            dcir_vals = [];
            if isfield(Lab_DC_DCIR_200cyc, chg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(chg_struct_name), soc_level) && isfield(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level), dc_profile)
                events = fieldnames(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level).(dc_profile).(evt_name);
                    field_name = sprintf('DCIR_%ds', dt_sec);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_vals = [dcir_vals; evt_data.(field_name).val];
                    end
                end
            end
            if ~isempty(dcir_vals)
                boxplot(dcir_vals, 'Labels', {sprintf('%s-%s - %ds', soc_level, dc_profile, dt_sec)});
            end
            title(sprintf('%s-%s - %ds DCIR', soc_level, dc_profile, dt_sec));
            ylabel('DCIR (mΩ)');
            set(gca, 'XTickLabel', []);
            grid on;
        end
    end
    sgtitle(sprintf('DCIR Boxplot by DC/SOC90/time - %s (Charging)', channel));
    
    % 파일로 저장
    saveas(gcf, sprintf('figures/DC_Events/%s_DCIR_Charging_Boxplot.fig', channel));
    close(gcf);
    fprintf('Saved DCIR charging boxplot for %s\n', channel);
end

%% 채널별 DCIR boxplot (방전 이벤트) - Visible off
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    
    figure('Name', sprintf('DCIR Boxplot by DC/SOC90/time - %s (Discharging)', channel), ...
           'Position', [200 + (ch_idx-1)*50, 100 + (ch_idx-1)*50, 1200, 800], 'Visible', 'off');
    dchg_struct_name = sprintf('%s_Drive_200cyc_DchEvent', channel);
    soc_level = 'SOC90';
    
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        for dt_col = 1:length(dt_box) % 1s, 5s, 10s, 30s
            dt_sec = dt_box(dt_col);
            subplot_idx = (dc_idx-1)*4 + dt_col;
            subplot(8,4,subplot_idx);
            dcir_vals = [];
            if isfield(Lab_DC_DCIR_200cyc, dchg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name), soc_level) && isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level), dc_profile)
                events = fieldnames(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level).(dc_profile));
                for evt_idx = 1:length(events)
                    evt_name = events{evt_idx};
                    evt_data = Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level).(dc_profile).(evt_name);
                    field_name = sprintf('DCIR_%ds', dt_sec);
                    if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                        dcir_vals = [dcir_vals; evt_data.(field_name).val];
                    end
                end
            end
            if ~isempty(dcir_vals)
                boxplot(dcir_vals, 'Labels', {sprintf('%s-%s - %ds', soc_level, dc_profile, dt_sec)});
            end
            title(sprintf('%s-%s - %ds DCIR', soc_level, dc_profile, dt_sec));
            ylabel('DCIR (mΩ)');
            set(gca, 'XTickLabel', []);
            grid on;
        end
    end
    sgtitle(sprintf('DCIR Boxplot by DC/SOC90/time - %s (Discharging)', channel));
    
    % 파일로 저장
    saveas(gcf, sprintf('figures/DC_Events/%s_DCIR_Discharging_Boxplot.fig', channel));
    close(gcf);
    fprintf('Saved DCIR discharging boxplot for %s\n', channel);
end

%% DCIR 계산값들을 하나의 figure에 표시 (모든 채널)
figure('Name', 'All Channels DCIR Comparison', 'Position', [100, 100, 1200, 800]);
soc_level = 'SOC90';

% DCIR 시간 구간별로 subplot 생성
for dt_idx = 1:length(dt_box)
    dt_sec = dt_box(dt_idx);
    subplot(2, 2, dt_idx);
    
    % 모든 채널의 DCIR 데이터 수집
    all_dcir_data = [];
    channel_labels = {};
    
    for ch_idx = 1:length(channels)
        channel = channels{ch_idx};
        chg_struct_name = sprintf('%s_Drive_200cyc_ChgEvent', channel);
        dchg_struct_name = sprintf('%s_Drive_200cyc_DchEvent', channel);
        
        % 충전 이벤트 DCIR 수집
        chg_dcir_vals = [];
        if isfield(Lab_DC_DCIR_200cyc, chg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(chg_struct_name), soc_level)
            for dc_idx = 1:length(dc_profiles)
                dc_profile = dc_profiles{dc_idx};
                if isfield(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level), dc_profile)
                    events = fieldnames(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level).(dc_profile));
                    for evt_idx = 1:length(events)
                        evt_name = events{evt_idx};
                        evt_data = Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level).(dc_profile).(evt_name);
                        field_name = sprintf('DCIR_%ds', dt_sec);
                        if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                            chg_dcir_vals = [chg_dcir_vals; evt_data.(field_name).val];
                        end
                    end
                end
            end
        end
        
        % 방전 이벤트 DCIR 수집
        dchg_dcir_vals = [];
        if isfield(Lab_DC_DCIR_200cyc, dchg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name), soc_level)
            for dc_idx = 1:length(dc_profiles)
                dc_profile = dc_profiles{dc_idx};
                if isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level), dc_profile)
                    events = fieldnames(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level).(dc_profile));
                    for evt_idx = 1:length(events)
                        evt_name = events{evt_idx};
                        evt_data = Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level).(dc_profile).(evt_name);
                        field_name = sprintf('DCIR_%ds', dt_sec);
                        if isfield(evt_data, field_name) && ~isnan(evt_data.(field_name).val)
                            dchg_dcir_vals = [dchg_dcir_vals; evt_data.(field_name).val];
                        end
                    end
                end
            end
        end
        
        % 모든 DCIR 값 결합
        all_vals = [chg_dcir_vals; dchg_dcir_vals];
        if ~isempty(all_vals)
            all_dcir_data = [all_dcir_data; all_vals];
            channel_labels = [channel_labels; repmat({channel}, length(all_vals), 1)];
        end
    end
    
    % Boxplot 생성
    if ~isempty(all_dcir_data)
        boxplot(all_dcir_data, channel_labels);
        title(sprintf('%ds DCIR - All Channels', dt_sec));
        ylabel('DCIR (mΩ)');
        xlabel('Channel');
        grid on;
    end
end

sgtitle('DCIR Comparison Across All Channels (SOC90)');

% 파일로 저장
saveas(gcf, 'figures/DC_Events/All_Channels_DCIR_Comparison.fig');

%% 모든 채널 요약 통계
fprintf('\n=== 모든 채널 이벤트 개수 요약 (SOC90) ===\n');
for ch_idx = 1:length(channels)
    channel = channels{ch_idx};
    chg_struct_name = sprintf('%s_Drive_200cyc_ChgEvent', channel);
    dchg_struct_name = sprintf('%s_Drive_200cyc_DchEvent', channel);
    soc_level = 'SOC90';
    
    fprintf('\n%s:\n', channel);
    fprintf('  충전: ');
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        n_evt = 0;
        if isfield(Lab_DC_DCIR_200cyc, chg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(chg_struct_name), soc_level) && isfield(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level), dc_profile)
            n_evt = length(fieldnames(Lab_DC_DCIR_200cyc.(chg_struct_name).(soc_level).(dc_profile)));
        end
        if dc_idx == 1
            fprintf('DC%d: %d', dc_idx, n_evt);
        else
            fprintf(', DC%d: %d', dc_idx, n_evt);
        end
    end
    fprintf('\n');
    
    fprintf('  방전: ');
    for dc_idx = 1:length(dc_profiles)
        dc_profile = dc_profiles{dc_idx};
        n_evt = 0;
        if isfield(Lab_DC_DCIR_200cyc, dchg_struct_name) && isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name), soc_level) && isfield(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level), dc_profile)
            n_evt = length(fieldnames(Lab_DC_DCIR_200cyc.(dchg_struct_name).(soc_level).(dc_profile)));
        end
        if dc_idx == 1
            fprintf('DC%d: %d', dc_idx, n_evt);
        else
            fprintf(', DC%d: %d', dc_idx, n_evt);
        end
    end
    fprintf('\n');
end

fprintf('\nVisualization completed for all channels\n');
fprintf('All figures saved to figures/DC_Events/\n');