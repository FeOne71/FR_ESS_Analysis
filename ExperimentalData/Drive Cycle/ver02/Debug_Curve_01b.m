clear; clc; close all;

%% 사용자 설정
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';
targetChannel = 'Ch16'; % 변화가 가장 뚜렷한 가혹 조건 채널
targetSOC = 'SOC90';    % 상관성이 높았던 SOC 구간
targetProfile = 'DC1';  % 임의의 프로파일 선택

% 비교할 사이클 (초기 vs 후기)
cycles = {'0cyc', '600cyc'}; % 800cyc 파일이 없으면 있는 것 중 가장 큰 것

figure('Name', 'Curve Shape Analysis', 'Position', [100, 100, 1200, 800]);

for i = 1:length(cycles)
    cycName = cycles{i};
    filePath = fullfile(dataDir, sprintf('parsedDriveCycle_%s_filtered.mat', cycName));
    
    if ~exist(filePath, 'file')
        fprintf('파일 없음: %s\n', filePath); continue;
    end
    
    fprintf('Loading %s...\n', cycName);
    load(filePath);
    
    % 구조체 이름 찾기
    varName = sprintf('parsedDriveCycle_%s', cycName);
    eval(sprintf('data = %s;', varName));
    
    % 사용 가능한 채널 이름 확인
    availableChannels = fieldnames(data);
    fprintf('Available channels: %s\n', strjoin(availableChannels, ', '));
    
    % targetChannel이 사용 가능한 채널 중에 있는지 확인
    % 채널 이름 형식이 Ch16인지 ch16_Drive_0cyc_ChgEvent인지 확인
    matchingChannel = [];
    for ch = 1:length(availableChannels)
        chName = availableChannels{ch};
        % Ch16 또는 ch16으로 시작하는 채널 찾기
        if strcmpi(chName, targetChannel) || startsWith(chName, lower(targetChannel)) || ...
           startsWith(chName, targetChannel) || contains(chName, targetChannel)
            matchingChannel = chName;
            break;
        end
    end
    
    % 숫자로 매칭 시도 (Ch16 -> ch16 또는 16)
    if isempty(matchingChannel)
        chNum = regexp(targetChannel, '\d+', 'match');
        if ~isempty(chNum)
            chNumStr = chNum{1};
            for ch = 1:length(availableChannels)
                chName = availableChannels{ch};
                if contains(chName, chNumStr) && (startsWith(chName, 'ch') || startsWith(chName, 'Ch'))
                    matchingChannel = chName;
                    break;
                end
            end
        end
    end
    
    if isempty(matchingChannel)
        fprintf('ERROR: Channel %s not found in data!\n', targetChannel);
        fprintf('Available channels: %s\n', strjoin(availableChannels, ', '));
        continue;
    end
    
    fprintf('Using channel: %s\n', matchingChannel);
    
    % 데이터 추출
    try
        profData = data.(matchingChannel).(targetSOC).(targetProfile);
        V = profData.V;
        I = profData.I;
        t = profData.t;
        
        % 방전 구간 찾기 (전류 < 0)
        % 가장 긴 방전 구간 하나만 뽑아서 분석
        is_dchg = I < -1.0; % 1A 이상 방전
        dchg_starts = find(diff([0; is_dchg]) == 1);
        dchg_ends = find(diff([is_dchg; 0]) == -1);
        
        [max_len, max_idx] = max(dchg_ends - dchg_starts);
        if isempty(max_idx), continue; end
        
        idx_s = dchg_starts(max_idx);
        idx_e = dchg_ends(max_idx);
        
        V_seg = V(idx_s:idx_e);
        I_seg = I(idx_s:idx_e);
        t_seg = t(idx_s:idx_e);
        
        % --- [논문 방식 전처리] ---
        % 1. Smoothing
        smooth_window = max(5, round(length(V_seg)/10));
        V_filt = smoothdata(V_seg, 'gaussian', smooth_window);
        
        % 2. 용량(Q) 계산 (Ah)
        % 시간 데이터 타입 확인 및 변환
        if isa(t_seg, 'datetime')
            t_rel = seconds(t_seg - t_seg(1));
        elseif isa(t_seg, 'duration')
            t_rel = seconds(t_seg - t_seg(1));
        else
            % 이미 숫자형(초 단위)인 경우
            t_rel = t_seg - t_seg(1);
        end
        Q_seg = cumtrapz(t_rel, I_seg) / 3600; 
        Q_filt = smoothdata(Q_seg, 'gaussian', smooth_window);
        
        % 3. 미분 계산
        dQ = diff(Q_filt);
        dV = diff(V_filt);
        
        % dQ/dV (Clipping 적용)
        dQdV = dQ ./ dV;
        dQdV(abs(dQdV) > 50) = NaN; % 이상치 제거 후 시각화
        
        % dV/dQ (Clipping 적용)
        dVdQ = dV ./ dQ;
        dVdQ(abs(dVdQ) > 50) = NaN;
        
        % X축 (Q 또는 V)
        x_axis_V = V_filt(2:end);
        x_axis_Q = Q_filt(2:end);
        
        % --- [Plotting] ---
        color = 'b'; 
        if i==2, color = 'r'; end % 0cyc=Blue, 600cyc=Red
        
        % 1. V vs Q (기본 방전 곡선)
        subplot(2,2,1); hold on;
        plot(abs(Q_filt), V_filt, 'Color', color, 'LineWidth', 2);
        xlabel('Capacity Discharged (Ah)'); ylabel('Voltage (V)');
        title('Discharge Curve (V vs Q)');
        legend(cycles); grid on;
        
        % 2. dQ/dV Curve (문제의 구간)
        subplot(2,2,2); hold on;
        plot(x_axis_V, dQdV, 'Color', color, 'LineWidth', 1.5);
        xlabel('Voltage (V)'); ylabel('dQ/dV (Ah/V)');
        title('dQ/dV Curve (High Noise?)');
        ylim([-50, 50]); grid on;
        
        % 3. dV/dQ Curve (우리가 찾은 해답)
        subplot(2,2,3); hold on;
        plot(abs(x_axis_Q), dVdQ, 'Color', color, 'LineWidth', 1.5);
        xlabel('Capacity (Ah)'); ylabel('dV/dQ (V/Ah)');
        title('dV/dQ Curve (Resistance-like)');
        ylim([-5, 5]); grid on;
        
        % 4. Voltage Slope Histogram (통계 분포 확인)
        subplot(2,2,4); hold on;
        histogram(dVdQ, 20, 'FaceColor', color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        xlabel('dV/dQ Value'); ylabel('Count');
        title('Distribution of dV/dQ');
        legend(cycles); grid on;
        
    catch ME
        fprintf('Error in %s: %s\n', cycName, ME.message);
    end
end