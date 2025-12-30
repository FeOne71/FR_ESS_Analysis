%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field_Operational_Impedance_Spectrum_NewData.m
%
% KEPCO FR ESS 필드 데이터로부터 SOH 및 열화 모드 분석을 위한
% 운용 데이터 기반의 고급 특징(HI)을 추출하는 스크립트.
%
% 추출 특징:
% 1. OIS (Operational Impedance Slope): 동적 임피던스 기울기
% 2. CDRAI (Charge-Discharge Resistance Asymmetry Index): 충/방전 저항 비대칭성
% 3. POR (Polarization-to-Ohmic Ratio): 펄스 전압 프로파일 형상 분석
%
% 최종 수정: 2025-10-02 (완성본)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% 1. 경로 및 파라미터 설정 (User Definition)
% =========================================================================
% --- 경로 설정 ---
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\New'; % 일별 .mat 파일이 있는 폴더
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_FFT','Operational_Impedance_Spectrum_OldData');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% --- 기본 파라미터 ---
Ns = 17*14; % 직렬 셀 개수 (랙 전압 계산용)
yearList = {'2023', '2024', '2025'}; % 분석할 연도
rackNames_all = {'Rack01'}; % 분석할 랙 이름
Cnom_cell = 64;

% --- OIS 파라미터 ---
OIS_SEGMENT_SEC = 900;           % 15분 (900초) 단위로 데이터 분할
OIS_ACTIVE_THR_A = 64*0.05;          % 세그먼트가 '활성'으로 판단될 전류 표준편차 임계값 (A)
OIS_F_HIGH_HZ = 0.1;             % OIS 계산용 고주파 (10초 주기)
OIS_F_LOW_HZ = 0.01;             % OIS 계산용 저주파 (100초 주기)

% --- CDRAI 파라미터 ---
CDRAI_SOC_BINS = (40:5:100)';     % 분석할 SOC 구간 (50-55%, 55-60%, ...)
CDRAI_REPRESENTATIVE_SOC = 65;   % 일별 대표값으로 사용할 중심 SOC (%)
CDRAI_MIN_POINTS_FOR_AVG = 10;   % 평균 저항 계산을 위한 최소 데이터 포인트 수

% --- PVPSA (POR) 파라미터 ---
POR_PULSE_MIN_DUR_SEC = 30;      % 안정적인 펄스로 판단할 최소 지속 시간 (초)
POR_PULSE_STABILITY_THR = 0.1;   % 전류 안정성 임계값 (변동폭 < 10%)
POR_OHMIC_DUR_SEC = 2;           % Ohmic 전압 강하 계산 시간 (초)

%% =========================================================================
total_results_table = table(); % 최종 결과를 저장할 테이블을 여기서 단 한 번만 생성!

if ~exist(saveDir, 'dir'), mkdir(saveDir); end

for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    yearPath = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));

    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);

        % 💡💡💡 [수정 2] monthly_table 초기화를 제거했습니다. 💡💡💡
        % 이제 total_results_table에 바로 누적합니다.

        matFiles = dir(fullfile(monthPath, '*.mat'));
        [~, idx] = sort({matFiles.name});
        matFiles = matFiles(idx);

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            fprintf('Processing: %s\n', matFiles(f).name);

            try
                dailyData = load(matFilePath);
            catch ME
                fprintf('  -> Error loading file: %s\n', ME.message);
                continue;
            end

            try
                % --- 데이터 준비 ---

                % 💡💡💡 [수정 1] 데이터 접근 경로 수정 💡💡💡
                % dailyData 구조체 안에 있는 'Raw' 구조체에 접근
                if isfield(dailyData, 'Raw')
                    rackData = dailyData.Raw;
                else
                    % 만약 'Raw' 필드가 없는 경우를 대비 (이전 데이터 형식 호환)
                    rackData = dailyData;
                end

                [~, name, ~] = fileparts(matFiles(f).name);
                date_match = regexp(name, '\d{8}', 'match');
                if ~isempty(date_match)
                    file_date = datetime(date_match{1}, 'InputFormat', 'yyyyMMdd');
                else
                    fprintf('  -> Could not parse date from filename: %s\n', name);
                    continue;
                end

                % 시간 변수가 duration 객체일 경우와 char일 경우 모두 처리
                if isduration(rackData.Date_Time)
                    time_of_day = rackData.Date_Time;
                else
                    time_of_day = duration(rackData.Date_Time, 'InputFormat', 'hh:mm:ss');
                end
                t_vec = file_date + time_of_day;

                % 💡💡💡 [수정 2] 'rackData' 변수를 사용하여 변수 할당 💡💡💡
                I = rackData.DCCurrent;
                V_cell_avg = rackData.CVavg;
                SOC = rackData.SOC_BMS;
                P = rackData.DCPower;
                T_batt = rackData.MTavg;

                % NaN 데이터 처리
                valid_idx = ~isnan(I) & ~isnan(V_cell_avg) & ~isnan(SOC) & ~isnat(t_vec) & ~isnan(P) & ~isnan(T_batt);

                t_vec = t_vec(valid_idx);
                I = I(valid_idx);
                V_cell_avg = V_cell_avg(valid_idx);
                SOC = SOC(valid_idx);
                P = P(valid_idx);
                T_batt = T_batt(valid_idx);

                if length(t_vec) < OIS_SEGMENT_SEC
                    fprintf('  -> Not enough valid data.\n');
                    continue;
                end

                fs = 1/mean(seconds(diff(t_vec)));
                if isnan(fs) || isinf(fs), fs = 1; end

                active_idx_temp = abs(I) > (Cnom_cell * 0.05);
                avg_temp = mean(T_batt(active_idx_temp));

                % --- 특징 추출 함수 호출 ---
                ois_val = calculate_ois(t_vec, I, V_cell_avg, fs, OIS_SEGMENT_SEC, OIS_ACTIVE_THR_A, OIS_F_HIGH_HZ, OIS_F_LOW_HZ);
                cdrai_val = calculate_cdrai(I, V_cell_avg, SOC, Cnom_cell, CDRAI_MIN_POINTS_FOR_AVG);
                por_val = calculate_por(t_vec, P, V_cell_avg, I, fs, Cnom_cell);

            catch ME
                fprintf('  -> Feature extraction failed: %s\n', ME.message);
                ois_val = NaN; cdrai_val = NaN; por_val = NaN; avg_temp = NaN;
            end

            % --- 결과 테이블 저장 ---
            rackName = "Rack01";
            daily_result = table(file_date, rackName, ois_val, cdrai_val, por_val, avg_temp, ...
                'VariableNames', {'Date', 'RackName', 'OIS', 'CDRAI', 'POR', 'Avg_Temp'});

            % 💡💡💡 [수정 3] monthly_table 대신 total_results_table에 바로 누적 💡💡💡
            total_results_table = [total_results_table; daily_result];

            fprintf('  -> Results: OIS=%.4f, CDRAI=%.4f, POR=%.4f\n', ois_val, cdrai_val, por_val);

        end % file loop�

    end % file loop

    disp('All processing complete.');

    %% 3. 최종 결과 저장 및 전체 기간 시각화
    % =========================================================================
    % 💡💡💡 [수정 4] 모든 루프가 끝난 후, 최종적으로 한 번만 저장하고 시각화 💡💡💡

    if ~isempty(total_results_table)
        % 날짜순으로 최종 정렬
        total_results_table = sortrows(total_results_table, 'Date');

        % 최종 결과 파일 저장
        finalSaveName = ['Advanced_HI_Results_Total_', strjoin(yearList, '_'), '.mat'];
        save(fullfile(saveDir, finalSaveName), 'total_results_table');
        fprintf('\nSuccessfully saved all results to: %s\n', fullfile(saveDir, finalSaveName));

        % --- 전체 기간에 대한 시계열 플롯 생성 ---
        fig_timeseries = figure('Name', 'Total Time Series Analysis', 'Position', [100 100 1400 900]);
        sgtitle('Total Health Indicator (HI) Trend Analysis', 'FontSize', 16, 'FontWeight', 'bold');

        ax1 = subplot(4,1,1);
        plot(total_results_table.Date, total_results_table.OIS, '-o', 'MarkerSize', 3);
        title('OIS (Operational Impedance Slope)'); ylabel('Slope (Ohm/log(Hz))'); grid on;

        ax2 = subplot(4,1,2);
        plot(total_results_table.Date, total_results_table.CDRAI, '-s', 'MarkerSize', 3, 'Color', 'r');
        title('CDRAI (Charge-Discharge Resistance Asymmetry Index)'); ylabel('Ratio (R_{ch}/R_{dch})'); grid on;

        ax3 = subplot(4,1,3);
        plot(total_results_table.Date, total_results_table.POR, '-d', 'MarkerSize', 3, 'Color', 'g');
        title('POR (Polarization-to-Ohmic Ratio)'); ylabel('Ratio (|\DeltaV_{pol}/\DeltaV_{ohm}|)'); grid on;

        ax4 = subplot(4,1,4);
        plot(total_results_table.Date, total_results_table.Avg_Temp, '.-', 'Color', [0.5 0.5 0.5]);
        title('Average Operating Temperature'); ylabel('Temperature (°C)'); grid on;

        % 모든 서브플롯의 x축을 연결하여 확대/축소 시 같이 움직이도록 함
        linkaxes([ax1, ax2, ax3, ax4], 'x');

        % 최종 플롯 저장
        finalPlotName = ['Total_TimeSeries_Plot_', strjoin(yearList, '_'), '.fig'];
        savefig(fig_timeseries, fullfile(saveDir, finalPlotName));

        disp('Total trend plot has been saved.');
    else
        disp('No data was processed. No results to save or plot.');
    end
end

% #### **2. OIS 계산 함수: `calculate_ois.m`**

function daily_ois = calculate_ois(t_vec, I, V, fs, segment_sec, active_thr, f_high, f_low)
% [수정] 디버깅을 위한 시각화 코드 추가

    % 💡💡💡 디버깅 플래그 💡💡💡
    DEBUG_PLOT = true; % true로 바꾸면 특정 날짜 처리 시 그래프가 나타남
    DEBUG_DATE = '2023-01-15'; % 확인하고 싶은 날짜를 'YYYY-MM-DD' 형식으로 입력

    % --- 디버깅용 Figure 생성 ---
    if DEBUG_PLOT && contains(string(t_vec(1)), DEBUG_DATE)
        fig_debug = figure('Name', ['Debug OIS - ' DEBUG_DATE], 'Position', [50 50 1500 800]);
        plot_count = 1;
    end

    if isempty(t_vec)
        daily_ois = NaN;
        return;
    end
    
    num_points_per_segment = floor(segment_sec * fs);
    if num_points_per_segment < 10, daily_ois = NaN; return; end
    num_segments = floor(length(I) / num_points_per_segment);
    
    ois_values = [];
    
    for i = 1:num_segments
        start_idx = (i-1) * num_points_per_segment + 1;
        end_idx = i * num_points_per_segment;
        
        I_seg = I(start_idx:end_idx);
        V_seg = V(start_idx:end_idx);
        
        if std(I_seg) < active_thr
            continue;
        end
        
        nfft = 2^nextpow2(length(I_seg));
        I_f = fft(I_seg, nfft);
        V_f = fft(V_seg, nfft);
        F_vec = (0:nfft-1)*(fs/nfft);
        
        I_f(1) = 0;
        Z_op_f = V_f ./ I_f;
        Z_op_f(abs(I_f) < 1e-3) = NaN;
        mag_Z = abs(Z_op_f);
        
        [~, idx_high] = min(abs(F_vec - f_high));
        [~, idx_low] = min(abs(F_vec - f_low));
        
        mag_Z_high = mag_Z(idx_high);
        mag_Z_low = mag_Z(idx_low);
        
        if isnan(mag_Z_high) || isnan(mag_Z_low) || isinf(mag_Z_high) || isinf(mag_Z_low)
            continue;
        end
        
        current_ois = (mag_Z_low - mag_Z_high) / (log10(f_high) - log10(f_low));
        ois_values = [ois_values; current_ois];
        
        % --- 💡💡💡 디버깅 플롯 로직 💡💡💡 ---
        if DEBUG_PLOT && contains(string(t_vec(1)), DEBUG_DATE)
            figure(fig_debug); % 생성해둔 Figure 활성화
            
            % 최대 12개까지만 플롯
            if plot_count <= 12
                subplot(3, 4, plot_count);
                loglog(F_vec, mag_Z, 'b-'); % 로그 스케일로 임피던스 스펙트럼 플롯
                hold on;
                % f_high, f_low 지점 표시
                loglog(F_vec(idx_high), mag_Z_high, 'ro', 'MarkerFaceColor', 'r');
                loglog(F_vec(idx_low), mag_Z_low, 'go', 'MarkerFaceColor', 'g');
                title(sprintf('Segment %d (OIS=%.4f)', i, current_ois));
                xlabel('Frequency (Hz)');
                ylabel('|Z_{op}| (Ohm)');
                grid on;
                xlim([0.001, fs/2]); % x축 범위 조절
                hold off;
                plot_count = plot_count + 1;
            end
        end
    end
    
    if isempty(ois_values)
        daily_ois = NaN;
    else
        ois_values = ois_values(~isoutlier(ois_values));
        daily_ois = median(ois_values, 'omitnan');
    end
    
    % --- 디버깅용 Figure 제목 추가 ---
    if DEBUG_PLOT && contains(string(t_vec(1)), DEBUG_DATE) && exist('fig_debug', 'var')
        sgtitle(sprintf('Impedance Spectrums for each Active Segment on %s (Final OIS = %.4f)', DEBUG_DATE, daily_ois), 'FontSize', 14);
    end
end