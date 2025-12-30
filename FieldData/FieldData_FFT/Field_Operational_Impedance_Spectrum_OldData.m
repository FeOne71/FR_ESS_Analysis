%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field_Operational_Impedance_Spectrum.m
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
dataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old'; % 일별 .mat 파일이 있는 폴더
saveDir = fullfile('D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_FFT','Operational_Impedance_Spectrum_OldData');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% --- 기본 파라미터 ---
Ns = 17*14; % 직렬 셀 개수 (랙 전압 계산용)
yearList = {'2021', '2022', '2023'}; % 분석할 연도
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

% =========================================================================
%% 2. 메인 처리 루프
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

for year_idx = 1:length(yearList)
    year = yearList{year_idx};
    yearPath = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));
    
    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, '*.mat'));
        [~, idx] = sort({matFiles.name});
        matFiles = matFiles(idx);
        
        % 월별 테이블 초기화
        monthly_table = table();
        
        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            fprintf('Processing: %s\n', matFiles(f).name);
            
            try
                 load(matFilePath, 'Raw'); % 'Raw' 구조체 로드
            catch ME
                 fprintf('  -> Error loading file: %s\n', ME.message);
                 continue;
            end
            
            for rack_idx = 1:length(rackNames_all)
                rackName = rackNames_all{rack_idx};
                
                if ~isfield(Raw, rackName)
                    continue;
                end
                
                % --- 데이터 준비 ---
                rackData = Raw.(rackName);                
                I = rackData.DCCurrent_A;       % 전류 데이터 (A)
                V_cell_avg = rackData.AverageCV_V; % 평균 셀 전압 (V)
                SOC = rackData.SOCPct;            % SOC 데이터 (%)            
                P = rackData.DCPower_kW; % 전력 데이터 추가 로드
                T_batt = rackData.AverageMT_degC;
                t_vec = datetime(rackData.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
                % t_vec = rackData.Time;

                % NaN 데이터 처리
                valid_idx = ~isnan(I) & ~isnan(V_cell_avg) & ~isnan(SOC);
                t_vec = t_vec(valid_idx);
                I = I(valid_idx);
                V_cell_avg = V_cell_avg(valid_idx);
                SOC = SOC(valid_idx);
                P = P(valid_idx);
                T_batt = T_batt(valid_idx);

                % 💡 일별 평균 운전 온도 계산
                active_idx_temp = abs(I) > (Cnom_cell * 0.05);
                T_batt_avgTemp = mean(T_batt(active_idx_temp));

                if length(t_vec) < OIS_SEGMENT_SEC
                    fprintf('  -> Not enough valid data for rack %s.\n', rackName);
                    continue;
                end
                
                fs = 1/mean(seconds(diff(t_vec)));

                % --- 특징 추출 함수 호출 ---
                try
                    % 1. OIS 추출
                    ois_val = calculate_ois(t_vec, I, V_cell_avg, fs, OIS_SEGMENT_SEC, OIS_ACTIVE_THR_A, OIS_F_HIGH_HZ, OIS_F_LOW_HZ);
                    
                    % 2. CDRAI 추출
                    cdrai_val = calculate_cdrai(I, V_cell_avg, SOC, Cnom_cell, CDRAI_MIN_POINTS_FOR_AVG);

                    % 3. POR 추출
                    por_val = calculate_por(t_vec, P, V_cell_avg, I, fs, Cnom_cell);

                catch ME
                    fprintf('  -> Feature extraction failed for rack %s: %s\n', rackName, ME.message);
                    % 에러 발생 시 NaN 값으로 처리하여 테이블 구조 유지
                    ois_val = NaN; cdrai_val = NaN; por_val = NaN;
                end

                % --- 결과 테이블 저장 ---
                [~, name, ~] = fileparts(matFiles(f).name);
                DateStr = name(5:12);
                
                daily_result = table(datetime(DateStr, 'InputFormat', 'yyyyMMdd'), string(rackName), ois_val, cdrai_val, por_val, T_batt_avgTemp, ...
                    'VariableNames', {'Date', 'RackName', 'OIS', 'CDRAI', 'POR', 'Avg_Temp'});
                
                monthly_table = [monthly_table; daily_result];
                fprintf('  -> Rack %s: OIS=%.4f, CDRAI=%.4f, POR=%.4f\n', rackName, ois_val, cdrai_val, por_val);

            end % rack loop
        end % file loop
        
        % --- 월별 결과 저장 및 시각화 ---
        if ~isempty(monthly_table)
            month_name = monthDirs(m).name;
            saveFileName = fullfile(saveDir, ['Advanced_HI_', year, '_', month_name, '.mat']);
            save(saveFileName, 'monthly_table');
            fprintf('Saved monthly results to: %s\n\n', saveFileName);
            
            % 월별 시계열 플롯 생성
            fig_timeseries = figure('Name', ['Advanced HI Time Series: ', month_name], 'Position', [100 100 1200 800], 'Visible', 'off');
            sgtitle(['Advanced Health Indicators: ' month_name], 'FontSize', 14, 'FontWeight', 'bold');
            
            subplot(3,1,1); 
            plot(monthly_table.Date, monthly_table.OIS, '-o', 'MarkerSize', 4); 
            title('OIS (Operational Impedance Slope)'); ylabel('Slope (Ohm/log(Hz))'); grid on;
            
            subplot(3,1,2); 
            plot(monthly_table.Date, monthly_table.CDRAI, '-s', 'MarkerSize', 4); 
            title('CDRAI (Charge-Discharge Resistance Asymmetry Index)'); ylabel('Ratio (R_{ch}/R_{dch})'); grid on;
            
            subplot(3,1,3); 
            plot(monthly_table.Date, monthly_table.POR, '-d', 'MarkerSize', 4); 
            title('POR (Polarization-to-Ohmic Ratio)'); ylabel('Ratio (|\DeltaV_{pol}/\DeltaV_{ohm}|)'); grid on;
            
            % 시계열 플롯 저장
            saveas(fig_timeseries, fullfile(saveDir, ['TimeSeries_Plot_', year, '_', month_name, '.fig']));
            
            % 월별 박스 플롯 생성
            fig_boxplot = figure('Name', ['Advanced HI Boxplot: ', month_name], 'Position', [200 200 800 500], 'Visible', 'off');
            boxplot_data = [monthly_table.OIS, monthly_table.CDRAI, monthly_table.POR];
            boxplot(boxplot_data, 'Labels', {'OIS', 'CDRAI', 'POR'});
            title(['Monthly Distribution of Advanced Health Indicators: ' month_name]);
            ylabel('Value'); grid on;
            
            % 박스 플롯 저장
            saveas(fig_boxplot, fullfile(saveDir, ['BoxPlot_', year, '_', month_name, '.fig']));
            
            close(fig_timeseries);
            close(fig_boxplot);
        end % isempty(monthly_table)
    end % month loop
end % year loop

disp('All processing complete.');

% #### **2. OIS 계산 함수: `calculate_ois.m`**

function daily_ois = calculate_ois(t_vec, I, V, fs, segment_sec, active_thr, f_high, f_low)
% 하루 데이터로부터 대표 OIS 값을 계산합니다.

    if isempty(t_vec)
        daily_ois = NaN;
        return;
    end
    
    num_points_per_segment = round(segment_sec * fs);
    num_segments = floor(length(I) / num_points_per_segment);
    
    ois_values = [];
    
    for i = 1:num_segments
        start_idx = (i-1) * num_points_per_segment + 1;
        end_idx = i * num_points_per_segment;
        
        I_seg = I(start_idx:end_idx);
        V_seg = V(start_idx:end_idx);
        
        % 유효 구간 필터링
        if std(I_seg) < active_thr
            continue;
        end
        
        % FFT 수행
        nfft = 2^nextpow2(length(I_seg));
        I_f = fft(I_seg, nfft);
        V_f = fft(V_seg, nfft);
        F_vec = (0:nfft-1)*(fs/nfft);
        
        % 0 주파수(DC) 성분 제외
        I_f(1) = 0; 
        
        % 임피던스 계산
        Z_op_f = V_f ./ I_f;
        
        % 전류 크기가 작은 지점의 임피던스 값은 불안정하므로 제외
        Z_op_f(abs(I_f) < 1e-3) = NaN;
        
        mag_Z = abs(Z_op_f);
        
        % 대표 주파수에서 임피던스 크기 추출
        [~, idx_high] = min(abs(F_vec - f_high));
        [~, idx_low] = min(abs(F_vec - f_low));
        
        mag_Z_high = mag_Z(idx_high);
        mag_Z_low = mag_Z(idx_low);
        
        if isnan(mag_Z_high) || isnan(mag_Z_low) || isinf(mag_Z_high) || isinf(mag_Z_low)
            continue;
        end
        
        % OIS 계산
        current_ois = (mag_Z_low - mag_Z_high) / (log10(f_high) - log10(f_low));
        ois_values = [ois_values; current_ois];
    end
    
    if isempty(ois_values)
        daily_ois = NaN;
    else
        % 이상치 제거 후 중앙값 계산
        ois_values = ois_values(~isoutlier(ois_values));
        daily_ois = median(ois_values, 'omitnan');
    end
end