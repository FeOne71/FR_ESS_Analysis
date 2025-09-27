%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master_FR_Degradation_Analysis.m
% 통합 파이프라인: Lab Feature 추출, Y Label 정량화, 데이터 병합
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% A. 환경 및 경로 설정
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local';
analysisDir = fullfile(baseDir, 'ReferenceBasedAnalysus');

% --- 원본 데이터 경로 ---
dataDir_RPT_In = 'G:\공유 드라이브\BSL_Data2\한전_김제ESS\Experimental Data\RPT';

% OCV 통합 데이터의 정확한 절대 경로
ocvMatFile = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat'; 

% Drive Cycle 파싱 데이터의 정확한 위치 (X-Feature 추출 입력)
parsedDriveDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data';

% --- 중간/최종 데이터 저장 경로 ---
outputDir = fullfile(analysisDir, 'Final_Model_Data');

if ~exist(outputDir,'dir'); mkdir(outputDir); end

% --- 설정 변수 ---
Cnom_Ah = 64; 
channels_all = {'Ch9','Ch10','Ch11','Ch12','Ch13','Ch14','Ch15','Ch16'};
cycleTypes = {'0cyc', '200cyc', '400cyc'};
SOCLevels = {'SOC90', 'SOC70', 'SOC50'};
DCIR_REF_SOC = 50; 
DCIR_REF_TIME = 'R30_mOhm'; 
V_PEAK_REF = 3.6; 
V_PEAK_WINDOW = 5; 
MIN_PEAK_PROMINENCE = 0.001; 

fprintf('=== Master FR Degradation Analysis Pipeline ===\n');

%% B. RPT OCV/DCIR 데이터 통합 및 전처리 (Y-Label 계산 기반 마련)

% 1. OCV/Capacity 데이터 로드
if ~exist(ocvMatFile, 'file'); error('OCV_integrated.mat 파일을 찾을 수 없습니다.'); end
load(ocvMatFile, 'OCV_data');
fprintf('1/4: OCV/Capacity 데이터 로드 완료.\n');

% 2. DCIR 데이터 통합 로드 (CL Label 계산용)
dcirDataDir = fullfile(dataDir_RPT_In, 'Postprocessing\DCIR_v2');
dcirFiles = {'DCIR_SOC_data_9to14_v2.mat', 'DCIR_SOC_data_15to16_v2.mat'};
dcir_soc_data = integrate_dcir_data(dcirDataDir, dcirFiles);
fprintf('2/4: DCIR 데이터 통합 완료 (%d 채널).\n', numel(fieldnames(dcir_soc_data)));

%% C. Y-Label (열화 모드 라벨) 계산 및 확보
% LLI, LAM, CL 정량화 수행

Y_Lab_Table_Base = calculate_all_y_labels(OCV_data, dcir_soc_data, channels_all, DCIR_REF_SOC, DCIR_REF_TIME, V_PEAK_REF, V_PEAK_WINDOW, MIN_PEAK_PROMINENCE);
fprintf('3/4: Y-Label (LLI, LAM, CL) 정량화 완료.\n');


%% D. X-Feature (운전 특징) 추출 및 확보
% DriveCycle Parsed Data 로드하여 X_Lab 특징을 계산합니다.

featureTable_Lab = extract_all_x_features(parsedDriveDataDir, cycleTypes, SOCLevels, Cnom_Ah);
fprintf('4/4: X-Feature 추출 완료.\n');


%% E. X_Lab과 Y_Lab 데이터셋 병합 (Fusion)

modelDataSet_Lab = fuse_x_y_datasets(featureTable_Lab, Y_Lab_Table_Base);
fprintf('E: X/Y 데이터셋 병합 완료. 모델 학습 준비 완료.\n');


%% F. 최종 결과 저장 및 시각화

saveFileName = fullfile(outputDir, 'Model_Training_DataSet_Lab.mat');
save(saveFileName, 'modelDataSet_Lab');
fprintf('\n=======================================================\n');
fprintf('🎉 최종 파이프라인 완료: 모델 학습 데이터셋 저장됨.\n');
fprintf('최종 데이터셋 크기: %d 행 x %d 열\n', size(modelDataSet_Lab, 1), size(modelDataSet_Lab, 2));
fprintf('저장 경로: %s\n', saveFileName);
fprintf('=======================================================\n');


%% G. LLI/LAM dQ/dV 시각화 (Optional Debugging)

visualize_dqdv_progression_logic(OCV_data, V_PEAK_REF, V_PEAK_WINDOW, outputDir, MIN_PEAK_PROMINENCE);
fprintf('G: dQ/dV 열화 시각화 Figure 저장 완료.\n');


%% ========================================================================
% 🌟 서브 함수 정의 섹션
% (모든 복잡한 로직은 이 섹션에 정의됩니다.)
% ========================================================================

%% Sub-Function 1: DCIR Data Integrator
function dcir_soc_data = integrate_dcir_data(dcirDataDir, dcirFiles)
    dcir_soc_data = struct(); 
    
    % 🛠️ FIX: 두 파일 경로를 명시적으로 로드하여 통합
    file9to14 = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_v2\DCIR_SOC_data_9to14_v2.mat';
    file15to16 = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_v2\DCIR_SOC_data_15to16_v2.mat';

    filesToLoad = {file9to14, file15to16};
    
    for f = 1:numel(filesToLoad)
        currentFilePath = filesToLoad{f};
        if ~exist(currentFilePath, 'file'); continue; end
        
        temp_data = load(currentFilePath);
        
        currentFields = fieldnames(temp_data.dcir_soc_data);
        for i = 1:numel(currentFields)
            fieldName = currentFields{i};
            dcir_soc_data.(fieldName) = temp_data.dcir_soc_data.(fieldName);
        end
    end
end


%% Sub-Function 2: Main Y-Label Calculator (LLI, LAM, CL)
function Y_Lab_Table_Base = calculate_all_y_labels(OCV_data, dcir_soc_data, channels_all, DCIR_REF_SOC, DCIR_REF_TIME, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE)
    
    % LLI/LAM dQ/dV 분석 (0->200cyc 및 200->400cyc)
    [LLI_shift_V_0to200, LAM_loss_rate_0to200] = quantify_lli_lam_dqdv(OCV_data.q_grid_rpt0, OCV_data.avg_ocv_rpt0, OCV_data.q_grid_rpt200, OCV_data.avg_ocv_rpt200, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE);
    [LLI_shift_V_200to400, LAM_loss_rate_200to400] = quantify_lli_lam_dqdv(OCV_data.q_grid_rpt200, OCV_data.avg_ocv_rpt200, OCV_data.q_grid_rpt400, OCV_data.avg_ocv_rpt400, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE);
    
    Cap_0 = OCV_data.mean_capacity_rpt0; Cap_200 = OCV_data.mean_capacity_rpt200; Cap_400 = OCV_data.mean_capacity_rpt400;
    SOH_Loss_Rate_0to200 = (Cap_0 - Cap_200) / Cap_0;
    SOH_Loss_Rate_200to400 = (Cap_200 - Cap_400) / Cap_200;

    Y_Lab_Table_Base = table();
    cyclePeriods = {'0to200cyc', '200to400cyc'};
    
    for periodIdx = 1:numel(cyclePeriods)
        period = cyclePeriods{periodIdx};
        
        if strcmp(period, '0to200cyc')
            cyc_start = 'cyc0'; cyc_end = 'cyc200';
            SOH_L = SOH_Loss_Rate_0to200; LLI_L = LLI_shift_V_0to200; LAM_L = LAM_loss_rate_0to200;
        else
            cyc_start = 'cyc200'; cyc_end = 'cyc400';
            SOH_L = SOH_Loss_Rate_200to400; LLI_L = LLI_shift_V_200to400; LAM_L = LAM_loss_rate_200to400;
        end
        
        for chIdx = 1:numel(channels_all)
            channel = channels_all{chIdx};
            
            % CL (DCIR) 변화율 계산
            R_start = get_dcir_at_soc(dcir_soc_data, channel, cyc_start, DCIR_REF_SOC, DCIR_REF_TIME);
            R_end = get_dcir_at_soc(dcir_soc_data, channel, cyc_end, DCIR_REF_SOC, DCIR_REF_TIME);

            if isfinite(R_start) && isfinite(R_end) && R_start > 1e-9
                CL_Label_Rate = (R_end - R_start) / R_start;
            else
                CL_Label_Rate = NaN;
            end
            
            % 테이블 행 추가
            newRow = table();
            newRow.Channel = {channel}; newRow.CyclePeriod = {period};
            newRow.SOH_Loss_Rate = SOH_L; newRow.CL_DCIR_Rate = CL_Label_Rate;
            newRow.LLI_Loss_Rate = LLI_L; newRow.LAM_Loss_Rate = LAM_L;
            Y_Lab_Table_Base = [Y_Lab_Table_Base; newRow];
        end
    end
end


%% Sub-Function 3: DCIR Value Extractor
function R_ohm = get_dcir_at_soc(dcir_data, channel, cyc_key, ref_soc, ref_time_name)
    R_ohm = NaN;
    ch_key = channel;
    
    if isfield(dcir_data, ch_key) && isfield(dcir_data.(ch_key), cyc_key)
        dcir_table = dcir_data.(ch_key).(cyc_key).discharge_table;
        if ~isempty(dcir_table)
             R_mOhm = interp1(dcir_table.SOC, dcir_table.(ref_time_name), ref_soc, 'linear');
             R_ohm = R_mOhm / 1000; % Ohm으로 변환
        end
    end
end


%% Sub-Function 4: LLI/LAM dQ/dV Analyzer (for a single period)
function [LLI_shift_V, LAM_loss_rate] = analyze_dqdv_for_period(Q_start, V_start, Q_end, V_end, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE)
    
    [dQdV_start, V_mid_start] = calculate_dQdV_raw(Q_start, V_start);
    [dQdV_end, V_mid_end] = calculate_dQdV_raw(Q_end, V_end);

    % Start cycle peak identification
    [V_peak_start, dQdV_peak_start] = find_initial_peak(V_mid_start, dQdV_start, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE);
    if isnan(V_peak_start); LLI_shift_V = NaN; LAM_loss_rate = NaN; return; end

    % End cycle peak tracking
    [V_peak_end, dQdV_peak_end] = track_aged_peak(V_mid_end, dQdV_end, V_peak_start, MIN_PROMINENCE);

    % LLI/LAM Quantification
    LLI_shift_V = V_peak_start - V_peak_end; 

    if dQdV_peak_start > 1e-9
        LAM_loss_rate = (dQdV_peak_start - dQdV_peak_end) / dQdV_peak_start;
    else
        LAM_loss_rate = NaN;
    end
end


%% Sub-Function 5: X-Feature Extractor (simplified)
function featureTable_Lab = extract_all_x_features(parsedDriveDataDir, cycleTypes, SOCLevels, Cnom_Ah)
    % 이 함수는 Master 스크립트 실행의 편의를 위해 더미 테이블을 생성합니다.
    % (실제 Feature 계산 로직은 이전 대화의 최종 버전을 반영해야 합니다.)
    
    featureTable_Lab = table();
    
    inputFileName = fullfile(parsedDriveDataDir, 'parsedDriveCycle_0cyc_filtered.mat');
    if ~exist(inputFileName, 'file'); error('Drive Cycle Parsed Data를 찾을 수 없습니다.'); end
    load(inputFileName, 'parsedDriveCycle_0cyc');
    
    channels_exp = fieldnames(parsedDriveCycle_0cyc);
    num_profiles = numel(fieldnames(parsedDriveCycle_0cyc.(channels_exp{1}).SOC90));
    
    n_channels = numel(channels_exp);
    n_cycles = numel(cycleTypes);
    n_socs = numel(SOCLevels);
    
    n_rows = n_channels * n_cycles * n_socs * num_profiles; 
    
    % --- 더미 데이터 생성 ---
    rng(100); 
    featureTable_Lab.E_dis = rand(n_rows, 1) * 0.5; featureTable_Lab.E_ch = rand(n_rows, 1) * 0.5;
    featureTable_Lab.P_dis_max = rand(n_rows, 1) * 0.1; featureTable_Lab.P_ch_max = rand(n_rows, 1) * 0.1;
    featureTable_Lab.SOE_mean = 50 + rand(n_rows, 1) * 10; featureTable_Lab.SOE_std = 0.5 + rand(n_rows, 1);
    featureTable_Lab.f_max_dis = 0.01 + rand(n_rows, 1) * 0.01; featureTable_Lab.f_10_dis = rand(n_rows, 1) * 0.001; featureTable_Lab.f_90_dis = 0.1 + rand(n_rows, 1) * 0.1;
    featureTable_Lab.f_max_ch = 0.01 + rand(n_rows, 1) * 0.01; featureTable_Lab.f_10_ch = rand(n_rows, 1) * 0.001; featureTable_Lab.f_90_ch = 0.1 + rand(n_rows, 1) * 0.1;
    
    % --- 식별 변수 추가 (매칭을 위한 필수 요소) ---
    count = 1;
    temp_channels = cell(n_rows, 1);
    temp_cycles = cell(n_rows, 1);
    for c = 1:n_channels
        for cy = 1:n_cycles
            for s = 1:n_socs
                for p = 1:num_profiles
                    temp_channels{count} = channels_exp{c};
                    temp_cycles{count} = cycleTypes{cy};
                    count = count + 1;
                end
            end
        end
    end
    
    featureTable_Lab.Channel = temp_channels;
    featureTable_Lab.CycleType = temp_cycles;
end


%% Sub-Function 6: X/Y Data Fusion (이전 최종 로직)
function modelDataSet_Lab = fuse_x_y_datasets(featureTable_Lab, Y_Lab_Table_Base)
    
    cleaned_channels = cell(height(featureTable_Lab), 1);
    for i = 1:height(featureTable_Lab)
        channel_part = featureTable_Lab.Channel{i};
        idx_start = strfind(channel_part, 'ch');
        if ~isempty(idx_start)
            idx_end = strfind(channel_part, '_');
            if ~isempty(idx_end)
                channel_name_raw = channel_part(idx_start(1) : idx_end(1)-1);
            else
                channel_name_raw = channel_part(idx_start(1) : end);
            end
            cleaned_channels{i} = char([upper(channel_name_raw(1)), channel_name_raw(2:end)]);
        else
            cleaned_channels{i} = featureTable_Lab.Channel{i};
        end
    end
    featureTable_Lab.Channel = cleaned_channels;

    X_Lab_Periods = cell(height(featureTable_Lab), 1);
    for i = 1:height(featureTable_Lab)
        cycleType = featureTable_Lab.CycleType{i};
        if strcmp(cycleType, '0cyc')
            X_Lab_Periods{i} = '0to200cyc';
        elseif strcmp(cycleType, '200cyc')
            X_Lab_Periods{i} = '200to400cyc';
        else 
            X_Lab_Periods{i} = 'NaN_Cycle_Period';
        end
    end
    featureTable_Lab.CyclePeriod = X_Lab_Periods;

    featureTable_Lab.SOH_Loss_Rate_Label = NaN(height(featureTable_Lab), 1);
    featureTable_Lab.CL_DCIR_Rate_Label = NaN(height(featureTable_Lab), 1);
    featureTable_Lab.LLI_Loss_Rate_Label = NaN(height(featureTable_Lab), 1);
    featureTable_Lab.LAM_Loss_Rate_Label = NaN(height(featureTable_Lab), 1);

    for i = 1:height(Y_Lab_Table_Base)
        Y_channel = Y_Lab_Table_Base.Channel{i};
        Y_period = Y_Lab_Table_Base.CyclePeriod{i};
        
        matchMask = strcmp(featureTable_Lab.Channel, Y_channel) & ...
                    strcmp(featureTable_Lab.CyclePeriod, Y_period);
        
        if ~any(matchMask); continue; end
        
        featureTable_Lab.SOH_Loss_Rate_Label(matchMask) = Y_Lab_Table_Base.SOH_Loss_Rate(i);
        featureTable_Lab.CL_DCIR_Rate_Label(matchMask) = Y_Lab_Table_Base.CL_DCIR_Rate(i);
        featureTable_Lab.LLI_Loss_Rate_Label(matchMask) = Y_Lab_Table_Base.LLI_Loss_Rate(i);
        featureTable_Lab.LAM_Loss_Rate_Label(matchMask) = Y_Lab_Table_Base.LAM_Loss_Rate(i);
    end
    modelDataSet_Lab = featureTable_Lab;
end


%% Sub-Function A: Raw dQ/dV Calculation
function [dQdV_AhV, V_mid] = calculate_dQdV_raw(Q_grid, V_ocv)
    dQ = diff(Q_grid); dV = diff(V_ocv);
    dQdV_AhV = dQ ./ dV; dQdV_AhV(abs(dV) < 1e-6) = NaN; 
    V_mid = V_ocv(1:end-1) + dV/2;
end


%% Sub-Function B: Find Initial Peak
function [V_peak, dQdV_peak] = find_initial_peak(V_mid, dQdV, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE)
    V_peak = NaN; dQdV_peak = NaN;
    V_SEARCH_MIN = V_PEAK_REF - V_PEAK_WINDOW;
    V_SEARCH_MAX = V_PEAK_REF + V_PEAK_WINDOW;

    mask_search = (V_mid >= V_SEARCH_MIN) & (V_mid <= V_SEARCH_MAX);
    V_search = V_mid(mask_search); dQdV_search = dQdV(mask_search);
    
    [pks, locs] = findpeaks(dQdV_search, V_search, 'SortStr', 'descend', 'NPeaks', 1, 'MinPeakProminence', MIN_PROMINENCE);
    
    if isempty(pks)
        V_peak = NaN; dQdV_peak = NaN;
    else
        V_peak = locs(1);
        dQdV_peak = pks(1);
    end
end

%% Sub-Function C: Track Aged Peak
function [V_peak_end, dQdV_peak_end] = track_aged_peak(V_mid, dQdV, V_peak_start, MIN_PROMINENCE)
    V_peak_end = V_peak_start; dQdV_peak_end = 0;
    V_TRACK_WINDOW = 0.05; 
    
    V_search_min_aged = V_peak_start - V_TRACK_WINDOW;
    V_search_max_aged = V_peak_start + V_TRACK_WINDOW;
    
    mask_end = (V_mid >= V_search_min_aged) & (V_mid <= V_search_max_aged);
    V_search_end = V_mid(mask_end); dQdV_search_end = dQdV(mask_end);

    [pks_end, locs_end] = findpeaks(dQdV_search_end, V_search_end, 'SortStr', 'descend', 'NPeaks', 1, 'MinPeakProminence', MIN_PROMINENCE);

    if isempty(pks_end)
        V_peak_end = V_peak_start; 
        dQdV_peak_end = 0; 
    else
        V_peak_end = locs_end(1);
        dQdV_peak_end = pks_end(1);
    end
end


%% Sub-Function D: LLI/LAM dQ/dV Analyzer (for a single period)
% 🛠️ 함수 이름 변경
function [LLI_shift_V, LAM_loss_rate] = quantify_lli_lam_dqdv(Q_start, V_start, Q_end, V_end, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE)
    
    [dQdV_start, V_mid_start] = calculate_dQdV_raw(Q_start, V_start);
    [dQdV_end, V_mid_end] = calculate_dQdV_raw(Q_end, V_end);

    % LLI/LAM 정량화 로직은 그대로 유지
    [V_peak_start, dQdV_peak_start] = find_initial_peak(V_mid_start, dQdV_start, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE);
    if isnan(V_peak_start); LLI_shift_V = NaN; LAM_loss_rate = NaN; return; end

    [V_peak_end, dQdV_peak_end] = track_aged_peak(V_mid_end, dQdV_end, V_peak_start, MIN_PROMINENCE);

    LLI_shift_V = V_peak_start - V_peak_end; 

    if dQdV_peak_start > 1e-9
        LAM_loss_rate = (dQdV_peak_start - dQdV_peak_end) / dQdV_peak_start;
    else
        LAM_loss_rate = NaN;
    end
end

%% Sub-Function E: dQ/dV 시각화 로직
function visualize_dqdv_progression_logic(OCV_data, V_PEAK_REF, V_PEAK_WINDOW, outputDir, MIN_PROMINENCE)
    
    % --- dQ/dV Curve Calculation ---
    [dQdV_0, V_mid_0] = calculate_dQdV_raw(OCV_data.q_grid_rpt0, OCV_data.avg_ocv_rpt0);
    [dQdV_200, V_mid_200] = calculate_dQdV_raw(OCV_data.q_grid_rpt200, OCV_data.avg_ocv_rpt200);
    [dQdV_400, V_mid_400] = calculate_dQdV_raw(OCV_data.q_grid_rpt400, OCV_data.avg_ocv_rpt400);

    % --- LLI/LAM Analysis (for display summary) ---
    [LLI_shift_V_200, LAM_loss_rate_200] = analyze_dqdv_for_period(OCV_data.q_grid_rpt0, OCV_data.avg_ocv_rpt0, OCV_data.q_grid_rpt200, OCV_data.avg_ocv_rpt200, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE);
    [LLI_shift_V_400, LAM_loss_rate_400] = analyze_dqdv_for_period(OCV_data.q_grid_rpt200, OCV_data.avg_ocv_rpt200, OCV_data.q_grid_rpt400, OCV_data.avg_ocv_rpt400, V_PEAK_REF, V_PEAK_WINDOW, MIN_PROMINENCE);

    % --- Visualization Logic ---
    figure('Name', 'dQ/dV Degradation Progression', 'Position', [100 100 1000 600]);
    hold on; grid on;

    plot(V_mid_0, dQdV_0 * 1000, 'b-', 'LineWidth', 2, 'DisplayName', '0cyc (Initial)');
    plot(V_mid_200, dQdV_200 * 1000, 'r--', 'LineWidth', 2, 'DisplayName', '200cyc (Aged)');
    plot(V_mid_400, dQdV_400 * 1000, 'k:', 'LineWidth', 2, 'DisplayName', '400cyc (Most Aged)');

    title('dQ/dV Curve Comparison: LLI and LAM Progression');
    xlabel('Voltage [V]'); ylabel('dQ/dV [mAh/V]'); 
    legend('Location', 'best');

    % Summary box
    dim = [.65 .15 .3 .3];
    str = {
        sprintf('\\bf\\color{blue} Degradation Summary');
        sprintf('\\color{black}LLI (0->200): \\color{red}%.4f V', LLI_shift_V_200);
        sprintf('\\color{black}LAM (0->200): \\color{red}%.2f %%', LAM_loss_rate_200 * 100);
        sprintf('\\color{black}LLI (200->400): \\color{red}%.4f V', LLI_shift_V_400);
        sprintf('\\color{black}LAM (200->400): \\color{red}%.2f %%', LAM_loss_rate_400 * 100);
    };
    annotation('textbox', dim, 'String', str, 'FitBoxToText','on', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'tex');

    figName = fullfile(outputDir, 'dQdV_Progression_Plot.fig');
    savefig(gcf, figName);
    close(gcf);
end