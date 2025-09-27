%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract_Lab_Y_Labels.m
% Y_Lab (SOH Loss, CL Loss, LLI, LAM) 라벨 계산 및 저장
% 🛠️ 200cyc to 400cyc 분석 로직 완벽 적용 (NaN 문제 해결)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% 경로 설정
rptParsedDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Parsed';
dcirDataDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\DCIR_v2';
ocvDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated';
outputFolder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Labels';

if ~exist(outputFolder,'dir'); mkdir(outputFolder); end

%% 설정 변수
cyclePeriods = {'0to200cyc', '200to400cyc'}; % 분석 기간
DCIR_REF_SOC = 50; % CL 라벨 추출을 위한 참조 SOC [%]
DCIR_REF_TIME = 'R30_mOhm'; % CL 라벨 추출을 위한 참조 저항 (30초 저항)
V_PEAK_REF = 3.6; % LLI/LAM 분석을 위한 피크 참조 전압 (셀 제조사 데이터시트 기반으로 조정 필요)
V_PEAK_WINDOW = 5; % 피크 탐색 윈도우 크기 (포인트)

% --- OCV 데이터 로드 ---
ocvMatFile = fullfile(ocvDataPath, 'OCV_integrated.mat');
if ~exist(ocvMatFile, 'file'); error('OCV_integrated.mat not found.'); end
load(ocvMatFile, 'OCV_data');
fprintf('OCV_integrated.mat 로드 완료.\n');

%% 1. 다중 DCIR MAT 파일 통합 로드
dcirFiles = {'DCIR_SOC_data_9to14_v2.mat', 'DCIR_SOC_data_15to16_v2.mat'};
dcir_soc_data = struct(); 

fprintf('다중 DCIR MAT 파일 통합 로드 시작...\n');
for f = 1:numel(dcirFiles)
    dcirMatFile = fullfile(dcirDataDir, dcirFiles{f});
    if ~exist(dcirMatFile, 'file'); fprintf('경고: %s 파일 부재. 스킵합니다.\n', dcirFiles{f}); continue; end
    
    temp_data = load(dcirMatFile);
    
    currentFields = fieldnames(temp_data.dcir_soc_data);
    for i = 1:numel(currentFields)
        fieldName = currentFields{i};
        dcir_soc_data.(fieldName) = temp_data.dcir_soc_data.(fieldName);
    end
    fprintf('%s 로드 및 통합 완료.\n', dcirFiles{f});
end

% 전체 채널 목록 재구성 (Ch9 ~ Ch16)
channels = fieldnames(dcir_soc_data);
Y_Lab_Table_Base = table();

%% 2. LLI/LAM 및 SOH Loss 계산 (주기별로 미리 계산)

% --- 2-1. 0cyc vs 200cyc 분석 ---
Q_grid_0 = OCV_data.q_grid_rpt0; V_ocv_0 = OCV_data.avg_ocv_rpt0; 
Q_grid_200 = OCV_data.q_grid_rpt200; V_ocv_200 = OCV_data.avg_ocv_rpt200; 

[dQdV_0, V_mid_0] = calculate_dQdV(Q_grid_0, V_ocv_0);
[dQdV_200, V_mid_200] = calculate_dQdV(Q_grid_200, V_ocv_200);

[LLI_shift_V_0to200, LAM_loss_rate_0to200] = analyze_dQdV_change(V_mid_0, dQdV_0, V_mid_200, dQdV_200, V_PEAK_REF, V_PEAK_WINDOW);

Cap_0 = OCV_data.mean_capacity_rpt0;
Cap_200 = OCV_data.mean_capacity_rpt200;
SOH_Loss_Rate_0to200 = (Cap_0 - Cap_200) / Cap_0;


% --- 2-2. 200cyc vs 400cyc 분석 ---
Q_grid_400 = OCV_data.q_grid_rpt400; V_ocv_400 = OCV_data.avg_ocv_rpt400; 

[dQdV_400, V_mid_400] = calculate_dQdV(Q_grid_400, V_ocv_400);

[LLI_shift_V_200to400, LAM_loss_rate_200to400] = analyze_dQdV_change(V_mid_200, dQdV_200, V_mid_400, dQdV_400, V_PEAK_REF, V_PEAK_WINDOW);

Cap_400 = OCV_data.mean_capacity_rpt400;
SOH_Loss_Rate_200to400 = (Cap_200 - Cap_400) / Cap_200; % 200cyc 용량 대비 400cyc 용량 손실


%% 3. 채널/주기별 데이터 테이블 구성 및 CL 계산
for chIdx = 1:numel(channels)
    channel = channels{chIdx};
    ch_key = channel;
    
    for periodIdx = 1:numel(cyclePeriods)
        period = cyclePeriods{periodIdx};
        
        % LLI/LAM/SOH 라벨 값 할당
        if strcmp(period, '0to200cyc')
            cyc_start = 'cyc0';
            cyc_end = 'cyc200';
            SOH_Loss_Rate_current = SOH_Loss_Rate_0to200;
            LLI_Loss_Rate_current = LLI_shift_V_0to200;
            LAM_Loss_Rate_current = LAM_loss_rate_0to200;
            
        elseif strcmp(period, '200to400cyc')
            cyc_start = 'cyc200';
            cyc_end = 'cyc400';
            SOH_Loss_Rate_current = SOH_Loss_Rate_200to400;
            LLI_Loss_Rate_current = LLI_shift_V_200to400; 
            LAM_Loss_Rate_current = LAM_loss_rate_200to400;
            
        else
            continue;
        end
        
        % DCIR 데이터 존재 여부 확인
        if ~isfield(dcir_soc_data, ch_key) || ...
           ~isfield(dcir_soc_data.(ch_key), cyc_start) || ...
           ~isfield(dcir_soc_data.(ch_key), cyc_end)
            % DCIR 데이터가 없으면 CL도 NaN으로 처리
            CL_Label_Rate = NaN; 
        else
            % CL Label (DCIR 증가율) 계산
            dcir_table_start = dcir_soc_data.(ch_key).(cyc_start).discharge_table;
            R_start = interp1(dcir_table_start.SOC, dcir_table_start.(DCIR_REF_TIME), DCIR_REF_SOC, 'linear');
            
            dcir_table_end = dcir_soc_data.(ch_key).(cyc_end).discharge_table;
            R_end = interp1(dcir_table_end.SOC, dcir_table_end.(DCIR_REF_TIME), DCIR_REF_SOC, 'linear');
            
            if isfinite(R_start) && isfinite(R_end) && R_start > 1e-6
                CL_Label_Rate = (R_end - R_start) / R_start;
            else
                CL_Label_Rate = NaN;
            end
        end

        %% 4. Y_Lab 테이블에 추가
        newRow = table();
        newRow.Channel = {channel};
        newRow.CyclePeriod = {period};
        newRow.SOH_Loss_Rate = SOH_Loss_Rate_current;
        newRow.CL_DCIR_Rate = CL_Label_Rate;
        newRow.LLI_Loss_Rate = LLI_Loss_Rate_current;
        newRow.LAM_Loss_Rate = LAM_Loss_Rate_current;
        
        Y_Lab_Table_Base = [Y_Lab_Table_Base; newRow];
    end
end


%% 5. Y_Lab 테이블 저장
saveFileName = fullfile(outputFolder, 'Lab_Degradation_Y_Labels_Final.mat');
save(saveFileName, 'Y_Lab_Table_Base');
fprintf('\n=== Y_Lab 테이블 추출 완료 (모든 주기 포함) ===\n');
fprintf('결과 파일: %s\n', saveFileName);


%% ========================================================================
% 🌟 서브 함수 1: dQ/dV 곡선 계산
% ========================================================================

function [dQdV_AhV, V_mid] = calculate_dQdV(Q_grid, V_ocv)
    % Q-V 데이터를 받아 dQ/dV 곡선을 계산
    dQ = diff(Q_grid);
    dV = diff(V_ocv);

    % dQ/dV 계산 (0으로 나누는 오류 방지)
    dQdV_AhV = dQ ./ dV;
    dQdV_AhV(abs(dV) < 1e-6) = NaN; % 전압 변화가 없으면 NaN 처리
    
    % V_mid 계산
    V_mid = V_ocv(1:end-1) + dV/2;
end

%% ========================================================================
% 🌟 서브 함수 2: LLI/LAM 정량화 (dQ/dV 변화 분석)
% ========================================================================

function [LLI_shift_V, LAM_loss_rate] = analyze_dQdV_change(V_mid_start, dQdV_start, V_mid_end, dQdV_end, V_PEAK_REF, V_PEAK_WINDOW)
    % 시작(start) 주기와 종료(end) 주기 dQ/dV 커브를 비교하여 LLI/LAM을 정량화합니다.
    
    % 1. 시작(start) 주기 피크 위치 찾기
    [~, idx_ref_start] = min(abs(V_mid_start - V_PEAK_REF));
    
    idx_start_win = max(1, idx_ref_start - V_PEAK_WINDOW);
    idx_end_win = min(length(V_mid_start), idx_ref_start + V_PEAK_WINDOW);
    
    V_window_start = V_mid_start(idx_start_win : idx_end_win);
    dQdV_window_start = dQdV_start(idx_start_win : idx_end_win);

    [dQdV_peak_start, peak_local_idx_start] = max(dQdV_window_start);
    V_peak_start = V_window_start(peak_local_idx_start);
    
    
    % 2. 종료(end) 주기 피크 위치 찾기 (start 피크 전압 근처에서 찾음)
    [~, idx_ref_end] = min(abs(V_mid_end - V_peak_start)); 
    
    idx_start_win_end = max(1, idx_ref_end - V_PEAK_WINDOW);
    idx_end_win_end = min(length(V_mid_end), idx_ref_end + V_PEAK_WINDOW);
    
    V_window_end = V_mid_end(idx_start_win_end : idx_end_win_end);
    dQdV_window_end = dQdV_end(idx_start_win_end : idx_end_win_end);

    [dQdV_peak_end, peak_local_idx_end] = max(dQdV_window_end);
    V_peak_end = V_window_end(peak_local_idx_end);

    % 3. LLI 정량화: 피크 전압의 수평 이동 (dV)
    % LLI는 피크 전압을 낮게 이동시킵니다 (V_peak_start > V_peak_end이면 LLI 발생)
    LLI_shift_V = V_peak_start - V_peak_end; 

    % 4. LAM 정량화: 피크 높이의 수직 감소 (dQ/dV 높이 변화)
    % LAM은 피크 높이를 줄입니다.
    if dQdV_peak_start > 1e-6
        LAM_loss_rate = (dQdV_peak_start - dQdV_peak_end) / dQdV_peak_start;
    else
        LAM_loss_rate = NaN;
    end
    
end