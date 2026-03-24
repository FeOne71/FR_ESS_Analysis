%% Field_CEff_Extractor.m
% =====================================================================
% Phase 4A: 필드 데이터 전력(CP) → 실효 전류율(C_eff) 변환 스크립트
% - 필드의 정전력(Constant Power) 제어 구간을 탐지
% - I_cell(t) = P_rack(t) / (V_rack(t) * Np) 로 변환
% - C_eff = mean(|I_cell|) / Q_nominal 로 산출하여 랩 CC 데이터와 동기화
% =====================================================================
clear; clc; close all;

%% 1. 파라미터 및 경로 설정
rawDir = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat';
files = {
    fullfile(rawDir,'Old','2021','202106','Raw_20210603.mat'), 'Y2021', 'old'
    fullfile(rawDir,'New','2023','202310','Raw_20231016.mat'), 'Y2023', 'new'
    fullfile(rawDir,'New','2024','202409','Raw_20240909.mat'), 'Y2024', 'new'
    fullfile(rawDir,'New','2025','202507','Raw_20250711.mat'), 'Y2025', 'new'
};

Np = 2; % 병렬 개수
Q_nominal = 64; % Ah (NCM 공칭 용량)
thr_A = Q_nominal * 0.05; % 노이즈 컷 (3.2A)

fprintf('=== Phase 4A: Field CP -> C_eff Conversion ===\n');

res_CEff = struct();

for k = 1:size(files,1)
    fpath = files{k,1}; yr = files{k,2}; dtype = files{k,3};
    if ~exist(fpath,'file'), fprintf('[Skip] %s not found\n', yr); continue; end
    
    S = load(fpath);
    
    %% 2. 데이터 형식에 따른 I, V, P 추출
    if strcmp(dtype,'old')
        I_rack = S.Raw.Rack01.DCCurrent_A(:);
        V_rack = S.Raw.Rack01.SumCV_V(:);       % Rack 전체 전압
        P_rack = S.Raw.Rack01.DCPower_kW(:) * 1000; % kW -> W
        time_arr = S.Raw.Rack01.Time;
        dt = datetime(time_arr, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        t_s = seconds(dt - dt(1));
    else
        I_rack = S.Raw.DCCurrent(:);
        V_rack = S.Raw.CVavg(:) * 14 * 14;      % CVavg(Cell) * 14 Cell/Module * 14 Module/Rack = Rack Voltage
        P_rack = S.Raw.DCPower(:) * 1000;       % kW -> W
        time_arr = S.Raw.Date_Time;
        t_s = seconds(time_arr) - seconds(time_arr(1));
    end
    
    if isempty(t_s), t_s = (0:length(I_rack)-1)'; end
    
    %% 3. 충전 구간 탐지 및 C_eff 변환 로직
    % 충전 전류(I > 0) 활성화 마스크
    chg_mask = I_rack/Np > thr_A;
    d = diff([0; chg_mask; 0]);
    starts = find(d==1); ends = find(d==-1)-1;
    lens = ends - starts + 1;
    
    if isempty(lens)
        fprintf('-- %s: No charging sequence found.\n', yr);
        continue;
    end
    
    % 가장 긴 충전 시퀀스 추출
    [max_len, max_i] = max(lens);
    idx = starts(max_i):ends(max_i);
    
    % 해당 시퀀스 내부 물리량 절단
    t_seg = t_s(idx);
    V_seg = V_rack(idx);
    P_seg = P_rack(idx);
    I_meas_seg = I_rack(idx) / Np; % 측정된 Cell 전류 (비교용)
    
    %% 4. 수식 기반 C_eff 변환 (Physical Mapping)
    % 4-1. Power에서 파생된 전류 I_derive(t) = P(t) / (V(t) * Np)
    % 여기서 P(t)가 충전일 때 양수/음수 표기법 혼재될 수 있으므로 절대값 취함
    I_derive_seg = abs(P_seg) ./ (V_seg * Np); 
    
    % 4-2. 적분을 통한 실효 전류 평균 (C_eff) 계산
    % C_eff = [ (1/T) * \int |I(t)| dt ] / Q_nominal
    T_total = t_seg(end) - t_seg(1);
    int_I_derived = trapz(t_seg, I_derive_seg);
    I_mean_derived = int_I_derived / T_total;
    C_eff_derived = I_mean_derived / Q_nominal;
    
    % 측정된 전류 기준 C_eff 계산 (Cross-check 용)
    int_I_meas = trapz(t_seg, I_meas_seg);
    I_mean_meas = int_I_meas / T_total;
    C_eff_meas = I_mean_meas / Q_nominal;
    
    % 5. 결과 저장 및 표출
    fprintf('-- %s --\n', yr);
    fprintf('  Seq Length : %d seconds\n', max_len);
    fprintf('  Mean Meas I: %.2f A -> C_eff (Meas): %.3f C\n', I_mean_meas, C_eff_meas);
    fprintf('  Mean Deri I: %.2f A -> C_eff (Deri): %.3f C\n', I_mean_derived, C_eff_derived);
    
    res_CEff.(yr).C_eff_Meas = C_eff_meas;
    res_CEff.(yr).C_eff_Deri = C_eff_derived;
    res_CEff.(yr).I_meas_seg = I_meas_seg;
    res_CEff.(yr).I_derive_seg = I_derive_seg;
    res_CEff.(yr).t_seg = t_seg;
end

%% 6. 시각화 (CP Profile vs CC Profile의 차이 증명)
% 대표적으로 가장 최신 연도(Y2025) 데이터의 전류 프로파일 굵게 그리기
if isfield(res_CEff, 'Y2025')
    yr_target = 'Y2025';
else
    yr_target = fieldnames(res_CEff);
    yr_target = yr_target{1};
end

t_plot = res_CEff.(yr_target).t_seg - res_CEff.(yr_target).t_seg(1);
I_m = res_CEff.(yr_target).I_meas_seg;
I_d = res_CEff.(yr_target).I_derive_seg;
avg_C = res_CEff.(yr_target).C_eff_Meas;

fig = figure('Position', [100, 200, 800, 400], 'Name', 'CP to C_eff Physical Mapping');
hold on; grid on;
plot(t_plot/60, I_m, 'b-', 'LineWidth', 2, 'DisplayName', 'Measured I_{cell}(t)');
% plot(t_plot/60, I_d, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Derived I_{cell}(t) from CP');
plot([t_plot(1)/60 t_plot(end)/60], [avg_C*Q_nominal avg_C*Q_nominal], 'k--', 'LineWidth', 2, 'DisplayName', sprintf('Equivalent C_{eff} (%.3fC = %.1fA)', avg_C, avg_C*Q_nominal));

title(sprintf('Field [%s] Constant Power (CP) Control Profile', yr_target), 'FontWeight', 'bold');
xlabel('Time (Minutes)');
ylabel('Current (A)');
legend('Location', 'best');
ylim([0 max(I_m)*1.2]);

saveDir = fileparts(mfilename('fullpath'));
saveas(fig, fullfile(saveDir, 'Field_CP_Profile.fig'));
fprintf('\n>> Saved: Field_CP_Profile.fig\n');
