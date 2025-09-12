clear; clc; close all

%% *_SOC.mat 폴더
data_folder = "G:\공유 드라이브\BSL_Data2\BYD_LFP\Data\Processed_data\4C40\cell2_DCIR";
% save_path   = data_folder;
save_path = "D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT";
files = dir(fullfile(data_folder, '*_SOC.mat'));

%% 파라미터
dI_min_A = 0.1;   % |dI|<이 값이면 스킵
markerList = {'o','s','^','d','v','>','<','p','h','x','+'}; % 파일별 마커

%% 누적 CSV 버퍼
soc_all=[]; R_all=[]; file_all={};
cumQ_all=[]; Q_all=[];

%% 피겨 준비
fig = figure('Position', [100 100 1600 900]);
hold on; grid on;
xlabel('SOC (%)');
ylabel('R (m\Omega)');
title('SOC vs. R(30s)');

for fi = 1:numel(files)
    fname = files(fi).name;
    S = load(fullfile(data_folder, fname), 'data');
    if ~isfield(S,'data'), continue; end
    data = S.data;

    soc_vals = []; R_vals = []; cumQ_vals = []; Q_vals = [];

    for k = 1:numel(data)
        % n1C step만
        if ~isfield(data(k),'n1C_flag') || data(k).n1C_flag ~= 1, continue; end
        if ~isfield(data(k),'I') || ~isfield(data(k),'V') || numel(data(k).I)<2, continue; end

        % R 계산
        V_start = data(k).V(1);
        V_end   = data(k).V(end);
        I_end   = data(k).I(end);

        dV = V_end - V_start;
        dI = I_end - 0;   % 시작 전류 0 가정
        if ~isfinite(dV) || ~isfinite(dI) || abs(dI) < dI_min_A, continue; end

        R_step = abs(dV/dI); % [Ohm]
        if isfield(data(k),'SOC_start') && ~isnan(data(k).SOC_start)
            soc_vals(end+1,1) = data(k).SOC_start;
            R_vals(end+1,1)   = R_step;
            
            % dQ_Ah 계산 (structure.m과 동일한 변수명)
            if isfield(data(k),'dQ_Ah')
                cumQ_vals(end+1,1) = data(k).dQ_Ah;  % dQ_Ah (각 step의 전하량 변화)
                Q_vals(end+1,1) = abs(data(k).dQ_Ah); % |dQ_Ah| (절댓값)
            else
                cumQ_vals(end+1,1) = NaN;
                Q_vals(end+1,1) = NaN;
            end
            
            % 누적 CSV
            soc_all(end+1,1)  = data(k).SOC_start;
            R_all(end+1,1)    = R_step;
            cumQ_all(end+1,1) = cumQ_vals(end);
            Q_all(end+1,1)    = Q_vals(end);
            file_all{end+1,1} = fname;
        end
    end

    if ~isempty(soc_vals)
        [soc_sorted, idx] = sort(soc_vals, 'descend');
        R_sorted = R_vals(idx);
        cumQ_sorted = cumQ_vals(idx);
        Q_sorted = Q_vals(idx);

        mk = markerList{mod(fi-1,numel(markerList))+1};
        plot(soc_sorted, 1e3*R_sorted, ['-' mk], ...
            'LineWidth',1.5,'MarkerSize',4, 'DisplayName', fname);
            
        % 디버깅 표 출력
        fprintf('\n=== %s ===\n', fname);
        fprintf('SOC[%%]    dQ_Ah[Ah]    |dQ_Ah|[Ah]    R30[mΩ]\n');
        for i = 1:length(soc_sorted)
            fprintf('%6.2f    %8.3f    %8.3f    %8.4f\n', ...
                soc_sorted(i), cumQ_sorted(i), Q_sorted(i), 1e3*R_sorted(i));
        end
    end
end

legend('Location','northeast','Interpreter','none'); % 범례 그림 안쪽에
grid on; box on;
hold off;
ylim([0 20]);
yticks(0:2:20);

% % 전체 데이터 요약 표 출력
% fprintf('\n=== 전체 데이터 요약 ===\n');
% fprintf('SOC[%%]    dQ_Ah[Ah]    |dQ_Ah|[Ah]    R30[mΩ]    File\n');
% for i = 1:length(soc_all)
%     fprintf('%6.2f    %8.3f    %8.3f    %8.4f    %s\n', ...
%         soc_all(i), cumQ_all(i), Q_all(i), 1e3*R_all(i), file_all{i});
% end