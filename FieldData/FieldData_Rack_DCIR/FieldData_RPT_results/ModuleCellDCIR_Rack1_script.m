%% Rack#01 × Module#01 개별 셀 DCIR (2행 헤더, 3행 데이터)
rack_file   = 'Rack2021_RPT.xlsx';
module_file = 'ModuleData.xlsx';
ModuleID      = 1;     % <- 1번 모듈
DI_abs_min_A  = 20;    % 펄스 탐지 최소 전류 스텝[A]
min_sep_s     = 60;    % 펄스 간 최소 분리[초]
tol_module_s  = 180;   % 모듈 최근접 허용[초] (순차 로깅 감안, 120~180 권장)
t0_list       = [];    % 수동 지정 없으면 자동 탐지

% 0) Import: 2행=변수명, 3행=데이터
optsR = detectImportOptions(rack_file,   'VariableNamingRule','preserve');
optsR.VariableNamesRange = 'A2'; optsR.DataRange = 'A3';
TR = readtable(rack_file, optsR);

optsM = detectImportOptions(module_file, 'VariableNamingRule','preserve');
optsM.VariableNamesRange = 'A2'; optsM.DataRange = 'A3';
TM = readtable(module_file, optsM);

% 1) 필수 컬럼
t_rack = datetime(TR.('Time'));            % Rack Time
I_rack = TR.('DC Current(A)');             % Rack Current
t_mod_all = datetime(TM.('Time'));         % Module Time
mod_id_all = TM.('Rack No.');              % Module index (1~8)
cellCols = TM.Properties.VariableNames(startsWith(TM.Properties.VariableNames,'Cell#'));
nCells = numel(cellCols);

% 선택 모듈만 필터
TM1    = TM(mod_id_all == ModuleID, :);
t_mod1 = t_mod_all(mod_id_all == ModuleID);

% 2) 펄스 탐지(자동) 또는 수동
if isempty(t0_list)
    dI  = [0; diff(I_rack)];
    idx = find(abs(dI) >= DI_abs_min_A);
    keep = true(size(idx));
    for k = 2:numel(idx)
        lastKept = find(keep,1,'last');
        if seconds(t_rack(idx(k)) - t_rack(idx(lastKept))) < min_sep_s
            keep(k) = false;
        end
    end
    idx = idx(keep);
    t0_list = t_rack(idx);
    dI_list = I_rack(idx+1) - I_rack(idx);
else
    dI_list = nan(size(t0_list));
    for k = 1:numel(t0_list)
        [~,ii] = min(abs(t_rack - t0_list(k)));
        dI_list(k) = I_rack(ii+1) - I_rack(ii);
    end
end

% 3) 개별 셀 DCIR 계산(최근접, 무보간/무스무딩)
Tol = seconds(tol_module_s);
res_all = table();
for k = 1:numel(t0_list)
    t0  = t0_list(k);
    dIk = dI_list(k);
    pre_idx  = find(t_mod1 <= t0, 1, 'last');
    post_idx = find(t_mod1 >= t0, 1, 'first');
    if isempty(pre_idx) || isempty(post_idx), continue; end
    if (t0 - t_mod1(pre_idx) > Tol) || (t_mod1(post_idx) - t0 > Tol), continue; end

    R = zeros(1,nCells);
    for c = 1:nCells
        Vpre  = TM1.(cellCols{c})(pre_idx);
        Vpost = TM1.(cellCols{c})(post_idx);
        R(c)  = (Vpost - Vpre) / dIk;  % [Ohm]
    end
    T = array2table(R, 'VariableNames', regexprep(cellCols,'\(V\)','_Rohm'));
    T.t0 = t0; T.dI_A = dIk; T.preTime = t_mod1(pre_idx); T.postTime = t_mod1(post_idx);
    res_all = [res_all; T]; %#ok<AGROW>
end

% 4) 시각화
figure('Name','Rack Current & Pulses');
plot(t_rack, I_rack, 'LineWidth',1); grid on; hold on;
for k = 1:height(res_all), xline(res_all.t0(k)); end
xlabel('Time'); ylabel('Rack DC Current [A]');
title('Rack#01 Current and Detected Pulses');

if ~isempty(res_all)
    lastRow  = res_all(end,:);
    cellNames = regexprep(cellCols,'\(V\)','');
    R_mOhm   = 1000 * table2array(lastRow(:, regexprep(cellCols,'\(V\)','_Rohm')));
    figure('Name','Cell DCIR (last pulse)');
    bar(1:nCells, R_mOhm, 'LineWidth',1);
    xticks(1:nCells); xticklabels(cellNames); xtickangle(0);
    ylabel('DCIR per Cell [m\Omega]');
    title(sprintf('Rack#01 - Module#%02d - DCIR per Cell (t0=%s)', ModuleID, datestr(lastRow.t0)));
    grid on;
end

if height(res_all) >= 2
    Rmat = zeros(height(res_all), nCells);
    for r = 1:height(res_all)
        Rmat(r,:) = table2array(res_all(r, regexprep(cellCols,'\(V\)','_Rohm')));
    end
    figure('Name','Cell DCIR Heatmap (all pulses)');
    imagesc(1000*Rmat); colorbar;
    xlabel('Cell index'); ylabel('Pulse index'); 
    title(sprintf('Rack#01 - Module#%02d - DCIR Heatmap [m\\Omega]', ModuleID));
    xticks(1:nCells); xticklabels(strrep(cellNames,'Cell#',''));
end

disp(res_all);
