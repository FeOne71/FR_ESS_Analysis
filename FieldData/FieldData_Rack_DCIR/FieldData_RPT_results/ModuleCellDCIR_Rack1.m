function res_all = ModuleCellDCIR_Rack1(module_file, rack_file, varargin)
% 목적: Rack#01 전류로 펄스(ΔI) 탐지 → Module#01의 개별 셀 전압으로 ΔV 산출 → DCIR(Ω) 계산 및 시각화
% 입력:
%   module_file : 'ModuleData.xlsx'
%   rack_file   : 'Rack2021_RPT.xlsx'
% 옵션(name-value):
%   'ModuleID'     (기본값 1)   % 모듈 번호
%   'DIminA'       (기본값 20)  % 펄스 탐지 최소 전류 스텝[A]
%   'MinSepSec'    (기본값 60)  % 펄스 간 최소 분리[초]
%   'TolModuleSec' (기본값 60)  % 모듈 최근접 허용 시간[초] (전/후 각각)
%   'T0List'       ([] 이면 자동탐지) % 수동 펄스 시각 목록: [datetime, datetime, ...]

p = inputParser;
addParameter(p,'ModuleID',1);
addParameter(p,'DIminA',20);
addParameter(p,'MinSepSec',60);
addParameter(p,'TolModuleSec',60);
addParameter(p,'T0List',datetime.empty);
parse(p,varargin{:});
ModuleID     = p.Results.ModuleID;
DI_abs_min_A = p.Results.DIminA;
min_sep_s    = p.Results.MinSepSec;
tol_module_s = p.Results.TolModuleSec;
t0_list      = p.Results.T0List;

%% 0) 데이터 로드 (원래 헤더 보존)
% --- Rack ---
TR = readtable(rack_file, 'VariableNamingRule','preserve');

% Time 컬럼 찾기 (예: 'Time')
timeCols = TR.Properties.VariableNames(contains(TR.Properties.VariableNames,'Time','IgnoreCase',true));
assert(~isempty(timeCols), 'Rack 파일에서 Time 컬럼을 찾지 못했습니다.');
% 문자열/엑셀 시리얼 모두 대응
try
    t_rack = datetime(TR.(timeCols{1}));
catch
    t_rack = datetime(TR.(timeCols{1}),'ConvertFrom','excel');
end

% 전류 컬럼 찾기 (예: 'DC Current(A)')
candI = TR.Properties.VariableNames(contains(TR.Properties.VariableNames,'Current','IgnoreCase',true) ...
    & contains(TR.Properties.VariableNames,'DC','IgnoreCase',true));
assert(~isempty(candI), '랙 전류 컬럼을 찾지 못했습니다. 예: ''DC Current(A)''');
I_rack = TR.(candI{1});

% --- Module ---
TM = readtable(module_file, 'VariableNamingRule','preserve');

% Time
assert(any(strcmpi(TM.Properties.VariableNames,'Time')), 'Module 파일의 Time 컬럼을 찾지 못했습니다.');
try
    t_mod  = datetime(TM.('Time'));
catch
    t_mod  = datetime(TM.('Time'),'ConvertFrom','excel');
end

% 모듈 번호 (파일에선 'Rack No.'로 표기됨)
if any(strcmpi(TM.Properties.VariableNames,'Rack No.'))
    mod_id = TM.('Rack No.');
elseif any(strcmpi(TM.Properties.VariableNames,'Rack No_'))
    mod_id = TM.('Rack No_');
else
    error('Module 파일에서 ''Rack No.'' 컬럼을 찾지 못했습니다.');
end

% 셀 전압 컬럼들 (예: 'Cell#01(V)' ~ 'Cell#14(V)')
cellCols = TM.Properties.VariableNames(startsWith(TM.Properties.VariableNames,'Cell#'));
assert(~isempty(cellCols), 'Module 파일에서 ''Cell#xx(V)'' 컬럼들을 찾지 못했습니다.');
nCells = numel(cellCols);

% 선택 모듈만 추출
TM1   = TM(mod_id == ModuleID, :);
t_mod1 = t_mod(mod_id == ModuleID);

%% 1) 펄스 시각 결정: 자동 탐지 또는 수동 지정
if isempty(t0_list)
    dI = [0; diff(I_rack)];
    idx = find(abs(dI) >= DI_abs_min_A);
    if isempty(idx)
        error('전류 step(>= %g A)을 찾지 못했습니다.', DI_abs_min_A);
    end
    % 너무 근접한 후보 제거
    keep = true(size(idx));
    for k = 2:numel(idx)
        if seconds(t_rack(idx(k)) - t_rack(idx(find(keep,1,'last')))) < min_sep_s
            keep(k) = false;
        end
    end
    idx = idx(keep);
    t0_list = t_rack(idx);
    dI_list = I_rack(idx+1) - I_rack(idx);
else
    % 수동 지정 시 ΔI 재계산
    dI_list = nan(size(t0_list));
    for k = 1:numel(t0_list)
        [~,ii] = min(abs(t_rack - t0_list(k)));
        assert(ii>=1 && ii<numel(I_rack), '수동 펄스 시각이 랙 데이터 범위를 벗어났습니다.');
        dI_list(k) = I_rack(ii+1) - I_rack(ii);
    end
end

%% 2) 개별 셀 DCIR 계산 (최근접, 보간/스무딩 없음)
Tol = seconds(tol_module_s);
res_all = table();
for k = 1:numel(t0_list)
    t0  = t0_list(k);
    dIk = dI_list(k);

    % 모듈 데이터에서 t0 직전/직후 최근접 샘플
    pre_idx  = find(t_mod1 <= t0, 1, 'last');
    post_idx = find(t_mod1 >= t0, 1, 'first');
    if isempty(pre_idx) || isempty(post_idx), continue; end
    if (t0 - t_mod1(pre_idx) > Tol) || (t_mod1(post_idx) - t0 > Tol), continue; end

    % 각 셀별 ΔV 및 R 계산
    dV    = zeros(1,nCells);
    Rcell = zeros(1,nCells);
    for c = 1:nCells
        Vpre  = TM1.(cellCols{c})(pre_idx);
        Vpost = TM1.(cellCols{c})(post_idx);
        dV(c)    = Vpost - Vpre;
        Rcell(c) = dV(c) / dIk;   % [Ohm]
    end

    % 결과 적재
    T = array2table(Rcell, 'VariableNames', regexprep(cellCols,'\(V\)','_Rohm'));
    T.t0       = t0;
    T.dI_A     = dIk;
    T.preTime  = t_mod1(pre_idx);
    T.postTime = t_mod1(post_idx);
    res_all = [res_all; T]; %#ok<AGROW>
end

%% 3) 시각화
% (a) 랙 전류 타임라인 + 펄스 시각
figure('Name','Rack Current & Pulses');
plot(t_rack, I_rack, 'LineWidth', 1); hold on; grid on;
for k = 1:height(res_all)
    xline(res_all.t0(k));
end
xlabel('Time'); ylabel('Rack DC Current [A]');
title('Rack#01 Current and Detected Pulses');

% (b) 마지막 펄스의 셀별 DCIR 막대그래프 (mΩ)
if ~isempty(res_all)
    lastRow = res_all(end,:);
    cellNames = regexprep(cellCols,'\(V\)',''); % ex) 'Cell#01'
    R_mOhm = 1000 * table2array(lastRow(:, regexprep(cellCols,'\(V\)','_Rohm')));
    figure('Name','Cell DCIR (last pulse)');
    bar(1:nCells, R_mOhm, 'LineWidth', 1);
    xticks(1:nCells); xticklabels(cellNames); xtickangle(0);
    ylabel('DCIR per Cell [m\Omega]');
    title(sprintf('Rack#01 - Module#%02d - DCIR per Cell (t0=%s)', ...
        ModuleID, datestr(lastRow.t0)));
    grid on;
end

% (c) 여러 펄스가 있으면 heatmap (셀 × 펄스)
if height(res_all) >= 2
    Rmat = zeros(height(res_all), nCells);
    for r = 1:height(res_all)
        Rmat(r,:) = table2array(res_all(r, regexprep(cellCols,'\(V\)','_Rohm')));
    end
    figure('Name','Cell DCIR Heatmap (all pulses)');
    imagesc(1000*Rmat); colorbar; % mΩ
    xlabel('Cell index'); ylabel('Pulse index (time order)');
    title(sprintf('Rack#01 - Module#%02d - DCIR Heatmap [m\\Omega]', ModuleID));
    xticks(1:nCells); xticklabels(strrep(cellNames,'Cell#',''));
end
end
