%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 03_EventVisualization.m
% 목적: 추출된 이벤트와 [모든 피쳐]를 시각적으로 검증
% 기능: 
%   1. V, I 프로파일 그래프
%   2. dV/dQ (저항성분) 그래프
%   3. 계산된 모든 피쳐 값(dQ, Slope, Var_dVdQ 등) 텍스트 출력
%
% 입력:
%   - Lab_DC_Events_Features_*cyc.mat (02번 스크립트 출력)
%
% 출력:
%   - figures/EventVisualization/*.fig (시각화 그래프)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

fprintf('=== Event Visualization (Comprehensive) ===\n');

%% 설정
inputDir = fullfile(pwd, 'Results');
outputDir = fullfile(pwd, 'Results');
figuresDir = fullfile(outputDir, 'figures', 'EventVisualization');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

filePattern = 'Lab_DC_Events_Features_*cyc.mat';

% 이벤트 타입 선택: 'Charge', 'Discharge', 'Both'
event_type_selection = 'Discharge';  % 'Charge': 충전만, 'Discharge': 방전만, 'Both': 모두

% 이벤트 타입 선택 검증
if ~ismember(event_type_selection, {'Charge', 'Discharge', 'Both'})
    error('event_type_selection must be ''Charge'', ''Discharge'', or ''Both''');
end

fprintf('Input directory: %s\n', inputDir);
fprintf('Figures directory: %s\n', figuresDir);
fprintf('Event type selection: %s\n', event_type_selection);

%% 파일 찾기
matFiles = dir(fullfile(inputDir, filePattern));
if isempty(matFiles)
    error('No feature files found in %s\nPlease run 02_FeatureExtraction.m first', inputDir);
end

fprintf('Found %d feature files\n', length(matFiles));

%% 파일 하나 로드 (예: 첫 번째 파일)
fileIdx = 1; 
filePath = fullfile(inputDir, matFiles(fileIdx).name);
fprintf('\nLoading: %s\n', matFiles(fileIdx).name);
load(filePath);

% 변수명 자동 탐지
vars = who('-file', filePath);
dataVarName = '';
for v = 1:length(vars)
    if contains(vars{v}, 'Lab_DC_DCIR_')
        dataVarName = vars{v};
        break;
    end
end

if isempty(dataVarName)
    error('Could not find data variable in file');
end

data = eval(dataVarName);
fprintf('Loaded data structure: %s\n', dataVarName);

%% 랜덤 채널 및 이벤트 선택 (선택된 타입만 필터링)
allChFields = fieldnames(data);
fprintf('Available channels in data: %s\n', strjoin(allChFields, ', '));

if strcmp(event_type_selection, 'Charge')
    chFields = allChFields(contains(allChFields, '_Charge'));
elseif strcmp(event_type_selection, 'Discharge')
    chFields = allChFields(contains(allChFields, '_Discharge'));
else % 'Both'
    chFields = allChFields;
end

if isempty(chFields)
    % 사용 가능한 채널 타입 분석
    chargeChannels = allChFields(contains(allChFields, '_Charge'));
    dischargeChannels = allChFields(contains(allChFields, '_Discharge'));
    otherChannels = allChFields(~contains(allChFields, '_Charge') & ~contains(allChFields, '_Discharge'));
    
    errorMsg = sprintf(['No channels found for selected event type (%s) in data structure.\n' ...
                       'Available channel types:\n' ...
                       '  Charge channels: %d\n' ...
                       '  Discharge channels: %d\n' ...
                       '  Other channels: %d\n' ...
                       'Please check:\n' ...
                       '  1. If you ran 01_Simple_EventDetection.m with event_type_selection=''%s''\n' ...
                       '  2. If 02_FeatureExtraction.m processed the correct event type\n' ...
                       '  3. Or change event_type_selection to ''Charge'' or ''Both'' in this script'], ...
                       event_type_selection, ...
                       length(chargeChannels), ...
                       length(dischargeChannels), ...
                       length(otherChannels), ...
                       event_type_selection);
    error(errorMsg);
end

fprintf('Filtered channels for type ''%s'': %d channels\n', event_type_selection, length(chFields));

selCh = chFields{randi(length(chFields))}; % 랜덤 채널
fprintf('Selected channel: %s\n', selCh);

if ~isfield(data, selCh) || ~isstruct(data.(selCh))
    error('Invalid channel structure');
end

socFields = fieldnames(data.(selCh));
if isempty(socFields)
    error('No SOC levels found');
end

selSoc = socFields{1}; % 첫번째 SOC
fprintf('Selected SOC: %s\n', selSoc);

if ~isfield(data.(selCh), selSoc) || ~isstruct(data.(selCh).(selSoc))
    error('Invalid SOC structure');
end

profFields = fieldnames(data.(selCh).(selSoc));
if isempty(profFields)
    error('No profiles found');
end

selProf = profFields{1}; % Nth 프로파일
fprintf('Selected profile: %s\n', selProf);

if ~isfield(data.(selCh).(selSoc), selProf) || ~isstruct(data.(selCh).(selSoc).(selProf))
    error('Invalid profile structure');
end

evts = fieldnames(data.(selCh).(selSoc).(selProf));
if isempty(evts)
    error('No events found');
end

selEvt = evts{randi(length(evts))}; % 랜덤 이벤트
fprintf('Selected event: %s\n', selEvt);

evtData = data.(selCh).(selSoc).(selProf).(selEvt);

fprintf('\nVisualizing: %s / %s / %s / %s\n', selCh, selSoc, selProf, selEvt);

%% 시각화
fig = figure('Name', 'Event Feature Check', 'Position', [100, 100, 1400, 900], 'Visible', 'on');

% 1. V & I Profile
subplot(2,2,1);
t_rel = evtData.t - evtData.t(1);
yyaxis left; 
plot(t_rel, evtData.V, 'b-', 'LineWidth', 1.5); 
ylabel('Voltage (V)', 'FontSize', 11, 'FontWeight', 'bold');
yyaxis right; 
plot(t_rel, evtData.I, 'r-', 'LineWidth', 1.5); 
ylabel('Current (A)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
title('Raw Profile (V & I)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('Voltage', 'Current', 'Location', 'best', 'FontSize', 9);

% 2. dV/dQ Curve
subplot(2,2,2);
Q = cumtrapz(t_rel, evtData.I) / 3600; % Ah
% 직접 다시 계산하여 그리기 (검증용)
if length(evtData.V) > 10 && length(Q) > 10
    win_size = max(3, floor(length(evtData.V)/10));
    win_size = min(win_size, floor(length(evtData.V)/3));
    V_s = smoothdata(evtData.V, 'gaussian', win_size);
    Q_s = smoothdata(Q, 'gaussian', win_size);
    
    dVdQ = diff(V_s) ./ diff(Q_s);
    dVdQ(abs(dVdQ) > 50) = NaN; % Clipping
    Q_valid = Q_s(2:end);
    
    plot(abs(Q_valid), dVdQ, 'k-', 'LineWidth', 1.5);
    xlabel('Capacity (Ah)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('dV/dQ (V/Ah)', 'FontSize', 11, 'FontWeight', 'bold');
    title('dV/dQ Curve (Resistance Proxy)', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
else
    text(0.5, 0.5, 'Insufficient data for dV/dQ', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Units', 'normalized', 'FontSize', 12);
    title('dV/dQ Curve', 'FontSize', 12, 'FontWeight', 'bold');
end

% 3. Feature Summary (Text Box) - 모든 피쳐 값 표시
subplot(2,2,[3 4]); 
axis off;

% DCIR 텍스트 포맷팅
dcir_str = '';
dcir_fields = {'DCIR_1s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s'};
for d = 1:length(dcir_fields)
    if isfield(evtData, dcir_fields{d})
        dcir_val = evtData.(dcir_fields{d});
        if isstruct(dcir_val) && isfield(dcir_val, 'val') && ~isnan(dcir_val.val)
            dcir_str = [dcir_str, sprintf('  %s: %.2f mΩ\n', dcir_fields{d}, dcir_val.val)];
        end
    end
end
if isempty(dcir_str)
    dcir_str = '  (No DCIR data)\n';
end

% 기본 통계값 포맷팅
v_skew_str = 'N/A';
v_kurt_str = 'N/A';
i_skew_str = 'N/A';
i_kurt_str = 'N/A';

if isfield(evtData, 'Feat_V_Skew') && ~isnan(evtData.Feat_V_Skew)
    v_skew_str = sprintf('%.4f', evtData.Feat_V_Skew);
end
if isfield(evtData, 'Feat_V_Kurt') && ~isnan(evtData.Feat_V_Kurt)
    v_kurt_str = sprintf('%.4f', evtData.Feat_V_Kurt);
end
if isfield(evtData, 'Feat_I_Skew') && ~isnan(evtData.Feat_I_Skew)
    i_skew_str = sprintf('%.4f', evtData.Feat_I_Skew);
end
if isfield(evtData, 'Feat_I_Kurt') && ~isnan(evtData.Feat_I_Kurt)
    i_kurt_str = sprintf('%.4f', evtData.Feat_I_Kurt);
end

% dV/dQ 통계값 포맷팅
dvddq_skew_str = 'N/A';
dvddq_kurt_str = 'N/A';
if isfield(evtData, 'Feat_Skew_dVdQ') && ~isnan(evtData.Feat_Skew_dVdQ)
    dvddq_skew_str = sprintf('%.4f', evtData.Feat_Skew_dVdQ);
end
if isfield(evtData, 'Feat_Kurt_dVdQ') && ~isnan(evtData.Feat_Kurt_dVdQ)
    dvddq_kurt_str = sprintf('%.4f', evtData.Feat_Kurt_dVdQ);
end

% Entropy 포맷팅
entropy_v_str = 'N/A';
entropy_i_str = 'N/A';
if isfield(evtData, 'Feat_Entropy_V') && ~isnan(evtData.Feat_Entropy_V)
    entropy_v_str = sprintf('%.4f', evtData.Feat_Entropy_V);
end
if isfield(evtData, 'Feat_Entropy_I') && ~isnan(evtData.Feat_Entropy_I)
    entropy_i_str = sprintf('%.4f', evtData.Feat_Entropy_I);
end

infoStr = sprintf([...
    '=== Feature Summary (Paper-based) ===\n\n' ...
    '1. Basic Statistics:\n' ...
    '   Voltage:  Mean=%.4f V,  Var=%.4f,  Skew=%s,  Kurt=%s\n' ...
    '   Current:  Mean=%.4f A,  Var=%.4f,  Skew=%s,  Kurt=%s\n' ...
    '   Duration: %.2f s\n\n' ...
    '2. Capacity & Shape Features:\n' ...
    '   Delta Q (Partial Capacity): %.6f Ah\n' ...
    '   Voltage Slope (dV/dt): %.4e V/s\n' ...
    '   Voltage Entropy (Shannon): %s\n' ...
    '   Current Entropy (Shannon): %s\n\n' ...
    '3. Differential Analysis (dV/dQ):\n' ...
    '   Mean dV/dQ: %.4f V/Ah\n' ...
    '   Var dV/dQ:  %.4f (V/Ah)²\n' ...
    '   Skew dV/dQ: %s\n' ...
    '   Kurt dV/dQ: %s\n\n' ...
    '4. Differential Analysis (dQ/dV):\n' ...
    '   Mean dQ/dV: %.4f Ah/V\n' ...
    '   Var dQ/dV:  %.4f (Ah/V)²\n\n' ...
    '5. DCIR (Dynamic Internal Resistance):\n%s'], ...
    evtData.Feat_V_Mean, evtData.Feat_V_Var, v_skew_str, v_kurt_str, ...
    evtData.Feat_I_Mean, evtData.Feat_I_Var, i_skew_str, i_kurt_str, ...
    evtData.duration, ...
    evtData.Feat_dQ, ...
    evtData.Feat_Slope_V, ...
    entropy_v_str, entropy_i_str, ...
    evtData.Feat_Mean_dVdQ, evtData.Feat_Var_dVdQ, ...
    dvddq_skew_str, dvddq_kurt_str, ...
    evtData.Feat_Mean_dQdV, evtData.Feat_Var_dQdV, ...
    dcir_str);

text(0.05, 0.98, infoStr, 'FontSize', 10, 'FontName', 'FixedWidth', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
     'Interpreter', 'none');
title('Calculated Features (Paper-based)', 'FontSize', 14, 'FontWeight', 'bold');

% 전체 제목
sgtitle(sprintf('Event Visualization: %s / %s / %s / %s', ...
    selCh, selSoc, selProf, selEvt), 'FontSize', 14, 'FontWeight', 'bold');

% 저장
saveFileName = sprintf('Viz_%s_%s_%s_%s', selCh, selSoc, selProf, selEvt);
saveFileName = regexprep(saveFileName, '[^a-zA-Z0-9_\.]', '_'); % 파일명 안전하게 만들기
savePath_fig = fullfile(figuresDir, [saveFileName, '.fig']);
savePath_png = fullfile(figuresDir, [saveFileName, '.png']);

saveas(fig, savePath_fig);
saveas(fig, savePath_png);
fprintf('\nSaved visualization to:\n  %s\n  %s\n', savePath_fig, savePath_png);

fprintf('\n=== Visualization Complete ===\n');
