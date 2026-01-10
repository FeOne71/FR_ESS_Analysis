%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 02_FeatureExtraction.m
% 목적: 저장된 Raw Event 파일에서 논문 기반의 SOH 예측용 Feature 추출
% 추가되는 피쳐:
%   1. 기본 통계: Mean, Var, Skew, Kurt (V, I)
%   2. 용량: Feat_dQ (Ah) - Delta Q (Partial Capacity)
%   3. 형태학적 피쳐: Voltage Slope, Entropy
%   4. 미분 피쳐: dV/dQ (Mean, Var, Skew, Kurt - 방전용), dQ/dV
%   5. DCIR: 1s, 5s, 10s, 30s
%
% 입력:
%   - Lab_DC_Events_Raw_*cyc.mat (01번 스크립트 출력)
%
% 출력:
%   - Lab_DC_Events_Features_*cyc.mat (피쳐가 추가된 이벤트 데이터)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

fprintf('=== Feature Extraction (Comprehensive) ===\n');

%% 설정
inputDir = fullfile(pwd, 'Results');
outputDir = fullfile(pwd, 'Results');

dt_list = [1, 5, 10, 30]; % DCIR 계산할 시간대

% 이벤트 타입 선택: 'Charge', 'Discharge', 'Both'
event_type_selection = 'Discharge';  % 'Charge': 충전만, 'Discharge': 방전만, 'Both': 모두

% 이벤트 타입 선택 검증
if ~ismember(event_type_selection, {'Charge', 'Discharge', 'Both'})
    error('event_type_selection must be ''Charge'', ''Discharge'', or ''Both''');
end

fprintf('Input directory: %s\n', inputDir);
fprintf('Output directory: %s\n', outputDir);
fprintf('Event type selection: %s\n', event_type_selection);

%% 파일 루프
matFiles = dir(fullfile(inputDir, 'Lab_DC_Events_Raw_*cyc.mat'));
if isempty(matFiles)
    error('No raw event files found in %s\nPlease run 01_Simple_EventDetection.m first', inputDir);
end

fprintf('Found %d raw event files\n', length(matFiles));

for i = 1:length(matFiles)
    fileName = matFiles(i).name;
    filePath = fullfile(inputDir, fileName);
    
    fprintf('\n--- Extracting features from %s ---\n', fileName);
    load(filePath, 'rawEvents'); 
    
    % 구조체 순회 (선택된 타입만 필터링)
    allChTypes = fieldnames(rawEvents);
    if strcmp(event_type_selection, 'Charge')
        chTypes = allChTypes(contains(allChTypes, '_Charge'));
    elseif strcmp(event_type_selection, 'Discharge')
        chTypes = allChTypes(contains(allChTypes, '_Discharge'));
    else % 'Both'
        chTypes = allChTypes;
    end
    
    total_features_added = 0;
    
    for c = 1:length(chTypes)
        ctName = chTypes{c};
        if ~isfield(rawEvents, ctName) || ~isstruct(rawEvents.(ctName))
            continue;
        end
        
        socs = fieldnames(rawEvents.(ctName));
        for s = 1:length(socs)
            sName = socs{s};
            if ~isfield(rawEvents.(ctName), sName) || ~isstruct(rawEvents.(ctName).(sName))
                continue;
            end
            
            profs = fieldnames(rawEvents.(ctName).(sName));
            for p = 1:length(profs)
                pName = profs{p};
                if ~isfield(rawEvents.(ctName).(sName), pName) || ~isstruct(rawEvents.(ctName).(sName).(pName))
                    continue;
                end
                
                evts = fieldnames(rawEvents.(ctName).(sName).(pName));
                
                for e = 1:length(evts)
                    evtName = evts{e};
                    if ~isfield(rawEvents.(ctName).(sName).(pName), evtName)
                        continue;
                    end
                    
                    D = rawEvents.(ctName).(sName).(pName).(evtName);
                    
                    if ~isfield(D, 't') || ~isfield(D, 'I') || ~isfield(D, 'V')
                        continue;
                    end
                    
                    % ----------------------------------------------------
                    % 1. 기본 통계 (V, I)
                    % ----------------------------------------------------
                    rawEvents.(ctName).(sName).(pName).(evtName).Feat_V_Mean = mean(D.V);
                    rawEvents.(ctName).(sName).(pName).(evtName).Feat_V_Var  = var(D.V);
                    rawEvents.(ctName).(sName).(pName).(evtName).Feat_V_Skew = skewness(D.V, 0); % bias correction
                    rawEvents.(ctName).(sName).(pName).(evtName).Feat_V_Kurt = kurtosis(D.V, 0); % bias correction
                    
                    rawEvents.(ctName).(sName).(pName).(evtName).Feat_I_Mean = mean(D.I);
                    rawEvents.(ctName).(sName).(pName).(evtName).Feat_I_Var  = var(D.I);
                    rawEvents.(ctName).(sName).(pName).(evtName).Feat_I_Skew = skewness(D.I, 0);
                    rawEvents.(ctName).(sName).(pName).(evtName).Feat_I_Kurt = kurtosis(D.I, 0);
                    
                    % ----------------------------------------------------
                    % 2. 용량 (Delta Q - Partial Capacity)
                    % ----------------------------------------------------
                    t_rel = D.t - D.t(1);
                    Q_profile = cumtrapz(t_rel, D.I) / 3600; % Ah
                    % 중요: 방전이면 Q가 감소하므로 절대값 사용
                    rawEvents.(ctName).(sName).(pName).(evtName).Feat_dQ = abs(Q_profile(end));
                    
                    % ----------------------------------------------------
                    % 3. Voltage Slope (Linear Fit)
                    % ----------------------------------------------------
                    if length(t_rel) > 1 && range(t_rel) > 0
                        P = polyfit(t_rel, D.V, 1);
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Slope_V = P(1); % dV/dt
                    else
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Slope_V = NaN;
                    end
                    
                    % ----------------------------------------------------
                    % 4. Entropy (Shannon Entropy of Voltage and Current)
                    % ----------------------------------------------------
                    % 데이터 분포(Histogram)의 엔트로피 계산
                    try
                        % Voltage Entropy
                        [counts_V, ~] = histcounts(D.V, 20); % 20 bins
                        probs_V = counts_V / sum(counts_V);
                        probs_V = probs_V(probs_V > 0); % log(0) 방지
                        if ~isempty(probs_V)
                            entropy_V = -sum(probs_V .* log2(probs_V));
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Entropy_V = entropy_V;
                        else
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Entropy_V = NaN;
                        end
                        
                        % Current Entropy
                        [counts_I, ~] = histcounts(D.I, 20);
                        probs_I = counts_I / sum(counts_I);
                        probs_I = probs_I(probs_I > 0);
                        if ~isempty(probs_I)
                            entropy_I = -sum(probs_I .* log2(probs_I));
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Entropy_I = entropy_I;
                        else
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Entropy_I = NaN;
                        end
                    catch
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Entropy_V = NaN;
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Entropy_I = NaN;
                    end
                    
                    % ----------------------------------------------------
                    % 5. 미분 피쳐 (dV/dQ, dQ/dV)
                    % ----------------------------------------------------
                    if length(D.V) > 10 && length(Q_profile) > 10
                        % 스무딩 (노이즈 감소)
                        win_size = max(3, floor(length(D.V)/10));
                        win_size = min(win_size, floor(length(D.V)/3)); % 너무 큰 윈도우 방지
                        V_smooth = smoothdata(D.V, 'gaussian', win_size);
                        Q_smooth = smoothdata(Q_profile, 'gaussian', win_size);
                        
                        dV = diff(V_smooth);
                        dQ = diff(Q_smooth);
                        
                        % --- dV/dQ (방전 데이터에서 강력함: 저항 유사 지표) ---
                        valid_dQ = abs(dQ) > 1e-6; % 0 나누기 방지
                        if sum(valid_dQ) > 5
                            dVdQ = dV(valid_dQ) ./ dQ(valid_dQ);
                            % Clipping (튀는 값 제거)
                            dVdQ(abs(dVdQ) > 100) = NaN; 
                            dVdQ_valid = dVdQ(~isnan(dVdQ));
                            
                            if length(dVdQ_valid) > 3
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Mean_dVdQ = nanmean(dVdQ);
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Var_dVdQ  = nanvar(dVdQ);
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Skew_dVdQ = skewness(dVdQ_valid, 0);
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Kurt_dVdQ = kurtosis(dVdQ_valid, 0);
                            else
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Mean_dVdQ = NaN;
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Var_dVdQ  = NaN;
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Skew_dVdQ = NaN;
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Kurt_dVdQ = NaN;
                            end
                        else
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Mean_dVdQ = NaN;
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Var_dVdQ  = NaN;
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Skew_dVdQ = NaN;
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Kurt_dVdQ = NaN;
                        end
                        
                        % --- dQ/dV (충전 데이터, 혹은 평탄하지 않은 구간에서 유효) ---
                        valid_dV = abs(dV) > 1e-6;
                        if sum(valid_dV) > 5
                            dQdV = dQ(valid_dV) ./ dV(valid_dV);
                            dQdV(abs(dQdV) > 1000) = NaN;
                            dQdV_valid = dQdV(~isnan(dQdV));
                            
                            if length(dQdV_valid) > 3
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Mean_dQdV = nanmean(dQdV);
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Var_dQdV  = nanvar(dQdV);
                            else
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Mean_dQdV = NaN;
                                rawEvents.(ctName).(sName).(pName).(evtName).Feat_Var_dQdV  = NaN;
                            end
                        else
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Mean_dQdV = NaN;
                            rawEvents.(ctName).(sName).(pName).(evtName).Feat_Var_dQdV  = NaN;
                        end
                    else
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Mean_dVdQ = NaN;
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Var_dVdQ  = NaN;
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Skew_dVdQ = NaN;
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Kurt_dVdQ = NaN;
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Mean_dQdV = NaN;
                        rawEvents.(ctName).(sName).(pName).(evtName).Feat_Var_dQdV  = NaN;
                    end
                    
                    % ----------------------------------------------------
                    % 6. DCIR (Dynamic Charge/Discharge Internal Resistance)
                    % ----------------------------------------------------
                    V1 = D.V(1); 
                    I1 = D.I(1);
                    for d_idx = 1:length(dt_list)
                        dt_target = dt_list(d_idx);
                        idx = find(t_rel >= dt_target, 1);
                        if ~isempty(idx) && idx <= length(D.V) && idx <= length(D.I)
                            dV_ir = D.V(idx) - V1;
                            dI_ir = D.I(idx) - I1;
                            if abs(dI_ir) > 0.1
                                R = abs(dV_ir / dI_ir) * 1000; % mOhm
                            else
                                R = NaN;
                            end
                        else
                            R = NaN;
                        end
                        fname = sprintf('DCIR_%ds', dt_target);
                        rawEvents.(ctName).(sName).(pName).(evtName).(fname).val = R;
                    end
                    
                    total_features_added = total_features_added + 1;
                end
            end
        end
    end
    
    % 저장 (변수명 변경: Lab_DC_DCIR_...)
    varNameNew = strrep(fileName, 'Lab_DC_Events_Raw_', 'Lab_DC_DCIR_');
    varNameNew = strrep(varNameNew, '.mat', '');
    varNameNew = regexprep(varNameNew, '[^a-zA-Z0-9_]', '_');
    
    eval(sprintf('%s = rawEvents;', varNameNew));
    
    newFileName = strrep(fileName, 'Raw', 'Features');
    savePath = fullfile(outputDir, newFileName);
    eval(sprintf('save(''%s'', ''%s'');', savePath, varNameNew));
    
    fprintf('  Processed %d events. Saved to %s\n', total_features_added, newFileName);
end

fprintf('\n=== Feature Extraction Complete ===\n');
fprintf('All results saved to: %s\n', outputDir);
