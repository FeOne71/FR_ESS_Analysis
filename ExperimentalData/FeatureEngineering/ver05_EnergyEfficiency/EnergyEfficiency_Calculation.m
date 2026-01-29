%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EnergyEfficiency_Calculation.m
% - Calculate energy efficiency for Ch09 drive cycle data
% - Use OCV-based SOC correction formula
% - Generate Long Format table for correlation analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% Configuration
% File paths
ocvPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\OCV_integrated\OCV_integrated.mat';
driveCycleDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\parsed_data'; % Reference: Simple_EventDetection_01.m
outputDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\ver05_EnergyEfficiency\Results';
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% Target channel
targetChannel = 'Ch09';

% RPT cycles
rptCycles = {'0cyc', '200cyc', '400cyc', '600cyc', '800cyc', '1000cyc'};
rptCycleNums = [0, 200, 400, 600, 800, 1000];

% SOC levels
socLevels = {'SOC90', 'SOC70', 'SOC50'};

% Drive cycle types
dcTypes = {'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6', 'DC7', 'DC8'};

% Current threshold for E_net calculation (A)
I_threshold = 0.64;

%% Load OCV data
fprintf('\n=== Loading OCV Data ===\n');
if ~exist(ocvPath, 'file')
    error('OCV file not found: %s', ocvPath);
end

ocvData = load(ocvPath, 'OCV_data');
OCV_data = ocvData.OCV_data;

% Verify OCV data structure
fprintf('OCV data loaded successfully\n');
ocvFields = fieldnames(OCV_data);
fprintf('Available OCV fields: %s\n', strjoin(ocvFields(1:min(10, length(ocvFields))), ', '));

%% Find drive cycle data files
fprintf('\n=== Finding Drive Cycle Data Files ===\n');
% Search for parsedDriveCycle files (Reference: Simple_EventDetection_01.m line 48)
driveCycleFiles = dir(fullfile(driveCycleDir, 'parsedDriveCycle_*cyc_filtered.mat'));

if isempty(driveCycleFiles)
    error('No parsedDriveCycle files found in %s. Please check the path.', driveCycleDir);
end

fprintf('Found %d drive cycle files in: %s\n', length(driveCycleFiles), driveCycleDir);

%% Initialize results table (Long Format)
% Columns: Cycle, SOC_Level, DC_Type, E_net, delta_E_stored, Efficiency
resultsTable = table();

%% Initialize parallel pool if available
fprintf('\n=== Checking Parallel Computing Availability ===\n');
try
    poolObj = gcp('nocreate');
    if isempty(poolObj)
        poolObj = parpool('local');
        fprintf('Parallel pool started with %d workers\n', poolObj.NumWorkers);
    else
        fprintf('Parallel pool already exists with %d workers\n', poolObj.NumWorkers);
    end
    useParallel = true;
catch ME
    fprintf('Parallel Computing Toolbox not available or error: %s\n', ME.message);
    fprintf('Continuing with sequential processing...\n');
    useParallel = false;
end

%% Process each drive cycle file
fprintf('\n=== Processing Drive Cycle Data ===\n');

% Pre-allocate cell array for storing results from each file
fileResults = cell(length(driveCycleFiles), 1);
fileVisData = cell(length(driveCycleFiles), 1);

% Process files in parallel
if useParallel
    fprintf('Processing %d files in parallel...\n', length(driveCycleFiles));
    % Broadcast OCV_data to all workers
    OCV_data_broadcast = OCV_data;  %#ok<NASGU>
    parfor fileIdx = 1:length(driveCycleFiles)
        % Use broadcasted OCV_data
        OCV_data = OCV_data_broadcast;  %#ok<PFBNS>
    fileName = driveCycleFiles(fileIdx).name;
    
    % Extract cycle number from filename
    token = regexp(fileName, 'parsedDriveCycle_(\d+cyc)_filtered', 'tokens');
    if isempty(token)
        continue;
    end
    cycleType = token{1}{1};
    
    % Extract numeric cycle value
    cycleNum = str2double(strrep(cycleType, 'cyc', ''));
    
    % Find matching RPT cycle
    rptIdx = find(rptCycleNums == cycleNum, 1);
    if isempty(rptIdx)
        continue;
    end
    
    % Load drive cycle data
    filePath = fullfile(driveCycleFiles(fileIdx).folder, fileName);
    try
        driveData = load(filePath);
    catch
        continue;
    end
    
    % Find the variable name
    varName = sprintf('parsedDriveCycle_%s', cycleType);
    if ~isfield(driveData, varName)
        continue;
    end
    
    parsedData = driveData.(varName);
    
    % Find Ch09 data
    ch9FieldName = sprintf('ch9_Drive_%s', cycleType);
    if ~isfield(parsedData, ch9FieldName)
        continue;
    end
    
    ch9Data = parsedData.(ch9FieldName);
    
    % Get OCV data for this cycle (need to load OCV_data in each worker)
    % Note: OCV_data needs to be broadcast to workers
    ocvFieldName = sprintf('OCV_integrated_%d', cycleNum);
    if ~isfield(OCV_data, ocvFieldName)
        continue;
    end
    
    ocvCycleData = OCV_data.(ocvFieldName);
    
    % Extract OCV curve data
    SOC_grid = ocvCycleData.SOC_grid;
    V_avg_SOC = ocvCycleData.V_avg_SOC;
    
    % Use Ch09's actual aged capacity (Q_aged) instead of mean capacity
    capFieldName = sprintf('static_capacity_ch09_rpt%d', cycleNum);
    if isfield(OCV_data, capFieldName)
        Q_aged = OCV_data.(capFieldName);  % Ch09's actual degraded capacity
    else
        Q_aged = ocvCycleData.mean_capacity;  % Fallback to mean capacity if Ch09 data not available
    end
    Q_nominal = Q_aged;  % Use Q_aged for energy calculations
    
    % Initialize storage for this cycle
    cycleResults = cell(length(socLevels) * (length(dcTypes) + 1), 1);
    cycleVisData = struct();
    resultIdx = 1;
    
    % Process each SOC level
    for socIdx = 1:length(socLevels)
        socLevel = socLevels{socIdx};
        
        if ~isfield(ch9Data, socLevel)
            continue;
        end
        
        socData = ch9Data.(socLevel);
        cycleVisData.(socLevel) = struct();
        
        % Initialize accumulators for Net Efficiency
        Total_E_stored_ini = 0; Total_E_stored_aft = 0;
        Total_E_input = 0; Total_E_output = 0;
        Total_E_actual_dis = 0; Total_E_actual_chg = 0;
        Total_E_stored_theoretical_dis = 0; Total_E_stored_theoretical_chg = 0;
        
        % Process each drive cycle type
        for dcIdx = 1:length(dcTypes)
            dcType = dcTypes{dcIdx};
            
            if ~isfield(socData, dcType)
                continue;
            end
            
            dcData = socData.(dcType);
            
            % Check required fields
            if ~isfield(dcData, 'V') || ~isfield(dcData, 'I') || ~isfield(dcData, 't')
                continue;
            end
            
            V = dcData.V;
            I = dcData.I;
            t = dcData.t;
            
            % Convert time to seconds if needed
            if isa(t, 'duration')
                t_sec = seconds(t);
            else
                t_sec = t;
            end
            
            % Remove NaN and ensure column vectors
            validIdx = ~isnan(V) & ~isnan(I) & ~isnan(t_sec);
            if sum(validIdx) < 2
                continue;
            end
            
            V = V(validIdx);
            I = I(validIdx);
            t_sec = t_sec(validIdx);
            
            % Normalize time to start from 0
            t_sec = t_sec - t_sec(1);
            
            %% Step 1: Determine SOC1 and SOC2 from initial and final voltages
            V_ini = V(1);
            V_aft = V(end);
            
            % Interpolate SOC from OCV curve (inverse lookup)
            % Note: V_avg_SOC is increasing with SOC, so we can use interp1
            SOC1 = interp1(V_avg_SOC, SOC_grid, V_ini, 'linear');
            SOC2 = interp1(V_avg_SOC, SOC_grid, V_aft, 'linear');
            
            % Check for extrapolation warnings (silent in parallel)
            
            %% Step 2: Calculate real-time SOC using correction formula
            % SOC(t) = SOC1 + (SOC2 - SOC1) × [∫₀ᵗ I(τ)dτ] / [∫₀ᵉⁿᵈ I(τ)dτ]
            
            % Calculate cumulative current integral
            % Use trapz for numerical integration
            cumQ = cumtrapz(t_sec, I);  % Cumulative charge (Ah)
            Q_total = trapz(t_sec, I);  % Total charge (Ah)
            
            % Apply SOC correction formula
            if abs(Q_total) > 1e-6  % Avoid division by zero
                SOC_t = SOC1 + (SOC2 - SOC1) * (cumQ / Q_total);
            else
                % If no net current, SOC remains constant
                SOC_t = SOC1 * ones(size(t_sec));
            end
            
            % Ensure SOC stays within [0, 100]
            SOC_t = max(0, min(100, SOC_t));
            
            %% Step 3: Calculate real-time OCV(t) from SOC(t)
            OCV_t = interp1(SOC_grid, V_avg_SOC, SOC_t, 'linear');
            
            %% Step 4: Calculate ΔE_stored (시점 기반 에너지 계산)
            % 주행 전후 휴지기 끝점(OCV_1, OCV_2)에서의 에너지 값만 사용
            % E_stored(SOC) = (Q_nom / 100) × ∫₀^SOC OCV(SOC) dSOC
            % ΔE_stored = E_stored(ini) - E_stored(aft) (방전 시 양수)
            
            % Create fine SOC grid for integration
            soc_fine = linspace(0, 100, 10001);  % Fine grid for accurate integration
            ocv_fine = interp1(SOC_grid, V_avg_SOC, soc_fine, 'linear');
            
            % Calculate E_stored as cumulative integral
            % E_stored(SOC) = (Q_nominal/100) × ∫₀^SOC OCV(SOC') dSOC'
            E_stored_cum = cumtrapz(soc_fine, ocv_fine) * (Q_nominal / 100);  % Wh
            
            % Get E_stored at initial and final points (휴지기 끝점)
            E_stored_ini = interp1(soc_fine, E_stored_cum, SOC1, 'linear');  % E_stored(SOC1)
            E_stored_aft = interp1(soc_fine, E_stored_cum, SOC2, 'linear');  % E_stored(SOC2)
            
            % ΔE_stored = E_stored(ini) - E_stored(aft) (방전 시 양수)
            delta_E_stored = E_stored_ini - E_stored_aft;  % Wh
            
            %% [수정] Step 5: 개별 에너지 적분 (Actual & Theoretical)
            % 주행 중 발생하는 모든 충전 이벤트와 방전 이벤트를 각각 적분
            % Apply current threshold: only integrate when |I| >= I_threshold
            th = I_threshold;  % Threshold for clarity
            
            % Find discharge and charge regions with threshold
            discharge_mask = I < -th;  % 방전 구간 (전류가 -th 미만)
            charge_mask = I > th;      % 충전 구간 (전류가 th 초과)
            
            % 방전 에너지 (Actual & Theoretical)
            if sum(discharge_mask) >= 2
                % E_dis,actual = ∫(I<-th) |V·I| dt
                E_actual_dis = abs(trapz(t_sec(discharge_mask), V(discharge_mask) .* I(discharge_mask)) / 3600);
                % E_dis,ocv = ∫(I<-th) |OCV·I| dt
                E_stored_theoretical_dis = abs(trapz(t_sec(discharge_mask), OCV_t(discharge_mask) .* I(discharge_mask)) / 3600);
            else
                E_actual_dis = 0;
                E_stored_theoretical_dis = 0;
            end
            
            % 충전 에너지 (Actual & Theoretical)
            if sum(charge_mask) >= 2
                % E_chg,actual = ∫(I>th) (V·I) dt
                E_actual_chg = trapz(t_sec(charge_mask), V(charge_mask) .* I(charge_mask)) / 3600;
                % E_chg,ocv = ∫(I>th) (OCV·I) dt
                E_stored_theoretical_chg = trapz(t_sec(charge_mask), OCV_t(charge_mask) .* I(charge_mask)) / 3600;
            else
                E_actual_chg = 0;
                E_stored_theoretical_chg = 0;
            end
            
            % Legacy variables for compatibility
            E_output = E_actual_dis;
            E_input = E_actual_chg;
            E_net = E_output - E_input;
            
            %% [수정] Step 6: 논문 기준 효율 산출 (이원화 완료)
            % 주행 중 발생하는 모든 충전 이벤트와 방전 이벤트를 각각 적분하여 효율 계산
            % 
            % 충전 효율 (η_chg): "넣어준 전기 에너지 대비 실제로 화학적으로 저장된 양"
            %   E_chg,actual = ∫(I>th) (V·I) dt
            %   E_chg,ocv = ∫(I>th) (OCV·I) dt
            %   η_chg = E_chg,ocv / E_chg,actual
            % 
            % 방전 효율 (η_dis): "내부에서 꺼낸 화학 에너지 대비 실제로 밖으로 나온 양"
            %   E_dis,actual = ∫(I<-th) |V·I| dt
            %   E_dis,ocv = ∫(I<-th) |OCV·I| dt
            %   η_dis = E_dis,actual / E_dis,ocv
            
            % 1. 방전 효율: 얻은 것(출력) / 지불한 것(이론적 소모)
            if E_stored_theoretical_dis > 1e-4
                Discharge_Efficiency = E_actual_dis / E_stored_theoretical_dis;
            else
                Discharge_Efficiency = NaN;  % 이론적 방전 에너지가 너무 작으면 계산 불가
            end
            
            % 2. 충전 효율: 얻은 것(이론적 저장) / 지불한 것(실제 투입)
            if E_actual_chg > 1e-4
                Charge_Efficiency = E_stored_theoretical_chg / E_actual_chg;
            else
                Charge_Efficiency = NaN;  % 실제 충전 에너지가 너무 작으면 계산 불가
            end
            
            % 효율 범위 제한 (0~1 또는 0~100%)
            % 물리적으로 효율은 0~1 사이여야 함
            if ~isnan(Discharge_Efficiency)
                if Discharge_Efficiency < 0
                    Discharge_Efficiency = NaN;  % 음수는 물리적으로 불가능
                elseif Discharge_Efficiency > 1
                    Discharge_Efficiency = NaN;  % 100% 초과는 계산 오류
                end
            end
            
            if ~isnan(Charge_Efficiency)
                if Charge_Efficiency < 0
                    Charge_Efficiency = NaN;  % 음수는 물리적으로 불가능
                elseif Charge_Efficiency > 1
                    Charge_Efficiency = NaN;  % 100% 초과는 계산 오류
                end
            end
            
            % Legacy variables for backward compatibility
            % E_stored_theoretical_dis, E_stored_theoretical_chg는 Step 5에서 계산됨
            
            % Legacy variables for backward compatibility
            % Determine dominant efficiency type based on delta_E_stored
            if delta_E_stored > 0
                efficiency_type = 'Discharge';
                efficiency = Discharge_Efficiency;  % Use discharge efficiency
            elseif delta_E_stored < 0
                efficiency_type = 'Charge';
                efficiency = Charge_Efficiency;  % Use charge efficiency
            else
                efficiency_type = 'Unknown';
                efficiency = NaN;
            end
            
            %% Step 7: Calculate Coulomb Counting SOC (for comparison)
            % Simple current integration without correction
            cumQ_cc = cumtrapz(t_sec, I);  % Cumulative charge (Ah)
            SOC_cc = SOC1 + (cumQ_cc / Q_nominal) * 100;  % Coulomb counting SOC
            SOC_cc = max(0, min(100, SOC_cc));  % Clamp to [0, 100]
            
            %% Step 8: Store visualization data
            cycleVisData.(socLevel).(dcType).t_sec = t_sec;
            cycleVisData.(socLevel).(dcType).V = V;
            cycleVisData.(socLevel).(dcType).I = I;
            cycleVisData.(socLevel).(dcType).SOC_t = SOC_t;  % Corrected SOC
            cycleVisData.(socLevel).(dcType).SOC_cc = SOC_cc;  % Coulomb counting SOC
            cycleVisData.(socLevel).(dcType).OCV_t = OCV_t;
            cycleVisData.(socLevel).(dcType).E_net = E_net;
            cycleVisData.(socLevel).(dcType).delta_E_stored = delta_E_stored;
            cycleVisData.(socLevel).(dcType).efficiency = efficiency;
            cycleVisData.(socLevel).(dcType).SOC1 = SOC1;
            cycleVisData.(socLevel).(dcType).SOC2 = SOC2;
            cycleVisData.(socLevel).(dcType).cycleNum = cycleNum;
            cycleVisData.(socLevel).(dcType).Q_nominal = Q_nominal;
            
            % Accumulate for Net Efficiency
            Total_E_stored_ini = Total_E_stored_ini + E_stored_ini;
            Total_E_stored_aft = Total_E_stored_aft + E_stored_aft;
            Total_E_input = Total_E_input + E_input;
            Total_E_output = Total_E_output + E_output;
            Total_E_actual_dis = Total_E_actual_dis + E_actual_dis;
            Total_E_actual_chg = Total_E_actual_chg + E_actual_chg;
            Total_E_stored_theoretical_dis = Total_E_stored_theoretical_dis + E_stored_theoretical_dis;
            Total_E_stored_theoretical_chg = Total_E_stored_theoretical_chg + E_stored_theoretical_chg;

            %% Step 9: Store results for this DC (with both Discharge and Charge efficiencies)
            cycleResults{resultIdx} = struct('Cycle', cycleType, 'SOC_Level', socLevel, ...
                'DC_Type', dcType, 'E_stored_ini', E_stored_ini, 'E_stored_aft', E_stored_aft, ...
                'delta_E_stored', delta_E_stored, 'E_input', E_input, 'E_output', E_output, ...
                'E_net', E_net, 'E_actual_dis', E_actual_dis, 'E_actual_chg', E_actual_chg, ...
                'E_stored_theoretical_dis', E_stored_theoretical_dis, 'E_stored_theoretical_chg', E_stored_theoretical_chg, ...
                'Discharge_Efficiency', Discharge_Efficiency, 'Charge_Efficiency', Charge_Efficiency, ...
                'Efficiency', efficiency, 'Efficiency_Type', efficiency_type);
            resultIdx = resultIdx + 1;
            
        end  % DC loop
        
        % Calculate Net Efficiency (Aggregation of 8 DCs)
        if Total_E_stored_theoretical_dis > 1e-4
            Net_Discharge_Efficiency = Total_E_actual_dis / Total_E_stored_theoretical_dis;
        else
            Net_Discharge_Efficiency = NaN;
        end
        
        if Total_E_actual_chg > 1e-4
            Net_Charge_Efficiency = Total_E_stored_theoretical_chg / Total_E_actual_chg;
        else
            Net_Charge_Efficiency = NaN;
        end
        
        % Boundary checks for Net Efficiency
        if ~isnan(Net_Discharge_Efficiency)
            if Net_Discharge_Efficiency < 0, Net_Discharge_Efficiency = NaN;
            elseif Net_Discharge_Efficiency > 1, Net_Discharge_Efficiency = NaN; end
        end
        if ~isnan(Net_Charge_Efficiency)
            if Net_Charge_Efficiency < 0, Net_Charge_Efficiency = NaN;
            elseif Net_Charge_Efficiency > 1, Net_Charge_Efficiency = NaN; end
        end
        
        % Net Energy Balance
        Total_E_net = Total_E_output - Total_E_input;
        Total_delta_E_stored = Total_E_stored_ini - Total_E_stored_aft;
        
        % Determine Net Efficiency Type
        if Total_delta_E_stored > 0
            net_eff_type = 'Discharge';
            net_efficiency = Net_Discharge_Efficiency;
        elseif Total_delta_E_stored < 0
            net_eff_type = 'Charge';
            net_efficiency = Net_Charge_Efficiency;
        else
            net_eff_type = 'Unknown';
            net_efficiency = NaN;
        end
        
        % Store Net Results
        cycleResults{resultIdx} = struct('Cycle', cycleType, 'SOC_Level', socLevel, ...
            'DC_Type', 'Net_Efficiency', 'E_stored_ini', Total_E_stored_ini, ...
            'E_stored_aft', Total_E_stored_aft, ...
            'delta_E_stored', Total_delta_E_stored, 'E_input', Total_E_input, 'E_output', Total_E_output, ...
            'E_net', Total_E_net, 'E_actual_dis', Total_E_actual_dis, 'E_actual_chg', Total_E_actual_chg, ...
            'E_stored_theoretical_dis', Total_E_stored_theoretical_dis, ...
            'E_stored_theoretical_chg', Total_E_stored_theoretical_chg, ...
            'Discharge_Efficiency', Net_Discharge_Efficiency, 'Charge_Efficiency', Net_Charge_Efficiency, ...
            'Efficiency', net_efficiency, 'Efficiency_Type', net_eff_type);
        resultIdx = resultIdx + 1;
    end  % SOC loop
    
    % Store results for this file
    fileResults{fileIdx} = cycleResults;
    fileVisData{fileIdx} = struct('cycleType', cycleType, 'cycleNum', cycleNum, ...
        'visData', cycleVisData, 'Q_nominal', Q_nominal, 'SOC_grid', SOC_grid, ...
        'V_avg_SOC', V_avg_SOC);
    
end  % File loop (parallel)

else
    % Sequential processing (fallback)
    fprintf('Processing %d files sequentially...\n', length(driveCycleFiles));
    fileResults = cell(length(driveCycleFiles), 1);
    fileVisData = cell(length(driveCycleFiles), 1);
    
    for fileIdx = 1:length(driveCycleFiles)
        fileName = driveCycleFiles(fileIdx).name;
        
        % Extract cycle number from filename
        token = regexp(fileName, 'parsedDriveCycle_(\d+cyc)_filtered', 'tokens');
        if isempty(token)
            continue;
        end
        cycleType = token{1}{1};
        cycleNum = str2double(strrep(cycleType, 'cyc', ''));
        
        % Find matching RPT cycle
        rptIdx = find(rptCycleNums == cycleNum, 1);
        if isempty(rptIdx)
            continue;
        end
        
        fprintf('  Processing %s...\n', cycleType);
        
        % Load drive cycle data
        filePath = fullfile(driveCycleFiles(fileIdx).folder, fileName);
        try
            driveData = load(filePath);
        catch
            continue;
        end
        
        % Find the variable name
        varName = sprintf('parsedDriveCycle_%s', cycleType);
        if ~isfield(driveData, varName)
            continue;
        end
        
        parsedData = driveData.(varName);
        
        % Find Ch09 data
        ch9FieldName = sprintf('ch9_Drive_%s', cycleType);
        if ~isfield(parsedData, ch9FieldName)
            continue;
        end
        
        ch9Data = parsedData.(ch9FieldName);
        
        % Get OCV data for this cycle
        ocvFieldName = sprintf('OCV_integrated_%d', cycleNum);
        if ~isfield(OCV_data, ocvFieldName)
            continue;
        end
        
        ocvCycleData = OCV_data.(ocvFieldName);
        SOC_grid = ocvCycleData.SOC_grid;
        V_avg_SOC = ocvCycleData.V_avg_SOC;
        
        % Use Ch09's actual aged capacity (Q_aged) instead of mean capacity
        capFieldName = sprintf('static_capacity_ch09_rpt%d', cycleNum);
        if isfield(OCV_data, capFieldName)
            Q_aged = OCV_data.(capFieldName);  % Ch09's actual degraded capacity
        else
            Q_aged = ocvCycleData.mean_capacity;  % Fallback to mean capacity if Ch09 data not available
        end
        Q_nominal = Q_aged;  % Use Q_aged for energy calculations
        
        % Initialize storage
        cycleResults = cell(length(socLevels) * (length(dcTypes) + 1), 1);
        cycleVisData = struct();
        resultIdx = 1;
        
        % Process each SOC level and DC (same logic as parallel version)
        for socIdx = 1:length(socLevels)
            socLevel = socLevels{socIdx};
            if ~isfield(ch9Data, socLevel)
                continue;
            end
            
            socData = ch9Data.(socLevel);
            cycleVisData.(socLevel) = struct();
            
            % Initialize accumulators for Net Efficiency
            Total_E_stored_ini = 0; Total_E_stored_aft = 0;
            Total_E_input = 0; Total_E_output = 0;
            Total_E_actual_dis = 0; Total_E_actual_chg = 0;
            Total_E_stored_theoretical_dis = 0; Total_E_stored_theoretical_chg = 0;
            
            for dcIdx = 1:length(dcTypes)
                dcType = dcTypes{dcIdx};
                if ~isfield(socData, dcType)
                    continue;
                end
                
                dcData = socData.(dcType);
                if ~isfield(dcData, 'V') || ~isfield(dcData, 'I') || ~isfield(dcData, 't')
                    continue;
                end
                
                V = dcData.V;
                I = dcData.I;
                t = dcData.t;
                
                if isa(t, 'duration')
                    t_sec = seconds(t);
                else
                    t_sec = t;
                end
                
                validIdx = ~isnan(V) & ~isnan(I) & ~isnan(t_sec);
                if sum(validIdx) < 2
                    continue;
                end
                
                V = V(validIdx);
                I = I(validIdx);
                t_sec = t_sec(validIdx);
                t_sec = t_sec - t_sec(1);
                
                % Same calculations as parallel version
                V_ini = V(1);
                V_aft = V(end);
                SOC1 = interp1(V_avg_SOC, SOC_grid, V_ini, 'linear');
                SOC2 = interp1(V_avg_SOC, SOC_grid, V_aft, 'linear');
                
                cumQ = cumtrapz(t_sec, I);
                Q_total = trapz(t_sec, I);
                
                if abs(Q_total) > 1e-6
                    SOC_t = SOC1 + (SOC2 - SOC1) * (cumQ / Q_total);
                else
                    SOC_t = SOC1 * ones(size(t_sec));
                end
                SOC_t = max(0, min(100, SOC_t));
                
                OCV_t = interp1(SOC_grid, V_avg_SOC, SOC_t, 'linear');
                
                soc_fine = linspace(0, 100, 10001);
                ocv_fine = interp1(SOC_grid, V_avg_SOC, soc_fine, 'linear');
                E_stored_cum = cumtrapz(soc_fine, ocv_fine) * (Q_nominal / 100);
                E_stored_ini = interp1(soc_fine, E_stored_cum, SOC1, 'linear');
                E_stored_aft = interp1(soc_fine, E_stored_cum, SOC2, 'linear');
                delta_E_stored = E_stored_ini - E_stored_aft;
                
                %% [수정] Step 5: 개별 에너지 적분 (Actual & Theoretical)
                % 주행 중 발생하는 모든 충전 이벤트와 방전 이벤트를 각각 적분
                % Apply current threshold: only integrate when |I| >= I_threshold
                th = I_threshold;  % Threshold for clarity
                
                % Find discharge and charge regions with threshold
                discharge_mask = I < -th;  % 방전 구간 (전류가 -th 미만)
                charge_mask = I > th;      % 충전 구간 (전류가 th 초과)
                
                % 방전 에너지 (Actual & Theoretical)
                if sum(discharge_mask) >= 2
                    % E_dis,actual = ∫(I<-th) |V·I| dt
                    E_actual_dis = abs(trapz(t_sec(discharge_mask), V(discharge_mask) .* I(discharge_mask)) / 3600);
                    % E_dis,ocv = ∫(I<-th) |OCV·I| dt
                    E_stored_theoretical_dis = abs(trapz(t_sec(discharge_mask), OCV_t(discharge_mask) .* I(discharge_mask)) / 3600);
                else
                    E_actual_dis = 0;
                    E_stored_theoretical_dis = 0;
                end
                
                % 충전 에너지 (Actual & Theoretical)
                if sum(charge_mask) >= 2
                    % E_chg,actual = ∫(I>th) (V·I) dt
                    E_actual_chg = trapz(t_sec(charge_mask), V(charge_mask) .* I(charge_mask)) / 3600;
                    % E_chg,ocv = ∫(I>th) (OCV·I) dt
                    E_stored_theoretical_chg = trapz(t_sec(charge_mask), OCV_t(charge_mask) .* I(charge_mask)) / 3600;
                else
                    E_actual_chg = 0;
                    E_stored_theoretical_chg = 0;
                end
                
                % Legacy variables for compatibility
                E_output = E_actual_dis;
                E_input = E_actual_chg;
                E_net = E_output - E_input;
                
                %% [수정] Step 6: 논문 기준 효율 산출 (이원화 완료)
                % 주행 중 발생하는 모든 충전 이벤트와 방전 이벤트를 각각 적분하여 효율 계산
                % 
                % 충전 효율 (η_chg): "넣어준 전기 에너지 대비 실제로 화학적으로 저장된 양"
                %   E_chg,actual = ∫(I>th) (V·I) dt
                %   E_chg,ocv = ∫(I>th) (OCV·I) dt
                %   η_chg = E_chg,ocv / E_chg,actual
                % 
                % 방전 효율 (η_dis): "내부에서 꺼낸 화학 에너지 대비 실제로 밖으로 나온 양"
                %   E_dis,actual = ∫(I<-th) |V·I| dt
                %   E_dis,ocv = ∫(I<-th) |OCV·I| dt
                %   η_dis = E_dis,actual / E_dis,ocv
                
                % 1. 방전 효율: 얻은 것(출력) / 지불한 것(이론적 소모)
                if E_stored_theoretical_dis > 1e-4
                    Discharge_Efficiency = E_actual_dis / E_stored_theoretical_dis;
                else
                    Discharge_Efficiency = NaN;  % 이론적 방전 에너지가 너무 작으면 계산 불가
                end
                
                % 2. 충전 효율: 얻은 것(이론적 저장) / 지불한 것(실제 투입)
                if E_actual_chg > 1e-4
                    Charge_Efficiency = E_stored_theoretical_chg / E_actual_chg;
                else
                    Charge_Efficiency = NaN;  % 실제 충전 에너지가 너무 작으면 계산 불가
                end
                
                % 효율 범위 제한 (0~1 또는 0~100%)
                % 물리적으로 효율은 0~1 사이여야 함
                if ~isnan(Discharge_Efficiency)
                    if Discharge_Efficiency < 0
                        Discharge_Efficiency = NaN;  % 음수는 물리적으로 불가능
                    elseif Discharge_Efficiency > 1
                        Discharge_Efficiency = NaN;  % 100% 초과는 계산 오류
                    end
                end
                
                if ~isnan(Charge_Efficiency)
                    if Charge_Efficiency < 0
                        Charge_Efficiency = NaN;  % 음수는 물리적으로 불가능
                    elseif Charge_Efficiency > 1
                        Charge_Efficiency = NaN;  % 100% 초과는 계산 오류
                    end
                end
                
                % Legacy variables for backward compatibility
                % E_stored_theoretical_dis, E_stored_theoretical_chg는 Step 5에서 계산됨
                
                % Legacy variables for backward compatibility
                % Determine dominant efficiency type based on delta_E_stored
                if delta_E_stored > 0
                    efficiency_type = 'Discharge';
                    efficiency = Discharge_Efficiency;  % Use discharge efficiency
                elseif delta_E_stored < 0
                    efficiency_type = 'Charge';
                    efficiency = Charge_Efficiency;  % Use charge efficiency
                else
                    efficiency_type = 'Unknown';
                    efficiency = NaN;
                end
                
                cumQ_cc = cumtrapz(t_sec, I);
                SOC_cc = SOC1 + (cumQ_cc / Q_nominal) * 100;
                SOC_cc = max(0, min(100, SOC_cc));
                
                cycleVisData.(socLevel).(dcType).t_sec = t_sec;
                cycleVisData.(socLevel).(dcType).V = V;
                cycleVisData.(socLevel).(dcType).I = I;
                cycleVisData.(socLevel).(dcType).SOC_t = SOC_t;
                cycleVisData.(socLevel).(dcType).SOC_cc = SOC_cc;
                cycleVisData.(socLevel).(dcType).OCV_t = OCV_t;
                cycleVisData.(socLevel).(dcType).E_net = E_net;
                cycleVisData.(socLevel).(dcType).delta_E_stored = delta_E_stored;
                cycleVisData.(socLevel).(dcType).efficiency = efficiency;
                cycleVisData.(socLevel).(dcType).SOC1 = SOC1;
                cycleVisData.(socLevel).(dcType).SOC2 = SOC2;
                cycleVisData.(socLevel).(dcType).cycleNum = cycleNum;
                cycleVisData.(socLevel).(dcType).Q_nominal = Q_nominal;
                
                % Accumulate for Net Efficiency
                Total_E_stored_ini = Total_E_stored_ini + E_stored_ini;
                Total_E_stored_aft = Total_E_stored_aft + E_stored_aft;
                Total_E_input = Total_E_input + E_input;
                Total_E_output = Total_E_output + E_output;
                Total_E_actual_dis = Total_E_actual_dis + E_actual_dis;
                Total_E_actual_chg = Total_E_actual_chg + E_actual_chg;
                Total_E_stored_theoretical_dis = Total_E_stored_theoretical_dis + E_stored_theoretical_dis;
                Total_E_stored_theoretical_chg = Total_E_stored_theoretical_chg + E_stored_theoretical_chg;

                cycleResults{resultIdx} = struct('Cycle', cycleType, 'SOC_Level', socLevel, ...
                    'DC_Type', dcType, 'E_stored_ini', E_stored_ini, 'E_stored_aft', E_stored_aft, ...
                    'delta_E_stored', delta_E_stored, 'E_input', E_input, 'E_output', E_output, ...
                    'E_net', E_net, 'E_actual_dis', E_actual_dis, 'E_actual_chg', E_actual_chg, ...
                    'E_stored_theoretical_dis', E_stored_theoretical_dis, 'E_stored_theoretical_chg', E_stored_theoretical_chg, ...
                    'Discharge_Efficiency', Discharge_Efficiency, 'Charge_Efficiency', Charge_Efficiency, ...
                    'Efficiency', efficiency, 'Efficiency_Type', efficiency_type);
                resultIdx = resultIdx + 1;
            end
        
            % Calculate Net Efficiency (Aggregation of 8 DCs)
            if Total_E_stored_theoretical_dis > 1e-4
                Net_Discharge_Efficiency = Total_E_actual_dis / Total_E_stored_theoretical_dis;
            else
                Net_Discharge_Efficiency = NaN;
            end
            
            if Total_E_actual_chg > 1e-4
                Net_Charge_Efficiency = Total_E_stored_theoretical_chg / Total_E_actual_chg;
            else
                Net_Charge_Efficiency = NaN;
            end
            
            % Boundary checks for Net Efficiency
            if ~isnan(Net_Discharge_Efficiency)
                if Net_Discharge_Efficiency < 0, Net_Discharge_Efficiency = NaN;
                elseif Net_Discharge_Efficiency > 1, Net_Discharge_Efficiency = NaN; end
            end
            if ~isnan(Net_Charge_Efficiency)
                if Net_Charge_Efficiency < 0, Net_Charge_Efficiency = NaN;
                elseif Net_Charge_Efficiency > 1, Net_Charge_Efficiency = NaN; end
            end
            
            % Net Energy Balance
            Total_E_net = Total_E_output - Total_E_input;
            Total_delta_E_stored = Total_E_stored_ini - Total_E_stored_aft;
            
            % Determine Net Efficiency Type
            if Total_delta_E_stored > 0
                net_eff_type = 'Discharge';
                net_efficiency = Net_Discharge_Efficiency;
            elseif Total_delta_E_stored < 0
                net_eff_type = 'Charge';
                net_efficiency = Net_Charge_Efficiency;
            else
                net_eff_type = 'Unknown';
                net_efficiency = NaN;
            end
            
            % Store Net Results
            cycleResults{resultIdx} = struct('Cycle', cycleType, 'SOC_Level', socLevel, ...
                'DC_Type', 'Net_Efficiency', 'E_stored_ini', Total_E_stored_ini, ...
                'E_stored_aft', Total_E_stored_aft, ...
                'delta_E_stored', Total_delta_E_stored, 'E_input', Total_E_input, 'E_output', Total_E_output, ...
                'E_net', Total_E_net, 'E_actual_dis', Total_E_actual_dis, 'E_actual_chg', Total_E_actual_chg, ...
                'E_stored_theoretical_dis', Total_E_stored_theoretical_dis, ...
                'E_stored_theoretical_chg', Total_E_stored_theoretical_chg, ...
                'Discharge_Efficiency', Net_Discharge_Efficiency, 'Charge_Efficiency', Net_Charge_Efficiency, ...
                'Efficiency', net_efficiency, 'Efficiency_Type', net_eff_type);
            resultIdx = resultIdx + 1;
        end
        
        fileResults{fileIdx} = cycleResults;
        fileVisData{fileIdx} = struct('cycleType', cycleType, 'cycleNum', cycleNum, ...
            'visData', cycleVisData, 'Q_nominal', Q_nominal, 'SOC_grid', SOC_grid, ...
            'V_avg_SOC', V_avg_SOC);
    end
end

% Combine results from all files
fprintf('\n=== Combining Results ===\n');
for fileIdx = 1:length(fileResults)
    cycleResults = fileResults{fileIdx};
    if ~isempty(cycleResults)
        for dcIdx = 1:length(cycleResults)
            if ~isempty(cycleResults{dcIdx})
                % Check if Efficiency_Type field exists, if not set to 'Unknown'
                if isfield(cycleResults{dcIdx}, 'Efficiency_Type')
                    effType = cycleResults{dcIdx}.Efficiency_Type;
                else
                    effType = 'Unknown';
                    fprintf('  Warning: Efficiency_Type field missing for %s %s %s, using Unknown\n', ...
                        cycleResults{dcIdx}.Cycle, cycleResults{dcIdx}.SOC_Level, cycleResults{dcIdx}.DC_Type);
                end
                
                % Check if new efficiency fields exist
                if isfield(cycleResults{dcIdx}, 'Discharge_Efficiency')
                    dischargeEff = cycleResults{dcIdx}.Discharge_Efficiency;
                else
                    dischargeEff = NaN;
                end
                
                if isfield(cycleResults{dcIdx}, 'Charge_Efficiency')
                    chargeEff = cycleResults{dcIdx}.Charge_Efficiency;
                else
                    chargeEff = NaN;
                end
                
                newRow = table({cycleResults{dcIdx}.Cycle}, {cycleResults{dcIdx}.SOC_Level}, ...
                    {cycleResults{dcIdx}.DC_Type}, cycleResults{dcIdx}.E_stored_ini, ...
                    cycleResults{dcIdx}.E_stored_aft, cycleResults{dcIdx}.delta_E_stored, ...
                    cycleResults{dcIdx}.E_input, cycleResults{dcIdx}.E_output, ...
                    cycleResults{dcIdx}.E_net, dischargeEff, chargeEff, ...
                    cycleResults{dcIdx}.Efficiency, {effType}, ...
                    'VariableNames', {'Cycle', 'SOC_Level', 'DC_Type', 'E_stored_ini', ...
                    'E_stored_aft', 'delta_E_stored', 'E_input', 'E_output', 'E_net', ...
                    'Discharge_Efficiency', 'Charge_Efficiency', 'Efficiency', 'Efficiency_Type'});
                resultsTable = [resultsTable; newRow];
            end
        end
    end
end

%% Create visualizations (sequential, after data collection)
fprintf('\n=== Creating Visualizations ===\n');
figDir = fullfile(outputDir, 'Figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

for fileIdx = 1:length(fileVisData)
    fileData = fileVisData{fileIdx};
    if isempty(fileData) || ~isfield(fileData, 'visData')
        continue;
    end
    
    cycleType = fileData.cycleType;
    cycleNum = fileData.cycleNum;
    visData = fileData.visData;
    Q_nominal = fileData.Q_nominal;
    SOC_grid = fileData.SOC_grid;
    V_avg_SOC = fileData.V_avg_SOC;
    
    fprintf('Creating visualizations for %s...\n', cycleType);
    
    % Process each SOC level
    for socIdx = 1:length(socLevels)
        socLevel = socLevels{socIdx};
        
        if ~isfield(visData, socLevel)
            continue;
        end
        
        %% Create visualizations for this SOC level
        if isfield(visData, socLevel)
            
            % Create Figures directory
            figDir = fullfile(outputDir, 'Figures');
            if ~exist(figDir, 'dir')
                mkdir(figDir);
            end
            
            %% Visualization 1: SOC Profile Verification (8 subplots: DC1~DC8)
            % fprintf('\n    === SOC Verification Debug Info - %s %s ===\n', cycleType, socLevel);
            
            fig1 = figure('Name', sprintf('SOC Verification - %s %s', cycleType, socLevel), ...
                'Position', [100, 100, 1800, 1000]);
            set(gcf, 'Visible', 'off');
            
            % Initialize statistics storage (commented out - debugging not needed)
            % verificationStats = struct();
            
            for dcIdx = 1:length(dcTypes)
                dcType = dcTypes{dcIdx};
                
                if ~isfield(visData.(socLevel), dcType)
                    subplot(2, 4, dcIdx);
                    text(0.5, 0.5, sprintf('%s\nNo Data', dcType), ...
                        'HorizontalAlignment', 'center', 'FontSize', 12);
                    axis off;
                    continue;
                end
                
                subplot(2, 4, dcIdx);
                
                t_vis = visData.(socLevel).(dcType).t_sec / 60;  % Convert to minutes
                V_vis = visData.(socLevel).(dcType).V;  % Terminal voltage
                OCV_vis = visData.(socLevel).(dcType).OCV_t;  % OCV
                
                % Calculate voltage difference (error)
                V_diff = V_vis - OCV_vis;  % Terminal voltage - OCV (should be positive during discharge, negative during charge)
                V_diff_abs = abs(V_diff);  % Absolute difference
                
                % Calculate statistics (commented out - debugging not needed)
                % V_diff_mean = mean(V_diff);
                % V_diff_std = std(V_diff);
                % V_diff_max = max(V_diff_abs);
                % V_diff_min = min(V_diff_abs);
                % V_diff_rmse = sqrt(mean(V_diff.^2));  % Root Mean Square Error
                
                % Store statistics (commented out - debugging not needed)
                % verificationStats.(dcType).mean_diff = V_diff_mean;
                % verificationStats.(dcType).std_diff = V_diff_std;
                % verificationStats.(dcType).max_abs_diff = V_diff_max;
                % verificationStats.(dcType).min_abs_diff = V_diff_min;
                % verificationStats.(dcType).rmse = V_diff_rmse;
                % verificationStats.(dcType).mean_abs_diff = mean(V_diff_abs);
                
                % Print debugging information (commented out - debugging not needed)
                % fprintf('      %s:\n', dcType);
                % fprintf('        Mean difference (V_term - OCV): %.4f V\n', V_diff_mean);
                % fprintf('        Mean absolute difference: %.4f V\n', mean(V_diff_abs));
                % fprintf('        Max absolute difference: %.4f V\n', V_diff_max);
                % fprintf('        Min absolute difference: %.4f V\n', V_diff_min);
                % fprintf('        Std deviation: %.4f V\n', V_diff_std);
                % fprintf('        RMSE: %.4f V\n', V_diff_rmse);
                % fprintf('        Relative error (RMSE/mean_V): %.2f%%\n', ...
                %     V_diff_rmse / mean(V_vis) * 100);
                
                % Calculate common y-axis range for both voltages
                v_min = min([min(V_vis), min(OCV_vis)]);
                v_max = max([max(V_vis), max(OCV_vis)]);
                v_range = v_max - v_min;
                ylim_common = [v_min - v_range*0.01, v_max + v_range*0.01];
                
                % Left y-axis: Terminal Voltage
                yyaxis left
                plot(t_vis, V_vis, 'b-', 'LineWidth', 2, 'DisplayName', 'Terminal Voltage');
                ylabel('Terminal Voltage [V]', 'FontSize', 10);
                ylim(ylim_common);
                
                % Right y-axis: OCV
                yyaxis right
                plot(t_vis, OCV_vis, 'r--', 'LineWidth', 1.5, 'DisplayName', 'OCV');
                ylabel('OCV [V]', 'FontSize', 10);
                ylim(ylim_common);
                
                % Add title
                xlabel('Time [min]', 'FontSize', 10);
                title(sprintf('%s', dcType), 'FontSize', 11, 'FontWeight', 'bold');
                legend('Location', 'best', 'FontSize', 9);
                grid on;
                grid minor;
            end
            
            % Print summary statistics (commented out - debugging not needed)
            % dcNames = fieldnames(verificationStats);
            % if ~isempty(dcNames)
            %     all_rmse = [];
            %     all_mean_abs = [];
            %     for i = 1:length(dcNames)
            %         all_rmse = [all_rmse, verificationStats.(dcNames{i}).rmse];
            %         all_mean_abs = [all_mean_abs, verificationStats.(dcNames{i}).mean_abs_diff];
            %     end
            %     fprintf('      Summary:\n');
            %     fprintf('        Average RMSE across all DCs: %.4f V\n', mean(all_rmse));
            %     fprintf('        Average mean absolute diff: %.4f V\n', mean(all_mean_abs));
            %     fprintf('        Max RMSE: %.4f V (%s)\n', max(all_rmse), ...
            %         dcNames{all_rmse == max(all_rmse)});
            %     fprintf('        Min RMSE: %.4f V (%s)\n', min(all_rmse), ...
            %         dcNames{all_rmse == min(all_rmse)});
            % end
            
            sgtitle(sprintf('SOC Profile Verification - %s %s\n(Terminal Voltage vs OCV)', ...
                cycleType, socLevel), 'FontSize', 14, 'FontWeight', 'bold');
            
            fig1Name = sprintf('01_SOC_Verification_%s_%s', cycleType, socLevel);
            savefig(fig1, fullfile(figDir, [fig1Name, '.fig']));
            close(fig1);
            fprintf('    Saved figure: %s\n', fig1Name);
            
            %% Visualization 2: V-SOC Hysteresis Loop (8 subplots: DC1~DC8)
            fig2 = figure('Name', sprintf('V-SOC Hysteresis - %s %s', cycleType, socLevel), ...
                'Position', [100, 100, 1800, 1000]);
            set(gcf, 'Visible', 'off');
            
            for dcIdx = 1:length(dcTypes)
                dcType = dcTypes{dcIdx};
                
                if ~isfield(visData.(socLevel), dcType)
                    subplot(2, 4, dcIdx);
                    text(0.5, 0.5, sprintf('%s\nNo Data', dcType), ...
                        'HorizontalAlignment', 'center', 'FontSize', 12);
                    axis off;
                    continue;
                end
                
                subplot(2, 4, dcIdx);
                
                SOC_vis = visData.(socLevel).(dcType).SOC_t;
                V_vis = visData.(socLevel).(dcType).V;
                OCV_vis = visData.(socLevel).(dcType).OCV_t;
                eff_vis = visData.(socLevel).(dcType).efficiency;
                
                % Plot OCV curve (equilibrium state)
                plot(SOC_grid, V_avg_SOC, 'k-', 'LineWidth', 2.5, 'DisplayName', 'OCV (Equilibrium)');
                hold on;
                
                % Plot actual terminal voltage trajectory
                plot(SOC_vis, V_vis, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Terminal Voltage');
                
                % Fill area between OCV and terminal voltage (energy loss)
                % Interpolate OCV to match SOC_vis points
                ocv_at_soc = interp1(SOC_grid, V_avg_SOC, SOC_vis, 'linear');
                
                % Remove duplicate or near-duplicate points to avoid fill warnings
                % Use unique with tolerance for SOC values (0.01% tolerance)
                [~, uniqueIdx] = unique(round(SOC_vis * 1000) / 1000, 'stable');
                if length(uniqueIdx) < length(SOC_vis)
                    SOC_vis_clean = SOC_vis(uniqueIdx);
                    V_vis_clean = V_vis(uniqueIdx);
                    ocv_at_soc_clean = ocv_at_soc(uniqueIdx);
                else
                    SOC_vis_clean = SOC_vis;
                    V_vis_clean = V_vis;
                    ocv_at_soc_clean = ocv_at_soc;
                end
                
                % Create fill polygon: [SOC_vis, reverse(SOC_vis)], [OCV, reverse(V)]
                soc_fill = [SOC_vis_clean(:); flipud(SOC_vis_clean(:))];
                v_fill = [ocv_at_soc_clean(:); flipud(V_vis_clean(:))];
                
                % Remove any NaN or Inf values
                validIdx = isfinite(soc_fill) & isfinite(v_fill);
                if sum(validIdx) > 2
                    soc_fill = soc_fill(validIdx);
                    v_fill = v_fill(validIdx);
                    
                    % Suppress fill warnings and use fill with explicit settings
                    warning('off', 'MATLAB:fill:EmptyPolygon');
                    try
                        fill(soc_fill, v_fill, 'r', 'FaceAlpha', 0.15, ...
                            'EdgeColor', 'r', 'EdgeAlpha', 0.3, ...
                            'DisplayName', 'Energy Loss Area');
                    catch
                        % If fill fails, use patch as fallback
                        patch(soc_fill, v_fill, 'r', 'FaceAlpha', 0.15, ...
                            'EdgeColor', 'none', 'DisplayName', 'Energy Loss Area');
                    end
                    warning('on', 'MATLAB:fill:EmptyPolygon');
                end
                
                xlabel('SOC [%]', 'FontSize', 10);
                ylabel('Voltage [V]', 'FontSize', 10);
                if ~isnan(eff_vis)
                    title(sprintf('%s (Eff=%.2f%%)', dcType, eff_vis*100), 'FontSize', 11, 'FontWeight', 'bold');
                else
                    title(sprintf('%s (Eff=NaN)', dcType), 'FontSize', 11, 'FontWeight', 'bold');
                end
                legend('Location', 'best', 'FontSize', 9);
                grid on;
                grid minor;
                hold off;
            end
            
            sgtitle(sprintf('V-SOC Hysteresis Loop - %s %s\n(Energy Loss = Area between OCV and Terminal Voltage)', ...
                cycleType, socLevel), 'FontSize', 14, 'FontWeight', 'bold');
            
            fig2Name = sprintf('02_VSOC_Hysteresis_%s_%s', cycleType, socLevel);
            savefig(fig2, fullfile(figDir, [fig2Name, '.fig']));
            close(fig2);
            fprintf('    Saved figure: %s\n', fig2Name);
        end
        
    end  % SOC loop
end  % File loop

%% Visualization 3: Efficiency vs SOH (Bar Graph for SOC90) - Separated by Discharge/Charge
fprintf('\n=== Creating Efficiency vs SOH Visualization (SOC90) ===\n');
figDir = fullfile(outputDir, 'Figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

% Get Ch09 C/3 discharge capacity for SOH calculation
% OCV_data contains: static_capacity_ch09_rpt0, static_capacity_ch09_rpt200, etc.
% Use cycle number as key (not cycle string) to avoid invalid field name issue
capacity_ch09 = containers.Map();  % Use Map instead of struct for numeric keys
capacity_0cyc = NaN;  % Initial capacity for SOH calculation

for rptIdx = 1:length(rptCycles)
    rptCycle = rptCycles{rptIdx};
    cycleNum = rptCycleNums(rptIdx);
    
    % Get Ch09 static capacity (C/3 discharge capacity)
    capFieldName = sprintf('static_capacity_ch09_rpt%d', cycleNum);
    if isfield(OCV_data, capFieldName)
        capacity_ch09(num2str(cycleNum)) = OCV_data.(capFieldName);
        if cycleNum == 0
            capacity_0cyc = OCV_data.(capFieldName);
        end
    else
        capacity_ch09(num2str(cycleNum)) = NaN;
        fprintf('  Warning: %s not found in OCV_data\n', capFieldName);
    end
end

if isnan(capacity_0cyc)
    fprintf('  Error: Initial capacity (0cyc) not found. Cannot calculate SOH.\n');
else
    fprintf('  Ch09 Initial Capacity (0cyc): %.3f Ah\n', capacity_0cyc);
    
    % Create visualization for SOC70 only (충전과 방전을 각각 따로 표시)
    socLevel = 'SOC70';
    
    % Figure 1: Discharge Efficiency
    figTitle_Dis = sprintf('Discharge Efficiency vs SOH - %s\n(Ch09 C/3 Discharge Capacity)', socLevel);
    figName_Dis = sprintf('03_Discharge_Efficiency_vs_SOH_%s', socLevel);
    
    fig_Dis = figure('Name', figTitle_Dis, 'Position', [100, 100, 1800, 1000]);
    set(gcf, 'Visible', 'off');
    
    % Figure 2: Charge Efficiency
    figTitle_Chg = sprintf('Charge Efficiency vs SOH - %s\n(Ch09 C/3 Discharge Capacity)', socLevel);
    figName_Chg = sprintf('03_Charge_Efficiency_vs_SOH_%s', socLevel);
    
    fig_Chg = figure('Name', figTitle_Chg, 'Position', [100, 100, 1800, 1000]);
    set(gcf, 'Visible', 'off');
    
    % Create 8 subplots for DC1~DC8 for each figure
    fprintf('  Processing %d DC types...\n', length(dcTypes));
    for dcIdx = 1:length(dcTypes)
        dcType = dcTypes{dcIdx};
        fprintf('    Processing %s...\n', dcType);
        
        % Filter results for this DC and SOC level
        mask = strcmp(resultsTable.SOC_Level, socLevel) & ...
               strcmp(resultsTable.DC_Type, dcType);
        filteredData = resultsTable(mask, :);
        
        if height(filteredData) == 0
            % Discharge figure
            figure(fig_Dis);
            subplot(2, 4, dcIdx);
            text(0.5, 0.5, sprintf('%s\nNo Data', dcType), ...
                'HorizontalAlignment', 'center', 'FontSize', 12);
            axis off;
            
            % Charge figure
            figure(fig_Chg);
            subplot(2, 4, dcIdx);
            text(0.5, 0.5, sprintf('%s\nNo Data', dcType), ...
                'HorizontalAlignment', 'center', 'FontSize', 12);
            axis off;
            continue;
        end
        
        % Extract cycle numbers, both efficiencies, and calculate SOH
        % Vectorized operations for better performance
        cycleStrs = filteredData.Cycle;
        cycleNums = cellfun(@(x) str2double(strrep(x, 'cyc', '')), cycleStrs);
        
        % Get Discharge and Charge efficiencies (vectorized)
        if ismember('Discharge_Efficiency', filteredData.Properties.VariableNames)
            discharge_eff = filteredData.Discharge_Efficiency;
        else
            discharge_eff = NaN(height(filteredData), 1);
        end
        
        if ismember('Charge_Efficiency', filteredData.Properties.VariableNames)
            charge_eff = filteredData.Charge_Efficiency;
        else
            charge_eff = NaN(height(filteredData), 1);
        end
        
        % Calculate SOH: 현재 사이클 용량값 / 초기 사이클 용량값 (0사이클) × 100
        soh_values = NaN(height(filteredData), 1);
        for i = 1:height(filteredData)
            cycleKey = num2str(cycleNums(i));
            if capacity_ch09.isKey(cycleKey)
                cap_value = capacity_ch09(cycleKey);
                if ~isnan(cap_value) && ~isnan(capacity_0cyc)
                    soh_values(i) = (cap_value / capacity_0cyc) * 100;
                end
            end
        end
        
        % Get unique cycle numbers and aggregate efficiencies
        uniqueCycles = unique(cycleNums);
        discharge_eff_agg = zeros(length(uniqueCycles), 1);
        charge_eff_agg = zeros(length(uniqueCycles), 1);
        soh_agg = zeros(length(uniqueCycles), 1);
        
        for i = 1:length(uniqueCycles)
            cycleMask = cycleNums == uniqueCycles(i);
            discharge_eff_agg(i) = mean(discharge_eff(cycleMask), 'omitnan');
            charge_eff_agg(i) = mean(charge_eff(cycleMask), 'omitnan');
            soh_agg(i) = mean(soh_values(cycleMask), 'omitnan');
        end
        
        % Sort by cycle number
        [cycleNums_sorted, sortIdx] = sort(uniqueCycles);
        discharge_eff_sorted = discharge_eff_agg(sortIdx) * 100;  % Convert to percentage
        charge_eff_sorted = charge_eff_agg(sortIdx) * 100;  % Convert to percentage
        soh_sorted = soh_agg(sortIdx);
        
        % Set x-axis labels with SOH values: "0 (98%)", "200 (96%)", etc.
        xTickLabels = cell(length(cycleNums_sorted), 1);
        for i = 1:length(cycleNums_sorted)
            xTickLabels{i} = sprintf('%d (%.0f%%)', cycleNums_sorted(i), soh_sorted(i));
        end
        
        % Plot Discharge Efficiency
        figure(fig_Dis);
        subplot(2, 4, dcIdx);
        
        discharge_valid = ~isnan(discharge_eff_sorted) & ~isnan(soh_sorted);
        if sum(discharge_valid) > 0
            x_pos = 1:length(cycleNums_sorted);
            bar(x_pos(discharge_valid), discharge_eff_sorted(discharge_valid), ...
                'FaceColor', [0 0 1], 'EdgeColor', 'k', 'LineWidth', 1);
            set(gca, 'XTick', x_pos, 'XTickLabel', xTickLabels);
            xlabel('Cycle (SOH)', 'FontSize', 11);
            ylabel('Discharge Efficiency [%]', 'FontSize', 11);
            title(sprintf('%s', dcType), 'FontSize', 12, 'FontWeight', 'bold');
            grid on;
            grid minor;
            if length(cycleNums_sorted) > 0
                xlim([0.5, length(cycleNums_sorted) + 0.5]);
            end
        else
            text(0.5, 0.5, sprintf('%s\nNo Valid Data', dcType), ...
                'HorizontalAlignment', 'center', 'FontSize', 12);
            axis off;
        end
        
        % Plot Charge Efficiency
        figure(fig_Chg);
        subplot(2, 4, dcIdx);
        
        charge_valid = ~isnan(charge_eff_sorted) & ~isnan(soh_sorted);
        if sum(charge_valid) > 0
            x_pos = 1:length(cycleNums_sorted);
            bar(x_pos(charge_valid), charge_eff_sorted(charge_valid), ...
                'FaceColor', [0.85 0.33 0.1], 'EdgeColor', 'k', 'LineWidth', 1);
            set(gca, 'XTick', x_pos, 'XTickLabel', xTickLabels);
            xlabel('Cycle (SOH)', 'FontSize', 11);
            ylabel('Charge Efficiency [%]', 'FontSize', 11);
            title(sprintf('%s', dcType), 'FontSize', 12, 'FontWeight', 'bold');
            grid on;
            grid minor;
            if length(cycleNums_sorted) > 0
                xlim([0.5, length(cycleNums_sorted) + 0.5]);
            end
        else
            text(0.5, 0.5, sprintf('%s\nNo Valid Data', dcType), ...
                'HorizontalAlignment', 'center', 'FontSize', 12);
            axis off;
        end
    end
    
    % Save figures
    fprintf('  Saving Discharge Efficiency figure...\n');
    figure(fig_Dis);
    sgtitle(figTitle_Dis, 'FontSize', 14, 'FontWeight', 'bold');
    savefig(fig_Dis, fullfile(figDir, [figName_Dis, '.fig']));
    close(fig_Dis);
    fprintf('  Saved figure: %s\n', figName_Dis);
    
    fprintf('  Saving Charge Efficiency figure...\n');
    figure(fig_Chg);
    sgtitle(figTitle_Chg, 'FontSize', 14, 'FontWeight', 'bold');
    savefig(fig_Chg, fullfile(figDir, [figName_Chg, '.fig']));
    close(fig_Chg);
    fprintf('  Saved figure: %s\n', figName_Chg);
    fprintf('  Saved figure: %s\n', figName_Chg);
    
    %% Visualization 3-2: Net Efficiency vs SOH (New)
    fprintf('  Creating Net Efficiency visualization...\n');
    
    % Filter for Net Efficiency
    netMask = strcmp(resultsTable.SOC_Level, socLevel) & ...
              strcmp(resultsTable.DC_Type, 'Net_Efficiency');
    netData = resultsTable(netMask, :);
    
    if height(netData) > 0
        % Extract cycle info
        cycleStrs = netData.Cycle;
        cycleNums = cellfun(@(x) str2double(strrep(x, 'cyc', '')), cycleStrs);
        
        % Get efficiencies
        if ismember('Discharge_Efficiency', netData.Properties.VariableNames)
            net_dis_eff = netData.Discharge_Efficiency * 100;
        else
            net_dis_eff = NaN(height(netData), 1);
        end
        
        if ismember('Charge_Efficiency', netData.Properties.VariableNames)
            net_chg_eff = netData.Charge_Efficiency * 100;
        else
            net_chg_eff = NaN(height(netData), 1);
        end
        
        % Calculate SOH
        soh_values = NaN(height(netData), 1);
        for i = 1:height(netData)
            cycleKey = num2str(cycleNums(i));
            if capacity_ch09.isKey(cycleKey)
                cap_value = capacity_ch09(cycleKey);
                if ~isnan(cap_value) && ~isnan(capacity_0cyc)
                    soh_values(i) = (cap_value / capacity_0cyc) * 100;
                end
            end
        end
        
        % Sort by cycle number
        [cycleNums_sorted, sortIdx] = sort(cycleNums);
        net_dis_eff_sorted = net_dis_eff(sortIdx);
        net_chg_eff_sorted = net_chg_eff(sortIdx);
        soh_sorted = soh_values(sortIdx);
        
        % Create x-axis labels
        xTickLabels = cell(length(cycleNums_sorted), 1);
        for i = 1:length(cycleNums_sorted)
            xTickLabels{i} = sprintf('%d (%.0f%%)', cycleNums_sorted(i), soh_sorted(i));
        end
        
        % Plot Net Efficiency
        figName_Net = sprintf('03_2_Net_Efficiency_vs_SOH_%s', socLevel);
        fig_Net = figure('Name', 'Net Efficiency vs SOH', 'Position', [150, 150, 1000, 600]);
        set(gcf, 'Visible', 'off');
        
        x_pos = 1:length(cycleNums_sorted);
        b = bar(x_pos, [net_dis_eff_sorted, net_chg_eff_sorted]);
        b(1).FaceColor = [0 0 1];      % Discharge: Blue
        b(2).FaceColor = [0.85 0.33 0.1]; % Charge: Orange
        
        legend('Net Discharge Eff.', 'Net Charge Eff.', 'Location', 'best', 'FontSize', 10);
        set(gca, 'XTick', x_pos, 'XTickLabel', xTickLabels);
        xlabel('Cycle (SOH)', 'FontSize', 11);
        ylabel('Efficiency [%]', 'FontSize', 11);
        title(sprintf('Net Energy Efficiency vs SOH - %s\n(Aggregated DC1~DC8)', socLevel), ...
              'FontSize', 14, 'FontWeight', 'bold');
        ylim([80 110]);
        grid on;
        grid minor;
        
        % Add data labels
        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels1 = string(round(b(1).YData, 2)) + "%";
        text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontSize', 8);
         
        xtips2 = b(2).XEndPoints;
        ytips2 = b(2).YEndPoints;
        labels2 = string(round(b(2).YData, 2)) + "%";
        text(xtips2, ytips2, labels2, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontSize', 8);
             
        % Save figure
        savefig(fig_Net, fullfile(figDir, [figName_Net, '.fig']));      
        close(fig_Net);
        fprintf('  Saved figure: %s\n', figName_Net);
    else
        fprintf('  No Net Efficiency data found for visualization.\n');
    end
end

%% Visualization 4: Resistance vs Efficiency Scatter
% This visualization requires resistance data to be loaded separately
fprintf('\n=== Resistance vs Efficiency Visualization ===\n');

% Check if resistance data file exists (user should provide path)
resistanceDataPath = fullfile(outputDir, '..', '..', 'ver03_OnlyChargeEvents', 'Results', 'Resistance_Data.mat');
% Alternative path
if ~exist(resistanceDataPath, 'file')
    resistanceDataPath = fullfile(outputDir, 'Resistance_Data.mat');
end

if exist(resistanceDataPath, 'file')
    fprintf('Loading resistance data from: %s\n', resistanceDataPath);
    resistanceData = load(resistanceDataPath);
    
    % Try to find resistance table (adjust field name as needed)
    resistanceTable = [];
    if isfield(resistanceData, 'resistanceTable')
        resistanceTable = resistanceData.resistanceTable;
    elseif isfield(resistanceData, 'resultsTable')
        resistanceTable = resistanceData.resultsTable;
    else
        % Try to find any table variable
        vars = fieldnames(resistanceData);
        for v = 1:length(vars)
            if istable(resistanceData.(vars{v}))
                resistanceTable = resistanceData.(vars{v});
                fprintf('Found resistance table: %s\n', vars{v});
                break;
            end
        end
    end
    
    if ~isempty(resistanceTable)
        % Merge with efficiency results
        % Expected merge keys: Cycle, SOC_Level, DC_Type
        mergedTable = outerjoin(resultsTable, resistanceTable, ...
            'Keys', {'Cycle', 'SOC_Level', 'DC_Type'}, ...
            'MergeKeys', true);
        
        % Find resistance columns (columns containing 'Resistance' or 'R_')
        allCols = mergedTable.Properties.VariableNames;
        resistanceCols = {};
        for c = 1:length(allCols)
            if contains(allCols{c}, 'Resistance', 'IgnoreCase', true) || ...
               contains(allCols{c}, 'R_', 'IgnoreCase', false) || ...
               contains(allCols{c}, 'R3s', 'IgnoreCase', true) || ...
               contains(allCols{c}, 'R10s', 'IgnoreCase', true)
                resistanceCols{end+1} = allCols{c};
            end
        end
        
        if ~isempty(resistanceCols)
            fprintf('Found resistance columns: %s\n', strjoin(resistanceCols, ', '));
            
            % Create scatter plots for each resistance type
            for rIdx = 1:length(resistanceCols)
                rCol = resistanceCols{rIdx};
                
                % Filter valid data
                validIdx = ~isnan(mergedTable.Efficiency) & ~isnan(mergedTable.(rCol));
                if sum(validIdx) < 3
                    fprintf('  Skipping %s (insufficient data points)\n', rCol);
                    continue;
                end
                
                x_data = mergedTable.(rCol)(validIdx);
                y_data = mergedTable.Efficiency(validIdx) * 100;  % Convert to percentage
                
                % Create scatter plot
                fig4 = figure('Name', sprintf('Resistance vs Efficiency - %s', rCol), ...
                    'Position', [100, 100, 1000, 800]);
                set(gcf, 'Visible', 'off');
                
                scatter(x_data, y_data, 80, 'filled', 'MarkerFaceAlpha', 0.6);
                hold on;
                
                % Fit linear regression
                p = polyfit(x_data, y_data, 1);
                x_fit = linspace(min(x_data), max(x_data), 100);
                y_fit = polyval(p, x_fit);
                plot(x_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
                
                % Calculate R² and p-value
                [R, P] = corrcoef(x_data, y_data);
                R_squared = R(1,2)^2;
                p_value = P(1,2);
                
                % Determine resistance unit from column name
                if contains(rCol, 'mOhm', 'IgnoreCase', true) || contains(rCol, 'mΩ', 'IgnoreCase', true)
                    rUnit = '[mΩ]';
                elseif contains(rCol, 'Ohm', 'IgnoreCase', true) || contains(rCol, 'Ω', 'IgnoreCase', true)
                    rUnit = '[Ω]';
                else
                    rUnit = '[Ω]';  % Default unit
                end
                xlabel(sprintf('%s %s', strrep(rCol, '_', ' '), rUnit), 'FontSize', 12);
                ylabel('Energy Efficiency [%]', 'FontSize', 12);
                title(sprintf('Resistance vs Efficiency\nR² = %.4f, p-value = %.2e', ...
                    R_squared, p_value), 'FontSize', 14, 'FontWeight', 'bold');
                legend('Data Points', 'Linear Fit', 'Location', 'best', 'FontSize', 11);
                grid on;
                grid minor;
                hold off;
                
                fig4Name = sprintf('04_Resistance_Efficiency_%s', rCol);
                savefig(fig4, fullfile(figDir, [fig4Name, '.fig']));
                close(fig4);
                fprintf('  Saved figure: %s (R²=%.4f, p=%.2e)\n', fig4Name, R_squared, p_value);
            end
        else
            fprintf('No resistance columns found in the data\n');
        end
    else
        fprintf('Resistance table not found in the loaded data\n');
    end
else
    fprintf('Resistance data file not found. Skipping resistance visualization.\n');
    fprintf('To create this visualization, provide resistance data in Long Format table.\n');
    fprintf('Expected path: %s\n', resistanceDataPath);
end

%% Save results
fprintf('\n=== Saving Results ===\n');
fprintf('Total rows in results table: %d\n', height(resultsTable));

% Sort resultsTable by Cycle number (numerically)
cycleStrs = resultsTable.Cycle;
cycleNums = cellfun(@(x) str2double(strrep(x, 'cyc', '')), cycleStrs);
[~, sortIdx] = sort(cycleNums);
resultsTable = resultsTable(sortIdx, :);
fprintf('Sorted results table by Cycle number.\n');

% Save as MAT file
matSavePath = fullfile(outputDir, 'EnergyEfficiency_Results.mat');
save(matSavePath, 'resultsTable', 'I_threshold', 'Q_nominal', '-v7.3');
fprintf('Saved: %s\n', matSavePath);

% Save as Excel file (with units in column names)
excelSavePath = fullfile(outputDir, 'EnergyEfficiency_Results.xlsx');
% Create a copy of resultsTable with units in column names
resultsTableExcel = resultsTable;
% Check which columns exist
varNames = resultsTableExcel.Properties.VariableNames;
newVarNames = varNames;
for i = 1:length(varNames)
    switch varNames{i}
        case 'E_stored_ini'
            newVarNames{i} = 'E_stored_ini [Wh]';
        case 'E_stored_aft'
            newVarNames{i} = 'E_stored_aft [Wh]';
        case 'delta_E_stored'
            newVarNames{i} = 'delta_E_stored [Wh]';
        case 'E_input'
            newVarNames{i} = 'E_input [Wh]';
        case 'E_output'
            newVarNames{i} = 'E_output [Wh]';
        case 'E_net'
            newVarNames{i} = 'E_net [Wh]';
        case 'Discharge_Efficiency'
            newVarNames{i} = 'Discharge_Efficiency [-]';
        case 'Charge_Efficiency'
            newVarNames{i} = 'Charge_Efficiency [-]';
        case 'Efficiency'
            newVarNames{i} = 'Efficiency [-]';
    end
end
resultsTableExcel.Properties.VariableNames = newVarNames;

% Try to delete existing file if it exists (might be locked)
if exist(excelSavePath, 'file')
    try
        delete(excelSavePath);
        pause(0.5);  % Brief pause to ensure file is released
    catch
        fprintf('Warning: Could not delete existing Excel file. It may be open in another application.\n');
    end
end

% Try to write Excel file
try
    writetable(resultsTableExcel, excelSavePath);
    fprintf('Saved: %s\n', excelSavePath);
catch ME
    fprintf('Warning: Could not write Excel file: %s\n', ME.message);
    fprintf('Attempting to save as CSV instead...\n');
    % Fallback: Save as CSV file
    csvSavePath = fullfile(outputDir, 'EnergyEfficiency_Results.csv');
    writetable(resultsTableExcel, csvSavePath);
    fprintf('Saved as CSV: %s\n', csvSavePath);
end

% Display summary statistics
fprintf('\n=== Summary Statistics ===\n');
fprintf('Efficiency statistics:\n');
fprintf('  Mean: %.4f (%.2f%%)\n', mean(resultsTable.Efficiency, 'omitnan'), mean(resultsTable.Efficiency, 'omitnan')*100);
fprintf('  Std:  %.4f (%.2f%%)\n', std(resultsTable.Efficiency, 'omitnan'), std(resultsTable.Efficiency, 'omitnan')*100);
fprintf('  Min:  %.4f (%.2f%%)\n', min(resultsTable.Efficiency, [], 'omitnan'), min(resultsTable.Efficiency, [], 'omitnan')*100);
fprintf('  Max:  %.4f (%.2f%%)\n', max(resultsTable.Efficiency, [], 'omitnan'), max(resultsTable.Efficiency, [], 'omitnan')*100);

fprintf('\n=== Energy Efficiency Calculation Completed ===\n');
