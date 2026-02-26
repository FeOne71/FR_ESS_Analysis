classdef app1 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        VisualizationCheckBox_3        matlab.ui.control.CheckBox
        VisualizationCheckBox_2        matlab.ui.control.CheckBox
        VisualizationCheckBox          matlab.ui.control.CheckBox
        Gauge_4                        matlab.ui.control.LinearGauge
        Gauge_3                        matlab.ui.control.LinearGauge
        DescriptionTextArea            matlab.ui.control.TextArea
        DescriptionTextAreaLabel       matlab.ui.control.Label
        ProgressTextArea               matlab.ui.control.TextArea
        ProgressTextAreaLabel          matlab.ui.control.Label
        NormalizationCheckBox          matlab.ui.control.CheckBox
        StandardizationCheckBox        matlab.ui.control.CheckBox
        FineTuningOnCheckBox           matlab.ui.control.CheckBox
        StatusLamp_3                   matlab.ui.control.Lamp
        StatusLamp_3Label              matlab.ui.control.Label
        Tree_3                         matlab.ui.container.CheckBoxTree
        DynamicNode                    matlab.ui.container.TreeNode
        R1secNode                      matlab.ui.container.TreeNode
        R3secNode                      matlab.ui.container.TreeNode
        R5secNode                      matlab.ui.container.TreeNode
        R10secNode                     matlab.ui.container.TreeNode
        R30secNode                     matlab.ui.container.TreeNode
        DCIRChgDchgNode_2              matlab.ui.container.TreeNode
        EquilibriumNode                matlab.ui.container.TreeNode
        EnergyEfficientyNode           matlab.ui.container.TreeNode
        DriveCycle2ShallowNode_2       matlab.ui.container.TreeNode
        DriveCycle3ShallowNode_3       matlab.ui.container.TreeNode
        DriveCycle4ShallowNode_2       matlab.ui.container.TreeNode
        DriveCycle5ShallowNode_2       matlab.ui.container.TreeNode
        DriveCycle6ShallowNode_2       matlab.ui.container.TreeNode
        DriveCycle7ShallowNode_2       matlab.ui.container.TreeNode
        DriveCycle8ShallowNode_2       matlab.ui.container.TreeNode
        LabelNode                      matlab.ui.container.TreeNode
        AvailableCapacityAhNode        matlab.ui.container.TreeNode
        EmergyEfficiencyNode           matlab.ui.container.TreeNode
        AgingModeNode                  matlab.ui.container.TreeNode
        LLINode                        matlab.ui.container.TreeNode
        LAMpNode                       matlab.ui.container.TreeNode
        Gauge_2                        matlab.ui.control.LinearGauge
        Gauge                          matlab.ui.control.LinearGauge
        Tree_2                         matlab.ui.container.CheckBoxTree
        RandomForestNode               matlab.ui.container.TreeNode
        SupportVectorMachineNode       matlab.ui.container.TreeNode
        LSTMNode                       matlab.ui.container.TreeNode
        GaussianProcessRegressionNode  matlab.ui.container.TreeNode
        TransferLearningNode_2         matlab.ui.container.TreeNode
        Tree                           matlab.ui.container.CheckBoxTree
        RPTNode                        matlab.ui.container.TreeNode
        DriveCycleNode                 matlab.ui.container.TreeNode
        StatusLamp_2                   matlab.ui.control.Lamp
        StatusLamp_2Label              matlab.ui.control.Label
        DeploymentButton               matlab.ui.control.StateButton
        ModelButton                    matlab.ui.control.StateButton
        FeaturesLabelsButton           matlab.ui.control.StateButton
        StatusLamp                     matlab.ui.control.Lamp
        StatusLampLabel                matlab.ui.control.Label
        DataLoadButton                 matlab.ui.control.StateButton
        
        % ==============================================
        % [NEW] 앱 데이터 저장용 커스텀 변수들
        % ==============================================
        AppPath_RPT = ''                 % RPT 폴더 경로
        AppPath_DC = ''                  % Drive Cycle 폴더 경로
        LoadRPT = false                  % RPT 로드 여부 (Tree 체크 결과)
        LoadDC = false                   % Drive Cycle 로드 여부 (Tree 체크 결과)
        NumRPTFiles = 0                  % 스캔된 RPT 파일 개수
        NumDCFiles = 0                   % 스캔된 DC 파일 개수
        ProcessData_VQ = struct()        % DataLoader가 넘겨줄 스무딩/보간된 데이터
        FeatureTable = table()           % FeatureExtractor가 만들어낸 피처/라벨 최종 테이블
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Callbacks that handle component events
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        % ============================================================
        % [Callback 1] DataLoadButton → 파일 스캔 전용
        %   - Tree 체크박스에서 선택된 항목 확인
        %   - 선택된 폴더의 파일 존재 여부 & 개수만 스캔
        %   - 실제 전처리는 수행하지 않음
        % ============================================================
        function DataLoadButtonValueChanged(app, event)
            addpath('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\App');

            % 0. 로그 초기화
            app.ProgressTextArea.Value = {'=== Data Load Started ==='};
            drawnow limitrate;
            
            % 1. Tree 체크박스에서 체크된 항목 읽기
            checkedNodes = app.Tree.CheckedNodes;
            if isempty(checkedNodes)
                uialert(app.UIFigure, 'Please check the data items to load (RPT / Drive Cycle).', 'No Selection');
                app.DataLoadButton.Value = 0;
                return;
            end
            
            app.LoadRPT = false;
            app.LoadDC = false;
            for i = 1:numel(checkedNodes)
                nodeName = checkedNodes(i).Text;
                if strcmp(nodeName, 'RPT')
                    app.LoadRPT = true;
                elseif strcmp(nodeName, 'Drive Cycle')
                    app.LoadDC = true;
                end
            end
            app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('  Selected: RPT=%s, DC=%s', mat2str(app.LoadRPT), mat2str(app.LoadDC))}]; drawnow limitrate;
            
            % 2. 폴더 경로 선택
            if app.LoadRPT
                app.AppPath_RPT = uigetdir('D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data', 'Select RPT Data Folder');
                if app.AppPath_RPT == 0
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [Canceled] RPT folder selection canceled'}]; drawnow limitrate;
                    app.LoadRPT = false;
                else
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('  RPT Path: %s', app.AppPath_RPT)}]; drawnow limitrate;
                end
            end
            
            if app.LoadDC
                app.AppPath_DC = uigetdir('D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data', 'Select Drive Cycle Data Folder');
                if app.AppPath_DC == 0
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [Canceled] DC folder selection canceled'}]; drawnow limitrate;
                    app.LoadDC = false;
                    app.AppPath_DC = '';
                else
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('  DC Path: %s', app.AppPath_DC)}]; drawnow limitrate;
                end
            else
                app.AppPath_DC = '';
            end
            
            if ~app.LoadRPT && ~app.LoadDC
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [Warning] No valid folders selected.'}]; drawnow limitrate;
                app.StatusLamp.Color = [1.00, 0.00, 0.00];
                return;
            end
            
            % 3. 상태등 노란색(스캔 중)
            app.StatusLamp.Color = [1.00, 0.84, 0.00];
            drawnow;
            
            % 4. 파일 스캔 (존재 여부 + 개수만 확인)
            app.NumRPTFiles = 0;
            app.NumDCFiles = 0;
            
            if app.LoadRPT && ~isempty(app.AppPath_RPT)
                rpt_files = dir(fullfile(app.AppPath_RPT, '**\*_RPT_*.csv'));
                app.NumRPTFiles = length(rpt_files);
                if app.NumRPTFiles > 0
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('  [OK] Found %d RPT files', app.NumRPTFiles)}]; drawnow limitrate;
                else
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [X] No RPT files found (*_RPT_*.csv)'}]; drawnow limitrate;
                end
            end
            
            if app.LoadDC && ~isempty(app.AppPath_DC)
                dc_mat = dir(fullfile(app.AppPath_DC, '**\*.mat'));
                dc_csv = dir(fullfile(app.AppPath_DC, '**\*Drive*.csv'));
                app.NumDCFiles = length(dc_mat) + length(dc_csv);
                if app.NumDCFiles > 0
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('  [OK] Found %d DC files (mat:%d, csv:%d)', app.NumDCFiles, length(dc_mat), length(dc_csv))}]; drawnow limitrate;
                else
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [X] No Drive Cycle files found'}]; drawnow limitrate;
                end
            end
            
            % 5. 결과 표시
            if (app.NumRPTFiles + app.NumDCFiles) > 0
                app.StatusLamp.Color = [0.00, 1.00, 0.00]; % Green
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('=== Scan Complete: Total %d files ===', app.NumRPTFiles + app.NumDCFiles)}];
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'>> Please click [Features & Labels] button.'}]; drawnow limitrate;
            else
                app.StatusLamp.Color = [1.00, 0.00, 0.00];
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'=== Files Not Found ==='}]; drawnow limitrate;
            end
        end

        % [Callback 2] FeaturesLabelsButton
        % ============================================================
        function FeaturesLabelsButtonValueChanged(app, event)
            addpath('D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\App');

            % HTML progress bar template
            barTpl = ['<style>*{margin:0;padding:0;box-sizing:border-box;font-family:Arial,sans-serif;}' ...
              '.title{text-align:center;font-size:10px;color:#999;letter-spacing:2px;margin-bottom:4px;}' ...
              '.row{display:flex;align-items:center;gap:12px;padding:0 4px;}' ...
              '.track{flex:1;background:#d9d9d9;height:30px;border-radius:2px;}' ...
              '.fill{height:100%%;background:#6cc0f5;border-radius:2px;}' ...
              '.val{font-size:32px;font-weight:bold;color:#333;min-width:60px;text-align:right;}' ...
              '.ticks{display:flex;justify-content:space-between;padding:4px 4px 0;font-size:10px;color:#999;}</style>' ...
              '<div class="title">PROGRESS</div>' ...
              '<div class="row"><div class="track"><div class="fill" style="width:%d%%"></div></div>' ...
              '<div class="val">%03d</div></div>' ...
              '<div class="ticks" style="margin-right:72px;">' ...
              '<span>0</span><span>20</span><span>40</span><span>60</span><span>80</span><span>100</span></div>'];

            if (app.NumRPTFiles + app.NumDCFiles) == 0
                uialert(app.UIFigure, 'No scanned files. Run [Data Load] first.', 'Error');
                app.FeaturesLabelsButton.Value = 0;
                return;
            end
            
            app.StatusLamp_3.Color = [1.00, 0.84, 0.00];
            app.HTML.HTMLSource = sprintf(barTpl, 0, 0); drawnow;
            app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {''}; {'=== Features & Labels Start ==='}]; drawnow;
            
            try
                % ====================================================
                % Step 1: Process RPT files one by one (0% ~ 50%)
                % ====================================================
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [1/2] Processing RPT files...'}]; drawnow;
                
                App_VQ_grid = struct();
                
                if app.LoadRPT && ~isempty(app.AppPath_RPT)
                    csv_files = dir(fullfile(app.AppPath_RPT, '**\*_RPT_*.csv'));
                    nFiles = length(csv_files);
                    
                    for k = 1:nFiles
                        fn = csv_files(k).name;
                        fp = fullfile(csv_files(k).folder, fn);
                        
                        tok = regexp(fn, 'Ch(\d+)_RPT_(\d+)cyc', 'tokens');
                        if isempty(tok), continue; end
                        
                        ch_key = sprintf('Ch%02d', str2double(tok{1}{1}));
                        cyc_key = sprintf('cyc%d', str2double(tok{1}{2}));
                        
                        ch_data = process_single_rpt_file(fp, 3, [8, 10], 20);
                        
                        if ~isfield(App_VQ_grid, cyc_key)
                            App_VQ_grid.(cyc_key) = struct();
                        end
                        App_VQ_grid.(cyc_key).(ch_key) = ch_data;
                        
                        % Update progress bar after each file (0%~50%)
                        pct = round((k / nFiles) * 50);
                        app.HTML.HTMLSource = sprintf(barTpl, pct, pct);
                        drawnow;
                    end
                end
                
                app.ProcessData_VQ = App_VQ_grid;
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [1/2] RPT processing done!'}]; drawnow;
                
                % ====================================================
                % Step 2: Feature Extraction (50% ~ 100%)
                % ====================================================
                app.HTML.HTMLSource = sprintf(barTpl, 50, 50); drawnow;
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [2/2] Extracting features (dQ/dV, Peak)...'}]; drawnow;
                
                master_ruler_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\Dataset\MasterRulers.mat';
                app.FeatureTable = App_FeatureExtractor(app.ProcessData_VQ, master_ruler_path);
                
                app.HTML.HTMLSource = sprintf(barTpl, 100, 100); drawnow;
                
                if isempty(app.FeatureTable)
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [Warning] No features extracted'}]; drawnow;
                    app.StatusLamp_3.Color = [1.00, 0.00, 0.00];
                    app.FeaturesLabelsButton.Value = 0;
                    return;
                end
                
                app.StatusLamp_3.Color = [0.00, 1.00, 0.00];
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('  [2/2] Done! (%dx%d)', size(app.FeatureTable,1), size(app.FeatureTable,2))}];
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'=== Features & Labels Complete ==='}]; drawnow;
                
                % ====================================================
                % Visualization (if checked)
                % ====================================================
                visDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
                
                % Data Load Visualization
                if app.VisualizationCheckBox.Value
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [VIS] Running Data Load visualization...'}]; drawnow;
                    try
                        App_Visualizer_DataLoad(app.ProcessData_VQ, visDir);
                    catch ME_vis
                        app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('  [VIS ERROR] %s', ME_vis.message)}]; drawnow;
                    end
                end
                
                % Features & Labels Visualization
                if app.VisualizationCheckBox_2.Value
                    app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {'  [VIS] Running Features & Labels visualization...'}]; drawnow;
                    try
                        % Collect checked Tree_3 node names
                        checkedNodes_T3 = {};
                        if ~isempty(app.Tree_3.CheckedNodes)
                            for nn = 1:numel(app.Tree_3.CheckedNodes)
                                checkedNodes_T3{end+1} = app.Tree_3.CheckedNodes(nn).Text;
                            end
                        end
                        App_Visualizer_Features(app.ProcessData_VQ, app.FeatureTable, checkedNodes_T3, visDir);
                    catch ME_vis
                        app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('  [VIS ERROR] %s', ME_vis.message)}]; drawnow;
                    end
                end
                
            catch ME
                app.StatusLamp_3.Color = [1.00, 0.00, 0.00];
                app.HTML.HTMLSource = sprintf(barTpl, 0, 0); drawnow;
                app.ProgressTextArea.Value = [app.ProgressTextArea.Value; {sprintf('  [ERROR] %s', ME.message)}]; drawnow;
                app.FeaturesLabelsButton.Value = 0;
            end
        end

    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Component initialization (UI 생성)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        % ============================================================
        % createComponents: 모든 UI 컴포넌트 생성
        % ============================================================
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 678 560];
            app.UIFigure.Name = 'MATLAB App';

            % Create DataLoadButton
            app.DataLoadButton = uibutton(app.UIFigure, 'state');
            app.DataLoadButton.ValueChangedFcn = createCallbackFcn(app, @DataLoadButtonValueChanged, true);
            app.DataLoadButton.Text = 'Data Load';
            app.DataLoadButton.FontWeight = 'bold';
            app.DataLoadButton.Position = [42 447 100 23];

            % Create StatusLampLabel
            app.StatusLampLabel = uilabel(app.UIFigure);
            app.StatusLampLabel.HorizontalAlignment = 'right';
            app.StatusLampLabel.Position = [403 449 39 22];
            app.StatusLampLabel.Text = 'Status';

            % Create StatusLamp
            app.StatusLamp = uilamp(app.UIFigure);
            app.StatusLamp.Position = [457 449 20 20];

            % Create FeaturesLabelsButton
            app.FeaturesLabelsButton = uibutton(app.UIFigure, 'state');
            % [MODIFIED] Added Callback for Features & Labels Button
            app.FeaturesLabelsButton.ValueChangedFcn = createCallbackFcn(app, @FeaturesLabelsButtonValueChanged, true);
            app.FeaturesLabelsButton.Text = {'Features &'; 'Labels'};
            app.FeaturesLabelsButton.FontWeight = 'bold';
            app.FeaturesLabelsButton.Position = [42 343 100 38];

            % Create ModelButton
            app.ModelButton = uibutton(app.UIFigure, 'state');
            app.ModelButton.Text = 'Model';
            app.ModelButton.FontWeight = 'bold';
            app.ModelButton.Position = [42 222 100 23];

            % Create DeploymentButton
            app.DeploymentButton = uibutton(app.UIFigure, 'state');
            app.DeploymentButton.Text = 'Deployment';
            app.DeploymentButton.FontWeight = 'bold';
            app.DeploymentButton.Position = [42 61 100 23];

            % Create StatusLamp_2Label
            app.StatusLamp_2Label = uilabel(app.UIFigure);
            app.StatusLamp_2Label.HorizontalAlignment = 'right';
            app.StatusLamp_2Label.Position = [181 61 39 22];
            app.StatusLamp_2Label.Text = 'Status';

            % Create StatusLamp_2
            app.StatusLamp_2 = uilamp(app.UIFigure);
            app.StatusLamp_2.Position = [235 61 20 20];

            % Create Tree (체크박스 트리 - RPT와 Drive Cycle 최상위 노드만)
            app.Tree = uitree(app.UIFigure, 'checkbox');
            app.Tree.Position = [181 430 204 40];

            % Create RPTNode
            app.RPTNode = uitreenode(app.Tree);
            app.RPTNode.Text = 'RPT';

            % Create DriveCycleNode
            app.DriveCycleNode = uitreenode(app.Tree);
            app.DriveCycleNode.Text = 'Drive Cycle';

            % Create Tree_2
            app.Tree_2 = uitree(app.UIFigure, 'checkbox');
            app.Tree_2.Position = [181 139 204 106];

            % Create RandomForestNode
            app.RandomForestNode = uitreenode(app.Tree_2);
            app.RandomForestNode.Text = 'Random Forest';

            % Create SupportVectorMachineNode
            app.SupportVectorMachineNode = uitreenode(app.Tree_2);
            app.SupportVectorMachineNode.Text = 'Support Vector Machine';

            % Create LSTMNode
            app.LSTMNode = uitreenode(app.Tree_2);
            app.LSTMNode.Text = 'LSTM';

            % Create GaussianProcessRegressionNode
            app.GaussianProcessRegressionNode = uitreenode(app.Tree_2);
            app.GaussianProcessRegressionNode.Text = 'Gaussian Process Regression';

            % Create TransferLearningNode_2
            app.TransferLearningNode_2 = uitreenode(app.Tree_2);
            app.TransferLearningNode_2.Text = 'Transfer Learning';

            % Create Gauge
            app.Gauge = uigauge(app.UIFigure, 'linear');
            app.Gauge.Position = [403 139 119 41];

            % Create Gauge_2
            app.Gauge_2 = uigauge(app.UIFigure, 'linear');
            app.Gauge_2.Position = [403 55 119 41];

            % Create Tree_3
            app.Tree_3 = uitree(app.UIFigure, 'checkbox');
            app.Tree_3.Position = [181 301 204 80];

            % Create DynamicNode
            app.DynamicNode = uitreenode(app.Tree_3);
            app.DynamicNode.Text = 'Dynamic';

            % Create R1secNode
            app.R1secNode = uitreenode(app.DynamicNode);
            app.R1secNode.Text = 'R1sec';

            % Create R3secNode
            app.R3secNode = uitreenode(app.DynamicNode);
            app.R3secNode.Text = 'R3sec';

            % Create R5secNode
            app.R5secNode = uitreenode(app.DynamicNode);
            app.R5secNode.Text = 'R5sec';

            % Create R10secNode
            app.R10secNode = uitreenode(app.DynamicNode);
            app.R10secNode.Text = 'R10sec';

            % Create R30secNode
            app.R30secNode = uitreenode(app.DynamicNode);
            app.R30secNode.Text = 'R30sec';

            % Create DCIRChgDchgNode_2
            app.DCIRChgDchgNode_2 = uitreenode(app.DynamicNode);
            app.DCIRChgDchgNode_2.Text = 'DCIR (Chg, Dchg)';

            % Create EquilibriumNode
            app.EquilibriumNode = uitreenode(app.Tree_3);
            app.EquilibriumNode.Text = 'Equilibrium';

            % Create EnergyEfficientyNode
            app.EnergyEfficientyNode = uitreenode(app.EquilibriumNode);
            app.EnergyEfficientyNode.Text = 'Energy Efficienty';

            % Create DriveCycle2ShallowNode_2
            app.DriveCycle2ShallowNode_2 = uitreenode(app.EquilibriumNode);
            app.DriveCycle2ShallowNode_2.Text = 'Drive Cycle 2 (Shallow)';

            % Create DriveCycle3ShallowNode_3
            app.DriveCycle3ShallowNode_3 = uitreenode(app.EquilibriumNode);
            app.DriveCycle3ShallowNode_3.Text = 'Drive Cycle 3 (Shallow)';

            % Create DriveCycle4ShallowNode_2
            app.DriveCycle4ShallowNode_2 = uitreenode(app.EquilibriumNode);
            app.DriveCycle4ShallowNode_2.Text = 'Drive Cycle 4 (Shallow)';

            % Create DriveCycle5ShallowNode_2
            app.DriveCycle5ShallowNode_2 = uitreenode(app.EquilibriumNode);
            app.DriveCycle5ShallowNode_2.Text = 'Drive Cycle 5 (Shallow)';

            % Create DriveCycle6ShallowNode_2
            app.DriveCycle6ShallowNode_2 = uitreenode(app.EquilibriumNode);
            app.DriveCycle6ShallowNode_2.Text = 'Drive Cycle 6 (Shallow)';

            % Create DriveCycle7ShallowNode_2
            app.DriveCycle7ShallowNode_2 = uitreenode(app.EquilibriumNode);
            app.DriveCycle7ShallowNode_2.Text = 'Drive Cycle 7 (Shallow)';

            % Create DriveCycle8ShallowNode_2
            app.DriveCycle8ShallowNode_2 = uitreenode(app.EquilibriumNode);
            app.DriveCycle8ShallowNode_2.Text = 'Drive Cycle 8 (Shallow)';

            % Create LabelNode
            app.LabelNode = uitreenode(app.Tree_3);
            app.LabelNode.Text = 'Label';

            % Create AvailableCapacityAhNode
            app.AvailableCapacityAhNode = uitreenode(app.LabelNode);
            app.AvailableCapacityAhNode.Text = 'Available Capacity (Ah)';

            % Create EmergyEfficiencyNode
            app.EmergyEfficiencyNode = uitreenode(app.LabelNode);
            app.EmergyEfficiencyNode.Text = 'Emergy Efficiency';

            % Create AgingModeNode
            app.AgingModeNode = uitreenode(app.LabelNode);
            app.AgingModeNode.Text = 'Aging Mode';

            % Create LLINode
            app.LLINode = uitreenode(app.AgingModeNode);
            app.LLINode.Text = 'LLI';

            % Create LAMpNode
            app.LAMpNode = uitreenode(app.AgingModeNode);
            app.LAMpNode.Text = 'LAMp';

            % Create StatusLamp_3Label
            app.StatusLamp_3Label = uilabel(app.UIFigure);
            app.StatusLamp_3Label.HorizontalAlignment = 'right';
            app.StatusLamp_3Label.Position = [401 349 39 22];
            app.StatusLamp_3Label.Text = 'Status';

            % Create StatusLamp_3
            app.StatusLamp_3 = uilamp(app.UIFigure);
            app.StatusLamp_3.Position = [455 349 20 20];

            % Create FineTuningOnCheckBox
            app.FineTuningOnCheckBox = uicheckbox(app.UIFigure);
            app.FineTuningOnCheckBox.Text = 'Fine Tuning On';
            app.FineTuningOnCheckBox.Position = [403 222 217 22];
            app.FineTuningOnCheckBox.Value = true;

            % Create StandardizationCheckBox
            app.StandardizationCheckBox = uicheckbox(app.UIFigure);
            app.StandardizationCheckBox.Text = 'Standardization';
            app.StandardizationCheckBox.Position = [403 193 105 22];

            % Create NormalizationCheckBox
            app.NormalizationCheckBox = uicheckbox(app.UIFigure);
            app.NormalizationCheckBox.Text = 'Normalization';
            app.NormalizationCheckBox.Position = [549 193 95 22];

            % Create DescriptionTextAreaLabel
            app.DescriptionTextAreaLabel = uilabel(app.UIFigure);
            app.DescriptionTextAreaLabel.HorizontalAlignment = 'right';
            app.DescriptionTextAreaLabel.FontWeight = 'bold';
            app.DescriptionTextAreaLabel.Position = [42 527 71 22];
            app.DescriptionTextAreaLabel.Text = 'Description';

            % Create DescriptionTextArea
            app.DescriptionTextArea = uitextarea(app.UIFigure);
            app.DescriptionTextArea.FontSize = 10;
            app.DescriptionTextArea.Position = [128 514 513 37];
            app.DescriptionTextArea.Value = {'This application is for state estimator for only using Lab Dataset (RPT and Drive Cycle). You can choose the dataset option using dropdown menu, features and labels, and also ML / DL options.As the model development has finished, you can delploy your model. '};

            % Create Gauge_3
            app.Gauge_3 = uigauge(app.UIFigure, 'linear');
            app.Gauge_3.Position = [401 400 119 41];

            % Create Gauge_4
            app.Gauge_4 = uigauge(app.UIFigure, 'linear');
            app.Gauge_4.Position = [401 300 119 41];

            % Create ProgressTextAreaLabel
            app.ProgressTextAreaLabel = uilabel(app.UIFigure);
            app.ProgressTextAreaLabel.FontWeight = 'bold';
            app.ProgressTextAreaLabel.Position = [490 475 80 22];
            app.ProgressTextAreaLabel.Text = 'Progress Status';

            % Create ProgressTextArea (진행 상황 실시간 로그)
            app.ProgressTextArea = uitextarea(app.UIFigure);
            app.ProgressTextArea.Editable = 'off';
            app.ProgressTextArea.FontName = 'Consolas';
            app.ProgressTextArea.FontSize = 10;
            app.ProgressTextArea.BackgroundColor = [0.15 0.15 0.18];
            app.ProgressTextArea.FontColor = [0.0 1.0 0.4];
            app.ProgressTextArea.Position = [490 300 160 175];
            app.ProgressTextArea.Value = {'Waiting...'};

            % Create VisualizationCheckBox
            app.VisualizationCheckBox = uicheckbox(app.UIFigure);
            app.VisualizationCheckBox.Text = 'Visualization';
            app.VisualizationCheckBox.Position = [549 222 89 22];

            % Create VisualizationCheckBox_2
            app.VisualizationCheckBox_2 = uicheckbox(app.UIFigure);
            app.VisualizationCheckBox_2.Text = 'Visualization';
            app.VisualizationCheckBox_2.Position = [549 349 89 22];

            % Create VisualizationCheckBox_3
            app.VisualizationCheckBox_3 = uicheckbox(app.UIFigure);
            app.VisualizationCheckBox_3.Text = 'Visualization';
            app.VisualizationCheckBox_3.Position = [549 449 89 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
