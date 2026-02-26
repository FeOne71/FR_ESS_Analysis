% =========================================================================
% 이 코드는 복사-붙여넣기 전용 텍스트 파일입니다.
% MATLAB App Designer에서 app1.mlapp를 여신 후, '코드 보기(Code View)' 화면에서
% 해당하는 버튼을 눌러 속성(Property) 및 콜백(Callback) 함수를 추가하고 
% 아래의 코드를 붙여넣으세요.
% =========================================================================

% -------------------------------------------------------------------------
% 1. 앱 속성(Property) 추가하기
% -------------------------------------------------------------------------
% 1) 앱 디자이너 상단 리본 메뉴에서 [Property] -> [Public Property] 클릭
% 2) 생성된 'properties (Access = public)' 블록 안에 아래 4줄을 붙여넣으세요.

        AppPath_RPT = ''                 % RPT 폴더 경로
        AppPath_DC = ''                  % Drive Cycle 폴더 경로
        ProcessData_VQ = struct()        % DataLoader가 넘겨줄 스무딩/보간된 데이터
        FeatureTable = table()           % FeatureExtractor가 만들어낸 피처/라벨 최종 테이블


% -------------------------------------------------------------------------
% 2. DataLoadButton 콜백 함수 (수정)
% -------------------------------------------------------------------------
% 1) 앱 디자이너 디자인 뷰에서 'Data Load' 버튼 우클릭
% 2) Callbacks -> 'DataLoadButtonValueChanged' (또는 'ButtonPushed') 클릭
% 3) 자동으로 생성된 함수 블록의 내용을 아래 코드로 모두 덮어쓰세요.
%    (만약 기존 함수 이름이 다르다면, 함수 내부 내용만 복사하세요)

        function DataLoadButtonValueChanged(app, event)
            % 1. 경로 선택 창 띄우기 (폴더 2개 선택)
            app.AppPath_RPT = uigetdir('D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data', 'Select RPT Data Folder');
            if app.AppPath_RPT == 0, return; end  % 취소 누르면 종료
            
            app.AppPath_DC = uigetdir('D:\JCW\Projects\KEPCO_ESS_Local\Experimental Data', 'Select Drive Cycle Data Folder (Cancel to skip)');
            if app.AppPath_DC == 0, app.AppPath_DC = ''; end % DC는 없을 수도 있으니 빈칸 허용
            
            % 2. 상태등 노란색(로딩 중)으로 변경
            app.StatusLamp.Color = [1.00, 0.84, 0.00]; % Yellow
            drawnow;
            
            try
                % 게이지 초기화
                app.Gauge_4.Value = 0;
                
                % 3. 방금 만든 백엔드 스크립트 실행 (RPT폴더, DC폴더 윈도우20, 병렬false, 진행률 게이지 전달)
                app.ProcessData_VQ = App_DataLoader(app.AppPath_RPT, app.AppPath_DC, 20, 20, false, app.Gauge_4);
                
                % 4. 데이터 로딩 성공 시 상태등 녹색으로 변경 및 알림
                app.Gauge_4.Value = 100;
                app.StatusLamp.Color = [0.00, 1.00, 0.00]; % Green
                uialert(app.UIFigure, 'Data interpolation and moving average (Window=20) completed successfully!', 'Data Load Success');
            catch ME
                % 에러 발생 시 상태등 빨간색
                app.Gauge_4.Value = 0;
                app.StatusLamp.Color = [1.00, 0.00, 0.00]; % Red
                uialert(app.UIFigure, ['An error occurred during data processing: ', ME.message], 'Error');
            end
        end


% -------------------------------------------------------------------------
% 3. Features & Labels 버튼 콜백 함수 (신규)
% -------------------------------------------------------------------------
% 1) 앱 디자이너 디자인 뷰에서 'Features & Labels' 버튼을 찾아 우클릭
% 2) Callbacks -> 'Add ButtonPushedFcn' 클릭
% 3) 자동으로 생성된 함수 블록 내부에 아래 코드를 덮어쓰세요.
%    (만약 이름이 'FeaturesLabelsButtonPushed' 와 다르면 내부 내용만 복사하세요)

        function FeaturesLabelsButtonPushed(app, event)
            % 데이터 로드가 안되어 있으면 경고
            if isempty(fieldnames(app.ProcessData_VQ))
                uialert(app.UIFigure, 'No interpolated data. Please run [Data Load] first.', 'Error');
                return;
            end
            
            % 상태등 노란색(계산 중)으로 변경
            app.StatusLamp_3.Color = [1.00, 0.84, 0.00]; % Yellow
            drawnow;
            
            try
                % 1. MasterRulers 경로 지정 (고정)
                master_ruler_path = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\Lab_RPT_Analysis\Dataset\MasterRulers.mat';
                
                % 2. 피처 추출 및 라벨링 스크립트 실행
                app.FeatureTable = App_FeatureExtractor(app.ProcessData_VQ, master_ruler_path);
                
                if isempty(app.FeatureTable)
                    uialert(app.UIFigure, 'No data extracted.', 'Warning');
                    app.StatusLamp_3.Color = [1.00, 0.00, 0.00];
                    return;
                end
                
                % 성공 시 상태등 녹색으로 변경
                app.StatusLamp_3.Color = [0.00, 1.00, 0.00]; % Green
                uialert(app.UIFigure, 'Features and labels extracted successfully based on Master Ruler.', 'Feature Extraction Success');
            catch ME
                app.StatusLamp_3.Color = [1.00, 0.00, 0.00]; % Red
                uialert(app.UIFigure, ['An error occurred during feature extraction: ', ME.message], 'Error');
            end
        end
