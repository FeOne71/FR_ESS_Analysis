%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab_Data_Fusion.m (FINAL FIX: Channel/CyclePeriod String Type Issue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
warning off;

%% 경로 설정
featureDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Drive Cycle\features_anchored_soc';
labelDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing\Labels';
outputFolder = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\Model_Data';

if ~exist(outputFolder,'dir'); mkdir(outputFolder); end

%% 데이터 로드
load(fullfile(featureDataPath, 'Lab_DutyCycle_Features_X_Lab.mat'), 'featureTable_Lab');
load(fullfile(labelDataPath, 'Lab_Degradation_Y_Labels_Final.mat'), 'Y_Lab_Table_Base');

fprintf('X_Lab (%d rows)와 Y_Lab (%d rows) 로드 완료.\n', height(featureTable_Lab), height(Y_Lab_Table_Base));

%% 1. X_Lab 테이블의 채널 이름 정리 (최종 수정)
% 'ch9_Drive_0cyc' 형태를 'Ch9' 형태로 변경하고 char 배열로 통일합니다.
original_channels = featureTable_Lab.Channel;
cleaned_channels = cell(height(featureTable_Lab), 1);

for i = 1:length(original_channels)
    channel_part = original_channels{i};
    idx_start = strfind(channel_part, 'ch');
    
    if ~isempty(idx_start)
        idx_end = strfind(channel_part, '_');
        if ~isempty(idx_end)
            channel_name_raw = channel_part(idx_start(1) : idx_end(1)-1);
        else
            channel_name_raw = channel_part(idx_start(1) : end);
        end
        % 🛠️ 최종 FIX: char 배열로 명시적 변환
        cleaned_channels{i} = char([upper(channel_name_raw(1)), channel_name_raw(2:end)]);
    else
        cleaned_channels{i} = original_channels{i};
    end
end
featureTable_Lab.Channel = cleaned_channels;
fprintf('X_Lab 채널 이름 정리 완료: 예) ch9_Drive_0cyc -> Ch9\n');


%% 2. CyclePeriod 칼럼 생성 및 정리 (CycleType 문자열 정돈)

% X_Lab의 CycleType 칼럼 자체를 정돈합니다.
featureTable_Lab.CycleType = strtrim(featureTable_Lab.CycleType);

% CyclePeriod 생성
X_Lab_Periods = repmat({''}, height(featureTable_Lab), 1);

for i = 1:height(featureTable_Lab)
    
    cycleType = featureTable_Lab.CycleType{i};
    
    if strcmp(cycleType, '0cyc')
        X_Lab_Periods{i} = '0to200cyc';
    elseif strcmp(cycleType, '200cyc')
        X_Lab_Periods{i} = '200to400cyc';
    else 
        X_Lab_Periods{i} = 'NaN_Cycle_Period';
    end
end

featureTable_Lab.CyclePeriod = X_Lab_Periods;
fprintf('X_Lab CyclePeriod 매핑 완료: CycleType 정돈 후 매칭\n');


%% 3. Y_Lab 칼럼을 X_Lab에 매칭하여 추가 (Fusion)

featureTable_Lab.SOH_Loss_Rate_Label = NaN(height(featureTable_Lab), 1);
featureTable_Lab.CL_DCIR_Rate_Label = NaN(height(featureTable_Lab), 1);
featureTable_Lab.LLI_Loss_Rate_Label = NaN(height(featureTable_Lab), 1);
featureTable_Lab.LAM_Loss_Rate_Label = NaN(height(featureTable_Lab), 1);

for i = 1:height(Y_Lab_Table_Base)
    
    Y_channel = Y_Lab_Table_Base.Channel{i};
    Y_period = Y_Lab_Table_Base.CyclePeriod{i};
    
    % 🛠️ 매칭 로직: Channel과 CyclePeriod가 모두 일치하는 행을 찾음
    matchMask = strcmp(featureTable_Lab.Channel, Y_channel) & ...
                strcmp(featureTable_Lab.CyclePeriod, Y_period);
    
    if ~any(matchMask)
        fprintf('경고: Y_Lab 행 %d (%s, %s)에 매칭되는 X_Lab 데이터가 없습니다. (데이터 불일치)\n', i, Y_channel, Y_period);
        continue;
    end
    
    % 매칭되는 모든 행에 Y_Lab 값 복사 (192행씩 복사되어야 함)
    featureTable_Lab.SOH_Loss_Rate_Label(matchMask) = Y_Lab_Table_Base.SOH_Loss_Rate(i);
    featureTable_Lab.CL_DCIR_Rate_Label(matchMask) = Y_Lab_Table_Base.CL_DCIR_Rate(i);
    featureTable_Lab.LLI_Loss_Rate_Label(matchMask) = Y_Lab_Table_Base.LLI_Loss_Rate(i);
    featureTable_Lab.LAM_Loss_Rate_Label(matchMask) = Y_Lab_Table_Base.LAM_Loss_Rate(i);
end

% 최종 모델 학습 데이터셋 (X+Y)
modelDataSet_Lab = featureTable_Lab;

%% 4. 결과 저장
saveFileName = fullfile(outputFolder, 'Model_Training_DataSet_Lab.mat');
save(saveFileName, 'modelDataSet_Lab');
fprintf('\n=== 모델 학습 데이터셋 (X_Lab + Y_Lab) 통합 완료 ===\n');
fprintf('최종 데이터셋 크기: %d 행 x %d 열\n', size(modelDataSet_Lab, 1), size(modelDataSet_Lab, 2));
fprintf('결과 파일: %s\n', saveFileName);