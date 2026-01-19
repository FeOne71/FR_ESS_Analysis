%% 배터리 데이터 구조 분석 스크립트
% 작성일: 2026-01-13
% 목적: 연도별/월별 폴더 구조의 mat 파일 분석

clear; clc; close all;

%% 경로 설정
baseDataPath = 'D:\JCW\Projects\KEPCO_ESS_Local\Rack_raw2mat\Old';
savePath = 'D:\JCW\Projects\KEPCO_ESS_Local\FieldData\FieldData_0113';

% 저장 경로가 없으면 생성
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

fprintf('=== 배터리 데이터 구조 분석 시작 ===\n\n');

%% 1단계: 폴더 구조 탐색
fprintf('[1단계] 폴더 구조 탐색 중...\n');

% 연도별 폴더 찾기
yearFolders = dir(baseDataPath);
yearFolders = yearFolders([yearFolders.isdir] & ~ismember({yearFolders.name}, {'.', '..'}));

fprintf('발견된 연도별 폴더 수: %d\n', length(yearFolders));

% 전체 구조 저장용
folderStructure = cell(length(yearFolders), 1);
allMatFiles = {};

for i = 1:length(yearFolders)
    yearName = yearFolders(i).name;
    yearPath = fullfile(baseDataPath, yearName);
    fprintf('  연도: %s\n', yearName);
    
    % 월별 폴더 찾기
    monthFolders = dir(yearPath);
    monthFolders = monthFolders([monthFolders.isdir] & ~ismember({monthFolders.name}, {'.', '..'}));
    
    yearData.name = yearName;
    yearData.months = cell(length(monthFolders), 1);
    
    for j = 1:length(monthFolders)
        monthName = monthFolders(j).name;
        monthPath = fullfile(yearPath, monthName);
        fprintf('    월: %s\n', monthName);
        
        % mat 파일 찾기
        matFiles = dir(fullfile(monthPath, '*.mat'));
        
        monthData.name = monthName;
        monthData.files = {matFiles.name};
        yearData.months{j} = monthData;
        
        % 전체 mat 파일 목록에 추가
        for k = 1:length(matFiles)
            allMatFiles{end+1} = fullfile(monthPath, matFiles(k).name);
        end
        
        fprintf('      mat 파일 수: %d\n', length(matFiles));
    end
    
    folderStructure{i} = yearData;
end

fprintf('\n총 mat 파일 수: %d\n\n', length(allMatFiles));

%% 2단계: 샘플 파일 분석
fprintf('[2단계] 샘플 파일 분석 중...\n');

if isempty(allMatFiles)
    error('mat 파일을 찾을 수 없습니다!');
end

% 첫 번째 파일을 샘플로 사용
sampleFile = allMatFiles{1};
fprintf('샘플 파일: %s\n\n', sampleFile);

% mat 파일 로드
try
    data = load(sampleFile);
    fprintf('파일 로드 성공!\n\n');
catch ME
    error('파일 로드 실패: %s', ME.message);
end

%% 3단계: 변수 정보 추출
fprintf('[3단계] 변수 정보 추출 중...\n');

varNames = fieldnames(data);
fprintf('발견된 변수 수: %d\n\n', length(varNames));

varInfo = struct();

for i = 1:length(varNames)
    varName = varNames{i};
    varData = data.(varName);
    
    varInfo(i).name = varName;
    varInfo(i).size = size(varData);
    varInfo(i).class = class(varData);
    varInfo(i).bytes = whos('varData').bytes;
    
    fprintf('변수 %d: %s\n', i, varName);
    fprintf('  - 크기: %s\n', mat2str(size(varData)));
    fprintf('  - 타입: %s\n', class(varData));
    fprintf('  - 메모리: %.2f KB\n', whos('varData').bytes / 1024);
    
    % 데이터 미리보기
    if isnumeric(varData)
        if numel(varData) <= 10
            fprintf('  - 데이터: %s\n', mat2str(varData));
        else
            fprintf('  - 데이터 범위: [%.4f ~ %.4f]\n', min(varData(:)), max(varData(:)));
            fprintf('  - 평균: %.4f\n', mean(varData(:), 'omitnan'));
        end
    elseif ischar(varData) || isstring(varData)
        if length(varData) < 100
            fprintf('  - 내용: %s\n', varData);
        else
            fprintf('  - 내용: %s... (잘림)\n', varData(1:100));
        end
    elseif iscell(varData)
        fprintf('  - 셀 배열 크기: %s\n', mat2str(size(varData)));
        if ~isempty(varData)
            fprintf('  - 첫 요소 타입: %s\n', class(varData{1}));
        end
    elseif isstruct(varData)
        subFields = fieldnames(varData);
        fprintf('  - 구조체 필드 수: %d\n', length(subFields));
        fprintf('  - 필드: %s\n', strjoin(subFields, ', '));
    end
    
    fprintf('\n');
end

%% 4단계: 결과 저장
fprintf('[4단계] 결과 저장 중...\n');

% 구조 정보 저장
save(fullfile(savePath, 'folder_structure.mat'), 'folderStructure', 'allMatFiles');
save(fullfile(savePath, 'variable_info.mat'), 'varInfo', 'sampleFile');

% 텍스트 리포트 생성
reportFile = fullfile(savePath, 'data_structure_report.txt');
fid2 = fopen(reportFile, 'w');

fprintf(fid2, '===== 배터리 데이터 구조 분석 리포트 =====\n');
fprintf(fid2, '분석 일시: %s\n\n', datestr(now));

fprintf(fid2, '[폴더 구조]\n');
fprintf(fid2, '기본 경로: %s\n', baseDataPath);
fprintf(fid2, '총 mat 파일 수: %d\n\n', length(allMatFiles));

fprintf(fid2, '[샘플 파일]\n');
fprintf(fid2, '%s\n\n', sampleFile);

fprintf(fid2, '[변수 정보]\n');
for i = 1:length(varInfo)
    fprintf(fid2, '%d. %s\n', i, varInfo(i).name);
    fprintf(fid2, '   크기: %s\n', mat2str(varInfo(i).size));
    fprintf(fid2, '   타입: %s\n', varInfo(i).class);
    fprintf(fid2, '   메모리: %.2f KB\n\n', varInfo(i).bytes / 1024);
end

fclose(fid2);

fprintf('결과 저장 완료!\n');
fprintf('  - folder_structure.mat\n');
fprintf('  - variable_info.mat\n');
fprintf('  - data_structure_report.txt\n\n');

fprintf('=== 분석 완료 ===\n');

% 작업 공간에 주요 변수 표시
fprintf('\n작업 공간에 로드된 주요 변수:\n');
fprintf('  - folderStructure: 폴더 구조 정보\n');
fprintf('  - allMatFiles: 모든 mat 파일 경로 목록\n');
fprintf('  - varInfo: 변수 상세 정보\n');
fprintf('  - data: 샘플 파일 데이터\n');
