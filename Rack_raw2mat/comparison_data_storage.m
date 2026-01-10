%% 데이터 저장 방식 비교 예제
% 이 파일은 참고용입니다

%% 현재 방식: 구조체 (Struct)
% 저장:
Raw_struct = struct();
Raw_struct.Rack01 = struct();
Raw_struct.Rack01.Time = T.Time;
Raw_struct.Rack01.DCCurrent_A = T.DCCurrent_A;
Raw_struct.Rack01.AverageCV_V = T.AverageCV_V;
save('data_struct.mat', 'Raw_struct', '-v7.3');

% 로드 및 사용:
load('data_struct.mat');
current = Raw_struct.Rack01.DCCurrent_A;
time = Raw_struct.Rack01.Time;

%% 대안 1: 테이블(Table)로 저장
% 저장:
Raw_table = struct();
Raw_table.Rack01 = T;  % 테이블 그대로 저장
save('data_table.mat', 'Raw_table', '-v7.3');

% 로드 및 사용:
load('data_table.mat');
current = Raw_table.Rack01.DCCurrent_A;  % 더 직관적
time = Raw_table.Rack01.Time;
% 추가 장점: Raw_table.Rack01(:, {'Time', 'DCCurrent_A'}) 같은 슬라이싱 가능

%% 대안 2: 타임테이블(Timetable)로 저장
% 저장:
Raw_timetable = struct();
Raw_timetable.Rack01 = table2timetable(T, 'RowTimes', T.Time);
save('data_timetable.mat', 'Raw_timetable', '-v7.3');

% 로드 및 사용:
load('data_timetable.mat');
% 시간 기반 인덱싱/리샘플링이 매우 편함
% Raw_timetable.Rack01('2021-06-07 10:00:00', :)  % 특정 시간 데이터
% retime(Raw_timetable.Rack01, 'hourly', 'mean')  % 시간 리샘플링

%% 대안 3: HDF5 직접 사용 (더 효율적, 하지만 복잡)
% h5create, h5write 사용
% 파일 크기 더 작고, 부분 로드 가능
% 하지만 MATLAB 코드가 복잡해짐

%% 대안 4: Parquet 형식 (Python 등과 호환)
% MATLAB R2021b+에서 지원
% 하지만 MATLAB에서만 사용한다면 MAT 파일이 더 편함
