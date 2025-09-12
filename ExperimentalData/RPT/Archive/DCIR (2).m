clc;clear;close all;

%RPT data load
original_mat_path ="G:\공유 드라이브\Battery Software Group (2025)\Members\정미정\Parsing data\250423 RPT-Aging_Li_1.mat";
load(original_mat_path);

% 데이터 파일들이 저장된 경로를 변수로 지정
data_folder = "G:\공유 드라이브\Battery Software Group (2025)\Members\정미정\데이터";

% 데이터 파일 경로에서 디렉토리 추출
splitPath = split(data_folder, filesep);

% 최종 경로가 위치의 순서를 찾아서 반환
index = find(strcmp('데이터', splitPath), 1);

% splitPath의 {N번째} 행에 "DCIR graph plot"라는 문자열을 반환
splitPath{index} = 'DCIR graph plot';

% "Parsing data"가 저장될 경로를 변수로 지정 
save_path = strjoin(splitPath, filesep);

% 저장될 경로에 디렉토리가 없을 시, 생성
if ~exist(save_path, 'dir')
    mkdir(save_path);
end
%% 데이터 전처리
% 그래프 이름을 파일 이름으로 저장 
[~,base_filename,~]= fileparts(original_mat_path);

%% Charge DCIR
% 'C','R'이 6번 이상 반복되는 구간을 Charge DCIR 이라 지정
% 최소 반복 횟수 6번
min_repeat_CR = 6;

i =1;
data_len= length(data);
all_types = {data.type};

charge_start_index = NaN;
charge_end_index = NaN;

while i <= data_len -1
    repeat_CR = 0;
    start_idx = i;

% 'C','R'로 반복되는 구간을 탐색
    while i +1 <= data_len
        if strcmp(all_types{i},'C')&& strcmp (all_types{i+1},'R')
            repeat_CR = repeat_CR +1;
            i= i+2 ;
        else 
            break
        end
    end
    %'C','R'이 반복되는 구간의 시작 index 설정
    if repeat_CR >= min_repeat_CR
        charge_start_index = start_idx;
        charge_end_index = i - 1;  % 마지막 R까지 포함 
        break;  % 조건 만족 시 종료
    else
        i = start_idx + 1;  % 불충족 시 한 칸 뒤로 이동해서 재탐색
    end
end

% 구조체 추출
if ~isnan(charge_start_index)
    data_DCIR_chg = data(charge_start_index : charge_end_index);
end

%% Dicharge DCIR
% 'D','R'이 6번 이상 반복되는 구간을 Discharge DCIR 이라 지정
% 최소 반복 횟수 6번
min_repeat_DR = 6;

i =1;
data_len= length(data);
all_types = {data.type};

Discharge_start_index = NaN;
Discharge_end_index = NaN;

while i <= data_len -1
    repeat_DR = 0;
    start_idx = i;
% 'D','R'로 반복되는 구간을 탐색
    while i +1 <= data_len
        if strcmp(all_types{i},'D')&& strcmp (all_types{i+1},'R')
            repeat_DR = repeat_DR +1;
            i= i+2 ;
        else 
            break
        end
    end
        %'D','R'이 반복되는 구간의 시작 index 설정
    if repeat_DR >= min_repeat_DR
        Discharge_start_index = start_idx;
        Discharge_end_index = i - 1;  % 마지막 R까지 포함
        break;  % 조건 만족 시 종료
    else
        i = start_idx + 1;  % 불충족 시 한 칸 뒤로 이동해서 재탐색
    end
end

% 구조체 추출
if ~isnan(Discharge_start_index)
    data_DCIR_dis = data(Discharge_start_index : Discharge_end_index);
end

%% Q 계산하기
%charge 
% charge 상태에서 t, I 계산하기 위해 반출할 공간 생성
t_chg_all = [];
I_chg_all = [];

for j = 1:length(data_DCIR_chg)
    t_tmp = data_DCIR_chg(j).t(:);
    I_tmp = data_DCIR_chg(j).I(:)*1000; %[mAh]

    t_chg_all = [t_chg_all; t_tmp];
    I_chg_all = [I_chg_all; I_tmp];
end
% cumQ, Q 계산
cumQ_chg = cumtrapz(t_chg_all, I_chg_all) / 3600;
Q_chg = trapz(t_chg_all, I_chg_all) / 3600;


% discharge 상태
t_dis_all = [];
I_dis_all = [];

% Discharge 상태에서 t, I 계산하기 위해 반출할 공간 생성
for j = 1:length(data_DCIR_dis)
    t_dis_tmp = data_DCIR_dis(j).t(:);
    I_dis_tmp = abs(data_DCIR_dis(j).I(:))*1000; %[mAh]

    t_dis_all = [t_dis_all; t_dis_tmp];
    I_dis_all = [I_dis_all; I_dis_tmp];
end

% cumQ, Q 계산
cumQ_dis = cumtrapz(t_dis_all, I_dis_all) / 3600;
Q_dis = trapz(t_dis_all, I_dis_all) / 3600;


%SOC 계산 (용량이 큰 것을 기준으로)
    if Q_dis > Q_chg
        Qref = Q_dis;
    else
        Qref= Q_chg;
    end

% Charge SOC
data_DCIR_chg_soc = cumQ_chg/Qref;

% Discharge SOC
data_DCIR_dis_soc = cumQ_dis/Qref;

%% step 시간별로 SOC 잘라 넣기
%SOC 삽입
% charge 
% 계산한 SOC 를 step 별로 잘라 넣음
start_idx = 1;
for j = 1:length(data_DCIR_chg)
    n = length(data_DCIR_chg(j).steptime);
    idx_range = start_idx : start_idx + n - 1;

    data_DCIR_chg(j).soc = data_DCIR_chg_soc(idx_range);  % row로 맞춤
    start_idx = start_idx + n;
end

%discharge
% 계산한 SOC 를 step 별로 잘라 넣음
start_idx = 1;
   
for j = 1:length(data_DCIR_dis)
    n = length(data_DCIR_dis(j).steptime);
    idx_range = start_idx : start_idx + n - 1;

    data_DCIR_dis(j).soc = data_DCIR_dis_soc(idx_range);  % row로 맞춤
    start_idx = start_idx + n;
end

%% DCIR 계산하기(charge)

% DCIR 경과 시간 계산
chg_data_DCIR_selected = data_DCIR_chg([]);  
count = 0;

% DCIR test C 만 추출 DCIR_time이 30초이상 31초 미만인 것
for i = 1:length(data_DCIR_chg)
    if isfield(data_DCIR_chg(i), 't') && ~isempty(data_DCIR_chg(i).t)
        t_vec = data_DCIR_chg(i).t;
        t_elapsed = t_vec(end) - t_vec(1);  % 경과 시간 계산

        if t_elapsed >=30 && t_elapsed <=31
            count = count + 1;
            chg_data_DCIR_selected(count) = data_DCIR_chg(i);
        end
    end
end

% step time을 double형태로 만들기
for i = 1:length(chg_data_DCIR_selected)
    if isfield(chg_data_DCIR_selected(i), 'steptime')
        chg_data_DCIR_selected(i).steptime_double = seconds(chg_data_DCIR_selected(i).steptime);
    end
end
%% R 계산하기 (Charge) 
target_times = [1, 5, 25];
s= chg_data_DCIR_selected;
% I_1C =0.0038

% target time에 따라 추출하기
for i = 1:length(s)
    t_vec = s(i).steptime_double;  % 이미 double
    V_vec = s(i).V;
    soc_vec = s(i).soc;
    I_1C_chg = mean(s(i).I);

    idx_list = zeros(size(target_times));

    for k = 1:length(target_times)
        idx_list(k) = find(t_vec == target_times(k), 1, 'first');  % '정확히' 일치
    end

    % 인덱스로 추출
    chg_data_DCIR_selected(i).V_pulse = V_vec(idx_list);
    chg_data_DCIR_selected(i).SOC = soc_vec(1)*100;
    chg_data_DCIR_selected(i).V0 = V_vec(3);
    
end

% 계산하기
for i = 1:length(chg_data_DCIR_selected)
    V = chg_data_DCIR_selected(i).V_pulse;
    V0 = chg_data_DCIR_selected(i).V0;
      
    % ΔV
    dV1 = V(1) - V0 ;  % 1s - 0s
    dV2 = V(2) - V(1);  % 5s - 1s
    dV3 = V(3) - V(2);  % 25s - 5s
    
    % 저항 계산: ΔV / ΔR
    chg_data_DCIR_selected(i).dv= [dV1, dV2, dV3];
    chg_data_DCIR_selected(i).dR = [dV1/I_1C_chg,dV2/I_1C_chg,dV3/I_1C_chg];
    chg_data_DCIR_selected(i).R = [dV1/I_1C_chg, (dV1+dV2)/I_1C_chg,(dV1+dV2+dV3)/I_1C_chg];
end


%% DICR 계산하기(Discharge)
% DCIR 경과 시간 계산
dis_data_DCIR_selected = data_DCIR_dis([]);  
count = 0;

% DCIR test D만 추출 DCIR_time이 30초 이상 31초 이하인 것
for i = 1:length(data_DCIR_dis)
    if isfield(data_DCIR_dis(i), 't') && ~isempty(data_DCIR_dis(i).t)
        t_vec = data_DCIR_dis(i).t;
        t_elapsed = t_vec(end) - t_vec(1);  % 경과 시간 계산

        if t_elapsed >=30 && t_elapsed <=31
            count = count + 1;
            dis_data_DCIR_selected(count) = data_DCIR_dis(i);
        end
    end
end

% step time을 double형태로 만들기
for i = 1:length(dis_data_DCIR_selected)
    if isfield(dis_data_DCIR_selected(i), 'steptime')
        dis_data_DCIR_selected(i).steptime_double = seconds(dis_data_DCIR_selected(i).steptime);
    end
end
%% R 계산하기 (Discharge) 
target_times = [1, 5, 25];
s= dis_data_DCIR_selected;

% target time에 따라 추출하기
for i = 1:length(s)
    t_vec = s(i).steptime_double;  % 이미 double
    V_vec = s(i).V;
    soc_vec = s(i).soc;
    I_1C_dis = mean(s(i).I);

    idx_list = zeros(size(target_times));

    for k = 1:length(target_times)
        idx_list(k) = find(t_vec == target_times(k), 1, 'first');  % target과 '정확히' 일치
    end

    % 인덱스로 추출
    dis_data_DCIR_selected(i).V_pulse = V_vec(idx_list);
    dis_data_DCIR_selected(i).SOC = soc_vec(1)*100;
    dis_data_DCIR_selected(i).V0 = V_vec(3);
end

% 계산하기
for i = 1:length(dis_data_DCIR_selected)
    V = dis_data_DCIR_selected(i).V_pulse;
    V0 = dis_data_DCIR_selected(i).V0;

    % ΔV
    dV1 = V(1) - V0;  % 1s - 0s
    dV2 = V(2) - V(1);  % 5s - 1s
    dV3 = V(3) - V(2);  % 30s - 5s

    % 저항 계산: ΔV / ΔI
    dis_data_DCIR_selected(i).dv= [dV1, dV2, dV3];
    dis_data_DCIR_selected(i).dR = [dV1/I_1C_dis,dV2/I_1C_dis, dV3/I_1C_dis];
    dis_data_DCIR_selected(i).R = [dV1/I_1C_dis, (dV1+dV2)/I_1C_dis,(dV1+dV2+dV3)/I_1C_dis];
end

%% figure 그리기
figure(1);
hold on;

n = length(chg_data_DCIR_selected);
SOC = zeros(n, 1);
dR = zeros(n, 3);

for i = 1:n
    SOC(i) = chg_data_DCIR_selected(i).SOC;  % 실제 SOC 값
    dR(i,:) = chg_data_DCIR_selected(i).dR;  % 1x3 vector
end

% 균일하게 배치된 가짜 x축 위치 (1, 2, 3, ...)
x_pos = 1:n;

% bar 그리기
bar(x_pos, dR, 'stacked');

% x축 눈금 설정 → SOC 값 표시
xticks(x_pos);
[~, file_name, ~] = fileparts(original_mat_path);  % 확장자 없이 파일 이름만 추출
title(file_name, 'Interpreter', 'none','FontSize',14)
xticklabels(arrayfun(@(x) sprintf('%.2f', x), SOC, 'UniformOutput', false));
xtickangle(45);  % x축 라벨 기울이기

set(gca)
xlabel('State of Charge (%)','FontSize',14);
ylabel('DCIR (Ω)','FontSize',14);

legend({'1s','5s','30s'}, 'Location', 'northeast','FontSize',14) 


%% figure 그리기
% Discharge
figure(2);
hold on;

n = length(dis_data_DCIR_selected);
SOC = zeros(n, 1);
dR = zeros(n, 3);

for i = 1:n
    SOC(i) = dis_data_DCIR_selected(i).SOC;  % 실제 SOC 값
    dR(i,:) = dis_data_DCIR_selected(i).dR;  % 1x3 vector
end

% 균일하게 배치된 가짜 x축 위치 (1, 2, 3, ...)
x_pos = 1:n;

% bar 그리기
bar(x_pos, dR, 'stacked');

% x축 눈금 설정 → SOC 값 표시
xticks(x_pos);[~, file_name, ~] = fileparts(original_mat_path);  % 확장자 없이 파일 이름만 추출
title(file_name, 'Interpreter', 'none','FontSize',14)
xticklabels(arrayfun(@(x) sprintf('%.2f', x), SOC, 'UniformOutput', false));
xtickangle(45);  % x축 라벨 기울이기

set(gca)
xlabel('Depth of Discharge (%)','FontSize',14);
ylabel('DCIR (Ω)','FontSize',14);

legend({'1s','5s','30s'}, 'Location', 'northwest','FontSize',14) 


%% 저장하기
for fig_num = 1:2
    figure(fig_num);  % 해당 figure로 포커스
    png_save_path = fullfile(save_path, sprintf('%s_fig%d.png', char(base_filename), fig_num));
    saveas(gcf, png_save_path);
end
mat_save_path = fullfile(save_path,[char(base_filename),'.mat']);
save(mat_save_path);
% clearvars -except data data_DCIR_chg data_DCIR_dis chg_data_DCIR_selected dis_data_DCIR_selected I_1C_chg I_1C_dis cumQ_chg cumQ_dis ;
