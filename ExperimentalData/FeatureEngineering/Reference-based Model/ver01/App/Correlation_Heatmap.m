%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPT Correlation Heatmap
% 피처 간 + 피처-라벨 피어슨 상관관계 히트맵
% 정사각형 전체 행렬 + p-value 유의성 표시
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Configuration
baseDir = 'D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData';
path_feature_mat = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'Dataset', 'Feature_Matrix_Final.mat');
saveDir = fullfile(baseDir, 'FeatureEngineering', 'Lab_RPT_Analysis', 'Pipeline_Visualizations');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

%% DataDir
fprintf('Loading Data...\n');
load(path_feature_mat, 'FeatureTable');

% 정규화된 피처 사용 (C-rate 스케일 차이 제거)
X = FeatureTable.X_Normalized;
Y = FeatureTable.Y_Labels;

% 전체 데이터 결합이 feat 14개 + label 3개 = 17개
all_data = [X, Y];

var_names = {'Chg dQ Seg1','Chg dQ Seg2','Chg dQ Seg3','Chg dQ Seg4','Chg dQ Seg5', ...
             'Dch dQ Seg1','Dch dQ Seg2','Dch dQ Seg3','Dch dQ Seg4','Dch dQ Seg5', ...
             'Chg dQdV PkH','Chg dQdV PkA','Dch dQdV PkH','Dch dQdV PkA', ...
             'SOH','LLI','LAM'};

n_vars = length(var_names);

%% correlation + p-value 
[R, P] = corr(all_data, 'type', 'Pearson');

% colormap: Blue → White → Red (논문 스타일)
n = 128;
blue = [0.15 0.25 0.55];
white = [1 1 1];
red = [0.7 0.15 0.15];

c1_r = linspace(blue(1), white(1), n)';
c1_g = linspace(blue(2), white(2), n)';
c1_b = linspace(blue(3), white(3), n)';
c2_r = linspace(white(1), red(1), n)';
c2_g = linspace(white(2), red(2), n)';
c2_b = linspace(white(3), red(3), n)';
cmap = [c1_r c1_g c1_b; c2_r c2_g c2_b];

%% P-value notation
% p<0.05: *, p<0.01: **, p<0.001: ***

%% === Figure 1: 충전 피처 상관 히트맵 (Charge) ===
chg_idx = [1:5, 11, 12, 15, 16, 17];
chg_labels = var_names(chg_idx);
n_chg = length(chg_idx);
R_chg = R(chg_idx, chg_idx);
P_chg = P(chg_idx, chg_idx);

fig1 = figure('Position', [50, 50, 800, 800], 'Color', 'w');
imagesc(R_chg);
colormap(cmap);
caxis([-1, 1]);
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Pearson correlation coefficients';
cb.Label.FontSize = 11;
cb.Label.FontWeight = 'bold';
axis square;

% 라벨: Y축(왼쪽) + X축(위쪽)
set(gca, 'XTick', 1:n_chg, 'XTickLabel', chg_labels, 'FontSize', 9, 'FontWeight', 'bold');
set(gca, 'YTick', 1:n_chg, 'YTickLabel', chg_labels, 'FontSize', 9, 'FontWeight', 'bold');
set(gca, 'XAxisLocation', 'top');
xtickangle(45);

% 대각선
hold on;
line([0.5, n_chg+0.5], [0.5, n_chg+0.5], 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

% 구분선: 피처/라벨 경계
line([7.5, 7.5], [0.5, n_chg+0.5], 'Color', 'k', 'LineWidth', 2);
line([0.5, n_chg+0.5], [7.5, 7.5], 'Color', 'k', 'LineWidth', 2);

% 숫자 + 유의성 표시
for i = 1:n_chg
    for j = 1:n_chg
        if i == j, continue; end
        val = R_chg(i,j);
        pval = P_chg(i,j);
        
        if pval < 0.001,     star = '***';
        elseif pval < 0.01,  star = '**';
        elseif pval < 0.05,  star = '*';
        else,                star = '';
        end
        
        if abs(val) > 0.6, clr = 'w'; else, clr = 'k'; end
        
        txt = sprintf('%.2f%s', val, star);
        text(j, i, txt, 'HorizontalAlignment', 'center', ...
            'FontSize', 7, 'FontWeight', 'bold', 'Color', clr);
    end
end
hold off;

title({'Charge Features Correlation Matrix', '(*p<0.05, **p<0.01, ***p<0.001)'}, ...
    'FontSize', 13, 'FontWeight', 'bold');

saveas(fig1, fullfile(saveDir, 'Heatmap_Charge_Matrix.fig'));
close(fig1);
fprintf('충전 피처 상관행렬 히트맵 완료\n');

%% === Figure 2: 방전 피처 상관 히트맵 (Discharge) ===
dch_idx = [6:10, 13, 14, 15, 16, 17];
dch_labels = var_names(dch_idx);
n_dch = length(dch_idx);
R_dch = R(dch_idx, dch_idx);
P_dch = P(dch_idx, dch_idx);

fig2 = figure('Position', [50, 50, 800, 800], 'Color', 'w');
imagesc(R_dch);
colormap(cmap);
caxis([-1, 1]);
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Pearson correlation coefficients';
cb.Label.FontSize = 11;
cb.Label.FontWeight = 'bold';
axis square;

set(gca, 'XTick', 1:n_dch, 'XTickLabel', dch_labels, 'FontSize', 9, 'FontWeight', 'bold');
set(gca, 'YTick', 1:n_dch, 'YTickLabel', dch_labels, 'FontSize', 9, 'FontWeight', 'bold');
set(gca, 'XAxisLocation', 'top');
xtickangle(45);

hold on;
line([0.5, n_dch+0.5], [0.5, n_dch+0.5], 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
line([7.5, 7.5], [0.5, n_dch+0.5], 'Color', 'k', 'LineWidth', 2);
line([0.5, n_dch+0.5], [7.5, 7.5], 'Color', 'k', 'LineWidth', 2);

for i = 1:n_dch
    for j = 1:n_dch
        if i == j, continue; end
        val = R_dch(i,j);
        pval = P_dch(i,j);
        
        if pval < 0.001,     star = '***';
        elseif pval < 0.01,  star = '**';
        elseif pval < 0.05,  star = '*';
        else,                star = '';
        end
        
        if abs(val) > 0.6, clr = 'w'; else, clr = 'k'; end
        
        txt = sprintf('%.2f%s', val, star);
        text(j, i, txt, 'HorizontalAlignment', 'center', ...
            'FontSize', 7, 'FontWeight', 'bold', 'Color', clr);
    end
end
hold off;

title({'Discharge Features Correlation Matrix', '(*p<0.05, **p<0.01, ***p<0.001)'}, ...
    'FontSize', 13, 'FontWeight', 'bold');

saveas(fig2, fullfile(saveDir, 'Heatmap_Discharge_Matrix.fig'));
close(fig2);
fprintf('방전 피처 상관행렬 히트맵 완료\n');

%% === Figure 3: 전체 피처 상관 히트맵 ===
fig3 = figure('Position', [50, 30, 1000, 1000], 'Color', 'w');
imagesc(R);
colormap(cmap);
caxis([-1, 1]);
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Pearson correlation coefficients';
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';
axis square;

set(gca, 'XTick', 1:n_vars, 'XTickLabel', var_names, 'FontSize', 8, 'FontWeight', 'bold');
set(gca, 'YTick', 1:n_vars, 'YTickLabel', var_names, 'FontSize', 8, 'FontWeight', 'bold');
set(gca, 'XAxisLocation', 'top');
xtickangle(45);

hold on;
line([0.5, n_vars+0.5], [0.5, n_vars+0.5], 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
% 충전dQ / 방전dQ / dQdV / 라벨 구분선
line([5.5, 5.5], [0.5, n_vars+0.5], 'Color', 'k', 'LineWidth', 2);
line([0.5, n_vars+0.5], [5.5, 5.5], 'Color', 'k', 'LineWidth', 2);
line([10.5, 10.5], [0.5, n_vars+0.5], 'Color', 'k', 'LineWidth', 2);
line([0.5, n_vars+0.5], [10.5, 10.5], 'Color', 'k', 'LineWidth', 2);
line([14.5, 14.5], [0.5, n_vars+0.5], 'Color', 'k', 'LineWidth', 2);
line([0.5, n_vars+0.5], [14.5, 14.5], 'Color', 'k', 'LineWidth', 2);

for i = 1:n_vars
    for j = 1:n_vars
        if i == j, continue; end
        val = R(i,j);
        pval = P(i,j);
        
        if pval < 0.001,     star = '***';
        elseif pval < 0.01,  star = '**';
        elseif pval < 0.05,  star = '*';
        else,                star = '';
        end
        
        if abs(val) > 0.6, clr = 'w'; else, clr = 'k'; end
        
        txt = sprintf('%.2f%s', val, star);
        text(j, i, txt, 'HorizontalAlignment', 'center', ...
            'FontSize', 5.5, 'FontWeight', 'bold', 'Color', clr);
    end
end
hold off;

title({'All Features & Labels Correlation Matrix', '(*p<0.05, **p<0.01, ***p<0.001)'}, ...
    'FontSize', 14, 'FontWeight', 'bold');

saveas(fig3, fullfile(saveDir, 'Heatmap_All_Matrix.fig'));
close(fig3);
fprintf('전체 피처 상관행렬 히트맵 완료\n');

fprintf('\n=== 모든 히트맵 생성 완료 ===\n');
fprintf('저장 위치: %s\n', saveDir);

%% === 엑셀 복사용 탭 구분 출력 ===
fprintf('\n\n========== 엑셀 복사용 (탭 구분) ==========\n');
fprintf('아래 내용을 복사하여 엑셀에 붙여넣으세요.\n');

var_names_raw = {'Chg_dQ_Seg1','Chg_dQ_Seg2','Chg_dQ_Seg3','Chg_dQ_Seg4','Chg_dQ_Seg5', ...
                 'Dch_dQ_Seg1','Dch_dQ_Seg2','Dch_dQ_Seg3','Dch_dQ_Seg4','Dch_dQ_Seg5', ...
                 'Chg_dQdV_PkH','Chg_dQdV_PkA','Dch_dQdV_PkH','Dch_dQdV_PkA', ...
                 'SOH','LLI','LAM'};

chg_names_raw = var_names_raw(chg_idx);
dch_names_raw = var_names_raw(dch_idx);

% =========================================================
% 1. 충전 (Charge) 상관계수 + 유의성
% =========================================================
fprintf('\n\n===== [Charge] Pearson r + significance =====\n');
fprintf('\t');
for j = 1:n_chg, fprintf('%s\t', chg_names_raw{j}); end
fprintf('\n');
for i = 1:n_chg
    fprintf('%s\t', chg_names_raw{i});
    for j = 1:n_chg
        if i == j
            fprintf('1.000\t');
        else
            pval = P_chg(i,j);
            if pval < 0.001,     star = '***';
            elseif pval < 0.01,  star = '**';
            elseif pval < 0.05,  star = '*';
            else,                star = '';
            end
            fprintf('%.3f%s\t', R_chg(i,j), star);
        end
    end
    fprintf('\n');
end

% 충전 p-value
fprintf('\n===== [Charge] p-value =====\n');
fprintf('\t');
for j = 1:n_chg, fprintf('%s\t', chg_names_raw{j}); end
fprintf('\n');
for i = 1:n_chg
    fprintf('%s\t', chg_names_raw{i});
    for j = 1:n_chg
        if i == j, fprintf('-\t');
        else, fprintf('%.2e\t', P_chg(i,j));
        end
    end
    fprintf('\n');
end

% =========================================================
% 2. 방전 (Discharge) 상관계수 + 유의성
% =========================================================
fprintf('\n\n===== [Discharge] Pearson r + significance =====\n');
fprintf('\t');
for j = 1:n_dch, fprintf('%s\t', dch_names_raw{j}); end
fprintf('\n');
for i = 1:n_dch
    fprintf('%s\t', dch_names_raw{i});
    for j = 1:n_dch
        if i == j
            fprintf('1.000\t');
        else
            pval = P_dch(i,j);
            if pval < 0.001,     star = '***';
            elseif pval < 0.01,  star = '**';
            elseif pval < 0.05,  star = '*';
            else,                star = '';
            end
            fprintf('%.3f%s\t', R_dch(i,j), star);
        end
    end
    fprintf('\n');
end

% 방전 p-value
fprintf('\n===== [Discharge] p-value =====\n');
fprintf('\t');
for j = 1:n_dch, fprintf('%s\t', dch_names_raw{j}); end
fprintf('\n');
for i = 1:n_dch
    fprintf('%s\t', dch_names_raw{i});
    for j = 1:n_dch
        if i == j, fprintf('-\t');
        else, fprintf('%.2e\t', P_dch(i,j));
        end
    end
    fprintf('\n');
end

% =========================================================
% 3. 전체 상관계수 + 유의성
% =========================================================
fprintf('\n\n===== [All] Pearson r + significance =====\n');
fprintf('\t');
for j = 1:n_vars, fprintf('%s\t', var_names_raw{j}); end
fprintf('\n');
for i = 1:n_vars
    fprintf('%s\t', var_names_raw{i});
    for j = 1:n_vars
        if i == j
            fprintf('1.000\t');
        else
            pval = P(i,j);
            if pval < 0.001,     star = '***';
            elseif pval < 0.01,  star = '**';
            elseif pval < 0.05,  star = '*';
            else,                star = '';
            end
            fprintf('%.3f%s\t', R(i,j), star);
        end
    end
    fprintf('\n');
end
