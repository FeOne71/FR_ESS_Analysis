%% 실험 전압 데이터 변동 폭 확인
clc; clear ; close all;

dataFolder = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\figures';

matFiles = {
    'Ch9_Drive_0cyc_real_load.mat';
    'Ch10_Drive_0cyc_real_load.mat';
    'Ch11_Drive_0cyc_real_load.mat';
    'Ch12_Drive_0cyc_real_load.mat';
    'Ch13_Drive_0cyc_real_load.mat';
    'Ch14_Drive_0cyc_real_load.mat';
    'Ch15_Drive_0cyc_real_load.mat';
    'Ch16_Drive_0cyc_real_load.mat';
};

ranges = {
    'Deep1',    1.5e4, 5.5e4;
    'Deep2',    0.8e5, 1.3e5;
    'Deep3',    1.65e5, 2.05e5;
    'Deep4',    2.35e5, 2.75e5;
    'Shallow1', 3.08e5, 3.74e5;
    'Shallow2', 3.52e5, 3.74e5;
    'Shallow3', 4.00e5, 4.30e5;
    'Shallow4', 4.50e5, 4.80e5;
};

for f = 1:length(matFiles)
    matPath = fullfile(dataFolder, matFiles{f});
    S = load(matPath);
    
    fieldName = fieldnames(S);
    data = S.(fieldName{1});  
    
    totalTime= data.totalTime;
    voltage  = data.voltage;
    current  = data.current;
    capacity = data.capacity;

    figure('Name', matFiles{f}, 'NumberTitle', 'off');
    
    for i = 1:size(ranges, 1)
        name = ranges{i, 1};
        idx_start = round(ranges{i, 2});
        idx_end   = round(ranges{i, 3});
        
        % if idx_end > length(totalTime)
        %     continue;
        % end

        % idx_end = min(idx_end, length(totalTime));

        t_seg = totalTime(idx_start:idx_end);
        V_seg = voltage(idx_start:idx_end);
        I_seg = current(idx_start:idx_end);
        Q_seg = capacity(idx_start:idx_end);

        [Vmin, idx_min] = min(V_seg);
        [Vmax, idx_max] = max(V_seg);
        Vdiff = Vmax - Vmin;

        Vinfo.(name).Vmin = Vmin;
        Vinfo.(name).Vmax = Vmax;
        Vinfo.(name).Vdiff = Vdiff;

        subplot(4,2,i);
        yyaxis left
        plot(t_seg, V_seg, 'b-', 'LineWidth', 1.2); hold on;
        plot(t_seg(idx_min), Vmin, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');  % Vmin
        plot(t_seg(idx_max), Vmax, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');  % Vmax
        ylabel('Voltage [V]');

        yyaxis right
        plot(t_seg, I_seg, 'r-', 'LineWidth', 1.2);
        ylabel('Current [A]');

        title(sprintf('%s (ΔV = %.3f V)', name, Vdiff), 'Interpreter', 'none');
        xlabel('Time');
        grid on;
    end

    disp(Vinfo);
end



selectedProfilePath = 'G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Experimental Data\Drive Cycle\selected_profiles.mat';

S = load(selectedProfilePath);  % 
profileNames = {'deep_1', 'deep_2', 'deep_3', 'deep_4', ...
                'shallow_1', 'shallow_2', 'shallow_3', 'shallow_4'};

Vinfo = struct();


figure('Name', 'Voltage/Current of Selected Profiles', 'NumberTitle', 'off');

for i = 1:length(profileNames)
    name = profileNames{i};
    profile = S.all_profiles.(name);  

    t = profile.time;
    V = profile.voltage;
    I = profile.current;

    [Vmin, idx_min] = min(V);
    [Vmax, idx_max] = max(V);
    Vdiff = (Vmax - Vmin)/(17*14);

    Vinfo.(name).Vmin = Vmin;
    Vinfo.(name).Vmax = Vmax;
    Vinfo.(name).Vdiff = Vdiff;

 
    subplot(4, 2, i);
    
    yyaxis left
    plot(t, V, 'b-', 'LineWidth', 1.3); hold on;
    plot(t(idx_min), Vmin, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');  % Vmin
    plot(t(idx_max), Vmax, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');  % Vmax
    ylabel('Voltage [V]');

    yyaxis right
    plot(t, I, 'r-', 'LineWidth', 1.3);
    ylabel('Current [A]');

    title(sprintf('%s (ΔV = %.3f V)', name, Vdiff), 'Interpreter', 'none');
    xlabel('Time');
    grid on;
end

% disp(Vinfo);