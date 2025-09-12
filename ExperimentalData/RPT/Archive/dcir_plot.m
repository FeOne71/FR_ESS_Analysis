clear; clc; close all

%% 폴더/파일
data_folder = "G:\공유 드라이브\BSL_Data4\HNE_reduced_RPT_0812_processed_plot";
% data_folder = "HNE_10degC_1C0.33C_RuducedRPT_2_postprocessing_HPPC";
% save_path   = data_folder;
save_path ="D:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\RPT\Postprocessing";
files = dir(fullfile(data_folder, '*.mat'));

%% 파라미터
t_targets   = [1];        
dI_min_A    = 0.1;       
markerList  = {'o','s','^','d','v','>','<','p','h','x','+'};

%% 누적 버퍼
soc_all = []; R_all = []; file_all = {};

%% 피겨
fig = figure('Position',[100 100 1600 900]); hold on; grid on
xlabel('SOC (%)'); ylabel('R (m\Omega)'); title('SOC vs. R(1s)');

for fi = 1:numel(files)
    fname = files(fi).name;

   
    S = load(fullfile(data_folder, fname));
    if ~isfield(S,'n1C_pulse'), continue; end
    T = S.n1C_pulse;              

    soc_vals = []; R_vals = [];

    for k = 1:height(T)

        if any(ismissing(T{k,{'V','I','t'}})), continue; end
        V = T.V{k,1};            
        I = T.I{k,1};           
        t = T.t{k,1};            

        if isempty(V) || isempty(I) || isempty(t), continue; end


        t_rel = t - t(1);

       
        I_step = I(1);              
        if ~isfinite(I_step) || abs(I_step) < dI_min_A, continue; end

       
        R_step = NaN(numel(t_targets),1);
        for ti = 1:numel(t_targets)
            idx = find(t_rel >= t_targets(ti), 1, 'first');
            if ~isempty(idx) && isfinite(V(idx)) && isfinite(V(1))
                dV = V(idx) - V(1);
                dI = I_step;                    
                if abs(dI) >= dI_min_A
                    R_step(ti) = abs(dV/dI);     
                end
            end
        end


        if any(isnan(R_step)) && all(ismember({'V_final','V_after'}, T.Properties.VariableNames))
            V_final = T.V_final(k); 
            V_after = T.V_after(k);  
            if isfinite(V_final) && isfinite(V_after) && abs(I_step) >= dI_min_A
                R_step = abs((V_after - V_final) / I_step);  % [Ohm]
            end
        else
            R_step = R_step(1);     
        end

        if ~isnan(R_step)
            soc_vals(end+1,1) = T.SOC0(k);   
            R_vals(end+1,1)   = R_step;

     
            soc_all(end+1,1)  = T.SOC0(k);
            R_all(end+1,1)    = R_step;
            file_all{end+1,1} = fname;
        end
    end

    if ~isempty(soc_vals)
        [soc_sorted, idx] = sort(soc_vals, 'descend');
        R_sorted = R_vals(idx);

        mk = markerList{mod(fi-1,numel(markerList))+1};
        plot(soc_sorted, 1e3*R_sorted, ['-' mk], ...
            'LineWidth',1.5,'MarkerSize',4,'DisplayName',fname);
    end
end

legend('Location','northeast','Interpreter','none');
grid on; box on; hold off

