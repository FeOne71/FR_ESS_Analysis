function [ocv_struct, static_struct] = process_one_channel_ocv(channel, rpt_cycles, ExperimentalDataPath, saveDir_OCV, ocv_conditions, ocv_steps, static_capacity_step)
% Process OCV and static capacity for one channel (for parfor/serial).
% - Reads each RPT CSV for the given channel and cycles
% - Extracts OCV (step 8,10 / substep 2) and Static (step 3 / substep 2)
% - Returns per-cycle structs and saves a per-channel OCV figure.

ocv_struct = struct();
static_struct = struct();

fig = figure('Name', sprintf('OCV - %s', channel), ...
    'Position', [100 100 1200 800], 'Visible', 'off');

for rpt_idx = 1:length(rpt_cycles)
    rpt_cycle = rpt_cycles{rpt_idx};
    field_name = sprintf('cyc%s', rpt_cycle(1:end-3));  % '0cyc' -> 'cyc0'

    % Ch14 RPT800cyc 스킵 (실험 중단)
    if strcmp(channel, 'Ch14') && strcmp(rpt_cycle, '800cyc')
        continue;
    end

    filename = sprintf('%s_RPT_%s.csv', channel, rpt_cycle);
    filepath = fullfile(ExperimentalDataPath, filename);
    if ~isfile(filepath)
        continue;
    end

    try
        T = readtable(filepath, 'VariableNamingRule', 'preserve');
    catch
        continue;
    end

    subplot(5,2,rpt_idx);
    hold on;

    charge_ocv = [];
    discharge_ocv = [];

    % Column indices (same as original script):
    % 2 = Step Index, 4 = Substep, 8 = Voltage(V), 9 = Capacity(Ah)
    for ocv_idx = 1:length(ocv_conditions)
        step_idx = ocv_steps(ocv_idx);
        ocv_data = (T{:,2} == step_idx) & (T{:,4} == 2);

        voltage_data = T{ocv_data, 8};

        % Reverse discharge OCV so voltage increases with SOC
        if strcmp(ocv_conditions{ocv_idx}, 'discharge')
            voltage_data = flipud(voltage_data);
        end

        % Moving average (window size = length / 200)
        window_size = max(1, round(length(voltage_data) / 200));
        voltage_smoothed = movmean(voltage_data, window_size);

        % SOC grid: 0~100%, 201 points
        soc_data = linspace(0, 100, 201);

        % Interpolate smoothed voltage onto SOC grid
        voltage_interp = interp1(linspace(0, 100, length(voltage_smoothed)), ...
                                 voltage_smoothed, soc_data, 'linear');

        if strcmp(ocv_conditions{ocv_idx}, 'charge')
            charge_ocv = voltage_interp;
            plot(soc_data, voltage_interp, 'b-', 'LineWidth', 2, 'DisplayName', 'Charge');
        else
            discharge_ocv = voltage_interp;
            plot(soc_data, voltage_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'Discharge');
        end
    end

    % Average of charge/discharge
    avg_ocv = (charge_ocv + discharge_ocv) / 2;
    plot(soc_data, avg_ocv, 'g-', 'LineWidth', 2, 'DisplayName', 'Average');

    % Save per-cycle OCV data into struct
    charge_mask    = (T{:,2} == ocv_steps(1)) & (T{:,4} == 2);
    discharge_mask = (T{:,2} == ocv_steps(2)) & (T{:,4} == 2);

    ocv_struct.(field_name).soc      = soc_data;
    ocv_struct.(field_name).avg_ocv  = avg_ocv;
    ocv_struct.(field_name).capacity = ...
        (T{find(charge_mask,1,'last'), 9} + T{find(discharge_mask,1,'last'), 9}) / 2;

    % Static Capacity (Step 3, Substep 2, discharge)
    static_dch_mask = (T{:,2} == static_capacity_step) & (T{:,4} == 2);
    if any(static_dch_mask)
        static_struct.(field_name).discharge = ...
            T{find(static_dch_mask,1,'last'), 9};
    else
        static_struct.(field_name).discharge = NaN;
    end

    xlabel('SOC [%]');
    ylabel('Voltage [V]');
    title(sprintf('%s - %s OCV', channel, rpt_cycle));
    legend('Location', 'best');
    grid on;
    xlim([0 100]);
end

% Save per-channel OCV figure
figName = fullfile(saveDir_OCV, sprintf('%s_OCV.fig', channel));
savefig(fig, figName);
close(fig);

end

