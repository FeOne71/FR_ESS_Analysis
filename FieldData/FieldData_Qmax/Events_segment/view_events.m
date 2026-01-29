% View events in the saved structure
clear; clc;

load('EventsResults/all_events_raw_cell_level.mat', 'all_events');

fprintf('=== Viewing Events Structure ===\n\n');

% Navigate to first available event
if isfield(all_events, 'Rack01')
    rack01 = all_events.Rack01;
    years = fieldnames(rack01);
    
    fprintf('Available years: %s\n\n', strjoin(years, ', '));
    
    % Find first year with events
    for y_idx = 1:length(years)
        year = years{y_idx};
        year_data = rack01.(year);
        months = fieldnames(year_data);
        
        if ~isempty(months)
            fprintf('Year %s:\n', year);
            for m_idx = 1:length(months)
                month = months{m_idx};
                month_events = year_data.(month);
                
                if iscell(month_events) && length(month_events) > 0
                    fprintf('  Month %s: %d events\n', month, length(month_events));
                    
                    % Show first event details
                    evt = month_events{1};
                    fprintf('\n  First event in %s.%s:\n', year, month);
                    fprintf('    Type: %s\n', evt.type);
                    fprintf('    File: %s\n', evt.file_info.filename);
                    fprintf('    Duration: %.1f sec\n', seconds(evt.time_seq_datetime(end) - evt.time_seq_datetime(1)));
                    fprintf('    Fields: %s\n', strjoin(fieldnames(evt), ', '));
                    
                    % Show data sizes
                    fprintf('    Data sizes:\n');
                    fprintf('      I_raw: %d points\n', length(evt.current_seq_cell_A_raw));
                    fprintf('      I_smooth: %d points\n', length(evt.current_seq_cell_A));
                    fprintf('      V_raw: %d points\n', length(evt.voltage_seq_cell_V_raw));
                    fprintf('      V_smooth: %d points\n', length(evt.voltage_seq_cell_V));
                    fprintf('      SOC: %d points\n', length(evt.soc_seq_pct));
                    fprintf('      Time: %d points\n', length(evt.time_seq_datetime));
                    
                    % Show some data values
                    fprintf('\n    Sample data (first 5 points):\n');
                    fprintf('      Time: %s ... %s\n', ...
                        datestr(evt.time_seq_datetime(1), 'HH:MM:SS'), ...
                        datestr(evt.time_seq_datetime(min(5, length(evt.time_seq_datetime))), 'HH:MM:SS'));
                    fprintf('      I_raw: %.2f ... %.2f A\n', ...
                        evt.current_seq_cell_A_raw(1), ...
                        evt.current_seq_cell_A_raw(min(5, length(evt.current_seq_cell_A_raw))));
                    fprintf('      V_raw: %.3f ... %.3f V\n', ...
                        evt.voltage_seq_cell_V_raw(1), ...
                        evt.voltage_seq_cell_V_raw(min(5, length(evt.voltage_seq_cell_V_raw))));
                    fprintf('      SOC: %.1f ... %.1f %%\n', ...
                        evt.soc_seq_pct(1), ...
                        evt.soc_seq_pct(min(5, length(evt.soc_seq_pct))));
                    
                    fprintf('\n  To access this event in workspace:\n');
                    fprintf('    evt = all_events.Rack01.%s.%s{1};\n', year, month);
                    fprintf('    plot(evt.time_seq_datetime, evt.current_seq_cell_A);\n');
                    fprintf('\n');
                    
                    % Only show first month with events
                    break;
                end
            end
            break;
        end
    end
    
    % Count total events
    total_events = 0;
    for y_idx = 1:length(years)
        year_data = rack01.(years{y_idx});
        months = fieldnames(year_data);
        for m_idx = 1:length(months)
            month_events = year_data.(months{m_idx});
            if iscell(month_events)
                total_events = total_events + length(month_events);
            end
        end
    end
    fprintf('=== Total events: %d ===\n', total_events);
else
    fprintf('Rack01 not found\n');
end
