% Check saved events structure
load('EventsResults/all_events_raw_cell_level.mat', 'all_events');

fprintf('=== Checking saved events structure ===\n');
fprintf('Fields in all_events: %s\n', strjoin(fieldnames(all_events), ', '));

if isfield(all_events, 'Rack01')
    fprintf('\nRack01 exists\n');
    rack01 = all_events.Rack01;
    rack01_fields = fieldnames(rack01);
    fprintf('Fields in Rack01: %s\n', strjoin(rack01_fields, ', '));
    
    if ~isempty(rack01_fields)
        year = rack01_fields{1};
        fprintf('\nFirst year: %s\n', year);
        y_data = rack01.(year);
        y_fields = fieldnames(y_data);
        fprintf('Fields in %s: %s\n', year, strjoin(y_fields, ', '));
        
        if ~isempty(y_fields)
            month = y_fields{1};
            fprintf('\nFirst month: %s\n', month);
            m_data = y_data.(month);
            fprintf('Type of month data: %s\n', class(m_data));
            
            if iscell(m_data)
                fprintf('Number of events: %d\n', length(m_data));
                if length(m_data) > 0
                    evt = m_data{1};
                    fprintf('\nFields in first event: %s\n', strjoin(fieldnames(evt), ', '));
                    fprintf('Event type: %s\n', evt.type);
                    fprintf('Event duration: %.1f sec\n', seconds(evt.time_seq_datetime(end) - evt.time_seq_datetime(1)));
                else
                    fprintf('No events in this month\n');
                end
            else
                fprintf('Month data is not a cell array\n');
            end
        else
            fprintf('No months in this year\n');
        end
    else
        fprintf('No years in Rack01\n');
    end
else
    fprintf('Rack01 does not exist\n');
end

% Count total events
total_events = 0;
if isfield(all_events, 'Rack01')
    rack01 = all_events.Rack01;
    years = fieldnames(rack01);
    for y = 1:length(years)
        year_data = rack01.(years{y});
        months = fieldnames(year_data);
        for m = 1:length(months)
            month_events = year_data.(months{m});
            if iscell(month_events)
                total_events = total_events + length(month_events);
            end
        end
    end
end
fprintf('\n=== Total events: %d ===\n', total_events);
