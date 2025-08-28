%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESS Charging Event Analysis Report Generator
% Automatic report generation for DCIR analysis and current clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DCIR_Report_Generator()
    
    % Load the analysis results
    load('G:\공유 드라이브\Battery Software Lab\0_Group Meeting\개인별_미팅자료\정철원\Materials\ESS_Data_Preprocessing\FieldData_DCIR_Charge\AutoResults\all_chg_events_current_clustering_auto_2021.mat');
    
    % Create report object
    import mlreportgen.report.*
    import mlreportgen.dom.*
    
    % Create report
    rpt = Report('ESS_Charging_Analysis_Report', 'pdf');
    
    % Add title page
    titlePage = TitlePage();
    titlePage.Title = 'ESS Charging Event Analysis Report';
    titlePage.Subtitle = 'DCIR Analysis and Current-based Clustering';
    titlePage.Author = 'Battery Software Lab';
    add(rpt, titlePage);
    
    % Add table of contents
    toc = TableOfContents();
    add(rpt, toc);
    
    % Executive Summary
    add(rpt, Heading1('Executive Summary'));
    summary = Paragraph();
    summary.append('This report presents the analysis of ESS charging events from 2021 field data. ');
    summary.append('The analysis includes automatic current-based clustering using K-means algorithm ');
    summary.append('and DCIR (Direct Current Internal Resistance) calculation for each cluster. ');
    summary.append('The results provide insights into charging patterns and resistance characteristics.');
    add(rpt, summary);
    
    % Analysis Parameters
    add(rpt, Heading1('Analysis Parameters'));
    params = Table();
    params.Style = {Border('solid'), Width('100%')};
    
    % Parameter table
    paramData = {
        'Parameter', 'Value', 'Description';
        'Battery Capacity', '1024 Ah', 'Nominal capacity';
        'Min Charge Duration', '30 s', 'Minimum charging time';
        'Max Power Std', '10 kW', 'Maximum power stability';
        'Max Current Std', '10.24 A', 'Maximum current stability';
        'Min Cluster Size', '3', 'Minimum events per cluster';
        'Max Clusters', '6', 'Maximum number of clusters';
        'Analysis Year', '2021', 'Data period'
    };
    
    for i = 1:size(paramData, 1)
        row = TableRow();
        for j = 1:size(paramData, 2)
            cell = TableEntry(paramData{i, j});
            row.append(cell);
        end
        params.append(row);
    end
    add(rpt, params);
    
    % Clustering Results
    add(rpt, Heading1('Current-based Clustering Results'));
    
    % Get cluster information
    clusterLabels = fieldnames(eventStruct);
    numClusters = length(clusterLabels);
    
    clusterSummary = Paragraph();
    clusterSummary.append(sprintf('A total of %d clusters were identified using K-means clustering. ', numClusters));
    clusterSummary.append('Each cluster represents a distinct charging current range. ');
    clusterSummary.append('The clustering was performed based on average charging current values.');
    add(rpt, clusterSummary);
    
    % Cluster details table
    clusterTable = Table();
    clusterTable.Style = {Border('solid'), Width('100%')};
    
    % Header
    headerRow = TableRow();
    headerRow.append(TableEntry('Cluster'));
    headerRow.append(TableEntry('C-rate'));
    headerRow.append(TableEntry('Events'));
    headerRow.append(TableEntry('Avg Current (A)'));
    headerRow.append(TableEntry('Description'));
    clusterTable.append(headerRow);
    
    % Data rows
    for i = 1:numClusters
        label = clusterLabels{i};
        events = fieldnames(eventStruct.(label));
        numEvents = length(events);
        
        % Calculate average current for this cluster
        totalCurrent = 0;
        validEvents = 0;
        for j = 1:numEvents
            if isfield(eventStruct.(label).(events{j}), 'avg_current')
                totalCurrent = totalCurrent + eventStruct.(label).(events{j}).avg_current;
                validEvents = validEvents + 1;
            end
        end
        avgCurrent = totalCurrent / validEvents;
        
        % Extract C-rate from label
        currentMatch = regexp(label, 'cluster_(\d+)A', 'tokens');
        if ~isempty(currentMatch)
            currentVal = str2double(currentMatch{1}{1});
            cRate = currentVal / 1024;
            cRateStr = sprintf('%.2fC', cRate);
        else
            cRateStr = 'N/A';
            cRate = NaN;  % Set default value
        end
        
        % Determine description based on C-rate
        if isnan(cRate) || cRate < 0.1
            desc = 'Low rate charging';
        elseif cRate < 0.2
            desc = 'Medium rate charging';
        else
            desc = 'High rate charging';
        end
        
        row = TableRow();
        row.append(TableEntry(label));
        row.append(TableEntry(cRateStr));
        row.append(TableEntry(sprintf('%d', numEvents)));
        row.append(TableEntry(sprintf('%.1f', avgCurrent)));
        row.append(TableEntry(desc));
        clusterTable.append(row);
    end
    add(rpt, clusterTable);
    
    % DCIR Analysis
    add(rpt, Heading1('DCIR Analysis Results'));
    
    dcirIntro = Paragraph();
    dcirIntro.append('DCIR (Direct Current Internal Resistance) analysis was performed for each cluster. ');
    dcirIntro.append('The analysis includes time-based DCIR calculations at 1s, 3s, 5s, 10s, 30s, and 50s intervals, ');
    dcirIntro.append('as well as DCIR difference calculations to assess resistance evolution during charging.');
    add(rpt, dcirIntro);
    
    % DCIR summary for each cluster
    for i = 1:numClusters
        label = clusterLabels{i};
        events = fieldnames(eventStruct.(label));
        numEvents = length(events);
        
        add(rpt, Heading2(sprintf('Cluster: %s (%d events)', label, numEvents)));
        
        % Calculate DCIR statistics
        dcirFields = {'DCIR_1s', 'DCIR_3s', 'DCIR_5s', 'DCIR_10s', 'DCIR_30s', 'DCIR_50s'};
        dcirLabels = {'1s', '3s', '5s', '10s', '30s', '50s'};
        
        dcirTable = Table();
        dcirTable.Style = {Border('solid'), Width('100%')};
        
        % Header
        dcirHeader = TableRow();
        dcirHeader.append(TableEntry('Time Interval'));
        dcirHeader.append(TableEntry('Mean DCIR (mΩ)'));
        dcirHeader.append(TableEntry('Std DCIR (mΩ)'));
        dcirHeader.append(TableEntry('Min DCIR (mΩ)'));
        dcirHeader.append(TableEntry('Max DCIR (mΩ)'));
        dcirTable.append(dcirHeader);
        
        % Data rows
        for j = 1:length(dcirFields)
            field = dcirFields{j};
            values = [];
            
            for k = 1:numEvents
                if isfield(eventStruct.(label).(events{k}), field) && ...
                   isfield(eventStruct.(label).(events{k}).(field), 'val')
                    val = eventStruct.(label).(events{k}).(field).val;
                    if ~isnan(val)
                        values = [values, val];
                    end
                end
            end
            
            if ~isempty(values)
                meanVal = mean(values);
                stdVal = std(values);
                minVal = min(values);
                maxVal = max(values);
                
                row = TableRow();
                row.append(TableEntry(dcirLabels{j}));
                row.append(TableEntry(sprintf('%.2f', meanVal)));
                row.append(TableEntry(sprintf('%.2f', stdVal)));
                row.append(TableEntry(sprintf('%.2f', minVal)));
                row.append(TableEntry(sprintf('%.2f', maxVal)));
                dcirTable.append(row);
            end
        end
        add(rpt, dcirTable);
    end
    
    % Key Findings
    add(rpt, Heading1('Key Findings'));
    
    findings = Paragraph();
    findings.append('1. **Current Distribution**: The analysis reveals distinct charging patterns with ');
    findings.append(sprintf('%d different current ranges identified.', numClusters));
    findings.append('\n\n');
    findings.append('2. **DCIR Characteristics**: Each cluster shows unique DCIR evolution patterns ');
    findings.append('during charging, indicating different resistance behaviors at various current levels.');
    findings.append('\n\n');
    findings.append('3. **Charging Efficiency**: Lower current clusters typically show more stable ');
    findings.append('DCIR values, while higher current clusters may exhibit more variation.');
    add(rpt, findings);
    
    % Conclusions
    add(rpt, Heading1('Conclusions'));
    
    conclusions = Paragraph();
    conclusions.append('The automatic clustering approach successfully identified distinct charging patterns ');
    conclusions.append('in the ESS field data. The DCIR analysis provides valuable insights into ');
    conclusions.append('resistance characteristics at different charging rates. This information ');
    conclusions.append('can be used for optimizing charging strategies and monitoring battery health.');
    add(rpt, conclusions);
    
    % Generate the report
    close(rpt);
    
    fprintf('Report generated successfully: ESS_Charging_Analysis_Report.pdf\n');
    
end 