% MATLAB Script to Plot LQN Optimization Results Over Experiment Sequence
% AND Compare results across different Performance Weights
% Reads CSVs from './data/' folder.
% Exports individual subplot plots to './plots/' folder as PDF files.

clearvars;
clc;
close all; % Chiude eventuali figure aperte

% --- Configuration ---
data_folder = 'data'; % Subfolder for input CSVs
plots_folder = 'plots'; % Subfolder for output PDFs
filename_pattern_base = 'variable_load_results_PW*.csv';
export_pdf = true; % Flag per abilitare/disabilitare l'export PDF
target_fontsize = 20;

% --- Ensure plots folder exists ---
if export_pdf && ~isfolder(plots_folder)
    fprintf('Creating output folder: %s\n', plots_folder);
    mkdir(plots_folder);
end

% --- Construct full file pattern ---
filename_pattern = fullfile(data_folder, filename_pattern_base);

% --- PART 1: Plot Results for a Single CSV (Oldest Run based on modification time) ---
files = dir(filename_pattern);
if isempty(files)
    error('No CSV files found matching pattern: %s in folder %s', filename_pattern_base, data_folder);
end
[~,idx] = sort([files.datenum], 'ascend');
files = files(idx);
oldest_filename_short = files(1).name;
oldest_filename_full = fullfile(data_folder, oldest_filename_short);

fprintf('--- Part 1: Plotting results for oldest run: %s ---\n', oldest_filename_short);

% --- Read Data (Oldest Run) ---
if ~isfile(oldest_filename_full) error('File not found: %s.', oldest_filename_full); end
try
    T = readtable(oldest_filename_full);
    fprintf('Data read successfully. %d rows found.\n', height(T));
catch ME
    error('Failed to read CSV: %s\nError: %s', oldest_filename_full, ME.message);
end

% --- Extract Data Columns (Oldest Run) ---
T_valid = T; 
axesHandles1 = gobjects(5, 1);

if isempty(T_valid)
    warning('No results found in the oldest file to plot.');
else
    timeIndex = 1:height(T_valid);
users = T_valid.Users;
totalCost = T_valid.TotalCost;
responseTime = T_valid.ResponseTime;
selectedRate2 = T_valid.SelectedRate_Station2;
    hasStation3 = ismember('SelectedRate_Station3', T_valid.Properties.VariableNames);
    if hasStation3 selectedRate3 = T_valid.SelectedRate_Station3; end

    minResponseTime = 1 + (1/128);
    maxResponseTime = 1 + 1;

    % --- Create Figure (Invisible) ---
    % Aumenta le dimensioni della figura per provare ad accomodare il font grande
    fig1 = figure('Name', ['LQN Sequence - ', oldest_filename_short], 'NumberTitle', 'off', 'Position', [50, 50, 1200, 1500], 'Visible', 'off'); 
    sgtitle(['LQN Results Over Sequence (File: ', strrep(oldest_filename_short, '_', '\_'), ')'], 'FontSize', target_fontsize); 

    % 1. Users vs Experiment Sequence
    axesHandles1(1) = subplot(5, 1, 1);
    plot(timeIndex, users, 'b-', 'LineWidth', 1.5); % Linee leggermente più spesse
    title('User Load Over Time', 'FontSize', 25); 
    xlabel('Experiment Sequence Index', 'FontSize', target_fontsize); 
    ylabel('Users', 'FontSize', target_fontsize); grid on; xlim([1, length(timeIndex)]);
    set(gca, 'FontSize', target_fontsize); % Imposta font size per tick e label

    % 2. Total Cost vs Experiment Sequence
    axesHandles1(2) = subplot(5, 1, 2);
    plot(timeIndex, totalCost, 'r-', 'LineWidth', 1.5);
    title('Total Cost Over Time', 'FontSize', 25); 
    xlabel('Experiment Sequence Index', 'FontSize', target_fontsize); 
    ylabel('Total Cost', 'FontSize', target_fontsize); grid on; xlim([1, length(timeIndex)]);
    set(gca, 'FontSize', target_fontsize);

    % 3. Response Time vs Experiment Sequence
    axesHandles1(3) = subplot(5, 1, 3);
    hold on;
    plot(timeIndex, responseTime, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Actual Response Time');
    yline(minResponseTime, 'g--', 'LineWidth', 1, 'DisplayName', sprintf('Min RT (%.3f)', minResponseTime));
    yline(maxResponseTime, 'r--', 'LineWidth', 1, 'DisplayName', sprintf('Max RT (%.3f)', maxResponseTime));
    hold off;
    title('Average Response Time Over Time', 'FontSize', 25); 
    xlabel('Experiment Sequence Index', 'FontSize', target_fontsize); 
    ylabel('Response Time', 'FontSize', target_fontsize); grid on; xlim([1, length(timeIndex)]);
    legend('Location', 'best', 'FontSize', target_fontsize*0.8); % Legenda leggermente più piccola del titolo
    set(gca, 'FontSize', target_fontsize);

    % 4. Selected Service Rate Index for Station 2 vs Experiment Sequence
    axesHandles1(4) = subplot(5, 1, 4);
    stairs(timeIndex, selectedRate2, 'c-', 'LineWidth', 1.5);
    title('Selected Service Level - Task T_{1} Over Time', 'FontSize', 25); 
    xlabel('Experiment Sequence Index', 'FontSize', target_fontsize); 
    ylabel('Service Level', 'FontSize', target_fontsize); grid on; ylim([0 9]); yticks(1:8); xlim([1, length(timeIndex)]);
    set(gca, 'FontSize', target_fontsize);

    % 5. Selected Service Rate Index for Station 3 vs Experiment Sequence
    axesHandles1(5) = subplot(5, 1, 5);
    if hasStation3
        stairs(timeIndex, selectedRate3, 'k-', 'LineWidth', 1.5);
        title('Selected Service Level - Task T_{2} Over Time', 'FontSize', 25); 
        ylabel('Service Level', 'FontSize', target_fontsize); ylim([0 9]); yticks(1:8);
    else
        title('Station 3 Data Not Found', 'FontSize', target_fontsize); axis off;
    end
    xlabel('Experiment Sequence Index', 'FontSize', target_fontsize); grid on; xlim([1, length(timeIndex)]);
    set(gca, 'FontSize', target_fontsize);

    % --- Export Individual Subplots from Part 1 ---
    if export_pdf
        base_pdf_name1 = oldest_filename_short(1:end-4); 
        subplot_names1 = {'Users', 'Cost', 'ResponseTime', 'RateStn2', 'RateStn3'};

        for i = 1:length(axesHandles1)
            ax = axesHandles1(i);
            if isgraphics(ax) 
                % Riduci i margini al minimo possibile
                outerpos = ax.OuterPosition;
                ti = ax.TightInset; 
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax.Position = [left bottom ax_width ax_height];
                 
                pdf_filename_sub = fullfile(plots_folder, sprintf('%s_subplot_%s.pdf', base_pdf_name1, subplot_names1{i}));
                try
                    fprintf('Exporting subplot %d to %s...\n', i, pdf_filename_sub);
                    exportgraphics(ax, pdf_filename_sub, 'ContentType', 'vector', 'Resolution', 300);
                catch ME_export
                    warning('Failed to export subplot %d to PDF: %s', i, ME_export.message);
                end
            end
        end
         fprintf('Part 1 subplots exported (if any).\n');
    end
    close(fig1); 
end
fprintf('Part 1 analysis complete.\n');


% --- PART 2: Compare Results Across Different Performance Weights ---
fprintf('\n--- Part 2: Comparing results across different Performance Weights ---\n');

files = dir(filename_pattern); 
if isempty(files) error('No CSV files found matching pattern: %s in folder %s', filename_pattern_base, data_folder); end

num_files = length(files);
performanceWeights = zeros(num_files, 1);
avgCosts = zeros(num_files, 1);
avgResponseTimes = zeros(num_files, 1);
avgRate2 = zeros(num_files, 1);
avgRate3 = zeros(num_files, 1);
axesHandles2 = gobjects(3, 1); 

fprintf('Found %d files to analyze.\n', num_files);

for i = 1:num_files
    current_filename_short = files(i).name;
    current_filename_full = fullfile(data_folder, current_filename_short);
    fprintf('Processing file: %s\n', current_filename_short);

    % Estrai PerformanceWeight (come prima)
    try
        pw_str = regexp(current_filename_short, 'PW(\d+p\d+).csv', 'tokens');
        if isempty(pw_str) pw_str = regexp(current_filename_short, 'PW(\d+).csv', 'tokens'); end
        if isempty(pw_str) error('Pattern not found'); end
        weight_str = strrep(pw_str{1}{1}, 'p', '.');
        performanceWeights(i) = str2double(weight_str);
    catch
        warning('Could not parse PW from filename: %s. Skipping.', current_filename_short);
        performanceWeights(i) = NaN; continue;
    end

    % Leggi e processa i dati (come prima)
    try
        T_current = readtable(current_filename_full);
        valid_statuses = {'OPTIMAL', 'LOCALLY_SOLVED'};
        T_current_valid = T_current(ismember(T_current.Status, valid_statuses), :);
        if isempty(T_current_valid)
            warning('No successful runs in %s.', current_filename_short);
             avgCosts(i) = NaN; avgResponseTimes(i) = NaN; avgRate2(i) = NaN; avgRate3(i) = NaN;
             continue;
        end
        avgCosts(i) = mean(T_current_valid.TotalCost, 'omitnan');
        avgResponseTimes(i) = mean(T_current_valid.ResponseTime, 'omitnan');
        avgRate2(i) = mean(T_current_valid.SelectedRate_Station2, 'omitnan');
        if ismember('SelectedRate_Station3', T_current.Properties.VariableNames)
            avgRate3(i) = mean(T_current_valid.SelectedRate_Station3, 'omitnan');
        else avgRate3(i) = NaN; end
    catch ME
        warning('Failed read/process %s: %s. Skipping.', current_filename_short, ME.message);
        performanceWeights(i) = NaN; 
        avgCosts(i) = NaN; avgResponseTimes(i) = NaN; avgRate2(i) = NaN; avgRate3(i) = NaN;
        continue;
    end
end

valid_idx = ~isnan(performanceWeights) & ~isnan(avgCosts); 
performanceWeights = performanceWeights(valid_idx);
avgCosts = avgCosts(valid_idx);
avgResponseTimes = avgResponseTimes(valid_idx);
avgRate2 = avgRate2(valid_idx);
avgRate3 = avgRate3(valid_idx);

if isempty(performanceWeights)
    warning('No valid data found across different Performance Weights to plot.');
else
    [performanceWeights, sort_idx] = sort(performanceWeights);
    avgCosts = avgCosts(sort_idx);
    avgResponseTimes = avgResponseTimes(sort_idx);
    avgRate2 = avgRate2(sort_idx);
    avgRate3 = avgRate3(sort_idx);

    % --- Crea Figura di Confronto (invisibile) ---
    fig2 = figure('Name', 'Comparison Across PW', 'NumberTitle', 'off', 'Position', [100, 100, 1400, 1000], 'Visible', 'off'); % Figura più grande
    sgtitle('Average Metrics vs Performance Weight', 'FontSize', target_fontsize);

    % 1. Average Cost vs Performance Weight
    axesHandles2(1) = subplot(2, 2, 1);
    plot(performanceWeights, avgCosts, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 6); % Linee e marker più grandi
    title('Average Total Cost', 'FontSize', target_fontsize); 
    xlabel('Performance Weight (w_p)', 'FontSize', target_fontsize); 
    ylabel('Average Cost', 'FontSize', target_fontsize); grid on;
    xlim([min(performanceWeights)-0.05, max(performanceWeights)+0.05]);
    set(gca, 'FontSize', target_fontsize);

    % 2. Average Response Time vs Performance Weight
    axesHandles2(2) = subplot(2, 2, 2);
    plot(performanceWeights, avgResponseTimes, 'ms-', 'LineWidth', 1.5, 'MarkerSize', 6);
    title('Average Response Time', 'FontSize', target_fontsize); 
    xlabel('Performance Weight (w_p)', 'FontSize', target_fontsize); 
    ylabel('Avg Response Time', 'FontSize', target_fontsize); grid on;
    xlim([min(performanceWeights)-0.05, max(performanceWeights)+0.05]);
    set(gca, 'FontSize', target_fontsize);

    % 3. Average Selected Rate Index vs Performance Weight
    axesHandles2(3) = subplot(2, 2, 3); 
    plot(performanceWeights, avgRate2, 'cs-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Avg Rate Idx Stn 2');
    hold on;
    if any(~isnan(avgRate3))
        plot(performanceWeights, avgRate3, 'ks-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Avg Rate Idx Stn 3');
    end
    hold off;
    title('Average Selected Rate Index', 'FontSize', target_fontsize); 
    xlabel('Performance Weight (w_p)', 'FontSize', target_fontsize); 
    ylabel('Avg Rate Index (1-8)', 'FontSize', target_fontsize); grid on;
    ylim([0 9]); yticks(1:8);
    xlim([min(performanceWeights)-0.05, max(performanceWeights)+0.05]);
    legend('Location', 'best', 'FontSize', target_fontsize*0.8);
    set(gca, 'FontSize', target_fontsize);

    % --- Export Individual Subplots from Part 2 ---
     if export_pdf
        subplot_names2 = {'AvgCost', 'AvgResponseTime', 'AvgRates'};
        for i = 1:length(axesHandles2)
            ax = axesHandles2(i);
            if isgraphics(ax)
                 % Riduci margini
                 outerpos = ax.OuterPosition;
                 ti = ax.TightInset; 
                 left = outerpos(1) + ti(1);
                 bottom = outerpos(2) + ti(2);
                 ax_width = outerpos(3) - ti(1) - ti(3);
                 ax_height = outerpos(4) - ti(2) - ti(4);
                 ax.Position = [left bottom ax_width ax_height];

                 pdf_filename_sub = fullfile(plots_folder, sprintf('comparison_%s.pdf', subplot_names2{i}));
                 try
                    fprintf('Exporting subplot %d to %s...\n', i, pdf_filename_sub);
                    exportgraphics(ax, pdf_filename_sub, 'ContentType', 'vector', 'Resolution', 300);
                 catch ME_export
                    warning('Failed to export subplot %d to PDF: %s', i, ME_export.message);
                 end
            end
        end
         fprintf('Part 2 subplots exported (if any).\n');
     end
     close(fig2); 
end
fprintf('Part 2 analysis complete.\n');
