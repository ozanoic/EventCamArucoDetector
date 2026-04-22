function viewAllResults(dataRoot)
%VIEWALLRESULTS Tabbed GUI to browse every dataset's detection results.
%
%  viewAllResults()              % defaults to 'Data'
%  viewAllResults('Data')
%
%  Scans dataRoot for <folder>/<folder>_results.mat (merged files produced
%  by mergeAllResults) and opens a single window with:
%    - Tab 1:    Summary comparing all datasets
%    - Tab 2..N: One detailed tab per dataset

if nargin < 1 || isempty(dataRoot)
    dataRoot = 'Data';
end

if ~isfolder(dataRoot)
    error('viewAllResults: folder not found: %s', dataRoot);
end

%% ---- Discover merged result files ----
d = dir(dataRoot);
subDirs = d([d.isdir] & ~ismember({d.name}, {'.', '..', 'archive', 'esim', 'old'}));

datasets = struct('name', {}, 'path', {}, 'data', {});
for i = 1:length(subDirs)
    folder = fullfile(dataRoot, subDirs(i).name);
    candidate = fullfile(folder, [subDirs(i).name '_results.mat']);
    if isfile(candidate)
        try
            datasets(end+1).name = subDirs(i).name; %#ok<AGROW>
            datasets(end).path = candidate;
            datasets(end).data = load(candidate);
        catch ME
            warning('Could not load %s: %s', candidate, ME.message);
        end
    end
end

if isempty(datasets)
    error('viewAllResults: no <name>_results.mat files found under %s', dataRoot);
end

fprintf('Found %d dataset(s) under %s\n', length(datasets), dataRoot);
for i = 1:length(datasets)
    fprintf('  [%2d] %s\n', i, datasets(i).name);
end

%% ---- Build GUI ----
fig = uifigure('Name', sprintf('Detection Results - %s', dataRoot), ...
    'Position', [60 60 1500 850]);
tg = uitabgroup(fig, 'Position', [0 0 fig.Position(3) fig.Position(4)]);
fig.SizeChangedFcn = @(src,~) set(tg, 'Position', [0 0 src.Position(3) src.Position(4)]);

%% ---- Tab 1: Summary ----
summaryTab = uitab(tg, 'Title', 'Summary');
buildSummaryTab(summaryTab, datasets);

%% ---- Tabs 2..N: Per-dataset details ----
for i = 1:length(datasets)
    tab = uitab(tg, 'Title', shortenName(datasets(i).name));
    buildDetailTab(tab, datasets(i));
end

fprintf('\nOpened window with %d tabs (1 summary + %d datasets).\n', ...
    1 + length(datasets), length(datasets));
end


%% =========================================================================
%                       SUMMARY TAB
%% =========================================================================
function buildSummaryTab(parent, datasets)
    n = length(datasets);

    % Build per-dataset summary rows
    names = {datasets.name}';
    overallRates = zeros(n, 1);
    totalTicks = zeros(n, 1);
    bestWin = cell(n, 1);
    bestRate = zeros(n, 1);
    allWinMs = [];

    % Collect union of all window durations for the heatmap
    for i = 1:n
        if isfield(datasets(i).data, 'windowDurations_ms')
            allWinMs = union(allWinMs, double(datasets(i).data.windowDurations_ms));
        end
    end
    allWinMs = sort(allWinMs);

    rateMat = nan(n, length(allWinMs));  % per-dataset × per-window rate (%)

    for i = 1:n
        d = datasets(i).data;
        nTicks = length(d.tNow_us);
        totalTicks(i) = nTicks;
        overallRates(i) = 100 * sum(d.anyDetected > 0) / max(nTicks, 1);

        wins = double(d.windowDurations_ms);
        rates = zeros(size(wins));
        for wi = 1:length(wins)
            col = d.(sprintf('win_%dms', wins(wi)));
            if isfield(d, 'attemptedPerWindow')
                denom = d.attemptedPerWindow(wi);
            else
                denom = nTicks;
            end
            rates(wi) = 100 * sum(col >= 0) / max(denom, 1);

            [~, colIdx] = ismember(wins(wi), allWinMs);
            rateMat(i, colIdx) = rates(wi);
        end

        [bestRate(i), bi] = max(rates);
        bestWin{i} = sprintf('%dms', wins(bi));
    end

    % ---- Layout: table on top-left, bars on top-right, heatmap below ----
    gl = uigridlayout(parent, [2 2]);
    gl.RowHeight = {'1x', '1.2x'};
    gl.ColumnWidth = {'1x', '1x'};

    % --- Table ---
    tblPanel = uipanel(gl, 'Title', 'Per-Dataset Summary');
    tblPanel.Layout.Row = 1; tblPanel.Layout.Column = 1;

    T = table(names, totalTicks, round(overallRates, 2), bestWin, round(bestRate, 2), ...
        'VariableNames', {'Dataset', 'Ticks', 'AnyDetectPct', 'BestWindow', 'BestRatePct'});
    uit = uitable(tblPanel, 'Data', T, 'Units', 'normalized', 'Position', [0 0 1 1]);
    uit.ColumnSortable = true(1, 5);

    % --- Overall detection bar chart ---
    barPanel = uipanel(gl, 'Title', 'Overall Detection Rate per Dataset');
    barPanel.Layout.Row = 1; barPanel.Layout.Column = 2;
    ax1 = uiaxes(barPanel);
    ax1.Units = 'normalized'; ax1.Position = [0.08 0.18 0.9 0.78];
    b = bar(ax1, 1:n, overallRates);
    b.FaceColor = 'flat';
    cmap = turbo(n);
    for i = 1:n
        b.CData(i,:) = cmap(i,:);
    end
    ax1.XTick = 1:n;
    ax1.XTickLabel = cellfun(@shortenName, names, 'UniformOutput', false);
    ax1.XTickLabelRotation = 35;
    ylabel(ax1, 'Detection Rate (%)');
    title(ax1, 'Any-window detection rate');
    ylim(ax1, [0 105]);
    grid(ax1, 'on');
    for i = 1:n
        text(ax1, i, overallRates(i) + 1.5, sprintf('%.1f%%', overallRates(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end

    % --- Heatmap: dataset × window ---
    hmPanel = uipanel(gl, 'Title', 'Detection Rate (%) -- Dataset × Window');
    hmPanel.Layout.Row = 2; hmPanel.Layout.Column = [1 2];
    ax2 = uiaxes(hmPanel);
    ax2.Units = 'normalized'; ax2.Position = [0.18 0.15 0.78 0.78];
    imagesc(ax2, rateMat, 'AlphaData', ~isnan(rateMat));
    ax2.Color = [0.25 0.25 0.25];
    colormap(ax2, parula);
    cb = colorbar(ax2);
    cb.Label.String = 'Rate (%)';
    caxis(ax2, [0 100]);
    ax2.XTick = 1:length(allWinMs);
    ax2.XTickLabel = arrayfun(@(x) sprintf('%dms', x), allWinMs, 'UniformOutput', false);
    ax2.XTickLabelRotation = 45;
    ax2.YTick = 1:n;
    ax2.YTickLabel = cellfun(@shortenName, names, 'UniformOutput', false);
    xlabel(ax2, 'Window Duration');
    ylabel(ax2, 'Dataset');
    % Overlay numeric values
    for i = 1:n
        for j = 1:length(allWinMs)
            if ~isnan(rateMat(i,j))
                txtColor = 'k';
                if rateMat(i,j) < 50, txtColor = 'w'; end
                text(ax2, j, i, sprintf('%.0f', rateMat(i,j)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', txtColor);
            end
        end
    end
end


%% =========================================================================
%                       DETAIL TAB (per dataset)
%% =========================================================================
function buildDetailTab(parent, dataset)
    d = dataset.data;
    name = dataset.name;

    winMs = double(d.windowDurations_ms);
    numWindows = length(winMs);
    nTicks = length(d.tNow_us);
    tSec = (double(d.tNow_us) - double(d.tNow_us(1))) / 1e6;

    % Build detection matrix
    detMat = false(nTicks, numWindows);
    for wi = 1:numWindows
        detMat(:, wi) = d.(sprintf('win_%dms', winMs(wi))) >= 0;
    end
    anyDet = d.anyDetected > 0;

    % 2x2 layout: raster | bar / heatmap | rolling rate
    gl = uigridlayout(parent, [3 2]);
    gl.RowHeight = {35, '1x', '1x'};
    gl.ColumnWidth = {'1.2x', '1x'};

    % --- Header label with key stats ---
    hdr = uilabel(gl, 'Text', buildHeaderText(name, d, anyDet, detMat, winMs), ...
        'FontWeight', 'bold', 'FontSize', 11);
    hdr.Layout.Row = 1; hdr.Layout.Column = [1 2];

    % --- Raster plot ---
    p1 = uipanel(gl, 'Title', 'Detection Timeline (per window)');
    p1.Layout.Row = 2; p1.Layout.Column = [1 2];
    ax1 = uiaxes(p1);
    ax1.Units = 'normalized'; ax1.Position = [0.07 0.18 0.92 0.77];
    hold(ax1, 'on');
    cmap = hsv(numWindows);
    for wi = 1:numWindows
        tHits = tSec(detMat(:, wi));
        if ~isempty(tHits)
            plot(ax1, tHits, wi*ones(size(tHits)), '|', ...
                'Color', cmap(wi,:), 'MarkerSize', 4);
        end
    end
    hold(ax1, 'off');
    ax1.YTick = 1:numWindows;
    ax1.YTickLabel = arrayfun(@(x) sprintf('%dms', x), winMs, 'UniformOutput', false);
    ylim(ax1, [0.5 numWindows+0.5]);
    if ~isempty(tSec), xlim(ax1, [tSec(1) tSec(end)]); end
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Window');
    grid(ax1, 'on');

    % --- Per-window bar ---
    p2 = uipanel(gl, 'Title', 'Detection Rate per Window');
    p2.Layout.Row = 3; p2.Layout.Column = 1;
    ax2 = uiaxes(p2);
    ax2.Units = 'normalized'; ax2.Position = [0.12 0.22 0.85 0.72];
    rates = zeros(1, numWindows);
    for wi = 1:numWindows
        if isfield(d, 'attemptedPerWindow')
            denom = d.attemptedPerWindow(wi);
        else
            denom = nTicks;
        end
        rates(wi) = 100 * sum(detMat(:, wi)) / max(denom, 1);
    end
    b = bar(ax2, 1:numWindows, rates, 'FaceColor', 'flat');
    for wi = 1:numWindows, b.CData(wi,:) = cmap(wi,:); end
    ax2.XTick = 1:numWindows;
    ax2.XTickLabel = arrayfun(@(x) sprintf('%dms', x), winMs, 'UniformOutput', false);
    ax2.XTickLabelRotation = 45;
    ylabel(ax2, 'Rate (%)');
    ylim(ax2, [0 105]);
    grid(ax2, 'on');
    for wi = 1:numWindows
        text(ax2, wi, rates(wi)+2, sprintf('%.0f%%', rates(wi)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end

    % --- Rolling rate ---
    p3 = uipanel(gl, 'Title', 'Rolling Detection Rate');
    p3.Layout.Row = 3; p3.Layout.Column = 2;
    ax3 = uiaxes(p3);
    ax3.Units = 'normalized'; ax3.Position = [0.1 0.22 0.87 0.72];
    if nTicks > 1
        durS = tSec(end) - tSec(1);
        winS = max(0.5, durS / 20);
        step = tSec(2) - tSec(1);
        if step <= 0, step = 1e-3; end
        winT = max(round(winS / step), 10);
        hold(ax3, 'on');
        for wi = 1:numWindows
            rr = movmean(double(detMat(:, wi)), winT) * 100;
            plot(ax3, tSec, rr, 'Color', [cmap(wi,:) 0.5], 'LineWidth', 0.8);
        end
        ra = movmean(double(anyDet), winT) * 100;
        plot(ax3, tSec, ra, 'k-', 'LineWidth', 2);
        hold(ax3, 'off');
        xlim(ax3, [tSec(1) tSec(end)]);
    end
    ylim(ax3, [0 105]);
    xlabel(ax3, 'Time (s)');
    ylabel(ax3, 'Rate (%)');
    grid(ax3, 'on');
end


%% =========================================================================
%                       HELPERS
%% =========================================================================
function s = shortenName(name)
    % Trim common prefix for readable tab / axis labels
    s = regexprep(name, '^marker_z2_', '');
    s = strrep(s, '_', ' ');
end

function txt = buildHeaderText(name, d, anyDet, detMat, winMs)
    nTicks = length(d.tNow_us);
    nAny = sum(anyDet);
    dur = (double(d.tNow_us(end)) - double(d.tNow_us(1))) / 1e6;

    % Best window
    rates = zeros(1, length(winMs));
    for wi = 1:length(winMs)
        if isfield(d, 'attemptedPerWindow')
            denom = d.attemptedPerWindow(wi);
        else
            denom = nTicks;
        end
        rates(wi) = 100 * sum(detMat(:, wi)) / max(denom, 1);
    end
    [bestRate, bi] = max(rates);

    txt = sprintf(['%s   |   %d ticks (%.2fs)   |   Any-window: %d (%.1f%%)   |   ' ...
        'Best: %dms at %.1f%%   |   Windows: %s ms'], ...
        name, nTicks, dur, nAny, 100*nAny/max(nTicks,1), ...
        winMs(bi), bestRate, mat2str(winMs));
end
