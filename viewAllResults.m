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

% Sort datasets by (motion, speed) so speeds appear as low -> med -> high
sortKeys = cell(length(datasets), 1);
for i = 1:length(datasets)
    sortKeys{i} = datasetSortKey(datasets(i).name);
end
[~, order] = sort(sortKeys);
datasets = datasets(order);

fprintf('Found %d dataset(s) under %s\n', length(datasets), dataRoot);
for i = 1:length(datasets)
    fprintf('  [%2d] %s\n', i, datasets(i).name);
end

%% ---- Build GUI ----
fig = uifigure('Name', sprintf('Detection Results - %s', dataRoot), ...
    'Position', [60 60 1500 850]);
% Use a grid layout so the tabgroup resizes with the figure automatically
% (avoids the SizeChangedFcn / AutoResizeChildren conflict warning).
figGL = uigridlayout(fig, [2 1]);
figGL.RowHeight   = {32, '1x'};
figGL.Padding     = [4 4 4 4];
figGL.RowSpacing  = 4;

% --- Toolbar (save current tab / full window) ---
tb = uigridlayout(figGL, [1 3]);
tb.Layout.Row     = 1;
tb.ColumnWidth    = {140, 160, '1x'};
tb.Padding        = [0 0 0 0];
tb.ColumnSpacing  = 6;
btnTab = uibutton(tb, 'Text', 'Save Current Tab...', ...
    'ButtonPushedFcn', @(~,~) saveCurrentTab(fig, tg));
btnTab.Layout.Column = 1;
btnAll = uibutton(tb, 'Text', 'Save Whole Window...', ...
    'ButtonPushedFcn', @(~,~) saveWholeWindow(fig));
btnAll.Layout.Column = 2;

tg = uitabgroup(figGL);
tg.Layout.Row = 2;

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
    barGL = uigridlayout(barPanel, [1 1]);
    barGL.Padding = [5 5 5 5];
    ax1 = uiaxes(barGL);
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
    hmGL = uigridlayout(hmPanel, [1 1]);
    hmGL.Padding = [5 5 5 5];
    ax2 = uiaxes(hmGL);
    imagesc(ax2, rateMat, 'AlphaData', ~isnan(rateMat));
    ax2.Color = [0.25 0.25 0.25];
    colormap(ax2, parula);
    cb = colorbar(ax2);
    cb.Label.String = 'Rate (%)';
    clim(ax2, [0 100]);
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

    % Layout: header, raster, per-window bar
    gl = uigridlayout(parent, [3 1]);
    gl.RowHeight = {35, '1.2x', '1x'};

    % --- Header label with key stats ---
    hdr = uilabel(gl, 'Text', buildHeaderText(name, d, anyDet, detMat, winMs), ...
        'FontWeight', 'bold', 'FontSize', 11);
    hdr.Layout.Row = 1;

    % --- Raster plot ---
    p1 = uipanel(gl, 'Title', 'Detection Timeline (per window)');
    p1.Layout.Row = 2;
    p1GL = uigridlayout(p1, [1 1]);
    p1GL.Padding = [5 5 5 5];
    ax1 = uiaxes(p1GL);
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
    p2 = uipanel(gl, 'Title', 'Detection Rate per Window (rate = detections / attempts for that window)');
    p2.Layout.Row = 3;
    p2GL = uigridlayout(p2, [1 1]);
    p2GL.Padding = [5 5 5 5];
    ax2 = uiaxes(p2GL);
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
end


%% =========================================================================
%                       HELPERS
%% =========================================================================
function s = shortenName(name)
    % Trim common prefix for readable tab / axis labels
    s = regexprep(name, '^marker_z2_', '');
    s = strrep(s, '_', ' ');
end

function k = datasetSortKey(name)
    % Build a sort key that orders datasets by (motion, speed)
    % where speed is forced to the sequence low -> med -> high.
    base = regexprep(name, '^marker_z2_', '');
    tok = regexp(base, '^(.*)_(low|med|high)$', 'tokens', 'once');
    if isempty(tok)
        k = ['z_' base];           % unknown pattern: sort to the end
        return;
    end
    motion = tok{1};
    speed  = tok{2};
    switch speed
        case 'low',  s = '1';
        case 'med',  s = '2';
        case 'high', s = '3';
        otherwise,   s = '9';
    end
    k = [motion '_' s];
end

function saveWholeWindow(fig)
    % Save the entire uifigure (all tabs composited as currently shown) to a file.
    filters = {'*.png','PNG image (*.png)'; ...
               '*.pdf','PDF document (*.pdf)'; ...
               '*.jpg','JPEG image (*.jpg)'};
    [file, path] = uiputfile(filters, 'Save window as', 'detection_results.png');
    if isequal(file, 0), return; end
    outPath = fullfile(path, file);
    try
        exportapp(fig, outPath);
        uialert(fig, sprintf('Saved to:\n%s', outPath), 'Saved', 'Icon', 'success');
    catch ME
        uialert(fig, sprintf('Could not save: %s', ME.message), 'Save failed');
    end
end

function saveCurrentTab(fig, tg)
    % Save only the currently selected tab content.
    currentTab = tg.SelectedTab;
    tabTitle   = matlab.lang.makeValidName(currentTab.Title);
    filters = {'*.png','PNG image (*.png)'; ...
               '*.pdf','PDF document (*.pdf)'; ...
               '*.jpg','JPEG image (*.jpg)'};
    defaultName = ['detection_' tabTitle '.png'];
    [file, path] = uiputfile(filters, 'Save current tab as', defaultName);
    if isequal(file, 0), return; end
    outPath = fullfile(path, file);
    try
        exportapp(currentTab, outPath);
        uialert(fig, sprintf('Saved to:\n%s', outPath), 'Saved', 'Icon', 'success');
    catch ME
        uialert(fig, sprintf('Could not save: %s', ME.message), 'Save failed');
    end
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
