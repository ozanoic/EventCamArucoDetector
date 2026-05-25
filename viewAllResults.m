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

% --- Toolbar (marker filter + save buttons) ---
tb = uigridlayout(figGL, [1 6]);
tb.Layout.Row     = 1;
tb.ColumnWidth    = {65, 130, 160, 170, 215, '1x'};
tb.Padding        = [0 0 0 0];
tb.ColumnSpacing  = 6;

% Discover every marker id any window ever decoded into.
markerIds = discoverMarkerIds(datasets);

% Tabgroup is created first so save-button closures can capture it.
tg = uitabgroup(figGL);
tg.Layout.Row = 2;

lblMarker = uilabel(tb, 'Text', 'Marker:', 'HorizontalAlignment', 'right');
lblMarker.Layout.Column = 1;
markerItems    = ['All markers', arrayfun(@(x) sprintf('ID %d', x), markerIds, ...
                                          'UniformOutput', false)];
markerItemData = [{[]}, num2cell(markerIds)];
ddMarker = uidropdown(tb, 'Items', markerItems, 'ItemsData', markerItemData);
ddMarker.Layout.Column = 2;

btnTab = uibutton(tb, 'Text', 'Save Current Tab (PDF)...', ...
    'ButtonPushedFcn', @(~,~) saveCurrentTab(fig, tg));
btnTab.Layout.Column = 3;
btnAll = uibutton(tb, 'Text', 'Save Whole Window (PDF)...', ...
    'ButtonPushedFcn', @(~,~) saveWholeWindow(fig));
btnAll.Layout.Column = 4;
btnEvery = uibutton(tb, 'Text', 'Save All Tabs (multi-page PDF)...', ...
    'ButtonPushedFcn', @(~,~) saveAllTabs(fig, tg));
btnEvery.Layout.Column = 5;

% (Re)build every tab for the currently selected marker filter.
ddMarker.ValueChangedFcn = @(~,~) buildAllTabs(tg, datasets, ddMarker.Value);
buildAllTabs(tg, datasets, ddMarker.Value);   % initial build

fprintf('\nOpened window with %d tabs (1 summary + %d datasets).\n', ...
    1 + length(datasets), length(datasets));
end


%% =========================================================================
%                       BUILD / REBUILD ALL TABS
%% =========================================================================
function buildAllTabs(tg, datasets, markerFilter)
    % markerFilter: [] -> all markers; otherwise scalar marker ID.
    delete(tg.Children);
    titleSuffix = markerLabelSuffix(markerFilter);

    summaryTab = uitab(tg, 'Title', ['Summary' titleSuffix]);
    buildSummaryTab(summaryTab, datasets, markerFilter);

    for i = 1:length(datasets)
        tab = uitab(tg, 'Title', shortenName(datasets(i).name));
        buildDetailTab(tab, datasets(i), markerFilter);
    end
end

function s = markerLabelSuffix(markerFilter)
    if isempty(markerFilter)
        s = '';
    else
        s = sprintf(' (ID %d)', markerFilter);
    end
end

function ids = discoverMarkerIds(datasets)
    % Scan every win_*ms column across every dataset and return the sorted
    % unique set of marker IDs that were actually decoded (>= 0).
    seen = [];
    for i = 1:length(datasets)
        d = datasets(i).data;
        if ~isfield(d, 'windowDurations_ms'), continue; end
        wins = double(d.windowDurations_ms);
        for wi = 1:length(wins)
            col = d.(sprintf('win_%dms', wins(wi)));
            v = double(col(col >= 0));
            if ~isempty(v), seen = [seen; v]; end %#ok<AGROW>
        end
    end
    ids = unique(seen(:))';
end

function mask = perWindowPerMarkerMask(d, windowMs, markerFilter)
    % Multi-marker-aware mask for a single (window, marker) cell.
    %
    %  markerFilter == [] : any marker -> read scalar win_<X>ms (>= 0)
    %  markerFilter == M  : prefer win_<X>ms_idM (logical per (tick,window))
    %                       which captures markers that the scalar column
    %                       would miss when several markers shared a tick.
    %                       Falls back to win_<X>ms == M when the new field
    %                       is absent (older result files).
    scalarField = sprintf('win_%dms', windowMs);
    if isempty(markerFilter)
        mask = d.(scalarField)(:) >= 0;
        return;
    end
    perField = sprintf('win_%dms_id%d', windowMs, markerFilter);
    if isfield(d, perField)
        mask = d.(perField)(:) > 0;
    else
        mask = d.(scalarField)(:) == markerFilter;
    end
end

function ids = collectReportIds(datasets)
    % Returns the set of marker IDs that should get their own column /
    % series in the summary view. Prefers the canonical fields written
    % by detectAruco/mergeResults; falls back to discoverMarkerIds.
    reported = [];
    requested = [];
    for i = 1:length(datasets)
        d = datasets(i).data;
        if isfield(d, 'markerIdsReported')
            reported = union(reported, double(d.markerIdsReported(:)'));
        end
        if isfield(d, 'requestedMarkerIds')
            requested = union(requested, double(d.requestedMarkerIds(:)'));
        end
    end
    if ~isempty(reported)
        ids = sort(reported(:)');
    elseif ~isempty(requested)
        ids = sort(requested(:)');
    else
        ids = discoverMarkerIds(datasets);
    end
end

function anyDet = perMarkerAnyDetect(d, mid)
    % Returns an Nx1 logical: was marker `mid` decoded by ANY window at
    % each tick?  Resolution order:
    %   (1) anyDetected_id<mid> -- written by detectAruco / mergeResults
    %   (2) win_<X>ms_id<mid>   -- per-window per-marker columns
    %   (3) win_<X>ms == mid    -- legacy scalar column (loses some hits)
    fname = sprintf('anyDetected_id%d', mid);
    if isfield(d, fname)
        anyDet = d.(fname) > 0;
        return;
    end
    nTicks = length(d.tNow_us);
    anyDet = false(nTicks, 1);
    wins = double(d.windowDurations_ms);
    for wi = 1:length(wins)
        perField = sprintf('win_%dms_id%d', wins(wi), mid);
        if isfield(d, perField)
            anyDet = anyDet | (d.(perField)(:) > 0);
        else
            col = d.(sprintf('win_%dms', wins(wi)));
            anyDet = anyDet | (col(:) == mid);
        end
    end
end


%% =========================================================================
%                       SUMMARY TAB
%% =========================================================================
function buildSummaryTab(parent, datasets, markerFilter)
    if nargin < 3, markerFilter = []; end
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

    % Discover the set of "interesting" marker IDs so we can show one
    % column per marker in the summary table.  Preference order:
    %   (1) markerIdsReported (set by detectAruco/mergeResults), if every
    %       dataset has it - that's the canonical answer.
    %   (2) requestedMarkerIds union across datasets.
    %   (3) every ID actually decoded in the data.
    perMarkerIds = collectReportIds(datasets);

    perMarkerPct = nan(n, length(perMarkerIds));

    for i = 1:n
        d = datasets(i).data;
        nTicks = length(d.tNow_us);
        totalTicks(i) = nTicks;

        % Any-window detection rate, filtered by marker
        wins = double(d.windowDurations_ms);
        winMask = false(nTicks, length(wins));
        for wi = 1:length(wins)
            winMask(:, wi) = perWindowPerMarkerMask(d, wins(wi), markerFilter);
        end
        overallRates(i) = 100 * sum(any(winMask, 2)) / max(nTicks, 1);

        rates = zeros(size(wins));
        for wi = 1:length(wins)
            if isfield(d, 'attemptedPerWindow')
                denom = d.attemptedPerWindow(wi);
            else
                denom = nTicks;
            end
            rates(wi) = 100 * sum(winMask(:, wi)) / max(denom, 1);

            [~, colIdx] = ismember(wins(wi), allWinMs);
            rateMat(i, colIdx) = rates(wi);
        end

        [bestRate(i), bi] = max(rates);
        bestWin{i} = sprintf('%dms', wins(bi));

        % Per-marker any-detect % (independent of the dropdown filter)
        for mi = 1:length(perMarkerIds)
            mid = perMarkerIds(mi);
            anyForId = perMarkerAnyDetect(d, mid);
            perMarkerPct(i, mi) = 100 * sum(anyForId) / max(nTicks, 1);
        end
    end

    % ---- Layout: table on top-left, bars on top-right, heatmap below ----
    gl = uigridlayout(parent, [2 2]);
    gl.RowHeight = {'1x', '1.2x'};
    gl.ColumnWidth = {'1x', '1x'};

    % --- Table ---
    tblPanel = uipanel(gl, 'Title', 'Per-Dataset Summary');
    tblPanel.Layout.Row = 1; tblPanel.Layout.Column = 1;

    varNames = {'Dataset', 'Ticks', 'AnyDetectPct', 'BestWindow', 'BestRatePct'};
    tableCols = {names, totalTicks, round(overallRates, 2), bestWin, round(bestRate, 2)};
    for mi = 1:length(perMarkerIds)
        varNames{end+1} = sprintf('Id%dPct', perMarkerIds(mi)); %#ok<AGROW>
        tableCols{end+1} = round(perMarkerPct(:, mi), 2);       %#ok<AGROW>
    end
    T = table(tableCols{:}, 'VariableNames', varNames);
    uit = uitable(tblPanel, 'Data', T, 'Units', 'normalized', 'Position', [0 0 1 1]);
    uit.ColumnSortable = true(1, width(T));

    % --- Overall detection bar chart -----------------------------------
    % If there are multiple reported markers, draw one grouped bar per
    % marker so you can read e.g. "marker 3 vs marker 8" side-by-side
    % at a glance. Otherwise fall back to the simple one-bar-per-dataset
    % view.
    barPanel = uipanel(gl, 'Title', 'Overall Detection Rate per Dataset');
    barPanel.Layout.Row = 1; barPanel.Layout.Column = 2;
    barGL = uigridlayout(barPanel, [1 1]);
    barGL.Padding = [5 5 5 5];
    ax1 = uiaxes(barGL);

    if length(perMarkerIds) >= 2 && isempty(markerFilter)
        % Grouped bars: rows = datasets, columns = markers
        b = bar(ax1, 1:n, perMarkerPct, 'grouped');
        legendLabels = arrayfun(@(x) sprintf('ID %d', x), perMarkerIds, ...
                                'UniformOutput', false);
        legend(ax1, legendLabels, 'Location', 'best');
        title(ax1, 'Any-window detection rate per marker');
    else
        b = bar(ax1, 1:n, overallRates);
        b.FaceColor = 'flat';
        cmap = turbo(n);
        for i = 1:n, b.CData(i,:) = cmap(i,:); end
        if isempty(markerFilter)
            title(ax1, 'Any-window detection rate');
        else
            title(ax1, sprintf('Any-window detection rate (marker ID %d)', markerFilter));
        end
        for i = 1:n
            text(ax1, i, overallRates(i) + 1.5, sprintf('%.1f%%', overallRates(i)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end
    ax1.XTick = 1:n;
    ax1.XTickLabel = cellfun(@shortenName, names, 'UniformOutput', false);
    ax1.XTickLabelRotation = 35;
    ylabel(ax1, 'Detection Rate (%)');
    ylim(ax1, [0 105]);
    grid(ax1, 'on');

    % --- Heatmap: dataset × window ---
    if isempty(markerFilter)
        hmTitle = 'Detection Rate (%) -- Dataset × Window  (any marker)';
    else
        hmTitle = sprintf('Detection Rate (%%) -- Dataset × Window  (marker ID %d only)', markerFilter);
    end
    hmPanel = uipanel(gl, 'Title', hmTitle);
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
function buildDetailTab(parent, dataset, markerFilter)
    if nargin < 3, markerFilter = []; end
    d = dataset.data;
    name = dataset.name;

    winMs = double(d.windowDurations_ms);
    numWindows = length(winMs);
    nTicks = length(d.tNow_us);
    tSec = (double(d.tNow_us) - double(d.tNow_us(1))) / 1e6;

    % Build detection matrix (filtered by marker selection)
    detMat = false(nTicks, numWindows);
    for wi = 1:numWindows
        detMat(:, wi) = perWindowPerMarkerMask(d, winMs(wi), markerFilter);
    end
    anyDet = any(detMat, 2);

    % Layout: header, raster, per-window bar
    gl = uigridlayout(parent, [3 1]);
    gl.RowHeight = {35, '1.2x', '1x'};

    % --- Header label with key stats ---
    hdr = uilabel(gl, 'Text', buildHeaderText(name, d, anyDet, detMat, winMs, markerFilter), ...
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
    % Save the entire uifigure (only the currently visible tab is captured;
    % use "Save All Tabs" to get every tab).  PDF is the default.
    filters = {'*.pdf','PDF document (*.pdf)'; ...
               '*.png','PNG image (*.png)'; ...
               '*.jpg','JPEG image (*.jpg)'};
    [file, path] = uiputfile(filters, 'Save window as', 'detection_results.pdf');
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
    % Save only the currently selected tab content.  PDF by default.
    currentTab = tg.SelectedTab;
    tabTitle   = matlab.lang.makeValidName(currentTab.Title);
    filters = {'*.pdf','PDF document (*.pdf)'; ...
               '*.png','PNG image (*.png)'; ...
               '*.jpg','JPEG image (*.jpg)'};
    defaultName = ['detection_' tabTitle '.pdf'];
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

function saveAllTabs(fig, tg)
    % Save every tab as one high-quality multi-page PDF.
    %
    % For each tab we briefly select it (so the layout is realised) and
    % then append it to the output PDF with exportgraphics(..., 'Append',
    % true), which produces vector output where possible. If the running
    % MATLAB cannot append a uitab container to a PDF, we fall back to
    % writing one PDF per tab into a sibling folder.
    [file, path] = uiputfile({'*.pdf','PDF document (*.pdf)'}, ...
        'Save all tabs as multi-page PDF', 'detection_results_all.pdf');
    if isequal(file, 0), return; end
    outPath = fullfile(path, file);
    if isfile(outPath), delete(outPath); end

    tabs    = tg.Children;
    nTabs   = numel(tabs);
    origTab = tg.SelectedTab;
    progressDlg = uiprogressdlg(fig, 'Title', 'Saving PDF', ...
        'Message', sprintf('Rendering tab 1 of %d...', nTabs), ...
        'Cancelable', 'on');
    restoreTab = onCleanup(@() restoreSelectedTab(tg, origTab)); %#ok<NASGU>
    closeProg  = onCleanup(@() delete(progressDlg));              %#ok<NASGU>

    appendOk = true;
    try
        for i = 1:nTabs
            if progressDlg.CancelRequested, return; end
            progressDlg.Value   = (i-1) / nTabs;
            progressDlg.Message = sprintf('Rendering tab %d of %d: %s', ...
                i, nTabs, tabs(i).Title);

            tg.SelectedTab = tabs(i);
            drawnow;

            try
                exportgraphics(tabs(i), outPath, ...
                    'ContentType', 'vector', ...
                    'Append',      i > 1);
            catch
                % Older MATLAB or unsupported container -> per-tab files.
                appendOk = false;
                break;
            end
        end

        if appendOk
            progressDlg.Value = 1;
            uialert(fig, sprintf('Saved %d tabs to:\n%s', nTabs, outPath), ...
                'Saved', 'Icon', 'success');
            return;
        end

        % --- Fallback: render each tab via exportapp into a sibling folder
        if isfile(outPath), delete(outPath); end
        [outDir, base] = fileparts(outPath);
        pagesDir = fullfile(outDir, [base '_pages']);
        if ~isfolder(pagesDir), mkdir(pagesDir); end
        for i = 1:nTabs
            if progressDlg.CancelRequested, return; end
            progressDlg.Value   = (i-1) / nTabs;
            progressDlg.Message = sprintf('Fallback: rendering tab %d of %d', i, nTabs);
            tg.SelectedTab = tabs(i);
            drawnow;
            safeTitle = matlab.lang.makeValidName(tabs(i).Title);
            dst = fullfile(pagesDir, sprintf('%02d_%s.pdf', i, safeTitle));
            exportapp(tabs(i), dst);
        end
        uialert(fig, sprintf(['Multi-page PDF append not supported here.\n' ...
            'Saved %d individual tab PDFs in:\n%s'], nTabs, pagesDir), ...
            'Saved (per-tab)', 'Icon', 'warning');
    catch ME
        uialert(fig, sprintf('Could not save: %s', ME.message), 'Save failed');
    end
end

function restoreSelectedTab(tg, tab)
    try
        tg.SelectedTab = tab;
    catch
        % ignore -- figure may already be closed
    end
end

function txt = buildHeaderText(name, d, anyDet, detMat, winMs, markerFilter)
    if nargin < 6, markerFilter = []; end
    nTicks = length(d.tNow_us);
    nAny = sum(anyDet);
    dur = (double(d.tNow_us(end)) - double(d.tNow_us(1))) / 1e6;

    % Best window (for this filter)
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

    if isempty(markerFilter)
        filterStr = 'any marker';
    else
        filterStr = sprintf('marker ID %d', markerFilter);
    end

    txt = sprintf(['%s   |   %s   |   %d ticks (%.2fs)   |   Any-window: %d (%.1f%%)   |   ' ...
        'Best: %dms at %.1f%%   |   Windows: %s ms'], ...
        name, filterStr, nTicks, dur, nAny, 100*nAny/max(nTicks,1), ...
        winMs(bi), bestRate, mat2str(winMs));
end
