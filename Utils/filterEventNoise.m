function outPath = filterEventNoise(matFile, opts)
%FILTEREVENTNOISE  Clean an event-camera recording before detection.
%
%  outPath = filterEventNoise(matFile)
%  outPath = filterEventNoise(matFile, opts)
%  outPaths = filterEventNoise(["a.mat","b.mat"], opts)        % batch
%
%  Pipeline (each stage opt-in via opts.method):
%
%    'hot'   - hot-pixel mask. Per-pixel event count over the whole
%              recording; pixels whose count exceeds hotMult * median
%              get dropped entirely. Cheap; removes broken pixels.
%
%    'bg'    - background-activity filter. For each event at (x,y,t),
%              look at the r x r box around (x,y); keep it only if
%              >= minNeighbors of those pixels fired within dt_ms of t.
%              Removes isolated thermal noise. This is your "no event
%              in 10 ms -> noise" idea, generalised.
%
%    'hot+bg' (default) - hot-pixel mask then BAF. Recommended.
%
%    'refractory' - drop events that fire less than tauRef_us after the
%              previous event at the same pixel. Kills hot-pixel bursts
%              that survive the hot-pixel mask.
%
%  Inputs:
%    matFile : path (string/char) or array of paths
%    opts    : struct with optional fields
%      .method        (default 'hot+bg')  'hot' | 'bg' | 'hot+bg' | 'refractory'
%      .dt_ms         (default 10)        BAF time window in ms
%      .radius        (default 1)         BAF spatial radius (1 = 3x3)
%      .minNeighbors  (default 1)         BAF: require >= N neighbours
%      .hotMult       (default 10)        hot-pixel threshold (rate >= mult*median)
%      .tauRef_us     (default 1000)      refractory window (us)
%      .sensorSize    (default inferred)  [H W]
%      .show          (default true)      side-by-side preview
%      .windowMs      (default 10)        accumulation window in preview
%      .numFrames     (default 6)         preview slices
%      .save          (default true)      write <name>_filtered.mat
%      .outPath       (default '')        explicit save path (single file)
%
%  Output:
%    outPath : path of the saved <name>_filtered.mat (or cell array if
%              matFile was an array).
%
%  Example:
%    filterEventNoise('Data/.../1.mat');                            % defaults
%    filterEventNoise('Data/.../1.mat', struct('dt_ms', 5));        % stricter
%    filterEventNoise(["Data/.../1.mat","Data/.../2.mat"], ...
%        struct('show', false, 'method', 'hot+bg'));                % batch

if nargin < 2, opts = struct(); end
if ~isfield(opts,'method'),       opts.method       = 'hot+bg'; end
if ~isfield(opts,'dt_ms'),        opts.dt_ms        = 10;       end
if ~isfield(opts,'radius'),       opts.radius       = 1;        end
if ~isfield(opts,'minNeighbors'), opts.minNeighbors = 1;        end
if ~isfield(opts,'hotMult'),      opts.hotMult      = 10;       end
if ~isfield(opts,'tauRef_us'),    opts.tauRef_us    = 1000;     end
if ~isfield(opts,'sensorSize'),   opts.sensorSize   = [];       end
if ~isfield(opts,'show'),         opts.show         = true;     end
if ~isfield(opts,'windowMs'),     opts.windowMs     = 10;       end
if ~isfield(opts,'numFrames'),    opts.numFrames    = 6;        end
if ~isfield(opts,'save'),         opts.save         = true;     end
if ~isfield(opts,'outPath'),      opts.outPath      = '';       end

%% ---- Batch dispatch ----
matFile = string(matFile);
matFile = matFile(:);
if numel(matFile) > 1
    if ~isempty(opts.outPath)
        error('filterEventNoise: opts.outPath is only valid for a single input.');
    end
    outPath = strings(numel(matFile), 1);
    for k = 1:numel(matFile)
        fprintf('\n=== [%d/%d] %s ===\n', k, numel(matFile), matFile(k));
        outPath(k) = filterEventNoise(matFile(k), opts);
    end
    outPath = cellstr(outPath);
    return;
end
matFile = char(matFile);

if ~isfile(matFile)
    error('filterEventNoise: file not found: %s', matFile);
end

%% ---- Load ----
fprintf('Loading %s...\n', matFile);
allVars = load(matFile);
if ~isfield(allVars, 'events')
    error('filterEventNoise: %s has no `events` variable.', matFile);
end
events = allVars.events;
N0 = size(events, 1);
fprintf('  %d events loaded\n', N0);
if size(events, 2) < 4
    error('events must be Nx4 (x, y, polarity, t)');
end

% Sensor size
if ~isempty(opts.sensorSize)
    H = opts.sensorSize(1); W = opts.sensorSize(2);
else
    H = double(max(events(:,2))) + 1;
    W = double(max(events(:,1))) + 1;
end
fprintf('  sensor: %d x %d (HxW)\n', H, W);

% Make sure events are sorted by time -- every filter assumes this.
if ~issorted(events(:,4))
    [~, si] = sort(events(:,4));
    events = events(si, :);
    fprintf('  sorted events by timestamp\n');
end

%% ---- Apply filters ----
method = lower(string(opts.method));
keep = true(N0, 1);

tFilter = tic;

if any(method == "hot") || method == "hot+bg"
    [keepHot, hotMask] = applyHotPixelMask(events, [H W], opts.hotMult);
    nDropped = sum(~keepHot);
    fprintf('  [hot ] dropped %d events  (%.2f%%)  on %d hot pixels\n', ...
        nDropped, 100*nDropped/N0, sum(hotMask(:)));
    keep = keep & keepHot;
end

if method == "bg" || method == "hot+bg"
    events_for_bg = events(keep, :);
    keepBgLocal = applyBackgroundActivityFilter( ...
        events_for_bg, [H W], opts.dt_ms*1000, opts.radius, opts.minNeighbors);
    keepBgFull = false(N0, 1);
    keepBgFull(find(keep, sum(keep))) = keepBgLocal;
    nDropped = sum(keep) - sum(keepBgFull & keep);
    fprintf('  [BAF ] dropped %d events  (%.2f%%)  dt=%dms  r=%d  minN=%d\n', ...
        nDropped, 100*nDropped/N0, opts.dt_ms, opts.radius, opts.minNeighbors);
    keep = keep & keepBgFull;
end

if method == "refractory"
    keepRef = applyRefractoryFilter(events, [H W], opts.tauRef_us);
    nDropped = sum(~keepRef);
    fprintf('  [refr] dropped %d events  (%.2f%%)  tauRef=%dus\n', ...
        nDropped, 100*nDropped/N0, opts.tauRef_us);
    keep = keep & keepRef;
end

filtered = events(keep, :);
fprintf('  ---- summary ----\n');
fprintf('  kept %d / %d events  (%.2f%%)  in %.1fs\n', ...
    size(filtered,1), N0, 100*size(filtered,1)/N0, toc(tFilter));

%% ---- Decide output path ----
if isempty(opts.outPath)
    [d, n] = fileparts(matFile);
    outPath = fullfile(d, [n '_filtered.mat']);
else
    outPath = opts.outPath;
end

%% ---- Save ----
if opts.save
    allVars.events    = filtered;
    allVars.noiseFilter = struct( ...
        'method',        char(method), ...
        'dt_ms',         opts.dt_ms, ...
        'radius',        opts.radius, ...
        'minNeighbors',  opts.minNeighbors, ...
        'hotMult',       opts.hotMult, ...
        'tauRef_us',     opts.tauRef_us, ...
        'eventsIn',      N0, ...
        'eventsOut',     size(filtered,1), ...
        'keepRatio',     size(filtered,1) / max(N0,1), ...
        'sensorSize',    [H W], ...
        'sourceFile',    matFile, ...
        'createdAt',     datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<DATST,TNOW1>

    % Force MAT v7.3 -- the events array easily exceeds the 2 GB limit
    % of older MAT versions on 480x640 recordings.
    save(outPath, '-struct', 'allVars', '-v7.3');

    info = dir(outPath);
    fprintf('  saved -> %s  (%s)\n', outPath, humanBytes(info.bytes));
end

%% ---- Visualise ----
if opts.show
    showComparison(events, filtered, H, W, opts.windowMs, opts.numFrames, matFile);
end
end


%% =========================================================================
%                       HOT-PIXEL MASK
%% =========================================================================
function [keep, hotMask] = applyHotPixelMask(events, sensorSize, mult)
    H = sensorSize(1); W = sensorSize(2);
    x = double(events(:,1)) + 1;
    y = double(events(:,2)) + 1;
    valid = (x>=1) & (x<=W) & (y>=1) & (y<=H);
    counts = accumarray([y(valid), x(valid)], 1, [H W]);

    activeCounts = counts(counts > 0);
    if isempty(activeCounts)
        keep = true(size(events,1), 1);
        hotMask = false(H, W);
        return;
    end
    threshold = mult * median(activeCounts);
    hotMask = counts > threshold;

    % Drop events that originate from any hot pixel
    pixIdx = sub2ind([H W], min(max(y,1),H), min(max(x,1),W));
    keep = ~hotMask(pixIdx);
    keep = keep & valid;
end


%% =========================================================================
%                       BACKGROUND ACTIVITY FILTER
%% =========================================================================
function keep = applyBackgroundActivityFilter(events, sensorSize, dt_us, r, minNeighbors)
    % events must already be sorted by time (column 4).
    H = sensorSize(1); W = sensorSize(2);
    N = size(events, 1);

    % Padded last-event-time map so we can read an (2r+1) x (2r+1) box
    % without bounds checks. Pad with -inf so unset pixels never count
    % as a neighbour.
    lastT = -inf(H + 2*r, W + 2*r);

    keep = false(N, 1);
    xs = double(events(:,1)) + 1 + r;   % padded 1-based column
    ys = double(events(:,2)) + 1 + r;   % padded 1-based row
    ts = double(events(:,4));

    threshold = double(dt_us);
    rr = r;

    for i = 1:N
        x = xs(i); y = ys(i); t = ts(i);
        nbhd = lastT(y-rr:y+rr, x-rr:x+rr);
        % We don't count the pixel itself as its own neighbour; mask
        % the centre so the next-event-at-same-pixel doesn't trivially
        % keep an isolated burst alive.
        nbhd(rr+1, rr+1) = -inf;
        if nnz(nbhd >= t - threshold) >= minNeighbors
            keep(i) = true;
        end
        lastT(y, x) = t;
    end
end


%% =========================================================================
%                       REFRACTORY FILTER
%% =========================================================================
function keep = applyRefractoryFilter(events, sensorSize, tauRef_us)
    H = sensorSize(1); W = sensorSize(2);
    N = size(events, 1);
    lastT = -inf(H, W);
    keep = true(N, 1);
    xs = double(events(:,1)) + 1;
    ys = double(events(:,2)) + 1;
    ts = double(events(:,4));
    for i = 1:N
        x = xs(i); y = ys(i); t = ts(i);
        if x < 1 || x > W || y < 1 || y > H, keep(i) = false; continue; end
        if t - lastT(y, x) < tauRef_us
            keep(i) = false;
        else
            lastT(y, x) = t;
        end
    end
end


%% =========================================================================
%                       VISUAL COMPARISON
%% =========================================================================
function showComparison(eventsRaw, eventsClean, H, W, winMs, nFrames, fileA)
    tRaw = eventsRaw(:,4);
    tMin = double(min(tRaw));
    tMax = double(max(tRaw));
    winUs = winMs * 1000;

    if tMax - tMin < winUs
        warning('Recording shorter than the window; using one slice.');
        tNowList = (tMin + tMax) / 2;
        nFrames = 1;
    else
        tNowList = linspace(tMin + winUs, tMax, nFrames);
    end

    [~, baseName] = fileparts(fileA);
    fig = figure('Name', sprintf('Noise filter preview: %s', baseName), ...
                 'Color', 'w', ...
                 'Position', [60 60 250*nFrames + 80, 620]);
    tl = tiledlayout(fig, 2, nFrames, 'TileSpacing','compact','Padding','compact');

    for i = 1:nFrames
        tNow = tNowList(i);

        nexttile(tl, i);
        imgA = accumulateWindow(eventsRaw, tNow, winUs, H, W);
        imshow(autoStretch(imgA), 'InitialMagnification', 'fit');
        title(sprintf('t = %.3f s', tNow / 1e6), 'FontSize', 10);
        if i == 1
            ylabel('raw', 'FontWeight','bold','Visible','on', ...
                   'Rotation', 0, 'HorizontalAlignment', 'right');
        end

        nexttile(tl, nFrames + i);
        imgB = accumulateWindow(eventsClean, tNow, winUs, H, W);
        imshow(autoStretch(imgB), 'InitialMagnification', 'fit');
        if i == 1
            ylabel('filtered', 'FontWeight','bold','Visible','on', ...
                   'Rotation', 0, 'HorizontalAlignment', 'right');
        end
    end
    title(tl, sprintf('Noise filter  -  raw vs filtered  -  %d ms window', winMs), ...
          'FontSize', 12, 'FontWeight', 'bold');
end


%% =========================================================================
function img = accumulateWindow(events, tNow, winUs, H, W)
    tCol = events(:, 4);
    mask = (tCol > (tNow - winUs)) & (tCol <= tNow);
    if ~any(mask)
        img = zeros(H, W); return;
    end
    x = double(events(mask, 1)) + 1;
    y = double(events(mask, 2)) + 1;
    valid = (x>=1) & (x<=W) & (y>=1) & (y<=H);
    if ~any(valid), img = zeros(H, W); return; end
    img = accumarray([y(valid), x(valid)], 1, [H, W]);
end


%% =========================================================================
function out = autoStretch(img)
    img = double(img);
    mx = max(img(:));
    if mx == 0
        out = uint8(zeros(size(img)));
    else
        out = uint8(img / mx * 255);
    end
end


%% =========================================================================
function s = humanBytes(n)
    units = {'B','KB','MB','GB','TB'};
    i = 1;
    while n >= 1024 && i < length(units)
        n = n / 1024; i = i + 1;
    end
    s = sprintf('%.2f %s', n, units{i});
end
