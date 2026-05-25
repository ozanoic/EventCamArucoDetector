function diagnoseFailures(resultFile, eventFile, opts)
%DIAGNOSEFAILURES  Visualise ticks where the detector missed.
%
%  diagnoseFailures(resultFile, eventFile)
%  diagnoseFailures(resultFile, eventFile, opts)
%
%  Loads a result struct + its source event file, picks a handful of
%  ticks where `anyDetected == 0`, and renders a grid where each panel
%  shows the accumulated event image at `opts.windowMs` ms and the
%  candidate quads the blob detector found there.
%
%  Use it to answer: "on a failed tick, is the real marker shape
%  outlined in red?"
%    - If yes  -> the blob detector saw the marker but the decoder
%                 couldn't read it (try shorter windows, polarity
%                 split, or sub-pixel corner refinement).
%    - If no   -> the marker simply isn't a quad-shaped blob in this
%                 frame (motion blur / occlusion / out-of-frame /
%                 too-small after resolution reduction).
%
%  Inputs:
%    resultFile : path to a .mat written by detectAruco
%    eventFile  : path to the .mat with `events` used to produce it
%    opts (optional struct):
%        .numTicks    (default 9)        ticks to show (NxN grid)
%        .windowMs    (default 8)        accumulation window
%        .mode        (default 'failed') 'failed' | 'success' | 'both'
%        .tickIndices (default [])       explicit tick indices to show
%        .sensorSize  (auto-inferred)    [H W]
%        .blobParams  (defaults)         struct .minArea/.maxArea/.maxAspect
%
%  Example:
%    diagnoseFailures( ...
%        'Data/OzanEventData_22.05.2026/3/3_reduced_results_v3.mat', ...
%        'Data/OzanEventData_22.05.2026/3/3_reduced.mat', ...
%        struct('numTicks', 9, 'windowMs', 8));

if nargin < 3, opts = struct(); end
if ~isfield(opts,'numTicks'),    opts.numTicks    = 9;        end
if ~isfield(opts,'windowMs'),    opts.windowMs    = 8;        end
if ~isfield(opts,'mode'),        opts.mode        = 'failed'; end
if ~isfield(opts,'tickIndices'), opts.tickIndices = [];       end
if ~isfield(opts,'sensorSize'),  opts.sensorSize  = [];       end
if ~isfield(opts,'blobParams')
    opts.blobParams = struct( ...
        'minArea',   100, ...
        'maxArea',   inf, ...
        'maxAspect', 3.0);
end

if ~isfile(resultFile), error('result file not found: %s', resultFile); end
if ~isfile(eventFile),  error('event file not found: %s',  eventFile);  end

fprintf('Loading result: %s\n', resultFile);
R = load(resultFile);
fprintf('Loading events: %s\n', eventFile);
E = load(eventFile);
events = E.events;
fprintf('  %d events\n', size(events, 1));

% Sensor size
if isempty(opts.sensorSize)
    H = double(max(events(:, 2))) + 1;
    W = double(max(events(:, 1))) + 1;
else
    H = opts.sensorSize(1); W = opts.sensorSize(2);
end
if isinf(opts.blobParams.maxArea)
    opts.blobParams.maxArea = 0.6 * H * W;
end
fprintf('  sensor: %d x %d  |  blob area in [%d, %d]\n', ...
    H, W, opts.blobParams.minArea, opts.blobParams.maxArea);

% Sort events by timestamp (binary search needs this)
[~, si] = sort(events(:, 4));
events  = events(si, :);
evT     = double(events(:, 4));
evX     = double(events(:, 1)) + 1;
evY     = double(events(:, 2)) + 1;

% Pick which ticks to show
anyDet = R.anyDetected > 0;
nDet   = sum(anyDet);
nMiss  = sum(~anyDet);
fprintf('  %d ticks total  |  %d detections  |  %d misses\n', ...
    length(anyDet), nDet, nMiss);

if ~isempty(opts.tickIndices)
    tickIndices = opts.tickIndices(:);
else
    switch lower(opts.mode)
        case 'failed'
            pool = find(~anyDet);
        case 'success'
            pool = find(anyDet);
        case 'both'
            half = ceil(opts.numTicks / 2);
            failPool    = find(~anyDet);
            successPool = find(anyDet);
            pool = [pickN(failPool, half); pickN(successPool, opts.numTicks - half)];
        otherwise
            error('Unknown mode: %s', opts.mode);
    end
    tickIndices = pickN(pool, opts.numTicks);
end

if isempty(tickIndices)
    error('No ticks matching mode=%s.', opts.mode);
end

nTicks = numel(tickIndices);
sz = ceil(sqrt(nTicks));
winUs = opts.windowMs * 1000;

% Window durations actually computed in the result, for the title line
if isfield(R, 'windowDurations_ms')
    runWins = double(R.windowDurations_ms);
else
    runWins = [];
end

[~, baseName] = fileparts(eventFile);
fig = figure('Name', sprintf('Diagnose: %s (%s ticks)', baseName, opts.mode), ...
             'Position', [60 60 1300 1100], 'Color', 'w');
tl = tiledlayout(fig, sz, sz, 'TileSpacing', 'compact', 'Padding', 'compact');

for ti = 1:nTicks
    tickIdx = tickIndices(ti);
    tNow    = double(R.tNow_us(tickIdx));
    tFrom   = tNow - winUs;

    iStart = binSearchLeft(evT,  tFrom);
    iEnd   = binSearchRight(evT, tNow);
    if iEnd >= iStart
        x = evX(iStart:iEnd);
        y = evY(iStart:iEnd);
        valid = (x >= 1) & (x <= W) & (y >= 1) & (y <= H);
        x = x(valid); y = y(valid);
    else
        x = []; y = [];
    end

    if isempty(x)
        countImg = zeros(H, W);
    else
        countImg = accumarray([y, x], 1, [H, W]);
    end
    activeMask = countImg > 0;

    quads = detectQuadBlob(activeMask, opts.blobParams);

    % decoded IDs at this tick across ALL run windows (from result file)
    decodedIDs = decodedIDsAtTick(R, tickIdx, runWins);

    nexttile(tl, ti);
    if max(countImg(:)) == 0
        imshow(uint8(zeros(H, W)));
    else
        imshow(uint8(countImg / max(countImg(:)) * 255), ...
               'InitialMagnification', 'fit');
    end
    hold on;
    for qi = 1:length(quads)
        q = quads{qi};
        plot([q(:,1); q(1,1)], [q(:,2); q(1,2)], '-', ...
             'Color', [1 0.25 0.25], 'LineWidth', 1.0);
    end
    hold off;

    if isempty(decodedIDs)
        status = 'no decode';
        color  = [0.85 0.1 0.1];
    else
        status = sprintf('IDs %s', mat2str(decodedIDs));
        color  = [0.1 0.55 0.1];
    end
    title(sprintf('tick %d  t=%.3fs\n%d quads  |  %s', ...
                  tickIdx, tNow / 1e6, length(quads), status), ...
          'FontSize', 8, 'Color', color);
end

title(tl, sprintf('Diagnose %s ticks  -  %d ms window  -  %s', ...
                  opts.mode, opts.windowMs, baseName), ...
      'FontSize', 11, 'FontWeight', 'bold');

end


%% =========================================================================
function idx = pickN(pool, n)
    if isempty(pool) || n <= 0
        idx = [];
        return;
    end
    if length(pool) <= n
        idx = pool(:);
        return;
    end
    pick = round(linspace(1, length(pool), n));
    idx = pool(pick);
end

function ids = decodedIDsAtTick(R, tickIdx, wins)
    ids = [];
    if isempty(wins), return; end
    for w = wins
        f = sprintf('win_%dms', w);
        if ~isfield(R, f), continue; end
        v = R.(f)(tickIdx);
        if v >= 0, ids = [ids, double(v)]; end %#ok<AGROW>
    end
    ids = unique(ids);
end

function idx = binSearchLeft(evT, tTarget)
    lo = 1; hi = length(evT);
    if hi == 0, idx = 1; return; end
    while lo < hi
        mid = floor((lo + hi) / 2);
        if evT(mid) < tTarget
            lo = mid + 1;
        else
            hi = mid;
        end
    end
    idx = lo;
end

function idx = binSearchRight(evT, tTarget)
    lo = 1; hi = length(evT);
    if hi == 0, idx = 0; return; end
    while lo < hi
        mid = ceil((lo + hi) / 2);
        if evT(mid) > tTarget
            hi = mid - 1;
        else
            lo = mid;
        end
    end
    idx = lo;
end
