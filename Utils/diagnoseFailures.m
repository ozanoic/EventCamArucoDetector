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
%        .failingId   (default [])       show ticks where this marker
%                                        was MISSED. Overrides .mode.
%                                        Requires anyDetected_id<N> in
%                                        the result file (re-run
%                                        detectAruco with this ID in
%                                        requestedMarkerIds).
%        .referenceId (default [])       used with .failingId. When
%                                        set, only show ticks where
%                                        the reference marker WAS
%                                        detected (= "we know the
%                                        camera could see something
%                                        useful here"). Empty -> any
%                                        detection counts as reference.
%        .sensorSize  (auto-inferred)    [H W]
%        .blobParams  (defaults)         struct .minArea/.maxArea/.maxAspect
%
%  Playback (video-style) mode:
%        .playback     (default false)   when true, play ticks as a
%                                        movie instead of drawing a
%                                        static grid. Quads are drawn
%                                        GREEN when the result file
%                                        says any marker was decoded
%                                        at this tick (for the chosen
%                                        window), RED otherwise.
%        .tickStart    (default 1)       first tick to play
%        .tickEnd      (default Inf)     last tick to play
%        .stepMs       (default 50)      wall-clock pause between
%                                        frames (ms). Lower = faster
%                                        playback.
%        .saveVideo    (default false)   also write the playback as
%                                        an MP4 alongside the result.
%        .videoPath    (default '')      explicit output path.
%
%        While playing, the figure listens for:
%            space  -> pause / resume
%            q      -> quit
%            -> / left arrow  -> step (only useful while paused)
%
%  Examples:
%    % All failed ticks
%    diagnoseFailures(resultFile, eventFile);
%
%    % Why did marker 3 fail when something else was visible?
%    diagnoseFailures(resultFile, eventFile, struct('failingId', 3));
%
%    % Same, but only when marker 8 specifically was found
%    diagnoseFailures(resultFile, eventFile, ...
%        struct('failingId', 3, 'referenceId', 8));
%
%    % Play the whole recording at 10 ms, green = decoded, red = miss
%    diagnoseFailures(resultFile, eventFile, ...
%        struct('playback', true, 'windowMs', 10, 'stepMs', 60));
%
%    % Play and also save to MP4
%    diagnoseFailures(resultFile, eventFile, ...
%        struct('playback', true, 'saveVideo', true));

if nargin < 3, opts = struct(); end
if ~isfield(opts,'numTicks'),    opts.numTicks    = 9;        end
if ~isfield(opts,'windowMs'),    opts.windowMs    = 8;        end
if ~isfield(opts,'mode'),        opts.mode        = 'failed'; end
if ~isfield(opts,'tickIndices'), opts.tickIndices = [];       end
if ~isfield(opts,'failingId'),   opts.failingId   = [];       end
if ~isfield(opts,'referenceId'), opts.referenceId = [];       end
if ~isfield(opts,'sensorSize'),  opts.sensorSize  = [];       end
if ~isfield(opts,'playback'),    opts.playback    = false;    end
if ~isfield(opts,'tickStart'),   opts.tickStart   = 1;        end
if ~isfield(opts,'tickEnd'),     opts.tickEnd     = Inf;      end
if ~isfield(opts,'stepMs'),      opts.stepMs      = 50;       end
if ~isfield(opts,'saveVideo'),   opts.saveVideo   = false;    end
if ~isfield(opts,'videoPath'),   opts.videoPath   = '';       end
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

% --- Playback branch: bail out of the grid path entirely ---
if opts.playback
    playbackImpl(R, opts, H, W, evT, evX, evY, eventFile, resultFile, anyDet);
    return;
end

selectorLabel = opts.mode;
if ~isempty(opts.tickIndices)
    tickIndices = opts.tickIndices(:);
    selectorLabel = sprintf('explicit %d ticks', numel(tickIndices));
elseif ~isempty(opts.failingId)
    % Per-marker miss diagnosis. Need the per-marker any-detect column
    % the detector writes when requestedMarkerIds is non-empty.
    fId    = opts.failingId;
    fField = sprintf('anyDetected_id%d', fId);
    if ~isfield(R, fField)
        error('diagnoseFailures:missingField', ...
              'result file has no %s field. Re-run detectAruco with %d in params.requestedMarkerIds.', ...
              fField, fId);
    end
    missed = R.(fField) == 0;

    if ~isempty(opts.referenceId)
        rId    = opts.referenceId;
        rField = sprintf('anyDetected_id%d', rId);
        if ~isfield(R, rField)
            error('diagnoseFailures:missingField', ...
                  'result file has no %s field. Re-run detectAruco with %d in params.requestedMarkerIds.', ...
                  rField, rId);
        end
        reference     = R.(rField) > 0;
        selectorLabel = sprintf('id %d missed, id %d found', fId, rId);
    else
        reference     = anyDet;
        selectorLabel = sprintf('id %d missed (anything else found)', fId);
    end

    pool = find(missed & reference);
    fprintf('  %d ticks where %s\n', length(pool), selectorLabel);
    if isempty(pool)
        error(['No ticks matching "%s". Either the marker is always ' ...
               'detected when something else is, or the reference ' ...
               'marker never coincides with this miss.'], selectorLabel);
    end
    tickIndices = pickN(pool, opts.numTicks);
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
    error('No ticks matching selector=%s.', selectorLabel);
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
fig = figure('Name', sprintf('Diagnose: %s (%s)', baseName, selectorLabel), ...
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

title(tl, sprintf('Diagnose  -  %s  -  %d ms window  -  %s', ...
                  selectorLabel, opts.windowMs, baseName), ...
      'FontSize', 11, 'FontWeight', 'bold');

end


%% =========================================================================
%                              PLAYBACK MODE
%% =========================================================================
function playbackImpl(R, opts, H, W, evT, evX, evY, eventFile, resultFile, anyDet)
    nTicks    = length(R.tNow_us);
    tickStart = max(1, opts.tickStart);
    tickEnd   = min(nTicks, opts.tickEnd);
    winUs     = opts.windowMs * 1000;

    % Decide which ticks to actually iterate over based on mode.
    switch lower(opts.mode)
        case 'failed'
            allowed = ~anyDet;
            modeStr = 'failed only';
        case 'success'
            allowed = anyDet;
            modeStr = 'successful only';
        case {'both', 'all'}
            allowed = true(size(anyDet));
            modeStr = 'all ticks';
        otherwise
            error('Unknown mode for playback: %s', opts.mode);
    end
    playableTicks = find(allowed);
    playableTicks = playableTicks(playableTicks >= tickStart & playableTicks <= tickEnd);
    if isempty(playableTicks)
        error('No ticks in range [%d, %d] matching mode=%s.', tickStart, tickEnd, opts.mode);
    end
    fprintf('Playback: %d ticks (%s), %d ms window\n', ...
        length(playableTicks), modeStr, opts.windowMs);

    % Pre-compute, for each playable tick, whether the chosen window
    % decoded any marker -- this is what flips quads green vs red.
    runWins = double(R.windowDurations_ms);
    closestWin = pickClosestWindow(runWins, opts.windowMs);
    if closestWin ~= opts.windowMs
        fprintf('Note: result file has no %d ms window, using nearest %d ms for green/red marking.\n', ...
            opts.windowMs, closestWin);
    end
    winField = sprintf('win_%dms', closestWin);
    winColRaw = R.(winField);

    % --- Optional video writer ---
    vw = [];
    if opts.saveVideo
        outPath = opts.videoPath;
        if isempty(outPath)
            [d, n] = fileparts(resultFile);
            outPath = fullfile(d, [n '_playback.mp4']);
        end
        try
            vw = VideoWriter(outPath, 'MPEG-4');
        catch
            vw = VideoWriter(outPath);                 % AVI fallback
        end
        vw.FrameRate = max(1, round(1000 / opts.stepMs));
        open(vw);
        fprintf('Recording video to %s @ %d fps\n', outPath, vw.FrameRate);
    end
    cleanupVid = onCleanup(@() closeIfOpen(vw));

    % --- Figure / axes / state ---
    [~, baseName] = fileparts(eventFile);
    state = struct('paused', false, 'quit', false, 'step', false);
    fig = figure('Name', sprintf('Playback: %s', baseName), ...
                 'Color', [0.08 0.08 0.08], ...
                 'Position', [60 60 1100 850]);
    set(fig, 'KeyPressFcn', @(s, e) onKey(s, e));
    ax = axes(fig, 'Position', [0.04 0.06 0.93 0.88], 'Color', 'k');

    function onKey(~, e)
        switch e.Key
            case 'space',  state.paused = ~state.paused;
            case 'q',      state.quit = true; state.paused = false;
            case {'rightarrow', 'leftarrow', 'd'},  state.step = true;
        end
    end

    fprintf('Controls: SPACE = pause/resume, q = quit, right-arrow = step (while paused)\n');

    for k = 1:length(playableTicks)
        if state.quit || ~isvalid(fig), break; end
        tick = playableTicks(k);

        % --- accumulate events ---
        tNow  = double(R.tNow_us(tick));
        tFrom = tNow - winUs;
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

        % --- candidate quads ---
        quads = detectQuadBlob(countImg > 0, opts.blobParams);

        % --- did this window decode anything at this tick?
        % Quads turn green when yes, red otherwise.
        windowDecoded = winColRaw(tick) >= 0;
        decodedIDs    = decodedIDsAtTick(R, tick, runWins);
        if windowDecoded
            quadColor = [0.2 1.0 0.2];
        else
            quadColor = [1.0 0.3 0.3];
        end

        % --- render ---
        if ~isvalid(fig), break; end
        cla(ax);
        if max(countImg(:)) == 0
            imshow(uint8(zeros(H, W)), 'Parent', ax);
        else
            imshow(uint8(countImg / max(countImg(:)) * 255), 'Parent', ax);
        end
        hold(ax, 'on');
        for qi = 1:length(quads)
            q = quads{qi};
            plot(ax, [q(:,1); q(1,1)], [q(:,2); q(1,2)], '-', ...
                 'Color', quadColor, 'LineWidth', 1.6);
        end
        hold(ax, 'off');

        if isempty(decodedIDs)
            statusStr = 'no decode';
        else
            statusStr = sprintf('decoded IDs: %s', mat2str(decodedIDs));
        end
        title(ax, sprintf(['tick %d / %d   t = %.3f s   %d ms window   ' ...
                          '%d quads   %s   (this window: %s)'], ...
            tick, nTicks, tNow / 1e6, opts.windowMs, length(quads), ...
            statusStr, tern(windowDecoded, 'HIT', 'miss')), ...
            'Color', quadColor, 'FontSize', 11, 'FontWeight', 'bold');
        drawnow;

        % --- video frame ---
        if ~isempty(vw)
            try
                writeVideo(vw, getframe(fig));
            catch
                % closed window or write failure; stop recording.
                close(vw);
                vw = [];
            end
        end

        % --- pacing / interactivity ---
        if state.paused
            while state.paused && ~state.quit && isvalid(fig)
                pause(0.05);
                if state.step, state.step = false; break; end
            end
        else
            pause(max(opts.stepMs, 1) / 1000);
        end
    end

    if ~isempty(vw)
        close(vw);
        fprintf('Video saved.\n');
    end
end

function closeIfOpen(vw)
    if ~isempty(vw)
        try, close(vw); catch, end %#ok<NOCOM,CTCH>
    end
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end

function w = pickClosestWindow(wins, target)
    if isempty(wins), w = target; return; end
    [~, idx] = min(abs(double(wins) - double(target)));
    w = double(wins(idx));
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
