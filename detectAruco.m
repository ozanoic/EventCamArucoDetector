function results = detectAruco(matFile, sensorSize, params)
%DETECTARUCO Multi-window sliding ArUco detection for event camera data.
%
%  results = detectAruco(matFile, sensorSize, params)
%
%  Inputs:
%    matFile    - path to .mat file containing 'events' (Nx4: x,y,pol,t)
%    sensorSize - [height, width]
%    params     - struct with fields:
%      .windowDurations_ms   - vector of lookback windows in ms (e.g. [5 10 20 50 100])
%      .tickStep_us          - tick step in microseconds (default 1000)
%      .showVis              - enable per-tick visualization (default false)
%      .useParallel          - true: use parfor if toolbox available (default);
%                              false: force sequential execution
%      .numCells             - marker grid size (default 8)
%      .codeSize             - inner code grid size (default 6)
%      .cellPx               - pixels per cell in unwarped image (default 20)
%      .blobParams           - struct with .minArea, .maxArea, .maxAspect
%      .requestedMarkerIds   - scalar or vector of ArUco IDs to accept
%                              (default [] = accept any decoded marker).
%                              Decodes whose ID is NOT in the set are
%                              treated as no-detection for that tick/window.
%      .earlyExitOnFirstHit  - true: stop a tick's window loop after the
%                              first successful decode (faster, but the
%                              per-window detection breakdown is no longer
%                              meaningful). Default false.
%      .detectScale          - scalar in (0, 1] (default 1). When < 1, blob
%                              detection runs on a downscaled mask; quad
%                              corners are scaled back up before the warp.
%                              0.5 is roughly 4x faster on the blob stage
%                              at the cost of ~1 px corner accuracy.
%      .useGPU               - true: push the per-quad imwarp to GPU via
%                              gpuArray (forces sequential mode because
%                              the GPU cannot be shared across parfor
%                              workers efficiently). Default false.
%      .decoderHammingDist   - >=0, default 0. Accept decodes within k
%                              bit-flips of any dict code (Hamming
%                              tolerance for noisy DVS). 0 = exact match.
%      .minEventsPerWindow   - skip windows whose event count is below
%                              this. Default 10. Lower it for faint or
%                              short-duration markers.
%      .refineCorners        - true: sub-pixel corner refinement using
%                              perpendicular gradient peaks (per-edge
%                              line fit + intersection). Default false.
%                              Big win when blob-detector corners are
%                              imprecise (typical for real DVS).
%      .refineSearchPx       - perpendicular search half-width (px).
%                              Default 5.
%
%  Output:
%    results - struct with fields:
%      .tNow_us             - timestamps (us) for each tick
%      .anyDetected         - 1 if any window detected a marker at that tick
%      .win_Xms             - detected marker ID per window (-1 = none / filtered out)
%      .windowDurations_ms  - copy of window durations used
%      .detectionsPerWindow - total detections per window
%      .requestedMarkerIds  - echo of the filter that was applied (double row)

%% ---- Default parameters ----
if ~isfield(params,'tickStep_us'),  params.tickStep_us  = 1000;  end
if ~isfield(params,'showVis'),      params.showVis      = false; end
if ~isfield(params,'useParallel'),  params.useParallel  = true;  end
if ~isfield(params,'numCells'),     params.numCells     = 8;     end
if ~isfield(params,'codeSize'),     params.codeSize     = 6;     end
if ~isfield(params,'cellPx'),       params.cellPx       = 20;    end
if ~isfield(params,'requestedMarkerIds'),   params.requestedMarkerIds   = []; end
if ~isfield(params,'earlyExitOnFirstHit'),  params.earlyExitOnFirstHit  = false; end
if ~isfield(params,'detectScale'),          params.detectScale          = 1.0;  end
if ~isfield(params,'useGPU'),               params.useGPU               = false; end
if ~isfield(params,'decoderHammingDist'),   params.decoderHammingDist   = 0;    end
if ~isfield(params,'minEventsPerWindow'),   params.minEventsPerWindow   = 10;   end
if ~isfield(params,'refineCorners'),        params.refineCorners        = false; end
if ~isfield(params,'refineSearchPx'),       params.refineSearchPx       = 5;     end

% Normalise requested marker IDs into a sorted int32 row vector.
requestedMarkerIds = int32(params.requestedMarkerIds(:)');
if isempty(requestedMarkerIds)
    fprintf('Requested markers: ANY (no filter)\n');
else
    fprintf('Requested markers: %s\n', mat2str(requestedMarkerIds));
end

windowDurations_ms = params.windowDurations_ms;
windowDurations_us = windowDurations_ms * 1000;
numWindows = length(windowDurations_us);
tickStep_us = params.tickStep_us;

numCells = params.numCells;
codeSize = params.codeSize;
cellPx   = params.cellPx;
sideSize = numCells * cellPx;
markerCoords = [0 0; sideSize-1 0; sideSize-1 sideSize-1; 0 sideSize-1];

H = sensorSize(1);
W = sensorSize(2);

%% ---- Load data ----
fprintf('Loading %s...\n', matFile);
tmp = load(matFile, 'events');
events = tmp.events;
numEvents = size(events, 1);
fprintf('Loaded %d events  |  sensor %dx%d\n', numEvents, W, H);

% Build dictionary (sorted arrays, parfor-compatible)
[dictCodes, dictIDs] = buildDictionaryArrays();
fprintf('Dictionary: %d markers\n', length(dictCodes));

%% ---- Precompute event data (sorted by time, compact types) ----
% Use the smallest types that still represent the data losslessly. This
% keeps the per-worker broadcast bundle small enough to deserialize on
% large recordings (a 480x640 sensor + long capture can easily exceed
% RAM if every coordinate is stored as a double).
evX = uint16(events(:,1) + 1);          % 1..W           (max 65535)
evY = uint16(events(:,2) + 1);          % 1..H
evT = int64(events(:,4));               % microseconds  (>= int64 range)

[evT, sortIdx] = sort(evT);
evX = evX(sortIdx);
evY = evY(sortIdx);
clear events tmp sortIdx;               % free the Nx4 double copy

fprintf('Event arrays:  evT %s (%.1f MB), evX/evY %s (%.1f MB each)\n', ...
    class(evT), numel(evT)*8/1e6, class(evX), numel(evX)*2/1e6);

tMin = double(evT(1));
tMax = double(evT(end));
tStart = tMin + max(windowDurations_us);
tEnd   = tMax;
numTicks = floor((tEnd - tStart) / tickStep_us) + 1;

fprintf('Time range: %.3fs to %.3fs\n', tMin/1e6, tMax/1e6);
fprintf('Ticks: %d  (every %d us = %.0f ms)\n', numTicks, tickStep_us, tickStep_us/1000);
fprintf('Windows: %d  (%s ms)\n', numWindows, mat2str(windowDurations_ms));

%% ---- Decide execution mode ----
if params.useParallel
    hasParallel = ~isempty(ver('parallel'));
    if hasParallel
        fprintf('useParallel=TRUE  |  Parallel Computing Toolbox: FOUND\n');
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool('local');
        end
        fprintf('Using %d workers\n', pool.NumWorkers);
    else
        fprintf('useParallel=TRUE  |  Parallel Computing Toolbox: NOT FOUND (falling back to sequential)\n');
    end
else
    hasParallel = false;
    fprintf('useParallel=FALSE |  Running sequential (parallel disabled by user)\n');
end

%% ---- Preallocate output ----
resultTable = zeros(numTicks, 2 + numWindows);
detectionsPerWindow = zeros(1, numWindows);

% Per-marker per-window tracking.  When requestedMarkerIds is non-empty
% we run multi-marker mode: each window scans EVERY quad and records,
% for each requested ID, whether that ID was decoded by some quad in
% this (tick, window). The result is a numTicks x numWindows x nMarkers
% logical tensor.  Stored flat (numTicks x (numWindows*nMarkers)) so
% parfor can write a row per iteration.
reportedIdsLocal = int32(requestedMarkerIds);
nReportedLocal   = numel(reportedIdsLocal);
if nReportedLocal > 0
    perMarkerFlat = false(numTicks, numWindows * nReportedLocal);
else
    perMarkerFlat = [];   % single-marker (legacy) mode
end

blobParams        = params.blobParams;
showVis           = params.showVis;
earlyExitOnHit    = logical(params.earlyExitOnFirstHit);
detectScale       = params.detectScale;
useGPU            = logical(params.useGPU);
hammingThresh     = max(0, round(params.decoderHammingDist));
minEventsPerWin   = max(1, round(params.minEventsPerWindow));
refineCorners     = logical(params.refineCorners);
refineSearchPx    = max(2, round(params.refineSearchPx));
if refineCorners
    fprintf('Corner refinement: ON  (perpendicular search half-width = %d px)\n', refineSearchPx);
end
if hammingThresh > 0
    fprintf('Decoder Hamming threshold: %d  (accept decodes within %d bit-flips of any dict code)\n', ...
        hammingThresh, hammingThresh);
end
if minEventsPerWin ~= 10
    fprintf('minEventsPerWindow = %d  (skip windows with fewer events)\n', minEventsPerWin);
end

if detectScale <= 0 || detectScale > 1
    error('detectAruco: params.detectScale must be in (0, 1]; got %g.', detectScale);
end
if useGPU && hasParallel
    fprintf('useGPU=TRUE  but useParallel=TRUE -> forcing sequential (GPU is shared per process).\n');
    hasParallel = false;
end
if useGPU
    try
        gd = gpuDevice;
        fprintf('GPU enabled: %s (CUDA %d.%d, %.1f GB)\n', gd.Name, ...
            gd.ComputeCapabilityMajor, gd.ComputeCapabilityMinor, ...
            gd.AvailableMemory / 1e9);
    catch
        fprintf('useGPU=TRUE but no usable GPU found -> falling back to CPU.\n');
        useGPU = false;
    end
end
if detectScale < 1
    fprintf('detectScale = %.2f  ->  blob stage at %dx%d, corners refined to %dx%d\n', ...
        detectScale, round(W*detectScale), round(H*detectScale), W, H);
end
if earlyExitOnHit
    fprintf('earlyExitOnFirstHit = TRUE  ->  per-tick window loop stops at first decode\n');
end

if showVis
    hFig = figure('Name', 'Multi-Window Detection', 'Position', [30 30 1800 900]);
end

%% ---- Process ----
fprintf('\nProcessing %d ticks x %d windows = %d attempts...\n', ...
    numTicks, numWindows, numTicks * numWindows);
tic;

if hasParallel
    % ==================================================================
    %  PARALLEL PATH
    % ==================================================================
    tickTimes = tStart + (0:(numTicks-1))' * tickStep_us;

    parResultTable = zeros(numTicks, 2 + numWindows);
    parDetPerWindow = zeros(numTicks, numWindows);

    dq = parallel.pool.DataQueue;
    parTicStart = tic;
    progressCount = containers.Map('count', 0);
    lastPct = containers.Map('pct', 0);
    afterEach(dq, @(~) parProgressCallback(progressCount, lastPct, numTicks, parTicStart));

    % --- Send the big event arrays to each worker as broadcast vars.
    % We used to wrap these in parallel.pool.Constant to avoid one copy
    % per iteration, but that proved fragile (Constant serialisation can
    % fail on OneDrive-locked temp dirs and on a stale parpool).  Since
    % we now store events as int64/uint16 instead of doubles, the
    % broadcast bundle is ~6x smaller than it used to be and fits
    % comfortably for typical recordings.
    %
    % If you hit "Out of Memory during deserialization" here:
    %   * reduce the parpool worker count: delete(gcp('nocreate'));
    %     parpool('local', 4)
    %   * or set params.useParallel = false to process sequentially.
    bytes = numel(evT)*8 + numel(evX)*2 + numel(evY)*2;
    fprintf('Broadcast bundle per worker: %.1f MB\n', bytes / 1e6);

    parPerMarkerFlat = false(numTicks, max(numWindows * nReportedLocal, 1));
    parStats = zeros(numTicks, 6);   % rolled up to a global funnel after parfor

    parfor tick = 1:numTicks
        tNow = tickTimes(tick);
        tickRow = -1 * ones(1, numWindows);
        tickPerMarker = false(numWindows, max(nReportedLocal, 1));
        tickStats = zeros(1, 6);

        for wi = 1:numWindows
            [bestID, detVec, statsRow] = processWindow_local( ...
                tNow, windowDurations_us(wi), ...
                evT, evX, evY, H, W, ...
                blobParams, detectScale, useGPU, ...
                markerCoords, sideSize, ...
                numCells, codeSize, cellPx, ...
                dictCodes, dictIDs, ...
                requestedMarkerIds, reportedIdsLocal, ...
                hammingThresh, minEventsPerWin, ...
                refineCorners, refineSearchPx);

            tickRow(wi) = bestID;
            if nReportedLocal > 0
                tickPerMarker(wi, :) = detVec;
            end
            tickStats = tickStats + statsRow;

            if earlyExitOnHit
                if nReportedLocal > 0
                    % stop only when EVERY reported marker has been seen
                    if all(any(tickPerMarker, 1)), break; end
                else
                    if bestID >= 0, break; end
                end
            end
        end

        anyDet = double(any(tickRow >= 0));
        parResultTable(tick, :) = [tNow, anyDet, tickRow];
        parDetPerWindow(tick, :) = double(tickRow >= 0);
        if nReportedLocal > 0
            parPerMarkerFlat(tick, :) = reshape(tickPerMarker, 1, []);
        end
        parStats(tick, :) = tickStats;

        send(dq, tick);
    end

    resultTable = parResultTable;
    detectionsPerWindow = sum(parDetPerWindow, 1);
    if nReportedLocal > 0
        perMarkerFlat = parPerMarkerFlat;
    end
    funnelStats = sum(parStats, 1);

else
    % ==================================================================
    %  SEQUENTIAL PATH
    % ==================================================================
    seqTicStart = tic;
    funnelStats = zeros(1, 6);
    for tick = 1:numTicks
        tNow = tStart + (tick - 1) * tickStep_us;
        tickRow = -1 * ones(1, numWindows);
        tickPerMarker = false(numWindows, max(nReportedLocal, 1));

        for wi = 1:numWindows
            [bestID, detVec, statsRow] = processWindow_local( ...
                tNow, windowDurations_us(wi), ...
                evT, evX, evY, H, W, ...
                blobParams, detectScale, useGPU, ...
                markerCoords, sideSize, ...
                numCells, codeSize, cellPx, ...
                dictCodes, dictIDs, requestedMarkerIds, reportedIdsLocal, ...
                hammingThresh, minEventsPerWin);

            tickRow(wi) = bestID;
            if nReportedLocal > 0
                tickPerMarker(wi, :) = detVec;
            end
            if bestID >= 0
                detectionsPerWindow(wi) = detectionsPerWindow(wi) + 1;
            end
            funnelStats = funnelStats + statsRow;

            if earlyExitOnHit
                if nReportedLocal > 0
                    if all(any(tickPerMarker, 1)), break; end
                else
                    if bestID >= 0, break; end
                end
            end
        end

        anyDet = double(any(tickRow >= 0));
        resultTable(tick, :) = [tNow, anyDet, tickRow];
        if nReportedLocal > 0
            perMarkerFlat(tick, :) = reshape(tickPerMarker, 1, []);
        end

        % Progress output (every 1%)
        pctNow = floor(100 * tick / numTicks);
        pctPrev = floor(100 * (tick-1) / numTicks);
        if pctNow > pctPrev
            elapsedSec = toc(seqTicStart);
            etaSec = elapsedSec / tick * (numTicks - tick);
            fprintf('  %3d%% processed  (%d/%d ticks)  elapsed: %.0fs  ETA: %.0fs\n', ...
                pctNow, tick, numTicks, elapsedSec, etaSec);
        end

        % Visualization
        if showVis
            figure(hFig); clf;
            subplot(2,1,1);
            colors = zeros(numWindows, 3);
            for wi = 1:numWindows
                if tickRow(wi) >= 0
                    colors(wi,:) = [0 0.8 0];
                else
                    colors(wi,:) = [0.8 0 0];
                end
            end
            b = bar(1:numWindows, ones(1,numWindows), 'FaceColor', 'flat');
            b.CData = colors;
            set(gca, 'XTick', 1:numWindows, 'XTickLabel', ...
                arrayfun(@(x) sprintf('%dms', x), windowDurations_ms, 'UniformOutput', false));
            xlabel('Window duration');
            title(sprintf('Tick %d/%d  |  t = %.4f s', tick, numTicks, tNow/1e6));
            ylabel('Detection'); ylim([0 1.5]);

            refWin = min(6, numWindows);
            dt = windowDurations_us(refWin);
            if isa(evT, 'int64')
                iS = bsearchLeft(evT, int64(tNow - dt));
                iE = bsearchRight(evT, int64(tNow));
            else
                iS = bsearchLeft(evT, tNow - dt);
                iE = bsearchRight(evT, tNow);
            end
            if iE >= iS
                rX = double(evX(iS:iE)); rY = double(evY(iS:iE));
                vm = (rY >= 1) & (rY <= H) & (rX >= 1) & (rX <= W);
                if any(vm)
                    refImg = accumarray([rY(vm), rX(vm)], 1, [H, W]);
                    subplot(2,1,2);
                    imshow(uint8(refImg / max(refImg(:)) * 255));
                    title(sprintf('Event image (%dms window)', windowDurations_ms(refWin)));
                end
            end
            drawnow;
        end
    end
end

elapsed = toc;

%% ---- Performance Summary ----
totalDetections = sum(detectionsPerWindow);
ticksWithAnyDetection = sum(resultTable(:, 2));

fprintf('\n========================================\n');
fprintf('       PERFORMANCE SUMMARY\n');
fprintf('========================================\n');
fprintf('Execution mode:          %s\n', ternary(hasParallel, 'PARALLEL', 'SEQUENTIAL'));
fprintf('Total ticks:             %d\n', numTicks);
fprintf('Ticks with detection:    %d (%.1f%%)\n', ...
    ticksWithAnyDetection, 100*ticksWithAnyDetection/max(numTicks,1));
fprintf('Total detections:        %d (across all windows)\n', totalDetections);
fprintf('Elapsed time:            %.1fs\n', elapsed);
% --- Pipeline funnel diagnostic ---
% Tells you where you're losing detections. The big drop between any
% two consecutive rows is the stage to look at.
fprintf('\n--- Pipeline funnel (sum over every tick x window attempt) ---\n');
fnNames = {'window attempts', ...
           'with >= minEvents events', ...
           'quads found (sum)', ...
           'decode attempts', ...
           'decode success (any ID)', ...
           'decode in requestedIds'};
fprintf('%-30s  %10s  %s\n', 'stage', 'count', 'survival %');
fprintf('%-30s  %10s  %s\n', '-----', '-----', '----------');
totalAttempts = max(funnelStats(1), 1);
for fi = 1:length(fnNames)
    fprintf('%-30s  %10d  %9.2f%%\n', fnNames{fi}, funnelStats(fi), ...
        100 * funnelStats(fi) / totalAttempts);
end

fprintf('\n--- Per-window breakdown (all markers combined) ---\n');
fprintf('%-10s  %10s  %10s\n', 'Window', 'Detections', 'Rate');
fprintf('%-10s  %10s  %10s\n', '------', '----------', '----');
for wi = 1:numWindows
    fprintf('%-10s  %10d  %9.1f%%\n', ...
        sprintf('%dms', windowDurations_ms(wi)), ...
        detectionsPerWindow(wi), ...
        100 * detectionsPerWindow(wi) / max(numTicks, 1));
end

%% ---- Per-marker breakdown ----
% In multi-marker mode (requestedMarkerIds non-empty) we have a
% full numTicks x numWindows x nMarkers logical tensor (perMarkerFlat),
% so we can answer "did marker 8 appear at tick T in window W?" even
% when another marker was decoded by an earlier quad in the same frame.
%
% In legacy single-marker mode we fall back to the scalar win_Xms
% columns and only see the FIRST decoded ID per (tick, window).
winColsAll = resultTable(:, 3:end);    % numTicks x numWindows scalar IDs
if ~isempty(requestedMarkerIds)
    reportIds = double(requestedMarkerIds);
else
    reportIds = unique(winColsAll(winColsAll >= 0))';
end

perMarker = struct();
useMultiMarker = ~isempty(perMarkerFlat) && ~isempty(reportIds);
if useMultiMarker
    perMarkerTensor = reshape(perMarkerFlat, numTicks, numWindows, length(reportIds));
end

if ~isempty(reportIds)
    fprintf('\n--- Per-marker breakdown ---\n');
    fprintf('%-8s  %8s  %12s  %8s\n', 'Marker', 'AnyHit', 'AnyHit %', 'TotalHit');
    fprintf('%-8s  %8s  %12s  %8s\n', '------', '------', '--------', '--------');
    for ri = 1:length(reportIds)
        mid = reportIds(ri);
        if useMultiMarker
            hits = perMarkerTensor(:, :, ri);     % numTicks x numWindows logical
        else
            hits = winColsAll == double(mid);
        end
        anyForId = any(hits, 2);
        countPerWin = sum(hits, 1);
        nAny  = sum(anyForId);
        nHit  = sum(countPerWin);
        fprintf('id=%-4d   %8d   %11.1f%%  %8d\n', ...
            mid, nAny, 100*nAny/max(numTicks,1), nHit);
        perMarker(ri).id                  = mid;
        perMarker(ri).anyDetected         = double(anyForId);
        perMarker(ri).detectionsPerWindow = countPerWin;
        perMarker(ri).hits                = hits;   % carried into results below
    end

    % Per-window table showing how each marker fares per window
    fprintf('\n--- Per-window, per-marker detection counts ---\n');
    hdr = sprintf('%-10s', 'Window');
    for ri = 1:length(reportIds)
        hdr = [hdr sprintf('  %10s', sprintf('id=%d', reportIds(ri)))]; %#ok<AGROW>
    end
    fprintf('%s\n', hdr);
    fprintf('%s\n', repmat('-', 1, length(hdr)));
    for wi = 1:numWindows
        row = sprintf('%-10s', sprintf('%dms', windowDurations_ms(wi)));
        for ri = 1:length(reportIds)
            row = [row sprintf('  %10d', perMarker(ri).detectionsPerWindow(wi))]; %#ok<AGROW>
        end
        fprintf('%s\n', row);
    end
end
fprintf('========================================\n');

%% ---- Build output struct ----
results.tNow_us = resultTable(:, 1);
results.anyDetected = resultTable(:, 2);
for wi = 1:numWindows
    varName = sprintf('win_%dms', windowDurations_ms(wi));
    results.(varName) = resultTable(:, 2 + wi);
end
results.windowDurations_ms = windowDurations_ms;
results.detectionsPerWindow = detectionsPerWindow;
results.requestedMarkerIds  = double(requestedMarkerIds);

% Per-marker fields. For each reported ID we write:
%   anyDetected_id<N>           : Nx1 double (0/1)  -- any window saw N
%   detectionsPerWindow_id<N>   : 1xW double         -- per-window counts
%   win_<X>ms_id<N>             : Nx1 double (0/1)  -- per (tick, window)
%                                                       per-marker hits.
% The win_<X>ms_id<N> columns are the per-marker analogue of the legacy
% scalar win_<X>ms columns and let downstream code answer "in this
% window at this tick, did marker N fire?" -- even if the same tick had
% several markers visible. They are only available in multi-marker mode
% (requestedMarkerIds non-empty).
results.markerIdsReported = reportIds;
for ri = 1:length(reportIds)
    mid = reportIds(ri);
    fname_any = sprintf('anyDetected_id%d', mid);
    fname_cnt = sprintf('detectionsPerWindow_id%d', mid);
    results.(fname_any) = perMarker(ri).anyDetected;
    results.(fname_cnt) = perMarker(ri).detectionsPerWindow;

    if useMultiMarker
        hits = perMarker(ri).hits;       % numTicks x numWindows logical
        for wi = 1:numWindows
            fname = sprintf('win_%dms_id%d', windowDurations_ms(wi), mid);
            results.(fname) = double(hits(:, wi));
        end
    end
end

end


%% =========================================================================
%                       LOCAL FUNCTIONS
%% =========================================================================

function parProgressCallback(progressCount, lastPct, numTicks, parTicStart)
    progressCount('count') = progressCount('count') + 1;
    done = progressCount('count');
    pctNow = floor(100 * done / numTicks);
    if pctNow > lastPct('pct')
        lastPct('pct') = pctNow;
        elapsedSec = toc(parTicStart);
        etaSec = elapsedSec / done * (numTicks - done);
        fprintf('  %3d%% processed  (%d/%d ticks)  elapsed: %.0fs  ETA: %.0fs\n', ...
            pctNow, done, numTicks, elapsedSec, etaSec);
    end
end

function idx = bsearchLeft(evT, tTarget)
    lo = 1; hi = length(evT);
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

function idx = bsearchRight(evT, tTarget)
    lo = 1; hi = length(evT);
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

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function [bestID, detVec, stats] = processWindow_local( ...
        tNow, dt, evT, evX, evY, H, W, ...
        blobParams, detectScale, useGPU, ...
        markerCoords, sideSize, ...
        numCells, codeSize, cellPx, ...
        dictCodes, dictIDs, requestedMarkerIds, reportedIds, ...
        hammingThresh, minEventsPerWindow, ...
        refineCorners, refineSearchPx)
    % stats : 1 x 6 double row recording the pipeline funnel for this
    %         single (tick, window) attempt. Columns are
    %           [ windowAttempted, eventsSliced, quadsFound,
    %             decodeAttempts, decodeSuccess, decodeInFilter ].
    %         The caller sums these to a global funnel report.
    stats = zeros(1, 6);
    stats(1) = 1;   % windowAttempted
    % Full per-window pipeline: slice events -> count image -> blob detect
    % -> per-quad unwarp -> decode.
    %
    % Returns:
    %   bestID  : the FIRST decoded marker ID that matched the filter
    %             (-1 if none).  Kept for backward compatibility with
    %             the existing win_<X>ms columns.
    %   detVec  : 1 x length(reportedIds) logical.  detVec(k) is true
    %             if ANY quad in this (tick, window) decoded to
    %             reportedIds(k). Empty when reportedIds is empty
    %             (legacy single-marker mode).
    %
    % We scan EVERY surviving quad so that a frame containing several
    % markers (e.g. ids 3 and 8 visible at once) is recorded as such.
    % Only break early if all reported markers have been seen.

    bestID = -1;
    nReported = numel(reportedIds);
    detVec = false(1, nReported);

    % Cast bsearch targets to the storage type of evT so mixed-type
    % comparisons don't fall back to double arithmetic on every step.
    if isa(evT, 'int64')
        tFrom = int64(tNow - dt);
        tNow_q = int64(tNow);
    else
        tFrom = tNow - dt;
        tNow_q = tNow;
    end
    iStart = bsearchLeft(evT, tFrom);
    iEnd   = bsearchRight(evT, tNow_q);
    nEv = iEnd - iStart + 1;
    stats(2) = double(nEv > 0);    % eventsSliced (this window had ANY events)
    if nEv < minEventsPerWindow, return; end

    % evX / evY may be uint16 (compact storage); cast to double for the
    % logical comparisons and accumarray subscripts.
    wX = double(evX(iStart:iEnd));
    wY = double(evY(iStart:iEnd));
    validMask = (wY >= 1) & (wY <= H) & (wX >= 1) & (wX <= W);
    wXv = wX(validMask);
    wYv = wY(validMask);
    if isempty(wXv), return; end

    % --- Always build the full-resolution count image (used for decoding) ---
    countImg = accumarray([wYv, wXv], 1, [H, W]);

    % --- Blob detection mask: full res, or downscaled when detectScale<1 ---
    if detectScale < 1
        sH = max(8, round(H * detectScale));
        sW = max(8, round(W * detectScale));
        smallMask = imresize(countImg, [sH sW], 'nearest') > 0;
        % adjust blob area thresholds to the downscaled image
        scaledParams = blobParams;
        scaleArea = detectScale^2;
        scaledParams.minArea = blobParams.minArea * scaleArea;
        scaledParams.maxArea = blobParams.maxArea * scaleArea;
        quadsSmall = detectQuadBlob(smallMask, scaledParams);
        stats(3) = length(quadsSmall);          % quadsFound
        if isempty(quadsSmall), return; end
        % Scale corners back to full resolution
        invScale = 1 / detectScale;
        quads = cell(1, length(quadsSmall));
        for qi = 1:length(quadsSmall)
            quads{qi} = quadsSmall{qi} * invScale;
        end
    else
        activeMask = countImg > 0;
        quads = detectQuadBlob(activeMask, blobParams);
        stats(3) = length(quads);                % quadsFound
        if isempty(quads), return; end
    end

    countU8 = uint8(countImg / max(countImg(:)) * 255);

    % --- Optional GPU offload of the per-quad warps -------------------
    % imwarp is the dominant per-quad cost; pushing the source image to
    % the GPU once and running all warps there saves the warp time on
    % large sensors (e.g. 480x640) when several quads survive blob
    % filtering. We gather only the small 160x160 warped tile.
    countSrc = countU8;
    if useGPU
        try
            countSrc = gpuArray(countU8);
        catch
            % keep CPU on any failure
            countSrc = countU8;
        end
    end

    for qi = 1:length(quads)
        corners = orderCornersForUnwarp_local(quads{qi});

        % Sub-pixel corner refinement: walk perpendicular to each rough
        % edge, find the gradient peak in countImg, fit a line, and
        % intersect adjacent lines. Without this, the min-area-rect
        % corners are typically 2-4 px off on noisy event blobs, which
        % shifts every cell-boundary sample in the unwarped image and
        % makes the decoder miss even with Hamming tolerance.
        if refineCorners
            refined = refineQuadCorners_local(countImg, corners, refineSearchPx);
            if ~any(isnan(refined(:))) && ~any(isinf(refined(:)))
                corners = orderCornersForUnwarp_local(refined);
            end
        end

        warpedImg = unwarpQuad_local(countSrc, corners, markerCoords, sideSize);
        if isempty(warpedImg), continue; end
        if isa(warpedImg, 'gpuArray')
            warpedImg = gather(warpedImg);
        end

        stats(4) = stats(4) + 1;          % decodeAttempts
        mid = decodeMarker_local( ...
            warpedImg, numCells, codeSize, cellPx, dictCodes, dictIDs, hammingThresh);
        if mid < 0, continue; end
        stats(5) = stats(5) + 1;          % decodeSuccess (any ID, before filter)

        % Filter: in legacy single-marker mode, accept any decoded ID
        % that passes requestedMarkerIds. In multi-marker mode the
        % filter is implicit in reportedIds (== requestedMarkerIds).
        if ~isempty(requestedMarkerIds) && ...
                ~any(requestedMarkerIds == int32(mid))
            continue;
        end
        stats(6) = stats(6) + 1;          % decodeInFilter

        % Backward-compat scalar: first match wins.
        if bestID < 0
            bestID = mid;
        end

        % Multi-marker mask: which reported ID did this quad match?
        if nReported > 0
            idx = find(reportedIds == int32(mid), 1);
            if ~isempty(idx)
                detVec(idx) = true;
            end
            % Stop early only when every reported marker is in.
            if all(detVec)
                return;
            end
        else
            % Legacy mode: original "first hit wins" early exit.
            return;
        end
    end
end


%% =========================================================================
%                    SUB-PIXEL CORNER REFINEMENT
%% =========================================================================
function refined = refineQuadCorners_local(countImg, roughCorners, searchHalfWidth)
    % Improve quad corners by fitting lines to the actual event-density
    % ridges that form each edge.
    %
    % For each of the 4 edges:
    %   1) walk along the rough edge in small steps,
    %   2) at each step sample a perpendicular line of countImg values,
    %   3) find the position of the peak (parabolic sub-pixel fit),
    %   4) fit a straight line through every collected peak.
    % Then intersect adjacent fitted lines to recover the refined corners.
    %
    % If anything goes wrong (degenerate line, near-parallel intersection,
    % too few valid samples), we fall back to the corresponding rough
    % corner for that vertex.

    [H, W] = size(countImg);
    edgeLines = zeros(4, 3);             % each row: [a b c]  ax+by+c=0
    countImgD = double(countImg);

    for ei = 1:4
        p1 = roughCorners(ei, :);
        p2 = roughCorners(mod(ei, 4) + 1, :);

        edgeVec = p2 - p1;
        edgeLen = norm(edgeVec);
        if edgeLen < 6
            edgeLines(ei, :) = pointsToLine_local(p1, p2);
            continue;
        end

        edgeUnit = edgeVec / edgeLen;
        perpUnit = [-edgeUnit(2), edgeUnit(1)];

        % Sample positions along the edge, stay clear of the corners.
        nSamples = max(8, floor(edgeLen / 2));
        ts = (0.15 + (0:nSamples-1)' * (0.70 / max(nSamples - 1, 1)));
        centers = p1 + ts .* edgeVec;

        offsets = (-searchHalfWidth:searchHalfWidth);
        nOff = length(offsets);

        % Build (nSamples x nOff) sample positions perpendicular to the edge
        xs = round(centers(:, 1) + offsets * perpUnit(1));
        ys = round(centers(:, 2) + offsets * perpUnit(2));

        inB = xs >= 1 & xs <= W & ys >= 1 & ys <= H;
        xsC = min(max(xs, 1), W);
        ysC = min(max(ys, 1), H);
        linIdx = sub2ind([H W], ysC, xsC);
        vals = countImgD(linIdx);
        vals(~inB) = 0;

        % Peak per row (per along-edge sample point)
        [maxVals, maxIdx] = max(vals, [], 2);

        peakPts = zeros(nSamples, 2);
        keep    = false(nSamples, 1);
        for s = 1:nSamples
            if maxVals(s) <= 0, continue; end
            mi = maxIdx(s);
            peakOff = offsets(mi);

            % Parabolic sub-pixel refinement when the peak isn't on the edge.
            if mi > 1 && mi < nOff && vals(s, mi-1) > 0 && vals(s, mi+1) > 0
                v1 = vals(s, mi - 1);
                v2 = vals(s, mi);
                v3 = vals(s, mi + 1);
                denom = v1 - 2*v2 + v3;
                if abs(denom) > 1e-6
                    sub = 0.5 * (v1 - v3) / denom;
                    if abs(sub) <= 1, peakOff = peakOff + sub; end
                end
            end

            peakPts(s, :) = centers(s, :) + peakOff * perpUnit;
            keep(s) = true;
        end

        peakPts = peakPts(keep, :);
        if size(peakPts, 1) < 3
            edgeLines(ei, :) = pointsToLine_local(p1, p2);
        else
            edgeLines(ei, :) = fitLine_local(peakPts);
        end
    end

    % Intersect adjacent edges to recover corners.
    refined = zeros(4, 2);
    for ei = 1:4
        prevEi = mod(ei - 2, 4) + 1;
        ipt = intersectLines_local(edgeLines(prevEi, :), edgeLines(ei, :));
        % Reject crazy results -- fall back to the rough corner.
        if any(isnan(ipt)) || any(isinf(ipt)) || ...
                ipt(1) < -100 || ipt(1) > W + 100 || ...
                ipt(2) < -100 || ipt(2) > H + 100 || ...
                norm(ipt - roughCorners(ei, :)) > 2 * searchHalfWidth + 5
            refined(ei, :) = roughCorners(ei, :);
        else
            refined(ei, :) = ipt;
        end
    end
end


function L = pointsToLine_local(p1, p2)
    a = p2(2) - p1(2);
    b = p1(1) - p2(1);
    c = p2(1) * p1(2) - p1(1) * p2(2);
    n = sqrt(a*a + b*b);
    if n > 1e-9
        L = [a/n, b/n, c/n];
    else
        L = [1 0 0];
    end
end


function L = fitLine_local(pts)
    % Total least squares line fit (orthogonal regression).
    centroid = mean(pts, 1);
    centered = pts - centroid;
    [~, ~, V] = svd(centered, 0);
    normal = V(:, 2);                    % smallest singular vector
    a = normal(1); b = normal(2);
    c = -(a * centroid(1) + b * centroid(2));
    n = sqrt(a*a + b*b);
    if n > 1e-9
        L = [a/n, b/n, c/n];
    else
        L = [1 0 0];
    end
end


function pt = intersectLines_local(L1, L2)
    a = L1(1); b = L1(2); c = L1(3);
    d = L2(1); e = L2(2); f = L2(3);
    denom = a * e - b * d;
    if abs(denom) < 1e-9
        pt = [NaN, NaN];
        return;
    end
    pt = [(b*f - c*e) / denom, (c*d - a*f) / denom];
end

function corners = orderCornersForUnwarp_local(corners)
    centroid = mean(corners, 1);
    angles = atan2(corners(:,2)-centroid(2), corners(:,1)-centroid(1));
    [~, si] = sort(angles);
    corners = corners(si, :);
    sums = corners(:,1) + corners(:,2);
    [~, tl] = min(sums);
    corners = circshift(corners, -(tl-1), 1);
    v1 = corners(2,:) - corners(1,:);
    v2 = corners(4,:) - corners(1,:);
    if v1(1)*v2(2) - v1(2)*v2(1) > 0
        corners = corners([1 4 3 2], :);
    end
end

function warped = unwarpQuad_local(img, srcCorners, dstCorners, sideSize)
    warped = [];
    try
        tform = fitgeotrans(srcCorners, dstCorners+1, 'projective');
        warped = imwarp(img, tform, 'OutputView', ...
            imref2d([sideSize sideSize], [1 sideSize], [1 sideSize]), ...
            'InterpolationMethod', 'bilinear');
    catch
        warped = [];
    end
end

function markerID = decodeMarker_local( ...
        warpedImg, numCells, codeSize, cellPx, dictCodes, dictIDs, hammingThresh)
    % hammingThresh : 0 (default) -> exact bit-for-bit dictionary match.
    %                 k > 0       -> accept the closest dictionary code if
    %                                its Hamming distance to the decoded
    %                                36-bit pattern is <= k. ARUCO_MIP_36h12
    %                                has minimum inter-marker Hamming
    %                                distance 12, so k up to ~4 still gives
    %                                unambiguous decoding.
    if nargin < 7, hammingThresh = 0; end

    markerID = -1;
    bestApproxID = -1;
    bestApproxDist = inf;

    warpedDbl = double(warpedImg);
    if max(warpedDbl(:)) == 0, return; end

    boundaryHalfWidth = 5;

    % ---- Compute transition threshold ----
    nHBnd = (numCells-1) * numCells;
    nVBnd = numCells * (numCells-1);
    boundaryVals = zeros(1, nHBnd + nVBnd);
    bi = 0;
    for row = 1:(numCells-1)
        bndRow = row * cellPx;
        r1 = max(1, bndRow - boundaryHalfWidth);
        r2 = min(size(warpedDbl,1), bndRow + boundaryHalfWidth);
        for col = 1:numCells
            cCenter = round((col-0.5) * cellPx);
            c1 = max(1, cCenter - 1);
            c2 = min(size(warpedDbl,2), cCenter + 1);
            bi = bi + 1;
            boundaryVals(bi) = mean(warpedDbl(r1:r2, c1:c2), 'all');
        end
    end
    for col = 1:(numCells-1)
        bndCol = col * cellPx;
        c1 = max(1, bndCol - boundaryHalfWidth);
        c2 = min(size(warpedDbl,2), bndCol + boundaryHalfWidth);
        for row = 1:numCells
            rCenter = round((row-0.5) * cellPx);
            r1 = max(1, rCenter - 1);
            r2 = min(size(warpedDbl,1), rCenter + 1);
            bi = bi + 1;
            boundaryVals(bi) = mean(warpedDbl(r1:r2, c1:c2), 'all');
        end
    end

    if max(boundaryVals) == 0, return; end
    transThresh = graythresh(uint8(boundaryVals / max(boundaryVals) * 255)) * max(boundaryVals);

    % ---- Vertical scan ----
    codeV = zeros(numCells);
    for col = 1:numCells
        cCenter = round((col-0.5) * cellPx);
        c1 = max(1, cCenter - 1);
        c2 = min(size(warpedDbl,2), cCenter + 1);
        currentColor = 0;
        codeV(1, col) = currentColor;
        for row = 2:numCells
            bndRow = (row-1) * cellPx;
            r1 = max(1, bndRow - boundaryHalfWidth);
            r2 = min(size(warpedDbl,1), bndRow + boundaryHalfWidth);
            bndVal = mean(warpedDbl(r1:r2, c1:c2), 'all');
            if bndVal > transThresh
                currentColor = 1 - currentColor;
            end
            codeV(row, col) = currentColor;
        end
    end

    % ---- Horizontal scan ----
    codeH = zeros(numCells);
    for row = 1:numCells
        rCenter = round((row-0.5) * cellPx);
        r1 = max(1, rCenter - 1);
        r2 = min(size(warpedDbl,1), rCenter + 1);
        currentColor = 0;
        codeH(row, 1) = currentColor;
        for col = 2:numCells
            bndCol = (col-1) * cellPx;
            c1 = max(1, bndCol - boundaryHalfWidth);
            c2 = min(size(warpedDbl,2), bndCol + boundaryHalfWidth);
            bndVal = mean(warpedDbl(r1:r2, c1:c2), 'all');
            if bndVal > transThresh
                currentColor = 1 - currentColor;
            end
            codeH(row, col) = currentColor;
        end
    end

    % ---- Combine: try vertical, horizontal, and majority vote ----
    candidates = {codeV, codeH, double((codeV + codeH) > 1)};

    for ci = 1:length(candidates)
        codeImg = candidates{ci};

        % --- Border-uniformity gate (cheap kill-switch for noise) -------
        % ARUCO_MIP_36h12 markers have a 1-cell-thick uniform border.
        % The outer ring is invariant under rotation and horizontal flip,
        % so we can check it ONCE per candidate (before trying any of the
        % 16 inv/flip/rot variants).
        %   * all 0 -> canonical orientation, only try inv=0
        %   * all 1 -> inverted orientation,  only try inv=1
        %   * mixed -> not a marker, skip the candidate entirely
        border = [codeImg(1, :), codeImg(end, :), ...
                  codeImg(2:end-1, 1)', codeImg(2:end-1, end)'];
        bSum = sum(border);
        if bSum == 0
            invList = 0;
        elseif bSum == numel(border)
            invList = 1;
        else
            continue;        % border not uniform -- noise, skip 16 variants
        end

        for inv = invList
            if inv == 1
                testCode = 1 - codeImg;
            else
                testCode = codeImg;
            end
            for doFlip = [0, 1]
                if doFlip == 1
                    testCodeF = fliplr(testCode);
                else
                    testCodeF = testCode;
                end
                inner = testCodeF(2:end-1, 2:end-1);
                for rot = 0:3
                    rotInner = rot90(inner, rot);
                    code = uint64(0);
                    for r = 1:codeSize
                        for c = 1:codeSize
                            bit = uint64(rotInner(r, c));
                            code = bitor(code, bitshift(bit, 36 - ((r-1)*codeSize + c)));
                        end
                    end
                    if hammingThresh == 0
                        idx = bsearchCode(dictCodes, code);
                        if idx > 0
                            markerID = dictIDs(idx);
                            return;
                        end
                    else
                        % Hamming-tolerant lookup: XOR with every dict
                        % code, count bits, take the minimum.
                        xors  = bitxor(uint64(code), dictCodes);
                        dists = popcount64Vec(xors);
                        [minDist, mi] = min(dists);
                        if minDist == 0
                            markerID = dictIDs(mi);
                            return;
                        end
                        if minDist <= hammingThresh && minDist < bestApproxDist
                            bestApproxDist = minDist;
                            bestApproxID   = dictIDs(mi);
                        end
                    end
                end
            end
        end
    end

    % After exhausting every variant: if we collected an approximate
    % match (Hamming threshold > 0 and best <= threshold), return it.
    if hammingThresh > 0 && bestApproxID >= 0
        markerID = bestApproxID;
    end
end

function n = popcount64Vec(x)
    % Population count (Hamming weight) for a vector/array of uint64 values.
    % Uses a 256-entry byte lookup table -- works around MATLAB's saturating
    % uint64 multiply, which makes the classic SWAR trick give wrong answers.
    persistent T8
    if isempty(T8)
        T8 = uint8(zeros(1, 256));
        for ii = 0:255
            T8(ii + 1) = sum(bitget(uint8(ii), 1:8));
        end
    end
    x = uint64(x);
    n = zeros(size(x));
    for k = 0:7
        byte = bitand(bitshift(x, -8*k), uint64(255));
        n = n + double(T8(double(byte) + 1));
    end
end

function idx = bsearchCode(sortedCodes, target)
    lo = 1; hi = length(sortedCodes);
    idx = 0;
    while lo <= hi
        mid = floor((lo + hi) / 2);
        if sortedCodes(mid) == target
            idx = mid;
            return;
        elseif sortedCodes(mid) < target
            lo = mid + 1;
        else
            hi = mid - 1;
        end
    end
end

function [sortedCodes, sortedIDs] = buildDictionaryArrays()
    codes = uint64([ ...
        hex2dec('d2b63a09d'), hex2dec('6001134e5'), hex2dec('1206fbe72'), hex2dec('ff8ad6cb4'), ...
        hex2dec('85da9bc49'), hex2dec('b461afe9c'), hex2dec('6db51fe13'), hex2dec('5248c541f'), ...
        hex2dec('8f34503'),   hex2dec('8ea462ece'), hex2dec('eac2be76d'), hex2dec('1af615c44'), ...
        hex2dec('b48a49f27'), hex2dec('2e4e1283b'), hex2dec('78b1f2fa8'), hex2dec('27d34f57e'), ...
        hex2dec('89222fff1'), hex2dec('4c1669406'), hex2dec('bf49b3511'), hex2dec('dc191cd5d'), ...
        hex2dec('11d7c3f85'), hex2dec('16a130e35'), hex2dec('e29f27eff'), hex2dec('428d8ae0c'), ...
        hex2dec('90d548477'), hex2dec('2319cbc93'), hex2dec('c3b0c3dfc'), hex2dec('424bccc9'),  ...
        hex2dec('2a081d630'), hex2dec('762743d96'), hex2dec('d0645bf19'), hex2dec('f38d7fd60'), ...
        hex2dec('c6cbf9a10'), hex2dec('3c1be7c65'), hex2dec('276f75e63'), hex2dec('4490a3f63'), ...
        hex2dec('da60acd52'), hex2dec('3cc68df59'), hex2dec('ab46f9dae'), hex2dec('88d533d78'), ...
        hex2dec('b6d62ec21'), hex2dec('b3c02b646'), hex2dec('22e56d408'), hex2dec('ac5f5770a'), ...
        hex2dec('aaa993f66'), hex2dec('4caa07c8d'), hex2dec('5c9b4f7b0'), hex2dec('aa9ef0e05'), ...
        hex2dec('705c5750'),  hex2dec('ac81f545e'), hex2dec('735b91e74'), hex2dec('8cc35cee4'), ...
        hex2dec('e44694d04'), hex2dec('b5e121de0'), hex2dec('261017d0f'), hex2dec('f1d439eb5'), ...
        hex2dec('a1a33ac96'), hex2dec('174c62c02'), hex2dec('1ee27f716'), hex2dec('8b1c5ece9'), ...
        hex2dec('6a05b0c6a'), hex2dec('d0568dfc'),  hex2dec('192d25e5f'), hex2dec('1adbeccc8'), ...
        hex2dec('cfec87f00'), hex2dec('d0b9dde7a'), hex2dec('88dcef81e'), hex2dec('445681cb9'), ...
        hex2dec('dbb2ffc83'), hex2dec('a48d96df1'), hex2dec('b72cc2e7d'), hex2dec('c295b53f'),  ...
        hex2dec('f49832704'), hex2dec('9968edc29'), hex2dec('9e4e1af85'), hex2dec('8683e2d1b'), ...
        hex2dec('810b45c04'), hex2dec('6ac44bfe2'), hex2dec('645346615'), hex2dec('3990bd598'), ...
        hex2dec('1c9ed0f6a'), hex2dec('c26729d65'), hex2dec('83993f795'), hex2dec('3ac05ac5d'), ...
        hex2dec('357adff3b'), hex2dec('d5c05565'),  hex2dec('2f547ef44'), hex2dec('86c115041'), ...
        hex2dec('640fd9e5f'), hex2dec('ce08bbcf7'), hex2dec('109bb343e'), hex2dec('c21435c92'), ...
        hex2dec('35b4dfce4'), hex2dec('459752cf2'), hex2dec('ec915b82c'), hex2dec('51881eed0'), ...
        hex2dec('2dda7dc97'), hex2dec('2e0142144'), hex2dec('42e890f99'), hex2dec('9a8856527'), ...
        hex2dec('8e80d9d80'), hex2dec('891cbcf34'), hex2dec('25dd82410'), hex2dec('239551d34'), ...
        hex2dec('8fe8f0c70'), hex2dec('94106a970'), hex2dec('82609b40c'), hex2dec('fc9caf36'),  ...
        hex2dec('688181d11'), hex2dec('718613c08'), hex2dec('f1ab7629'),  hex2dec('a357bfc18'), ...
        hex2dec('4c03b7a46'), hex2dec('204dedce6'), hex2dec('ad6300d37'), hex2dec('84cc4cd09'), ...
        hex2dec('42160e5c4'), hex2dec('87d2adfa8'), hex2dec('7850e7749'), hex2dec('4e750fc7c'), ...
        hex2dec('bf2e5dfda'), hex2dec('d88324da5'), hex2dec('234b52f80'), hex2dec('378204514'), ...
        hex2dec('abdf2ad53'), hex2dec('365e78ef9'), hex2dec('49caa6ca2'), hex2dec('3c39ddf3'),  ...
        hex2dec('c68c5385d'), hex2dec('5bfcbbf67'), hex2dec('623241e21'), hex2dec('abc90d5cc'), ...
        hex2dec('388c6fe85'), hex2dec('da0e2d62d'), hex2dec('10855dfe9'), hex2dec('4d46efd6b'), ...
        hex2dec('76ea12d61'), hex2dec('9db377d3d'), hex2dec('eed0efa71'), hex2dec('e6ec3ae2f'), ...
        hex2dec('441faee83'), hex2dec('ba19c8ff5'), hex2dec('313035eab'), hex2dec('6ce8f7625'), ...
        hex2dec('880dab58d'), hex2dec('8d3409e0d'), hex2dec('2be92ee21'), hex2dec('d60302c6c'), ...
        hex2dec('469ffc724'), hex2dec('87eebeed3'), hex2dec('42587ef7a'), hex2dec('7a8cc4e52'), ...
        hex2dec('76a437650'), hex2dec('999e41ef4'), hex2dec('7d0969e42'), hex2dec('c02baf46b'), ...
        hex2dec('9259f3e47'), hex2dec('2116a1dc0'), hex2dec('9f2de4d84'), hex2dec('effac29'),   ...
        hex2dec('7b371ff8c'), hex2dec('668339da9'), hex2dec('d010aee3f'), hex2dec('1cd00b4c0'), ...
        hex2dec('95070fc3b'), hex2dec('f84c9a770'), hex2dec('38f863d76'), hex2dec('3646ff045'), ...
        hex2dec('ce1b96412'), hex2dec('7a5d45da8'), hex2dec('14e00ef6c'), hex2dec('5e95abfd8'), ...
        hex2dec('b2e9cb729'), hex2dec('36c47dd7'),  hex2dec('b8ee97c6b'), hex2dec('e9e8f657'),  ...
        hex2dec('d4ad2ef1a'), hex2dec('8811c7f32'), hex2dec('47bde7c31'), hex2dec('3adadfb64'), ...
        hex2dec('6e5b28574'), hex2dec('33e67cd91'), hex2dec('2ab9fdd2d'), hex2dec('8afa67f2b'), ...
        hex2dec('e6a28fc5e'), hex2dec('72049cdbd'), hex2dec('ae65dac12'), hex2dec('1251a4526'), ...
        hex2dec('1089ab841'), hex2dec('e2f096ee0'), hex2dec('b0caee573'), hex2dec('fd6677e86'), ...
        hex2dec('444b3f518'), hex2dec('be8b3a56a'), hex2dec('680a75cfc'), hex2dec('ac02baea8'), ...
        hex2dec('97d815e1c'), hex2dec('1d4386e08'), hex2dec('1a14f5b0e'), hex2dec('e658a8d81'), ...
        hex2dec('a3868efa7'), hex2dec('3668a9673'), hex2dec('e8fc53d85'), hex2dec('2e2b7edd5'), ...
        hex2dec('8b2470f13'), hex2dec('f69795f32'), hex2dec('4589ffc8e'), hex2dec('2e2080c9c'), ...
        hex2dec('64265f7d'),  hex2dec('3d714dd10'), hex2dec('1692c6ef1'), hex2dec('3e67f2f49'), ...
        hex2dec('5041dad63'), hex2dec('1a1503415'), hex2dec('64c18c742'), hex2dec('a72eec35'),  ...
        hex2dec('1f0f9dc60'), hex2dec('a9559bc67'), hex2dec('f32911d0d'), hex2dec('21c0d4ffc'), ...
        hex2dec('e01cef5b0'), hex2dec('4e23a3520'), hex2dec('aa4f04e49'), hex2dec('e1c4fcc43'), ...
        hex2dec('208e8f6e8'), hex2dec('8486774a5'), hex2dec('9e98c7558'), hex2dec('2c59fb7dc'), ...
        hex2dec('9446a4613'), hex2dec('8292dcc2e'), hex2dec('4d61631'),   hex2dec('d05527809'), ...
        hex2dec('a0163852d'), hex2dec('8f657f639'), hex2dec('cca6c3e37'), hex2dec('cb136bc7a'), ...
        hex2dec('fc5a83e53'), hex2dec('9aa44fc30'), hex2dec('bdec1bd3c'), hex2dec('e020b9f7c'), ...
        hex2dec('4b8f35fb0'), hex2dec('b8165f637'), hex2dec('33dc88d69'), hex2dec('10a2f7e4d'), ...
        hex2dec('c8cb5ff53'), hex2dec('de259ff6b'), hex2dec('46d070dd4'), hex2dec('32d3b9741'), ...
        hex2dec('7075f1c04'), hex2dec('4d58dbea0') ...
    ]);
    ids = int32(0:(length(codes)-1));
    [sortedCodes, si] = sort(codes);
    sortedIDs = ids(si);
end
