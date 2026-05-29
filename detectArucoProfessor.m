function results = detectArucoProfessor(matFile, sensorSize, params)
%DETECTARUCOPROFESSOR Real-data-focused event-camera ArUco detector.
%
% Key differences from the original detector:
%   - uses a correct millisecond tick step by default
%   - fuses event count, time surface, ON-only, and OFF-only evidence
%   - scores only the requested marker IDs when provided
%   - accepts controlled Hamming error with a second-best margin
%   - uses simple temporal tracks as predicted fallback quads
%
% The output keeps the familiar win_<X>ms fields so existing analysis code
% can still inspect per-window detections.

params = fillDefaultsProfessor(params, sensorSize);

sensorH = sensorSize(1);
sensorW = sensorSize(2);
H = sensorH;
W = sensorW;

requestedMarkerIds = int32(params.requestedMarkerIds(:)');
if isempty(requestedMarkerIds)
    fprintf('Professor requested markers: ANY\n');
else
    fprintf('Professor requested markers: %s\n', mat2str(requestedMarkerIds));
end

[dictCodesById, dictIdsById] = loadMIP36h12CodesProfessor();
if isempty(requestedMarkerIds)
    targetIds = dictIdsById;
else
    bad = requestedMarkerIds < 0 | requestedMarkerIds >= numel(dictCodesById);
    if any(bad)
        error('detectArucoProfessor: requested IDs outside dictionary range: %s', ...
            mat2str(requestedMarkerIds(bad)));
    end
    targetIds = requestedMarkerIds;
end
targetCodes = dictCodesById(double(targetIds) + 1);

% Per-marker output is kept only for explicitly requested markers. In ANY
% mode the detector can still return scalar win_<X>ms IDs, but it avoids a
% huge numTicks x numWindows x 250 tensor.
reportedIds = requestedMarkerIds;
nReported = numel(reportedIds);

windowDurations_ms = double(params.windowDurations_ms(:)');
windowDurations_us = windowDurations_ms * 1000;
numWindows = numel(windowDurations_us);
tickStep_us = double(params.tickStep_us);

sideSize = params.numCells * params.cellPx;
markerCoords = [0 0; sideSize-1 0; sideSize-1 sideSize-1; 0 sideSize-1];

%% ---- Load data -----------------------------------------------------------
fprintf('Loading %s...\n', matFile);
tmp = load(matFile, 'events');
events = tmp.events;
if size(events, 2) < 4
    error('detectArucoProfessor: events must be Nx4 [x y polarity timestamp].');
end
fprintf('Loaded %d events | sensor %dx%d\n', size(events, 1), W, H);

evX = uint16(events(:, 1) + 1);
evY = uint16(events(:, 2) + 1);
evPol = int8(events(:, 3));
evT = int64(events(:, 4));
[evT, sortIdx] = sort(evT);
evX = evX(sortIdx);
evY = evY(sortIdx);
evPol = evPol(sortIdx);
clear events tmp sortIdx;

processingScale = double(params.processingScale);
if processingScale <= 0 || processingScale > 1
    error('detectArucoProfessor: processingScale must be in the interval (0, 1].');
end
if processingScale < 1
    procH = max(1, round(sensorH * processingScale));
    procW = max(1, round(sensorW * processingScale));
    evX = uint16(min(procW, max(1, floor((double(evX) - 1) * processingScale) + 1)));
    evY = uint16(min(procH, max(1, floor((double(evY) - 1) * processingScale) + 1)));
    H = procH;
    W = procW;

    areaScale = processingScale * processingScale;
    params.blobParams.minArea = max(20, round(params.blobParams.minArea * areaScale));
    params.blobParams.maxArea = max(params.blobParams.minArea + 1, ...
        round(params.blobParams.maxArea * areaScale));
    params.blobParams.minPixels = max(5, round(params.blobParams.minPixels * areaScale));
end

tMin = double(evT(1));
tMax = double(evT(end));
timelineStartWindow_ms = double(params.timelineStartWindow_ms);
if isempty(timelineStartWindow_ms) || isnan(timelineStartWindow_ms)
    timelineStartWindow_ms = max(windowDurations_ms);
end
timelineStart_us = max(max(windowDurations_us), timelineStartWindow_ms * 1000);
tStart = tMin + timelineStart_us;
tEnd = tMax;
numTicks = floor((tEnd - tStart) / tickStep_us) + 1;
if numTicks < 1
    error('detectArucoProfessor: recording shorter than largest window.');
end

fprintf('Time range: %.3fs to %.3fs\n', tMin / 1e6, tMax / 1e6);
if processingScale < 1
    fprintf('Processing scale: %.2f -> sensor %dx%d\n', processingScale, W, H);
end
fprintf('Ticks: %d (every %.1f ms)\n', numTicks, tickStep_us / 1000);
fprintf('Windows: %s ms\n', mat2str(windowDurations_ms));
if abs(timelineStart_us / 1000 - max(windowDurations_ms)) > eps
    fprintf('Timeline alignment start: %.1f ms (not processed as a window)\n', ...
        timelineStart_us / 1000);
end
fprintf('Time surface tau: %.1f ms | max Hamming: %d | track max age: %.0f ms\n', ...
    params.timeSurfaceTau_ms, params.maxHammingDist, params.trackMaxAge_ms);
fprintf('Runtime: %s | early exit: %s | window order: %s | masks: %s | quad finder: %s | max quads/window: %g\n', ...
    char(params.speedProfile), char(params.earlyExitMode), ...
    char(params.windowOrderMode), char(params.candidateMaskMode), ...
    char(params.quadFinderMode), params.maxQuadsPerWindow);
if isfinite(params.maxEventsPerWindow)
    fprintf('Event cap/window: %d (%s)\n', round(params.maxEventsPerWindow), ...
        char(params.eventLimitMode));
end

statNames = {'attempts', 'events_in_window', 'quads_after_merge', ...
             'predicted_quads', 'warp_decode_attempts', 'decode_hits', ...
             'accepted_hits', 'track_assisted_hits'};

%% ---- Output buffers ------------------------------------------------------
resultTable = -1 * ones(numTicks, 2 + numWindows);
perMarkerFlat = false(numTicks, max(numWindows * nReported, 1));
detectionsPerWindow = zeros(1, numWindows);
attemptedPerWindow = zeros(1, numWindows);
profStats = zeros(1, 8);
tracks = initTracksProfessor(reportedIds);
lastDetectedWindowIdx = NaN(1, max(nReported, 1));
previousTickHadDetection = false;
saveWindowStats = params.saveWindowStats;
if saveWindowStats
    windowStatsLog = zeros(numTicks, numWindows, numel(statNames));
else
    windowStatsLog = [];
end
saveCornerLog = params.saveCorners && nReported > 0;
if saveCornerLog
    cornerLog = NaN(numTicks, nReported, 4, 2);
    cornerWindowMs = NaN(numTicks, nReported);
    cornerScore = NaN(numTicks, nReported);
    cornerHamming = NaN(numTicks, nReported);
    cornerConfidence = NaN(numTicks, nReported);
    cornerBorderError = NaN(numTicks, nReported);
    cornerBoundaryScore = NaN(numTicks, nReported);
    cornerTrackAssisted = false(numTicks, nReported);
else
    cornerLog = [];
    cornerWindowMs = [];
    cornerScore = [];
    cornerHamming = [];
    cornerConfidence = [];
    cornerBorderError = [];
    cornerBoundaryScore = [];
    cornerTrackAssisted = [];
end

fprintf('\nProfessor detector processing %d ticks x %d windows = up to %d attempts...\n', ...
    numTicks, numWindows, numTicks * numWindows);

nextProgressPct = params.showProgressEveryPct;
nextCheckpointPct = params.checkpointEveryPct;
lastProgressTick = 0;
tStartRun = tic;

if params.saveCheckpoints && strlength(string(params.checkpointFile)) > 0
    fprintf('Checkpoint file: %s (every %d%%)\n', ...
        char(params.checkpointFile), params.checkpointEveryPct);
end

for tick = 1:numTicks
    tNow = tStart + (tick - 1) * tickStep_us;
    tickRow = -1 * ones(1, numWindows);
    tickPerMarker = false(numWindows, max(nReported, 1));
    tickDetections = emptyDetectionsProfessor();
    windowOrder = windowOrderForTickProfessor(numWindows, detectionsPerWindow, ...
        lastDetectedWindowIdx, previousTickHadDetection, params);

    for oi = 1:numel(windowOrder)
        wi = windowOrder(oi);
        attemptedPerWindow(wi) = attemptedPerWindow(wi) + 1;
        [bestID, detVec, detections, statsRow] = processWindowProfessor( ...
            tNow, windowDurations_us(wi), ...
            evT, evX, evY, evPol, H, W, ...
            params, markerCoords, sideSize, ...
            targetIds, targetCodes, reportedIds, tracks);

        tickRow(wi) = bestID;
        if bestID >= 0
            detectionsPerWindow(wi) = detectionsPerWindow(wi) + 1;
        end
        if nReported > 0
            tickPerMarker(wi, :) = detVec;
        end
        if ~isempty(detections)
            for di = 1:numel(detections)
                detections(di).windowIdx = wi;
                detections(di).windowMs = windowDurations_ms(wi);
            end
            tickDetections = appendDetectionsProfessor(tickDetections, detections);
        end
        profStats = profStats + statsRow;
        if saveWindowStats
            windowStatsLog(tick, wi, :) = statsRow;
        end

        if shouldStopWindowSearchProfessor(tickRow, tickPerMarker, nReported, params)
            break;
        end
    end

    [lastDetectedWindowIdx, previousTickHadDetection] = ...
        updateAdaptiveWindowMemoryProfessor(tickDetections, reportedIds, ...
        lastDetectedWindowIdx);

    tracks = updateTracksProfessor(tracks, tickDetections, tNow, params);

    resultTable(tick, :) = [tNow, double(any(tickRow >= 0)), tickRow];
    if nReported > 0
        perMarkerFlat(tick, :) = reshape(tickPerMarker, 1, []);
    end
    if saveCornerLog && ~isempty(tickDetections)
        [cornerLog, cornerWindowMs, cornerScore, cornerHamming, ...
            cornerConfidence, cornerBorderError, cornerBoundaryScore, ...
            cornerTrackAssisted] = ...
            updateCornerLogProfessor(cornerLog, cornerWindowMs, cornerScore, ...
            cornerHamming, cornerConfidence, cornerBorderError, ...
            cornerBoundaryScore, cornerTrackAssisted, tick, reportedIds, ...
            tickDetections);
    end

    pct = floor(100 * tick / numTicks);
    if pct >= nextProgressPct || tick == numTicks
        elapsed = toc(tStartRun);
        eta = elapsed / tick * (numTicks - tick);
        detectedTicksSoFar = sum(resultTable(1:tick, 2) > 0);
        detectionRateSoFar = 100 * detectedTicksSoFar / max(tick, 1);
        recentStart = lastProgressTick + 1;
        recentDetected = sum(resultTable(recentStart:tick, 2) > 0);
        recentTicks = tick - lastProgressTick;
        recentRate = 100 * recentDetected / max(recentTicks, 1);
        fprintf('  %3d%% processed (%d/%d ticks, t=%.2fs) elapsed %.0fs ETA %.0fs detection %.1f%% total, %.1f%% recent (%d/%d), attempts %d\n', ...
            pct, tick, numTicks, tNow / 1e6, elapsed, eta, ...
            detectionRateSoFar, recentRate, recentDetected, recentTicks, ...
            sum(attemptedPerWindow));

        if params.saveCheckpoints && strlength(string(params.checkpointFile)) > 0 && ...
                (pct >= nextCheckpointPct || tick == numTicks)
            saveCheckpointProfessor(params.checkpointFile, resultTable, perMarkerFlat, ...
                tick, numTicks, numWindows, nReported, windowDurations_ms, ...
                detectionsPerWindow, attemptedPerWindow, requestedMarkerIds, reportedIds, ...
                profStats, statNames, params, elapsed, saveCornerLog, ...
                cornerLog, cornerWindowMs, cornerScore, cornerHamming, ...
                cornerConfidence, cornerBorderError, cornerBoundaryScore, ...
                cornerTrackAssisted, saveWindowStats, windowStatsLog);
            nextCheckpointPct = nextCheckpointPct + params.checkpointEveryPct;
        end

        lastProgressTick = tick;
        nextProgressPct = nextProgressPct + params.showProgressEveryPct;
    end
end

elapsed = toc(tStartRun);

%% ---- Summary -------------------------------------------------------------
fprintf('\n========================================\n');
fprintf(' PROFESSOR DETECTOR SUMMARY\n');
fprintf('========================================\n');
fprintf('Total ticks:          %d\n', numTicks);
fprintf('Ticks with detection: %d (%.1f%%)\n', sum(resultTable(:, 2) > 0), ...
    100 * sum(resultTable(:, 2) > 0) / max(numTicks, 1));
fprintf('Elapsed time:         %.1fs\n', elapsed);

fprintf('\n--- Professor funnel ---\n');
for si = 1:numel(statNames)
    fprintf('%-24s %10d\n', statNames{si}, round(profStats(si)));
end

fprintf('\n--- Per-window detections ---\n');
for wi = 1:numWindows
    fprintf('%6d ms: %6d / %6d reached (%.1f%%)\n', windowDurations_ms(wi), ...
        detectionsPerWindow(wi), attemptedPerWindow(wi), ...
        100 * detectionsPerWindow(wi) / max(attemptedPerWindow(wi), 1));
end

%% ---- Build result struct -------------------------------------------------
results = struct();
results.algorithmName = 'detectArucoProfessor';
results.tNow_us = resultTable(:, 1);
results.anyDetected = resultTable(:, 2);
for wi = 1:numWindows
    fname = sprintf('win_%dms', windowDurations_ms(wi));
    results.(fname) = resultTable(:, 2 + wi);
end
results.windowDurations_ms = windowDurations_ms;
results.detectionsPerWindow = detectionsPerWindow;
results.attemptedPerWindow = attemptedPerWindow;
results.requestedMarkerIds = double(requestedMarkerIds);
results.markerIdsReported = double(reportedIds);
results.professorStats = profStats;
results.professorStatNames = statNames;
results.paramsUsed = params;
if saveWindowStats
    results.windowStats = windowStatsLog;
    results.windowStatNames = statNames;
end
if saveCornerLog
    results.cornerLog = cornerLog;
    results.cornerLogMarkerIds = double(reportedIds);
    results.cornerWindowMs = cornerWindowMs;
    results.cornerScore = cornerScore;
    results.cornerHamming = cornerHamming;
    results.cornerConfidence = cornerConfidence;
    results.cornerBorderError = cornerBorderError;
    results.cornerBoundaryScore = cornerBoundaryScore;
    results.cornerTrackAssisted = cornerTrackAssisted;
end

if nReported > 0
    perMarkerTensor = reshape(perMarkerFlat, numTicks, numWindows, nReported);
    fprintf('\n--- Per-marker detections ---\n');
    for ri = 1:nReported
        mid = double(reportedIds(ri));
        hits = perMarkerTensor(:, :, ri);
        anyForId = any(hits, 2);
        countPerWin = sum(hits, 1);
        fprintf('id=%d any-window: %d/%d (%.1f%%)\n', mid, sum(anyForId), ...
            numTicks, 100 * sum(anyForId) / max(numTicks, 1));

        results.(sprintf('anyDetected_id%d', mid)) = double(anyForId);
        results.(sprintf('detectionsPerWindow_id%d', mid)) = countPerWin;
        for wi = 1:numWindows
            results.(sprintf('win_%dms_id%d', windowDurations_ms(wi), mid)) = ...
                double(hits(:, wi));
        end
    end
end

fprintf('========================================\n');
end


function saveCheckpointProfessor(checkpointFile, resultTable, perMarkerFlat, ...
        processedTick, totalTicks, numWindows, nReported, windowDurations_ms, ...
        detectionsPerWindow, attemptedPerWindow, requestedMarkerIds, reportedIds, ...
        profStats, statNames, params, elapsed, saveCornerLog, ...
        cornerLog, cornerWindowMs, cornerScore, cornerHamming, ...
        cornerConfidence, cornerBorderError, cornerBoundaryScore, ...
        cornerTrackAssisted, saveWindowStats, windowStatsLog)

processedRows = 1:processedTick;
partial = struct();
partial.algorithmName = 'detectArucoProfessor';
partial.isPartial = true;
partial.processedTick = processedTick;
partial.totalTicks = totalTicks;
partial.elapsedSeconds = elapsed;
partial.tNow_us = resultTable(processedRows, 1);
partial.anyDetected = resultTable(processedRows, 2);
partial.windowDurations_ms = windowDurations_ms;
for wi = 1:numWindows
    fname = sprintf('win_%dms', windowDurations_ms(wi));
    partial.(fname) = resultTable(processedRows, 2 + wi);
end

partial.detectionsPerWindow = detectionsPerWindow;
partial.attemptedPerWindow = attemptedPerWindow;
partial.requestedMarkerIds = double(requestedMarkerIds);
partial.markerIdsReported = double(reportedIds);
partial.professorStats = profStats;
partial.professorStatNames = statNames;
partial.paramsUsed = params;
if saveWindowStats
    partial.windowStats = windowStatsLog(processedRows, :, :);
    partial.windowStatNames = statNames;
end
if saveCornerLog
    partial.cornerLog = cornerLog(processedRows, :, :, :);
    partial.cornerLogMarkerIds = double(reportedIds);
    partial.cornerWindowMs = cornerWindowMs(processedRows, :);
    partial.cornerScore = cornerScore(processedRows, :);
    partial.cornerHamming = cornerHamming(processedRows, :);
    partial.cornerConfidence = cornerConfidence(processedRows, :);
    partial.cornerBorderError = cornerBorderError(processedRows, :);
    partial.cornerBoundaryScore = cornerBoundaryScore(processedRows, :);
    partial.cornerTrackAssisted = cornerTrackAssisted(processedRows, :);
end

if nReported > 0
    perMarkerTensor = reshape(perMarkerFlat(processedRows, :), ...
        processedTick, numWindows, nReported);
    for ri = 1:nReported
        mid = double(reportedIds(ri));
        hits = perMarkerTensor(:, :, ri);
        partial.(sprintf('anyDetected_id%d', mid)) = double(any(hits, 2));
        partial.(sprintf('detectionsPerWindow_id%d', mid)) = sum(hits, 1);
        for wi = 1:numWindows
            partial.(sprintf('win_%dms_id%d', windowDurations_ms(wi), mid)) = ...
                double(hits(:, wi));
        end
    end
end

try
    save(checkpointFile, '-struct', 'partial', '-v7.3');
    fprintf('      checkpoint saved: %s\n', char(checkpointFile));
catch ME
    warning('detectArucoProfessor:checkpointFailed', ...
        'Could not save checkpoint %s: %s', char(checkpointFile), ME.message);
end
end


function [cornerLog, cornerWindowMs, cornerScore, cornerHamming, ...
        cornerConfidence, cornerBorderError, cornerBoundaryScore, ...
        cornerTrackAssisted] = updateCornerLogProfessor( ...
        cornerLog, cornerWindowMs, cornerScore, cornerHamming, ...
        cornerConfidence, cornerBorderError, cornerBoundaryScore, ...
        cornerTrackAssisted, tick, reportedIds, detections)

ids = [detections.id];
for ri = 1:numel(reportedIds)
    markerId = double(reportedIds(ri));
    idxs = find(ids == markerId);
    if isempty(idxs)
        continue;
    end

    [~, bestLocal] = min([detections(idxs).score]);
    det = detections(idxs(bestLocal));
    if isempty(det.corners) || ~isequal(size(det.corners), [4 2])
        continue;
    end

    cornerLog(tick, ri, :, :) = det.corners;
    cornerWindowMs(tick, ri) = det.windowMs;
    cornerScore(tick, ri) = det.score;
    cornerHamming(tick, ri) = det.hamming;
    cornerConfidence(tick, ri) = det.confidence;
    cornerBorderError(tick, ri) = det.borderError;
    cornerBoundaryScore(tick, ri) = det.boundaryScore;
    cornerTrackAssisted(tick, ri) = det.trackAssisted;
end
end


%% =========================================================================
%                              WINDOW PROCESSING
%% =========================================================================
function [bestID, detVec, detections, stats] = processWindowProfessor( ...
        tNow, dt, evT, evX, evY, evPol, H, W, ...
        params, markerCoords, sideSize, targetIds, targetCodes, reportedIds, tracks)

stats = zeros(1, 8);
stats(1) = 1;
bestID = -1;
detections = emptyDetectionsProfessor();
detVec = false(1, numel(reportedIds));

if isa(evT, 'int64')
    tFrom = int64(tNow - dt);
    tNowQ = int64(tNow);
else
    tFrom = tNow - dt;
    tNowQ = tNow;
end

iStart = bsearchLeftProfessor(evT, tFrom);
iEnd = bsearchRightProfessor(evT, tNowQ);
nEv = iEnd - iStart + 1;
stats(2) = max(nEv, 0);
if nEv < params.minEventsPerWindow
    return;
end

winIdx = selectWindowEventIndicesProfessor(iStart, iEnd, ...
    params.maxEventsPerWindow, params.eventLimitMode);
stats(2) = numel(winIdx);
imgs = accumulateWindowProfessor( ...
    evX(winIdx), evY(winIdx), evPol(winIdx), evT(winIdx), ...
    tNow, H, W, params);
if imgs.validEvents < params.minEventsPerWindow
    return;
end

predQuads = predictedTrackQuadsProfessor(tracks, tNow, H, W, params);
stats(4) = numel(predQuads);
if params.tryTrackFirst && ~isempty(predQuads)
    trackParams = params;
    trackParams.decodeImageMode = params.trackFirstDecodeImageMode;
    trackParams.usePolarityDecode = params.trackFirstUsePolarityDecode;
    [bestID, detVec, detections, evalStats] = evaluateQuadsProfessor( ...
        predQuads, imgs, trackParams, markerCoords, sideSize, targetIds, ...
        targetCodes, reportedIds, tracks, tNow, H, W, params.trackFirstRefineCorners);
    stats = stats + evalStats;
    if bestID >= 0
        return;
    end
end

quads = proposeQuadsProfessor(imgs, params);
if ~isempty(predQuads)
    quads = [quads, predQuads]; %#ok<AGROW>
end
quads = deduplicateQuadsProfessor(quads);
quads = limitQuadsProfessor(quads, params.maxQuadsPerWindow);
stats(3) = numel(quads);
if isempty(quads)
    return;
end

[bestID, detVec, detections, evalStats] = evaluateQuadsProfessor( ...
    quads, imgs, params, markerCoords, sideSize, targetIds, targetCodes, ...
    reportedIds, tracks, tNow, H, W, params.refineCorners);
stats = stats + evalStats;
end


function [bestID, detVec, detections, stats] = evaluateQuadsProfessor( ...
        quads, imgs, params, markerCoords, sideSize, targetIds, ...
        targetCodes, reportedIds, tracks, tNow, H, W, doRefine)

stats = zeros(1, 8);
bestID = -1;
bestScore = inf;
detections = emptyDetectionsProfessor();
detVec = false(1, numel(reportedIds));

for qi = 1:numel(quads)
    rough = orderCornersProfessor(quads{qi});
    cornerSets = {rough};
    if doRefine
        refined = refineQuadCornersProfessor(imgs.fusedImg, rough, params.refineSearchPx);
        if ~isempty(refined) && all(isfinite(refined(:)))
            cornerSets = {orderCornersProfessor(refined), rough};
        end
    end

    for ci = 1:numel(cornerSets)
        corners = cornerSets{ci};
        if ~quadInReasonableBoundsProfessor(corners, H, W)
            continue;
        end

        warped.count = unwarpProfessor(imgs.countNorm, corners, markerCoords, sideSize);
        if isempty(warped.count), continue; end
        warped.time = unwarpProfessor(imgs.timeNorm, corners, markerCoords, sideSize);
        if params.usePolarityDecode
            warped.on = unwarpProfessor(imgs.onNorm, corners, markerCoords, sideSize);
            warped.off = unwarpProfessor(imgs.offNorm, corners, markerCoords, sideSize);
        end

        stats(5) = stats(5) + 1;
        hit = scoreWarpProfessor(warped, targetIds, targetCodes, params);
        if ~hit.accepted
            continue;
        end
        stats(6) = stats(6) + 1;

        hit.corners = corners;
        hit.score = hit.score + 0.5 * quadGeometryPenaltyProfessor(corners);
        [trackPenalty, trackAssisted] = trackPenaltyProfessor(hit.id, corners, tracks, tNow, params);
        hit.score = hit.score + trackPenalty;
        hit.trackAssisted = trackAssisted;
        if trackAssisted
            stats(8) = stats(8) + 1;
        end

        detections = storeDetectionProfessor(detections, hit);
        stats(7) = stats(7) + 1;
        if ~isempty(reportedIds)
            idx = find(reportedIds == int32(hit.id), 1);
            if ~isempty(idx)
                detVec(idx) = true;
            end
        end

        if hit.score < bestScore
            bestScore = hit.score;
            bestID = double(hit.id);
        end

        if params.stopAfterAllHitsInWindow && ~isempty(reportedIds) && all(detVec)
            return;
        end

        if params.stopAfterFirstHitInWindow
            return;
        end
    end
end

if ~isempty(reportedIds) && ~isempty(detections)
    for di = 1:numel(detections)
        idx = find(reportedIds == int32(detections(di).id), 1);
        if ~isempty(idx)
            detVec(idx) = true;
        end
    end
end
end


function stop = shouldStopWindowSearchProfessor(tickRow, tickPerMarker, nReported, params)
mode = char(lower(string(params.earlyExitMode)));
switch mode
    case 'none'
        stop = false;
    case 'any'
        stop = any(tickRow >= 0);
    case 'all'
        stop = nReported > 0 && all(any(tickPerMarker, 1));
    otherwise
        error('detectArucoProfessor: unknown earlyExitMode "%s".', mode);
end
end


function order = windowOrderForTickProfessor(numWindows, detectionsPerWindow, ...
        lastDetectedWindowIdx, previousTickHadDetection, params)
mode = char(lower(string(params.windowOrderMode)));
switch mode
    case 'fixed'
        order = 1:numWindows;
    case 'adaptive'
        fallbackIdx = mostSuccessfulWindowIdxProfessor(detectionsPerWindow);
        if previousTickHadDetection
            seeds = lastDetectedWindowIdx(isfinite(lastDetectedWindowIdx));
            seeds = unique(round(seeds), 'stable');
            if isempty(seeds)
                seeds = fallbackIdx;
            end
        else
            seeds = fallbackIdx;
        end
        order = neighborExpandedWindowOrderProfessor(numWindows, seeds);
    otherwise
        error('detectArucoProfessor: unknown windowOrderMode "%s".', mode);
end
end


function bestIdx = mostSuccessfulWindowIdxProfessor(detectionsPerWindow)
if isempty(detectionsPerWindow) || all(detectionsPerWindow <= 0)
    bestIdx = 1;
    return;
end
[~, bestIdx] = max(detectionsPerWindow);
end


function order = neighborExpandedWindowOrderProfessor(numWindows, seeds)
seeds = seeds(seeds >= 1 & seeds <= numWindows);
if isempty(seeds)
    seeds = 1;
end

order = zeros(1, numWindows);
used = false(1, numWindows);
count = 0;
for radius = 0:(numWindows - 1)
    for si = 1:numel(seeds)
        seed = seeds(si);
        if radius == 0
            candidates = seed;
        else
            candidates = [seed - radius, seed + radius];
        end
        for ci = 1:numel(candidates)
            idx = candidates(ci);
            if idx < 1 || idx > numWindows || used(idx)
                continue;
            end
            count = count + 1;
            order(count) = idx;
            used(idx) = true;
        end
    end
end

if count < numWindows
    missing = find(~used);
    order(count + (1:numel(missing))) = missing;
end
end


function [lastDetectedWindowIdx, hadDetection] = ...
        updateAdaptiveWindowMemoryProfessor(detections, reportedIds, ...
        lastDetectedWindowIdx)
hadDetection = ~isempty(detections);
if ~hadDetection
    return;
end

if isempty(reportedIds)
    [~, bestLocal] = min([detections.score]);
    wi = detections(bestLocal).windowIdx;
    if isfinite(wi)
        lastDetectedWindowIdx(1) = wi;
    end
    return;
end

ids = double([detections.id]);
for ri = 1:numel(reportedIds)
    markerId = double(reportedIds(ri));
    idxs = find(ids == markerId);
    if isempty(idxs)
        continue;
    end
    [~, bestLocal] = min([detections(idxs).score]);
    wi = detections(idxs(bestLocal)).windowIdx;
    if isfinite(wi)
        lastDetectedWindowIdx(ri) = wi;
    end
end
end


function idx = selectWindowEventIndicesProfessor(iStart, iEnd, maxEvents, mode)
n = iEnd - iStart + 1;
if n <= 0
    idx = [];
    return;
end
if ~isfinite(maxEvents) || n <= maxEvents
    idx = iStart:iEnd;
    return;
end

maxEvents = max(1, floor(maxEvents));
mode = char(lower(string(mode)));
switch mode
    case 'recent'
        idx = (iEnd - maxEvents + 1):iEnd;
    case 'uniform'
        idx = round(linspace(iStart, iEnd, maxEvents));
    otherwise
        error('detectArucoProfessor: unknown eventLimitMode "%s".', mode);
end
end


%% =========================================================================
%                              ACCUMULATION
%% =========================================================================
function imgs = accumulateWindowProfessor(wXRaw, wYRaw, wPolRaw, wTRaw, tNow, H, W, params)
wX = double(wXRaw);
wY = double(wYRaw);
wT = double(wTRaw);
wPol = double(wPolRaw);

valid = wX >= 1 & wX <= W & wY >= 1 & wY <= H;
wX = wX(valid);
wY = wY(valid);
wT = wT(valid);
wPol = wPol(valid);

imgs.validEvents = numel(wX);
if isempty(wX)
    imgs.countImg = zeros(H, W);
    imgs.timeImg = zeros(H, W);
    imgs.onImg = zeros(H, W);
    imgs.offImg = zeros(H, W);
    imgs.countNorm = zeros(H, W);
    imgs.timeNorm = zeros(H, W);
    imgs.onNorm = zeros(H, W);
    imgs.offNorm = zeros(H, W);
    imgs.fusedImg = zeros(H, W);
    return;
end

lin = sub2ind([H W], wY, wX);
countVec = accumarray(lin, 1, [H * W, 1], @sum, 0);
lastVec = accumarray(lin, wT, [H * W, 1], @max, -inf);

isOn = wPol > 0;
onVec = accumarray(lin(isOn), 1, [H * W, 1], @sum, 0);
offVec = accumarray(lin(~isOn), 1, [H * W, 1], @sum, 0);

imgs.countImg = reshape(countVec, [H W]);
imgs.onImg = reshape(onVec, [H W]);
imgs.offImg = reshape(offVec, [H W]);

lastImg = reshape(lastVec, [H W]);
tauUs = max(1, params.timeSurfaceTau_ms * 1000);
timeImg = exp(-(double(tNow) - lastImg) / tauUs);
timeImg(~isfinite(timeImg)) = 0;
timeImg(lastImg == -inf) = 0;
imgs.timeImg = timeImg;

imgs.countNorm = normalize01Professor(log1p(imgs.countImg));
imgs.timeNorm = normalize01Professor(imgs.timeImg);
imgs.onNorm = normalize01Professor(log1p(imgs.onImg));
imgs.offNorm = normalize01Professor(log1p(imgs.offImg));

polContrast = normalize01Professor(abs(imgs.onNorm - imgs.offNorm));
imgs.fusedImg = max(imgs.countNorm, params.timeSurfaceWeight * imgs.timeNorm);
imgs.fusedImg = max(imgs.fusedImg, 0.55 * polContrast);
end


function quads = proposeQuadsProfessor(imgs, params)
quads = {};

maskFused = maskFromEvidenceProfessor(imgs.fusedImg, params);
maskTime = imgs.timeNorm > params.timeSurfaceThreshold;
maskOn = imgs.onImg > 0;
maskOff = imgs.offImg > 0;
maskCount = imgs.countImg > 0;

maskMode = char(lower(string(params.candidateMaskMode)));
switch maskMode
    case 'fast'
        masks = {maskFused, cleanupMaskProfessor(maskTime)};
    case 'balanced'
        masks = {maskFused, cleanupMaskProfessor(maskTime), ...
                 cleanupMaskProfessor(maskOn | maskOff)};
    case 'full'
        masks = {maskFused, cleanupMaskProfessor(maskTime), ...
                 cleanupMaskProfessor(maskOn), cleanupMaskProfessor(maskOff), ...
                 cleanupMaskProfessor(maskCount)};
    otherwise
        error('detectArucoProfessor: unknown candidateMaskMode "%s".', maskMode);
end

finderMode = char(lower(string(params.quadFinderMode)));

for mi = 1:numel(masks)
    bw = masks{mi};
    if ~any(bw(:)), continue; end
    if any(strcmp(finderMode, {'blob', 'both'}))
        try
            q1 = detectQuadBlob(bw, params.blobParams);
            quads = [quads, q1]; %#ok<AGROW>
        catch
            % Keep going; some masks may be pathological.
        end
    end
    if any(strcmp(finderMode, {'connected', 'both'}))
        try
            qp = params.blobParams;
            qp.connectivity = 4;
            qp.minPixels = max(8, round(params.blobParams.minPixels));
            qp.minRectangularity = params.blobParams.minRectangularity;
            q2 = findQuadCandidates(bw, qp);
            quads = [quads, q2]; %#ok<AGROW>
        catch
        end
    end
end
quads = deduplicateQuadsProfessor(quads);
end


function mask = maskFromEvidenceProfessor(img, params)
if max(img(:)) <= 0
    mask = false(size(img));
    return;
end
level = graythresh(img);
thr = max(0.035, min(0.35, level * 0.70));
mask = img >= thr;
mask = cleanupMaskProfessor(mask);

% Keep the strongest sparse pixels too. This helps faint real edges that
% Otsu sometimes suppresses when the image is mostly empty.
vals = img(img > 0);
if ~isempty(vals)
    pct = prctile(vals, 82);
    mask = mask | cleanupMaskProfessor(img >= pct);
end
end


function mask = cleanupMaskProfessor(mask)
if ~any(mask(:)), return; end
mask = bwareaopen(mask, 5);
try
    mask = imclose(mask, strel('disk', 1));
catch
end
end


%% =========================================================================
%                              DECODING
%% =========================================================================
function hit = scoreWarpProfessor(warped, targetIds, targetCodes, params)
hit = emptyHitProfessor();

fused = max(normalize01Professor(warped.count), ...
            params.timeSurfaceWeight * normalize01Professor(warped.time));

if params.usePolarityDecode
    polAbs = normalize01Professor(abs(normalize01Professor(warped.on) - ...
                                      normalize01Professor(warped.off)));
    polAny = max(normalize01Professor(warped.on), normalize01Professor(warped.off));
else
    polAbs = zeros(size(fused));
    polAny = zeros(size(fused));
end

decodeMode = char(lower(string(params.decodeImageMode)));
switch decodeMode
    case 'fast'
        decodeImages = {fused};
    case 'balanced'
        decodeImages = {fused, max(fused, 0.65 * polAbs)};
    case 'full'
        decodeImages = {fused, max(fused, 0.65 * polAbs), ...
                        max(fused, 0.45 * polAny), polAbs};
    otherwise
        error('detectArucoProfessor: unknown decodeImageMode "%s".', decodeMode);
end

for ii = 1:numel(decodeImages)
    [codeCandidates, boundaryScore] = transitionDecodeCandidatesProfessor( ...
        decodeImages{ii}, params);
    if isempty(codeCandidates)
        continue;
    end
    candHit = scoreCodeCandidatesProfessor( ...
        codeCandidates, targetIds, targetCodes, boundaryScore, params);
    if candHit.accepted && candHit.score < hit.score
        hit = candHit;
    end
end
end


function [candidates, boundaryScore] = transitionDecodeCandidatesProfessor(img, params)
candidates = {};
boundaryScore = 0;
img = normalize01Professor(double(img));
if max(img(:)) <= 0
    return;
end

numCells = params.numCells;
cellPx = params.cellPx;
bw = params.boundaryHalfWidth;

nVals = 2 * (numCells - 1) * numCells;
boundaryVals = zeros(1, nVals);
bi = 0;

for row = 1:(numCells - 1)
    bndRow = row * cellPx;
    r1 = max(1, bndRow - bw);
    r2 = min(size(img, 1), bndRow + bw);
    for col = 1:numCells
        cCenter = round((col - 0.5) * cellPx);
        c1 = max(1, cCenter - 1);
        c2 = min(size(img, 2), cCenter + 1);
        bi = bi + 1;
        boundaryVals(bi) = mean(img(r1:r2, c1:c2), 'all');
    end
end

for col = 1:(numCells - 1)
    bndCol = col * cellPx;
    c1 = max(1, bndCol - bw);
    c2 = min(size(img, 2), bndCol + bw);
    for row = 1:numCells
        rCenter = round((row - 0.5) * cellPx);
        r1 = max(1, rCenter - 1);
        r2 = min(size(img, 1), rCenter + 1);
        bi = bi + 1;
        boundaryVals(bi) = mean(img(r1:r2, c1:c2), 'all');
    end
end

mx = max(boundaryVals);
if mx <= 0
    return;
end
level = graythresh(uint8(boundaryVals / mx * 255));
transThresh = max(0.05 * mx, level * mx * params.transitionThresholdScale);
strongVals = boundaryVals(boundaryVals > transThresh);
if ~isempty(strongVals)
    boundaryScore = min(1, mean(strongVals) / (mx + eps));
end

codeV = zeros(numCells);
for col = 1:numCells
    cCenter = round((col - 0.5) * cellPx);
    c1 = max(1, cCenter - 1);
    c2 = min(size(img, 2), cCenter + 1);
    currentColor = 0;
    codeV(1, col) = currentColor;
    for row = 2:numCells
        bndRow = (row - 1) * cellPx;
        r1 = max(1, bndRow - bw);
        r2 = min(size(img, 1), bndRow + bw);
        bndVal = mean(img(r1:r2, c1:c2), 'all');
        if bndVal > transThresh
            currentColor = 1 - currentColor;
        end
        codeV(row, col) = currentColor;
    end
end

codeH = zeros(numCells);
for row = 1:numCells
    rCenter = round((row - 0.5) * cellPx);
    r1 = max(1, rCenter - 1);
    r2 = min(size(img, 1), rCenter + 1);
    currentColor = 0;
    codeH(row, 1) = currentColor;
    for col = 2:numCells
        bndCol = (col - 1) * cellPx;
        c1 = max(1, bndCol - bw);
        c2 = min(size(img, 2), bndCol + bw);
        bndVal = mean(img(r1:r2, c1:c2), 'all');
        if bndVal > transThresh
            currentColor = 1 - currentColor;
        end
        codeH(row, col) = currentColor;
    end
end

candidates = {codeV, codeH, double((codeV + codeH) >= 1), double((codeV + codeH) > 1)};
end


function hit = scoreCodeCandidatesProfessor(codeCandidates, targetIds, targetCodes, boundaryScore, params)
hit = emptyHitProfessor();
numTargets = numel(targetCodes);
if numTargets == 0
    return;
end

for ci = 1:numel(codeCandidates)
    codeImg = codeCandidates{ci};
    for inv = [0 1]
        if inv == 1
            testCode = 1 - codeImg;
        else
            testCode = codeImg;
        end
        for doFlip = [0 1]
            if doFlip
                testCodeF = fliplr(testCode);
            else
                testCodeF = testCode;
            end

            border = [testCodeF(1, :), testCodeF(end, :), ...
                      testCodeF(2:end-1, 1)', testCodeF(2:end-1, end)'];
            borderError = sum(border) / numel(border);
            if borderError > params.maxBorderErrorFrac
                continue;
            end

            inner = testCodeF(2:end-1, 2:end-1);
            for rot = 0:3
                rotInner = rot90(inner, rot);
                code = packCodeProfessor(rotInner, params.codeSize);
                dists = popcount64VecProfessor(bitxor(uint64(code), targetCodes));
                [minDist, bestIdx] = min(dists);
                if numTargets > 1
                    sortedD = sort(dists);
                    secondDist = sortedD(2);
                else
                    secondDist = inf;
                end

                marginOk = numTargets == 1 || ...
                    (secondDist - minDist) >= params.secondBestMargin || minDist <= 1;
                if minDist <= params.maxHammingDist && marginOk
                    score = double(minDist) + ...
                        params.borderPenaltyWeight * borderError - ...
                        0.50 * boundaryScore;
                    if score < hit.score
                        hit.accepted = true;
                        hit.id = double(targetIds(bestIdx));
                        hit.score = score;
                        hit.hamming = double(minDist);
                        hit.secondBestDist = double(secondDist);
                        hit.borderError = borderError;
                        hit.boundaryScore = boundaryScore;
                        hit.confidence = max(0, 1 - double(minDist) / max(params.maxHammingDist, 1));
                        hit.codeImg = testCodeF;
                    end
                end
            end
        end
    end
end
end


function code = packCodeProfessor(bits, codeSize)
code = uint64(0);
for r = 1:codeSize
    for c = 1:codeSize
        bit = uint64(bits(r, c) > 0);
        shift = 36 - ((r - 1) * codeSize + c);
        code = bitor(code, bitshift(bit, shift));
    end
end
end


%% =========================================================================
%                              TRACKING
%% =========================================================================
function tracks = initTracksProfessor(reportedIds)
tracks = struct('id', {}, 'active', {}, 'corners', {}, 'center', {}, ...
    'velocity', {}, 'lastT', {}, 'confidence', {}, 'misses', {});
for i = 1:numel(reportedIds)
    tracks(i).id = double(reportedIds(i));
    tracks(i).active = false;
    tracks(i).corners = zeros(4, 2);
    tracks(i).center = [NaN NaN];
    tracks(i).velocity = [0 0];  % pixels / second
    tracks(i).lastT = -inf;
    tracks(i).confidence = 0;
    tracks(i).misses = 0;
end
end


function predQuads = predictedTrackQuadsProfessor(tracks, tNow, H, W, params)
predQuads = {};
for i = 1:numel(tracks)
    if ~tracks(i).active, continue; end
    ageUs = double(tNow) - double(tracks(i).lastT);
    if ageUs < 0 || ageUs > params.trackMaxAge_ms * 1000
        continue;
    end
    dtSec = ageUs / 1e6;
    shift = tracks(i).velocity * dtSec;
    q = tracks(i).corners + shift;
    if quadInReasonableBoundsProfessor(q, H, W)
        predQuads{end+1} = q; %#ok<AGROW>
    end
end
end


function tracks = updateTracksProfessor(tracks, detections, tNow, params)
for ti = 1:numel(tracks)
    if isempty(detections)
        tracks(ti).misses = tracks(ti).misses + 1;
        continue;
    end
    ids = [detections.id];
    idxs = find(ids == tracks(ti).id);
    if isempty(idxs)
        tracks(ti).misses = tracks(ti).misses + 1;
        continue;
    end
    [~, bestLocal] = min([detections(idxs).score]);
    det = detections(idxs(bestLocal));
    newCenter = mean(det.corners, 1);

    if tracks(ti).active && isfinite(tracks(ti).lastT)
        dtSec = max((double(tNow) - double(tracks(ti).lastT)) / 1e6, 1e-6);
        measVel = (newCenter - tracks(ti).center) / dtSec;
        a = params.trackSmoothing;
        tracks(ti).velocity = a * tracks(ti).velocity + (1 - a) * measVel;
    else
        tracks(ti).velocity = [0 0];
    end

    tracks(ti).active = true;
    tracks(ti).corners = det.corners;
    tracks(ti).center = newCenter;
    tracks(ti).lastT = tNow;
    tracks(ti).confidence = det.confidence;
    tracks(ti).misses = 0;
end

for ti = 1:numel(tracks)
    ageUs = double(tNow) - double(tracks(ti).lastT);
    if tracks(ti).active && ageUs > params.trackMaxAge_ms * 1000
        tracks(ti).active = false;
    end
end
end


function [penalty, assisted] = trackPenaltyProfessor(id, corners, tracks, tNow, params)
penalty = 0;
assisted = false;
if isempty(tracks), return; end
idx = find([tracks.id] == double(id), 1);
if isempty(idx) || ~tracks(idx).active
    return;
end
ageUs = double(tNow) - double(tracks(idx).lastT);
if ageUs < 0 || ageUs > params.trackMaxAge_ms * 1000
    return;
end

dtSec = ageUs / 1e6;
predCenter = tracks(idx).center + tracks(idx).velocity * dtSec;
center = mean(corners, 1);
sides = quadSideLengthsProfessor(corners);
scale = max(mean(sides), 1);
d = norm(center - predCenter) / scale;
penalty = params.trackPenaltyWeight * min(d, 3);
assisted = d < 0.65;
end


%% =========================================================================
%                              GEOMETRY
%% =========================================================================
function refined = refineQuadCornersProfessor(weightImg, roughCorners, searchHalfWidth)
refined = [];
[H, W] = size(weightImg);
img = normalize01Professor(weightImg);
edgeLines = zeros(4, 3);

for ei = 1:4
    p1 = roughCorners(ei, :);
    p2 = roughCorners(mod(ei, 4) + 1, :);
    edgeVec = p2 - p1;
    edgeLen = norm(edgeVec);
    if edgeLen < 6
        edgeLines(ei, :) = pointsToLineProfessor(p1, p2);
        continue;
    end
    edgeUnit = edgeVec / edgeLen;
    perpUnit = [-edgeUnit(2), edgeUnit(1)];
    nSamples = max(8, floor(edgeLen / 2));
    ts = linspace(0.12, 0.88, nSamples)';
    centers = p1 + ts .* edgeVec;
    offsets = -searchHalfWidth:searchHalfWidth;

    peakPts = zeros(nSamples, 2);
    keep = false(nSamples, 1);
    for si = 1:nSamples
        vals = zeros(1, numel(offsets));
        for oi = 1:numel(offsets)
            pt = centers(si, :) + offsets(oi) * perpUnit;
            x = round(pt(1));
            y = round(pt(2));
            if x >= 1 && x <= W && y >= 1 && y <= H
                vals(oi) = img(y, x);
            end
        end
        [mx, mi] = max(vals);
        if mx <= 0, continue; end
        peakPts(si, :) = centers(si, :) + offsets(mi) * perpUnit;
        keep(si) = true;
    end

    peakPts = peakPts(keep, :);
    if size(peakPts, 1) >= 3
        edgeLines(ei, :) = fitLineProfessor(peakPts);
    else
        edgeLines(ei, :) = pointsToLineProfessor(p1, p2);
    end
end

sides = quadSideLengthsProfessor(roughCorners);
refined = zeros(4, 2);
for ei = 1:4
    prev = mod(ei - 2, 4) + 1;
    pt = intersectLinesProfessor(edgeLines(prev, :), edgeLines(ei, :));
    cap = min(0.15 * min(sides(ei), sides(prev)), 2 * searchHalfWidth + 3);
    if any(~isfinite(pt)) || norm(pt - roughCorners(ei, :)) > max(cap, 1)
        refined(ei, :) = roughCorners(ei, :);
    else
        refined(ei, :) = pt;
    end
end
end


function warped = unwarpProfessor(img, srcCorners, dstCorners, sideSize)
warped = [];
try
    tform = fitgeotrans(srcCorners, dstCorners + 1, 'projective');
    warped = imwarp(img, tform, 'OutputView', ...
        imref2d([sideSize sideSize], [1 sideSize], [1 sideSize]), ...
        'InterpolationMethod', 'bilinear');
catch
    warped = [];
end
end


function corners = orderCornersProfessor(corners)
centroid = mean(corners, 1);
angles = atan2(corners(:, 2) - centroid(2), corners(:, 1) - centroid(1));
[~, si] = sort(angles);
corners = corners(si, :);
sums = corners(:, 1) + corners(:, 2);
[~, tl] = min(sums);
corners = circshift(corners, -(tl - 1), 1);
v1 = corners(2, :) - corners(1, :);
v2 = corners(4, :) - corners(1, :);
if v1(1) * v2(2) - v1(2) * v2(1) > 0
    corners = corners([1 4 3 2], :);
end
end


function ok = quadInReasonableBoundsProfessor(corners, H, W)
margin = 40;
ok = all(isfinite(corners(:))) && ...
    all(corners(:, 1) > -margin) && all(corners(:, 1) < W + margin) && ...
    all(corners(:, 2) > -margin) && all(corners(:, 2) < H + margin) && ...
    polyarea(corners(:, 1), corners(:, 2)) > 20;
end


function penalty = quadGeometryPenaltyProfessor(corners)
sides = quadSideLengthsProfessor(corners);
aspectPenalty = max(sides) / max(min(sides), 1) - 1;
anglePenalty = 0;
for i = 1:4
    a = corners(mod(i - 2, 4) + 1, :) - corners(i, :);
    b = corners(mod(i, 4) + 1, :) - corners(i, :);
    anglePenalty = anglePenalty + abs(dot(a, b)) / max(norm(a) * norm(b), eps);
end
penalty = aspectPenalty + 0.25 * anglePenalty;
end


function sides = quadSideLengthsProfessor(corners)
sides = zeros(4, 1);
for i = 1:4
    sides(i) = norm(corners(i, :) - corners(mod(i, 4) + 1, :));
end
end


function L = pointsToLineProfessor(p1, p2)
a = p2(2) - p1(2);
b = p1(1) - p2(1);
c = p2(1) * p1(2) - p1(1) * p2(2);
n = sqrt(a * a + b * b);
if n > 1e-9
    L = [a / n, b / n, c / n];
else
    L = [1 0 0];
end
end


function L = fitLineProfessor(pts)
centroid = mean(pts, 1);
centered = pts - centroid;
[~, ~, V] = svd(centered, 0);
normal = V(:, 2);
a = normal(1);
b = normal(2);
c = -(a * centroid(1) + b * centroid(2));
n = sqrt(a * a + b * b);
if n > 1e-9
    L = [a / n, b / n, c / n];
else
    L = [1 0 0];
end
end


function pt = intersectLinesProfessor(L1, L2)
a = L1(1); b = L1(2); c = L1(3);
d = L2(1); e = L2(2); f = L2(3);
denom = a * e - b * d;
if abs(denom) < 1e-9
    pt = [NaN NaN];
else
    pt = [(b * f - c * e) / denom, (c * d - a * f) / denom];
end
end


function merged = deduplicateQuadsProfessor(quads)
if isempty(quads)
    merged = {};
    return;
end
centers = zeros(numel(quads), 2);
sizes = zeros(numel(quads), 1);
for k = 1:numel(quads)
    c = quads{k};
    centers(k, :) = mean(c, 1);
    sizes(k) = mean(quadSideLengthsProfessor(c));
end
keep = true(numel(quads), 1);
for k = 1:numel(quads)
    if ~keep(k), continue; end
    for m = k + 1:numel(quads)
        if ~keep(m), continue; end
        if norm(centers(k, :) - centers(m, :)) < 0.45 * min(sizes(k), sizes(m))
            keep(m) = false;
        end
    end
end
merged = quads(keep);
end


function limited = limitQuadsProfessor(quads, maxQuads)
if isempty(quads) || ~isfinite(maxQuads) || numel(quads) <= maxQuads
    limited = quads;
    return;
end

areas = zeros(numel(quads), 1);
for k = 1:numel(quads)
    q = quads{k};
    if size(q, 1) >= 3
        areas(k) = abs(polyarea(q(:, 1), q(:, 2)));
    end
end

[~, order] = sort(areas, 'descend');
keep = order(1:max(1, round(maxQuads)));
limited = quads(sort(keep));
end


%% =========================================================================
%                              STRUCT HELPERS
%% =========================================================================
function detections = emptyDetectionsProfessor()
detections = struct('accepted', {}, 'id', {}, 'score', {}, 'hamming', {}, ...
    'secondBestDist', {}, 'borderError', {}, 'boundaryScore', {}, ...
    'confidence', {}, 'codeImg', {}, 'corners', {}, 'windowIdx', {}, ...
    'windowMs', {}, 'trackAssisted', {});
end


function hit = emptyHitProfessor()
hit = struct('accepted', false, 'id', -1, 'score', inf, 'hamming', inf, ...
    'secondBestDist', inf, 'borderError', inf, 'boundaryScore', 0, ...
    'confidence', 0, 'codeImg', [], 'corners', [], 'windowIdx', NaN, ...
    'windowMs', NaN, 'trackAssisted', false);
end


function detections = storeDetectionProfessor(detections, hit)
if ~hit.accepted
    return;
end
if isempty(detections)
    detections = hit;
    return;
end
idx = find([detections.id] == hit.id, 1);
if isempty(idx)
    detections(end + 1) = hit;
elseif hit.score < detections(idx).score
    detections(idx) = hit;
end
end


function out = appendDetectionsProfessor(out, in)
if isempty(in), return; end
if isempty(out)
    out = in;
else
    out = [out, in]; %#ok<AGROW>
end
end


%% =========================================================================
%                              MISC HELPERS
%% =========================================================================
function params = fillDefaultsProfessor(params, sensorSize)
if nargin < 1 || isempty(params), params = struct(); end
if ~isfield(params, 'requestedMarkerIds'), params.requestedMarkerIds = []; end
if ~isfield(params, 'speedProfile'), params.speedProfile = "custom"; end
if ~isfield(params, 'processingScale'), params.processingScale = 1.0; end
if ~isfield(params, 'maxEventsPerWindow'), params.maxEventsPerWindow = inf; end
if ~isfield(params, 'eventLimitMode'), params.eventLimitMode = "uniform"; end
if ~isfield(params, 'timelineStartWindow_ms'), params.timelineStartWindow_ms = []; end
if ~isfield(params, 'saveCheckpoints'), params.saveCheckpoints = false; end
if ~isfield(params, 'checkpointEveryPct'), params.checkpointEveryPct = 2; end
if ~isfield(params, 'checkpointFile'), params.checkpointFile = ""; end
params.checkpointEveryPct = max(1, min(100, round(params.checkpointEveryPct)));
if ~isfield(params, 'saveCorners'), params.saveCorners = false; end
if ~isfield(params, 'saveWindowStats'), params.saveWindowStats = false; end
if ~isfield(params, 'windowDurations_ms'), params.windowDurations_ms = [5 10 20 50 100 250]; end
if ~isfield(params, 'tickStep_us'), params.tickStep_us = 5000; end
if ~isfield(params, 'numCells'), params.numCells = 8; end
if ~isfield(params, 'codeSize'), params.codeSize = 6; end
if ~isfield(params, 'cellPx'), params.cellPx = 20; end
if ~isfield(params, 'minEventsPerWindow'), params.minEventsPerWindow = 20; end
if ~isfield(params, 'timeSurfaceTau_ms'), params.timeSurfaceTau_ms = 25; end
if ~isfield(params, 'timeSurfaceWeight'), params.timeSurfaceWeight = 0.75; end
if ~isfield(params, 'timeSurfaceThreshold'), params.timeSurfaceThreshold = 0.08; end
if ~isfield(params, 'maxHammingDist'), params.maxHammingDist = 5; end
if ~isfield(params, 'secondBestMargin'), params.secondBestMargin = 2; end
if ~isfield(params, 'transitionThresholdScale'), params.transitionThresholdScale = 0.75; end
if ~isfield(params, 'boundaryHalfWidth'), params.boundaryHalfWidth = 5; end
if ~isfield(params, 'borderPenaltyWeight'), params.borderPenaltyWeight = 3.0; end
if ~isfield(params, 'maxBorderErrorFrac'), params.maxBorderErrorFrac = 0.45; end
if ~isfield(params, 'refineCorners'), params.refineCorners = true; end
if ~isfield(params, 'refineSearchPx'), params.refineSearchPx = 7; end
if ~isfield(params, 'trackMaxAge_ms'), params.trackMaxAge_ms = 250; end
if ~isfield(params, 'trackPenaltyWeight'), params.trackPenaltyWeight = 2.5; end
if ~isfield(params, 'trackSmoothing'), params.trackSmoothing = 0.65; end
if ~isfield(params, 'earlyExitOnAllMarkers'), params.earlyExitOnAllMarkers = true; end
if ~isfield(params, 'earlyExitMode')
    if params.earlyExitOnAllMarkers
        params.earlyExitMode = "all";
    else
        params.earlyExitMode = "none";
    end
end
if ~isfield(params, 'stopAfterFirstHitInWindow'), params.stopAfterFirstHitInWindow = false; end
if ~isfield(params, 'stopAfterAllHitsInWindow'), params.stopAfterAllHitsInWindow = false; end
if ~isfield(params, 'windowOrderMode'), params.windowOrderMode = "fixed"; end
if ~isfield(params, 'tryTrackFirst'), params.tryTrackFirst = false; end
if ~isfield(params, 'trackFirstRefineCorners'), params.trackFirstRefineCorners = false; end
if ~isfield(params, 'candidateMaskMode'), params.candidateMaskMode = "full"; end
if ~isfield(params, 'quadFinderMode'), params.quadFinderMode = "both"; end
if ~isfield(params, 'maxQuadsPerWindow'), params.maxQuadsPerWindow = inf; end
if ~isfield(params, 'decodeImageMode'), params.decodeImageMode = "full"; end
if ~isfield(params, 'usePolarityDecode'), params.usePolarityDecode = true; end
if ~isfield(params, 'trackFirstDecodeImageMode'), params.trackFirstDecodeImageMode = params.decodeImageMode; end
if ~isfield(params, 'trackFirstUsePolarityDecode'), params.trackFirstUsePolarityDecode = params.usePolarityDecode; end
if ~isfield(params, 'showProgressEveryPct'), params.showProgressEveryPct = 2; end
if ~isfield(params, 'saveDebug'), params.saveDebug = false; end
if ~isfield(params, 'blobParams') || isempty(params.blobParams)
    params.blobParams.minArea = 250;
    params.blobParams.maxArea = sensorSize(1) * sensorSize(2) * 0.35;
    params.blobParams.maxAspect = 3.0;
end
if ~isfield(params.blobParams, 'minArea'), params.blobParams.minArea = 250; end
if ~isfield(params.blobParams, 'maxArea'), params.blobParams.maxArea = sensorSize(1) * sensorSize(2) * 0.35; end
if ~isfield(params.blobParams, 'maxAspect'), params.blobParams.maxAspect = 3.0; end
if ~isfield(params.blobParams, 'minRectangularity'), params.blobParams.minRectangularity = 0.45; end
if ~isfield(params.blobParams, 'minPixels'), params.blobParams.minPixels = 15; end
end


function out = normalize01Professor(img)
img = double(img);
mx = max(img(:));
if mx <= 0
    out = zeros(size(img));
else
    out = img / mx;
end
end


function idx = bsearchLeftProfessor(evT, tTarget)
lo = 1;
hi = length(evT);
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


function idx = bsearchRightProfessor(evT, tTarget)
lo = 1;
hi = length(evT);
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


function n = popcount64VecProfessor(x)
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
    byte = bitand(bitshift(x, -8 * k), uint64(255));
    n = n + double(T8(double(byte) + 1));
end
end


function [codesById, idsById] = loadMIP36h12CodesProfessor()
% Keep one dictionary source of truth by reading the existing detector's
% embedded ARUCO_MIP_36h12 table. This avoids silently creating a second,
% possibly inconsistent copy in the experimental detector.
persistent cachedCodes cachedIds
if ~isempty(cachedCodes)
    codesById = cachedCodes;
    idsById = cachedIds;
    return;
end

rootDir = fileparts(mfilename('fullpath'));
src = fullfile(rootDir, 'detectAruco.m');
txt = fileread(src);
startIdx = regexp(txt, 'function \[sortedCodes, sortedIDs\] = buildDictionaryArrays', 'once');
if isempty(startIdx)
    error('detectArucoProfessor: could not find dictionary in detectAruco.m');
end
block = txt(startIdx:end);
tokens = regexp(block, 'hex2dec\(''([0-9a-fA-F]+)''\)', 'tokens');
if numel(tokens) < 250
    error('detectArucoProfessor: expected at least 250 ARUCO_MIP_36h12 codes, found %d.', ...
        numel(tokens));
end

codesById = zeros(1, 250, 'uint64');
for i = 1:250
    codesById(i) = uint64(hex2dec(tokens{i}{1}));
end
idsById = int32(0:249);

cachedCodes = codesById;
cachedIds = idsById;
end
