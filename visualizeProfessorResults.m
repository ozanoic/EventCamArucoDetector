function visualizeProfessorResults(eventMatFile, resultMatFile, opts)
%VISUALIZEPROFESSORRESULTS Replay event frames with Professor detections.
%
% Existing result files can show the event image, detected marker ID, and
% winning window for each tick. Rectangles are drawn only when the result
% file contains saved corner logs, because older Professor results did not
% store marker corner coordinates.
%
% Example:
%   visualizeProfessorResults( ...
%       "Data/OzanEventData_22.05.2026/4/4.mat", ...
%       "Data/OzanEventData_22.05.2026/4/4_professor_results.mat");
%
% Example video:
%   opts = struct('step', 5, 'saveVideo', true, 'videoFile', "review.mp4");
%   visualizeProfessorResults(eventFile, resultFile, opts);

if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = fillVisualizerOptions(opts);

results = load(resultMatFile);
tmp = load(eventMatFile, 'events');
events = tmp.events;
clear tmp;

if ~isfield(results, 'tNow_us')
    error('visualizeProfessorResults:MissingTime', ...
        'Result file must contain tNow_us.');
end
if size(events, 2) < 4
    error('visualizeProfessorResults:BadEvents', ...
        'Event file must contain events as Nx4 [x y polarity timestamp].');
end

[evT, order] = sort(double(events(:, 4)));
evX = double(events(order, 1)) + 1;
evY = double(events(order, 2)) + 1;
evPol = double(events(order, 3));
clear events order;

sensorSize = opts.sensorSize;
if isempty(sensorSize)
    sensorSize = [max(evY), max(evX)];
end
H = sensorSize(1);
W = sensorSize(2);

tickIdx = selectTicksForVisualization(results, opts);
if isempty(tickIdx)
    warning('visualizeProfessorResults:NoTicks', 'No ticks selected to visualize.');
    return;
end

hasCorners = resultHasCorners(results, opts.markerId);
if ~hasCorners
    warning('visualizeProfessorResults:NoCorners', ...
        ['This result file has marker IDs/windows, but no saved corner log. ' ...
         'Event frames will be shown without rectangles. Re-run a detector ' ...
         'with saveCorners=true to draw marker boxes.']);
end

fig = figure('Name', 'Professor ArUco Result Visualizer', 'Color', 'w');
ax = axes(fig);

writer = [];
if opts.saveVideo
    videoFile = string(opts.videoFile);
    if strlength(videoFile) == 0
        [resultDir, resultName, ~] = fileparts(resultMatFile);
        videoFile = fullfile(resultDir, resultName + "_visualization.mp4");
    end
    writer = VideoWriter(char(videoFile), 'MPEG-4');
    writer.FrameRate = opts.fps;
    open(writer);
end

cleanupObj = onCleanup(@() closeVideoWriter(writer));

for k = 1:numel(tickIdx)
    ti = tickIdx(k);
    tNow = double(results.tNow_us(ti));
    detections = detectionsAtTick(results, ti, opts.markerId);
    [~, detectedWindowMs] = primaryDetection(detections);
    eventWindowMs = opts.eventWindowMs;
    if isempty(eventWindowMs)
        if ~isnan(detectedWindowMs)
            eventWindowMs = detectedWindowMs;
        else
            eventWindowMs = opts.defaultEventWindowMs;
        end
    end

    tFrom = tNow - 1000 * eventWindowMs;
    iStart = bsearchLeftVisualizer(evT, tFrom);
    iEnd = bsearchRightVisualizer(evT, tNow);

    rgb = buildEventRgbFrame(evX(iStart:iEnd), evY(iStart:iEnd), ...
        evPol(iStart:iEnd), H, W, opts);

    image(ax, rgb);
    axis(ax, 'image');
    axis(ax, 'off');
    hold(ax, 'on');

    scale = getProcessingScale(results);
    palette = markerPalette();
    numDrawn = 0;
    for di = 1:numel(detections)
        corners = detections(di).corners;
        if isempty(corners)
            continue;
        end
        corners = corners ./ scale;
        color = palette(mod(di - 1, size(palette, 1)) + 1, :);
        drawMarkerBox(ax, corners, detections(di), color);
        numDrawn = numDrawn + 1;
    end

    if isempty(detections)
        detText = 'no marker';
    else
        detText = detectionSummaryText(detections);
    end

    title(ax, sprintf('tick %d/%d | t=%.3fs | %s', ...
        ti, numel(results.tNow_us), tNow / 1e6, detText), ...
        'Interpreter', 'none');

    if hasCorners && numDrawn == 0 && ~isempty(detections)
        text(ax, 10, 22, 'detection has no saved corners at this tick', ...
            'Color', 'y', 'FontWeight', 'bold', 'BackgroundColor', 'k', ...
            'Margin', 4);
    elseif ~hasCorners
        text(ax, 10, 22, 'corners not saved in this result file', ...
            'Color', 'y', 'FontWeight', 'bold', 'BackgroundColor', 'k', ...
            'Margin', 4);
    end

    hold(ax, 'off');
    drawnow;

    if opts.saveVideo
        writeVideo(writer, getframe(fig));
    end
    if opts.pauseSeconds > 0
        pause(opts.pauseSeconds);
    end
end
end


function opts = fillVisualizerOptions(opts)
defaults.sensorSize = [];
defaults.markerId = [];
defaults.tickIdx = [];
defaults.step = 1;
defaults.maxFrames = inf;
defaults.onlyDetections = false;
defaults.startTime_s = [];
defaults.endTime_s = [];
defaults.eventWindowMs = [];
defaults.defaultEventWindowMs = 10;
defaults.pauseSeconds = 0.01;
defaults.saveVideo = false;
defaults.videoFile = "";
defaults.fps = 20;
defaults.contrastPercentile = 99;

names = fieldnames(defaults);
for i = 1:numel(names)
    name = names{i};
    if ~isfield(opts, name)
        opts.(name) = defaults.(name);
    end
end
end


function tickIdx = selectTicksForVisualization(results, opts)
if ~isempty(opts.tickIdx)
    tickIdx = opts.tickIdx(:)';
else
    tickIdx = 1:numel(results.tNow_us);
end

if ~isempty(opts.startTime_s)
    tickIdx = tickIdx(results.tNow_us(tickIdx) >= opts.startTime_s * 1e6);
end
if ~isempty(opts.endTime_s)
    tickIdx = tickIdx(results.tNow_us(tickIdx) <= opts.endTime_s * 1e6);
end
if opts.onlyDetections && isfield(results, 'anyDetected')
    tickIdx = tickIdx(results.anyDetected(tickIdx) > 0);
end
tickIdx = tickIdx(1:max(1, round(opts.step)):end);
if isfinite(opts.maxFrames) && numel(tickIdx) > opts.maxFrames
    tickIdx = tickIdx(1:opts.maxFrames);
end
end


function rgb = buildEventRgbFrame(x, y, pol, H, W, opts)
valid = x >= 1 & x <= W & y >= 1 & y <= H;
x = x(valid);
y = y(valid);
pol = pol(valid);

on = pol > 0;
off = ~on;
onImg = accumarray([y(on), x(on)], 1, [H, W], @sum, 0);
offImg = accumarray([y(off), x(off)], 1, [H, W], @sum, 0);

mx = prctile([onImg(onImg > 0); offImg(offImg > 0)], opts.contrastPercentile);
if isempty(mx) || mx <= 0
    mx = 1;
end
onImg = min(1, log1p(onImg) / log1p(mx));
offImg = min(1, log1p(offImg) / log1p(mx));

rgb = zeros(H, W, 3);
rgb(:, :, 1) = offImg;
rgb(:, :, 2) = onImg;
rgb(:, :, 3) = 0.35 * offImg;
end


function [markerId, windowMs] = detectedMarkerAtTick(results, tickIdx, requestedId)
markerId = NaN;
windowMs = NaN;
if ~isfield(results, 'windowDurations_ms')
    return;
end

durations = double(results.windowDurations_ms(:)');
for i = 1:numel(durations)
    field = sprintf('win_%dms', durations(i));
    if ~isfield(results, field)
        continue;
    end
    id = results.(field)(tickIdx);
    if id < 0
        continue;
    end
    if ~isempty(requestedId) && id ~= requestedId
        continue;
    end
    markerId = id;
    windowMs = durations(i);
    return;
end
end


function detections = detectionsAtTick(results, tickIdx, requestedId)
detections = struct('id', {}, 'windowMs', {}, 'corners', {}, ...
    'score', {}, 'hamming', {});

ids = markerIdsForVisualization(results);
if isempty(ids)
    [markerId, windowMs] = detectedMarkerAtTick(results, tickIdx, requestedId);
    if ~isnan(markerId)
        detections(1).id = markerId;
        detections(1).windowMs = windowMs;
        detections(1).corners = getCornersForTick(results, tickIdx, markerId);
        detections(1).score = NaN;
        detections(1).hamming = NaN;
    end
    return;
end

for i = 1:numel(ids)
    markerId = ids(i);
    if ~isempty(requestedId) && markerId ~= requestedId
        continue;
    end

    [isDetected, windowMs] = markerDetectedAtTick(results, tickIdx, markerId);
    corners = getCornersForTick(results, tickIdx, markerId);
    if ~isDetected && isempty(corners)
        continue;
    end

    det.id = markerId;
    det.windowMs = windowMs;
    det.corners = corners;
    det.score = getMarkerMetricAtTick(results, 'cornerScore', tickIdx, markerId);
    det.hamming = getMarkerMetricAtTick(results, 'cornerHamming', tickIdx, markerId);
    detections(end + 1) = det; %#ok<AGROW>
end

if isempty(detections)
    [markerId, windowMs] = detectedMarkerAtTick(results, tickIdx, requestedId);
    if ~isnan(markerId)
        detections(1).id = markerId;
        detections(1).windowMs = windowMs;
        detections(1).corners = getCornersForTick(results, tickIdx, markerId);
        detections(1).score = NaN;
        detections(1).hamming = NaN;
    end
end
end


function ids = markerIdsForVisualization(results)
ids = [];
if isfield(results, 'cornerLogMarkerIds')
    ids = double(results.cornerLogMarkerIds(:)');
elseif isfield(results, 'markerIdsReported')
    ids = double(results.markerIdsReported(:)');
elseif isfield(results, 'requestedMarkerIds')
    ids = double(results.requestedMarkerIds(:)');
end
ids = ids(isfinite(ids));
end


function [detected, windowMs] = markerDetectedAtTick(results, tickIdx, markerId)
detected = false;
windowMs = NaN;

field = sprintf('anyDetected_id%d', markerId);
if isfield(results, field)
    detected = results.(field)(tickIdx) > 0;
end

if isfield(results, 'windowDurations_ms')
    durations = double(results.windowDurations_ms(:)');
    for i = 1:numel(durations)
        field = sprintf('win_%dms_id%d', durations(i), markerId);
        if isfield(results, field) && results.(field)(tickIdx) > 0
            detected = true;
            windowMs = durations(i);
            return;
        end
        field = sprintf('win_%dms', durations(i));
        if isfield(results, field) && results.(field)(tickIdx) == markerId
            detected = true;
            windowMs = durations(i);
            return;
        end
    end
end

if detected
    windowMs = getMarkerMetricAtTick(results, 'cornerWindowMs', tickIdx, markerId);
end
end


function [markerId, windowMs] = primaryDetection(detections)
markerId = NaN;
windowMs = NaN;
if isempty(detections)
    return;
end
markerId = detections(1).id;
windowMs = detections(1).windowMs;
end


function txt = detectionSummaryText(detections)
parts = strings(1, numel(detections));
for i = 1:numel(detections)
    if isnan(detections(i).windowMs)
        parts(i) = sprintf('id=%d', detections(i).id);
    else
        parts(i) = sprintf('id=%d/%gms', detections(i).id, detections(i).windowMs);
    end
end
txt = char(strjoin(parts, ', '));
end


function value = getMarkerMetricAtTick(results, fieldName, tickIdx, markerId)
value = NaN;
if ~isfield(results, fieldName)
    return;
end
ids = markerIdsForVisualization(results);
slot = find(ids == markerId, 1);
if isempty(slot)
    return;
end
data = results.(fieldName);
if size(data, 1) >= tickIdx && size(data, 2) >= slot
    value = data(tickIdx, slot);
end
end


function drawMarkerBox(ax, corners, det, color)
xy = [corners; corners(1, :)];
plot(ax, xy(:, 1), xy(:, 2), '-', 'Color', color, 'LineWidth', 2.5);
plot(ax, corners(:, 1), corners(:, 2), 'o', ...
    'Color', color, 'MarkerFaceColor', 'k', 'MarkerSize', 5);

labelPos = mean(corners, 1);
if isnan(det.windowMs)
    label = sprintf('id %d', det.id);
else
    label = sprintf('id %d | %g ms', det.id, det.windowMs);
end
if ~isnan(det.hamming)
    label = sprintf('%s | h=%g', label, det.hamming);
end

text(ax, labelPos(1), labelPos(2), label, ...
    'Color', color, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'k', ...
    'Margin', 4);
end


function colors = markerPalette()
colors = [
    1.00 0.92 0.10
    0.10 0.75 1.00
    0.20 1.00 0.35
    1.00 0.35 0.15
    0.85 0.35 1.00
    1.00 0.55 0.80
    0.35 1.00 0.85
    1.00 1.00 1.00
];
end


function tf = resultHasCorners(results, markerId)
tf = isfield(results, 'cornerLog') || isfield(results, 'corners');
if ~tf && ~isempty(markerId)
    tf = isfield(results, sprintf('corners_id%d', markerId));
end
end


function corners = getCornersForTick(results, tickIdx, markerId)
corners = [];
if isnan(markerId)
    return;
end

if isfield(results, 'cornerLog')
    ids = [];
    if isfield(results, 'cornerLogMarkerIds')
        ids = double(results.cornerLogMarkerIds(:)');
    elseif isfield(results, 'markerIdsReported')
        ids = double(results.markerIdsReported(:)');
    end
    if isempty(ids)
        markerSlot = 1;
    else
        markerSlot = find(ids == markerId, 1);
    end
    if ~isempty(markerSlot)
        c = squeeze(results.cornerLog(tickIdx, markerSlot, :, :));
        if isequal(size(c), [4 2]) && all(isfinite(c(:)))
            corners = c;
            return;
        end
    end
end

field = sprintf('corners_id%d', markerId);
if isfield(results, field)
    c = squeeze(results.(field)(tickIdx, :, :));
    if isequal(size(c), [4 2]) && all(isfinite(c(:)))
        corners = c;
        return;
    end
end

if isfield(results, 'corners')
    c = squeeze(results.corners(tickIdx, :, :));
    if isequal(size(c), [4 2]) && all(isfinite(c(:)))
        corners = c;
    end
end
end


function scale = getProcessingScale(results)
scale = 1.0;
if isfield(results, 'paramsUsed') && isfield(results.paramsUsed, 'processingScale')
    scale = double(results.paramsUsed.processingScale);
end
if ~isfinite(scale) || scale <= 0
    scale = 1.0;
end
end


function idx = bsearchLeftVisualizer(vals, target)
lo = 1;
hi = numel(vals);
if hi == 0
    idx = 1;
    return;
end
while lo < hi
    mid = floor((lo + hi) / 2);
    if vals(mid) < target
        lo = mid + 1;
    else
        hi = mid;
    end
end
idx = lo;
end


function idx = bsearchRightVisualizer(vals, target)
lo = 1;
hi = numel(vals);
if hi == 0
    idx = 0;
    return;
end
while lo < hi
    mid = ceil((lo + hi) / 2);
    if vals(mid) > target
        hi = mid - 1;
    else
        lo = mid;
    end
end
if vals(lo) > target
    idx = lo - 1;
else
    idx = lo;
end
end


function closeVideoWriter(writer)
if ~isempty(writer)
    close(writer);
end
end
