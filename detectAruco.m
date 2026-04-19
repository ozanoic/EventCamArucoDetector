function results = detectAruco(matFile, sensorSize, params)
%DETECTARUCO Multi-window sliding ArUco detection for event camera data.
%
%  results = detectAruco(matFile, sensorSize, params)
%
%  Inputs:
%    matFile    - path to .mat file containing 'events' (Nx4: x,y,pol,t)
%    sensorSize - [height, width]
%    params     - struct with fields:
%      .windowDurations_ms  - vector of lookback windows in ms (e.g. [5 10 20 50 100])
%      .tickStep_us         - tick step in microseconds (default 1000)
%      .showVis             - enable per-tick visualization (default false)
%      .useParallel         - true: use parfor if toolbox available (default);
%                             false: force sequential execution
%      .numCells            - marker grid size (default 8)
%      .codeSize            - inner code grid size (default 6)
%      .cellPx              - pixels per cell in unwarped image (default 20)
%      .blobParams          - struct with .minArea, .maxArea, .maxAspect
%
%  Output:
%    results - struct with fields:
%      .tNow_us             - timestamps (us) for each tick
%      .anyDetected         - 1 if any window detected a marker at that tick
%      .win_Xms             - detected marker ID per window (-1 = none)
%      .windowDurations_ms  - copy of window durations used
%      .detectionsPerWindow - total detections per window

%% ---- Default parameters ----
if ~isfield(params,'tickStep_us'),  params.tickStep_us  = 1000;  end
if ~isfield(params,'showVis'),      params.showVis      = false; end
if ~isfield(params,'useParallel'),  params.useParallel  = true;  end
if ~isfield(params,'numCells'),     params.numCells     = 8;     end
if ~isfield(params,'codeSize'),     params.codeSize     = 6;     end
if ~isfield(params,'cellPx'),       params.cellPx       = 20;    end

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

%% ---- Precompute event data (sorted by time) ----
evX = events(:,1) + 1;
evY = events(:,2) + 1;
evT = events(:,4);

[evT, sortIdx] = sort(evT);
evX = evX(sortIdx);
evY = evY(sortIdx);

tMin = evT(1);
tMax = evT(end);
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

blobParams = params.blobParams;
showVis = params.showVis;

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

    parfor tick = 1:numTicks
        tNow = tickTimes(tick);
        tickRow = -1 * ones(1, numWindows);

        for wi = 1:numWindows
            dt = windowDurations_us(wi);
            tFrom = tNow - dt;

            iStart = bsearchLeft(evT, tFrom);
            iEnd   = bsearchRight(evT, tNow);

            nEv = iEnd - iStart + 1;
            if nEv < 10, continue; end

            wX = evX(iStart:iEnd);
            wY = evY(iStart:iEnd);

            validMask = (wY >= 1) & (wY <= H) & (wX >= 1) & (wX <= W);
            wXv = wX(validMask);
            wYv = wY(validMask);
            if isempty(wXv), continue; end

            countImg = accumarray([wYv, wXv], 1, [H, W]);
            activeMask = countImg > 0;

            quads = detectQuadBlob(activeMask, blobParams);
            if isempty(quads), continue; end

            countU8 = uint8(countImg / max(countImg(:)) * 255);
            bestID = -1;

            for qi = 1:length(quads)
                corners = orderCornersForUnwarp_local(quads{qi});
                warpedImg = unwarpQuad_local(countU8, corners, markerCoords, sideSize);
                if isempty(warpedImg), continue; end

                mid = decodeMarker_local( ...
                    warpedImg, numCells, codeSize, cellPx, dictCodes, dictIDs);

                if mid >= 0
                    bestID = mid;
                    break;
                end
            end

            tickRow(wi) = bestID;
        end

        anyDet = double(any(tickRow >= 0));
        parResultTable(tick, :) = [tNow, anyDet, tickRow];
        parDetPerWindow(tick, :) = double(tickRow >= 0);

        send(dq, tick);
    end

    resultTable = parResultTable;
    detectionsPerWindow = sum(parDetPerWindow, 1);

else
    % ==================================================================
    %  SEQUENTIAL PATH
    % ==================================================================
    seqTicStart = tic;
    for tick = 1:numTicks
        tNow = tStart + (tick - 1) * tickStep_us;
        tickRow = -1 * ones(1, numWindows);

        for wi = 1:numWindows
            dt = windowDurations_us(wi);
            tFrom = tNow - dt;

            iStart = bsearchLeft(evT, tFrom);
            iEnd   = bsearchRight(evT, tNow);

            nEv = iEnd - iStart + 1;
            if nEv < 10, continue; end

            wX = evX(iStart:iEnd);
            wY = evY(iStart:iEnd);

            validMask = (wY >= 1) & (wY <= H) & (wX >= 1) & (wX <= W);
            wXv = wX(validMask);
            wYv = wY(validMask);
            if isempty(wXv), continue; end

            countImg = accumarray([wYv, wXv], 1, [H, W]);
            activeMask = countImg > 0;

            quads = detectQuadBlob(activeMask, blobParams);
            if isempty(quads), continue; end

            countU8 = uint8(countImg / max(countImg(:)) * 255);
            bestID = -1;

            for qi = 1:length(quads)
                corners = orderCornersForUnwarp_local(quads{qi});
                warpedImg = unwarpQuad_local(countU8, corners, markerCoords, sideSize);
                if isempty(warpedImg), continue; end

                mid = decodeMarker_local( ...
                    warpedImg, numCells, codeSize, cellPx, dictCodes, dictIDs);

                if mid >= 0
                    bestID = mid;
                    break;
                end
            end

            tickRow(wi) = bestID;
            if bestID >= 0
                detectionsPerWindow(wi) = detectionsPerWindow(wi) + 1;
            end
        end

        anyDet = double(any(tickRow >= 0));
        resultTable(tick, :) = [tNow, anyDet, tickRow];

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

            refWin = 6;
            dt = windowDurations_us(refWin);
            tFrom = tNow - dt;
            iS = bsearchLeft(evT, tFrom);
            iE = bsearchRight(evT, tNow);
            if iE >= iS
                rX = evX(iS:iE); rY = evY(iS:iE);
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
fprintf('\n--- Per-window breakdown ---\n');
fprintf('%-10s  %10s  %10s\n', 'Window', 'Detections', 'Rate');
fprintf('%-10s  %10s  %10s\n', '------', '----------', '----');
for wi = 1:numWindows
    fprintf('%-10s  %10d  %9.1f%%\n', ...
        sprintf('%dms', windowDurations_ms(wi)), ...
        detectionsPerWindow(wi), ...
        100 * detectionsPerWindow(wi) / max(numTicks, 1));
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
        warpedImg, numCells, codeSize, cellPx, dictCodes, dictIDs)

    markerID = -1;

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
        for inv = [0, 1]
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
                    idx = bsearchCode(dictCodes, code);
                    if idx > 0
                        markerID = dictIDs(idx);
                        return;
                    end
                end
            end
        end
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
