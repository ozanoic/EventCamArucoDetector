%% OpenCV-Style ArUco Detection for Event Camera Data
%  Adapts the classic OpenCV ArUco pipeline to event camera data:
%    1. Accumulate events into a binary image (SAE active mask)
%    2. Adaptive thresholding at multiple levels
%    3. Find contours → approximate to quads (4-corner polygons)
%    4. Filter quads by area, convexity, aspect ratio
%    5. Perspective unwarp each quad to canonical square
%    6. Decode bits: Otsu threshold on cells → binary grid
%    7. Dictionary lookup (ARUCO_MIP_36h12, 250 markers, 4 rotations)
%
%  Works with both Sarmadi's DVS128 data and user's 640×480 data.

clear; close all; clc;

%% ---- User settings ---------------------------------------------------------
visualizeAll = 1;   % true: visualize every packet | false: only on detection

% Choose dataset: 'sarmadi' or 'user'
datasetMode = 'sarmadi';

%% ---- Load events -----------------------------------------------------------
switch datasetMode
    case 'sarmadi'
        binFile = '../Data/Sarmadi/side2side/side2side/packets.bin';
        fprintf('Loading Sarmadi .bin data...\n');
        events = convertSarmadiBin(binFile);
        sensorSize = [128, 128];
    case 'user'
        matFile = '../Data/zoom_1000_100_10ms_g1/events_out_start_cropped.mat';
        fprintf('Loading user .mat data...\n');
        events = double(loadEvents(matFile));
        sensorSize = [480, 640];
end
numEvents = size(events, 1);
H = sensorSize(1);
W = sensorSize(2);

fprintf('Loaded %d events  |  sensor %dx%d\n', numEvents, W, H);

%% ---- Parameters ------------------------------------------------------------
packetDt = 10000;           % 10 ms packets

% Quad detection
minPerimeter = 40;          % minimum contour perimeter (pixels)
maxPerimeter = 2*(H+W);    % maximum contour perimeter
polyEpsFrac  = 0.03;       % approxPolyDP epsilon as fraction of perimeter
minArea      = 200;         % minimum quad area
maxArea      = H*W*0.5;    % maximum quad area (half the image)
minSideLen   = 10;          % minimum side length of quad
maxAspect    = 4.0;         % maximum aspect ratio of bounding box

% Adaptive threshold window sizes (multiple scales)
adaptiveWindows = [7, 13, 23, 41];
adaptiveC       = 5;        % constant subtracted from mean

% Marker decoding
numCells     = 7;           % 7×7 grid (5×5 inner + 1-cell border)
codeSize     = 5;           % inner grid = 5×5
cellPx       = 20;          % pixels per cell in unwarped image
sideSize     = numCells * cellPx;  % 140 pixels
borderBlackThresh = 0.6;    % fraction of border cells that must be black

% Standard marker corners (0-indexed)
markerCoords = [0 0; sideSize-1 0; sideSize-1 sideSize-1; 0 sideSize-1];

% Build dictionary
dictionary = buildDictionary_4x4_50();
fprintf('Dictionary loaded: %d markers\n', dictionary.Count);

%% ---- Group events into packets by time ------------------------------------
tAll = events(:, 4);
tMin = min(tAll);
tMax = max(tAll);
packetEdges = tMin : packetDt : (tMax + packetDt);
numPackets = length(packetEdges) - 1;
[~, ~, packetIdx] = histcounts(tAll, packetEdges);

fprintf('Total events: %d  |  Packets: %d  (%.0f us each)\n', ...
    numEvents, numPackets, packetDt);

%% ---- Prepare figure --------------------------------------------------------
hFig = figure('Name','OpenCV-Style ArUco Detection', ...
              'Position',[50 50 1600 800]);

totalDetections = 0;

%% ---- Open output file ------------------------------------------------------
logFid = fopen('opencv_output.txt', 'w');
fprintf(logFid, 'OpenCV-Style ArUco Detection Log\n');
fprintf(logFid, 'Total events: %d  |  Packets: %d  |  packetDt: %d us\n', ...
    numEvents, numPackets, packetDt);
fprintf(logFid, '----------------------------------------------------------------------\n');
fprintf(logFid, '%6s  %6s  %5s  %5s  %4s  %s\n', ...
    'Packet', 'Events', 'Conts', 'Quads', 'Det', 'MarkerIDs');
fprintf(logFid, '----------------------------------------------------------------------\n');

%% ---- Process each packet ---------------------------------------------------
fprintf('\nProcessing packets...\n');
tic;

for p = 1:numPackets
    % ---- Get events in this packet ----
    idx = (packetIdx == p);
    if sum(idx) < 10, continue; end

    pktEvents = events(idx, :);
    nEvt = size(pktEvents, 1);

    px  = pktEvents(:,1) + 1;      % 0→1 indexed
    py  = pktEvents(:,2) + 1;
    pol = pktEvents(:,3);

    % ---- Build event image ----
    % Accumulate event count per pixel (both polarities)
    evtImg = zeros(H, W);
    for e = 1:nEvt
        r = py(e);  c = px(e);
        if r < 1 || r > H || c < 1 || c > W, continue; end
        evtImg(r, c) = evtImg(r, c) + 1;
    end

    % Also build separate polarity images for decoding
    onImg  = zeros(H, W);
    offImg = zeros(H, W);
    for e = 1:nEvt
        r = py(e);  c = px(e);
        if r < 1 || r > H || c < 1 || c > W, continue; end
        if pol(e) > 0
            onImg(r, c) = onImg(r, c) + 1;
        else
            offImg(r, c) = offImg(r, c) + 1;
        end
    end

    % ---- Step 1-3: Adaptive threshold + contour + quad detection ----
    % Normalize event image to 0-255
    maxVal = max(evtImg(:));
    if maxVal == 0, continue; end
    img8 = uint8(evtImg / maxVal * 255);

    allQuads = {};
    totalContours = 0;

    for wi = 1:length(adaptiveWindows)
        ws = adaptiveWindows(wi);

        % Adaptive threshold (mean-based)
        bw = imbinarize(img8, 'adaptive', 'ForegroundPolarity', 'dark', ...
            'Sensitivity', 0.5);
        % Also try with inverted polarity
        bwInv = imbinarize(img8, 'adaptive', 'ForegroundPolarity', 'bright', ...
            'Sensitivity', 0.5);

        % Try both
        for bwIdx = 1:2
            if bwIdx == 1
                bwCurr = bw;
            else
                bwCurr = bwInv;
            end

            % Find contours
            [B, ~] = bwboundaries(bwCurr, 'noholes');
            totalContours = totalContours + length(B);

            for bi = 1:length(B)
                contour = B{bi};  % Nx2 [row, col]
                perim = size(contour, 1);

                % Filter by perimeter
                if perim < minPerimeter || perim > maxPerimeter, continue; end

                % Approximate polygon (Douglas-Peucker)
                epsilon = polyEpsFrac * perim;
                poly = approxPoly(contour, epsilon);

                % Keep only quadrilaterals
                if size(poly, 1) ~= 4, continue; end

                % Convert to [x, y] (col, row) format
                corners = [poly(:,2), poly(:,1)];

                % Filter: area
                a = polyarea(corners(:,1), corners(:,2));
                if a < minArea || a > maxArea, continue; end

                % Filter: convexity
                if ~isConvex(corners), continue; end

                % Filter: minimum side length
                sides = zeros(4,1);
                for s = 1:4
                    sides(s) = norm(corners(s,:) - corners(mod(s,4)+1,:));
                end
                if min(sides) < minSideLen, continue; end

                % Filter: aspect ratio
                if max(sides)/min(sides) > maxAspect, continue; end

                % Order corners consistently (TL, TR, BR, BL)
                corners = orderCorners(corners);

                allQuads{end+1} = corners; %#ok<AGROW>
            end
        end
    end

    % Deduplicate overlapping quads
    if length(allQuads) > 1
        allQuads = deduplicateQuads(allQuads);
    end

    % ---- Step 5-7: Unwarp, decode, dictionary lookup ----
    detectedMarkers = [];

    for qi = 1:length(allQuads)
        corners = allQuads{qi};

        % Perspective unwarp to canonical square
        warpedImg = unwarpQuad(img8, corners, markerCoords, sideSize);
        if isempty(warpedImg), continue; end

        % Decode the marker
        [grid, codes4rot] = decodeMarker(warpedImg, numCells, codeSize, ...
            cellPx, borderBlackThresh);
        if isempty(grid), continue; end

        % Dictionary lookup (try all 4 rotations)
        markerIdx = -1;
        orientation = -1;
        for rot = 1:4
            key = codes4rot(rot);
            if dictionary.isKey(key)
                markerIdx = dictionary(key);
                orientation = rot;
                break;
            end
        end

        if markerIdx >= 0
            det.id          = markerIdx;
            det.corners     = corners;
            det.grid        = grid;
            det.warpedImg   = warpedImg;
            det.orientation = orientation;
            detectedMarkers = [detectedMarkers, det]; %#ok<AGROW>
            totalDetections = totalDetections + 1;
        end
    end

    nDet = length(detectedMarkers);
    nQuads = length(allQuads);

    % ---- Write to log file ----
    idStr = '';
    if nDet > 0
        ids = arrayfun(@(d) d.id, detectedMarkers);
        idStr = strjoin(arrayfun(@(x) sprintf('%d', x), ids, ...
            'UniformOutput', false), ',');
    end
    fprintf(logFid, '%6d  %6d  %5d  %5d  %4d  %s\n', ...
        p, nEvt, totalContours, nQuads, nDet, idStr);

    % ---- Console log ----
    if mod(p, 100) == 0 || nDet > 0
        fprintf('[pkt %4d/%d]  events:%d  contours:%d  quads:%d  detected:%d', ...
            p, numPackets, nEvt, totalContours, nQuads, nDet);
        if nDet > 0
            for di = 1:nDet
                fprintf('  ID=%d', detectedMarkers(di).id);
            end
        end
        fprintf('\n');
    end

    % ---- Visualize (conditional) ----
    if visualizeAll || nDet > 0
        figure(hFig); clf;

        % (a) Event image
        subplot(2,3,1);
        imagesc(evtImg); colormap(gca,'hot'); colorbar; axis image;
        title(sprintf('Event image | pkt %d | %d evts', p, nEvt));

        % (b) Thresholded + quads
        subplot(2,3,2);
        imshow(img8); hold on;
        % Draw all quads in yellow
        for qi2 = 1:nQuads
            cc = allQuads{qi2};
            plot([cc(:,1); cc(1,1)], [cc(:,2); cc(1,2)], 'y-', 'LineWidth', 1.5);
        end
        % Draw detected markers in green
        for di = 1:nDet
            cc = detectedMarkers(di).corners;
            plot([cc(:,1); cc(1,1)], [cc(:,2); cc(1,2)], 'g-', 'LineWidth', 2.5);
            plot(cc(:,1), cc(:,2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
            text(mean(cc(:,1)), mean(cc(:,2)), ...
                sprintf('ID:%d', detectedMarkers(di).id), ...
                'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', 'k');
        end
        title(sprintf('Quads: %d | Detected: %d', nQuads, nDet));
        axis image; hold off;

        % (c) ON/OFF overlay
        subplot(2,3,3);
        combRGB = zeros(H, W, 3);
        if max(onImg(:)) > 0
            combRGB(:,:,1) = onImg / max(onImg(:));
        end
        if max(offImg(:)) > 0
            combRGB(:,:,3) = offImg / max(offImg(:));
        end
        imshow(combRGB * 3); axis image;
        title('ON (red) / OFF (blue)');

        % (d-f) Show detected marker details
        if nDet > 0
            det1 = detectedMarkers(1);

            subplot(2,3,4);
            imshow(det1.warpedImg); axis image;
            title(sprintf('Unwarped (ID:%d)', det1.id));

            subplot(2,3,5);
            imagesc(det1.grid); colormap(gca, 'gray'); axis image;
            title('Decoded grid');
            hold on;
            for gi = 0:numCells
                plot([gi+0.5, gi+0.5], [0.5, numCells+0.5], 'r-', 'LineWidth', 0.5);
                plot([0.5, numCells+0.5], [gi+0.5, gi+0.5], 'r-', 'LineWidth', 0.5);
            end
            hold off;

            subplot(2,3,6);
            markerImg = zeros(numCells);
            markerImg(2:end-1, 2:end-1) = det1.grid(2:end-1, 2:end-1);
            imagesc(1 - markerImg); colormap(gca, 'gray'); axis image;
            title(sprintf('Marker ID: %d', det1.id));
            hold on;
            for gi = 0:numCells
                plot([gi+0.5, gi+0.5], [0.5, numCells+0.5], 'r-', 'LineWidth', 0.5);
                plot([0.5, numCells+0.5], [gi+0.5, gi+0.5], 'r-', 'LineWidth', 0.5);
            end
            hold off;
        else
            subplot(2,3,4); cla; title('No detection');
            subplot(2,3,5); cla; title('No detection');
            subplot(2,3,6); cla; title('No detection');
        end

        sgtitle(sprintf('OpenCV-Style | Packet %d/%d | %d events | Total det: %d', ...
            p, numPackets, nEvt, totalDetections));
        drawnow;
    end
end

elapsed = toc;
fprintf('\n=== Done ===\n');
fprintf('Total packets     : %d\n', numPackets);
fprintf('Total events      : %d\n', numEvents);
fprintf('Total detections  : %d\n', totalDetections);
fprintf('Elapsed time      : %.1f s\n', elapsed);

% Close output file
fprintf(logFid, '----------------------------------------------------------------------\n');
fprintf(logFid, 'Total packets    : %d\n', numPackets);
fprintf(logFid, 'Total events     : %d\n', numEvents);
fprintf(logFid, 'Total detections : %d\n', totalDetections);
fprintf(logFid, 'Elapsed time     : %.1f s\n', elapsed);
fclose(logFid);
fprintf('Debug log saved to: opencv_output.txt\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%  =========================================================================

%% ---- approxPoly: Douglas-Peucker polygon approximation --------------------
function poly = approxPoly(contour, epsilon)
%  Simplify a contour (Nx2 [row,col]) to fewer vertices.
%  This is MATLAB's equivalent of cv::approxPolyDP.
    N = size(contour, 1);
    if N < 4
        poly = contour;
        return;
    end

    % Use reducepoly if available, otherwise manual Douglas-Peucker
    try
        % reducepoly expects [x y] and tolerance as fraction
        pts = [contour(:,2), contour(:,1)];  % [col, row] → [x, y]
        totalLen = 0;
        for i = 1:size(pts,1)
            j = mod(i, size(pts,1)) + 1;
            totalLen = totalLen + norm(pts(i,:) - pts(j,:));
        end
        reduced = reducepoly(pts, epsilon / totalLen);
        % Remove last point if it duplicates first (closed polygon)
        if size(reduced,1) > 1 && norm(reduced(1,:) - reduced(end,:)) < 1
            reduced = reduced(1:end-1, :);
        end
        poly = [reduced(:,2), reduced(:,1)];  % back to [row, col]
    catch
        % Fallback: manual Douglas-Peucker
        poly = douglasPeucker(contour, epsilon);
    end
end

%% ---- douglasPeucker: manual implementation --------------------------------
function result = douglasPeucker(points, epsilon)
    N = size(points, 1);
    if N <= 2
        result = points;
        return;
    end

    % Find the point with maximum distance from the line (first→last)
    dmax = 0;
    idx = 1;
    p1 = points(1,:);
    p2 = points(end,:);
    for i = 2:N-1
        d = pointToLineDistance(points(i,:), p1, p2);
        if d > dmax
            dmax = d;
            idx = i;
        end
    end

    % If max distance > epsilon, recursively simplify
    if dmax > epsilon
        r1 = douglasPeucker(points(1:idx,:), epsilon);
        r2 = douglasPeucker(points(idx:end,:), epsilon);
        result = [r1(1:end-1,:); r2];
    else
        result = [points(1,:); points(end,:)];
    end
end

%% ---- pointToLineDistance ---------------------------------------------------
function d = pointToLineDistance(p, a, b)
    ab = b - a;
    ap = p - a;
    len2 = dot(ab, ab);
    if len2 < 1e-10
        d = norm(ap);
    else
        d = abs(ab(1)*ap(2) - ab(2)*ap(1)) / sqrt(len2);
    end
end

%% ---- isConvex: check if polygon is convex ---------------------------------
function ok = isConvex(corners)
    n = size(corners, 1);
    ok = true;
    sign = 0;
    for i = 1:n
        j = mod(i, n) + 1;
        k = mod(j, n) + 1;
        d1 = corners(j,:) - corners(i,:);
        d2 = corners(k,:) - corners(j,:);
        cross = d1(1)*d2(2) - d1(2)*d2(1);
        if sign == 0
            sign = (cross > 0) - (cross < 0);
        elseif sign * cross < 0
            ok = false;
            return;
        end
    end
end

%% ---- orderCorners: order 4 points as TL, TR, BR, BL ----------------------
function c = orderCorners(pts)
    % Sort by sum (x+y) → TL has smallest, BR has largest
    % Sort by diff (y-x) → TR has smallest, BL has largest
    center = mean(pts, 1);
    angles = atan2(pts(:,2) - center(2), pts(:,1) - center(1));
    [~, ord] = sort(angles);
    pts = pts(ord, :);

    % Now sorted CCW. Find TL (smallest x+y)
    sums = pts(:,1) + pts(:,2);
    [~, tlIdx] = min(sums);
    % Rotate so TL is first
    pts = circshift(pts, -(tlIdx-1), 1);

    % Check winding: should be TL→TR→BR→BL (clockwise)
    % If cross product says CCW, reverse the last 3
    d1 = pts(2,:) - pts(1,:);
    d2 = pts(4,:) - pts(1,:);
    cross = d1(1)*d2(2) - d1(2)*d2(1);
    if cross > 0
        % Currently CCW (TL→BL→BR→TR), flip to CW
        pts = pts([1, 4, 3, 2], :);
    end
    c = pts;  % TL, TR, BR, BL
end

%% ---- unwarpQuad: perspective transform to canonical square ----------------
function warped = unwarpQuad(img, corners, markerCoords, sideSize)
    warped = [];
    try
        % fitgeotrans: moving (image corners) → fixed (marker square)
        tform = fitgeotrans(corners, markerCoords + 1, 'projective');
        outputView = imref2d([sideSize, sideSize]);
        warped = imwarp(img, tform, 'OutputView', outputView, 'Interp', 'bilinear');
    catch
        warped = [];
    end
end

%% ---- decodeMarker: extract bit grid from unwarped image -------------------
function [grid, codes4rot] = decodeMarker(warpedImg, numCells, codeSize, cellPx, borderBlackThresh)
    grid = [];
    codes4rot = [];

    wImg = double(warpedImg);
    sideSize = numCells * cellPx;

    % Otsu threshold on the whole unwarped image
    level = graythresh(uint8(wImg));
    thresh = level * 255;

    % Build cell grid: compute mean intensity per cell
    cellMeans = zeros(numCells);
    for r = 1:numCells
        for c = 1:numCells
            rStart = (r-1)*cellPx + 1;
            rEnd   = r*cellPx;
            cStart = (c-1)*cellPx + 1;
            cEnd   = c*cellPx;
            % Use inner 60% of each cell to avoid border effects
            margin = round(cellPx * 0.2);
            rStart = rStart + margin;
            rEnd   = rEnd - margin;
            cStart = cStart + margin;
            cEnd   = cEnd - margin;
            rStart = max(1, rStart);  rEnd = min(sideSize, rEnd);
            cStart = max(1, cStart);  cEnd = min(sideSize, cEnd);
            cellMeans(r, c) = mean(wImg(rStart:rEnd, cStart:cEnd), 'all');
        end
    end

    % Binary grid: 1 = white (marker bit 1), 0 = black (marker bit 0)
    grid = double(cellMeans > thresh);

    % Check border: all border cells should be black (0)
    borderCells = [grid(1,:), grid(end,:), grid(2:end-1,1)', grid(2:end-1,end)'];
    nBorderBlack = sum(borderCells == 0);
    nBorderTotal = length(borderCells);
    if nBorderBlack < borderBlackThresh * nBorderTotal
        grid = [];
        codes4rot = [];
        return;
    end

    % Extract inner code (codeSize × codeSize)
    innerGrid = grid(2:end-1, 2:end-1);

    % Generate 4 rotation codes
    codes4rot = zeros(1, 4, 'uint64');
    for rot = 0:3
        rotGrid = rot90(innerGrid, rot);
        code = uint64(0);
        for r = 1:codeSize
            for c = 1:codeSize
                code = bitshift(code, 1);
                code = code + uint64(rotGrid(r, c));
            end
        end
        codes4rot(rot+1) = code;
    end
end

%% ---- deduplicateQuads: merge overlapping quad detections ------------------
function merged = deduplicateQuads(allQ)
    n = length(allQ);
    centers = zeros(n, 2);
    sizes = zeros(n, 1);
    for k = 1:n
        c = allQ{k};
        centers(k,:) = mean(c, 1);
        s = zeros(4,1);
        for j = 1:4
            s(j) = norm(c(j,:) - c(mod(j,4)+1,:));
        end
        sizes(k) = mean(s);
    end
    keep = true(n, 1);
    for k = 1:n
        if ~keep(k), continue; end
        for m = k+1:n
            if ~keep(m), continue; end
            if norm(centers(k,:) - centers(m,:)) < 0.5 * min(sizes(k), sizes(m))
                keep(m) = false;
            end
        end
    end
    merged = allQ(keep);
end

%% ---- buildDictionary_4x4_50: ARUCO_MIP_36h12 dictionary ------------------
%  250 markers, 36-bit codes (6×6 inner grid), 4 rotations checked.
%  Note: This is the same dictionary as Sarmadi's. For standard OpenCV
%  DICT_4X4_50 (25-bit, 5×5 inner), a different dictionary would be needed.
%  Using 36h12 for compatibility with Sarmadi's data.
function dict = buildDictionary_4x4_50()
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
    dict = containers.Map('KeyType','uint64','ValueType','int32');
    for i = 1:length(codes)
        dict(codes(i)) = int32(i - 1);   % 0-indexed marker ID
    end
end
