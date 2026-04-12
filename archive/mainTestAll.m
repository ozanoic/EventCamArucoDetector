%% Test All Border Detection Methods on Synthetic ESIM Data
%  Runs blob, line-pair, and OpenCV-contour methods side by side
%  on the same event packets to compare which finds border candidates.
%
%  Methods:
%    1. Blob (detectQuadBlob)    — fill + watershed + min-area rect
%    2. Lines (detectQuadLines)  — Hough parallel pairs → perpendicular quad
%    3. Contour (OpenCV-style)   — adaptive threshold → contours → approxPoly

clear; close all; clc;

%% ---- Settings ----
matFile    = '../Data/Synthetic/MovingCam/moving_events/moving_events.mat';
sensorSize = [240, 320];   % [H, W]
packetDt   = 20000;         % 10 ms

%% ---- Load data ----
fprintf('Loading %s...\n', matFile);
tmp = load(matFile, 'events');
events = tmp.events;
numEvents = size(events, 1);
H = sensorSize(1);
W = sensorSize(2);
fprintf('Loaded %d events  |  sensor %dx%d\n', numEvents, W, H);

%% ---- Parameters for each method ----
% Blob params
blobParams.minArea   = 500;
blobParams.maxArea   = H*W*0.5;
blobParams.maxAspect = 3.0;

% Line params
lineParams.minLen      = 30;
lineParams.maxLen      = 300;
lineParams.angleTol    = 15;
lineParams.lengthRatio = 0.5;
lineParams.gapRatio    = [0.2 4.0];
lineParams.numPeaks    = 60;

% Contour params (OpenCV-style)
adaptiveWindows = [7, 15, 25, 41];
adaptiveC       = 5;
minPerimeter    = 60;
maxPerimeter    = 2*(H+W);
polyEpsFrac     = 0.03;
contourMinArea  = 500;
contourMaxArea  = H*W*0.5;
contourMinSide  = 15;
contourMaxAspect= 4.0;

%% ---- Group into packets ----
tAll = events(:, 4);
tMin = min(tAll);
tMax = max(tAll);
packetEdges = tMin : packetDt : (tMax + packetDt);
numPackets = length(packetEdges) - 1;
[~, ~, packetIdx] = histcounts(tAll, packetEdges);
fprintf('Packets: %d  (%.0f us each)\n', numPackets, packetDt);

%% ---- Process ----
hFig = figure('Name', 'All Methods Comparison', 'Position', [30 30 1800 900]);

for p = 1:numPackets
    idx = (packetIdx == p);
    if sum(idx) < 10, continue; end

    pktEvents = events(idx, :);
    px  = pktEvents(:,1) + 1;   % 1-indexed
    py  = pktEvents(:,2) + 1;
    pol = pktEvents(:,3);        % 0=OFF, 1=ON

    % ---- Build event images ----
    % Combined binary mask (all events)
    activeMask = false(H, W);
    onImg  = zeros(H, W);
    offImg = zeros(H, W);
    for e = 1:size(pktEvents, 1)
        r = py(e); c = px(e);
        if r < 1 || r > H || c < 1 || c > W, continue; end
        activeMask(r, c) = true;
        if pol(e) == 1
            onImg(r, c) = onImg(r, c) + 1;
        else
            offImg(r, c) = offImg(r, c) + 1;
        end
    end

    % Count image for intensity-based thresholding
    countImg = onImg + offImg;
    % Normalize to 0-255 for adaptive threshold
    if max(countImg(:)) > 0
        grayImg = uint8(countImg / max(countImg(:)) * 255);
    else
        continue;
    end

    % ---- Method 1: Blob ----
    quadsBlob = detectQuadBlob(activeMask, blobParams);

    % ---- Method 2: Lines ----
    quadsLines = detectQuadLines(activeMask, lineParams);

    % ---- Method 3: Contour (OpenCV-style) ----
    quadsContour = {};
    for wi = 1:length(adaptiveWindows)
        ws = adaptiveWindows(wi);
        % Adaptive threshold
        localMean = imfilter(double(grayImg), ones(ws)/ws^2, 'replicate');
        bw = double(grayImg) > (localMean - adaptiveC);
        bw = ~bw;  % invert: we want dark regions (marker border) as foreground

        % Also try on the binary active mask directly
        if wi == 1
            bwDirect = activeMask;
            % Morphological close to fill small gaps
            bwDirect = imclose(bwDirect, strel('disk', 2));
            bwDirect = imfill(bwDirect, 'holes');
        end

        % Find contours via boundary tracing
        [B, L] = bwboundaries(bw, 'noholes');
        for bi = 1:length(B)
            boundary = B{bi};
            perim = size(boundary, 1);
            if perim < minPerimeter || perim > maxPerimeter, continue; end

            % Approximate polygon (Douglas-Peucker)
            eps_val = polyEpsFrac * perim;
            poly = approxPoly(boundary(:,[2 1]), eps_val);  % boundary is [row,col] → [x,y]

            if size(poly, 1) ~= 4, continue; end

            % Convexity check
            if ~isConvex(poly), continue; end

            % Area check
            area = polyarea(poly(:,1), poly(:,2));
            if area < contourMinArea || area > contourMaxArea, continue; end

            % Side length check
            sides = zeros(4,1);
            for si = 1:4
                sides(si) = norm(poly(si,:) - poly(mod(si,4)+1,:));
            end
            if min(sides) < contourMinSide, continue; end
            if max(sides)/min(sides) > contourMaxAspect, continue; end

            quadsContour{end+1} = poly; %#ok<AGROW>
        end
    end

    % Also find contours on the direct active mask (morphologically closed)
    [B2, ~] = bwboundaries(bwDirect, 'noholes');
    for bi = 1:length(B2)
        boundary = B2{bi};
        perim = size(boundary, 1);
        if perim < minPerimeter || perim > maxPerimeter, continue; end
        eps_val = polyEpsFrac * perim;
        poly = approxPoly(boundary(:,[2 1]), eps_val);
        if size(poly, 1) ~= 4, continue; end
        if ~isConvex(poly), continue; end
        area = polyarea(poly(:,1), poly(:,2));
        if area < contourMinArea || area > contourMaxArea, continue; end
        sides = zeros(4,1);
        for si = 1:4
            sides(si) = norm(poly(si,:) - poly(mod(si,4)+1,:));
        end
        if min(sides) < contourMinSide, continue; end
        if max(sides)/min(sides) > contourMaxAspect, continue; end
        quadsContour{end+1} = poly; %#ok<AGROW>
    end

    % Deduplicate contour quads
    quadsContour = deduplicateQuads(quadsContour);

    % ---- Visualize ----
    figure(hFig); clf;

    % Row 1: Raw data
    % (1) Combined event image
    subplot(2,4,1);
    combRGB = zeros(H, W, 3);
    combRGB(:,:,1) = min(onImg / max(onImg(:)+eps), 1);
    combRGB(:,:,3) = min(offImg / max(offImg(:)+eps), 1);
    imshow(combRGB); axis image;
    title(sprintf('Events (R=ON B=OFF) | %d evts', sum(idx)));

    % (2) Active mask
    subplot(2,4,2);
    imshow(activeMask); axis image;
    title('Active mask');

    % (3) Gray image (count-based)
    subplot(2,4,3);
    imshow(grayImg); axis image;
    title('Count image (gray)');

    % (4) Morphologically closed
    subplot(2,4,4);
    imshow(bwDirect); axis image;
    title('Closed + filled mask');

    % Row 2: Detection results
    % (5) Blob method
    subplot(2,4,5);
    imshow(activeMask); hold on;
    for qi = 1:length(quadsBlob)
        c = quadsBlob{qi};
        plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], 'g-', 'LineWidth', 2);
        plot(c(:,1), c(:,2), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    end
    title(sprintf('Blob: %d quads', length(quadsBlob)));
    axis image; hold off;

    % (6) Lines method
    subplot(2,4,6);
    imshow(activeMask); hold on;
    for qi = 1:length(quadsLines)
        c = quadsLines{qi};
        plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], 'y-', 'LineWidth', 2);
        plot(c(:,1), c(:,2), 'yo', 'MarkerSize', 6, 'MarkerFaceColor', 'y');
    end
    title(sprintf('Lines: %d quads', length(quadsLines)));
    axis image; hold off;

    % (7) Contour method
    subplot(2,4,7);
    imshow(activeMask); hold on;
    for qi = 1:length(quadsContour)
        c = quadsContour{qi};
        plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], 'c-', 'LineWidth', 2);
        plot(c(:,1), c(:,2), 'co', 'MarkerSize', 6, 'MarkerFaceColor', 'c');
    end
    title(sprintf('Contour: %d quads', length(quadsContour)));
    axis image; hold off;

    % (8) All overlaid
    subplot(2,4,8);
    imshow(combRGB); hold on;
    for qi = 1:length(quadsBlob)
        c = quadsBlob{qi};
        plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], 'g-', 'LineWidth', 2);
    end
    for qi = 1:length(quadsLines)
        c = quadsLines{qi};
        plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], 'y-', 'LineWidth', 1.5);
    end
    for qi = 1:length(quadsContour)
        c = quadsContour{qi};
        plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], 'c-', 'LineWidth', 1.5);
    end
    title(sprintf('All: B=%d L=%d C=%d', ...
        length(quadsBlob), length(quadsLines), length(quadsContour)));
    axis image; hold off;

    tStart = packetEdges(p);
    sgtitle(sprintf('Packet %d/%d  |  t=%.3f s  |  Blob:%d  Lines:%d  Contour:%d', ...
        p, numPackets, tStart/1e6, ...
        length(quadsBlob), length(quadsLines), length(quadsContour)));
    drawnow;

    % Console log
    if mod(p, 50) == 0 || ~isempty(quadsBlob) || ~isempty(quadsLines) || ~isempty(quadsContour)
        fprintf('[pkt %4d/%d]  Blob:%d  Lines:%d  Contour:%d\n', ...
            p, numPackets, length(quadsBlob), length(quadsLines), length(quadsContour));
    end
end

fprintf('\nDone.\n');


%% =========================================================================
%                       LOCAL FUNCTIONS
%% =========================================================================

function poly = approxPoly(pts, epsilon)
%APPROXPOLY Douglas-Peucker polygon approximation (like cv::approxPolyDP)
%  pts: Nx2 [x,y], epsilon: max distance threshold
%  Returns simplified polygon vertices.
    if size(pts,1) < 3, poly = pts; return; end

    % Close the polygon if not already closed
    if norm(pts(1,:) - pts(end,:)) > 1
        pts = [pts; pts(1,:)];
    end

    % Douglas-Peucker
    keep = dpSimplify(pts, epsilon);
    poly = pts(keep, :);

    % Remove last point if it duplicates the first (closed polygon)
    if size(poly,1) > 1 && norm(poly(1,:) - poly(end,:)) < 1
        poly = poly(1:end-1, :);
    end
end

function keep = dpSimplify(pts, epsilon)
    n = size(pts, 1);
    if n <= 2
        keep = true(n, 1);
        return;
    end

    % Find the point with maximum distance from line(first, last)
    dmax = 0;
    imax = 1;
    p1 = pts(1,:);
    p2 = pts(end,:);
    for i = 2:n-1
        d = pointLineDistance(pts(i,:), p1, p2);
        if d > dmax
            dmax = d;
            imax = i;
        end
    end

    keep = false(n, 1);
    if dmax > epsilon
        k1 = dpSimplify(pts(1:imax,:), epsilon);
        k2 = dpSimplify(pts(imax:end,:), epsilon);
        keep(1:imax) = k1;
        keep(imax:end) = keep(imax:end) | k2;
    else
        keep(1) = true;
        keep(end) = true;
    end
end

function d = pointLineDistance(p, a, b)
    ab = b - a;
    ap = p - a;
    len2 = dot(ab, ab);
    if len2 < 1e-10
        d = norm(ap);
    else
        d = abs(ab(1)*ap(2) - ab(2)*ap(1)) / sqrt(len2);
    end
end

function ok = isConvex(poly)
%ISCONVEX Check if polygon vertices form a convex shape
    n = size(poly, 1);
    if n < 3, ok = false; return; end
    signs = zeros(n, 1);
    for i = 1:n
        j = mod(i, n) + 1;
        k = mod(j, n) + 1;
        v1 = poly(j,:) - poly(i,:);
        v2 = poly(k,:) - poly(j,:);
        signs(i) = v1(1)*v2(2) - v1(2)*v2(1);
    end
    ok = all(signs > 0) || all(signs < 0);
end

function merged = deduplicateQuads(quadsIn)
    if isempty(quadsIn), merged = {}; return; end
    centers = zeros(length(quadsIn), 2);
    sizes   = zeros(length(quadsIn), 1);
    for k = 1:length(quadsIn)
        c = quadsIn{k};
        centers(k,:) = mean(c, 1);
        s = zeros(size(c,1), 1);
        for j = 1:size(c,1)
            s(j) = norm(c(j,:) - c(mod(j,size(c,1))+1,:));
        end
        sizes(k) = mean(s);
    end
    keep = true(length(quadsIn), 1);
    for k = 1:length(quadsIn)
        if ~keep(k), continue; end
        for m = k+1:length(quadsIn)
            if ~keep(m), continue; end
            if norm(centers(k,:)-centers(m,:)) < 0.5*min(sizes(k),sizes(m))
                keep(m) = false;
            end
        end
    end
    merged = quadsIn(keep);
end
