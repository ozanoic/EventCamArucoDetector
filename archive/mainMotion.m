%% Motion-Aware ArUco Detection — Step 1 & 2
%  Step 1: Find quad border candidates from event blobs (polarity-agnostic)
%  Step 2: Estimate marker motion from polarity distribution along edges
%
%  Key insight: ArUco markers have black borders on white background.
%  When moving, leading edges produce OFF events (white→black),
%  trailing edges produce ON events (black→white).
%  The polarity pattern reveals the motion direction.

clear; close all; clc;

%% ---- User settings ---------------------------------------------------------
visualizeAll = true;    % true: show every packet | false: only when quads found

% Choose dataset
datasetMode = 'user';   % 'sarmadi' or 'user'

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
packetDt = 10000;           % 10 ms per packet

% Morphological preprocessing
dilateRadius = 1;           % dilation radius to connect nearby events
closeRadius  = 1;           % closing radius to fill small gaps

% Quad detection
minPerimeter = 30;          % minimum contour perimeter
polyEpsFrac  = 0.04;       % Douglas-Peucker epsilon (fraction of perimeter)
minArea      = 100;         % minimum quad area (px^2)
maxArea      = H*W*0.5;    % maximum quad area
minSideLen   = 8;           % minimum side length
maxAspect    = 5.0;         % max side ratio

% Motion estimation
motionBandWidth = 5;        % how many pixels on each side of edge to sample

%% ---- Group events into packets ---------------------------------------------
tAll = events(:, 4);
tMin = min(tAll);
tMax = max(tAll);
packetEdges = tMin : packetDt : (tMax + packetDt);
numPackets = length(packetEdges) - 1;
[~, ~, packetIdx] = histcounts(tAll, packetEdges);

fprintf('Packets: %d  (%.0f us each)\n', numPackets, packetDt);

%% ---- Open log file ---------------------------------------------------------
logFid = fopen('motion_output.txt', 'w');
fprintf(logFid, 'Motion-Aware ArUco Detection — Step 1 & 2\n');
fprintf(logFid, 'Packets: %d | packetDt: %d us\n', numPackets, packetDt);
fprintf(logFid, '-----------------------------------------------------------\n');
fprintf(logFid, '%6s  %6s  %5s  %s\n', 'Packet', 'Events', 'Quads', 'MotionVectors');
fprintf(logFid, '-----------------------------------------------------------\n');

%% ---- Prepare figure --------------------------------------------------------
hFig = figure('Name','Motion-Aware ArUco — Step 1 & 2', ...
              'Position',[50 50 1600 800]);

%% ---- Process each packet ---------------------------------------------------
fprintf('\nProcessing packets...\n');
tic;
totalQuads = 0;

for p = 1:numPackets
    % ---- Get events in this packet ----
    mask = (packetIdx == p);
    if sum(mask) < 10, continue; end

    pktEvents = events(mask, :);
    nEvt = size(pktEvents, 1);

    px  = pktEvents(:,1) + 1;   % 0→1 indexed
    py  = pktEvents(:,2) + 1;
    pol = pktEvents(:,3);        % -1 or +1

    % ---- Build images ----
    % Combined event count (polarity-agnostic) for blob detection
    evtMask = false(H, W);
    % Separate polarity images for motion estimation
    onImg   = zeros(H, W);
    offImg  = zeros(H, W);
    polImg  = zeros(H, W);   % signed polarity accumulator

    for e = 1:nEvt
        r = py(e);  c = px(e);
        if r < 1 || r > H || c < 1 || c > W, continue; end
        evtMask(r, c) = true;
        polImg(r, c) = polImg(r, c) + pol(e);
        if pol(e) > 0
            onImg(r, c) = onImg(r, c) + 1;
        else
            offImg(r, c) = offImg(r, c) + 1;
        end
    end

    % ================================================================
    % STEP 1: Find quad border candidates
    % ================================================================

    % Morphological preprocessing: connect nearby events
    se_dilate = strel('disk', dilateRadius);
    se_close  = strel('disk', closeRadius);
    bwClean = imdilate(evtMask, se_dilate);
    bwClean = imclose(bwClean, se_close);

    % Fill holes (the inside of a marker border may have no events)
    bwFilled = imfill(bwClean, 'holes');

    % Find contours on the filled image
    [B, ~] = bwboundaries(bwFilled, 'noholes');

    % Also find contours on the unfilled (cleaned) image for visualization
    [B_unfilled, ~] = bwboundaries(bwClean, 'noholes');

    % Use regionprops on the filled image for blob stats
    blobStats = regionprops(bwFilled, 'BoundingBox', 'Area', 'Centroid', ...
        'Eccentricity', 'Solidity');

    % Approximate ALL contours to polygons (keep for visualization)
    allPolys = {};       % all polygon approximations (any vertex count)
    polyNVerts = [];     % number of vertices per polygon
    polyAreas = [];      % area of each polygon
    rejectReason = {};   % why each polygon was rejected (empty = kept)

    quads = {};
    for bi = 1:length(B)
        contour = B{bi};   % Nx2 [row, col]
        perim = size(contour, 1);
        if perim < minPerimeter
            continue;
        end

        % Douglas-Peucker polygon approximation
        epsilon = polyEpsFrac * perim;
        poly = approxPoly(contour, epsilon);
        nVerts = size(poly, 1);

        % Convert to [x, y] (col, row)
        corners = [poly(:,2), poly(:,1)];

        % Store for visualization
        allPolys{end+1} = corners; %#ok<AGROW>
        polyNVerts(end+1) = nVerts; %#ok<AGROW>

        a = polyarea(corners(:,1), corners(:,2));
        polyAreas(end+1) = a; %#ok<AGROW>

        % ---- Filtering with rejection reason ----
        if nVerts ~= 4
            rejectReason{end+1} = sprintf('%dv', nVerts); %#ok<AGROW>
            continue;
        end

        if a < minArea
            rejectReason{end+1} = 'small'; %#ok<AGROW>
            continue;
        end
        if a > maxArea
            rejectReason{end+1} = 'large'; %#ok<AGROW>
            continue;
        end

        if ~isConvex(corners)
            rejectReason{end+1} = 'concave'; %#ok<AGROW>
            continue;
        end

        sides = zeros(4,1);
        for s = 1:4
            sides(s) = norm(corners(s,:) - corners(mod(s,4)+1,:));
        end
        if min(sides) < minSideLen
            rejectReason{end+1} = 'short'; %#ok<AGROW>
            continue;
        end

        if max(sides)/min(sides) > maxAspect
            rejectReason{end+1} = 'aspect'; %#ok<AGROW>
            continue;
        end

        rejectReason{end+1} = ''; %#ok<AGROW>  % accepted

        % Order corners: TL, TR, BR, BL
        corners = orderCorners(corners);
        quads{end+1} = corners; %#ok<AGROW>
    end

    % Deduplicate overlapping quads
    if length(quads) > 1
        quads = deduplicateQuads(quads);
    end
    nQuads = length(quads);
    totalQuads = totalQuads + nQuads;

    % ================================================================
    % STEP 2: Estimate motion for each quad from polarity
    % ================================================================
    motionVecs = zeros(nQuads, 2);   % [vx, vy] per quad

    for qi = 1:nQuads
        corners = quads{qi};
        motionVecs(qi,:) = estimateMotion(corners, onImg, offImg, ...
            motionBandWidth, H, W);
    end

    % ---- Log ----
    mvStr = '';
    for qi = 1:nQuads
        mvStr = [mvStr, sprintf('(%.1f,%.1f) ', motionVecs(qi,1), motionVecs(qi,2))]; %#ok<AGROW>
    end
    fprintf(logFid, '%6d  %6d  %5d  %s\n', p, nEvt, nQuads, mvStr);

    if mod(p, 100) == 0 || nQuads > 0
        fprintf('[pkt %4d/%d]  events:%d  quads:%d', ...
            p, numPackets, nEvt, nQuads);
        for qi = 1:nQuads
            fprintf('  motion=(%.1f,%.1f)', motionVecs(qi,1), motionVecs(qi,2));
        end
        fprintf('\n');
    end

    % ---- Visualize (always show) ----
    figure(hFig); clf;

    % (a) Raw event mask
    subplot(2,4,1);
    imshow(evtMask);
    title(sprintf('Raw events | pkt %d | %d evts', p, nEvt));
    axis image;

    % (b) After morphology + blob bounding boxes
    subplot(2,4,2);
    imshow(bwFilled); hold on;
    for bi2 = 1:length(blobStats)
        bb = blobStats(bi2).BoundingBox;
        rectangle('Position', bb, 'EdgeColor', 'c', 'LineWidth', 1);
        ct = blobStats(bi2).Centroid;
        plot(ct(1), ct(2), 'c+', 'MarkerSize', 6);
    end
    title(sprintf('Blobs: %d | dilate+close+fill', length(blobStats)));
    axis image; hold off;

    % (c) All contours (raw boundaries)
    subplot(2,4,3);
    imshow(evtMask); hold on;
    contColors = lines(max(length(B_unfilled), 1));
    for bi2 = 1:length(B_unfilled)
        ct = B_unfilled{bi2};
        ci = mod(bi2-1, size(contColors,1)) + 1;
        plot(ct(:,2), ct(:,1), '-', 'Color', contColors(ci,:), 'LineWidth', 1);
    end
    title(sprintf('Contours: %d (on cleaned mask)', length(B_unfilled)));
    axis image; hold off;

    % (d) Polygon approximations (all, color-coded by vertex count)
    subplot(2,4,4);
    imshow(evtMask); hold on;
    for pi2 = 1:length(allPolys)
        pp = allPolys{pi2};
        nv = polyNVerts(pi2);
        rej = rejectReason{pi2};

        % Color by status: green=accepted quad, yellow=rejected quad, red=not-4-vertex
        if nv == 4 && isempty(rej)
            col = [0 1 0];  lw = 2.5;   % accepted quad: green thick
        elseif nv == 4
            col = [1 1 0];  lw = 1.5;   % rejected quad: yellow
        else
            col = [1 0.3 0.3]; lw = 1;  % not quad: red thin
        end

        plot([pp(:,1); pp(1,1)], [pp(:,2); pp(1,2)], '-', ...
            'Color', col, 'LineWidth', lw);
        center_p = mean(pp, 1);
        if nv == 4 && ~isempty(rej)
            text(center_p(1), center_p(2), rej, ...
                'Color', 'y', 'FontSize', 7, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0 0.5]);
        elseif nv ~= 4
            text(center_p(1), center_p(2), sprintf('%dv', nv), ...
                'Color', [1 0.3 0.3], 'FontSize', 7, ...
                'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0 0.5]);
        end
    end
    title(sprintf('Polygons: %d (green=quad, yellow=rejected, red=non-quad)', ...
        length(allPolys)));
    axis image; hold off;

    % (e) ON/OFF polarity overlay with accepted quads
    subplot(2,4,5);
    combRGB = zeros(H, W, 3);
    if max(onImg(:))  > 0, combRGB(:,:,1) = onImg  / max(onImg(:));  end
    if max(offImg(:)) > 0, combRGB(:,:,3) = offImg / max(offImg(:)); end
    imshow(combRGB * 3); hold on;
    for qi = 1:nQuads
        cc = quads{qi};
        plot([cc(:,1); cc(1,1)], [cc(:,2); cc(1,2)], 'g-', 'LineWidth', 2);
    end
    title(sprintf('ON(red) OFF(blue) | Quads: %d', nQuads));
    axis image; hold off;

    % (f) Accepted quads with motion vectors
    subplot(2,4,6);
    imshow(evtMask); hold on;
    if nQuads > 0
        qColors = lines(nQuads);
        for qi = 1:nQuads
            cc = quads{qi};
            plot([cc(:,1); cc(1,1)], [cc(:,2); cc(1,2)], '-', ...
                'Color', qColors(qi,:), 'LineWidth', 2);
            plot(cc(:,1), cc(:,2), 'o', 'Color', qColors(qi,:), ...
                'MarkerSize', 8, 'MarkerFaceColor', qColors(qi,:));

            center_q = mean(cc, 1);
            mv = motionVecs(qi,:);
            scale = 15;
            quiver(center_q(1), center_q(2), mv(1)*scale, mv(2)*scale, 0, ...
                'Color', 'y', 'LineWidth', 3, 'MaxHeadSize', 2);
            text(center_q(1)+5, center_q(2)-5, ...
                sprintf('v=(%.1f,%.1f)', mv(1), mv(2)), ...
                'Color', 'y', 'FontSize', 9, 'FontWeight', 'bold', ...
                'BackgroundColor', 'k');
        end
    end
    title(sprintf('Quads + motion | %d found', nQuads));
    axis image; hold off;

    % (g) Per-edge polarity (first quad) or blob details
    subplot(2,4,7);
    if nQuads > 0
        showEdgePolarity(quads{1}, onImg, offImg, motionBandWidth, H, W);
        title('Edge polarity (quad 1)');
    else
        % Show blob statistics table
        if ~isempty(blobStats)
            nBlobs = min(length(blobStats), 10);
            txt = sprintf('Blob stats (top %d):\n', nBlobs);
            % Sort by area descending
            [~, sortIdx] = sort([blobStats.Area], 'descend');
            for bi2 = 1:nBlobs
                bs = blobStats(sortIdx(bi2));
                txt = [txt, sprintf('#%d A=%.0f Sol=%.2f Ecc=%.2f\n', ...
                    bi2, bs.Area, bs.Solidity, bs.Eccentricity)]; %#ok<AGROW>
            end
            text(0.05, 0.5, txt, 'Units', 'normalized', ...
                'FontName', 'FixedWidth', 'FontSize', 8, ...
                'VerticalAlignment', 'middle');
        end
        title('Blob stats (no quads)');
        axis off;
    end

    % (h) Motion detail or polygon summary
    subplot(2,4,8);
    if nQuads > 0
        showMotionDetail(quads{1}, polImg, motionVecs(1,:), H, W);
        title('Polarity field + motion');
    else
        % Show polygon vertex distribution
        if ~isempty(polyNVerts)
            vertCounts = histcounts(polyNVerts, 0.5:1:max(polyNVerts)+0.5);
            bar(1:length(vertCounts), vertCounts, 'FaceColor', [0.3 0.6 1]);
            xlabel('Vertices'); ylabel('Count');
            title(sprintf('Polygon vertex distribution (n=%d)', length(polyNVerts)));
        else
            title('No polygons found');
            axis off;
        end
    end

    sgtitle(sprintf('Motion-Aware | Pkt %d/%d | %d events | Blobs: %d | Contours: %d | Polys: %d | Quads: %d', ...
        p, numPackets, nEvt, length(blobStats), length(B_unfilled), length(allPolys), nQuads));
    drawnow;
end

elapsed = toc;
fprintf('\n=== Done ===\n');
fprintf('Total packets : %d\n', numPackets);
fprintf('Total quads   : %d\n', totalQuads);
fprintf('Elapsed time  : %.1f s\n', elapsed);

fprintf(logFid, '-----------------------------------------------------------\n');
fprintf(logFid, 'Total quads  : %d\n', totalQuads);
fprintf(logFid, 'Elapsed time : %.1f s\n', elapsed);
fclose(logFid);
fprintf('Log saved to: motion_output.txt\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%  =========================================================================

%% ---- estimateMotion: polarity-based motion estimation for a quad --------
function mv = estimateMotion(corners, onImg, offImg, bandW, H, W)
%  For each edge of the quad, sample events on both sides (inside/outside).
%  OFF events on one side = black arriving = leading edge → motion towards that side.
%  ON events on one side  = black leaving  = trailing edge → motion away from that side.
%
%  For each edge, the motion component perpendicular to that edge is:
%    positive towards the OFF-dominant side.
%
%  We combine all 4 edges via least-squares to get 2D motion [vx, vy].

    nEdges = 4;
    % Each edge gives one constraint: dot(motion, edgeNormal) = polarity_sign
    A = zeros(nEdges, 2);
    b = zeros(nEdges, 1);

    for ei = 1:4
        p1 = corners(ei, :);
        p2 = corners(mod(ei, 4) + 1, :);

        edgeVec = p2 - p1;
        edgeLen = norm(edgeVec);
        if edgeLen < 1, continue; end

        % Outward-pointing normal (assumes CW winding: TL→TR→BR→BL)
        % For CW, outward normal is (dy, -dx) / len
        normal = [edgeVec(2), -edgeVec(1)] / edgeLen;

        % Sample points along the edge
        nSamples = max(round(edgeLen), 5);
        ts = linspace(0, 1, nSamples);

        onInside  = 0;  offInside  = 0;
        onOutside = 0;  offOutside = 0;

        for si = 1:nSamples
            pt = p1 + ts(si) * edgeVec;

            % Sample on inside (negative normal direction) and outside (positive normal)
            for d = 1:bandW
                % Outside point
                px_out = round(pt(1) + d * normal(1));
                py_out = round(pt(2) + d * normal(2));
                if px_out >= 1 && px_out <= W && py_out >= 1 && py_out <= H
                    onOutside  = onOutside  + onImg(py_out, px_out);
                    offOutside = offOutside + offImg(py_out, px_out);
                end

                % Inside point
                px_in = round(pt(1) - d * normal(1));
                py_in = round(pt(2) - d * normal(2));
                if px_in >= 1 && px_in <= W && py_in >= 1 && py_in <= H
                    onInside  = onInside  + onImg(py_in, px_in);
                    offInside = offInside + offImg(py_in, px_in);
                end
            end
        end

        % Motion signal along this edge's normal:
        % If border is moving in the +normal direction (outward):
        %   Outside: white→black = OFF events dominate
        %   Inside:  black→white = ON events dominate
        % If border is moving in the -normal direction (inward):
        %   Outside: black→white = ON events dominate
        %   Inside:  white→black = OFF events dominate
        %
        % Score: positive = motion in +normal direction
        %   offOutside + onInside → motion outward (+normal)
        %   onOutside + offInside → motion inward (-normal)
        totalEvents = onInside + offInside + onOutside + offOutside;
        if totalEvents < 1
            continue;
        end

        outwardScore = (offOutside + onInside) - (onOutside + offInside);
        motionSign = outwardScore / totalEvents;   % normalized [-1, 1]

        A(ei, :) = normal;
        b(ei)    = motionSign;
    end

    % Least-squares solve: A * [vx; vy] = b
    validRows = any(A ~= 0, 2);
    if sum(validRows) >= 2
        mv = (A(validRows,:) \ b(validRows))';
    else
        mv = [0, 0];
    end
end

%% ---- showEdgePolarity: visualize ON/OFF distribution per edge -----------
function showEdgePolarity(corners, onImg, offImg, bandW, H, W)
    combRGB = zeros(H, W, 3);
    if max(onImg(:))  > 0, combRGB(:,:,1) = onImg  / max(onImg(:));  end
    if max(offImg(:)) > 0, combRGB(:,:,3) = offImg / max(offImg(:)); end
    imshow(combRGB * 3); hold on;

    % Draw quad
    plot([corners(:,1); corners(1,1)], [corners(:,2); corners(1,2)], ...
        'g-', 'LineWidth', 2);

    % For each edge, show ON/OFF counts
    edgeLabels = {'Top', 'Right', 'Bottom', 'Left'};
    for ei = 1:4
        p1 = corners(ei, :);
        p2 = corners(mod(ei,4)+1, :);
        mid = (p1 + p2) / 2;

        edgeVec = p2 - p1;
        edgeLen = norm(edgeVec);
        if edgeLen < 1, continue; end
        normal = [edgeVec(2), -edgeVec(1)] / edgeLen;

        nSamples = max(round(edgeLen), 5);
        ts = linspace(0, 1, nSamples);

        onOut = 0; offOut = 0; onIn = 0; offIn = 0;
        for si = 1:nSamples
            pt = p1 + ts(si) * edgeVec;
            for d = 1:bandW
                px_o = round(pt(1) + d*normal(1));
                py_o = round(pt(2) + d*normal(2));
                if px_o>=1 && px_o<=W && py_o>=1 && py_o<=H
                    onOut  = onOut  + onImg(py_o, px_o);
                    offOut = offOut + offImg(py_o, px_o);
                end
                px_i = round(pt(1) - d*normal(1));
                py_i = round(pt(2) - d*normal(2));
                if px_i>=1 && px_i<=W && py_i>=1 && py_i<=H
                    onIn  = onIn  + onImg(py_i, px_i);
                    offIn = offIn + offImg(py_i, px_i);
                end
            end
        end

        % Label
        text(mid(1), mid(2), ...
            sprintf('%s\nout:+%d/-%d\nin:+%d/-%d', ...
            edgeLabels{ei}, onOut, offOut, onIn, offIn), ...
            'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold', ...
            'BackgroundColor', [0 0 0 0.7], ...
            'HorizontalAlignment', 'center');
    end
    axis image; hold off;
end

%% ---- showMotionDetail: polarity field with motion arrow -----------------
function showMotionDetail(corners, polImg, mv, H, W)
    % Show signed polarity image: red = ON, blue = OFF
    maxP = max(abs(polImg(:)));
    if maxP == 0, maxP = 1; end
    rgbPol = zeros(H, W, 3);
    posP = max(0, polImg) / maxP;     % ON → red
    negP = max(0, -polImg) / maxP;    % OFF → blue
    rgbPol(:,:,1) = posP;
    rgbPol(:,:,3) = negP;
    imshow(rgbPol * 3); hold on;

    % Draw quad
    plot([corners(:,1); corners(1,1)], [corners(:,2); corners(1,2)], ...
        'g-', 'LineWidth', 2);

    % Draw motion vector
    center = mean(corners, 1);
    scale = 20;
    quiver(center(1), center(2), mv(1)*scale, mv(2)*scale, 0, ...
        'Color', 'y', 'LineWidth', 3, 'MaxHeadSize', 2);

    % Draw normal directions for each edge
    for ei = 1:4
        p1 = corners(ei,:);
        p2 = corners(mod(ei,4)+1,:);
        mid = (p1+p2)/2;
        edgeVec = p2 - p1;
        edgeLen = norm(edgeVec);
        if edgeLen < 1, continue; end
        normal = [edgeVec(2), -edgeVec(1)] / edgeLen;
        quiver(mid(1), mid(2), normal(1)*8, normal(2)*8, 0, ...
            'Color', 'c', 'LineWidth', 1.5, 'MaxHeadSize', 2);
    end

    axis image; hold off;
end

%% ---- approxPoly: Douglas-Peucker polygon approximation -----------------
function poly = approxPoly(contour, epsilon)
    N = size(contour, 1);
    if N < 4, poly = contour; return; end

    try
        pts = [contour(:,2), contour(:,1)];  % [col, row] → [x, y]
        totalLen = 0;
        for i = 1:size(pts,1)
            j = mod(i, size(pts,1)) + 1;
            totalLen = totalLen + norm(pts(i,:) - pts(j,:));
        end
        reduced = reducepoly(pts, epsilon / totalLen);
        if size(reduced,1) > 1 && norm(reduced(1,:) - reduced(end,:)) < 1
            reduced = reduced(1:end-1, :);
        end
        poly = [reduced(:,2), reduced(:,1)];
    catch
        poly = douglasPeucker(contour, epsilon);
    end
end

%% ---- douglasPeucker -----------------------------------------------------
function result = douglasPeucker(points, epsilon)
    N = size(points, 1);
    if N <= 2, result = points; return; end
    dmax = 0; idx = 1;
    p1 = points(1,:); p2 = points(end,:);
    for i = 2:N-1
        d = ptLineDist(points(i,:), p1, p2);
        if d > dmax, dmax = d; idx = i; end
    end
    if dmax > epsilon
        r1 = douglasPeucker(points(1:idx,:), epsilon);
        r2 = douglasPeucker(points(idx:end,:), epsilon);
        result = [r1(1:end-1,:); r2];
    else
        result = [points(1,:); points(end,:)];
    end
end

%% ---- ptLineDist ----------------------------------------------------------
function d = ptLineDist(p, a, b)
    ab = b - a; ap = p - a;
    len2 = dot(ab, ab);
    if len2 < 1e-10, d = norm(ap);
    else, d = abs(ab(1)*ap(2) - ab(2)*ap(1)) / sqrt(len2);
    end
end

%% ---- isConvex -----------------------------------------------------------
function ok = isConvex(corners)
    n = size(corners, 1); ok = true; s = 0;
    for i = 1:n
        j = mod(i,n)+1; k = mod(j,n)+1;
        d1 = corners(j,:) - corners(i,:);
        d2 = corners(k,:) - corners(j,:);
        cr = d1(1)*d2(2) - d1(2)*d2(1);
        if s == 0, s = sign(cr);
        elseif s * cr < 0, ok = false; return;
        end
    end
end

%% ---- orderCorners: TL, TR, BR, BL --------------------------------------
function c = orderCorners(pts)
    center = mean(pts, 1);
    angles = atan2(pts(:,2)-center(2), pts(:,1)-center(1));
    [~, ord] = sort(angles);
    pts = pts(ord,:);
    sums = pts(:,1) + pts(:,2);
    [~, tlIdx] = min(sums);
    pts = circshift(pts, -(tlIdx-1), 1);
    d1 = pts(2,:) - pts(1,:);
    d2 = pts(4,:) - pts(1,:);
    cr = d1(1)*d2(2) - d1(2)*d2(1);
    if cr > 0, pts = pts([1,4,3,2],:); end
    c = pts;
end

%% ---- deduplicateQuads ---------------------------------------------------
function merged = deduplicateQuads(allQ)
    n = length(allQ);
    centers = zeros(n,2); sizes = zeros(n,1);
    for k = 1:n
        c = allQ{k};
        centers(k,:) = mean(c,1);
        s = zeros(4,1);
        for j = 1:4, s(j)=norm(c(j,:)-c(mod(j,4)+1,:)); end
        sizes(k) = mean(s);
    end
    keep = true(n,1);
    for k = 1:n
        if ~keep(k), continue; end
        for m = k+1:n
            if ~keep(m), continue; end
            if norm(centers(k,:)-centers(m,:)) < 0.5*min(sizes(k),sizes(m))
                keep(m) = false;
            end
        end
    end
    merged = allQ(keep);
end
