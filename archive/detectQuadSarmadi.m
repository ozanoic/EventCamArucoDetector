function [quads, dbg] = detectQuadSarmadi(sae_on, sae_off, currentTime, params)
%DETECTQUADSARMADI Detect square marker borders using Sarmadi et al. (2021).
%
%  Algorithm (adapted from "Detection of Binary Square Fiducial Markers
%  Using an Event Camera"):
%    1. Build separate ON / OFF normalized-timestamp images from SAE
%    2. Preprocess: majority-vote neighbourhood filter + masked Gaussian blur
%    3. Detect line segments in each image (Hough-based)
%    4. Match ON-OFF segment pairs (angle < 30°, length ratio < 2x,
%       endpoint-projection overlap) → 4-corner quad candidates
%    5. Filter candidates by area and aspect ratio
%
%  Inputs
%    sae_on      [H x W]  SAE for ON  (polarity = +1) events
%    sae_off     [H x W]  SAE for OFF (polarity = -1) events
%    currentTime scalar   current timestamp (µs)
%    params      struct   (all fields optional, sensible defaults below)
%
%  Outputs
%    quads   cell array of 4×2 [x y] corner matrices
%    dbg     struct with intermediate images for visualization:
%            .on_img, .off_img, .on_mask, .off_mask,
%            .on_segs, .off_segs   (struct arrays with .p1 .p2 .len .theta)

%% ---- defaults -------------------------------------------------------------
if nargin < 4, params = struct(); end
def = struct('activeWindow',100000, 'minSegLen',25, 'maxSegLen',250, ...
             'angleTol',30, 'lengthRatio',0.5, 'numPeaks',50, ...
             'minArea',400, 'maxArea',22500, 'maxAspect',3.0, ...
             'houghFillGap',5, 'houghThreshFrac',0.15);
flds = fieldnames(def);
for k = 1:numel(flds)
    if ~isfield(params, flds{k}), params.(flds{k}) = def.(flds{k}); end
end

[H, W] = size(sae_on);
aw     = params.activeWindow;

%% ---- 1. Build ON / OFF event images --------------------------------------
% Active masks (recent events in the time window)
on_mask  = (sae_on  > 0) & ((currentTime - sae_on)  <= aw);
off_mask = (sae_off > 0) & ((currentTime - sae_off) <= aw);

% Normalized timestamp images  → "1 - I_norm" (paper Eq. 6)
on_img  = buildNormImg(sae_on,  on_mask);
off_img = buildNormImg(sae_off, off_mask);

%% ---- 2. Preprocess --------------------------------------------------------
[on_img,  on_mask]  = preprocessFrame(on_img,  on_mask);
[off_img, off_mask] = preprocessFrame(off_img, off_mask);

on_img  = maskedGaussBlur(on_img,  on_mask,  3);
off_img = maskedGaussBlur(off_img, off_mask, 3);

%% ---- 3. Detect line segments in each image --------------------------------
on_segs  = detectSegmentsHough(on_img,  on_mask,  params);
off_segs = detectSegmentsHough(off_img, off_mask, params);

%% ---- 4. Match ON↔OFF pairs → quad candidates ----------------------------
quads = matchSegmentPairs(on_segs, off_segs, params, [H, W]);

%% ---- debug output --------------------------------------------------------
if nargout > 1
    dbg.on_img   = on_img;
    dbg.off_img  = off_img;
    dbg.on_mask  = on_mask;
    dbg.off_mask = off_mask;
    dbg.on_segs  = on_segs;
    dbg.off_segs = off_segs;
end
end

%% =========================================================================
%  Helper: build 1-I_norm image from SAE
%  =========================================================================
function img = buildNormImg(sae, mask)
    img = zeros(size(sae));
    ts  = sae(mask);
    if isempty(ts), return; end
    tMin = min(ts);
    tMax = max(ts);
    if tMax > tMin
        img(mask) = 1 - (sae(mask) - tMin) / (tMax - tMin);
    else
        img(mask) = 1;
    end
end

%% =========================================================================
%  Helper: majority-vote neighbourhood filter  (paper Eq. 7-11)
%  =========================================================================
function [imgOut, maskOut] = preprocessFrame(img, mask)
    [H, W]  = size(img);
    imgOut  = img;
    maskOut = mask;
    for y = 1:H
        yn_lo = max(1, y-1);  yn_hi = min(H, y+1);
        for x = 1:W
            xn_lo = max(1, x-1);  xn_hi = min(W, x+1);
            nbPatch  = mask(yn_lo:yn_hi, xn_lo:xn_hi);
            nTotal   = numel(nbPatch) - 1;          % exclude self
            cy = y - yn_lo + 1;  cx = x - xn_lo + 1;
            selfMask = nbPatch(cy, cx);
            nbPatch(cy, cx) = false;
            nValid = sum(nbPatch(:));
            halfN  = nTotal / 2;

            if nValid > halfN && ~selfMask
                % Hole-fill: most neighbours active → turn on
                valPatch = img(yn_lo:yn_hi, xn_lo:xn_hi);
                valPatch(cy, cx) = 0;
                nbActive = valPatch(nbPatch == 1);
                imgOut(y, x)  = mean(nbActive);
                maskOut(y, x) = true;
            elseif nValid < halfN && selfMask
                % Noise remove: most neighbours inactive → turn off
                imgOut(y, x)  = 0;
                maskOut(y, x) = false;
            end
        end
    end
end

%% =========================================================================
%  Helper: masked Gaussian blur  (paper Eq. 12)
%  =========================================================================
function out = maskedGaussBlur(img, mask, ksize)
    sigma   = 0.3*((ksize-1)*0.5-1) + 0.8;
    maskDbl = double(mask);
    imgBlur  = imgaussfilt(img,    sigma, 'FilterSize', ksize, 'Padding', 0);
    maskBlur = imgaussfilt(maskDbl, sigma, 'FilterSize', ksize, 'Padding', 0);
    out = zeros(size(img));
    valid = maskBlur > 0;
    out(valid) = imgBlur(valid) ./ maskBlur(valid);
end

%% =========================================================================
%  Helper: detect line segments using Hough transform
%  =========================================================================
function segs = detectSegmentsHough(img, mask, params)
    segs = struct('p1',{}, 'p2',{}, 'len',{}, 'theta',{});

    % Convert to uint8 for edge/Hough processing
    img8 = uint8(img * 255);

    % Thin the mask to get 1-pixel-wide edges
    bw   = mask;
    thin = bwmorph(bw, 'thin', Inf);

    % Hough transform
    [Ht, theta, rho] = hough(thin, 'RhoResolution', 1, 'Theta', -90:0.5:89.5);
    maxH = max(Ht(:));
    if maxH == 0, return; end

    peaks = houghpeaks(Ht, params.numPeaks, ...
        'Threshold', ceil(params.houghThreshFrac * maxH));
    if isempty(peaks), return; end

    lines = houghlines(thin, theta, rho, peaks, ...
        'FillGap', params.houghFillGap, 'MinLength', params.minSegLen);

    % Convert to our segment struct, filter by length
    for k = 1:length(lines)
        p1  = lines(k).point1;   % [col, row]
        p2  = lines(k).point2;
        dxy = p2 - p1;
        len = norm(dxy);
        if len < params.minSegLen || len > params.maxSegLen, continue; end
        th  = atan2(dxy(2), dxy(1));       % radians
        segs(end+1) = struct('p1',p1, 'p2',p2, 'len',len, 'theta',th); %#ok<AGROW>
    end
end

%% =========================================================================
%  Helper: match ON↔OFF segment pairs  (paper Section III-D, Eq. 21)
%  =========================================================================
function quads = matchSegmentPairs(on_segs, off_segs, params, imSize)
    quads = {};
    if isempty(on_segs) || isempty(off_segs), return; end

    angleTolRad = params.angleTol * pi / 180;

    for i = 1:length(on_segs)
        on = on_segs(i);
        if on.len < params.minSegLen, continue; end

        for j = 1:length(off_segs)
            off = off_segs(j);
            if off.len < params.minSegLen, continue; end

            % ---- Length compatibility (within 2x) ----
            if on.len > 2 * off.len, continue; end
            if off.len > 2 * on.len, continue; end

            % ---- Angle compatibility (< angleTol degrees) ----
            dth = abs(on.theta - off.theta);
            if dth > pi, dth = 2*pi - dth; end
            if dth > pi/2, dth = pi - dth; end          % lines are undirected
            if dth > angleTolRad, continue; end

            % ---- Projection overlap ----
            %  At least one endpoint of one segment must project within
            %  the span of the other segment.
            proj = projectsIn(on.p1, off) || projectsIn(on.p2, off) || ...
                   projectsIn(off.p1, on) || projectsIn(off.p2, on);
            if ~proj, continue; end

            % ---- Build 4-corner quad from the two segments ----
            corners = orderCorners([on.p1; on.p2; off.p1; off.p2]);

            % ---- Area check ----
            a = polyarea(corners(:,1), corners(:,2));
            if a < params.minArea || a > params.maxArea, continue; end

            % ---- Aspect ratio check ----
            sides = zeros(4,1);
            for s = 1:4
                sides(s) = norm(corners(s,:) - corners(mod(s,4)+1,:));
            end
            ar = max(sides) / max(min(sides), 1);
            if ar > params.maxAspect, continue; end

            quads{end+1} = corners; %#ok<AGROW>
        end
    end

    % Deduplicate overlapping detections
    if length(quads) > 1
        quads = deduplicateQuads(quads);
    end
end

%% =========================================================================
%  Helper: does point p project within segment ls?
%  =========================================================================
function ok = projectsIn(p, ls)
    % Project p onto the infinite line through ls.p1→ls.p2
    v  = ls.p2 - ls.p1;         % segment direction
    w  = p     - ls.p1;
    t  = dot(w, v) / dot(v, v); % parametric position along segment
    ok = (t >= 0) && (t <= 1);
    if ~ok, return; end
    % Also check perpendicular distance is reasonable (< segment length)
    pp   = ls.p1 + t * v;       % projected point
    dist = norm(pp - p);
    ok   = dist <= ls.len;      % perpendicular gap ≤ segment length
end

%% =========================================================================
%  Helper: order 4 points into a proper quadrilateral  (TL TR BR BL)
%  =========================================================================
function c = orderCorners(pts)
    % Sort by x, then assign TL/BL from left pair, TR/BR from right pair
    [~, ix] = sort(pts(:,1));
    left  = pts(ix(1:2), :);
    right = pts(ix(3:4), :);
    % Within each pair: lower y = top
    if left(1,2)  > left(2,2),  left  = left([2 1],:);  end
    if right(1,2) > right(2,2), right = right([2 1],:); end
    c = [left(1,:); right(1,:); right(2,:); left(2,:)];  % TL TR BR BL
end

%% =========================================================================
%  Helper: deduplicate overlapping quads
%  =========================================================================
function merged = deduplicateQuads(allQ)
    n = length(allQ);
    centers = zeros(n, 2);
    sizes   = zeros(n, 1);
    for k = 1:n
        c = allQ{k};
        centers(k,:) = mean(c, 1);
        s = zeros(4,1);
        for j = 1:4, s(j) = norm(c(j,:) - c(mod(j,4)+1,:)); end
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
