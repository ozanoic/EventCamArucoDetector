function quads = detectQuadLines(activeMask, params)
%DETECTQUADLINES Detect square marker borders using Hough line pairs.
%
% Pipeline:
%   1. Detect line segments via Hough transform
%   2. Filter by length [minLen, maxLen]
%   3. Find parallel line pairs (similar angle, similar length)
%   4. Find perpendicular pair-pairs that form a closed quadrilateral
%   5. Return the 4 intersection corners
%
% Input:
%   activeMask - logical image (white = active event pixels)
%   params     - struct with optional fields:
%     .minLen       (default 20)   - minimum line segment length (px)
%     .maxLen       (default 150)  - maximum line segment length (px)
%     .angleTol     (default 15)   - parallel angle tolerance (degrees)
%     .lengthRatio  (default 0.4)  - max allowed |L1-L2|/max(L1,L2)
%     .gapRatio     (default [0.3, 3.0]) - gap/length ratio range
%     .numPeaks     (default 30)   - max Hough peaks to consider
%
% Output:
%   quads - cell array of 4x2 [x y] matrices (TL, TR, BR, BL)

if nargin < 2, params = struct(); end
if ~isfield(params,'minLen'),      params.minLen      = 20;        end
if ~isfield(params,'maxLen'),      params.maxLen       = 150;       end
if ~isfield(params,'angleTol'),    params.angleTol     = 15;        end
if ~isfield(params,'lengthRatio'), params.lengthRatio  = 0.4;       end
if ~isfield(params,'gapRatio'),    params.gapRatio     = [0.3 3.0]; end
if ~isfield(params,'numPeaks'),    params.numPeaks     = 30;        end

quads = {};

% ---- 1. Edge detection on active mask ----
%   The active mask is already edge-like, but Hough works best with
%   thin single-pixel edges.  Use bwmorph to thin.
edges = bwmorph(activeMask, 'thin', Inf);

% ---- 2. Hough transform ----
[H, theta, rho] = hough(edges, 'RhoResolution', 1, 'Theta', -90:0.5:89.5);
peaks = houghpeaks(H, params.numPeaks, 'Threshold', ceil(0.2 * max(H(:))));
if isempty(peaks), return; end

segments = houghlines(edges, theta, rho, peaks, ...
    'FillGap', 8, 'MinLength', params.minLen);
if isempty(segments), return; end

% ---- 3. Compute length & midpoint for each segment ----
nSeg = length(segments);
seg_data = struct('p1',{},'p2',{},'len',{},'angle',{},'mid',{});
for k = 1:nSeg
    p1  = segments(k).point1;       % [col, row] = [x, y]
    p2  = segments(k).point2;
    d   = p2 - p1;
    len = norm(d);
    if len < params.minLen || len > params.maxLen, continue; end
    seg_data(end+1).p1    = p1;
    seg_data(end).p2      = p2;
    seg_data(end).len     = len;
    seg_data(end).angle   = atan2d(d(2), d(1));   % -180..180
    seg_data(end).mid     = (p1 + p2) / 2;
end
nSeg = length(seg_data);
if nSeg < 4, return; end

% ---- 4. Find parallel line pairs ----
%   Two lines are a "parallel pair" if:
%     - angles differ by < angleTol  (mod 180)
%     - lengths are similar
%     - perpendicular gap between them is in [gapRatio * avgLen]
pairs = [];   % Nx2 indices into seg_data
for a = 1:nSeg
    for b = a+1:nSeg
        % Angle check (parallel)
        dAng = abs(angleDiff(seg_data(a).angle, seg_data(b).angle));
        if dAng > params.angleTol, continue; end

        % Length similarity
        La = seg_data(a).len;
        Lb = seg_data(b).len;
        if abs(La-Lb)/max(La,Lb) > params.lengthRatio, continue; end

        % Perpendicular gap between the two lines
        gap = perpGap(seg_data(a), seg_data(b));
        avgLen = (La + Lb) / 2;
        ratio  = gap / avgLen;
        if ratio < params.gapRatio(1) || ratio > params.gapRatio(2), continue; end

        % Overlap check: the midpoints projected onto the line direction
        % should overlap (the lines should be "beside" each other)
        overlap = projectionOverlap(seg_data(a), seg_data(b));
        if overlap < 0.3, continue; end

        pairs = [pairs; a b]; %#ok<AGROW>
    end
end
if isempty(pairs), return; end

% ---- 5. Find two perpendicular pairs that form a square ----
nPairs = size(pairs, 1);
for p1 = 1:nPairs
    for p2 = p1+1:nPairs
        idxA = pairs(p1,:);
        idxB = pairs(p2,:);

        % All 4 segments must be distinct
        if any(ismember(idxA, idxB)), continue; end

        % The two pairs must be roughly perpendicular
        angA = meanAngle(seg_data(idxA(1)).angle, seg_data(idxA(2)).angle);
        angB = meanAngle(seg_data(idxB(1)).angle, seg_data(idxB(2)).angle);
        dPerp = abs(angleDiff(angA, angB));
        if abs(dPerp - 90) > params.angleTol, continue; end

        % Compute the 4 intersection corners of the 4 lines
        corners = quadFromPairs(seg_data, idxA, idxB);
        if isempty(corners), continue; end

        % Validate: area in range, aspect ratio OK
        area = polyarea(corners(:,1), corners(:,2));
        minA = params.minLen^2 * 0.3;
        maxA = params.maxLen^2 * 3;
        if area < minA || area > maxA, continue; end

        sides = zeros(4,1);
        for s = 1:4
            sides(s) = norm(corners(s,:) - corners(mod(s,4)+1,:));
        end
        if max(sides)/min(sides) > 3, continue; end

        corners = orderQuadCorners(corners);
        quads{end+1} = corners; %#ok<AGROW>
    end
end

% Deduplicate quads with close centers
quads = deduplicateQuads(quads);
end

%% =========================================================================
%                       LOCAL  HELPERS
%% =========================================================================

function d = angleDiff(a, b)
%ANGLEDIFF Smallest angle between two directions (mod 180).
    d = mod(a - b + 180, 360) - 180;
    d = min(abs(d), 180 - abs(d));
end

function a = meanAngle(a1, a2)
%MEANANGLE Average of two angles (mod 180).
    a = a1 + angleDiff(a1, a2) / 2;
end

function g = perpGap(s1, s2)
%PERPGAP Perpendicular distance between two parallel-ish segments.
%   Uses distance from s2's midpoint to the infinite line through s1.
    d  = s1.p2 - s1.p1;
    n  = [-d(2), d(1)];             % normal
    n  = n / (norm(n) + eps);
    g  = abs(dot(s2.mid - s1.mid, n));
end

function ov = projectionOverlap(s1, s2)
%PROJECTIONOVERLAP Fractional overlap when both segments are projected
%   onto their mean direction.
    dir = s1.p2 - s1.p1;
    dir = dir / (norm(dir) + eps);

    % Project all 4 endpoints onto the direction
    t = [dot(s1.p1, dir), dot(s1.p2, dir), ...
         dot(s2.p1, dir), dot(s2.p2, dir)];
    r1 = sort(t(1:2));
    r2 = sort(t(3:4));

    lo = max(r1(1), r2(1));
    hi = min(r1(2), r2(2));
    overlap = max(0, hi - lo);
    span    = max(r1(2) - r1(1), r2(2) - r2(1));
    ov = overlap / (span + eps);
end

function corners = quadFromPairs(seg, idxA, idxB)
%QUADFROMPAIRS Compute 4 corners by intersecting two line pairs.
%   Pair A has lines idxA(1), idxA(2);  Pair B has idxB(1), idxB(2).
%   We intersect each line of pair A with each line of pair B → 4 points.

    lines = [idxA(1) idxB(1);   % corner 1: A1 ∩ B1
             idxA(1) idxB(2);   % corner 2: A1 ∩ B2
             idxA(2) idxB(2);   % corner 3: A2 ∩ B2
             idxA(2) idxB(1)];  % corner 4: A2 ∩ B1

    corners = zeros(4, 2);
    for k = 1:4
        pt = lineIntersect(seg(lines(k,1)), seg(lines(k,2)));
        if isempty(pt), corners = []; return; end
        corners(k,:) = pt;
    end
end

function pt = lineIntersect(s1, s2)
%LINEINTERSECT Intersection of two infinite lines defined by segments.
    d1 = s1.p2 - s1.p1;
    d2 = s2.p2 - s2.p1;
    cross_d = d1(1)*d2(2) - d1(2)*d2(1);
    if abs(cross_d) < 1e-6, pt = []; return; end   % parallel
    dp = s2.p1 - s1.p1;
    t  = (dp(1)*d2(2) - dp(2)*d2(1)) / cross_d;
    pt = s1.p1 + t * d1;    % [x, y]
end

function ordered = orderQuadCorners(corners)
    centroid = mean(corners,1);
    angles = atan2(corners(:,2)-centroid(2), corners(:,1)-centroid(1));
    [~,si] = sort(angles);  sorted = corners(si,:);
    sums = sorted(:,1)+sorted(:,2);  [~,tl] = min(sums);
    ordered = circshift(sorted, -(tl-1), 1);
    v1 = ordered(2,:)-ordered(1,:); v2 = ordered(4,:)-ordered(1,:);
    if v1(1)*v2(2)-v1(2)*v2(1) > 0, ordered = ordered([1 4 3 2],:); end
end

function merged = deduplicateQuads(quadsIn)
    if isempty(quadsIn), merged = {}; return; end
    centers = zeros(length(quadsIn),2);
    sizes   = zeros(length(quadsIn),1);
    for k = 1:length(quadsIn)
        c = quadsIn{k};
        centers(k,:) = mean(c,1);
        s = zeros(4,1);
        for j=1:4, s(j)=norm(c(j,:)-c(mod(j,4)+1,:)); end
        sizes(k) = mean(s);
    end
    keep = true(length(quadsIn),1);
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
