function quads = findQuadCandidates(binaryImg, params)
%FINDQUADCANDIDATES Find square-like regions in a binary image.
%
% Works for BOTH filled blobs AND thin edge patterns by using the
% convex hull area (not blob pixel area) for size filtering.
%
% Input:
%   binaryImg - logical image (white = object)
%   params    - struct with optional fields:
%     .minArea          (default 100)  - minimum CONVEX HULL area
%     .maxArea          (default 50%)  - maximum CONVEX HULL area
%     .maxAspect        (default 2.5)  - max side ratio of min-area rect
%     .minRectangularity(default 0.6)  - convexArea / minRectArea
%     .minPixels        (default 30)   - minimum number of component pixels
%     .connectivity     (default 8)    - 4 or 8 connectivity
%
% Output:
%   quads - cell array of 4x2 [x y] corner matrices (TL, TR, BR, BL)

if nargin < 2, params = struct(); end
if ~isfield(params,'minArea'),          params.minArea          = 100;                  end
if ~isfield(params,'maxArea'),          params.maxArea          = 0.5*numel(binaryImg); end
if ~isfield(params,'maxAspect'),        params.maxAspect        = 2.5;                  end
if ~isfield(params,'minRectangularity'),params.minRectangularity= 0.6;                  end
if ~isfield(params,'minPixels'),        params.minPixels        = 30;                   end
if ~isfield(params,'connectivity'),     params.connectivity     = 8;                    end

cc    = bwconncomp(binaryImg, params.connectivity);
props = regionprops(cc, 'Area', 'ConvexHull', 'ConvexArea');

quads = {};
for i = 1:cc.NumObjects
    % ---- Minimum pixel count (skip tiny noise) ----
    if props(i).Area < params.minPixels, continue; end

    % ---- Area filter on CONVEX HULL (works for edges AND filled blobs) ----
    cxArea = props(i).ConvexArea;
    if cxArea < params.minArea || cxArea > params.maxArea, continue; end

    % ---- Fit minimum-area rectangle to convex hull ----
    hull = props(i).ConvexHull;
    [rect, rectArea, rectSides] = minAreaRect(hull);
    if isempty(rect), continue; end

    % ---- Rectangularity: how square is the convex hull? ----
    rectangularity = cxArea / (rectArea + eps);
    if rectangularity < params.minRectangularity, continue; end

    % ---- Aspect ratio of min-area rectangle ----
    aspect = max(rectSides) / (min(rectSides) + eps);
    if aspect > params.maxAspect, continue; end

    % ---- Order corners TL, TR, BR, BL ----
    rect = orderQuadCorners(rect);
    quads{end+1} = rect; %#ok<AGROW>
end
end

%% =========================================================================
function [rect, rectArea, sides] = minAreaRect(hull)
if size(hull,1) > 1 && norm(hull(1,:) - hull(end,:)) < 1
    hull = hull(1:end-1,:);
end
n = size(hull,1);
if n < 3, rect = []; rectArea = 0; sides = [0 0]; return; end

minA  = inf;  bestR = [];  bestW = 0;  bestH = 0;
for i = 1:n
    j    = mod(i, n) + 1;
    edge = hull(j,:) - hull(i,:);
    ang  = atan2(edge(2), edge(1));
    ca = cos(-ang);  sa = sin(-ang);
    rotH = ([ca -sa; sa ca] * hull')';
    xmn = min(rotH(:,1)); xmx = max(rotH(:,1));
    ymn = min(rotH(:,2)); ymx = max(rotH(:,2));
    w = xmx-xmn;  h = ymx-ymn;  area = w*h;
    if area < minA
        minA = area; bestW = w; bestH = h;
        cRot = [xmn ymn; xmx ymn; xmx ymx; xmn ymx];
        bestR = ([ca sa; -sa ca] * cRot')';
    end
end
rect = bestR;  rectArea = minA;  sides = [bestW, bestH];
end

%% =========================================================================
function ordered = orderQuadCorners(corners)
centroid = mean(corners,1);
angles = atan2(corners(:,2)-centroid(2), corners(:,1)-centroid(1));
[~,si] = sort(angles);  sorted = corners(si,:);
sums = sorted(:,1)+sorted(:,2);  [~,tl] = min(sums);
ordered = circshift(sorted, -(tl-1), 1);
v1 = ordered(2,:)-ordered(1,:);  v2 = ordered(4,:)-ordered(1,:);
if v1(1)*v2(2)-v1(2)*v2(1) > 0, ordered = ordered([1 4 3 2],:); end
end
