function quads = detectQuadBlob(activeMask, params)
%DETECTQUADBLOB Detect square marker borders using blob analysis.
%
% Two sub-methods combined:
%   A) imfill + watershed + min-area rect  (closed-edge markers)
%   B) convex-hull on raw edges            (markers with gaps)
%
% Input:
%   activeMask - logical image (white = active event pixels)
%   params     - struct with fields:
%     .minArea   (default 625)    - min convex-hull area
%     .maxArea   (default 22500)  - max convex-hull area
%     .maxAspect (default 2.5)    - max side ratio
%
% Output:
%   quads - cell array of 4x2 [x y] matrices (TL, TR, BR, BL)

if nargin < 2, params = struct(); end
if ~isfield(params,'minArea'),   params.minArea   = 625;   end
if ~isfield(params,'maxArea'),   params.maxArea   = 22500; end
if ~isfield(params,'maxAspect'), params.maxAspect = 2.5;   end

qp.minArea           = params.minArea;
qp.maxArea           = params.maxArea;
qp.maxAspect         = params.maxAspect;
qp.minRectangularity = 0.6;

% ---- Method A: fill-based ----
bw = imfill(activeMask, 'holes');
bw = bwareaopen(bw, 100);

% Watershed to split touching blobs
D      = bwdist(~bw);
peaks  = imextendedmax(D, 2);
peaks  = imdilate(peaks, strel('disk',1));
imp    = imimposemin(-D, peaks | ~bw);
W      = watershed(imp);
bw(W == 0) = false;
bw     = bwareaopen(bw, 100);

quads_A = findQuadCandidates(bw, qp);

% ---- Method B: convex-hull on raw edges (4-connectivity) ----
qp_e = qp;
qp_e.connectivity = 4;
qp_e.minPixels    = 30;
quads_B = findQuadCandidates(activeMask, qp_e);

% ---- Combine & deduplicate ----
quads = deduplicateQuads(quads_A, quads_B);
end

%% =========================================================================
function merged = deduplicateQuads(quadsA, quadsB)
    all = [quadsA, quadsB];
    if isempty(all), merged = {}; return; end

    centers = zeros(length(all), 2);
    sizes   = zeros(length(all), 1);
    for k = 1:length(all)
        c = all{k};
        centers(k,:) = mean(c, 1);
        s = zeros(4,1);
        for j = 1:4, s(j) = norm(c(j,:) - c(mod(j,4)+1,:)); end
        sizes(k) = mean(s);
    end

    keep = true(length(all),1);
    for k = 1:length(all)
        if ~keep(k), continue; end
        for m = k+1:length(all)
            if ~keep(m), continue; end
            if norm(centers(k,:)-centers(m,:)) < 0.5*min(sizes(k),sizes(m))
                keep(m) = false;
            end
        end
    end
    merged = all(keep);
end
