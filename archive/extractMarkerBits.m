function [bits, confidence] = extractMarkerBits(frame, corners, gridSize)
%EXTRACTMARKERBITS Extract binary bit pattern from a quadrilateral region
%
% [bits, confidence] = extractMarkerBits(frame, corners, gridSize)
%
% Input:
%   frame    - 2D image (integrated polarity or event-count frame)
%   corners  - 4x2 matrix [x,y] ordered TL, TR, BR, BL
%   gridSize - inner marker grid size (e.g. 4 for a 4x4 ArUco marker)
%
% Output:
%   bits       - gridSize x gridSize binary matrix (0=white, 1=black)
%   confidence - scalar in [0,1]; higher = more reliable detection
%
% The ArUco marker has a totalSize = gridSize+2 cell grid (1-cell black
% border on each side).  The function warps the quadrilateral to a
% canonical square, samples each cell, thresholds, checks the border,
% and returns the inner gridSize x gridSize pattern.

totalSize = gridSize + 2;          % grid cells including border
cellPx    = 50;                    % pixels per cell in canonical view
canonSz   = totalSize * cellPx;    % canonical square side length

% --- Destination corners of the canonical square ---
dst = [1, 1; canonSz, 1; canonSz, canonSz; 1, canonSz];

% --- Compute homography (src corners -> dst) via DLT ---
H = computeHomography(corners, dst);

% --- Warp the source frame into the canonical square (backward map) ---
warped = warpImageBilinear(frame, H, [canonSz, canonSz]);

% --- Sample cell values ---
cellVals = zeros(totalSize);
margin   = round(cellPx * 0.25);   % stay away from cell edges

for r = 1:totalSize
    for c = 1:totalSize
        cy = round((r - 0.5) * cellPx);
        cx = round((c - 0.5) * cellPx);

        r1 = max(1,      cy - margin);
        r2 = min(canonSz, cy + margin);
        c1 = max(1,      cx - margin);
        c2 = min(canonSz, cx + margin);

        region = warped(r1:r2, c1:c2);
        cellVals(r,c) = mean(region(:));
    end
end

% --- Normalize to [0,1] ---
mn = min(cellVals(:));
mx = max(cellVals(:));
if mx - mn < eps
    bits       = zeros(gridSize);
    confidence = 0;
    return;
end
cellVals = (cellVals - mn) / (mx - mn);

% --- Adaptive threshold (Otsu on the cell values) ---
level     = graythresh(cellVals);
binarized = cellVals < level;   % lower value = darker = black = 1

% --- Border check (border cells should all be black = 1) ---
borderMask = false(totalSize);
borderMask(1,:)   = true;  borderMask(end,:) = true;
borderMask(:,1)   = true;  borderMask(:,end) = true;
borderCorrect = mean(binarized(borderMask));

% If border is mostly white, the polarity may be inverted -> flip
if borderCorrect < 0.5
    binarized     = ~binarized;
    borderCorrect = 1 - borderCorrect;
end

% --- Extract inner bits ---
bits = double(binarized(2:end-1, 2:end-1));

% --- Confidence metric ---
%  1) Border correctness  (weight 0.4)
%  2) Bimodality: separation between mean of black and white cells (0.4)
%  3) Warped-image contrast (0.2)
innerVals = cellVals(2:end-1, 2:end-1);
blacks    = innerVals(bits == 1);
whites    = innerVals(bits == 0);
if isempty(blacks) || isempty(whites)
    separation = 0;
else
    separation = abs(mean(whites) - mean(blacks));
end
contrast   = mx - mn;
confidence = 0.4*borderCorrect + 0.4*min(1,separation/0.3) + 0.2*min(1,contrast/0.2);
end

%% ---- Local helpers ---------------------------------------------------------

function H = computeHomography(src, dst)
%COMPUTEHOMOGRAPHY 3x3 projective homography via Direct Linear Transform.
%   src, dst : 4x2  [x y]

n = size(src,1);
A = zeros(2*n, 9);
for i = 1:n
    x  = src(i,1);  y  = src(i,2);
    xp = dst(i,1);  yp = dst(i,2);
    A(2*i-1,:) = [-x -y -1  0  0  0  xp*x  xp*y  xp];
    A(2*i  ,:) = [ 0  0  0 -x -y -1  yp*x  yp*y  yp];
end
[~,~,V] = svd(A, 'econ');
h = V(:,end);
H = reshape(h, [3 3])';
end

function warped = warpImageBilinear(img, H, outSize)
%WARPIMAGEBILINEAR Backward-mapping warp with bilinear interpolation.
%   H maps source -> destination, so we invert for backward mapping.

Hinv = H \ eye(3);             % = inv(H)

[xG, yG] = meshgrid(1:outSize(2), 1:outSize(1));
coords   = [xG(:)'; yG(:)'; ones(1, numel(xG))];

srcC = Hinv * coords;
srcC = bsxfun(@rdivide, srcC, srcC(3,:));

srcX = reshape(srcC(1,:), outSize);
srcY = reshape(srcC(2,:), outSize);

[rows, cols] = size(img);
[X, Y]       = meshgrid(1:cols, 1:rows);
warped       = interp2(X, Y, double(img), srcX, srcY, 'linear', 0);
end
