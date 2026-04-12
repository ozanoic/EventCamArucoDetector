function dict = getArucoDict(dictType)
%GETARUCODICT Get or generate an ArUco marker dictionary
%
% dict = getArucoDict(dictType)
%
% Input:
%   dictType - 'DICT_4X4_50' (default), 'DICT_5X5_50', 'DICT_6X6_50'
%
% Output:
%   dict.gridSize   - inner grid dimension (e.g. 4)
%   dict.numMarkers - number of markers in the dictionary
%   dict.markers    - cell array of gridSize x gridSize binary matrices
%                     Convention: 0 = white cell, 1 = black cell
%
% The dictionary is generated once using a greedy Hamming-distance
% algorithm with a fixed random seed for reproducibility, then cached
% as a .mat file next to this function.
%
% NOTE: These markers are NOT identical to OpenCV's DICT_4X4_50.
% To use the exact OpenCV dictionary, generate the byte data in Python:
%
%   import cv2, numpy as np
%   d = cv2.aruco.getPredefinedDictionary(cv2.aruco.DICT_4X4_50)
%   for i in range(50):
%       row_bytes = []
%       for r in range(4):
%           row_bytes.append(d.bytesList[i, r, 0])
%       print(row_bytes)
%
% Then replace the cached .mat file or modify this function accordingly.

if nargin < 1, dictType = 'DICT_4X4_50'; end

switch upper(dictType)
    case 'DICT_4X4_50',  gridSize = 4; numMarkers = 50; minDist = 4;
    case 'DICT_5X5_50',  gridSize = 5; numMarkers = 50; minDist = 5;
    case 'DICT_6X6_50',  gridSize = 6; numMarkers = 50; minDist = 6;
    otherwise, error('getArucoDict:unknown', 'Unknown dictionary: %s', dictType);
end

nBits = gridSize^2;

% --- Try loading from cache ---
cacheFile = fullfile(fileparts(mfilename('fullpath')), ...
    sprintf('aruco_dict_%dx%d_%d.mat', gridSize, gridSize, numMarkers));

if exist(cacheFile, 'file')
    data = load(cacheFile);
    dict = data.dict;
    return;
end

% --- Generate dictionary ---
fprintf('Generating ArUco dictionary (%dx%d, %d markers)...\n', ...
    gridSize, gridSize, numMarkers);

rng_prev = rng;
rng(42, 'twister');

nTotal    = 2^nBits;
order     = randperm(nTotal) - 1;   % deterministic shuffle
selected  = zeros(0, nBits);

for idx = 1:nTotal
    val  = order(idx);
    bits = bitget(val, nBits:-1:1);   % MSB-first

    % Skip too-uniform patterns
    s = sum(bits);
    if s < 3 || s > nBits - 3, continue; end

    mat = reshape(bits, [gridSize, gridSize]);

    % Compute all 4 rotations
    rots = zeros(4, nBits);
    tmp  = mat;
    for r = 1:4
        rots(r,:) = tmp(:)';
        tmp = rot90(tmp);
    end

    % Self-symmetry check: each rotation must differ enough
    selfOk = true;
    for r = 2:4
        if sum(rots(1,:) ~= rots(r,:)) < minDist
            selfOk = false;
            break;
        end
    end
    if ~selfOk, continue; end

    % Hamming distance against all already-selected markers
    if ~isempty(selected)
        dists = zeros(4, size(selected,1));
        for r = 1:4
            dists(r,:) = sum(bsxfun(@ne, rots(r,:), selected), 2)';
        end
        if min(dists(:)) < minDist, continue; end
    end

    selected = [selected; bits]; %#ok<AGROW>
    if size(selected,1) >= numMarkers, break; end
end

rng(rng_prev);   % restore RNG state

% --- Build output struct ---
dict.gridSize   = gridSize;
dict.numMarkers = size(selected, 1);
dict.markers    = cell(dict.numMarkers, 1);
for i = 1:dict.numMarkers
    dict.markers{i} = reshape(selected(i,:), [gridSize, gridSize]);
end

save(cacheFile, 'dict');
fprintf('Dictionary generated: %d markers  ->  %s\n', dict.numMarkers, cacheFile);
end
