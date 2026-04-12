function results = detectArucoFromEvents(events, varargin)
%DETECTARUCOFROMEVENTS Detect ArUco markers from event camera data
%
% results = detectArucoFromEvents(events)
% results = detectArucoFromEvents(events, 'Param', value, ...)
%
% Pipeline:
%   1. Build event-count frame (edge map) and integrated-polarity frame
%   2. Preprocess: blur, threshold, morphological closing, hole filling
%   3. Detect quadrilateral candidates (Douglas-Peucker approximation)
%   4. For each quad: perspective-warp, sample cells, threshold bits
%   5. Match extracted bits against ArUco dictionary (4 rotations)
%   6. Return detected marker IDs, corners, bits, and confidence
%
% Input:
%   events - Nx4 double  [x, y, polarity, timestamp]
%            x,y 0-indexed; polarity {0,1}
%
% Name-Value Parameters:
%   SensorSize  - [height, width]  (default: auto-detect from data)
%   TimeWindow  - scalar; use only the last TimeWindow time-units of
%                 events  (default: Inf = use all)
%   DictType    - string  (default: 'DICT_4X4_50')
%   MinArea     - minimum quad area in px  (default: 100)
%   MaxArea     - maximum quad area in px  (default: 50% of image)
%   MorphRadius - radius of morphological structuring element (default: 3)
%   Visualize   - logical  (default: false)
%
% Output:
%   results  struct with fields:
%     .ids          - 1xM  detected marker IDs (0-based; -1 = unmatched)
%     .corners      - 1xM  cell of 4x2 [x,y] corner matrices
%     .bits         - 1xM  cell of gridSize x gridSize binary matrices
%     .confidences  - 1xM  confidence scores
%     .numDetected  - scalar
%     .edgeFrame    - event-count image  (useful for debugging)
%     .intFrame     - integrated-polarity image
%     .binaryImg    - preprocessed binary image used for quad search
%     .allQuads     - all quad candidates before bit matching

%% ---- Parse inputs ----------------------------------------------------------
p = inputParser;
addParameter(p, 'SensorSize',  [],    @isnumeric);
addParameter(p, 'TimeWindow',  inf,   @isnumeric);
addParameter(p, 'DictType',    'DICT_4X4_50', @ischar);
addParameter(p, 'MinArea',     100,   @isnumeric);
addParameter(p, 'MaxArea',     [],    @isnumeric);
addParameter(p, 'MorphRadius', 3,     @isnumeric);
addParameter(p, 'Visualize',   false, @islogical);
parse(p, varargin{:});
opts = p.Results;

%% ---- Sensor size -----------------------------------------------------------
if isempty(opts.SensorSize)
    opts.SensorSize = [ceil(max(events(:,2)))+1, ceil(max(events(:,1)))+1];
end
if isempty(opts.MaxArea)
    opts.MaxArea = 0.5 * prod(opts.SensorSize);
end

%% ---- Time windowing --------------------------------------------------------
if isfinite(opts.TimeWindow)
    tEnd   = max(events(:,4));
    tStart = tEnd - opts.TimeWindow;
    events = events(events(:,4) >= tStart, :);
end

%% ---- Load dictionary -------------------------------------------------------
dict = getArucoDict(opts.DictType);

%% ---- Build frames ----------------------------------------------------------
edgeFrame = buildEventFrame(events, opts.SensorSize, 'count');
intFrame  = buildEventFrame(events, opts.SensorSize, 'integrated');

%% ---- Preprocess edge frame for quad detection ------------------------------
bw = preprocessEdgeFrame(edgeFrame, opts.MorphRadius);

%% ---- Find quadrilateral candidates ----------------------------------------
qp.minArea = opts.MinArea;
qp.maxArea = opts.MaxArea;
quads = findQuadCandidates(bw, qp);

%% ---- Extract & match each candidate ---------------------------------------
ids   = [];
corns = {};
bitsC = {};
confs = [];

for i = 1:length(quads)
    corners = quads{i};

    % Try integrated frame (better for intensity) and edge frame
    [bits1, c1] = extractMarkerBits(intFrame,  corners, dict.gridSize);
    [bits2, c2] = extractMarkerBits(edgeFrame, corners, dict.gridSize);

    if c1 >= c2
        bits = bits1;  conf = c1;
    else
        bits = bits2;  conf = c2;
    end

    if conf < 0.3, continue; end

    [id, matchedBits] = matchDictionary(bits, dict);

    ids(end+1)    = id;         %#ok<AGROW>
    corns{end+1}  = corners;    %#ok<AGROW>
    bitsC{end+1}  = matchedBits;%#ok<AGROW>
    confs(end+1)  = conf;       %#ok<AGROW>
end

%% ---- Assemble output -------------------------------------------------------
results.ids          = ids;
results.corners      = corns;
results.bits         = bitsC;
results.confidences  = confs;
results.numDetected  = numel(ids);
results.edgeFrame    = edgeFrame;
results.intFrame     = intFrame;
results.binaryImg    = bw;
results.allQuads     = quads;

%% ---- Visualize -------------------------------------------------------------
if opts.Visualize
    visualizeDetection(results, quads);
end
end

%% ============================================================================
%                        LOCAL  FUNCTIONS
%% ============================================================================

function bw = preprocessEdgeFrame(edgeFrame, morphRadius)
%PREPROCESSEDGEFRAME Convert event-count frame to a binary image suitable
%   for quadrilateral contour detection.

    if max(edgeFrame(:)) == 0
        bw = false(size(edgeFrame));
        return;
    end

    % Normalize
    frame = double(edgeFrame) / max(double(edgeFrame(:)));

    % Gaussian blur to smooth noise
    frame = imgaussfilt(frame, 1.5);

    % Otsu threshold
    level = graythresh(frame);
    bw    = imbinarize(frame, max(level, 0.05));

    % Morphological closing to bridge small gaps in edges
    se = strel('disk', morphRadius);
    bw = imclose(bw, se);

    % Fill enclosed regions so that marker interior becomes one blob
    bw = imfill(bw, 'holes');

    % Remove tiny components
    bw = bwareaopen(bw, 50);
end

function [id, matchedBits] = matchDictionary(bits, dict)
%MATCHDICTIONARY Try all 4 rotations against every dictionary entry.
%   Allows up to 1-bit error correction.

    id          = -1;
    matchedBits = bits;
    minDist     = inf;

    for rot = 0:3
        rb = rot90(bits, rot);
        for m = 1:dict.numMarkers
            d = sum(rb(:) ~= dict.markers{m}(:));
            if d < minDist
                minDist = d;
                if d <= 1
                    id          = m - 1;   % 0-indexed
                    matchedBits = rb;
                end
            end
        end
    end
end

function visualizeDetection(results, quads)
%VISUALIZEDETECTION Six-panel figure showing intermediate and final results.

    figure('Name','ArUco Detection from Events','Position',[80 80 1400 800]);

    % 1 - Event count frame
    subplot(2,3,1);
    imagesc(results.edgeFrame); colormap(gca,'hot'); colorbar;
    title('Event Count Frame'); axis image;

    % 2 - Integrated polarity frame
    subplot(2,3,2);
    imagesc(results.intFrame); colormap(gca,'gray'); colorbar;
    title('Integrated Polarity Frame'); axis image;

    % 3 - Binary image after preprocessing
    subplot(2,3,3);
    imshow(results.binaryImg);
    title('Preprocessed Binary'); axis image;

    % 4 - All quad candidates (green)
    subplot(2,3,4);
    imagesc(results.edgeFrame); colormap(gca,'hot'); hold on;
    for i = 1:length(quads)
        c = quads{i};
        plot([c(:,1);c(1,1)], [c(:,2);c(1,2)], 'g-', 'LineWidth',2);
    end
    title(sprintf('Quad Candidates (%d)', length(quads))); axis image;

    % 5 - Detected markers (blue) with IDs
    subplot(2,3,5);
    imagesc(results.edgeFrame); colormap(gca,'hot'); hold on;
    for i = 1:results.numDetected
        c = results.corners{i};
        plot([c(:,1);c(1,1)], [c(:,2);c(1,2)], 'b-','LineWidth',2);
        plot(c(:,1), c(:,2), 'bs','MarkerSize',8,'MarkerFaceColor','c');
        text(mean(c(:,1)), mean(c(:,2)), ...
            sprintf('ID %d', results.ids(i)), ...
            'Color','w','FontSize',12,'FontWeight','bold', ...
            'HorizontalAlignment','center');
    end
    title(sprintf('Detected Markers (%d)', results.numDetected)); axis image;

    % 6 - Extracted bit grid of first detected marker
    subplot(2,3,6);
    if results.numDetected > 0
        gs = size(results.bits{1},1);
        full = ones(gs+2);          % border = black = 1
        full(2:end-1,2:end-1) = results.bits{1};
        imagesc(1-full); colormap(gca,'gray'); axis image;
        title(sprintf('Bits  ID %d  conf %.2f', ...
            results.ids(1), results.confidences(1)));
    else
        title('No markers detected');
    end

    sgtitle('Event-based ArUco Marker Detection');
end
