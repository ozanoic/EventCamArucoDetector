%% Sarmadi et al. (2021) — Complete MATLAB Port (with ID decoding)
%  "Detection of Binary Square Fiducial Markers Using an Event Camera"
%
%  Full pipeline per packet (~10 ms batch):
%    1. Separate ON / OFF events → min-timestamp images
%    2. Preprocess (majority filter + masked Gaussian blur)
%    3. Detect line segments (Hough — he used LSD)
%    4. Age correction (move segments to age=0.5)
%    5. Match ON↔OFF segment pairs → marker candidates
%    6. Unwarp candidates to 160×160 standard square
%    7. Decode: Gaussian convolution → threshold → binary code
%    8. Dictionary lookup (ARUCO_MIP_36h12, 250 markers)
%    9. Erase everything → next packet
%
%  DVS128 sensor: 128×128 pixels.

clear; close all; clc;

%% ---- User settings ---------------------------------------------------------
visualizeAll = 0;   % true: visualize every packet | false: only on detection

% Choose dataset: 'sarmadi' or 'esim'
dataset = 'esim';

%% ---- Load events -----------------------------------------------------------
switch dataset
    case 'sarmadi'
        binFile = '../Data/Sarmadi/side2side/side2side/packets.bin';
        fprintf('Loading Sarmadi .bin data...\n');
        events = convertSarmadiBin(binFile);
        sensorSize = [128, 128];
        dataLabel = binFile;

    case 'esim'
        matFile = '../Data/Synthetic/MovingCam/moving_events/moving_events.mat';
        fprintf('Loading ESIM .mat data...\n');
        tmp = load(matFile, 'events');
        events = tmp.events;
        % polarity: 0→-1 for internal use (pol>0 = ON, pol<=0 = OFF)
        events(events(:,3) == 0, 3) = -1;
        sensorSize = [240, 320];   % [H, W] from calib.yaml
        dataLabel = matFile;

    otherwise
        error('Unknown dataset: %s', dataset);
end

numEvents = size(events, 1);
H = sensorSize(1);
W = sensorSize(2);

%% ---- Parameters (from his code / paper Table 1) ---------------------------
% Scale parameters for sensor resolution
scaleFactor = max(H, W) / 128;   % relative to DVS128

packetDt    = 10000;              % 10 ms = 10000 µs per packet
minSegLen   = round(25 * scaleFactor);  % scale with resolution
fillGap     = round(5 * scaleFactor);   % Hough FillGap, scale with resolution
kernelSize  = 3;                  % Gaussian kernel for preprocessing
kernelSigma = 0.3*((kernelSize-1)*0.5-1) + 0.8;
numPeaks    = round(50 * scaleFactor);  % more peaks for larger image
angleTolDeg = 30;                 % angle tolerance for segment matching
angleTolRad = angleTolDeg * pi / 180;

% Marker decoding parameters (from markercandidate.h)
cellSize    = 20;                 % pixels per cell in unwarped image
numCells    = 8;                  % 8×8 grid (6×6 inner + border)
sideSize    = cellSize * numCells;% 160 pixels
threshOn    = 55;                 % threshold percentage for ON scores
threshOff   = 55;                 % threshold percentage for OFF scores
codeSize    = 6;                  % inner grid = 6×6

% Standard marker corners (0-indexed like OpenCV)
markerCoords = [0 0; sideSize-1 0; sideSize-1 sideSize-1; 0 sideSize-1];

% Build dictionary
dictionary = buildDictionary();
fprintf('Dictionary loaded: %d markers\n', dictionary.Count);

%% ---- Group events into packets by time ------------------------------------
tAll = events(:, 4);
tMin = min(tAll);
tMax = max(tAll);
packetEdges = tMin : packetDt : (tMax + packetDt);
numPackets = length(packetEdges) - 1;
[~, ~, packetIdx] = histcounts(tAll, packetEdges);

fprintf('Total events: %d  |  Packets: %d  (%.0f µs each)\n', ...
    numEvents, numPackets, packetDt);

%% ---- Prepare figure --------------------------------------------------------
hFig = figure('Name','Sarmadi Method — Full Pipeline', ...
              'Position',[50 50 1600 800]);

totalDetections = 0;

%% ---- Open output file ------------------------------------------------------
logFid = fopen('sarmadi_output.txt', 'w');
fprintf(logFid, 'Sarmadi et al. (2021) — Detection Log\n');
fprintf(logFid, 'Dataset: %s\n', dataLabel);
fprintf(logFid, 'Total events: %d  |  Packets: %d  |  packetDt: %d us\n', numEvents, numPackets, packetDt);
fprintf(logFid, '----------------------------------------------------------------------\n');
fprintf(logFid, '%6s  %6s  %5s  %5s  %5s  %4s  %s\n', ...
    'Packet', 'Events', 'OnSeg', 'OffSeg', 'Cand', 'Det', 'MarkerIDs');
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
    pol = pktEvents(:,3);           % -1 or +1
    ts  = pktEvents(:,4);

    % ---- Timestamp statistics ----
    ts_min = min(ts);
    ts_max = max(ts);
    ts_interval = ts_max - ts_min;
    if ts_interval == 0, continue; end

    ratio = (ts - ts_min) / ts_interval;

    % ---- Build ON / OFF images (his fill_in_frames) ----
    on_img   = zeros(H, W);   on_mask  = zeros(H, W);   on_min  = inf(H, W);
    off_img  = zeros(H, W);   off_mask = zeros(H, W);   off_min = inf(H, W);

    for e = 1:nEvt
        r = py(e);  c = px(e);
        if r < 1 || r > H || c < 1 || c > W, continue; end
        if pol(e) > 0
            on_img(r,c)  = on_img(r,c) + ratio(e);
            on_mask(r,c) = on_mask(r,c) + 1;
            on_min(r,c)  = min(on_min(r,c), ratio(e));
        else
            off_img(r,c)  = off_img(r,c) + ratio(e);
            off_mask(r,c) = off_mask(r,c) + 1;
            off_min(r,c)  = min(off_min(r,c), ratio(e));
        end
    end

    % ---- Compute on_one_min / off_one_min ----
    on_one_min_img  = zeros(H,W);   on_one_min_mask = false(H,W);
    on_min_img      = zeros(H,W);   on_min_mask     = false(H,W);
    off_one_min_img = zeros(H,W);   off_one_min_mask= false(H,W);
    off_min_img     = zeros(H,W);   off_min_mask    = false(H,W);

    for r = 1:H
        for c = 1:W
            if on_mask(r,c) > 0
                on_one_min_mask(r,c) = true;
                on_min_mask(r,c)     = true;
                on_one_min_img(r,c)  = 1.0 - on_min(r,c);
                on_min_img(r,c)      = on_min(r,c);
                on_img(r,c)          = on_img(r,c) / on_mask(r,c);
            end
            if off_mask(r,c) > 0
                off_one_min_mask(r,c) = true;
                off_min_mask(r,c)     = true;
                off_one_min_img(r,c)  = 1.0 - off_min(r,c);
                off_min_img(r,c)      = off_min(r,c);
                off_img(r,c)          = off_img(r,c) / off_mask(r,c);
            end
        end
    end

    % ---- Preprocess ----
    [on_one_min_img, on_one_min_mask]   = preprocessMajority(on_one_min_img,  on_one_min_mask, H, W);
    [off_one_min_img, off_one_min_mask] = preprocessMajority(off_one_min_img, off_one_min_mask, H, W);
    on_one_min_img  = maskedGaussBlur(on_one_min_img,  double(on_one_min_mask),  kernelSize, kernelSigma);
    off_one_min_img = maskedGaussBlur(off_one_min_img, double(off_one_min_mask), kernelSize, kernelSigma);

    % ---- Detect line segments ----
    on_segs  = detectLineSegments(on_one_min_img,  on_one_min_mask,  minSegLen, numPeaks, fillGap);
    off_segs = detectLineSegments(off_one_min_img, off_one_min_mask, minSegLen, numPeaks, fillGap);

    % ---- Age correction ----
    on_segs  = moveSegments(on_min_img,  on_min_mask,  on_segs,  0.5, H, W);
    off_segs = moveSegments(off_min_img, off_min_mask, off_segs, 0.5, H, W);

    % ---- Match ON↔OFF pairs → candidates ----
    [candidates, candOnSegs, candOffSegs] = findMarkerCandidates( ...
        on_segs, off_segs, minSegLen, angleTolRad);

    % ---- Unwarp, decode, dictionary lookup ----
    on_min_u8  = uint8(on_min_img * 255);
    off_min_u8 = uint8(off_min_img * 255);

    detectedMarkers = [];   % struct array: .id, .corners, .codeImg

    for ci = 1:length(candidates)
        corners = candidates{ci};

        % Unwarp both ON and OFF min images
        [warpedOn, warpedOff, tform] = unwarpCandidate( ...
            on_min_u8, off_min_u8, corners, markerCoords, sideSize);
        if isempty(warpedOn), continue; end

        % Decode the candidate
        [codeImg, codes4rot] = decodeCandidate( ...
            warpedOn, warpedOff, cellSize, numCells, codeSize, threshOn, threshOff);

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
            det.id       = markerIdx;
            det.corners  = corners;
            det.codeImg  = codeImg;
            det.warpedOn = warpedOn;
            det.warpedOff= warpedOff;
            det.orientation = orientation;
            detectedMarkers = [detectedMarkers, det]; %#ok<AGROW>
            totalDetections = totalDetections + 1;
        end
    end

    % ---- Log ----
    nDet = length(detectedMarkers);
    nCand = length(candidates);

    % Write to output file (every packet)
    idStr = '';
    if nDet > 0
        ids = arrayfun(@(d) d.id, detectedMarkers);
        idStr = strjoin(arrayfun(@(x) sprintf('%d', x), ids, 'UniformOutput', false), ',');
    end
    fprintf(logFid, '%6d  %6d  %5d  %5d  %5d  %4d  %s\n', ...
        p, nEvt, length(on_segs), length(off_segs), nCand, nDet, idStr);

    % Console log
    if mod(p, 100) == 0 || nDet > 0
        fprintf('[pkt %4d/%d]  events:%d  on:%d  off:%d  cand:%d  detected:%d', ...
            p, numPackets, nEvt, length(on_segs), length(off_segs), ...
            nCand, nDet);
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

    % (a) ON events preprocessed
    subplot(2,4,1);
    imagesc(on_one_min_img, [0 1]); colormap(gca,'hot'); axis image;
    title(sprintf('ON (1-min) | pkt %d', p));

    % (b) OFF events preprocessed
    subplot(2,4,2);
    imagesc(off_one_min_img, [0 1]); colormap(gca,'hot'); axis image;
    title(sprintf('OFF (1-min) | %d evts', nEvt));

    % (c) Combined + segments
    subplot(2,4,3);
    combRGB = zeros(H, W, 3);
    combRGB(:,:,1) = on_one_min_img;
    combRGB(:,:,3) = off_one_min_img;
    imshow(combRGB * 2); hold on;
    for si = 1:length(on_segs)
        s = on_segs(si);
        plot([s.p1(1) s.p2(1)], [s.p1(2) s.p2(2)], 'r-', 'LineWidth', 2);
    end
    for si = 1:length(off_segs)
        s = off_segs(si);
        plot([s.p1(1) s.p2(1)], [s.p1(2) s.p2(2)], 'b-', 'LineWidth', 2);
    end
    title(sprintf('Segs: ON=%d OFF=%d', length(on_segs), length(off_segs)));
    axis image; hold off;

    % (d) Candidates + detections
    subplot(2,4,4);
    combRGB2 = combRGB;
    imshow(combRGB2 * 2); hold on;
    % Draw all candidates in gray
    for ci = 1:length(candidates)
        cc = candidates{ci};
        plot([cc(:,1); cc(1,1)], [cc(:,2); cc(1,2)], '-', ...
            'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    end
    % Draw detected markers in green with ID
    for di = 1:nDet
        cc = detectedMarkers(di).corners;
        plot([cc(:,1); cc(1,1)], [cc(:,2); cc(1,2)], 'g-', 'LineWidth', 2.5);
        plot(cc(:,1), cc(:,2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        text(mean(cc(:,1)), mean(cc(:,2)), sprintf('ID:%d', detectedMarkers(di).id), ...
            'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'BackgroundColor', 'k');
    end
    title(sprintf('Cand:%d  Detected:%d', length(candidates), nDet));
    axis image; hold off;

    % (e-h) Show decoded marker details if any detected
    if nDet > 0
        det1 = detectedMarkers(1);

        subplot(2,4,5);
        imshow(det1.warpedOn); axis image;
        title(sprintf('Unwarped ON (ID:%d)', det1.id));

        subplot(2,4,6);
        imshow(det1.warpedOff); axis image;
        title('Unwarped OFF');

        subplot(2,4,7);
        imagesc(det1.codeImg); colormap(gca, 'gray'); axis image;
        title('Decoded grid');
        % Draw grid lines
        hold on;
        for gi = 0:numCells
            plot([gi+0.5, gi+0.5], [0.5, numCells+0.5], 'r-', 'LineWidth', 0.5);
            plot([0.5, numCells+0.5], [gi+0.5, gi+0.5], 'r-', 'LineWidth', 0.5);
        end
        hold off;

        subplot(2,4,8);
        % Show reconstructed marker image
        markerImg = zeros(numCells);
        markerImg(1,:) = 0; markerImg(end,:) = 0;
        markerImg(:,1) = 0; markerImg(:,end) = 0;
        markerImg(2:end-1, 2:end-1) = det1.codeImg(2:end-1, 2:end-1);
        imagesc(1 - markerImg); colormap(gca, 'gray'); axis image;
        title(sprintf('Marker ID: %d', det1.id));
        hold on;
        for gi = 0:numCells
            plot([gi+0.5, gi+0.5], [0.5, numCells+0.5], 'r-', 'LineWidth', 0.5);
            plot([0.5, numCells+0.5], [gi+0.5, gi+0.5], 'r-', 'LineWidth', 0.5);
        end
        hold off;
    else
        subplot(2,4,5); cla; title('No detection');
        subplot(2,4,6); cla; title('No detection');
        subplot(2,4,7); cla; title('No detection');
        subplot(2,4,8); cla; title('No detection');
    end

    sgtitle(sprintf('Sarmadi Full Pipeline  |  Packet %d/%d  |  %d events  |  Total detections: %d', ...
        p, numPackets, nEvt, totalDetections));
    drawnow;
    end  % if visualizeAll || nDet > 0
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
fprintf('Debug log saved to: sarmadi_output.txt\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%  =========================================================================

%% ---- buildDictionary: ARUCO_MIP_36h12 (250 markers) -----------------------
function dict = buildDictionary()
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

%% ---- preprocessMajority ----------------------------------------------------
function [imgOut, maskOut] = preprocessMajority(img, mask, H, W)
    imgOut  = img;
    maskOut = mask;
    for y = 1:H
        yn_lo = max(1, y-1);  yn_hi = min(H, y+1);
        for x = 1:W
            xn_lo = max(1, x-1);  xn_hi = min(W, x+1);
            halfTotal = ((yn_hi-yn_lo+1)*(xn_hi-xn_lo+1) - 1) / 2;
            nValid = 0;  valSum = 0;
            for yn = yn_lo:yn_hi
                for xn = xn_lo:xn_hi
                    if yn == y && xn == x, continue; end
                    if mask(yn, xn)
                        nValid = nValid + 1;
                        valSum = valSum + img(yn, xn);
                    end
                end
            end
            if nValid > halfTotal && ~mask(y, x)
                imgOut(y, x) = valSum / nValid;
                maskOut(y, x) = true;
            elseif nValid < halfTotal && mask(y, x)
                imgOut(y, x) = 0;
                maskOut(y, x) = false;
            end
        end
    end
end

%% ---- maskedGaussBlur --------------------------------------------------------
function out = maskedGaussBlur(img, maskDbl, ksize, sigma)
    imgBlur  = imgaussfilt(img,     sigma, 'FilterSize', ksize, 'Padding', 0);
    maskBlur = imgaussfilt(maskDbl, sigma, 'FilterSize', ksize, 'Padding', 0);
    out = zeros(size(img));            % C++ zeros everything first
    valid = maskBlur > 0;
    out(valid) = imgBlur(valid) ./ maskBlur(valid);
end

%% ---- detectLineSegments -----------------------------------------------------
function segs = detectLineSegments(img, mask, minLen, numPeaks, fillGap)
    segs = struct('p1',{}, 'p2',{}, 'len',{}, 'theta',{}, 'age',{});

    % C++ uses LSD on the float64 image scaled to 255.
    % MATLAB has no built-in LSD; use skeleton of the mask as edge map.
    % This matches the C++ candidate_detector_hough path which uses
    % ximgproc::thinning on the binary mask before Hough.
    if ~any(mask(:)), return; end
    edges = bwmorph(mask > 0, 'thin', Inf);

    % Hough on the edge image
    [Ht, theta, rho] = hough(edges, 'RhoResolution', 1, 'Theta', -90:0.5:89.5);
    maxH = max(Ht(:));
    if maxH == 0, return; end
    peaks = houghpeaks(Ht, numPeaks, 'Threshold', ceil(0.10 * maxH));
    if isempty(peaks), return; end
    hLines = houghlines(edges, theta, rho, peaks, 'FillGap', fillGap, 'MinLength', minLen);
    for k = 1:length(hLines)
        p1 = hLines(k).point1;  p2 = hLines(k).point2;
        dxy = p2 - p1;  len = norm(dxy);
        if len < minLen, continue; end
        th  = atan2(dxy(2), dxy(1));
        age = calcLineAge(img, mask, p1, p2);
        segs(end+1) = struct('p1',p1, 'p2',p2, 'len',len, 'theta',th, 'age',age); %#ok<AGROW>
    end
end

%% ---- calcLineAge ------------------------------------------------------------
function age = calcLineAge(img, mask, p1, p2)
    dx = p2(1)-p1(1);  dy = p2(2)-p1(2);
    steps = round(max(abs(dx), abs(dy)));
    if steps == 0, age = 0; return; end
    xStep = dx/steps;  yStep = dy/steps;
    total = 0;  count = 0;
    [H, W] = size(img);
    for s = 0:steps
        x = round(p1(1) + xStep*s);
        y = round(p1(2) + yStep*s);
        if x>=1 && x<=W && y>=1 && y<=H && mask(y,x)
            total = total + img(y,x);
            count = count + 1;
        end
    end
    if count > 0, age = total/count; else, age = 0; end
end

%% ---- moveSegments (matches C++ find_line_offset + move_segments) -----------
function segs = moveSegments(minImg, minMask, segs, goalAge, H, W)
    for i = 1:length(segs)
        ls = segs(i);
        x1 = ls.p1(1);  y1 = ls.p1(2);
        x2 = ls.p2(1);  y2 = ls.p2(2);
        x_diff = x2 - x1;
        neg_y_diff = y1 - y2;

        alphas = 0;
        betas  = calcLineAge(minImg, minMask, ls.p1, ls.p2);

        if abs(x_diff) > abs(neg_y_diff)
            % x = a*y + b  → step along y (perpendicular is mostly in y)
            a = neg_y_diff / x_diff;
            for dir = [-1, 1]
                for abs_dy = 1:200
                    dy_step = abs_dy * dir;
                    dx_step = a * dy_step;
                    nP1 = [x1+dx_step, y1+dy_step];
                    nP2 = [x2+dx_step, y2+dy_step];
                    if round(nP1(1))<1||round(nP1(1))>W||round(nP1(2))<1||round(nP1(2))>H||...
                       round(nP2(1))<1||round(nP2(1))>W||round(nP2(2))<1||round(nP2(2))>H
                        break;
                    end
                    if ~isEventful(minMask, nP1, nP2), break; end
                    alphas(end+1) = dy_step; %#ok<AGROW>
                    betas(end+1)  = calcLineAge(minImg, minMask, nP1, nP2); %#ok<AGROW>
                end
            end
            if length(alphas) >= 2
                [c_lr, d_lr] = linReg(alphas, betas);
                if abs(c_lr) > 1e-6
                    y_off = (goalAge - d_lr) / c_lr;
                    x_off = a * y_off;
                else
                    x_off = 0;  y_off = 0;
                end
            else
                x_off = 0;  y_off = 0;
            end
        else
            % y = a*x + b  → step along x
            if abs(neg_y_diff) < 1e-10, continue; end
            a = x_diff / neg_y_diff;
            for dir = [-1, 1]
                for abs_dx = 1:200
                    dx_step = abs_dx * dir;
                    dy_step = a * dx_step;
                    nP1 = [x1+dx_step, y1+dy_step];
                    nP2 = [x2+dx_step, y2+dy_step];
                    if round(nP1(1))<1||round(nP1(1))>W||round(nP1(2))<1||round(nP1(2))>H||...
                       round(nP2(1))<1||round(nP2(1))>W||round(nP2(2))<1||round(nP2(2))>H
                        break;
                    end
                    if ~isEventful(minMask, nP1, nP2), break; end
                    alphas(end+1) = dx_step; %#ok<AGROW>
                    betas(end+1)  = calcLineAge(minImg, minMask, nP1, nP2); %#ok<AGROW>
                end
            end
            if length(alphas) >= 2
                [c_lr, d_lr] = linReg(alphas, betas);
                if abs(c_lr) > 1e-6
                    x_off = (goalAge - d_lr) / c_lr;
                    y_off = a * x_off;
                else
                    x_off = 0;  y_off = 0;
                end
            else
                x_off = 0;  y_off = 0;
            end
        end

        % Clamp offsets to keep segment in bounds (matching C++ move_segments)
        if x_off >= 0
            x_off = min(x_off, min(W - x2, W - x1));
        else
            x_off = max(x_off, max(-x2+1, -x1+1));
        end
        if y_off >= 0
            y_off = min(y_off, min(H - y2, H - y1));
        else
            y_off = max(y_off, max(-y2+1, -y1+1));
        end

        segs(i).p1 = [x1 + x_off, y1 + y_off];
        segs(i).p2 = [x2 + x_off, y2 + y_off];
        segs(i).age = calcLineAge(minImg, minMask, segs(i).p1, segs(i).p2);
    end
end

%% ---- isEventful -------------------------------------------------------------
function ok = isEventful(mask, p1, p2)
    dx=p2(1)-p1(1); dy=p2(2)-p1(2);
    steps=round(max(abs(dx),abs(dy)));
    if steps==0, ok=false; return; end
    xS=dx/steps; yS=dy/steps;
    [H,W]=size(mask); count=0;
    for s=0:steps
        x=round(p1(1)+xS*s); y=round(p1(2)+yS*s);
        if x>=1&&x<=W&&y>=1&&y<=H&&mask(y,x), count=count+1; end
    end
    ok = count > steps/2;
end

%% ---- linReg -----------------------------------------------------------------
function [c, d] = linReg(alphas, betas)
    n=length(alphas); sA=sum(alphas); sB=sum(betas);
    sAB=sum(alphas.*betas); sAA=sum(alphas.*alphas);
    den=n*sAA-sA*sA;
    if abs(den)<1e-10, c=0; d=sB/n; else
        c=(n*sAB-sA*sB)/den; d=(sB-c*sA)/n;
    end
end

%% ---- findMarkerCandidates ---------------------------------------------------
function [candidates, candOn, candOff] = findMarkerCandidates(on_segs, off_segs, minLen, angleTolRad)
    candidates = {};  candOn = {};  candOff = {};
    if isempty(on_segs) || isempty(off_segs), return; end
    for i = 1:length(on_segs)
        on = on_segs(i);
        if on.len < minLen, continue; end
        for j = 1:length(off_segs)
            off = off_segs(j);
            if off.len < minLen, continue; end
            if off.len*2 < on.len || off.len > on.len*2, continue; end
            dth = abs(on.theta - off.theta);
            if dth > pi, dth = 2*pi-dth; end
            if dth > pi/2, dth = pi-dth; end
            if dth > angleTolRad, continue; end
            if ~(projectsIn(on.p1,off)||projectsIn(on.p2,off)||...
                 projectsIn(off.p1,on)||projectsIn(off.p2,on))
                continue;
            end
            [corners, ~] = orderCornersSarmadi(on, off);
            candidates{end+1} = corners; %#ok<AGROW>
            candOn{end+1} = on; candOff{end+1} = off; %#ok<AGROW>
        end
    end
    if length(candidates) > 1
        candidates = deduplicateQuads(candidates);
    end
end

%% ---- projectsIn -------------------------------------------------------------
function ok = projectsIn(p, ls)
    v=ls.p2-ls.p1; w=p-ls.p1;
    t=dot(w,v)/dot(v,v);
    ok=(t>=0)&&(t<=1);
end

%% ---- orderCornersSarmadi: exact replica of his MarkerCandidate constructor --
function [corners, offSide] = orderCornersSarmadi(onSeg, offSeg)
    pts = [onSeg.p1; onSeg.p2; offSeg.p1; offSeg.p2];  % 4×2

    % Sort by x
    [~, xOrd] = sort(pts(:,1));
    % Left pair (indices 1,2 in sorted order), right pair (3,4)
    xi1 = xOrd(1);  xi2 = xOrd(2);  xi3 = xOrd(3);  xi4 = xOrd(4);

    % Within left pair: lower y first (top)
    if pts(xi1,2) > pts(xi2,2)
        tmp=xi1; xi1=xi2; xi2=tmp;
    end
    % Within right pair: lower y first
    if pts(xi3,2) > pts(xi4,2)
        tmp=xi3; xi3=xi4; xi4=tmp;
    end

    % input_indices: [0]=TL, [1]=TR, [2]=BR, [3]=BL  (0-indexed in C++)
    inputIdx = [xi1, xi3, xi4, xi2];  % TL TR BR BL

    % output_indices: which output slot does each original point go to?
    outputIdx = zeros(1,4);
    outputIdx(xi1) = 0;  % TL
    outputIdx(xi2) = 3;  % BL
    outputIdx(xi3) = 1;  % TR
    outputIdx(xi4) = 2;  % BR

    % Determine off_side: which side has the OFF segment endpoints (indices 3,4)
    offP1out = outputIdx(3);  % offSeg.p1 → which corner slot
    offP2out = outputIdx(4);  % offSeg.p2 → which corner slot
    op1 = min(offP1out, offP2out);
    op2 = max(offP1out, offP2out);
    if op1 == 0 && op2 == 3
        offSide = 3;
    else
        offSide = op1;
    end

    % Rotate corners so OFF side is on the left (his p[i] = points[input_indices[(i+off_side+1)%4]])
    corners = zeros(4, 2);
    for i = 0:3
        srcIdx = mod(i + offSide + 1, 4) + 1;  % 1-indexed
        corners(i+1, :) = pts(inputIdx(srcIdx), :);
    end
end

%% ---- unwarpCandidate --------------------------------------------------------
function [warpedOn, warpedOff, tform] = unwarpCandidate(onMinU8, offMinU8, corners, markerCoords, sideSize)
%  Perspective warp from candidate corners to standard square.
%  C++ code: H = getPerspectiveTransform(markerCoords, corners)
%            warpPerspective(image, output, H, markerSize, WARP_INVERSE_MAP)
%  OpenCV coords are 0-indexed. WARP_INVERSE_MAP means: for each output
%  pixel p_out, compute p_src = H * p_out to find the source pixel.
%  So H maps FROM markerCoords TO corners (standard square → image).
    warpedOn = [];  warpedOff = [];  tform = [];

    try
        % C++ getPerspectiveTransform(markerCoords, corners) computes H
        % that maps markerCoords → corners.  With WARP_INVERSE_MAP,
        % OpenCV applies H directly to each output pixel to get source
        % coords: p_src = H * p_out.  Output is in marker space, so
        % H(marker) = image_coords, which is marker → image.
        %
        % MATLAB imwarp(img, srcRef, tform, 'OutputView', outView):
        %   tform maps INPUT world → OUTPUT world.
        %   For each output pixel p_out, it computes tform.inverse(p_out)
        %   to get the input coord to sample.
        %
        % We want input=image, output=marker.  So tform must map
        % image → marker.  fitgeotrans(movingPts, fixedPts) maps
        % moving → fixed.  So: moving=corners (image), fixed=markerCoords.
        movingPts = corners - 1;                  % image coords, 0-indexed
        fixedPts  = markerCoords;                  % marker coords, 0-indexed

        tform = fitgeotrans(movingPts, fixedPts, 'projective');

        % Source image spatial ref: 0-indexed pixel centers at 0..W-1, 0..H-1
        % imref2d world limits must span from -0.5 to N-0.5 so that
        % pixel centers land exactly at integer coords 0, 1, ..., N-1
        [srcH, srcW] = size(onMinU8);
        srcRef = imref2d([srcH, srcW], [-0.5, srcW-0.5], [-0.5, srcH-0.5]);
        outputView = imref2d([sideSize, sideSize], [-0.5, sideSize-0.5], [-0.5, sideSize-0.5]);

        warpedOn  = imwarp(onMinU8,  srcRef, tform, 'OutputView', outputView, 'Interp', 'bilinear');
        warpedOff = imwarp(offMinU8, srcRef, tform, 'OutputView', outputView, 'Interp', 'bilinear');
    catch
        % fitgeotrans can fail with degenerate points
        warpedOn = [];  warpedOff = [];
    end
end

%% ---- decodeCandidate: his decode_candidate() --------------------------------
function [codeImg, codes4rot] = decodeCandidate(warpedOn, warpedOff, cellSize, numCells, codeSize, threshOn, threshOff)
%  Gaussian convolution decoding (paper Section III-E, Eq. 25-26)
    gaussW = cellSize;
    gaussH = cellSize;

    % 2D Gaussian kernel (his: getGaussianKernel(cellSize,-1) → outer product)
    sig = 0.3*((cellSize-1)/2.0 - 1) + 0.8;
    gk1 = fspecial('gaussian', [gaussH, 1], sig);
    gk2 = fspecial('gaussian', [1, gaussW], sig);
    gKernel = gk1 * gk2;
    gKernel = gKernel / sum(gKernel(:));  % normalize

    onScores  = zeros(numCells);
    offScores = zeros(numCells);

    wOn  = double(warpedOn);
    wOff = double(warpedOff);

    maxOnScore  = -1;
    maxOffScore = -1;

    % Compute convolution scores for inner cells
    % C++ uses 0-indexed: r=1..6, c=1..6 (inner cells)
    % C++ ROI: y_min = r*cellSize + cellSize/2 - gaussH/2
    %          x_min = c*cellSize - gaussW/2
    % MATLAB is 1-indexed, so add +1 to convert from 0-indexed pixel coords
    for r = 2:(numCells-1)       % 1-indexed cell index (= C++ r+1)
        for c = 2:(numCells-1)
            % Match C++ exactly: r_cpp = r-1, c_cpp = c-1
            r_cpp = r - 1;
            c_cpp = c - 1;
            y_min = r_cpp*cellSize + cellSize/2 - gaussH/2 + 1;  % +1 for 1-indexing
            y_max = y_min + gaussH - 1;
            x_min = c_cpp*cellSize - gaussW/2 + 1;               % +1 for 1-indexing
            x_max = x_min + gaussW - 1;

            % ON: shifted left by gaussW/4 (his code: x_min - gaussian_width/4)
            shiftOn = floor(gaussW / 4);
            x_min_on = x_min - shiftOn;
            x_max_on = x_max - shiftOn;

            % Clamp to image bounds
            x_min_on = max(1, x_min_on);  x_max_on = min(size(wOn,2), x_max_on);
            y_min    = max(1, y_min);      y_max    = min(size(wOn,1), y_max);
            x_min    = max(1, x_min);      x_max    = min(size(wOff,2), x_max);

            % Extract ROIs
            roiOn  = wOn(y_min:y_max, x_min_on:x_max_on);
            roiOff = wOff(y_min:y_max, x_min:x_max);

            % Convolve with Gaussian (element-wise multiply and sum)
            kH = min(size(roiOn,1), size(gKernel,1));
            kW_on = min(size(roiOn,2), size(gKernel,2));
            kW_off = min(size(roiOff,2), size(gKernel,2));

            if kH > 0 && kW_on > 0
                onScore = sum(sum(roiOn(1:kH, 1:kW_on) .* gKernel(1:kH, 1:kW_on)));
            else
                onScore = 0;
            end
            if kH > 0 && kW_off > 0
                offScore = sum(sum(roiOff(1:kH, 1:kW_off) .* gKernel(1:kH, 1:kW_off)));
            else
                offScore = 0;
            end

            onScores(r, c)  = onScore;
            offScores(r, c) = offScore;
            maxOnScore  = max(maxOnScore, onScore);
            maxOffScore = max(maxOffScore, offScore);
        end
    end

    % Threshold scores (paper Eq. 25)
    onBin  = zeros(numCells);
    offBin = zeros(numCells);
    if maxOnScore > 0
        for r = 2:(numCells-1)
            for c = 2:(numCells-1)
                if onScores(r,c) * 100.0 / maxOnScore > threshOn
                    onBin(r,c) = 1;
                end
            end
        end
    end
    if maxOffScore > 0
        for r = 2:(numCells-1)
            for c = 2:(numCells-1)
                if offScores(r,c) * 100.0 / maxOffScore > threshOff
                    offBin(r,c) = 1;
                end
            end
        end
    end

    % Decode each row left-to-right (paper Eq. 26)
    codeImg = zeros(numCells);
    for r = 2:(numCells-1)
        codeImg(r, 1) = 0;   % first cell assumed black (border)
        for c = 2:(numCells-1)
            if codeImg(r, c-1) == 0
                codeImg(r, c) = 0;
                if onBin(r, c) > 0
                    codeImg(r, c) = 1;
                end
            else
                codeImg(r, c) = 1;
                if offBin(r, c) > 0
                    codeImg(r, c) = 0;
                end
            end
        end
    end

    % Extract 36-bit code in 4 rotations (his code lines 181-220)
    codes4rot = zeros(1, 4, 'uint64');
    for rot = 0:3
        if mod(rot, 3) == 0
            iRange = 0:(codeSize-1);
        else
            iRange = (codeSize-1):-1:0;
        end
        if rot < 2
            jRange = 0:(codeSize-1);
        else
            jRange = (codeSize-1):-1:0;
        end

        code = uint64(0);
        for ii = iRange
            for jj = jRange
                code = bitshift(code, 1);
                if mod(rot, 2) == 0
                    code = code + uint64(codeImg(ii+2, jj+2));  % +2: 1-indexed + skip border
                else
                    code = code + uint64(codeImg(jj+2, ii+2));  % transposed
                end
            end
        end
        codes4rot(rot+1) = code;
    end
end

%% ---- deduplicateQuads -------------------------------------------------------
function merged = deduplicateQuads(allQ)
    n = length(allQ);
    centers = zeros(n,2);  sizes = zeros(n,1);
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
