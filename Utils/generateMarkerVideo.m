%% Generate Synthetic ArUco Marker Video for ESIM
%  Single ARUCO_MIP_36h12 marker on white background.
%  Starts centered, then progresses through increasingly challenging motions.
%
%  Motion sequence (all times approximate):
%    0-2s   : Static (centered)
%    2-5s   : Slow translate right
%    5-8s   : Slow translate left (back)
%    8-11s  : Slow translate down
%    11-14s : Slow translate up (back to center)
%    14-18s : Diagonal movement (bottom-right then back)
%    18-22s : Faster horizontal oscillation
%    22-26s : Zoom in (marker grows)
%    26-30s : Zoom out (marker shrinks back)
%    30-34s : Slow rotation (CW)
%    34-38s : Faster rotation (CCW)
%    38-42s : Translate + rotate combined
%    42-46s : Fast diagonal + zoom
%    46-50s : Fast rotation + zoom + translate (hardest)
%
%  Output: ../Data/Synthetic/marker_motion_test.avi (1000 FPS)

clear; close all; clc;

%% ---- Settings --------------------------------------------------------------
frameW     = 320;
frameH     = 240;
fps        = 1000;
bgColor    = 255;            % white background

markerID   = 3;              % which marker to render
baseMarkerSize = 80;        % marker size in pixels at scale=1
cellPx     = baseMarkerSize / 8;

% Output
outDir = '../Data/Synthetic';
if ~exist(outDir, 'dir'), mkdir(outDir); end
videoFile = fullfile(outDir, 'marker_motion_test.avi');

%% ---- Build marker image ----------------------------------------------------
fprintf('Building marker ID %d...\n', markerID);
dictCodes = buildDictionary36h12();
baseMarkerImg = renderMarker(markerID, dictCodes, cellPx);
fprintf('  Marker image: %dx%d px\n', size(baseMarkerImg));

%% ---- Define motion sequence ------------------------------------------------
% Each segment: [startTime, endTime, description]
% Motion is computed per-frame based on which segment we're in.

segments = {
%   tStart  tEnd   name
    0       2      'static'
    2       5      'slow_right'
    5       8      'slow_left'
    8       11     'slow_down'
    11      14     'slow_up'
    14      18     'diagonal'
    18      22     'fast_oscillate_x'
    22      26     'zoom_in'
    26      30     'zoom_out'
    30      34     'slow_rotate'
    34      38     'fast_rotate'
    38      42     'translate_rotate'
    42      46     'fast_diagonal_zoom'
    46      50     'everything'
};

totalDuration = segments{end, 2};
nFrames = round(fps * totalDuration);

fprintf('Total duration: %.0f s  |  Frames: %d  |  FPS: %d\n', ...
    totalDuration, nFrames, fps);

%% ---- Create frames output directory ----------------------------------------
framesDir = fullfile(outDir, 'frames');
if ~exist(framesDir, 'dir'), mkdir(framesDir); end

%% ---- Generate video --------------------------------------------------------
vw = VideoWriter(videoFile, 'Grayscale AVI');
vw.FrameRate = fps;
open(vw);

% Open images.csv inside the frames folder (ESIM format)
csvFid = fopen(fullfile(framesDir, 'images.csv'), 'w');

fprintf('Generating frames...\n');

% Center position
cx0 = frameW / 2;
cy0 = frameH / 2;

for f = 1:nFrames
    t = (f-1) / fps;

    % ---- Determine current segment ----
    segIdx = 1;
    for si = 1:size(segments, 1)
        if t >= segments{si, 1} && t < segments{si, 2}
            segIdx = si;
            break;
        end
    end
    tStart = segments{segIdx, 1};
    tEnd   = segments{segIdx, 2};
    segDur = tEnd - tStart;
    tLocal = (t - tStart) / segDur;   % 0..1 within segment
    segName = segments{segIdx, 3};

    % ---- Compute motion parameters: dx, dy (from center), angle, scale ----
    % All amplitudes are relative to frame size so marker stays visible.
    % maxShift: how far the marker center can move from frame center
    % (accounting for marker size is done later by clamping)
    maxShiftX = frameW/2 - baseMarkerSize/2 - 10;  % leave 10px margin
    maxShiftY = frameH/2 - baseMarkerSize/2 - 10;

    switch segName
        case 'static'
            dx = 0; dy = 0; angle = 0; sc = 1.0;

        case 'slow_right'
            dx = lerp(0, maxShiftX*0.8, tLocal);
            dy = 0; angle = 0; sc = 1.0;

        case 'slow_left'
            dx = lerp(maxShiftX*0.8, -maxShiftX*0.8, tLocal);
            dy = 0; angle = 0; sc = 1.0;

        case 'slow_down'
            dx = lerp(-maxShiftX*0.8, 0, tLocal);
            dy = lerp(0, maxShiftY*0.8, tLocal);
            angle = 0; sc = 1.0;

        case 'slow_up'
            dx = 0;
            dy = lerp(maxShiftY*0.8, 0, tLocal);
            angle = 0; sc = 1.0;

        case 'diagonal'
            if tLocal < 0.5
                t2 = tLocal * 2;
                dx = lerp(0, maxShiftX*0.7, t2);
                dy = lerp(0, maxShiftY*0.7, t2);
            else
                t2 = (tLocal - 0.5) * 2;
                dx = lerp(maxShiftX*0.7, 0, t2);
                dy = lerp(maxShiftY*0.7, 0, t2);
            end
            angle = 0; sc = 1.0;

        case 'fast_oscillate_x'
            dx = maxShiftX*0.7 * sin(tLocal * 3 * 2 * pi);
            dy = 0; angle = 0; sc = 1.0;

        case 'zoom_in'
            dx = 0; dy = 0; angle = 0;
            sc = lerp(1.0, 2.0, tLocal);

        case 'zoom_out'
            dx = 0; dy = 0; angle = 0;
            sc = lerp(2.0, 0.6, tLocal);

        case 'slow_rotate'
            dx = 0; dy = 0;
            sc = lerp(0.6, 1.0, min(tLocal*2, 1));
            angle = lerp(0, 90, tLocal);

        case 'fast_rotate'
            dx = 0; dy = 0; sc = 1.0;
            angle = lerp(90, -180, tLocal);

        case 'translate_rotate'
            circleR = min(maxShiftX, maxShiftY) * 0.4;
            circleAngle = tLocal * 2 * pi;
            dx = circleR * cos(circleAngle);
            dy = circleR * sin(circleAngle);
            angle = tLocal * 360;
            sc = 1.0;

        case 'fast_diagonal_zoom'
            dx = maxShiftX*0.4 * sin(tLocal * 4 * 2 * pi);
            dy = maxShiftY*0.4 * cos(tLocal * 3 * 2 * pi);
            angle = 0;
            sc = 1.0 + 0.4 * sin(tLocal * 2 * 2 * pi);

        case 'everything'
            dx = maxShiftX*0.3 * sin(tLocal * 5 * 2 * pi);
            dy = maxShiftY*0.3 * cos(tLocal * 4 * 2 * pi);
            angle = tLocal * 720;
            sc = 1.0 + 0.3 * sin(tLocal * 3 * 2 * pi);

        otherwise
            dx = 0; dy = 0; angle = 0; sc = 1.0;
    end

    % ---- Render marker image (scale + rotate) ----
    scaledSize = round(baseMarkerSize * sc);
    if scaledSize < 8, scaledSize = 8; end
    scaledImg = imresize(baseMarkerImg, [scaledSize, scaledSize], 'nearest');

    if angle ~= 0
        scaledImg = imrotateWhiteBg(scaledImg, -angle, bgColor);
    end

    % ---- Clamp dx, dy so marker stays fully inside the frame ----
    [mH, mW] = size(scaledImg);
    % Marker top-left = (cx0 + dx - mW/2, cy0 + dy - mH/2)
    % Must satisfy: px >= 1  and  px + mW - 1 <= frameW
    %               py >= 1  and  py + mH - 1 <= frameH
    maxDx = frameW - mW/2 - cx0;    % rightmost dx
    minDx = mW/2 - cx0 + 1;         % leftmost dx
    maxDy = frameH - mH/2 - cy0;
    minDy = mH/2 - cy0 + 1;
    dx = max(minDx, min(maxDx, dx));
    dy = max(minDy, min(maxDy, dy));

    % ---- Place on frame ----
    frame = uint8(bgColor * ones(frameH, frameW));
    px = round(cx0 + dx - mW/2);
    py = round(cy0 + dy - mH/2);
    frame = placeMarker(frame, scaledImg, px, py, bgColor);

    % ---- Write video frame ----
    writeVideo(vw, frame);

    % ---- Write image file (ESIM naming: 10-digit zero-padded) ----
    frameName = sprintf('frames_%010d.png', f-1);
    imwrite(frame, fullfile(framesDir, frameName));

    % ---- Write images.csv line: timestamp_nanoseconds,filename ----
    tsNano = round(t * 1e9);
    fprintf(csvFid, '%d,%s\n', tsNano, frameName);

    % Progress
    if mod(f, 5000) == 0 || f == nFrames
        fprintf('  Frame %d / %d  (t=%.1fs, seg: %s)\n', ...
            f, nFrames, t, segName);
    end
end

close(vw);
fclose(csvFid);

fprintf('\nDone!\n');
fprintf('Video: %s\n', videoFile);
fprintf('Frames: %s/frames_XXXXXXXXXX.png\n', framesDir);
fprintf('CSV:    %s/images.csv\n', framesDir);
fprintf('  %d frames, %d FPS, %.0f seconds\n', nFrames, fps, totalDuration);

%% ---- Show preview of key moments ------------------------------------------
figure('Name', 'Video Preview', 'Position', [50 50 1800 600]);
previewTimes = [0, 3, 7, 10, 16, 20, 24, 28, 32, 36, 40, 44, 48, 49.5];
nPrev = length(previewTimes);
nCols = 7;
nRows = ceil(nPrev / nCols);

vid = VideoReader(videoFile);
for pi = 1:nPrev
    tPrev = previewTimes(pi);
    fIdx = round(tPrev * fps) + 1;
    fIdx = min(fIdx, nFrames);
    vid.CurrentTime = (fIdx-1) / fps;
    fr = readFrame(vid);

    subplot(nRows, nCols, pi);
    imshow(fr); axis image;

    % Find segment name
    sn = 'static';
    for si = 1:size(segments, 1)
        if tPrev >= segments{si,1} && tPrev < segments{si,2}
            sn = segments{si,3};
            break;
        end
    end
    title(sprintf('t=%.1fs\n%s', tPrev, strrep(sn,'_',' ')), 'FontSize', 8);
end
sgtitle('Video Preview — Key Moments');


%% =========================================================================
%  LOCAL FUNCTIONS
%  =========================================================================

%% ---- lerp: linear interpolation -------------------------------------------
function v = lerp(a, b, t)
    v = a + (b - a) * t;
end

%% ---- imrotateWhiteBg: rotate image with white background ------------------
function out = imrotateWhiteBg(img, angleDeg, bgVal)
    % imrotate with 'loose' to avoid clipping, fill with bgVal
    out = imrotate(img, angleDeg, 'bilinear', 'loose');
    % imrotate fills new pixels with 0 (black). Replace with bgVal.
    % Find the rotated mask
    mask = imrotate(ones(size(img)), angleDeg, 'bilinear', 'loose');
    out(mask < 0.5) = bgVal;
end

%% ---- renderMarker: create marker image from dictionary code ----------------
function img = renderMarker(markerID, dictCodes, cellPx)
    code = dictCodes(markerID + 1);
    innerGrid = zeros(6, 6);
    for bit = 0:35
        bitPos = 35 - bit;
        r = floor(bit / 6) + 1;
        c = mod(bit, 6) + 1;
        innerGrid(r, c) = bitand(bitshift(code, -bitPos), uint64(1));
    end
    fullGrid = zeros(8, 8);
    fullGrid(2:7, 2:7) = innerGrid;
    cellPx = round(cellPx);
    img = uint8(zeros(8 * cellPx, 8 * cellPx));
    for r = 1:8
        for c = 1:8
            rS = (r-1)*cellPx + 1;  rE = r*cellPx;
            cS = (c-1)*cellPx + 1;  cE = c*cellPx;
            if fullGrid(r, c) > 0
                img(rS:rE, cS:cE) = 255;
            else
                img(rS:rE, cS:cE) = 0;
            end
        end
    end
end

%% ---- placeMarker: paste marker onto frame with bounds checking -------------
function frame = placeMarker(frame, mImg, px, py, bgColor)
    [fH, fW] = size(frame);
    [mH, mW] = size(mImg);
    srcR1 = 1; srcR2 = mH; srcC1 = 1; srcC2 = mW;
    dstR1 = py; dstR2 = py+mH-1; dstC1 = px; dstC2 = px+mW-1;
    if dstR1 < 1, srcR1 = srcR1+(1-dstR1); dstR1 = 1; end
    if dstC1 < 1, srcC1 = srcC1+(1-dstC1); dstC1 = 1; end
    if dstR2 > fH, srcR2 = srcR2-(dstR2-fH); dstR2 = fH; end
    if dstC2 > fW, srcC2 = srcC2-(dstC2-fW); dstC2 = fW; end
    if dstR1 > dstR2 || dstC1 > dstC2, return; end
    if srcR1 > srcR2 || srcC1 > srcC2, return; end
    frame(dstR1:dstR2, dstC1:dstC2) = mImg(srcR1:srcR2, srcC1:srcC2);
end

%% ---- buildDictionary36h12 --------------------------------------------------
function codes = buildDictionary36h12()
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
end
