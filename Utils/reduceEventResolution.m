function outPath = reduceEventResolution(matFile, scale, opts)
%REDUCEEVENTRESOLUTION  Downscale event-camera (x, y) and save / visualise.
%
%  outPath = reduceEventResolution(matFile)             % scale = 2, show = true
%  outPath = reduceEventResolution(matFile, scale)
%  outPath = reduceEventResolution(matFile, scale, opts)
%
%  Loads `events` from matFile (Nx4 -- x, y, polarity, timestamp),
%  downsamples the (x, y) coordinates by an integer factor `scale`
%  (default 2), and saves a sibling .mat file named `<name>_reduced.mat`.
%  Other variables in the source file are preserved as-is.
%
%  Optionally shows a side-by-side comparison: original sensor resolution
%  on top, reduced resolution on the bottom, several time slices across
%  the recording, each accumulating a configurable window (default 10 ms).
%
%  Inputs:
%    matFile : path to a .mat with an `events` (Nx4) variable
%    scale   : positive integer downscale factor (default 2)
%    opts    : optional struct with fields
%        .show       (default true)        show the side-by-side figure
%        .windowMs   (default 10)          accumulation window per slice (ms)
%        .numFrames  (default 6)           how many time slices to draw
%        .outPath    (default '')          explicit output path
%        .origSize   (default inferred)    [H W] of source sensor
%
%  Output:
%    outPath : path to the saved <name>_reduced.mat
%
%  Example:
%    reduceEventResolution('Data/OzanEventData_22.05.2026/1/1.mat');
%    reduceEventResolution('Data/.../1.mat', 2, struct('windowMs', 20, 'numFrames', 8));

if nargin < 2 || isempty(scale), scale = 2;           end
if nargin < 3,                   opts  = struct();    end
if ~isfield(opts,'show'),        opts.show      = true; end
if ~isfield(opts,'windowMs'),    opts.windowMs  = 10;   end
if ~isfield(opts,'numFrames'),   opts.numFrames = 6;    end
if ~isfield(opts,'outPath'),     opts.outPath   = '';   end
if ~isfield(opts,'origSize'),    opts.origSize  = [];   end

if ~isfile(matFile)
    error('reduceEventResolution: file not found: %s', matFile);
end

scale = max(1, round(scale));
if scale == 1
    warning('reduceEventResolution: scale == 1 is a no-op; nothing to reduce.');
end

%% ---- Load --------------------------------------------------------------
fprintf('Loading %s...\n', matFile);
allVars = load(matFile);                    % keep every variable in the file
if ~isfield(allVars, 'events')
    error('reduceEventResolution: %s has no `events` variable.', matFile);
end
events = allVars.events;
if size(events, 2) < 4
    error('reduceEventResolution: events must be Nx4 (x, y, polarity, t), got %dx%d.', ...
        size(events,1), size(events,2));
end
fprintf('  %d events loaded\n', size(events, 1));

%% ---- Infer original / new sensor size ---------------------------------
if ~isempty(opts.origSize)
    origH = opts.origSize(1);
    origW = opts.origSize(2);
else
    origH = double(max(events(:, 2))) + 1;  % +1 because coords are 0-based
    origW = double(max(events(:, 1))) + 1;
end
newH = ceil(origH / scale);
newW = ceil(origW / scale);

fprintf('  sensor size: %dx%d (HxW)  ->  %dx%d  (scale = %d)\n', ...
    origH, origW, newH, newW, scale);

%% ---- Downscale coordinates --------------------------------------------
% floor(x/scale) keeps the 0-based convention so detectAruco's `+1`
% shift still works on the reduced file.
reducedEvents      = events;
reducedEvents(:,1) = floor(events(:,1) / scale);
reducedEvents(:,2) = floor(events(:,2) / scale);

%% ---- Decide output path -----------------------------------------------
if isempty(opts.outPath)
    [d, n] = fileparts(matFile);
    outPath = fullfile(d, [n '_reduced.mat']);
else
    outPath = opts.outPath;
end

%% ---- Save (preserve any other variables) ------------------------------
allVars.events    = reducedEvents;
allVars.reduction = struct( ...
    'scale',     scale, ...
    'origSize',  [origH origW], ...
    'newSize',   [newH  newW], ...
    'sourceFile', matFile, ...
    'createdAt', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<DATST,TNOW1>
save(outPath, '-struct', 'allVars');
fprintf('  saved -> %s\n', outPath);
fprintf('  REMINDER: set sensorSize = [%d %d] in main.m when running detectAruco on the reduced file.\n', ...
    newH, newW);

%% ---- Visual comparison ------------------------------------------------
if opts.show
    showComparison(events, reducedEvents, ...
        origH, origW, newH, newW, ...
        opts.windowMs, opts.numFrames, matFile);
end
end


%% =========================================================================
function showComparison(eventsA, eventsB, hA, wA, hB, wB, winMs, nFrames, fileA)
    tAll  = eventsA(:, 4);
    tMin  = double(min(tAll));
    tMax  = double(max(tAll));
    winUs = winMs * 1000;

    % Evenly-spaced sample times across the recording. We start at
    % tMin + winUs so each slice has a full window of history.
    if tMax - tMin < winUs
        warning('Recording (%.2f ms) shorter than window (%d ms); using one slice.', ...
            (tMax - tMin) / 1000, winMs);
        tNowList = (tMin + tMax) / 2;
        nFrames  = 1;
    else
        tNowList = linspace(tMin + winUs, tMax, nFrames);
    end

    [~, baseName] = fileparts(fileA);
    fig = figure('Name', sprintf('Reduction preview: %s', baseName), ...
                 'Color', 'w', ...
                 'Position', [80 80 250*nFrames + 80, 620]);
    tl = tiledlayout(fig, 2, nFrames, ...
                     'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:nFrames
        tNow = tNowList(i);

        % top row: original resolution
        nexttile(tl, i);
        imgA = accumulateWindow(eventsA, tNow, winUs, hA, wA);
        imshow(autoStretch(imgA), 'InitialMagnification', 'fit');
        title(sprintf('t = %.3f s', tNow / 1e6), 'FontSize', 10);
        if i == 1
            ylabel(sprintf('original\n%dx%d', hA, wA), ...
                   'FontWeight', 'bold', 'Visible', 'on', ...
                   'Rotation', 0, 'HorizontalAlignment', 'right');
        end

        % bottom row: reduced resolution
        nexttile(tl, nFrames + i);
        imgB = accumulateWindow(eventsB, tNow, winUs, hB, wB);
        imshow(autoStretch(imgB), 'InitialMagnification', 'fit');
        if i == 1
            ylabel(sprintf('reduced\n%dx%d', hB, wB), ...
                   'FontWeight', 'bold', 'Visible', 'on', ...
                   'Rotation', 0, 'HorizontalAlignment', 'right');
        end
    end

    title(tl, sprintf('Side-by-side cumulative events  ·  %d ms window', winMs), ...
          'FontSize', 12, 'FontWeight', 'bold');
end


%% =========================================================================
function img = accumulateWindow(events, tNow, winUs, H, W)
    tCol = events(:, 4);
    mask = (tCol > (tNow - winUs)) & (tCol <= tNow);
    if ~any(mask)
        img = zeros(H, W);
        return;
    end
    x = double(events(mask, 1)) + 1;   % to 1-based for accumarray
    y = double(events(mask, 2)) + 1;
    valid = (x >= 1) & (x <= W) & (y >= 1) & (y <= H);
    if ~any(valid)
        img = zeros(H, W);
        return;
    end
    img = accumarray([y(valid), x(valid)], 1, [H, W]);
end


%% =========================================================================
function out = autoStretch(img)
    img = double(img);
    mx  = max(img(:));
    if mx == 0
        out = uint8(zeros(size(img)));
    else
        out = uint8(img / mx * 255);
    end
end
