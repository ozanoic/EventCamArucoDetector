function outPath = reduceEventResolution(matFile, scale, opts)
%REDUCEEVENTRESOLUTION  Downscale event-camera (x, y) and save / visualise.
%
%  outPath = reduceEventResolution(matFile)             % scale = 2, show = true
%  outPath = reduceEventResolution(matFile, scale)
%  outPath = reduceEventResolution(matFile, scale, opts)
%  outPaths = reduceEventResolution(["a.mat","b.mat"], ...)    % batch
%
%  Loads `events` from matFile (Nx4 -- x, y, polarity, timestamp),
%  downsamples the (x, y) coordinates by an integer factor `scale`
%  (default 2), and saves a sibling .mat file named `<name>_reduced.mat`.
%  Other variables in the source file are preserved as-is. The save uses
%  MAT v7.3 so the events array can exceed 2 GB (older MAT versions
%  silently truncate such variables, which writes a near-empty file).
%
%  Optionally shows a side-by-side comparison: original sensor resolution
%  on top, reduced resolution on the bottom, several time slices across
%  the recording, each accumulating a configurable window (default 10 ms).
%
%  Inputs:
%    matFile : path (string or char), OR string/cell array of paths
%              (each is processed and saved to its own sibling).
%    scale   : positive integer downscale factor (default 2)
%    opts    : optional struct with fields
%        .show       (default true)        show the side-by-side figure
%        .windowMs   (default 10)          accumulation window per slice (ms)
%        .numFrames  (default 6)           how many time slices to draw
%        .outPath    (default '')          explicit output path (single file only)
%        .origSize   (default inferred)    [H W] of source sensor
%
%  Output:
%    outPath : path to the saved <name>_reduced.mat
%              (cell array when matFile is an array).
%
%  Examples:
%    reduceEventResolution('Data/OzanEventData_22.05.2026/1/1.mat');
%    reduceEventResolution(["Data/.../1.mat","Data/.../2.mat"], 2, ...
%                          struct('show', false));   % batch, no figure

if nargin < 2 || isempty(scale), scale = 2;           end
if nargin < 3,                   opts  = struct();    end
if ~isfield(opts,'show'),        opts.show      = true; end
if ~isfield(opts,'windowMs'),    opts.windowMs  = 10;   end
if ~isfield(opts,'numFrames'),   opts.numFrames = 6;    end
if ~isfield(opts,'outPath'),     opts.outPath   = '';   end
if ~isfield(opts,'origSize'),    opts.origSize  = [];   end

% -- Batch dispatch: if the caller passes more than one path, loop. --
matFile = string(matFile);
matFile = matFile(:);                  % column vector of strings
if numel(matFile) > 1
    if ~isempty(opts.outPath)
        error('reduceEventResolution: opts.outPath is only valid for a single input.');
    end
    outPath = strings(numel(matFile), 1);
    for k = 1:numel(matFile)
        fprintf('\n=== [%d/%d] %s ===\n', k, numel(matFile), matFile(k));
        outPath(k) = reduceEventResolution(matFile(k), scale, opts);
    end
    outPath = cellstr(outPath);
    return;
end
matFile = char(matFile);               % single path from here on

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

% -- Force MAT v7.3 (HDF5). v7 silently truncates variables bigger than
% 2 GB, which is exactly the "tiny output file" failure mode that
% appears with long 480x640 recordings. v7.3 has no per-variable limit.
save(outPath, '-struct', 'allVars', '-v7.3');

% -- Sanity check the file actually contains the events we just wrote.
info = dir(outPath);
expectedBytes = numel(reducedEvents) * 8;   % rough lower bound (double)
if isempty(info)
    error('reduceEventResolution: save failed -- no file at %s', outPath);
end
if info.bytes < 0.1 * expectedBytes
    warning(['reduceEventResolution: saved file is %s but we wrote ' ...
             '%s of events. Something is wrong with the save.'], ...
             humanBytes(info.bytes), humanBytes(expectedBytes));
end

fprintf('  saved -> %s  (%s on disk)\n', outPath, humanBytes(info.bytes));
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


%% =========================================================================
function s = humanBytes(n)
    units = {'B','KB','MB','GB','TB'};
    i = 1;
    while n >= 1024 && i < length(units)
        n = n / 1024;
        i = i + 1;
    end
    s = sprintf('%.2f %s', n, units{i});
end
