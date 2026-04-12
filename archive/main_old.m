%% Event-by-Event ArUco Border Detection
% Processes events one-by-one, applies noise filtering, maintains a
% Surface of Active Events (SAE), and periodically tries to detect
% square (marker) borders in the active event pattern.
%
% Uses Sarmadi et al. (2021) method: separate ON/OFF polarity images,
% detect line segments in each, match ON-OFF pairs to form quad candidates.
%
% Step 1: Just detect marker borders (squares). No ID decoding yet.

clear; close all; clc;

%% ---- Load events -----------------------------------------------------------
% dataset = '../Data/roll_90dps_10ms_g1/events_out_start_cropped.mat';
% dataset = '../Data/roll_360dps_10ms_g1/events_out_start_cropped.mat';
% dataset = '../Data/zoom_1000_10_10ms_g1_v3/events_out_start_cropped.mat';
dataset = '../Data/zoom_1000_100_10ms_g1/events_out_start_cropped.mat';
sarmadiDataset = '../Data/Sarmadi/side2side/side2side/packets.bin';
events = convertSarmadiBin(sarmadiDataset);
% events  = double(loadEvents(dataset));
numEvents = size(events, 1);

fprintf('Loaded %d events  |  sensor 640x480\n', numEvents);
fprintf('  x: [%d, %d]   y: [%d, %d]\n', ...
    min(events(:,1)), max(events(:,1)), min(events(:,2)), max(events(:,2)));
fprintf('  t: [%.0f, %.0f] us   duration: %.1f ms\n', ...
    min(events(:,4)), max(events(:,4)), ...
    (max(events(:,4))-min(events(:,4)))/1000);

%% ---- Parameters (tune these) -----------------------------------------------
% sensorSize   = [480, 640];       % [height, width]
sensorSize   = [128, 128];       % [height, width]

% Noise filtering
refractoryDt = 100;             % us - suppress same-pixel repeats within this
corrDt       = 100;             % us - spatial correlation window (3x3 neighborhood)

% Visualization & detection
activeWindow = 10000;            % us - time window for "current" active events
vizInterval  = 10000;            % visualize & detect every N events

% Sarmadi detection parameters
sarmadiParams.activeWindow    = activeWindow;
sarmadiParams.minSegLen       = 25;     % minimum line segment length (px)
sarmadiParams.maxSegLen       = 250;    % maximum line segment length (px)
sarmadiParams.angleTol        = 30;     % max angle difference ON↔OFF (degrees)
sarmadiParams.lengthRatio     = 0.5;    % not used directly, paper uses 2x rule
sarmadiParams.numPeaks        = 50;     % Hough peaks to detect
sarmadiParams.minArea         = 25*25;  % min quad area (pixels)
sarmadiParams.maxArea         = 150*150;% max quad area (pixels)
sarmadiParams.maxAspect       = 3.0;    % max side ratio
sarmadiParams.houghFillGap    = 5;      % max gap in Hough line segments
sarmadiParams.houghThreshFrac = 0.15;   % Hough peak threshold as fraction of max

%% ---- Initialize ------------------------------------------------------------
height = sensorSize(1);
width  = sensorSize(2);

sae     = zeros(sensorSize);    % Combined SAE (for visualization)
sae_on  = zeros(sensorSize);    % SAE for ON  events (polarity = +1)
sae_off = zeros(sensorSize);    % SAE for OFF events (polarity = -1)

accFrame = zeros(sensorSize);   % Accumulated filtered event count
nFiltered = 0;
nRefractory = 0;
nNoCorr    = 0;

hFig = figure('Name','Event-by-Event ArUco Border Detection (Sarmadi)', ...
              'Position',[50 50 1500 750]);

%% ---- Main event loop -------------------------------------------------------
fprintf('\nProcessing events...\n');
tic;

for i = 1:numEvents
    x   = events(i,1) + 1;        % 0-indexed -> 1-indexed
    y   = events(i,2) + 1;
    pol = events(i,3);             % polarity: +1 (ON) or -1 (OFF)
    t   = events(i,4);

    % Bounds check
    if x < 1 || x > width || y < 1 || y > height, continue; end

    % ---- (1) Refractory period filter ----
    lastTime = sae(y, x);
    if lastTime > 0 && (t - lastTime) < refractoryDt
        sae(y, x) = t;
        if pol > 0, sae_on(y,x) = t; else, sae_off(y,x) = t; end
        nRefractory = nRefractory + 1;
        continue;
    end

    % ---- (2) Spatial correlation filter ----
    x_lo = max(1, x-1);  x_hi = min(width,  x+1);
    y_lo = max(1, y-1);  y_hi = min(height, y+1);

    nbr = sae(y_lo:y_hi, x_lo:x_hi);
    cy = y - y_lo + 1;
    cx = x - x_lo + 1;
    nbr(cy, cx) = 0;
    isSignal = any(nbr(:) > 0 & (t - nbr(:)) < corrDt);

    % ---- Update SAE (always) ----
    sae(y, x) = t;
    if pol > 0
        sae_on(y, x) = t;
    else
        sae_off(y, x) = t;
    end

    if isSignal
        nFiltered = nFiltered + 1;
        accFrame(y, x) = accFrame(y, x) + 1;
    else
        nNoCorr = nNoCorr + 1;
    end

    % ---- (3) Periodic visualization & square detection ----
    if mod(i, vizInterval) == 0
        % Combined active mask for visualization
        activeMask = (sae > 0) & ((t - sae) <= activeWindow);

        % ---- Sarmadi et al. detection ----
        [quads, dbg] = detectQuadSarmadi(sae_on, sae_off, t, sarmadiParams);

        fprintf('[event %6d]  quads:%d  on_segs:%d  off_segs:%d\n', ...
            i, length(quads), length(dbg.on_segs), length(dbg.off_segs));

        % ---- Visualize (2x3 layout) ----
        figure(hFig); clf;

        % (a) Accumulated filtered events
        subplot(2,3,1);
        imagesc(accFrame); colormap(gca,'hot'); colorbar;
        title(sprintf('Accumulated filtered events\n%d / %d  (%.0f%%)', ...
            i, numEvents, 100*i/numEvents));
        axis image;

        % (b) ON events image (preprocessed)
        subplot(2,3,2);
        onRGB = zeros(height, width, 3);
        onRGB(:,:,1) = dbg.on_img;   % red channel = ON
        imshow(onRGB * 3); hold on;
        for si = 1:length(dbg.on_segs)
            s = dbg.on_segs(si);
            plot([s.p1(1) s.p2(1)], [s.p1(2) s.p2(2)], 'r-', 'LineWidth', 2);
        end
        title(sprintf('ON events + segments (%d)', length(dbg.on_segs)));
        axis image; hold off;

        % (c) OFF events image (preprocessed)
        subplot(2,3,3);
        offRGB = zeros(height, width, 3);
        offRGB(:,:,3) = dbg.off_img;  % blue channel = OFF
        imshow(offRGB * 3); hold on;
        for si = 1:length(dbg.off_segs)
            s = dbg.off_segs(si);
            plot([s.p1(1) s.p2(1)], [s.p1(2) s.p2(2)], 'b-', 'LineWidth', 2);
        end
        title(sprintf('OFF events + segments (%d)', length(dbg.off_segs)));
        axis image; hold off;

        % (d) Combined ON+OFF overlay
        subplot(2,3,4);
        combRGB = zeros(height, width, 3);
        combRGB(:,:,1) = dbg.on_img;     % red = ON
        combRGB(:,:,3) = dbg.off_img;    % blue = OFF
        imshow(combRGB * 3); hold on;
        % Draw all segments
        for si = 1:length(dbg.on_segs)
            s = dbg.on_segs(si);
            plot([s.p1(1) s.p2(1)], [s.p1(2) s.p2(2)], 'r-', 'LineWidth', 1.5);
        end
        for si = 1:length(dbg.off_segs)
            s = dbg.off_segs(si);
            plot([s.p1(1) s.p2(1)], [s.p1(2) s.p2(2)], 'b-', 'LineWidth', 1.5);
        end
        title('ON (red) + OFF (blue) overlay');
        axis image; hold off;

        % (e) Active events
        subplot(2,3,5);
        imshow(activeMask);
        title(sprintf('Active events (window=%d us)\nt = %.0f us', ...
            activeWindow, t));
        axis image;

        % (f) Detected quad borders
        subplot(2,3,6);
        imshow(activeMask); hold on;
        qColors = lines(max(length(quads),1));
        for q = 1:length(quads)
            c = quads{q};
            plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], '-', ...
                 'Color', qColors(q,:), 'LineWidth', 2);
            plot(c(:,1), c(:,2), 'o', 'Color', qColors(q,:), ...
                 'MarkerSize', 8, 'MarkerFaceColor', qColors(q,:));
        end
        title(sprintf('Detected borders: %d', length(quads)));
        axis image; hold off;

        sgtitle(sprintf('Sarmadi Method  |  Event %d / %d  |  filtered: %d  |  quads: %d', ...
            i, numEvents, nFiltered, length(quads)));
        drawnow;

        % Log when quads are found
        if ~isempty(quads)
            fprintf('[event %6d]  %d quad(s) detected!\n', i, length(quads));
            for q = 1:length(quads)
                c = quads{q};
                fprintf('   quad %d  area=%.0f  corners: (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f) (%.0f,%.0f)\n', ...
                    q, polyarea(c(:,1),c(:,2)), ...
                    c(1,1),c(1,2), c(2,1),c(2,2), c(3,1),c(3,2), c(4,1),c(4,2));
            end
        end
    end
end

elapsed = toc;
fprintf('\n=== Done ===\n');
fprintf('Total events     : %d\n', numEvents);
fprintf('Filtered (signal): %d  (%.1f%%)\n', nFiltered, 100*nFiltered/numEvents);
fprintf('Refractory reject: %d  (%.1f%%)\n', nRefractory, 100*nRefractory/numEvents);
fprintf('No-correlation   : %d  (%.1f%%)\n', nNoCorr, 100*nNoCorr/numEvents);
fprintf('Elapsed time     : %.1f s\n', elapsed);
