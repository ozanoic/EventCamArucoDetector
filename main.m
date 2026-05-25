%% Event Camera ArUco Marker Detection
% Configure input file and parameters, then run detection.
clear; close all; clc;

% Make utility functions (converters, loaders, generators, visualizers) available
addpath("Data");
addpath("Utils");

% filterEventNoise('Data/OzanEventData_22.05.2026/1/1.mat', struct('method', 'hot+bg', 'dt_ms', 1));

%% ---- Input ----
matFiles    = ["Data/OzanEventData_22.05.2026/3/3_reduced.mat"];
% reduceEventResolution('Data/OzanEventData_22.05.2026/4/4.mat');

% matFiles    = ["Data/OzanEventData_22.05.2026/1/1_reduced.mat",...
%     "Data/OzanEventData_22.05.2026/2/2_reduced.mat",...
%     "Data/OzanEventData_22.05.2026/3/3_reduced.mat",...
%     "Data/OzanEventData_22.05.2026/4/4_reduced.mat",...
%     "Data/OzanEventData_22.05.2026/5/5_reduced.mat"];

% matFiles    = ["Data/marker_z2_cross_low/marker_z2_cross_low.mat", ...
%     "Data/marker_z2_cross_med/marker_z2_cross_med.mat", ...
%     "Data/marker_z2_cross_high/marker_z2_cross_high.mat", ...
%     "Data/marker_z2_linear_low/marker_z2_linear_low.mat", ...
%     "Data/marker_z2_linear_med/marker_z2_linear_med.mat", ...
%     "Data/marker_z2_linear_high/marker_z2_linear_high.mat", ...
%     "Data/marker_z2_rotation_low/marker_z2_rotation_low.mat", ...
%     "Data/marker_z2_rotation_med/marker_z2_rotation_med.mat", ...
%     "Data/marker_z2_rotation_high/marker_z2_rotation_high.mat", ...
%     "Data/marker_z2_zoom_low/marker_z2_zoom_low.mat", ...
%     "Data/marker_z2_zoom_med/marker_z2_zoom_med.mat", ...
%     "Data/marker_z2_zoom_high/marker_z2_zoom_high.mat"];

sensorSize = [240, 320];   % [height, width]
% sensorSize = [480, 640];   % [height, width]

%% ---- Parameters ----
params.windowDurations_ms = [3, 5, 8, 10, 15, 20, 30, 50, 70];
% params.windowDurations_ms = [150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750];
params.tickStep_us        = 1000000;       % 1 ms tick step
params.showVis            = false;
params.useParallel        = true;       % true = use parfor if toolbox is available; false = always sequential

% Which ArUco IDs are we hunting for?
%   []          -> accept ANY decoded marker
%   3           -> only count detections of marker id 3
%   [3 7 12]    -> accept any of these IDs
params.requestedMarkerIds = [3, 8];

% ---- Speed knobs (recommended for 480x640 sensors) -----------------------
% earlyExitOnFirstHit: skip remaining windows once one decodes at a tick.
%   Big speedup; the per-window breakdown is no longer meaningful when on.
params.earlyExitOnFirstHit = true;

% detectScale: run blob detection at this fraction of full resolution.
%   1.0 = no downscale (default).  0.5 = ~4x faster blob stage.
%   Corners are scaled back up before the perspective warp.
params.detectScale = 1.0;

% useGPU: push the per-quad imwarp to GPU. Requires an NVIDIA GPU and
%   forces sequential execution (parfor can't share a GPU efficiently).
params.useGPU = false;

% Tip: the cheapest speedup is a coarser tick step. Doubling tickStep_us
%   (1000 -> 2000) halves the work with almost no detection-quality loss
%   for typical marker motion.

% ---- Detection-quality knobs (for noisy / real DVS data) -----------------
% decoderHammingDist: 0 = exact dictionary match (default). On real DVS,
%   set this to 1 or 2 to accept decodes that are 1-2 bits off. Big
%   detection-rate improvement on noisy data because a single misread
%   cell no longer kills the entire decode. ARUCO_MIP_36h12 has
%   minimum inter-marker Hamming distance 12, so values up to ~4 stay
%   unambiguous.
params.decoderHammingDist = 4;

% minEventsPerWindow: skip windows with fewer than N events in their
%   [tNow - dt, tNow] range. Was hard-coded to 10. Faint markers can
%   land just under that on 480x640 sensors -- lower it to 5 if the
%   short windows are returning zero detections in your runs.
params.minEventsPerWindow = 5;

% refineCorners: sub-pixel corner refinement. The blob detector's
%   min-area-rect corners are typically 2-4 px off for noisy event
%   blobs, which shifts every cell-boundary sample in the unwarped
%   image and makes the decoder miss the marker. Turn this on if
%   diagnoseFailures shows "marker is outlined in red, but the
%   outline corners are clearly off the real marker corners".
params.refineCorners  = true;
params.refineSearchPx = 5;     % perpendicular search half-width (px)

% Marker grid (ARUCO_MIP_36h12: 8x8 grid, 6x6 inner code)
params.numCells = 8;
params.codeSize = 6;
params.cellPx   = 20;      % pixels per cell in unwarped image

% Blob detection
params.blobParams.minArea   = 100;
params.blobParams.maxArea   = sensorSize(1) * sensorSize(2) * 0.6;
params.blobParams.maxAspect = 3.0;

for i = 1:length(matFiles)
    %% ---- Run detection ----
    results = detectAruco(matFiles(i), sensorSize, params);

    %% ---- Save results ----
    [inputDir, inputName, ~] = fileparts(matFiles(i));
    outputFile = fullfile(inputDir, inputName + "_results_v3.mat");
    save(outputFile, '-struct', 'results');
    fprintf('Results saved to %s\n', outputFile);

    diagnoseFailures(outputFile, matFiles(i), struct('mode', 'both'));
end

%%
% mergeAllResults('Data');

%%
% viewAllResults('Data')