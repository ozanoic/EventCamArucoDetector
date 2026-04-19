%% Event Camera ArUco Marker Detection
% Configure input file and parameters, then run detection.

clear; close all; clc;

% Make utility functions (converters, loaders, generators, visualizers) available
addpath('Utils');

%% ---- Input ----
matFile    = 'Data/marker_z2_zoom_high/marker_z2_zoom_high.mat';
sensorSize = [240, 320];   % [height, width]

%% ---- Parameters ----
% params.windowDurations_ms = [5, 10, 15, 20, 30, 50, 70, 100, 125, 150];
% params.windowDurations_ms = [150, 200, 250, 300, 350, 400, 450, 500, 600, 700];
params.windowDurations_ms = [100];
params.tickStep_us        = 1000;       % 1 ms tick step
params.showVis            = 1;
params.useParallel        = true;       % true = use parfor if toolbox is available; false = always sequential

% Marker grid (ARUCO_MIP_36h12: 8x8 grid, 6x6 inner code)
params.numCells = 8;
params.codeSize = 6;
params.cellPx   = 20;      % pixels per cell in unwarped image

% Blob detection
params.blobParams.minArea   = 500;
params.blobParams.maxArea   = sensorSize(1) * sensorSize(2) * 0.4;
params.blobParams.maxAspect = 3.0;

%% ---- Run detection ----
results = detectAruco(matFile, sensorSize, params);

%% ---- Save results ----
[inputDir, inputName, ~] = fileparts(matFile);
outputFile = fullfile(inputDir, [inputName '_results.mat']);
save(outputFile, '-struct', 'results');
fprintf('Results saved to %s\n', outputFile);
