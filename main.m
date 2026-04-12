%% Event Camera ArUco Marker Detection
% Configure input file and parameters, then run detection.

clear; close all; clc;

%% ---- Input ----
matFile    = '../Data/Synthetic/MovingCam/moving_events_fast/moving_events_fast.mat';
sensorSize = [240, 320];   % [height, width]

%% ---- Parameters ----
params.windowDurations_ms = [5, 10, 15, 20, 30, 50, 70, 100, 125, 150];
params.tickStep_us        = 1000;       % 1 ms tick step
params.showVis            = false;

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
save('resultTable.mat', '-struct', 'results');
fprintf('Results saved to resultTable.mat\n');
