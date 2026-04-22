%% Event Camera ArUco Marker Detection
% Configure input file and parameters, then run detection.
clear; close all; clc;

% Make utility functions (converters, loaders, generators, visualizers) available
addpath("Data");
addpath("Utils");

%% ---- Input ----
matFiles    = ["Data/marker_z2_cross_low/marker_z2_cross_low.mat", ...
    "Data/marker_z2_cross_med/marker_z2_cross_med.mat", ...
    "Data/marker_z2_linear_high/marker_z2_linear_high.mat", ...
    "Data/marker_z2_linear_low/marker_z2_linear_low.mat", ...
    "Data/marker_z2_linear_med/marker_z2_linear_med.mat", ...
    "Data/marker_z2_rotation_high/marker_z2_rotation_high.mat", ...
    "Data/marker_z2_rotation_low/marker_z2_rotation_low.mat", ...
    "Data/marker_z2_rotation_med/marker_z2_rotation_med.mat", ...
    "Data/marker_z2_zoom_high/marker_z2_zoom_high.mat", ...
    "Data/marker_z2_zoom_low/marker_z2_zoom_low.mat", ...
    "Data/marker_z2_zoom_med/marker_z2_zoom_med.mat"];

sensorSize = [240, 320];   % [height, width]

%% ---- Parameters ----
% params.windowDurations_ms = [5, 10, 15, 20, 30, 50, 70, 100, 125, 150];
params.windowDurations_ms = [150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750];
params.tickStep_us        = 1000;       % 1 ms tick step
params.showVis            = false;
params.useParallel        = true;       % true = use parfor if toolbox is available; false = always sequential

% Marker grid (ARUCO_MIP_36h12: 8x8 grid, 6x6 inner code)
params.numCells = 8;
params.codeSize = 6;
params.cellPx   = 20;      % pixels per cell in unwarped image

% Blob detection
params.blobParams.minArea   = 500;
params.blobParams.maxArea   = sensorSize(1) * sensorSize(2) * 0.4;
params.blobParams.maxAspect = 3.0;

for i = 1:length(matFiles)
    %% ---- Run detection ----
    results = detectAruco(matFiles(i), sensorSize, params);

    %% ---- Save results ----
    [inputDir, inputName, ~] = fileparts(matFiles(i));
    outputFile = fullfile(inputDir, inputName + "_results_v2.mat");
    save(outputFile, '-struct', 'results');
    fprintf('Results saved to %s\n', outputFile);
end

%%
mergeAllResults('Data');

%%
viewAllResults           % defaults to 'Data'
viewAllResults('Data')