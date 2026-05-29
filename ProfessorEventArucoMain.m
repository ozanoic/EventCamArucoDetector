function ProfessorEventArucoMain()
%PROFESSOREVENTARUCOMAIN Entry point for the real-data-focused detector.
%
% This wrapper intentionally writes new result files only:
%   <input_name>_professor_results.mat
%
% The original detector files are left intact. Edit the matFiles list below
% when you want to run a different recording.

close all; clc;

addpath("Data");
addpath("Utils");

%% ---- Input ---------------------------------------------------------------
matFiles = "Data/OzanEventData_22.05.2026/4/4.mat";
% matFiles = "Data/marker_z2_cross_low/marker_z2_cross_low.mat";
sensorSize = [480, 640];   % [height, width]
% sensorSize = [240, 320];

%% ---- Real-data detector parameters --------------------------------------
params = ProfessorEventArucoParams(sensorSize);

% Markers in the real recordings.
% params.requestedMarkerIds = [3, 8];
params.requestedMarkerIds = 3;

% Start with one file and inspect the console funnel before batch runs.
for i = 1:numel(matFiles)
    fprintf('\n=== Professor detector: %s ===\n', matFiles(i));
    [inputDir, inputName, ~] = fileparts(matFiles(i));
    paramsRun = params;
    paramsRun.checkpointFile = fullfile(inputDir, inputName + params.outputTag + "_partial.mat");

    results = detectArucoProfessor(matFiles(i), sensorSize, paramsRun);

    outputFile = fullfile(inputDir, inputName + params.outputTag + "_results.mat");
    save(outputFile, '-struct', 'results', '-v7.3');
    fprintf('Professor results saved to %s\n', outputFile);

    previousFile = fullfile(inputDir, inputName + "_professor_results.mat");
    mergedFile = fullfile(inputDir, inputName + paramsRun.outputTag + "_merged_results.mat");
    if isfile(previousFile) && ~strcmp(outputFile, previousFile)
        mergeProfessorArucoResults(previousFile, outputFile, mergedFile, paramsRun.speedProfile);
    end
end
end


function params = ProfessorEventArucoParams(sensorSize)
%PROFESSOREVENTARUCOPARAMS Defaults tuned for noisy real DVS recordings.

params.windowDurations_ms = [3 5 8 10 15 20 30 50 70 100 150 250 400 600];

% The old main had 1000000 us, which is 1 second. Real marker visibility
% changes much faster, so this detector samples every 5 ms by default.
params.tickStep_us = 5000;

% Marker grid: ARUCO_MIP_36h12 uses an 8x8 marker with 6x6 data bits.
params.numCells = 8;
params.codeSize = 6;
params.cellPx = 20;

% Candidate generation. The detector fuses count, time-surface, ON, and
% OFF masks, so minArea can be a little more permissive than the old run.
params.blobParams.minArea = 250;
params.blobParams.maxArea = sensorSize(1) * sensorSize(2) * 0.35;
params.blobParams.maxAspect = 3.0;
params.blobParams.minRectangularity = 0.45;
params.blobParams.minPixels = 15;

% Time surface. Small tau favors very recent edges; larger tau keeps faint
% slow marker edges alive. The fixed windows still gate the event history.
params.timeSurfaceTau_ms = 25;
params.timeSurfaceWeight = 0.75;
params.timeSurfaceThreshold = 0.08;

% Decode thresholds. ARUCO_MIP_36h12 has minimum distance 12, so accepting
% a few bit errors is still normally unambiguous.
params.maxHammingDist = 5;
params.secondBestMargin = 2;
params.transitionThresholdScale = 0.75;
params.boundaryHalfWidth = 5;
params.borderPenaltyWeight = 3.0;

% Corner refinement and temporal tracking.
params.refineCorners = true;
params.refineSearchPx = 7;
params.trackMaxAge_ms = 250;
params.trackPenaltyWeight = 2.5;
params.trackSmoothing = 0.65;

% Runtime behavior.
params.minEventsPerWindow = 20;
params.earlyExitOnAllMarkers = true;
params.showProgressEveryPct = 2;
params.saveDebug = false;
params.processingScale = 1.0;
params.maxEventsPerWindow = inf;
params.eventLimitMode = "uniform";
params.timelineStartWindow_ms = [];
params.saveCheckpoints = true;
params.checkpointEveryPct = 2;
params.checkpointFile = "";
params.saveCorners = true;
params.saveWindowStats = true;
params.outputTag = "_professor";

% Runtime profile:
%   "fast"         = quick inspection / synthetic checks
%   "balanced"     = high-quality synthetic/smaller recordings
%   "realbalanced" = practical pass for large real recordings
%   "realquality"  = slower high-quality pass for large real recordings
%   "microrescue"  = focused short-window rescue pass for merging
%   "cascadefull"  = full start-to-end cascade using the selected windows
%   "thorough"     = original Professor-style rescue mode
params.speedProfile = "cascadefull";
params = applyProfessorSpeedProfile(params);
end


function params = applyProfessorSpeedProfile(params)
switch char(lower(string(params.speedProfile)))
    case 'fast'
        params.windowDurations_ms = [50 100];
        params.tickStep_us = 10000;
        params.cellPx = 12;
        params.refineCorners = false;
        params.earlyExitMode = "any";
        params.stopAfterFirstHitInWindow = true;
        params.tryTrackFirst = true;
        params.trackFirstRefineCorners = false;
        params.trackFirstDecodeImageMode = "fast";
        params.trackFirstUsePolarityDecode = false;
        params.candidateMaskMode = "fast";
        params.quadFinderMode = "connected";
        params.maxQuadsPerWindow = 8;
        params.decodeImageMode = "fast";
        params.usePolarityDecode = false;
    case 'balanced'
        params.windowDurations_ms = [50 100 250];
        params.tickStep_us = 5000;
        params.cellPx = 16;
        params.refineCorners = true;
        params.earlyExitMode = "any";
        params.stopAfterFirstHitInWindow = true;
        params.tryTrackFirst = true;
        params.trackFirstRefineCorners = false;
        params.trackFirstDecodeImageMode = "fast";
        params.trackFirstUsePolarityDecode = false;
        params.candidateMaskMode = "balanced";
        params.quadFinderMode = "both";
        params.maxQuadsPerWindow = 12;
        params.decodeImageMode = "balanced";
        params.usePolarityDecode = true;
    case 'realbalanced'
        params.windowDurations_ms = [50 100 250];
        params.tickStep_us = 50000;
        params.cellPx = 16;
        params.processingScale = 1.0;
        params.maxEventsPerWindow = inf;
        params.eventLimitMode = "uniform";
        params.refineCorners = true;
        params.earlyExitMode = "any";
        params.stopAfterFirstHitInWindow = true;
        params.tryTrackFirst = true;
        params.trackFirstRefineCorners = false;
        params.trackFirstDecodeImageMode = "fast";
        params.trackFirstUsePolarityDecode = false;
        params.candidateMaskMode = "balanced";
        params.quadFinderMode = "both";
        params.maxQuadsPerWindow = 12;
        params.decodeImageMode = "balanced";
        params.usePolarityDecode = true;
    case 'realquality'
        params.windowDurations_ms = [50 100 250];
        params.tickStep_us = 20000;
        params.cellPx = 16;
        params.processingScale = 1.0;
        params.maxEventsPerWindow = inf;
        params.eventLimitMode = "uniform";
        params.refineCorners = true;
        params.earlyExitMode = "any";
        params.stopAfterFirstHitInWindow = true;
        params.tryTrackFirst = true;
        params.trackFirstRefineCorners = false;
        params.trackFirstDecodeImageMode = "fast";
        params.trackFirstUsePolarityDecode = false;
        params.candidateMaskMode = "balanced";
        params.quadFinderMode = "both";
        params.maxQuadsPerWindow = 12;
        params.decodeImageMode = "balanced";
        params.usePolarityDecode = true;
    case 'microrescue'
        params.windowDurations_ms = [1 2 4 6];
        params.timelineStartWindow_ms = 600;
        params.tickStep_us = 5000;
        params.cellPx = 20;
        params.outputTag = "_professor_micro";
        params.processingScale = 1.0;
        params.maxEventsPerWindow = inf;
        params.eventLimitMode = "uniform";
        params.refineCorners = true;
        params.earlyExitMode = "all";
        params.stopAfterFirstHitInWindow = false;
        params.tryTrackFirst = false;
        params.trackFirstRefineCorners = true;
        params.trackFirstDecodeImageMode = "full";
        params.trackFirstUsePolarityDecode = true;
        params.candidateMaskMode = "full";
        params.quadFinderMode = "both";
        params.maxQuadsPerWindow = inf;
        params.decodeImageMode = "full";
        params.usePolarityDecode = true;
    case 'cascadefull'
        params.windowDurations_ms = [1 2 3 4 5 6 8 10 15 20 30 50 70 100 150 250];
        params.timelineStartWindow_ms = 600;
        params.tickStep_us = 5000;
        params.cellPx = 20;
        params.outputTag = "_professor_cascade";
        params.processingScale = 1.0;
        params.maxEventsPerWindow = inf;
        params.eventLimitMode = "uniform";
        params.refineCorners = true;
        params.earlyExitMode = "all";
        params.windowOrderMode = "adaptive";
        params.stopAfterFirstHitInWindow = false;
        params.stopAfterAllHitsInWindow = true;
        params.tryTrackFirst = false;
        params.trackFirstRefineCorners = true;
        params.trackFirstDecodeImageMode = "full";
        params.trackFirstUsePolarityDecode = true;
        params.candidateMaskMode = "full";
        params.quadFinderMode = "both";
        params.maxQuadsPerWindow = inf;
        params.decodeImageMode = "full";
        params.usePolarityDecode = true;
    case 'thorough'
        params.earlyExitMode = "all";
        params.stopAfterFirstHitInWindow = false;
        params.stopAfterAllHitsInWindow = true;
        params.tryTrackFirst = false;
        params.trackFirstRefineCorners = true;
        params.trackFirstDecodeImageMode = "full";
        params.trackFirstUsePolarityDecode = true;
        params.candidateMaskMode = "full";
        params.quadFinderMode = "both";
        params.maxQuadsPerWindow = inf;
        params.decodeImageMode = "full";
        params.usePolarityDecode = true;
    otherwise
        error('Unknown Professor speedProfile: %s', char(params.speedProfile));
end
end
