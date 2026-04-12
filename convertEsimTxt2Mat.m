%% Convert ESIM event .txt files to .mat format
%  Reads events.txt (x y pol timestamp_ns) and frames.txt (timestamp_ns filename)
%  Saves events as [x, y, polarity, timestamp_us] with polarity {-1,+1}
%  Also saves frame info (timestamps + filenames) in the same .mat
%
%  Usage:
%    convertEsimTxt2Mat(dataDir)
%    convertEsimTxt2Mat(dataDir, outputFile)
%
%  dataDir should contain events.txt and frames.txt

clear; close all; clc;

%% ---- Settings ----
% Choose which dataset to convert (comment/uncomment)
% dataDir = '../Data/Synthetic/MovingCam/moving_events';
dataDir = '../Data/Synthetic/MovingCam/moving_events_fast';

%% ---- Read events ----
eventsFile = fullfile(dataDir, 'events.txt');
fprintf('Reading events from: %s\n', eventsFile);
fprintf('  (this may take a while for large files...)\n');

% Format: x y pol timestamp_ns (space-separated)
% Use textscan for speed on large files
fid = fopen(eventsFile, 'r');
if fid == -1, error('Cannot open %s', eventsFile); end
raw = textscan(fid, '%d %d %d %d64', 'Delimiter', ' \t', 'MultipleDelimsAsOne', true);
fclose(fid);

x   = double(raw{1});     % 0-indexed pixel x
y   = double(raw{2});     % 0-indexed pixel y
pol = double(raw{3});     % 0 or 1
ts  = double(raw{4});     % already in microseconds

nEvents = length(x);
fprintf('  Loaded %d events\n', nEvents);
fprintf('  Time range: %.3f - %.3f s\n', ts(1)/1e6, ts(end)/1e6);
fprintf('  X range: %d - %d\n', min(x), max(x));
fprintf('  Y range: %d - %d\n', min(y), max(y));

% Polarity stays as-is: 0 (OFF) and 1 (ON)

% Make timestamps relative (first event at 0)
ts = ts - ts(1);

% Build events matrix: [x, y, polarity, timestamp_us]
events = [x, y, pol, ts];

fprintf('  Polarity values: {%d, %d}\n', min(events(:,3)), max(events(:,3)));
fprintf('  Timestamp range: 0 - %.0f us (%.3f s)\n', ts(end), ts(end)/1e6);

%% ---- Read frames info ----
framesFile = fullfile(dataDir, 'frames.txt');
if exist(framesFile, 'file')
    fprintf('\nReading frames from: %s\n', framesFile);
    fid = fopen(framesFile, 'r');
    rawFrames = textscan(fid, '%d64 %s', 'Delimiter', ' \t', 'MultipleDelimsAsOne', true);
    fclose(fid);

    frameTimestamps_us = double(rawFrames{1});
    frameFilenames     = rawFrames{2};

    % Make relative to first event
    frameTimestamps_us = frameTimestamps_us - double(raw{4}(1));

    nFrames = length(frameFilenames);
    fprintf('  Loaded %d frames\n', nFrames);
    fprintf('  Frame time range: %.3f - %.3f s\n', ...
        frameTimestamps_us(1)/1e6, frameTimestamps_us(end)/1e6);
else
    frameTimestamps_us = [];
    frameFilenames = {};
    fprintf('\nNo frames.txt found — saving events only.\n');
end

%% ---- Save .mat ----
[~, dirName] = fileparts(dataDir);
outFile = fullfile(dataDir, [dirName '.mat']);

fprintf('\nSaving to: %s\n', outFile);
save(outFile, 'events', 'frameTimestamps_us', 'frameFilenames', '-v7.3');

info = dir(outFile);
fprintf('Done! File size: %.1f MB\n', info.bytes / 1e6);
fprintf('\nVariables saved:\n');
fprintf('  events             [%d x 4] — [x, y, polarity, timestamp_us]\n', size(events,1));
fprintf('  frameTimestamps_us [%d x 1] — frame capture times in us\n', length(frameTimestamps_us));
fprintf('  frameFilenames     [%d x 1] — frame image filenames\n', length(frameFilenames));
