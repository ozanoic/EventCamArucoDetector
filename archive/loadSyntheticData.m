%% Load parsed ESIM synthetic data (events + frames)
%  Loads events.txt, frames.txt, and frame images from parse_bag.py output.
%
%  Usage:
%    [events, frames, frameTimes] = loadSyntheticData('../Data/Synthetic/esim_output');
%
%  Returns:
%    events     — [N x 4] double: [x, y, polarity, timestamp_us]
%                  polarity: -1 (OFF) or +1 (ON)
%    frames     — {M x 1} cell array of uint8 grayscale images
%    frameTimes — [M x 1] double: frame timestamps in microseconds

function [events, frames, frameTimes] = loadSyntheticData(dataDir)

    %% ---- Load events ----
    evFile = fullfile(dataDir, 'events.txt');
    fprintf('Loading events from %s ...\n', evFile);
    raw = dlmread(evFile);
    % raw columns: x, y, polarity(0/1), timestamp_us
    events = raw;
    % Convert polarity: 0 → -1, 1 → +1
    events(:,3) = events(:,3) * 2 - 1;
    fprintf('  %d events loaded. Time range: %.3f - %.3f s\n', ...
        size(events,1), events(1,4)/1e6, events(end,4)/1e6);

    %% ---- Load frames ----
    frFile = fullfile(dataDir, 'frames.txt');
    if ~exist(frFile, 'file')
        fprintf('  No frames.txt found, skipping frames.\n');
        frames = {};
        frameTimes = [];
        return;
    end

    fprintf('Loading frames from %s ...\n', frFile);
    fid = fopen(frFile, 'r');
    lines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line)
            lines{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);

    nFrames = length(lines);
    frameTimes = zeros(nFrames, 1);
    frames = cell(nFrames, 1);

    framesDir = fullfile(dataDir, 'frames');
    for i = 1:nFrames
        parts = strsplit(lines{i});
        frameTimes(i) = str2double(parts{1});
        imgPath = fullfile(framesDir, parts{2});
        frames{i} = imread(imgPath);
    end

    fprintf('  %d frames loaded. Time range: %.3f - %.3f s\n', ...
        nFrames, frameTimes(1)/1e6, frameTimes(end)/1e6);
end
