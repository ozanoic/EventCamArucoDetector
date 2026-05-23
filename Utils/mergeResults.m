function merged = mergeResults(resultFiles, outputFile)
%MERGERESULTS Merge multiple result files produced from the same input.
%
%  merged = mergeResults(resultFiles)
%  merged = mergeResults(resultFiles, outputFile)
%
%  Inputs:
%    resultFiles - cell array of paths to result .mat files
%                  (all must come from the same input, but can have
%                  different timestamp ranges -- e.g. different window sets
%                  produce different tStart, so tick counts may differ)
%    outputFile  - optional path to save the merged result
%
%  Output:
%    merged - combined struct with union of all windows, aligned to the
%             union of all timestamps, with recomputed anyDetected and
%             detectionsPerWindow.
%
%  Behavior:
%    - Takes the UNION of timestamps across all files
%      (each window is aligned onto this master timeline; timestamps where
%       a window was not computed are marked -1)
%    - Validates that overlapping timestamps are identical across files
%    - Takes the UNION of window durations across all files
%    - If a window appears in multiple files, keeps the first occurrence
%      (with a warning showing which file was kept)
%    - Sorts windows by duration in the output
%    - Recomputes anyDetected = any(win_* >= 0, 2) across all merged windows
%    - Recomputes detectionsPerWindow from the merged win_* columns

if nargin < 2, outputFile = ''; end
if ischar(resultFiles) || isstring(resultFiles)
    resultFiles = {char(resultFiles)};
end
if length(resultFiles) < 2
    warning('mergeResults: at least 2 files recommended (got %d).', length(resultFiles));
end

fprintf('\n=== Merging %d result files ===\n', length(resultFiles));

%% ---- Load all files ----
allData = cell(length(resultFiles), 1);
for k = 1:length(resultFiles)
    fprintf('  [%d] Loading %s\n', k, resultFiles{k});
    allData{k} = load(resultFiles{k});
end

%% ---- Build union timeline ----
allT = [];
for k = 1:length(allData)
    allT = [allT; double(allData{k}.tNow_us(:))]; %#ok<AGROW>
end
mergedT = unique(allT);           % sorted ascending, deduplicated
refN = length(mergedT);

fprintf('  Union timeline: %d ticks\n', refN);
for k = 1:length(allData)
    nk = length(allData{k}.tNow_us);
    fprintf('    file [%d]: %d ticks  (%.1f%% of union)\n', k, nk, 100*nk/refN);
end

%% ---- Collect window columns from all files ----
% Map window_ms -> entry with source file info and data aligned to its file's timeline
winMap = containers.Map('KeyType', 'double', 'ValueType', 'any');

for k = 1:length(allData)
    d = allData{k};
    if ~isfield(d, 'windowDurations_ms')
        warning('File %d has no windowDurations_ms field, skipping.', k);
        continue;
    end
    wins = d.windowDurations_ms;
    for wi = 1:length(wins)
        winKey = double(wins(wi));
        colName = sprintf('win_%dms', winKey);
        if ~isfield(d, colName)
            warning('File %d: missing field %s', k, colName);
            continue;
        end
        if isKey(winMap, winKey)
            existing = winMap(winKey);
            fprintf('  Duplicate window %dms: keeping file [%d], ignoring file [%d]\n', ...
                winKey, existing.fileIdx, k);
        else
            entry.fileIdx   = k;
            entry.data      = d.(colName);
            entry.sourceT   = double(d.tNow_us);
            entry.source    = resultFiles{k};
            winMap(winKey) = entry;
        end
    end
end

%% ---- Sort windows by duration ----
allWins = sort(cell2mat(keys(winMap)));
numWindows = length(allWins);
fprintf('  Merged windows (%d): %s ms\n', numWindows, mat2str(allWins));

%% ---- Build merged struct on the union timeline ----
merged.tNow_us = mergedT;

winMatrix = -1 * ones(refN, numWindows);
for wi = 1:numWindows
    w = allWins(wi);
    entry = winMap(w);
    col = entry.data(:);

    % Map source timestamps to indices in the union timeline
    [isIn, idx] = ismember(entry.sourceT, mergedT);
    if any(~isIn)
        error('mergeResults: win_%dms has timestamps not in union (internal bug).', w);
    end

    if length(col) ~= length(entry.sourceT)
        error('mergeResults: win_%dms data length (%d) != tNow_us length (%d) in its file.', ...
            w, length(col), length(entry.sourceT));
    end

    % Place values at the right rows; rows with no data stay -1
    winMatrix(idx, wi) = col;
    merged.(sprintf('win_%dms', w)) = winMatrix(:, wi);
end

% Recompute overall detection
merged.anyDetected = double(any(winMatrix >= 0, 2));

% Recompute per-window detection counts (use >= 0, which excludes -1 = not-computed)
merged.windowDurations_ms = allWins;
merged.detectionsPerWindow = sum(winMatrix >= 0, 1);

% Also record how many ticks each window was actually attempted on
% (useful since windows may not cover the full union timeline)
attempted = zeros(1, numWindows);
for wi = 1:numWindows
    entry = winMap(allWins(wi));
    attempted(wi) = length(entry.sourceT);
end
merged.attemptedPerWindow = attempted;

% Carry the requested-marker-IDs through (union across files, if present)
reqSet = [];
for k = 1:length(allData)
    if isfield(allData{k}, 'requestedMarkerIds')
        reqSet = union(reqSet, double(allData{k}.requestedMarkerIds(:)'));
    end
end
merged.requestedMarkerIds = reqSet;

% Recompute per-marker fields from the merged win_*ms columns.
% Report IDs = requested set if any, otherwise every ID that actually
% appears in the merged data.
if ~isempty(reqSet)
    reportIds = reqSet;
else
    vals = winMatrix(winMatrix >= 0);
    reportIds = unique(vals)';
end
merged.markerIdsReported = reportIds;
for ri = 1:length(reportIds)
    mid = reportIds(ri);
    hits = winMatrix == mid;
    merged.(sprintf('anyDetected_id%d', mid))         = double(any(hits, 2));
    merged.(sprintf('detectionsPerWindow_id%d', mid)) = sum(hits, 1);
end

%% ---- Summary ----
nAny = sum(merged.anyDetected);
fprintf('\n--- Merged Summary ---\n');
fprintf('Ticks:                 %d\n', refN);
fprintf('Ticks w/ any detect:   %d  (%.1f%%)\n', nAny, 100*nAny/refN);
fprintf('%-10s  %10s  %10s  %10s  %10s\n', 'Window', 'Detect', 'Attempt', 'Rate/attempt', 'Rate/total');
fprintf('%-10s  %10s  %10s  %10s  %10s\n', '------', '------', '-------', '------------', '----------');
for wi = 1:numWindows
    nd = merged.detectionsPerWindow(wi);
    na = merged.attemptedPerWindow(wi);
    fprintf('%-10s  %10d  %10d  %11.1f%%  %9.1f%%\n', ...
        sprintf('%dms', allWins(wi)), ...
        nd, na, ...
        100 * nd / max(na, 1), ...
        100 * nd / refN);
end

%% ---- Save if requested ----
if ~isempty(outputFile)
    save(outputFile, '-struct', 'merged');
    fprintf('\nMerged results saved to %s\n', outputFile);
end

end
