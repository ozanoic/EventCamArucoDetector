function merged = mergeResults(resultFiles, outputFile)
%MERGERESULTS Merge multiple result files produced from the same input.
%
%  merged = mergeResults(resultFiles)
%  merged = mergeResults(resultFiles, outputFile)
%
%  Inputs:
%    resultFiles - cell array of paths to result .mat files
%                  (all must come from the same input, i.e. same timestamps)
%    outputFile  - optional path to save the merged result
%
%  Output:
%    merged - combined struct with union of all windows, recomputed
%             anyDetected and detectionsPerWindow
%
%  Behavior:
%    - Validates that all files share the same tNow_us vector
%    - Takes the union of window durations across all files
%    - If a window appears in multiple files, keeps the first occurrence
%      (with a warning showing which file was kept)
%    - Sorts windows by duration in the output
%    - Recomputes anyDetected = any(win_* >= 0) across all merged windows
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

%% ---- Validate matching timestamps ----
refT = allData{1}.tNow_us;
refN = length(refT);
for k = 2:length(allData)
    d = allData{k};
    if length(d.tNow_us) ~= refN
        error('mergeResults: file %d has %d ticks, file 1 has %d.', ...
            k, length(d.tNow_us), refN);
    end
    if any(abs(double(d.tNow_us) - double(refT)) > 0.5)
        error('mergeResults: file %d has different timestamps than file 1.', k);
    end
end
fprintf('  Timestamps match: %d ticks\n', refN);

%% ---- Collect window columns from all files ----
% Map window_ms -> (file index, column data, source file path)
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
            entry.fileIdx = k;
            entry.data    = d.(colName);
            entry.source  = resultFiles{k};
            winMap(winKey) = entry;
        end
    end
end

%% ---- Sort windows by duration ----
allWins = sort(cell2mat(keys(winMap)));
numWindows = length(allWins);
fprintf('  Merged windows (%d): %s ms\n', numWindows, mat2str(allWins));

%% ---- Build merged struct ----
merged.tNow_us = refT;

% Fill in win_Xms fields in sorted order
winMatrix = -1 * ones(refN, numWindows);
for wi = 1:numWindows
    w = allWins(wi);
    entry = winMap(w);
    col = entry.data(:);
    if length(col) ~= refN
        error('mergeResults: win_%dms has %d rows, expected %d.', w, length(col), refN);
    end
    merged.(sprintf('win_%dms', w)) = col;
    winMatrix(:, wi) = col;
end

% Recompute overall detection
merged.anyDetected = double(any(winMatrix >= 0, 2));

% Recompute per-window detection counts
merged.windowDurations_ms = allWins;
merged.detectionsPerWindow = sum(winMatrix >= 0, 1);

%% ---- Summary ----
nAny = sum(merged.anyDetected);
fprintf('\n--- Merged Summary ---\n');
fprintf('Ticks:                 %d\n', refN);
fprintf('Ticks w/ any detect:   %d  (%.1f%%)\n', nAny, 100*nAny/refN);
fprintf('%-10s  %10s  %10s\n', 'Window', 'Detections', 'Rate');
fprintf('%-10s  %10s  %10s\n', '------', '----------', '----');
for wi = 1:numWindows
    fprintf('%-10s  %10d  %9.1f%%\n', ...
        sprintf('%dms', allWins(wi)), ...
        merged.detectionsPerWindow(wi), ...
        100 * merged.detectionsPerWindow(wi) / refN);
end

%% ---- Save if requested ----
if ~isempty(outputFile)
    save(outputFile, '-struct', 'merged');
    fprintf('\nMerged results saved to %s\n', outputFile);
end

end
