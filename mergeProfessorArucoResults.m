function merged = mergeProfessorArucoResults(previousFile, newFile, outputFile, sourceTag)
%MERGEPROFESSORARUCORESULTS Merge two Professor detector result files.
%
% The previous result is kept as the base. The new result contributes extra
% any-detection hits and any new window fields. Duplicate window fields from
% the new file are preserved with a source suffix, for example
% win_600ms_micro, so old per-window statistics are not silently changed.

if nargin < 4 || strlength(string(sourceTag)) == 0
    sourceTag = "extra";
end

sourceTag = matlab.lang.makeValidName(char(sourceTag));

previous = load(previousFile);
extra = load(newFile);

if ~isfield(previous, 'tNow_us') || ~isfield(extra, 'tNow_us')
    error('mergeProfessorArucoResults:MissingTime', ...
        'Both result files must contain tNow_us.');
end

if numel(previous.tNow_us) ~= numel(extra.tNow_us) || ...
        any(previous.tNow_us(:) ~= extra.tNow_us(:))
    error('mergeProfessorArucoResults:TimeMismatch', ...
        ['The result time grids do not match. Use the same tickStep_us and ' ...
         'largest window duration before merging.']);
end

merged = previous;
merged.algorithmName = 'mergeProfessorArucoResults';
merged.mergeSources = {char(previousFile), char(newFile)};
merged.mergeSourceTag = sourceTag;

prevAny = double(previous.anyDetected(:) > 0);
extraAny = double(extra.anyDetected(:) > 0);
merged.(sprintf('anyDetected_%s', sourceTag)) = extraAny;
merged.anyDetected_previous = prevAny;
merged.anyDetected = double(prevAny | extraAny);

merged.previousDetectionCount = sum(prevAny);
merged.extraDetectionCount = sum(extraAny);
merged.mergedDetectionCount = sum(merged.anyDetected);
merged.mergedDetectionGain = merged.mergedDetectionCount - merged.previousDetectionCount;

merged = copyExtraWindowFields(merged, extra, sourceTag);
merged = mergePerMarkerAnyFields(merged, previous, extra, sourceTag);
merged = refreshWindowStats(merged);

merged.previousParamsUsed = getFieldOrEmpty(previous, 'paramsUsed');
merged.extraParamsUsed = getFieldOrEmpty(extra, 'paramsUsed');

save(outputFile, '-struct', 'merged', '-v7.3');

fprintf('\nMerged Professor results saved to %s\n', outputFile);
fprintf('Previous detections: %d / %d (%.1f%%)\n', ...
    merged.previousDetectionCount, numel(merged.anyDetected), ...
    100 * merged.previousDetectionCount / max(numel(merged.anyDetected), 1));
fprintf('New-run detections:  %d / %d (%.1f%%)\n', ...
    merged.extraDetectionCount, numel(merged.anyDetected), ...
    100 * merged.extraDetectionCount / max(numel(merged.anyDetected), 1));
fprintf('Merged detections:   %d / %d (%.1f%%), gain %+d ticks\n\n', ...
    merged.mergedDetectionCount, numel(merged.anyDetected), ...
    100 * merged.mergedDetectionCount / max(numel(merged.anyDetected), 1), ...
    merged.mergedDetectionGain);
end


function merged = copyExtraWindowFields(merged, extra, sourceTag)
names = fieldnames(extra);
for i = 1:numel(names)
    name = names{i};
    if isempty(regexp(name, '^win_\d+ms(_id\d+)?$', 'once'))
        continue;
    end

    if isfield(merged, name)
        merged.(sprintf('%s_%s', name, sourceTag)) = extra.(name);
    else
        merged.(name) = extra.(name);
    end
end

if isfield(merged, 'windowDurations_ms') && isfield(extra, 'windowDurations_ms')
    merged.previousWindowDurations_ms = merged.windowDurations_ms;
    merged.extraWindowDurations_ms = extra.windowDurations_ms;
    merged.windowDurations_ms = unique([merged.windowDurations_ms(:); extra.windowDurations_ms(:)])';
end
end


function merged = mergePerMarkerAnyFields(merged, previous, extra, sourceTag)
if ~isfield(previous, 'markerIdsReported') || ~isfield(extra, 'markerIdsReported')
    return;
end

ids = intersect(double(previous.markerIdsReported(:)'), double(extra.markerIdsReported(:)'));
for id = ids
    field = sprintf('anyDetected_id%d', id);
    if ~isfield(previous, field) || ~isfield(extra, field)
        continue;
    end

    prevAny = double(previous.(field)(:) > 0);
    extraAny = double(extra.(field)(:) > 0);
    merged.(sprintf('%s_previous', field)) = prevAny;
    merged.(sprintf('%s_%s', field, sourceTag)) = extraAny;
    merged.(field) = double(prevAny | extraAny);
end
end


function merged = refreshWindowStats(merged)
if ~isfield(merged, 'windowDurations_ms')
    return;
end

durations = double(merged.windowDurations_ms(:)');
counts = zeros(1, numel(durations));
attempts = numel(merged.anyDetected) * ones(1, numel(durations));
for i = 1:numel(durations)
    field = sprintf('win_%dms', durations(i));
    if isfield(merged, field)
        counts(i) = sum(merged.(field)(:) >= 0);
    else
        counts(i) = NaN;
        attempts(i) = NaN;
    end
end

merged.detectionsPerWindow = counts;
merged.attemptedPerWindow = attempts;
end


function value = getFieldOrEmpty(s, field)
if isfield(s, field)
    value = s.(field);
else
    value = [];
end
end
