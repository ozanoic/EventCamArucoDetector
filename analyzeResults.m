function analyzeResults(resultFile)
%ANALYZERESULTS Detailed analysis and visualization of detection results.
%
%  analyzeResults('Data/moving_events_fast/moving_events_fast_results.mat')
%
%  Produces:
%    - Console summary (overall + per-window stats)
%    - Figure 1: Detection timeline (per-window raster + combined)
%    - Figure 2: Per-window detection rate bar chart
%    - Figure 3: Detection heatmap (time x window)
%    - Figure 4: Gap analysis (consecutive misses)
%    - Figure 5: Rolling detection rate over time
%    - Figure 6: Window agreement analysis

%% ---- Load ----
r = load(resultFile);
winMs = r.windowDurations_ms;
numWindows = length(winMs);
numTicks = length(r.tNow_us);

tSec = (r.tNow_us - r.tNow_us(1)) / 1e6;  % relative time in seconds
duration_s = tSec(end) - tSec(1);

% Build detection matrix: numTicks x numWindows (logical)
detMat = false(numTicks, numWindows);
idMat  = -1 * ones(numTicks, numWindows);
for wi = 1:numWindows
    col = r.(sprintf('win_%dms', winMs(wi)));
    idMat(:, wi) = col;
    detMat(:, wi) = col >= 0;
end

anyDet = r.anyDetected > 0;

%% ========================================================================
%  CONSOLE SUMMARY
%% ========================================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('  RESULT ANALYSIS: %s\n', resultFile);
fprintf('================================================================\n');
fprintf('Total ticks:            %d\n', numTicks);
fprintf('Duration:               %.3f s\n', duration_s);
fprintf('Tick step:              %.0f us\n', r.tNow_us(2) - r.tNow_us(1));
fprintf('Windows:                %s ms\n', mat2str(winMs));

% --- Overall ---
nAny = sum(anyDet);
fprintf('\n--- Overall ---\n');
fprintf('Ticks with any detect:  %d / %d  (%.1f%%)\n', nAny, numTicks, 100*nAny/numTicks);
fprintf('Ticks with NO detect:   %d / %d  (%.1f%%)\n', numTicks-nAny, numTicks, 100*(numTicks-nAny)/numTicks);

% --- Per-window ---
fprintf('\n--- Per-Window Detection Rates ---\n');
fprintf('%-10s  %8s  %8s  %10s  %10s  %10s\n', ...
    'Window', 'Detect', 'Miss', 'Rate', 'Unique', 'UniqueRate');
fprintf('%-10s  %8s  %8s  %10s  %10s  %10s\n', ...
    '------', '------', '----', '----', '------', '----------');
for wi = 1:numWindows
    nDet = sum(detMat(:, wi));
    % Unique: ticks where ONLY this window detected
    onlyThis = detMat(:, wi) & sum(detMat, 2) == 1;
    nUnique = sum(onlyThis);
    fprintf('%-10s  %8d  %8d  %9.1f%%  %10d  %9.1f%%\n', ...
        sprintf('%dms', winMs(wi)), nDet, numTicks - nDet, ...
        100*nDet/numTicks, nUnique, 100*nUnique/numTicks);
end

% --- All-windows agreement ---
allAgree = sum(all(detMat, 2));
noneAgree = sum(all(~detMat, 2));
fprintf('\n--- Agreement ---\n');
fprintf('All windows detect:     %d  (%.1f%%)\n', allAgree, 100*allAgree/numTicks);
fprintf('No window detects:      %d  (%.1f%%)\n', noneAgree, 100*noneAgree/numTicks);
fprintf('Partial agreement:      %d  (%.1f%%)\n', ...
    numTicks - allAgree - noneAgree, 100*(numTicks - allAgree - noneAgree)/numTicks);

% --- How many windows detect per tick ---
numWinPerTick = sum(detMat, 2);
fprintf('\n--- Windows Detecting Per Tick ---\n');
fprintf('%-12s  %8s  %8s\n', '#Windows', 'Ticks', 'Pct');
fprintf('%-12s  %8s  %8s\n', '--------', '-----', '---');
for n = 0:numWindows
    cnt = sum(numWinPerTick == n);
    if cnt > 0
        fprintf('%-12d  %8d  %7.1f%%\n', n, cnt, 100*cnt/numTicks);
    end
end

% --- Detected marker IDs ---
allIDs = idMat(idMat >= 0);
if ~isempty(allIDs)
    uniqueIDs = unique(allIDs);
    fprintf('\n--- Detected Marker IDs ---\n');
    fprintf('%-10s  %8s  %8s\n', 'ID', 'Count', 'Pct');
    fprintf('%-10s  %8s  %8s\n', '--', '-----', '---');
    for k = 1:length(uniqueIDs)
        cnt = sum(allIDs == uniqueIDs(k));
        fprintf('%-10d  %8d  %7.1f%%\n', uniqueIDs(k), cnt, 100*cnt/length(allIDs));
    end
end

% --- Gap analysis ---
missRuns = getMissRuns(anyDet);
if ~isempty(missRuns)
    fprintf('\n--- Miss Gap Analysis (consecutive ticks with no detection) ---\n');
    fprintf('Total gap segments:     %d\n', length(missRuns));
    fprintf('Shortest gap:           %d ticks  (%.1f ms)\n', min(missRuns), min(missRuns));
    fprintf('Longest gap:            %d ticks  (%.1f ms)\n', max(missRuns), max(missRuns));
    fprintf('Mean gap:               %.1f ticks  (%.1f ms)\n', mean(missRuns), mean(missRuns));
    fprintf('Median gap:             %.1f ticks  (%.1f ms)\n', median(missRuns), median(missRuns));
end

% --- Detection streak analysis ---
hitRuns = getMissRuns(~anyDet);  % invert: runs of detections
if ~isempty(hitRuns)
    fprintf('\n--- Detection Streak Analysis (consecutive detections) ---\n');
    fprintf('Total streak segments:  %d\n', length(hitRuns));
    fprintf('Shortest streak:        %d ticks\n', min(hitRuns));
    fprintf('Longest streak:         %d ticks  (%.1f ms)\n', max(hitRuns), max(hitRuns));
    fprintf('Mean streak:            %.1f ticks\n', mean(hitRuns));
end

% --- Time-region analysis (split into 10 segments) ---
nSegments = 10;
segLen = floor(numTicks / nSegments);
fprintf('\n--- Detection Rate Over Time (10 segments) ---\n');
fprintf('%-15s  %-15s  %8s  %8s\n', 'Time Start', 'Time End', 'Detect', 'Rate');
fprintf('%-15s  %-15s  %8s  %8s\n', '----------', '--------', '------', '----');
for s = 1:nSegments
    i1 = (s-1)*segLen + 1;
    i2 = min(s*segLen, numTicks);
    if s == nSegments, i2 = numTicks; end
    segDet = sum(anyDet(i1:i2));
    segN = i2 - i1 + 1;
    fprintf('%-15s  %-15s  %8d  %7.1f%%\n', ...
        sprintf('%.3fs', tSec(i1)), sprintf('%.3fs', tSec(i2)), ...
        segDet, 100*segDet/segN);
end

fprintf('================================================================\n\n');

%% ========================================================================
%  FIGURE 1: Detection Timeline (raster plot)
%% ========================================================================
figure('Name', 'Detection Timeline', 'Position', [50 50 1400 600]);

% Top: per-window raster
subplot(4,1,1:3);
hold on;
for wi = 1:numWindows
    tDet = tSec(detMat(:, wi));
    if ~isempty(tDet)
        plot(tDet, wi * ones(size(tDet)), '|', 'Color', getWinColor(wi, numWindows), ...
            'MarkerSize', 4);
    end
end
set(gca, 'YTick', 1:numWindows, 'YTickLabel', ...
    arrayfun(@(x) sprintf('%dms', x), winMs, 'UniformOutput', false));
ylim([0.5 numWindows+0.5]);
xlim([tSec(1) tSec(end)]);
ylabel('Window');
title('Detection Timeline (per window)');
grid on; hold off;

% Bottom: combined any-detection
subplot(4,1,4);
tHit = tSec(anyDet);
tMiss = tSec(~anyDet);
hold on;
bar(tSec, double(anyDet), 1, 'FaceColor', [0.2 0.7 0.2], 'EdgeColor', 'none');
xlim([tSec(1) tSec(end)]);
ylim([0 1.5]);
xlabel('Time (s)');
ylabel('Any');
title(sprintf('Combined: %d/%d (%.1f%%)', nAny, numTicks, 100*nAny/numTicks));
hold off;

%% ========================================================================
%  FIGURE 2: Per-Window Detection Rate
%% ========================================================================
figure('Name', 'Per-Window Detection Rate', 'Position', [50 700 700 400]);

rates = zeros(1, numWindows);
for wi = 1:numWindows
    rates(wi) = 100 * sum(detMat(:, wi)) / numTicks;
end

b = bar(1:numWindows, rates, 'FaceColor', 'flat');
for wi = 1:numWindows
    b.CData(wi,:) = getWinColor(wi, numWindows);
end
set(gca, 'XTick', 1:numWindows, 'XTickLabel', ...
    arrayfun(@(x) sprintf('%dms', x), winMs, 'UniformOutput', false));
ylabel('Detection Rate (%)');
xlabel('Window Duration');
title('Detection Rate Per Window');
ylim([0 105]);
grid on;

% Add percentage labels on bars
for wi = 1:numWindows
    text(wi, rates(wi) + 1.5, sprintf('%.1f%%', rates(wi)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

%% ========================================================================
%  FIGURE 3: Detection Heatmap
%% ========================================================================
figure('Name', 'Detection Heatmap', 'Position', [800 700 700 400]);

% Downsample for readability if too many ticks
maxBins = 500;
if numTicks > maxBins
    binSize = ceil(numTicks / maxBins);
    nBins = ceil(numTicks / binSize);
    heatData = zeros(nBins, numWindows);
    heatTime = zeros(nBins, 1);
    for b2 = 1:nBins
        i1 = (b2-1)*binSize + 1;
        i2 = min(b2*binSize, numTicks);
        heatData(b2, :) = mean(detMat(i1:i2, :), 1);
        heatTime(b2) = tSec(i1);
    end
else
    heatData = double(detMat);
    heatTime = tSec;
end

imagesc(heatTime, 1:numWindows, heatData');
colormap(gca, [0.15 0.15 0.15; parula(255)]);
colorbar;
set(gca, 'YTick', 1:numWindows, 'YTickLabel', ...
    arrayfun(@(x) sprintf('%dms', x), winMs, 'UniformOutput', false));
xlabel('Time (s)');
ylabel('Window Duration');
title('Detection Heatmap (yellow = detected)');
axis xy;

%% ========================================================================
%  FIGURE 4: Gap Distribution
%% ========================================================================
figure('Name', 'Gap Analysis', 'Position', [50 50 1200 400]);

% Per-window gap histograms
subplot(1,2,1);
hold on;
maxGap = 0;
for wi = 1:numWindows
    gaps = getMissRuns(detMat(:, wi));
    if ~isempty(gaps)
        maxGap = max(maxGap, max(gaps));
    end
end
if maxGap > 0
    edges = 0:max(1, floor(maxGap/50)):maxGap+1;
    for wi = numWindows:-1:1
        gaps = getMissRuns(detMat(:, wi));
        if ~isempty(gaps)
            histogram(gaps, edges, 'FaceColor', getWinColor(wi, numWindows), ...
                'FaceAlpha', 0.4, 'DisplayName', sprintf('%dms', winMs(wi)));
        end
    end
    xlabel('Gap Length (ticks)');
    ylabel('Count');
    title('Miss Gap Distribution Per Window');
    legend('Location', 'northeast');
end
hold off;

% Combined gap distribution
subplot(1,2,2);
if ~isempty(missRuns)
    histogram(missRuns, 30, 'FaceColor', [0.8 0.2 0.2]);
    xlabel('Gap Length (ticks)');
    ylabel('Count');
    title(sprintf('Combined Miss Gaps (n=%d, max=%d)', length(missRuns), max(missRuns)));
    xline(mean(missRuns), 'b--', sprintf('mean=%.1f', mean(missRuns)), 'LineWidth', 1.5);
    xline(median(missRuns), 'g--', sprintf('median=%.1f', median(missRuns)), 'LineWidth', 1.5);
end

%% ========================================================================
%  FIGURE 5: Rolling Detection Rate
%% ========================================================================
figure('Name', 'Rolling Detection Rate', 'Position', [800 50 700 400]);

rollingWinSec = max(0.5, duration_s / 20);
rollingWinTicks = round(rollingWinSec / (tSec(2) - tSec(1)));
rollingWinTicks = max(rollingWinTicks, 10);

hold on;
% Per-window rolling rate
for wi = 1:numWindows
    rollingRate = movmean(double(detMat(:, wi)), rollingWinTicks) * 100;
    plot(tSec, rollingRate, 'Color', [getWinColor(wi, numWindows) 0.5], 'LineWidth', 0.8);
end
% Combined rolling rate (bold)
rollingAny = movmean(double(anyDet), rollingWinTicks) * 100;
plot(tSec, rollingAny, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Any window');
hold off;

xlabel('Time (s)');
ylabel('Detection Rate (%)');
title(sprintf('Rolling Detection Rate (window = %.1fs)', rollingWinSec));
ylim([0 105]);
grid on;

% Add legend for a subset of windows + combined
legendEntries = cell(1, numWindows + 1);
for wi = 1:numWindows
    legendEntries{wi} = sprintf('%dms', winMs(wi));
end
legendEntries{end} = 'Any (combined)';
legend(legendEntries, 'Location', 'best', 'FontSize', 7);

%% ========================================================================
%  FIGURE 6: Window Agreement
%% ========================================================================
figure('Name', 'Window Agreement', 'Position', [800 500 700 500]);

% Correlation matrix: how often do two windows agree?
corrMat = zeros(numWindows);
for a = 1:numWindows
    for b2 = 1:numWindows
        corrMat(a, b2) = sum(detMat(:,a) & detMat(:,b2)) / max(1, sum(detMat(:,a) | detMat(:,b2)));
    end
end

imagesc(corrMat);
colormap(gca, parula);
colorbar;
caxis([0 1]);
set(gca, 'XTick', 1:numWindows, 'XTickLabel', ...
    arrayfun(@(x) sprintf('%dms', x), winMs, 'UniformOutput', false));
set(gca, 'YTick', 1:numWindows, 'YTickLabel', ...
    arrayfun(@(x) sprintf('%dms', x), winMs, 'UniformOutput', false));
title('Window Agreement (IoU: intersection / union)');
xlabel('Window'); ylabel('Window');

% Add text values
for a = 1:numWindows
    for b2 = 1:numWindows
        text(b2, a, sprintf('%.2f', corrMat(a,b2)), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, ...
            'Color', ternaryVal(corrMat(a,b2) > 0.5, 'k', 'w'));
    end
end

fprintf('Analysis complete. 6 figures generated.\n');
end

%% =========================================================================
%                       LOCAL HELPERS
%% =========================================================================

function runs = getMissRuns(detected)
%GETMISSRUNS Lengths of consecutive false runs in a logical vector.
    runs = [];
    inRun = false;
    runLen = 0;
    for k = 1:length(detected)
        if ~detected(k)
            inRun = true;
            runLen = runLen + 1;
        else
            if inRun
                runs(end+1) = runLen; %#ok<AGROW>
                inRun = false;
                runLen = 0;
            end
        end
    end
    if inRun
        runs(end+1) = runLen; %#ok<AGROW>
    end
end

function c = getWinColor(idx, total)
%GETWINCOLOR Distinct color for each window index using hsv colormap.
    cmap = hsv(total);
    c = cmap(idx, :);
end

function v = ternaryVal(cond, a, b)
    if cond, v = a; else, v = b; end
end
