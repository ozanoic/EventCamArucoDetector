function mergedFiles = mergeAllResults(dataRoot)
%MERGEALLRESULTS Batch-merge v1+v2 result pairs found under dataRoot.
%
%  mergeAllResults()              % defaults to 'Data'
%  mergeAllResults('Data')
%  mergeAllResults('C:/.../Data')
%
%  Scans every immediate subfolder of dataRoot for files matching
%  <name>_results_v1.mat  AND  <name>_results_v2.mat.
%  For each matched pair, calls mergeResults and saves the combined
%  output as <name>_results.mat (no _vN suffix) in the same folder.
%
%  Returns a cell array of the merged output paths.

if nargin < 1 || isempty(dataRoot)
    dataRoot = 'Data';
end

if ~isfolder(dataRoot)
    error('mergeAllResults: folder not found: %s', dataRoot);
end

fprintf('\n### Scanning %s for v1/v2 result pairs ###\n', dataRoot);

d = dir(dataRoot);
subDirs = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));

mergedFiles = {};
pairCount = 0;
skipCount = 0;

for i = 1:length(subDirs)
    folder = fullfile(dataRoot, subDirs(i).name);

    v1Files = dir(fullfile(folder, '*_results_v1.mat'));
    if isempty(v1Files), continue; end

    for k = 1:length(v1Files)
        v1Path = fullfile(folder, v1Files(k).name);
        % Derive the base name (strip "_v1.mat")
        base = regexprep(v1Files(k).name, '_v1\.mat$', '');
        v2Name = [base '_v2.mat'];
        v2Path = fullfile(folder, v2Name);

        if ~isfile(v2Path)
            fprintf('  [skip] %s  (no matching v2)\n', v1Files(k).name);
            skipCount = skipCount + 1;
            continue;
        end

        outName = [base '.mat'];         % e.g. marker_z2_cross_high_results.mat
        outPath = fullfile(folder, outName);

        fprintf('\n>>> Merging pair in %s:\n', folder);
        fprintf('    v1: %s\n', v1Files(k).name);
        fprintf('    v2: %s\n', v2Name);
        fprintf('    -> %s\n', outName);

        try
            mergeResults({v1Path, v2Path}, outPath);
            mergedFiles{end+1} = outPath; %#ok<AGROW>
            pairCount = pairCount + 1;
        catch ME
            fprintf('  [error] %s\n', ME.message);
            skipCount = skipCount + 1;
        end
    end
end

fprintf('\n### Done: %d pair(s) merged, %d skipped ###\n', pairCount, skipCount);
end
