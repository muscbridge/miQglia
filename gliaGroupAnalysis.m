%% Skeleton and Fractal Analysis of Micrographs (Group Based)
% This script runs ramification analysis on a group. Your group should be a
% folder containing images of the same ROI acquired with the same imaging
% parameters.
%
% Author: Siddhartha Dhiman
%
% Parameters
% ----------
%   groupDir : str
%       Path to input image file (.tif)
%   outDir : str
%       Path to output directory write out files
% Returns
% -------
%   props : struct
%       Struct containing quantified parameters
%--------------------------------------------------------------------------
%% Enter Paths
groupDir =  'D:\Datasets\Falangola_Morphology\data\Benchmark';
outDir =    'D:\Datasets\Falangola_Morphology\data\Benchmark\Outputs\Version_2';

%% Begin
if ~isdir(groupDir)
    error('Input directory %s does not exist', groupDir);
end
if ~isdir(outDir)
    error('Output directory %s is not a directory. Please specify a valid directory.', outpath);
end
if ~exist(outDir, 'dir')
    error('Output directory %s does not exist.', outDir)
end

%% Search for Image Files
files = dir(fullfile(groupDir, ['**', filesep, '*', '.tif']));

%% Recursively Run ramificationStats
f = waitbar(0, 'Starting');
groupStats = table;
fnum = length(files);
for i = 1:fnum
    try
        fpath = fullfile(files(i).folder, files(i).name);
        stats = struct2table(ramificationStatsTest(fpath, outDir));
        groupStats = vertcat(groupStats, stats);
    catch
        continue
    end
    waitbar(i/fnum, f, strrep(sprintf('Progress: %d %%\n File: %s', floor(i/fnum*100), files(i).name), '_', '\_'));
end
writetable(groupStats, fullfile(outDir, 'Group_Stats.csv'));
