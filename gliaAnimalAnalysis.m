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
groupDir =  'D:\Datasets\Falangola_Morphology\data\Sub\Older_mice\TG';
outDir =    'D:\Datasets\Falangola_Morphology\data\Sub\By_Animal\Final_Run\TG';

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
fname = extractfield(files, 'name');
aname = cell(length(fname), 1);
for i = 1:length(fname)
    tmp = strsplit(fname{i}, '_');
    aname{i} = tmp{1};
end
aname_ = unique(aname);
% mid = cell(length(aname), 1);
% for i = 1:length(aname)
%     tmp = strsplit(aname{i}, 'm');
%     mid{i} = ['m' tmp{2}];
% end
% mid_ = unique(mid);
for i = 1:length(aname_)
    idx = logical(zeros(length(aname),1));
    for j = 1:length(aname)
        idx(j) = contains(aname{j}, aname_{i});
    end
    files_ = files(idx);
    odir = fullfile(outDir, aname_{i});
    if ~exist(odir, 'dir')
       mkdir(odir)
    end
    %% Recursively Run ramificationStats
    groupStats = table;
    for k = 1:length(files_)
        try
        fpath = fullfile(files_(k).folder, files_(k).name);
        stats = struct2table(ramificationStats(fpath, odir));
        groupStats = vertcat(groupStats, stats);
        catch
            continue
        end
    end
    writetable(groupStats, fullfile(outDir, [aname_{i} '_stats.csv']));
end
