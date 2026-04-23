function [fileNameLeft, fileNameRight, imageLeft, imageRight] = loadStereoImages(leftFolder, rightFolder)
%LOADSTEREOIMAGES  Non-interactive stereo image loader.
%
%   [fileNameLeft, fileNameRight, imageLeft, imageRight] = loadStereoImages(leftFolder, rightFolder)
%
%   Sorts image files in each folder by name and loads them as a cell
%   array of double grayscale matrices (transposed to match the DIC
%   pipeline's [x, y] = [row, col] convention).
%
%   Supported extensions: jpg, jpeg, tif, tiff, bmp, png, jp2.
%
%   Output format matches ReadImage3DStereo_STAQ.m so downstream code
%   (StereoMatch_STAQ, TemporalMatch_quadtree_ST1) can use the results
%   directly.
%
%       fileNameLeft{1,i}  = filename  (row 1 = name)
%       fileNameLeft{2,i}  = folder    (row 2 = directory path)
%       imageLeft{i}       = double(H x W) grayscale, transposed

if ~exist(leftFolder, 'dir')
    error('Left images folder does not exist: %s', leftFolder);
end
if ~exist(rightFolder, 'dir')
    error('Right images folder does not exist: %s', rightFolder);
end

addpath(leftFolder);
if ~strcmp(leftFolder, rightFolder)
    addpath(rightFolder);
end

fileNameLeft  = listImageFiles(leftFolder);
fileNameRight = listImageFiles(rightFolder);

numLeft  = size(fileNameLeft, 2);
numRight = size(fileNameRight, 2);
if numLeft == 0
    error('No image files found in %s', leftFolder);
end
if numLeft ~= numRight
    warning('Left has %d images, right has %d — pipeline expects same count.', numLeft, numRight);
end

numImages = numLeft;
imageLeft  = cell(numImages, 1);
imageRight = cell(numImages, 1);
for i = 1:numImages
    imageLeft{i}  = loadOneImage(fullfile(fileNameLeft{2,i},  fileNameLeft{1,i}));
    imageRight{i} = loadOneImage(fullfile(fileNameRight{2,i}, fileNameRight{1,i}));
end
end


function fileCell = listImageFiles(folder)
% Returns a 2xN cell like struct2cell(dir(...)), with filename (row 1) +
% directory (row 2). Concatenates all supported extensions, case-insensitive.
exts = {'jpg','jpeg','tif','tiff','bmp','png','jp2'};
entries = [];
for k = 1:length(exts)
    entries = [entries; dir(fullfile(folder, ['*.' exts{k}]))];  %#ok<AGROW>
end
% Sort by name for deterministic ordering
if ~isempty(entries)
    [~, sortIdx] = sort({entries.name});
    entries = entries(sortIdx);
end
fileCell = struct2cell(entries);
end


function img = loadOneImage(path)
raw = imread(path);
[~, ~, nCh] = size(raw);
if nCh == 4
    raw = rgb2gray(raw(:,:,1:3));  % drop alpha
elseif nCh == 3
    raw = rgb2gray(raw);
end
img = double(raw)';  % transpose for pipeline convention
end
