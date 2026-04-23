function [maskLeft, maskRight] = loadStereoMasks(leftFolderOrSingle, rightFolderOrMode, varargin)
%LOADSTEREOMASKS  Non-interactive stereo mask loader.
%
%   Two calling conventions:
%
%   1) Dual-folder (preferred):
%       [maskLeft, maskRight] = loadStereoMasks(leftFolder, rightFolder, isIncremental, numFrames)
%      leftFolder  holds left-camera masks (any name)
%      rightFolder holds right-camera masks (any name); sorted by filename
%
%   2) Single-folder with suffix convention:
%       [maskLeft, maskRight] = loadStereoMasks(singleFolder, isIncremental, numFrames)
%      singleFolder contains *_0.ext (left) and *_1.ext (right) files.
%
%   Output: cell arrays of logical matrices, transposed to match the DIC
%   pipeline's (x,y)=(row,col) convention.

% Dispatch on argument pattern
if nargin >= 3 && isnumeric(rightFolderOrMode) == false && ischar(rightFolderOrMode)
    % Dual-folder mode
    leftFolder  = leftFolderOrSingle;
    rightFolder = rightFolderOrMode;
    isIncremental = varargin{1};
    numFrames     = varargin{2};
    [maskLeft, maskRight] = loadDualFolder(leftFolder, rightFolder, isIncremental, numFrames);
else
    % Single-folder mode
    singleFolder = leftFolderOrSingle;
    isIncremental = rightFolderOrMode;
    numFrames     = varargin{1};
    [maskLeft, maskRight] = loadSingleFolder(singleFolder, isIncremental, numFrames);
end
end


function [maskLeft, maskRight] = loadDualFolder(leftFolder, rightFolder, isIncremental, numFrames)
if ~exist(leftFolder, 'dir'),  error('Left mask folder does not exist: %s', leftFolder);  end
if ~exist(rightFolder, 'dir'), error('Right mask folder does not exist: %s', rightFolder); end
addpath(leftFolder);
if ~strcmp(leftFolder, rightFolder), addpath(rightFolder); end

leftFiles  = listAllImageFiles(leftFolder);
rightFiles = listAllImageFiles(rightFolder);

nLeft  = size(leftFiles, 2);
nRight = size(rightFiles, 2);
if nLeft == 0,  error('No masks found in %s', leftFolder);  end
if nRight == 0, error('No masks found in %s', rightFolder); end

[maskLeft, maskRight] = loadAndTrim(leftFiles, rightFiles, isIncremental, numFrames);
end


function [maskLeft, maskRight] = loadSingleFolder(singleFolder, isIncremental, numFrames)
if ~exist(singleFolder, 'dir')
    error('Mask folder does not exist: %s', singleFolder);
end
addpath(singleFolder);

leftFiles  = listMaskFiles(singleFolder, '_0');
rightFiles = listMaskFiles(singleFolder, '_1');

nLeft  = size(leftFiles, 2);
nRight = size(rightFiles, 2);

if nLeft == 0,  error('No left masks (*_0.*) found in %s', singleFolder);  end
if nRight == 0, error('No right masks (*_1.*) found in %s', singleFolder); end

[maskLeft, maskRight] = loadAndTrim(leftFiles, rightFiles, isIncremental, numFrames);
end


function [maskLeft, maskRight] = loadAndTrim(leftFiles, rightFiles, isIncremental, numFrames)
nLeft  = size(leftFiles, 2);
nRight = size(rightFiles, 2);

if isIncremental
    if nLeft < numFrames || nRight < numFrames
        warning('Incremental mode needs one mask per frame (%d expected); got %d left, %d right.', ...
            numFrames, nLeft, nRight);
    end
    numToLoad = min([nLeft, nRight, numFrames]);
else
    % Accumulative mode: first mask only, but load all present for safety.
    numToLoad = min(nLeft, nRight);
end

maskLeft  = cell(numToLoad, 1);
maskRight = cell(numToLoad, 1);
for i = 1:numToLoad
    maskLeft{i}  = loadOneMask(fullfile(leftFiles{2,i},  leftFiles{1,i}));
    maskRight{i} = loadOneMask(fullfile(rightFiles{2,i}, rightFiles{1,i}));
end
end


function fileCell = listAllImageFiles(folder)
% Any supported extension, sorted by name.
exts = {'jpg','jpeg','tif','tiff','bmp','png','jp2'};
entries = [];
for k = 1:length(exts)
    entries = [entries; dir(fullfile(folder, ['*.' exts{k}]))];  %#ok<AGROW>
end
if ~isempty(entries)
    [~, sortIdx] = sort({entries.name});
    entries = entries(sortIdx);
end
fileCell = struct2cell(entries);
end


function fileCell = listMaskFiles(folder, suffixTag)
% Lists files matching *<suffixTag>.<ext> for any supported image ext.
exts = {'jpg','jpeg','tif','tiff','bmp','png','jp2'};
entries = [];
for k = 1:length(exts)
    entries = [entries; dir(fullfile(folder, ['*' suffixTag '.' exts{k}]))];  %#ok<AGROW>
end
if ~isempty(entries)
    [~, sortIdx] = sort({entries.name});
    entries = entries(sortIdx);
end
fileCell = struct2cell(entries);
end


function m = loadOneMask(path)
raw = imread(path);
[~, ~, nCh] = size(raw);
if nCh == 4
    raw = rgb2gray(raw(:,:,1:3));
elseif nCh == 3
    raw = rgb2gray(raw);
end
m = logical(raw)';
end
