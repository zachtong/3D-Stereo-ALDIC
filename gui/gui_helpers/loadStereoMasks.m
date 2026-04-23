function [maskLeft, maskRight] = loadStereoMasks(maskFolder, isIncremental, numFrames)
%LOADSTEREOMASKS  Non-interactive stereo mask loader.
%
%   [maskLeft, maskRight] = loadStereoMasks(maskFolder, isIncremental, numFrames)
%
%   Assumes masks are named with a '_0' suffix for left and '_1' for right
%   (matches the original ReadImage3DStereo_STAQ_mask_* convention), e.g.
%     maskFolder/frame_0000_0.tif   <- left mask, frame 1
%     maskFolder/frame_0000_1.tif   <- right mask, frame 1
%     maskFolder/frame_0001_0.tif   <- left mask, frame 2
%     ...
%
%   isIncremental : true  -> require numFrames masks per camera
%                   false -> only first mask is used (but more is OK)
%   numFrames     : expected count (sanity-check only)
%
%   Output: cell arrays of logical matrices (transposed).

if ~exist(maskFolder, 'dir')
    error('Mask folder does not exist: %s', maskFolder);
end
addpath(maskFolder);

leftFiles  = listMaskFiles(maskFolder, '_0');
rightFiles = listMaskFiles(maskFolder, '_1');

nLeft  = size(leftFiles, 2);
nRight = size(rightFiles, 2);

if nLeft == 0
    error('No left masks (*_0.*) found in %s', maskFolder);
end
if nRight == 0
    error('No right masks (*_1.*) found in %s', maskFolder);
end

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
