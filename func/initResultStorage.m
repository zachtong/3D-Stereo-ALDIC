function RD = initResultStorage(numImagePairs, isIncremental)
%INITRESULTSTORAGE  Pre-allocate result containers for one camera.
%
%   RD = initResultStorage(numImagePairs, isIncremental)
%
%   Inputs:
%       numImagePairs - length(imgNormalized_*) (i.e. total frame count)
%       isIncremental - true for DICpara.DICIncOrNot == 1, false for acc
%
%   Output:
%       RD - struct with ResultDisp, ResultDefGrad, ResultFEMeshEachFrame,
%            ResultFEMesh. Incremental mode additionally adds
%            ResultDisp_inc.

numFrames = numImagePairs - 1;
RD.ResultDisp    = cell(numFrames, 1);
RD.ResultDefGrad = cell(numFrames, 1);

if isIncremental
    RD.ResultFEMeshEachFrame = cell(numFrames, 1);
    RD.ResultFEMesh          = cell(numFrames, 1);
    RD.ResultDisp_inc        = cell(numFrames, 1);
else
    RD.ResultFEMeshEachFrame = cell(1, 1);
    RD.ResultFEMesh          = cell(1, 1);
end
end
