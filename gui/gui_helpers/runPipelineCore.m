function results = runPipelineCore(DICpara, logFcn, progressFcn)
%RUNPIPELINECORE  Run the full 3D-Stereo-ALDIC pipeline non-interactively.
%
%   results = runPipelineCore(DICpara, logFcn, progressFcn)
%
%   Inputs:
%       DICpara     - fully populated DICpara struct (all fields set by GUI
%                     or script caller)
%       logFcn      - @(msg) function to append msg to the GUI log (fallback
%                     to fprintf if empty)
%       progressFcn - @(frac, msg) function to update progress; empty = noop
%
%   Returns a struct with:
%       results.FinalResult, results.RD_L, results.RD_R, results.StereoInfo
%       results.imageLeft, results.imageRight, results.maskLeft, results.maskRight
%       results.fileNameLeft, results.fileNameRight
%       results.coefficientsPerFrame   (cell array, indexed by ImgSeqNum)
%       results.strainPerFrame         (cell array, indexed by ImgSeqNum)
%
%   When DICpara.PlotMode == 'none', per-frame plotting is skipped (saves
%   time); strain is still computed via computeStrain3D directly.

if nargin < 2 || isempty(logFcn),      logFcn = @(m) fprintf('%s\n', m); end
if nargin < 3 || isempty(progressFcn), progressFcn = @(~,~) []; end

%% 1. Load images + masks
progressFcn(0.02, 'Loading images...');
logFcn(sprintf('Loading images from %s', DICpara.imgFolderLeft));
[fileNameLeft, fileNameRight, imageLeft, imageRight] = ...
    loadStereoImages(DICpara.imgFolderLeft, DICpara.imgFolderRight);
numFrames = length(imageLeft);
logFcn(sprintf('Loaded %d frame pairs, size = %dx%d', numFrames, size(imageLeft{1},1), size(imageLeft{1},2)));

progressFcn(0.04, 'Loading masks...');
isIncremental = (DICpara.DICIncOrNot == 1);
[maskLeft, maskRight] = loadMasksFromDICpara(DICpara, isIncremental, numFrames);
logFcn(sprintf('Loaded %d masks per camera', length(maskLeft)));

% ROI derivation
if isempty(DICpara.gridxyROIRange.gridx) || isempty(DICpara.gridxyROIRange.gridy)
    if DICpara.DICIncOrNot == 0
        [rowsInMask, colsInMask] = find(maskLeft{1});
        DICpara.gridxyROIRange.gridx = [min(rowsInMask), max(rowsInMask)];
        DICpara.gridxyROIRange.gridy = [min(colsInMask), max(colsInMask)];
    else
        DICpara.gridxyROIRange.gridx = [1, size(maskLeft{1}, 1)];
        DICpara.gridxyROIRange.gridy = [1, size(maskLeft{1}, 2)];
    end
end

DICpara.ImgSize                 = size(imageLeft{1});
DICpara.ImgRefMask              = double(maskLeft{1});
DICpara.LoadImgMethod           = 1;
DICpara.ImgSeqIncUnit           = 1;
DICpara.ImgSeqIncROIUpdateOrNot = ~isIncremental;
DICpara.Subpb2FDOrFEM           = 1;

%% 2. Normalize + init storage
progressFcn(0.08, 'Normalizing images...');
[imgNormalized_L, DICpara.gridxyROIRange] = funNormalizeImg(imageLeft,  DICpara.gridxyROIRange);
[imgNormalized_R, DICpara.gridxyROIRange] = funNormalizeImg(imageRight, DICpara.gridxyROIRange);

RD_L = initResultStorage(length(imgNormalized_L), isIncremental);
RD_R = initResultStorage(length(imgNormalized_R), isIncremental);

%% 3. Calibration
progressFcn(0.10, 'Loading stereo calibration...');
[calibFolder, calibName, calibExt] = fileparts(DICpara.calibrationFile);
calibFileName = [calibName calibExt];
StereoInfo = loadCalibration(DICpara.calibrationMethod, calibFolder, calibFileName);
logFcn('Stereo calibration loaded');

%% 4. Stereo match
progressFcn(0.15, 'Stereo matching (frame 1)...');
t0 = tic;
DICpara.showImgOrNot = 0;  % force silent for GUI
[StereoInfo, RD_L, RD_R] = StereoMatch_STAQ( ...
    RD_L, RD_R, imgNormalized_L{1}, imgNormalized_R{1}, ...
    fileNameLeft, maskLeft{1}, maskRight{1}, DICpara, StereoInfo, 1);
logFcn(sprintf('Stereo match: %.1fs', toc(t0)));

%% 5. Temporal match
progressFcn(0.25, 'Temporal match (left camera)...');
t0 = tic;
RD_L = TemporalMatch_quadtree_ST1(DICpara, fileNameLeft, maskLeft, ...
    imgNormalized_L, RD_L, StereoInfo, 'camera0', 1);
logFcn(sprintf('Temporal (L): %.1fs', toc(t0)));

progressFcn(0.55, 'Temporal match (right camera)...');
t0 = tic;
RD_R = TemporalMatch_quadtree_ST1(DICpara, fileNameRight, maskRight, ...
    imgNormalized_R, RD_R, StereoInfo, 'notCamera0', 1);
logFcn(sprintf('Temporal (R): %.1fs', toc(t0)));

%% 6. 3D reconstruction
progressFcn(0.82, '3D reconstruction...');
t0 = tic;
matchedPairs = organizeMatchedPairds_quadtree_ST1(RD_L, RD_R, DICpara);
[FinalResult, reprojErr] = stereoReconstruction_quadtree(matchedPairs, StereoInfo.cameraParams);
logFcn(sprintf('3D reconstruction: %.1fs | mean reproj err frame 1: %.3f px', ...
    toc(t0), mean(reprojErr{1}, 'omitnan')));

%% 7. Strain (always computed; plotting gated by PlotMode)
progressFcn(0.88, 'Computing strain...');

if ~exist(DICpara.outputFilePath, 'dir'), mkdir(DICpara.outputFilePath); end

% Optional specimen-frame transform
RotationMatrix = eye(3); TranslationMatrix = [0 0 0];
if DICpara.transformDisp == 1 && ~isfield(DICpara, 'RforStrainCal')
    Base_Points2D = getBasePoints(imageLeft{1,1}', maskLeft{1,1}');
    [RotationMatrix, TranslationMatrix] = GetRTMatrix( ...
        RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM, Base_Points2D, FinalResult.Coordinates(1,:));
    DICpara.RforStrainCal = RotationMatrix;
elseif isfield(DICpara, 'RforStrainCal') && ~isempty(DICpara.RforStrainCal)
    RotationMatrix = DICpara.RforStrainCal;
end

strain_size = DICpara.strain_size;
strain_length = (strain_size - 1) * DICpara.winstepsize + 1;
VSG_length = strain_length + DICpara.winsize;
logFcn(sprintf('strain_size=%dx%d pts | VSG=%dx%d px', strain_size, strain_size, VSG_length, VSG_length));

FinalResult.ResultStrainWorld{1} = 0;
FinalResult.Displacement_smooth(1,:) = FinalResult.Displacement(1,:);

coefficientsPerFrame = cell(length(imgNormalized_L), 1);
strainPerFrame       = cell(length(imgNormalized_L), 1);

nFramesMinus1 = length(imgNormalized_L) - 1;
silent = strcmpi(DICpara.PlotMode, 'none');

for ImgSeqNum = 2 : length(imgNormalized_L)
    frac_strain = 0.88 + 0.10 * (ImgSeqNum - 1) / nFramesMinus1;
    progressFcn(frac_strain, sprintf('Strain frame %d/%d', ImgSeqNum, length(imgNormalized_L)));
    t0 = tic;

    % 3D smoothing
    try
        for SmoothTimes = 1:DICpara.Smooth3DTimes
            FinalResult.Displacement_smooth(ImgSeqNum,:) = funSmoothDisp_Quadtree( ...
                FinalResult.Displacement(ImgSeqNum,:), ...
                RD_L.ResultFEMeshEachFrame{1}, DICpara);
        end
    catch
    end
    if DICpara.Smooth3DTimes == 0
        FinalResult.Displacement_smooth(ImgSeqNum,:) = FinalResult.Displacement(ImgSeqNum,:);
    end

    FinalResult = ConvertCoorAndDisp(FinalResult, RotationMatrix, TranslationMatrix);

    [DICpara, coefficients, voidIndex] = PlaneFit3_Quadtree( ...
        strain_size, DICpara, RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM, ...
        FinalResult.Coordinates(1,:), FinalResult.Coordinates(ImgSeqNum,:), ...
        FinalResult.Displacement_smooth(ImgSeqNum,:), '', ...
        ImgSeqNum, 'Local');
    coefficientsPerFrame{ImgSeqNum} = coefficients;

    % Always compute strain components (cheap, no figure)
    [exx, eyy, exy, e1, e2, maxShear, vonMises, dwdx, dwdy] = computeStrain3D(coefficients);
    strainPerFrame{ImgSeqNum} = struct( ...
        'strain_exx', exx, 'strain_eyy', eyy, 'strain_exy', exy, ...
        'strain_principal_max', e1, 'strain_principal_min', e2, ...
        'strain_maxshear', maxShear, 'strain_vonMises', vonMises, ...
        'dwdx', dwdx, 'dwdy', dwdy, 'voidIndex', voidIndex);
    FinalResult.ResultStrainWorld{ImgSeqNum,1} = strainPerFrame{ImgSeqNum};

    if ~silent
        % Legacy path: still call plotters (produces figures per frame)
        fullFirst = fullfile(fileNameLeft{2,1}, fileNameLeft{1,1});
        fullCurr  = fullfile(fileNameLeft{2,ImgSeqNum}, fileNameLeft{1,ImgSeqNum});
        DICpara = setPlottingParameters(DICpara, 'auto');
        if DICpara.DICIncOrNot == 0
            if DICpara.transformDisp == 0
                PlotdispQuadtreeMasks3D_acc_ST1(FinalResult.DisplacementNew(ImgSeqNum,:), ...
                    RD_L.ResultDisp{ImgSeqNum-1,1}.U, RD_L.ResultFEMeshEachFrame{1,1}, ...
                    fullFirst, fullCurr, DICpara, voidIndex);
            else
                PlotdispQuadtreeMasks3D_acc_ST1(FinalResult.Displacement(ImgSeqNum,:), ...
                    RD_L.ResultDisp{ImgSeqNum-1,1}.U, RD_L.ResultFEMeshEachFrame{1,1}, ...
                    fullFirst, fullCurr, DICpara, voidIndex);
            end
        else
            if DICpara.transformDisp == 0
                PlotdispQuadtreeMasks3D_inc_ST1(FinalResult.DisplacementNew(ImgSeqNum,:), ...
                    RD_L.ResultDisp_inc{ImgSeqNum-1,1}.U, RD_L.ResultDisp{ImgSeqNum-1,1}.U, ...
                    RD_L.ResultFEMeshEachFrame{ImgSeqNum-1,1}, RD_L.ResultFEMeshEachFrame{1,1}, ...
                    fullFirst, fullCurr, maskLeft{ImgSeqNum}, DICpara, voidIndex);
            else
                PlotdispQuadtreeMasks3D_inc_ST1(FinalResult.Displacement(ImgSeqNum,:), ...
                    RD_L.ResultDisp_inc{ImgSeqNum-1,1}.U, RD_L.ResultDisp{ImgSeqNum-1,1}.U, ...
                    RD_L.ResultFEMeshEachFrame{ImgSeqNum-1,1}, RD_L.ResultFEMeshEachFrame{1,1}, ...
                    fullFirst, fullCurr, maskLeft{ImgSeqNum}, DICpara, voidIndex);
            end
        end
    end

    logFcn(sprintf('Frame %d/%d: strain computed (%.1fs)', ...
        ImgSeqNum, length(imgNormalized_L), toc(t0)));
end

%% 8. Save
progressFcn(0.99, 'Saving FinalResult.mat...');
save(fullfile(DICpara.outputFilePath, 'FinalResult.mat'), ...
     'FinalResult', 'DICpara', 'RD_L', 'RD_R', 'StereoInfo', '-v7.3');
logFcn(sprintf('Saved FinalResult.mat to %s', DICpara.outputFilePath));

progressFcn(1.0, 'Done');

%% Build return struct
results = struct();
results.FinalResult         = FinalResult;
results.RD_L                = RD_L;
results.RD_R                = RD_R;
results.StereoInfo          = StereoInfo;
results.DICpara             = DICpara;
results.imageLeft           = imageLeft;
results.imageRight          = imageRight;
results.maskLeft            = maskLeft;
results.maskRight           = maskRight;
results.fileNameLeft        = fileNameLeft;
results.fileNameRight       = fileNameRight;
results.coefficientsPerFrame = coefficientsPerFrame;
results.strainPerFrame       = strainPerFrame;
end


% ===== Local helpers =====

function [maskLeft, maskRight] = loadMasksFromDICpara(DICpara, isIncremental, numFrames)
% Support both (a) single-folder with _0/_1 suffix and (b) dual folders.
hasDual = isfield(DICpara, 'maskFolderLeft') && ~isempty(DICpara.maskFolderLeft) && ...
          isfield(DICpara, 'maskFolderRight') && ~isempty(DICpara.maskFolderRight);
if hasDual
    [maskLeft, maskRight] = loadStereoMasks( ...
        DICpara.maskFolderLeft, DICpara.maskFolderRight, isIncremental, numFrames);
elseif isfield(DICpara, 'maskFolder') && ~isempty(DICpara.maskFolder)
    [maskLeft, maskRight] = loadStereoMasks(DICpara.maskFolder, isIncremental, numFrames);
else
    error('DICpara must specify maskFolder (single) or maskFolderLeft+maskFolderRight (dual).');
end
end

function StereoInfo = loadCalibration(method, folder, fileName)
StereoInfo = struct();
switch method
    case 0, StereoInfo = cameraParamsFormatConvertFromMatlabCV;
    case 1, StereoInfo.cameraParams = cameraParamsFormatConvertFromMatchID(folder, fileName);
    case 2, StereoInfo.cameraParams = cameraParamsFormatConvertFromMMC(folder, fileName);
    case 3, StereoInfo.cameraParams = cameraParamsFormatConvertFromDICe(folder, fileName);
    case 4, StereoInfo.cameraParams = cameraParamsFormatConvertFromOpenCorrFormat(folder, fileName);
    otherwise, error('Unknown calibration method: %d', method);
end
end
