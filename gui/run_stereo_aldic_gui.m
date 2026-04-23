function run_stereo_aldic_gui()
%RUN_STEREO_ALDIC_GUI  GUI-driven 3D-Stereo-ALDIC pipeline.
%
%   Top-level entry point. Opens configureDICpara() to collect settings,
%   then runs the full pipeline (stereo match + temporal match L/R + 3D
%   reconstruction + strain + plot) with a progress dialog.
%
%   Algorithm and numerical behavior are identical to main_3D_ALDIC.m;
%   only the parameter-collection front-end differs.

%% ------------------------------------------------------------
%  Environment: addpath + mex compile (same as main_3D_ALDIC.m Section 1)
%  ------------------------------------------------------------
session_start = tic;

if isempty(getenv('MW_MINGW64_LOC'))
    default_mingw = 'C:\TDM-GCC-64';
    if exist(default_mingw, 'dir')
        setenv('MW_MINGW64_LOC', default_mingw);
    end
end

mex_src = 'ba_interp2_spline.cpp';
mex_bin = ['ba_interp2_spline.', mexext];
if ~exist(mex_bin, 'file') || (dir(mex_src).datenum > dir(mex_bin).datenum)
    try
        mex('-O', mex_src);
    catch ME
        error('Mex compile failed: %s', ME.message);
    end
end

addpath('./examples', './func', './func_quadtree/rbfinterp/', ...
        './plotFiles/', './func_quadtree', './func_quadtree/refinement', ...
        './plotFiles/export_fig-d966721/', './gui', './gui/gui_helpers');

%% ------------------------------------------------------------
%  1. Configuration GUI
%  ------------------------------------------------------------
DICpara = configureDICpara();
if isempty(DICpara)
    fprintf('Cancelled.\n');
    return;
end
validateDICpara(DICpara);

%% ------------------------------------------------------------
%  2. Progress dialog
%  ------------------------------------------------------------
fig = uifigure('Name', '3D-Stereo-ALDIC', 'Position', [300, 350, 500, 140]);
pd = uiprogressdlg(fig, 'Title', '3D-Stereo-ALDIC running...', ...
                   'Message', 'Starting...', 'Cancelable', 'on', ...
                   'Value', 0);
cleaner = onCleanup(@() cleanupProgress(fig, pd));

    function advance(frac, msg)
        if ~isvalid(pd), error('User cancelled.'); end
        pd.Value = frac;
        pd.Message = msg;
        drawnow;
    end

%% ------------------------------------------------------------
%  3. Load images + masks
%  ------------------------------------------------------------
advance(0.02, 'Loading images...');
[fileNameLeft, fileNameRight, imageLeft, imageRight] = ...
    loadStereoImages(DICpara.imgFolderLeft, DICpara.imgFolderRight);
numFrames = length(imageLeft);

advance(0.04, sprintf('Loading masks (%d frames)...', numFrames));
isIncremental = (DICpara.DICIncOrNot == 1);
[maskLeft, maskRight] = loadStereoMasks(DICpara.maskFolder, isIncremental, numFrames);

% Derive ROI from first left mask, if not preset
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

DICpara.ImgSize     = size(imageLeft{1});
DICpara.ImgRefMask  = double(maskLeft{1});
DICpara.LoadImgMethod = 1;
DICpara.ImgSeqIncUnit = 1;
DICpara.ImgSeqIncROIUpdateOrNot = isIncremental * 0 + (1-isIncremental) * 1;
DICpara.Subpb2FDOrFEM = 1;

%% ------------------------------------------------------------
%  4. Normalize + init storage
%  ------------------------------------------------------------
advance(0.08, 'Normalizing images...');
[imgNormalized_L, DICpara.gridxyROIRange] = funNormalizeImg(imageLeft,  DICpara.gridxyROIRange);
[imgNormalized_R, DICpara.gridxyROIRange] = funNormalizeImg(imageRight, DICpara.gridxyROIRange);

RD_L = initResultStorage(length(imgNormalized_L), isIncremental);
RD_R = initResultStorage(length(imgNormalized_R), isIncremental);

%% ------------------------------------------------------------
%  5. Stereo calibration
%  ------------------------------------------------------------
advance(0.10, 'Loading stereo calibration...');
[calibFolder, calibName, calibExt] = fileparts(DICpara.calibrationFile);
calibFileName = [calibName calibExt];
StereoInfo = loadCalibration(DICpara.calibrationMethod, calibFolder, calibFileName);

%% ------------------------------------------------------------
%  6. Stereo match
%  ------------------------------------------------------------
advance(0.15, 'Stereo matching (frame 1)...');
stereoMatchShapeOrder = 1;
[StereoInfo, RD_L, RD_R] = StereoMatch_STAQ( ...
    RD_L, RD_R, imgNormalized_L{1}, imgNormalized_R{1}, ...
    fileNameLeft, maskLeft{1}, maskRight{1}, DICpara, StereoInfo, stereoMatchShapeOrder);

%% ------------------------------------------------------------
%  7. Temporal match (left + right)
%  ------------------------------------------------------------
advance(0.25, 'Temporal match (left camera)...');
shapeFuncOrder = 1;
RD_L = TemporalMatch_quadtree_ST1(DICpara, fileNameLeft, maskLeft, ...
    imgNormalized_L, RD_L, StereoInfo, 'camera0', shapeFuncOrder);

advance(0.55, 'Temporal match (right camera)...');
RD_R = TemporalMatch_quadtree_ST1(DICpara, fileNameRight, maskRight, ...
    imgNormalized_R, RD_R, StereoInfo, 'notCamera0', shapeFuncOrder);

%% ------------------------------------------------------------
%  8. Organize + 3D reconstruction
%  ------------------------------------------------------------
advance(0.82, 'Organizing matched pairs + 3D reconstruction...');
matchedPairs = organizeMatchedPairds_quadtree_ST1(RD_L, RD_R, DICpara);
[FinalResult, ~] = stereoReconstruction_quadtree(matchedPairs, StereoInfo.cameraParams);

%% ------------------------------------------------------------
%  9. Section 6 — strain + plot
%  ------------------------------------------------------------
advance(0.88, 'Computing strain + plotting...');

% Ensure output folder exists
if ~exist(DICpara.outputFilePath, 'dir')
    mkdir(DICpara.outputFilePath);
end

% Optional specimen-frame transform
RotationMatrix = eye(3); TranslationMatrix = [0 0 0];
if DICpara.transformDisp == 1
    Base_Points2D = getBasePoints(imageLeft{1,1}', maskLeft{1,1}');
    [RotationMatrix, TranslationMatrix] = GetRTMatrix( ...
        RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM, Base_Points2D, FinalResult.Coordinates(1,:));
    DICpara.RforStrainCal = RotationMatrix;
end

strain_size = DICpara.strain_size;
strain_length = (strain_size - 1) * DICpara.winstepsize + 1;
VSG_length = strain_length + DICpara.winsize;
fprintf('strain_size=%dx%d pts | VSG=%dx%d px\n', strain_size, strain_size, VSG_length, VSG_length);

FinalResult.ResultStrainWorld{1} = 0;
FinalResult.Displacement_smooth(1,:) = FinalResult.Displacement(1,:);

nFramesMinus1 = length(imgNormalized_L) - 1;
for ImgSeqNum = 2 : length(imgNormalized_L)
    frac_strain = 0.88 + 0.10 * (ImgSeqNum - 1) / nFramesMinus1;
    advance(frac_strain, sprintf('Strain/plot frame %d/%d', ImgSeqNum, length(imgNormalized_L)));

    FullImageName_current = fullfile(fileNameLeft{2,ImgSeqNum}, fileNameLeft{1,ImgSeqNum});
    FullImageName_first   = fullfile(fileNameLeft{2,1}, fileNameLeft{1,1});

    % Smoothing
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
        FinalResult.Displacement_smooth(ImgSeqNum,:), FullImageName_first, ...
        ImgSeqNum, 'Local');

    DICpara = setPlottingParameters(DICpara, 'auto');
    if DICpara.DICIncOrNot == 0
        if DICpara.transformDisp == 0
            PlotdispQuadtreeMasks3D_acc_ST1(FinalResult.DisplacementNew(ImgSeqNum,:), ...
                RD_L.ResultDisp{ImgSeqNum-1,1}.U, RD_L.ResultFEMeshEachFrame{1,1}, ...
                FullImageName_first, FullImageName_current, DICpara, voidIndex);
        else
            PlotdispQuadtreeMasks3D_acc_ST1(FinalResult.Displacement(ImgSeqNum,:), ...
                RD_L.ResultDisp{ImgSeqNum-1,1}.U, RD_L.ResultFEMeshEachFrame{1,1}, ...
                FullImageName_first, FullImageName_current, DICpara, voidIndex);
        end
        [exx,exy,eyy,epMax,epMin,maxShear,vonMises,~,~] = PlotstrainQuadtreeMasks3D_acc_ST1( ...
            RD_L.ResultDisp{ImgSeqNum-1,1}.U, coefficients, voidIndex, ...
            RD_L.ResultFEMeshEachFrame{1,1}, FullImageName_first, FullImageName_current, DICpara);
    else
        if DICpara.transformDisp == 0
            PlotdispQuadtreeMasks3D_inc_ST1(FinalResult.DisplacementNew(ImgSeqNum,:), ...
                RD_L.ResultDisp_inc{ImgSeqNum-1,1}.U, RD_L.ResultDisp{ImgSeqNum-1,1}.U, ...
                RD_L.ResultFEMeshEachFrame{ImgSeqNum-1,1}, RD_L.ResultFEMeshEachFrame{1,1}, ...
                FullImageName_first, FullImageName_current, maskLeft{ImgSeqNum}, DICpara, voidIndex);
        else
            PlotdispQuadtreeMasks3D_inc_ST1(FinalResult.Displacement(ImgSeqNum,:), ...
                RD_L.ResultDisp_inc{ImgSeqNum-1,1}.U, RD_L.ResultDisp{ImgSeqNum-1,1}.U, ...
                RD_L.ResultFEMeshEachFrame{ImgSeqNum-1,1}, RD_L.ResultFEMeshEachFrame{1,1}, ...
                FullImageName_first, FullImageName_current, maskLeft{ImgSeqNum}, DICpara, voidIndex);
        end
        [exx,exy,eyy,epMax,epMin,maxShear,vonMises] = PlotstrainQuadtreeMasks3D_inc_ST1( ...
            RD_L.ResultDisp_inc{ImgSeqNum-1,1}.U, RD_L.ResultDisp{ImgSeqNum-1,1}.U, ...
            coefficients, voidIndex, RD_L.ResultFEMeshEachFrame{1,1}, ...
            RD_L.ResultFEMeshEachFrame{ImgSeqNum-1,1}, FullImageName_first, ...
            FullImageName_current, maskLeft{ImgSeqNum}, DICpara);
    end

    FinalResult.ResultStrainWorld{ImgSeqNum,1} = struct( ...
        'strain_exx', exx, 'strain_eyy', eyy, 'strain_exy', exy, ...
        'strain_principal_max', epMax, 'strain_principal_min', epMin, ...
        'strain_maxshear', maxShear, 'strain_vonMises', vonMises);

    if strcmp(DICpara.PlotMode, 'none') || strcmp(DICpara.PlotMode, 'save')
        close all;
    end
end

%% ------------------------------------------------------------
%  10. Save FinalResult + DICpara
%  ------------------------------------------------------------
advance(0.99, 'Saving results...');
save(fullfile(DICpara.outputFilePath, 'FinalResult.mat'), ...
     'FinalResult', 'DICpara', 'RD_L', 'RD_R', 'StereoInfo', '-v7.3');

advance(1.0, 'Done');
fprintf('Done. Total elapsed: %.1fs. Output: %s\n', ...
        toc(session_start), DICpara.outputFilePath);

end


% ================== Helpers ==================

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

function cleanupProgress(fig, pd)
if isvalid(pd), close(pd); end
if isvalid(fig), delete(fig); end
end
