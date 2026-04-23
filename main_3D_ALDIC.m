%% ===================================================================
% 3D Stereo Adaptive Mesh Augmented Lagrangian Digital Image Correlation (3D-Stereo-ALDIC)
% ===================================================================
%
% DESCRIPTION:
%   This program implements a 3D stereo DIC method using adaptive quadtree mesh
%   for full-field deformation measurement. The algorithm includes stereo
%   calibration, temporal matching, 3D reconstruction, and strain/stress computation.
%
% FEATURES:
%   - Stereo camera calibration and matching
%   - Adaptive quadtree mesh refinement
%   - Temporal DIC tracking
%   - 3D reconstruction and displacement field computation
%   - Strain and stress field computation and visualization
%
% MAIN WORKFLOW:
%   1. Environment setup and image loading
%   2. DIC parameter initialization
%   3. Stereo calibration and matching
%   4. Temporal matching using quadtree mesh
%   5. 3D reconstruction
%   6. Strain/stress computation and visualization
%
% AUTHOR:
%   Zixiang (Zach) Tong (zachtong@utexas.edu)
%   Jin Yang (jin.yang@austin.utexas.edu)
%   University of Texas at Austin
%
% VERSION:
%   1.0 - June 2025
%
% REFERENCES:
%   [Add relevant papers/publications]
%
% LICENSE:
%   [Add license information]
%
% ===================================================================
%% Section 1: Clear MATLAB environment & mex set up Spline interpolation
close all; clear; clc; clearvars -global
fprintf('------------ Section 1 Start ------------ \n')
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
try
    mex -O ba_interp2_spline.cpp;
    warning('off');
    fprintf('Mex compilation successful.\n');
catch ME
    fprintf('Mex compilation failed: %s\n', ME.message);
    fprintf('Please check complier installation and path.\n');
end

addpath("./examples","./func",'./func_quadtree/rbfinterp/','./plotFiles/','./func_quadtree','./func_quadtree/refinement','./plotFiles/export_fig-d966721/');
% TODO: addpath("./YOUR IMAGE FOLDER");
fprintf('------------ Section 1 Done ------------ \n \n')

%% Section 2: Load DIC parameters and set up DIC parameters
%-------- Notes ----------------
% Strategy 1 needs Left and Right masks
% Strategy 2 needs Left masks only (not provided here)
% Incremental mode needs all updated masks
% Accumulative mode needs only first mask
%-------------------------------
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images and masks ======
% Load DIC raw images
[fileNameLeft, fileNameRight, imageLeft,imageRight, LoadImgMethod] = ReadImage3DStereo_STAQ;
DICpara = setDICParas_IncOrNot(size(fileNameRight,2));

% Load DIC masks
if DICpara.DICIncOrNot == 0
    [~, ~, maskLeft, maskRight] = ReadImage3DStereo_STAQ_mask_acc(imageLeft{1},imageRight{1});
elseif DICpara.DICIncOrNot == 1
    [~, ~, maskLeft, maskRight] = ReadImage3DStereo_STAQ_mask_inc(imageLeft{1},imageRight{1});
end

% ====== Set up DIC paras =====
DICpara = setDICParas_STAQ(DICpara, fileNameLeft,imageLeft,maskLeft,maskRight,LoadImgMethod);
try
    fprintf('The finest element size in the adaptive quadtree mesh is %d.\n', DICpara.winsizeMin);
catch
    DICpara.winsizeMin = 8; % Assign the finest element size in the adaptive quadtree mesh
    fprintf('The finest element size in the adaptive quadtree mesh is set to 8 by default.\n');
end

% ====== Normalize images: fNormalized = (f-f_avg)/(f_std) ======
[imgNormalized_L,DICpara.gridxyROIRange] = funNormalizeImg(imageLeft,DICpara.gridxyROIRange);
[imgNormalized_R,DICpara.gridxyROIRange] = funNormalizeImg(imageRight,DICpara.gridxyROIRange);

% ====== Initialize variable storage ======
RD_L.ResultDisp = cell(length(imgNormalized_L)-1,1);    RD_L.ResultDefGrad = cell(length(imgNormalized_L)-1,1);
% RD_L.ResultStrainWorld = cell(length(imgNormalized_L)-1,1);  RD_L.ResultStressWorld = cell(length(imgNormalized_L)-1,1);
if DICpara.DICIncOrNot == 0
    RD_L.ResultFEMeshEachFrame = cell(1,1); % Quadtree mesh
    RD_L.ResultFEMesh = cell(1,1); % Non-Quadtree mesh.
else
    RD_L.ResultFEMeshEachFrame = cell(length(imgNormalized_L)-1,1); % Quadtree mesh
    RD_L.ResultFEMesh = cell(length(imgNormalized_L)-1,1); % Non-Quadtree mesh.
    RD_L.ResultDisp_inc = cell(length(imgNormalized_L)-1,1);
end

RD_R.ResultDisp = cell(length(imgNormalized_R)-1,1);    RD_R.ResultDefGrad = cell(length(imgNormalized_R)-1,1);
%RD_R.ResultStrainWorld = cell(length(imgNormalized_R)-1,1);  RD_R.ResultStressWorld = cell(length(imgNormalized_R)-1,1);
if DICpara.DICIncOrNot == 0
    RD_R.ResultFEMeshEachFrame = cell(1,1); % Quadtree mesh
    RD_R.ResultFEMesh = cell(1,1); % Non-Quadtree mesh.
else
    RD_R.ResultFEMeshEachFrame = cell(length(imgNormalized_R)-1,1);
    RD_R.ResultFEMesh = cell(length(imgNormalized_L)-1,1);
    RD_R.ResultDisp_inc = cell(length(imgNormalized_L)-1,1);
end

DICpara.ImgRefMask = double(maskLeft{1});

fprintf('------------ Section 2 Done ------------ \n \n')

%% Section 3.1 Stereo Calibration
calib_method = funParaInput('CalibrationMethod');
switch calib_method
    case 0
        % Install Computer Vision Toolbox and use the Stereo Camera
        % Calibrator app to calibrate your cameras; then export via
        % cameraParamsFormatConvertFromMatlabCV below.
        StereoInfo = cameraParamsFormatConvertFromMatlabCV;

    case 1
        % Import Calib. results from MatchID Calibrator
        [CalibrationFile, CalibrationFilepath]  = uigetfile({'*.caldat'}, 'choose a *.caldat file');
        StereoInfo.cameraParams = cameraParamsFormatConvertFromMatchID(CalibrationFilepath, CalibrationFile);
        clear CalibrationFile CalibrationFilepath
    case 2
        % Import Calib. results from MCC (Zhuoyi Yin's Calibrator)
        [CalibrationFile, CalibrationFilepath]  = uigetfile({'*.mat'}, 'choose a *.mat file');
        StereoInfo.cameraParams = cameraParamsFormatConvertFromMMC(CalibrationFilepath, CalibrationFile);
        clear CalibrationFile CalibrationFilepath
    case 3
        % Import Calib. results from DICe Calibrator
        [CalibrationFile, CalibrationFilepath]  = uigetfile({'*.xml'}, 'choose a *.xml file');
        StereoInfo.cameraParams = cameraParamsFormatConvertFromDICe(CalibrationFilepath, CalibrationFile);
        clear CalibrationFile CalibrationFilepath
    case 4
        % Import Calib. results with the same format as OpenCorr
        [CalibrationFile, CalibrationFilepath]  = uigetfile({'*.csv'}, 'choose a *.csv file');
        StereoInfo.cameraParams = cameraParamsFormatConvertFromOpenCorrFormat(CalibrationFilepath, CalibrationFile);
        clear CalibrationFile CalibrationFilepath
end
%% Section 3.2 Stereo Matching
stereoMatchShapeOrder = 1; % Currently, we only support 1st shape function

[StereoInfo, RD_L, RD_R] = StereoMatch_STAQ(RD_L,RD_R,imgNormalized_L{1},imgNormalized_R{1},...
    fileNameLeft,maskLeft{1},maskRight{1} ,DICpara,StereoInfo,stereoMatchShapeOrder);

% Optional: sanity-check stereo calibration on frame 1. Extract to its
% own helper (verify_stereo_calibration.m) in the C3.5 workstream.
if isfield(DICpara, 'verifyStereoReconstruction') && DICpara.verifyStereoReconstruction
    verify_stereo_calibration(StereoInfo);
end

%% Section 4 Temporal Matching
% Only 1st-order shape function is implemented; 2nd-order path not supported.
shapeFuncOrder = 1;

% Left image series
RD_L = TemporalMatch_quadtree_ST1(DICpara, fileNameLeft,maskLeft,imgNormalized_L,RD_L,StereoInfo, 'camera0',shapeFuncOrder);

%% Right image series
RD_R = TemporalMatch_quadtree_ST1(DICpara,fileNameRight,maskRight,imgNormalized_R,RD_R,StereoInfo, 'notCamera0',shapeFuncOrder);

%% Section 5 3D-Result
% Obtain the matched 2D point pairs from the left and right images
matchedPairs = organizeMatchedPairds_quadtree_ST1(RD_L,RD_R,DICpara);
% Calculate the 3D coordinates
[FinalResult,reprojectionErrors] = stereoReconstruction_quadtree(matchedPairs,StereoInfo.cameraParams);
                                                                                      
% Optional:Check the 3D reconstruction results
check3DReconstructionResults(reprojectionErrors, FinalResult, RD_L,13);

%% Section 6: Compute strains/ Plot disp. and strains
close all;
fprintf('------------ Section 6 Start ------------ \n')
% -------------------------------------------------------------------
% This section computes strain fields and plots disp and strain results.
% -------------------------------------------------------------------
% ------ Convert units from pixels to the physical world ------
DICpara.um2px = 1;

% Image save path
DICpara.outputFilePath = [];

% ------ Smooth displacements ------
DICpara.DoYouWantToSmoothOnceMore = 0;

% ------ Strain type (only Green-Lagrangian is supported) ------
DICpara.StrainType = 0;

% ------ Choose image to plot (first only, second and next images) ------
DICpara.Image2PlotResults = funParaInput('Image2PlotResults');

% ------ Save fig format ------
DICpara.MethodToSaveFig = funParaInput('SaveFigFormat');

% ------ Choose overlay image transparency ------
DICpara.OrigDICImgTransparency = 1;
if DICpara.MethodToSaveFig == 1
    DICpara.OrigDICImgTransparency = 0.8;
end

%--------- Transform the disp. to other coor. sys.?---------
RotationMatrix = [1 0 0; 0 1 0; 0 0 1]; TranslationMatrix = [0 0 0];
DICpara.transformDisp = funParaInput('TransformDispOrNot');
if DICpara.transformDisp == 1
    Base_Points2D = getBasePoints(imageLeft{1,1}',maskLeft{1,1}');
    [RotationMatrix,TranslationMatrix] = GetRTMatrix( RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM, Base_Points2D, FinalResult.Coordinates(1,:));
    DICpara.RforStrainCal = RotationMatrix;
    close all;
end


% ------ Strain window size ------
prompt = 'What is your strain size? e.g. 3,5,7...\nInput: ';
strain_size = input(prompt);
strain_length = (strain_size-1) * DICpara.winstepsize + 1;
VSG_length = strain_length + DICpara.winsize;
fprintf('Your strain size is %d * %d (Unit: Calculated Points) \nYour VSG size is %d * %d (Unit: Pixel) \n', strain_size,strain_size,VSG_length,VSG_length);


%% Initialization
coefficients = cell(3,1);
FinalResult.ResultStrainWorld{1} = 0;
FinalResult.Displacement_smooth(1,:) = FinalResult.Displacement(1,:); % All zeros

for ImgSeqNum = 2: length(imgNormalized_L)
    fprintf('  Frame %d/%d\n', ImgSeqNum, length(imgNormalized_L));
    close all;
    ImageName = fileNameLeft{1,ImgSeqNum};
    FullImageName_current = fullfile(fileNameLeft{2,ImgSeqNum}, fileNameLeft{1,ImgSeqNum});
    FullImageName_first   = fullfile(fileNameLeft{2,1}, fileNameLeft{1,1});

    % ------ Smooth displacements (post-ALDIC, 3D) ------
    % DICpara.Smooth3DTimes: number of outer-loop passes of
    % funSmoothDisp_Quadtree. 0 disables smoothing. Default is 3.
    if ~isfield(DICpara,'Smooth3DTimes') || isempty(DICpara.Smooth3DTimes)
        DICpara.Smooth3DTimes = 3;
    end
    try
        for SmoothTimes = 1:DICpara.Smooth3DTimes
            FinalResult.Displacement_smooth(ImgSeqNum,:) = funSmoothDisp_Quadtree( ...
                FinalResult.Displacement(ImgSeqNum,:), ...
                RD_L.ResultFEMeshEachFrame{1}, DICpara);
        end
    catch
    end
    if DICpara.Smooth3DTimes == 0
        % Not smoothed: mirror raw into the _smooth slot so downstream
        % consumers (PlaneFit3, plot) see a valid field.
        FinalResult.Displacement_smooth(ImgSeqNum,:) = FinalResult.Displacement(ImgSeqNum,:);
    end


    FinalResult = ConvertCoorAndDisp(FinalResult,RotationMatrix,TranslationMatrix);

    % ----- Compute strain field ------
    %  coefficients = [U,x V,x 0; U,y V,y 0; U,z V,z 0]
    [DICpara,coefficients,voidIndex] = PlaneFit3_Quadtree( strain_size, DICpara , RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM, FinalResult.Coordinates(1,:) ,...
        FinalResult.Coordinates(ImgSeqNum,:), FinalResult.Displacement_smooth(ImgSeqNum,:),...
        FullImageName_first, ImgSeqNum, 'Local'); % 'Local' or 'Specific Coordinate'

    % ---------------------------------------------------------
    % Plotting Configuration - Easy Mode Switching
    % ---------------------------------------------------------

    % Define what to plot
    DICpara.plots_disp_to_generate = {'u', 'v', 'w', 'magnitude'};
    DICpara.plots_strain_to_generate = {'exx', 'eyy', 'exy', 'e1', 'e2'};
    % === QUICK MODE SELECTION ===
    % Choose one of: 'auto', 'custom', 'interactive'
    plotting_mode = 'auto';  % Change this line to switch modes instantly

    DICpara = setPlottingParameters(DICpara, plotting_mode);

    % -----------------------------------------------------

    % % ----- Plot disp and strain ------
    if DICpara.DICIncOrNot == 0
        % Displacement
        if DICpara.transformDisp == 0
            PlotdispQuadtreeMasks3D_acc_ST1(FinalResult.DisplacementNew(ImgSeqNum,:),RD_L.ResultDisp{ImgSeqNum-1,1}.U,...
                RD_L.ResultFEMeshEachFrame{1,1}, FullImageName_first,FullImageName_current,DICpara,voidIndex);
        else
            PlotdispQuadtreeMasks3D_acc_ST1(FinalResult.Displacement(ImgSeqNum,:),RD_L.ResultDisp{ImgSeqNum-1,1}.U,...
                RD_L.ResultFEMeshEachFrame{1,1}, FullImageName_first,FullImageName_current,DICpara,voidIndex);
        end
        % Strain
        [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
            strain_maxshear,strain_vonMises,dwdx,dwdy] = PlotstrainQuadtreeMasks3D_acc_ST1(RD_L.ResultDisp{ImgSeqNum-1,1}.U,coefficients,voidIndex, ...
            RD_L.ResultFEMeshEachFrame{1,1},FullImageName_first,FullImageName_current, DICpara);
    elseif DICpara.DICIncOrNot == 1
        % Displacement
        if DICpara.transformDisp == 0
            PlotdispQuadtreeMasks3D_inc_ST1(FinalResult.DisplacementNew(ImgSeqNum,:),RD_L.ResultDisp_inc{ImgSeqNum-1,1}.U,RD_L.ResultDisp{ImgSeqNum-1,1}.U,RD_L.ResultFEMeshEachFrame{ImgSeqNum-1,1},...
                RD_L.ResultFEMeshEachFrame{1,1}, FullImageName_first,FullImageName_current, maskLeft{ ImgSeqNum },DICpara,voidIndex);
        else
            PlotdispQuadtreeMasks3D_inc_ST1(FinalResult.Displacement(ImgSeqNum,:),RD_L.ResultDisp_inc{ImgSeqNum-1,1}.U,RD_L.ResultDisp{ImgSeqNum-1,1}.U,RD_L.ResultFEMeshEachFrame{ImgSeqNum-1,1},...
                RD_L.ResultFEMeshEachFrame{1,1}, FullImageName_first,FullImageName_current, maskLeft{ ImgSeqNum },DICpara,voidIndex);
        end
        % Strain
        [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
            strain_maxshear,strain_vonMises] = PlotstrainQuadtreeMasks3D_inc_ST1(RD_L.ResultDisp_inc{ImgSeqNum-1,1}.U,RD_L.ResultDisp{ImgSeqNum-1,1}.U,coefficients,voidIndex, ...
            RD_L.ResultFEMeshEachFrame{1,1},RD_L.ResultFEMeshEachFrame{ImgSeqNum-1,1},FullImageName_first,FullImageName_current, ...
            maskLeft{ ImgSeqNum },DICpara);
    end

    % ----- Save strain results ------
    FinalResult.ResultStrainWorld{ImgSeqNum,1}= struct('strain_exx',strain_exx,...
        'strain_eyy',strain_eyy,'strain_exy',strain_exy,'strain_principal_max',strain_principal_max,'strain_principal_min',strain_principal_min, ...
        'strain_maxshear',strain_maxshear,'strain_vonMises',strain_vonMises);
    clear strain_exx  strain_eyy  strain_ezz  strain_exy  strain_eyz  strain_exz strain_maxshear strain_principal_max strain_principal_min strain_vonMises
    close all;

end

fprintf('------------ Section 6 Done ------------ \n \n')
