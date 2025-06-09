% ---------------------------------------------
% Augmented Lagrangian Digital Image Correlation (AL-DIC)
% working with an adaptive quadtree mesh
%
% Author: Jin Yang, PhD @Caltech
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Date: 2015.04,06,07; 2016.03,04; 2020.11;2021.12
% ---------------------------------------------

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
%--------   Notes --------------
% Strategy 1 needs Left and Right masks
% Strategy 2 needs Left masks only
% Inc mode needs all updated masks
% Acc mode needs only first mask
% NOTE: the determination of min. num. of masks to be developed
%-------------------------------
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images and masks ======
[fileNameLeft, fileNameRight, imageLeft,imageRight, LoadImgMethod] = ReadImage3DStereo_STAQ; % Load DIC raw images
[~, ~, maskLeft, maskRight, ~] = ReadImage3DStereo_STAQ_mask(imageLeft{1},imageRight{1}); % Load DIC masks

% ====== Set up DIC paras =====
DICpara = setDICParas_STAQ(fileNameLeft,imageLeft,maskLeft,maskRight,LoadImgMethod);
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
% calib_method = funParaInput('CalibrationMethod');
calib_method = 3;
switch calib_method
    case 0
        % Calibration in MATLAB
        StereoInfo = StereoCameraCalibration;
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

%%%%%%%%%%%%%%%%%%%%% Test 3D construction: START %%%%%%%%%%%%%%%%%%%%%
RD0_L_Pts = StereoInfo.ResultFEMeshEachFrame.coordinatesFEM;
RD0_R_Pts = StereoInfo.ResultFEMesh_corr;
figure, plot(RD0_L_Pts(:,1),RD0_L_Pts(:,2),'o')
hold on; plot(RD0_R_Pts(:,1),RD0_R_Pts(:,2),'o')

matchedPairs{1,1} = [RD0_L_Pts, RD0_R_Pts];
    
% Try to calculate the reprojection errors using the first two frames to reconstruct 3D coordinates
% [FinalResult.Coordinates,reprojectionErrors] = stereoReconstruction(StereoInfo.cameraParams, RD0_L_Pts, RD0_R_Pts);

cameraParams = StereoInfo.cameraParams;
%cameraParams = StereoInfo.cameraParams_optimized;

K_left = cameraParams.cameraParamsLeft.K;
K_right = cameraParams.cameraParamsRight.K;
R_left = [1 0 0; 0 1 0; 0 0 1]; % Set left Camera Coordinate as World Coordinate
T_left = [0 0 0]'; % Set left Camera Coordinate as World Coordinate
R_right = cameraParams.rotationMatrix;
T_right = cameraParams.translationVector';

% Undistort points before doing 3D reconstruction
[matchedPairs_undistort]= funUndistortPoints(matchedPairs,cameraParams);

P_left = K_left * [R_left, T_left];
P_right = K_right * [R_right, T_right];
reconstructedPoints = cell(size(matchedPairs_undistort,1),1);
reprojectionErrors = cell(size(matchedPairs_undistort,1),1);

for i = 1:size(matchedPairs_undistort,1)
    [reconstructedPoints{i,1},reprojectionErrors{i}]= triangulate(matchedPairs_undistort{i,1}(:, 1:2), matchedPairs_undistort{i,1}(:, 3:4), P_left, P_right);
end
figure; scatter3(reconstructedPoints{1,1}(:,1),reconstructedPoints{1,1}(:,2),reconstructedPoints{1,1}(:,3));
disp(['mean_repo = ',num2str(mean(reprojectionErrors{1}))]);
% [reconstructedPoints,reprojectionErrors] = triangulatePoints(matchedPairs_undistort, K_left, K_right, R_left, T_left, R_right, T_right);
%%%%%%%%%%%%%%%%%%%%% Test end %%%%%%%%%%%%%%%%%%%%%

%% Section 4 Temporal Matching
% Left image series
shapeFuncOrder = 1; % Curre1ntly, we only support 1st shape function
RD_L = TemporalMatch_quadtree_ST1(DICpara, fileNameLeft,maskLeft,imgNormalized_L,RD_L,StereoInfo, 'camera0',shapeFuncOrder);

%% Right image series
shapeFuncOrder = 1; % Curre1ntly, we only support 1st shape function
RD_R = TemporalMatch_quadtree_ST1(DICpara,fileNameRight,maskRight,imgNormalized_R,RD_R,StereoInfo, 'notCamera0',shapeFuncOrder);

%% Section 5 3D-Result
% Obtain the matched 2D point pairs from the left and right images
matchedPairs = organizeMatchedPairds_quadtree_ST1(RD_L,RD_R,DICpara);
% Calculate the 3D coordinates
[FinalResult,reprojectionErrors] = stereoReconstruction_quadtree(matchedPairs,StereoInfo.cameraParams);

%------------------ Check the 3D reconstruction results---------------
close all
disp(['3D_reconstruction_error_mean = ', num2str(mean(reprojectionErrors{2,1}(:)))]);
disp(['3D_reconstruction_error_std = ', num2str(std(reprojectionErrors{2,1}(:)))]);

figure, subplot(3,2,1);
[h1] = histogram(reprojectionErrors{1,1}, 100); title('reprojectError Frame 1')
subplot(3,2,2);
[h2] = histogram(reprojectionErrors{2,1}, 100); title('reprojectError Frame 2')

% For debug
subplot(3,2,3);
scatter3(FinalResult.Coordinates{1,1},FinalResult.Coordinates{1,2},FinalResult.Coordinates{1,3});
title('Coordiantes at Frame 1')

subplot(3,2,4);
scatter3(FinalResult.Coordinates{2,1},FinalResult.Coordinates{2,2},FinalResult.Coordinates{2,3});
title('Coordiantes at Frame 2')

subplot(3,2,5);
quiver3(FinalResult.Coordinates{1,1},FinalResult.Coordinates{1,2},FinalResult.Coordinates{1,3},...
    FinalResult.Displacement{2,1},FinalResult.Displacement{2,2},FinalResult.Displacement{2,2})
title('Tracked disp vector')

subplot(3,2,6);
scatter3(FinalResult.Coordinates{2,1}-FinalResult.Coordinates{1,1},...
    FinalResult.Coordinates{2,2}-FinalResult.Coordinates{1,2},...
    FinalResult.Coordinates{2,3}-FinalResult.Coordinates{1,3});
title('Tracked disp components')

% Display and check the first deformed image results
checkNum = 2;
U3Dstereo = [FinalResult.Displacement{checkNum,1}, FinalResult.Displacement{checkNum,2}, FinalResult.Displacement{checkNum,3}]';
Magnitude = sqrt(U3Dstereo(1,:).^2 + U3Dstereo(2,:).^2 + U3Dstereo(3,:).^2);

Plotdisp_show3D_STAQ(U3Dstereo(:),RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM,...
    RD_L.ResultFEMeshEachFrame{1,1}.elementsFEM, 'NoEdgeColor');

%% Section 6: Compute strains/ Plot disp. and strains
close all;
fprintf('------------ Section 8 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain fields and plot disp and strain results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Convert units from pixels to the physical world ------
% DICpara.um2px = funParaInput('ConvertUnit');
DICpara.um2px = 1;

% Image save path
DICpara.outputFilePath = [];

% ------ Smooth displacements ------
%DICpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot');
DICpara.DoYouWantToSmoothOnceMore = 0;

% ------ Choose strain computation method ------
%DICpara.MethodToComputeStrain = funParaInput('StrainMethodOp');

% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
%DICpara.StrainType = funParaInput('StrainType');
DICpara.StrainType = 0; % Currently, only support Green-Lagrangian strain

% ------ Choose image to plot (first only, second and next images) ------
DICpara.Image2PlotResults = funParaInput('Image2PlotResults');

% ------ Save fig format ------
DICpara.MethodToSaveFig = funParaInput('SaveFigFormat');
%DICpara.MethodToSaveFig = 1;

% ------ Choose overlay image transparency ------
DICpara.OrigDICImgTransparency = 1;
if DICpara.MethodToSaveFig == 1
    % DICpara.OrigDICImgTransparency = funParaInput('OrigDICImgTransparency');
    DICpara.OrigDICImgTransparency = 0.5;
end

%--------- Transform the disp. to other coor. sys.?---------
RotationMatrix = [1 0 0; 0 1 0; 0 0 1]; TranslationMatrix = [0 0 0];
if isempty(DICpara.transformDisp) || DICpara.transformDisp == 1
    DICpara.transformDisp = funParaInput('TransformDispOrNot');
    if DICpara.transformDisp == 0
        Base_Points2D = getBasePoints(imageLeft{1,1});
        [RotationMatrix,TranslationMatrix] = GetRTMatrix( RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM, Base_Points2D, FinalResult.Coordinates(1,:));
        DICpara.RforStrainCal = RotationMatrix;
        close all;
    end
end

% ------ Start main part ------
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
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(imgNormalized_L))]);
    close all;
    ImageName = fileNameLeft{1,ImgSeqNum};
    FullImageName_current = [fileNameLeft{2,ImgSeqNum} ,'\', fileNameLeft{1,ImgSeqNum}];
    FullImageName_first = [fileNameLeft{2,1} ,'\', fileNameLeft{1,1}];

    % x0World = DICpara.um2px*x0;
    % y0World = DICpara.um2px*y0; % Ignore this: (size(ImgNormalized_L{1},2)+1-y0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt);
    % ------ Smooth displacements ------

    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            % if DICpara.transformDisp == 0
            %     FinalResult.Displacement_smooth(ImgSeqNum,:) = funSmoothDisp_Quadtree(FinalResult.DisplacementNew(ImgSeqNum,:),RD_L.ResultFEMeshEachFrame{1},DICpara);
            % else
                FinalResult.Displacement_smooth(ImgSeqNum,:) = funSmoothDisp_Quadtree(FinalResult.Displacement(ImgSeqNum,:),RD_L.ResultFEMeshEachFrame{1},DICpara);
            % end
            %%DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end


    FinalResult = ConvertCoorAndDisp(FinalResult,RotationMatrix,TranslationMatrix);



    % ----- Compute strain field ------
    %  coefficients = [U,x V,x 0; U,y V,y 0; U,z V,z 0]
    [DICpara,coefficients,voidIndex] = PlaneFit3_Quadtree( strain_size, DICpara , RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM, FinalResult.Coordinates(1,:) ,...
        FinalResult.Coordinates(ImgSeqNum,:), FinalResult.Displacement_smooth(ImgSeqNum,:),...
        FullImageName_first, ImgSeqNum, 'Local'); % 'Local' or 'Specific Coordinate'

    % Pre-display displacement results
    % close all; Plotdisp_show_3D(FinalResult.Displacement(ImgSeqNum,:),RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM,...
    %     RD_L.ResultFEMeshEachFrame{1,1}.elementsFEM,DICpara,'NoEdgeColor');

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


%%%---------------------Plot 3D scatter figure------------------------------
%     close all;
%     XXX = FinalResult.CoordinatesNew{ImgSeqNum,1}(~voidIndex);YYY = FinalResult.CoordinatesNew{ImgSeqNum,2}(~voidIndex);ZZZ= FinalResult.CoordinatesNew{ImgSeqNum,3}(~voidIndex);
%     temp{1} = FinalResult.DisplacementNew{ImgSeqNum,1}(~voidIndex);
%     temp{2} = FinalResult.DisplacementNew{ImgSeqNum,2}(~voidIndex);
%     temp{3} = FinalResult.DisplacementNew{ImgSeqNum,3}(~voidIndex);
%     temp{4} = strain_exx(~voidIndex);    
%     temp{5} = strain_exy(~voidIndex);    
%     temp{6} = strain_eyy(~voidIndex);  
%     temp{7} = dwdx(~voidIndex);  
%     temp{8} = dwdy(~voidIndex); 
%     for iii = 1:8%1:3
%     figure('Position',[500,500,550,350]);
%     scatter3(XXX,YYY,ZZZ,20,temp{iii},'filled') ;
%     ax = gca; set(ax, 'View', [-51.0493,77.4470],'FontSize', 22,'FontName', 'Arial','LineWidth', 0.5,'GridColor', 'k',  'GridAlpha', 0.5); 
%     axis tight;
%     xlabel('X [mm]', 'FontSize', 18, 'FontName', 'Arial','Rotation',45); 
%     ylabel('Y [mm]', 'FontSize', 18, 'FontName', 'Arial','Rotation',-35); 
%     zlabel('Z [mm]', 'FontSize', 18, 'FontName', 'Arial','Rotation',90); 
%     figure(iii);
%     colormap("turbo");
%     colorbar; 
%     if iii == 3
%         caxis([0,1.2]);
%     elseif iii == 2
%         caxis([-0.1,0.05]);
%     elseif iii == 1
%         caxis([-0.08 0.03])
%     elseif iii ==4
% caxis([-0.04 0.025]) %bulge
%     elseif iii ==5
% caxis([-0.017 0.017]) % bulge
% 
%     elseif iii ==6 
% caxis([-0.02 0.025]) % bulge
% 
%     elseif iii ==7 
%  caxis([-0.13 0.12])  % bulge
% 
%     elseif iii ==8 
% caxis([-0.15 0.15]) %bulge
% 
%     end
%     zlim([-0.9 1.6])
%     [~,imgname,imgext] = fileparts([fileNameLeft{2,ImgSeqNum},'\',fileNameLeft{1,ImgSeqNum}]);
%     print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_',num2str(iii)],'-djpeg','-r0');
%     end
%     XY2D = RD_L.ResultFEMeshEachFrame{1, 1}.coordinatesFEM(~voidIndex,:);
%     imagetemp = im2uint8(im2double(imageLeft{1})/ 65535);
%     for iiiii = 1:size(XY2D,1)
%         temp{9}(iiiii,1) = imagetemp(XY2D(iiiii,1),XY2D(iiiii,2));
%     end
%     figure(9); scatter3(XXX,YYY,ZZZ,40,temp{9},'filled') ; colormap gray
%     ax = gca; set(ax, 'View', [-51.0493,77.4470],'FontSize', 28,'FontName', 'Arial','LineWidth', 0.5,'GridColor', 'k',  'GridAlpha', 0.5); 
%     print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_',num2str(9)],'-djpeg','-r0');

%%%--------------------------------------------------------------

    % debug
    %strain_exx = strain_exx(50:end,:);
    % disp(['Exx_mean = ', num2str(mean(strain_exx(:)))]);
    % disp(['Exx_std = ', num2str(std(strain_exx(:)))]);
    % %strain_eyy = strain_eyy(50:end,:);
    % disp(['Eyy_mean = ', num2str(mean(strain_eyy(:)))]);
    % disp(['Eyy_std = ', num2str(std(strain_eyy(:)))]);
    % close all

    % ----- Save strain results ------
    % FinalResult.ResultStrainWorld{ImgSeqNum,1}= struct('strain_exx',strain_exx,...
    %     'strain_eyy',strain_eyy,'strain_exy',strain_exy);
    % 'strain_principal_max',strain_principal_max,'strain_principal_min',strain_principal_min, ...
    %     'strain_maxshear',strain_maxshear,'strain_vonMises',strain_vonMises
    clear strain_exx  strain_eyy  strain_ezz  strain_exy  strain_eyz  strain_exz strain_maxshear strain_principal_max strain_principal_min strain_vonMises
    % ------ Save figures for tracked displacement and strain fields ------
    % SaveFigFilesDispAndStrain;
    % SaveFigFilesDispAndStrain_v2; % For Ch 1.0 Sample 3
    SaveFigFilesDispAndStrain_v3; % For pig heart
    close all;

end

% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized_L)} ------
fprintf('------------ Section 8 Done ------------ \n \n')

% ------ Save data again including solved strain fields ------
% selectedFolder = uigetdir('', 'Select save folder');
% results_name = [selectedFolder,'\','results_',ImageName,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
% save(results_name, 'FullImageName_current','DICpara','FinalResult');


%% Zach: Section 6.5: Plot quadtree mesh

shapeFuncOrder = 1;
for ImgSeqNum = [2]%[11 16 21 26]
    TemporalMatch_STAQ_GenerateQuadtreeMeshOnly(ImgSeqNum,DICpara, fileNameLeft,maskLeft,imgNormalized_L,RD_L,StereoInfo, 'camera0',shapeFuncOrder);
end


%% Section 7: Compute stress
fprintf('------------ Section 9 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute stress fields and plot stress fields
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Choose material model ------
DICpara.MaterialModel = funParaInput('MaterialModel');
% ------ Define parameters in material models ------
if (DICpara.MaterialModel == 1) || (DICpara.MaterialModel == 2) % Linear elasticity
    fprintf('Define Linear elasticity parameters \n')
    fprintf("Young's modulus (unit: Pa): \n"); prompt = 'Input here (e.g., 69e9): ';
    DICpara.MaterialModelPara.YoungsModulus = input(prompt);
    fprintf("Poisson's ratio: \n"); prompt = 'Input here (e.g., 0.3): ';
    DICpara.MaterialModelPara.PoissonsRatio = input(prompt);
    fprintf('------------------------------------- \n');
end

% ------ Start main part ------
for ImgSeqNum = 2 : length(imgNormalized_L)

    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(imgNormalized_L))]); close all;

    % ------ Plot stress ------
    if DICpara.OrigDICImgTransparency == 1
        [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
            stress_principal_min_xyplane, stress_maxshear_xyplane, ...
            stress_maxshear_xyz3d, stress_vonMises]  =  Plotstress0( ...
            DICpara,RD_L.ResultStrainWorld{ImgSeqNum-1},size(imgNormalized_L{1}));

    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name_L{1,1}" corresponds to the first image
            [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = Plotstress( ...
                DICpara,RD_L.ResultStrainWorld{ImgSeqNum-1},size(imgNormalized_L{1}),fileNameLeft{1,1});

        else % Plot over second or next deformed images
            [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = Plotstress( ...
                DICpara,RD_L.ResultStrainWorld{ImgSeqNum-1},size(imgNormalized_L{1}),fileNameLeft{1,ImgSeqNum});

        end
    end


    % ------ Save figures for computed stress fields ------
    SaveFigFilesStress;

    % ----- Save strain results ------
    RD_L.ResultStressWorld{ImgSeqNum-1} = struct('stressxCoord',RD_L.ResultStrainWorld{ImgSeqNum-1}.strainxCoord,'stressyCoord',RD_L.ResultStrainWorld{ImgSeqNum-1}.strainyCoord, ...
        'stress_sxx',stress_sxx,'stress_sxy',stress_sxy,'stress_syy',stress_syy, ...
        'stress_principal_max_xyplane',stress_principal_max_xyplane, 'stress_principal_min_xyplane',stress_principal_min_xyplane, ...
        'stress_maxshear_xyplane',stress_maxshear_xyplane,'stress_maxshear_xyz3d',stress_maxshear_xyz3d, ...
        'stress_vonMises',stress_vonMises);

end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized_L)} ------
fprintf('------------ Section 9 Done ------------ \n \n')

% ------ Save data again including solved stress fields ------
results_name = ['results_',ImageName,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'fileNameLeft','DICpara','DICmesh','RD_L.ResultDisp','RD_L.ResultDefGrad','RD_L.ResultFEMesh','RD_L.ResultFEMeshEachFrame',...
    'ALSub1Time','ALSub2Time','ALSolveStep','RD_L.ResultStrainWorld','RD_L.ResultStressWorld');




