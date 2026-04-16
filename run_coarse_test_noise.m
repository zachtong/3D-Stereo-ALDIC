%% Non-Interactive Coarse Test for 3D-Stereo-ALDIC
% Runs the full DIC pipeline with large step size for quick validation.
% Stops after displacement field (no strain computation).
%
% Usage: Run this script from the 3D-Stereo-ALDIC root directory.

close all; clear; clc; clearvars -global
fprintf('========================================\n');
fprintf('  3D-Stereo-ALDIC: Noise Analysis Run\n');
fprintf('========================================\n\n');

%% =============== USER CONFIGURATION ===============
% Data paths (relative to project root)
data_base = './Stereo_DIC_Challenge_2.1_Bespoke/000_zach_modified_Bespoke_corrected';

img_left_folder  = fullfile(data_base, 'Noise', 'left');
img_right_folder = fullfile(data_base, 'Noise', 'right');
mask_left_folder  = fullfile(data_base, 'Noise_masks', 'left');
mask_right_folder = fullfile(data_base, 'Noise_masks', 'right');
calib_xml_path = fullfile(data_base, 'Calibration_DICe.xml');

% DIC parameters for coarse test
WINSIZE      = 24;   % Subset size (pixels) — MUST be even (code uses winsize/2 for integer indexing)
WINSTEPSIZE  = 8;   % Step size (pixels, power of 2)
WINSIZE_MIN  = 8;   % Same as WINSTEPSIZE → uniform grid (no quadtree refinement)
SEARCH_DIST  = 40;  % FFT search distance (pixels)

% ===================================================

%% Section 1: Environment setup
fprintf('------------ Section 1 Start ------------ \n')
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
try mex -O ba_interp2_spline.cpp; warning('off'); catch; end

addpath('./examples','./func','./func_quadtree/rbfinterp/', ...
    './plotFiles/','./func_quadtree','./func_quadtree/refinement', ...
    './plotFiles/export_fig-d966721/');
fprintf('------------ Section 1 Done ------------ \n\n')

%% Section 2: Load images, masks, and set parameters
fprintf('------------ Section 2 Start ------------ \n')

% ====== Load images (replaces ReadImage3DStereo_STAQ) ======
fprintf('Loading images from:\n  L: %s\n  R: %s\n', img_left_folder, img_right_folder);

img_L_files = dir(fullfile(img_left_folder, '*.tif'));
img_R_files = dir(fullfile(img_right_folder, '*.tif'));

numImages = length(img_L_files);
fprintf('Found %d LEFT and %d RIGHT images.\n', numImages, length(img_R_files));
assert(numImages == length(img_R_files), 'Mismatch in left/right image count!');

% Build file_name cell arrays (same format as ReadImage3DStereo_STAQ output)
fileNameLeft  = cell(6, numImages);
fileNameRight = cell(6, numImages);
imageLeft  = cell(numImages, 1);
imageRight = cell(numImages, 1);

for i = 1:numImages
    fileNameLeft{1,i}  = img_L_files(i).name;
    fileNameLeft{2,i}  = img_L_files(i).folder;
    fileNameRight{1,i} = img_R_files(i).name;
    fileNameRight{2,i} = img_R_files(i).folder;

    img = imread(fullfile(img_L_files(i).folder, img_L_files(i).name));
    if size(img,3) > 1, img = rgb2gray(img); end
    imageLeft{i} = double(img)';

    img = imread(fullfile(img_R_files(i).folder, img_R_files(i).name));
    if size(img,3) > 1, img = rgb2gray(img); end
    imageRight{i} = double(img)';
end
fprintf('All images loaded and transposed.\n');

% ====== Set DIC parameters (replaces setDICParas_IncOrNot + setDICParas_STAQ) ======
DICpara.DICIncOrNot = 1;  % 1 = incremental mode
DICpara.ImgSeqIncUnit = 1;
DICpara.ImgSeqIncROIUpdateOrNot = 0;
DICpara.NewFFTSearch = 1;

DICpara.winsize = WINSIZE;
DICpara.winstepsize = WINSTEPSIZE;
DICpara.winsizeMin = WINSIZE_MIN;
DICpara.Subpb2FDOrFEM = 1;  % Finite difference
DICpara.ClusterNo = 0;      % Auto-detect in TemporalMatch
DICpara.ImgSize = size(imageLeft{1});
DICpara.transformDisp = [];
DICpara.LoadImgMethod = 1;  % Separate L/R folders

% ROI: full image for incremental mode (per-frame mask applied inside TemporalMatch)
DICpara.gridxyROIRange.gridx = [1, size(imageLeft{1}, 1)];
DICpara.gridxyROIRange.gridy = [1, size(imageLeft{1}, 2)];

% Non-interactive FFT search
DICpara.InitFFTSearchMethod = 1;       % Whole-field FFT search
DICpara.NewFFTSearchDistance = SEARCH_DIST;
DICpara.fixSearchDistanceOrNot = 0;    % 0 = fix distance for all frames

% Non-interactive: skip data-driven mode prompt
DICpara.dataDrivenOrNot = 1;  % 1 = do NOT use data-driven

% Plotting off for batch
DICpara.showImgOrNot = 0;

% ====== Load masks (replaces ReadImage3DStereo_STAQ_mask_inc) ======
fprintf('Loading masks from:\n  L: %s\n  R: %s\n', mask_left_folder, mask_right_folder);

mask_L_files = dir(fullfile(mask_left_folder, '*.tiff'));
mask_R_files = dir(fullfile(mask_right_folder, '*.tiff'));
% Also check .tif extension
if isempty(mask_L_files), mask_L_files = dir(fullfile(mask_left_folder, '*.tif')); end
if isempty(mask_R_files), mask_R_files = dir(fullfile(mask_right_folder, '*.tif')); end

fprintf('Found %d LEFT and %d RIGHT masks.\n', length(mask_L_files), length(mask_R_files));
assert(length(mask_L_files) == numImages, 'Mask count (%d) != image count (%d)!', length(mask_L_files), numImages);

maskLeft  = cell(numImages, 1);
maskRight = cell(numImages, 1);
for i = 1:numImages
    m = imread(fullfile(mask_L_files(i).folder, mask_L_files(i).name));
    if size(m,3) > 1, m = rgb2gray(m); end
    maskLeft{i} = logical(m)';

    m = imread(fullfile(mask_R_files(i).folder, mask_R_files(i).name));
    if size(m,3) > 1, m = rgb2gray(m); end
    maskRight{i} = logical(m)';
end
fprintf('All masks loaded.\n');

% ====== Normalize images ======
[imgNormalized_L, DICpara.gridxyROIRange] = funNormalizeImg(imageLeft, DICpara.gridxyROIRange);
[imgNormalized_R, DICpara.gridxyROIRange] = funNormalizeImg(imageRight, DICpara.gridxyROIRange);

% ====== Initialize result storage ======
nFrames = length(imgNormalized_L);

RD_L.ResultDisp = cell(nFrames-1, 1);
RD_L.ResultDefGrad = cell(nFrames-1, 1);
RD_L.ResultFEMeshEachFrame = cell(nFrames-1, 1);
RD_L.ResultFEMesh = cell(nFrames-1, 1);
RD_L.ResultDisp_inc = cell(nFrames-1, 1);

RD_R.ResultDisp = cell(nFrames-1, 1);
RD_R.ResultDefGrad = cell(nFrames-1, 1);
RD_R.ResultFEMeshEachFrame = cell(nFrames-1, 1);
RD_R.ResultFEMesh = cell(nFrames-1, 1);
RD_R.ResultDisp_inc = cell(nFrames-1, 1);

DICpara.ImgRefMask = double(maskLeft{1});

fprintf('------------ Section 2 Done ------------ \n\n')
fprintf('Config: winsize=%d, winstepsize=%d, winsizeMin=%d, searchDist=%d\n', ...
    WINSIZE, WINSTEPSIZE, WINSIZE_MIN, SEARCH_DIST);
fprintf('Images: %d frames, incremental mode\n\n', nFrames);

%% Section 3.1: Stereo Calibration (DICe)
fprintf('------------ Section 3.1 Start ------------ \n')
[calib_folder, calib_file, ~] = fileparts(calib_xml_path);
% cameraParamsFormatConvertFromDICe expects folder WITH trailing separator (like uigetfile)
if ~endsWith(calib_folder, filesep)
    calib_folder = [calib_folder, filesep];
end
StereoInfo.cameraParams = cameraParamsFormatConvertFromDICe(calib_folder, [calib_file, '.xml']);
fprintf('Calibration loaded from: %s\n', calib_xml_path);
fprintf('------------ Section 3.1 Done ------------ \n\n')

%% Section 3.2: Stereo Matching
fprintf('------------ Section 3.2 Start ------------ \n')
stereoMatchShapeOrder = 1;
[StereoInfo, RD_L, RD_R] = StereoMatch_STAQ(RD_L, RD_R, ...
    imgNormalized_L{1}, imgNormalized_R{1}, ...
    fileNameLeft, maskLeft{1}, maskRight{1}, ...
    DICpara, StereoInfo, stereoMatchShapeOrder);

% Quick 3D check with stereo match points
RD0_L_Pts = StereoInfo.ResultFEMeshEachFrame.coordinatesFEM;
RD0_R_Pts = StereoInfo.ResultFEMesh_corr;
matchedPairs{1,1} = [RD0_L_Pts, RD0_R_Pts];

cameraParams = StereoInfo.cameraParams;
K_left = cameraParams.cameraParamsLeft.K;
K_right = cameraParams.cameraParamsRight.K;
R_left = eye(3); T_left = [0 0 0]';
R_right = cameraParams.rotationMatrix;
T_right = cameraParams.translationVector';

[matchedPairs_undistort] = funUndistortPoints(matchedPairs, cameraParams);
P_left = K_left * [R_left, T_left];
P_right = K_right * [R_right, T_right];

[reconstructedPoints_check, reprojectionErrors_check] = triangulate( ...
    matchedPairs_undistort{1}(:,1:2), matchedPairs_undistort{1}(:,3:4), P_left, P_right);
fprintf('Stereo match reprojection error: mean=%.3f px\n', mean(reprojectionErrors_check));
fprintf('------------ Section 3.2 Done ------------ \n\n')

%% Section 4: Temporal Matching
fprintf('============================================\n');
fprintf('  Section 4: Temporal Matching\n');
fprintf('============================================\n\n');

% Left camera
fprintf('--- LEFT camera temporal matching ---\n');
tic;
RD_L = TemporalMatch_quadtree_ST1(DICpara, fileNameLeft, maskLeft, ...
    imgNormalized_L, RD_L, StereoInfo, 'camera0', 1);
t_left = toc;
fprintf('Left camera done in %.1f seconds (%.1f s/frame)\n\n', t_left, t_left/(nFrames-1));

% Right camera
fprintf('--- RIGHT camera temporal matching ---\n');
tic;
RD_R = TemporalMatch_quadtree_ST1(DICpara, fileNameRight, maskRight, ...
    imgNormalized_R, RD_R, StereoInfo, 'notCamera0', 1);
t_right = toc;
fprintf('Right camera done in %.1f seconds (%.1f s/frame)\n\n', t_right, t_right/(nFrames-1));

%% Section 5: 3D Reconstruction
fprintf('------------ Section 5 Start ------------ \n')
matchedPairs = organizeMatchedPairds_quadtree_ST1(RD_L, RD_R, DICpara);
[FinalResult, reprojectionErrors] = stereoReconstruction_quadtree(matchedPairs, StereoInfo.cameraParams);

fprintf('3D reconstruction done.\n');
fprintf('Mean reprojection errors per frame:\n');
for i = 1:length(reprojectionErrors)
    nFinite = nnz(isfinite(reprojectionErrors{i}));
    fprintf('  Frame %2d: mean=%.3f px (%d/%d pts)\n', i, ...
        mean(reprojectionErrors{i}, 'omitnan'), nFinite, numel(reprojectionErrors{i}));
end
fprintf('------------ Section 5 Done ------------ \n\n')

%% Quick Visualization: Displacement Field
fprintf('============================================\n');
fprintf('  Visualizing Displacement Fields\n');
fprintf('============================================\n\n');

coords2D = RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM;
maskRef = maskLeft{1};  % Reference frame mask (transposed like images)

% Build mask index: true = on specimen, false = inside hole or outside
onSpecimen = false(size(coords2D, 1), 1);
for k = 1:size(coords2D, 1)
    cx = round(coords2D(k, 1));
    cy = round(coords2D(k, 2));
    if cx >= 1 && cx <= size(maskRef, 1) && cy >= 1 && cy <= size(maskRef, 2)
        onSpecimen(k) = maskRef(cx, cy);
    end
end
fprintf('Mask filter: %d/%d quadtree nodes on specimen\n', nnz(onSpecimen), length(onSpecimen));

% Plot first, middle, and last frames
frames_to_plot = unique([2, round(nFrames/2), nFrames]);

for idx = 1:length(frames_to_plot)
    f = frames_to_plot(idx);
    if f < 2 || f > nFrames, continue; end

    try
        % FinalResult uses cell arrays: {frame, component}
        U = FinalResult.Displacement{f,1};
        V = FinalResult.Displacement{f,2};
        W = FinalResult.Displacement{f,3};
        X = FinalResult.Coordinates{f,1};
        Y = FinalResult.Coordinates{f,2};
        Z = FinalResult.Coordinates{f,3};

        mag = sqrt(U.^2 + V.^2 + W.^2);

        % Only plot points that are on specimen AND have finite values
        plotIdx = onSpecimen & isfinite(mag);

        figure('Name', sprintf('Frame %d: 3D Displacement', f), 'Position', [50+200*(idx-1), 100, 600, 500]);
        scatter3(X(plotIdx), Y(plotIdx), Z(plotIdx), 10, mag(plotIdx), 'filled');
        colorbar; colormap jet;
        title(sprintf('Frame %d/%d: Displacement Magnitude (px)', f, nFrames));
        xlabel('X'); ylabel('Y'); zlabel('Z');
        axis equal; view(2);
        drawnow;

        fprintf('Frame %d: %d plotted pts | |U|=%.1f, |V|=%.1f, |W|=%.1f (mean mag=%.1f px)\n', ...
            f, nnz(plotIdx), mean(abs(U(plotIdx))), mean(abs(V(plotIdx))), mean(abs(W(plotIdx))), mean(mag(plotIdx)));
    catch ME
        fprintf('Frame %d: visualization error - %s\n', f, ME.message);
    end
end

%% Save workspace (displacement only)
save_name = sprintf('results/noise/coarse_test_noise_ws%d_ss%d.mat', WINSIZE, WINSTEPSIZE);
save(save_name, 'DICpara', 'FinalResult', 'RD_L', 'RD_R', 'StereoInfo', ...
    'reprojectionErrors', 'matchedPairs', '-v7.3');
fprintf('\nResults saved to: %s\n', save_name);
fprintf('\n========================================\n');
fprintf('  Coarse test complete!\n');
fprintf('  Total time: LEFT=%.0fs + RIGHT=%.0fs\n', t_left, t_right);
fprintf('  Review displacement plots above.\n');
fprintf('  If results look reasonable, reduce winstepsize for final run.\n');
fprintf('========================================\n');

%% Section 7: Strain Computation (Green-Lagrange)
% Can be run independently after coarse test completes (Run Section).
% Uses PlaneFit3_Quadtree for KNN-based local plane fitting of deformation gradient,
% then computes Green-Lagrange strain E = 0.5*(F'*F - I).
fprintf('\n============================================\n');
fprintf('  Section 7: Strain Computation\n');
fprintf('============================================\n\n');

close all;
% --- Strain parameters ---
STRAIN_SIZE = 10;  % VSG size in calculated points (3=small/noisy, 5=balanced, 7=smooth)
strain_length = (STRAIN_SIZE - 1) * WINSTEPSIZE + 1;
VSG_length = strain_length + WINSIZE;
fprintf('Strain gauge: %dx%d pts = %dx%d px VSG (winstepsize=%d)\n', ...
    STRAIN_SIZE, STRAIN_SIZE, VSG_length, VSG_length, WINSTEPSIZE);

DICpara.um2px = 1;
DICpara.StrainType = 0;  % Green-Lagrangian

% --- Specimen coordinate system (O, X-direction, Y-direction in 2D pixel coords) ---
% O = origin, X = point defining +X axis, Y = point defining +Y axis
% Table values are 0-based (top-left = (0,0)), MATLAB is 1-based → add 1
Base_Points2D = [1024+1, 1224+1; 1509.69+1, 1219.38+1; 1020.10+1, 470.48+1];
fprintf('Base points (2D): O=[%.0f,%.0f], X=[%.1f,%.1f], Y=[%.1f,%.1f]\n', ...
    Base_Points2D(1,:), Base_Points2D(2,:), Base_Points2D(3,:));

coords2D_s = RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM;

% Compute rotation matrix R and translation T for specimen coordinate system
[RotationMatrix, TranslationMatrix] = GetRTMatrix( ...
    coords2D_s, Base_Points2D, FinalResult.Coordinates(1,:));
DICpara.RforStrainCal = RotationMatrix;
fprintf('Specimen coordinate system established (R, T computed)\n');

% Transform ALL frames to specimen coordinates
FinalResult = ConvertCoorAndDisp(FinalResult, RotationMatrix, TranslationMatrix);
fprintf('3D coordinates transformed to specimen coordinate system\n');

% --- Build mask filter for visualization ---
maskRef_s = maskLeft{1};
onSpecimen_s = false(size(coords2D_s, 1), 1);
for k = 1:size(coords2D_s, 1)
    cx = round(coords2D_s(k, 1));
    cy = round(coords2D_s(k, 2));
    if cx >= 1 && cx <= size(maskRef_s, 1) && cy >= 1 && cy <= size(maskRef_s, 2)
        onSpecimen_s(k) = maskRef_s(cx, cy);
    end
end
fprintf('Mask: %d/%d quadtree nodes on specimen\n\n', nnz(onSpecimen_s), length(onSpecimen_s));

% --- Displacement smoothing parameters ---
SMOOTH_SIGMA = WINSTEPSIZE;  % Gaussian sigma (px). 0=no smoothing, step/2=mild, step=moderate
fprintf('Displacement smoothing: sigma=%.0f px (0=off)\n', SMOOTH_SIGMA);

% Precompute KNN for smoothing (once, reused for all frames)
if SMOOTH_SIGMA > 0
    K_smooth = min(50, size(coords2D_s, 1));
    [sm_idx, sm_dist] = knnsearch(coords2D_s, coords2D_s, 'K', K_smooth, 'Distance', 'chebychev');
    sm_weights = exp(-sm_dist.^2 / (2 * SMOOTH_SIGMA^2));
    sm_weights(sm_dist > 3 * SMOOTH_SIGMA) = 0;  % Cutoff at 3*sigma
end

% --- Select frames for strain ---
frames_strain = 2:nFrames;  % ALL noise frames: noise strain IS required for noise baseline evaluation
FinalResult.ResultStrainWorld = cell(nFrames, 1);

for idx = 1:length(frames_strain)
    f = frames_strain(idx);
    if f < 2 || f > nFrames, continue; end
    fprintf('--- Frame %d/%d: computing strain ---\n', f, nFrames);

    % Smooth displacement before strain computation (camera coords)
    if SMOOTH_SIGMA > 0
        DispSmoothed = FinalResult.Displacement(f, :);
        for comp = 1:3
            raw = FinalResult.Displacement{f, comp};
            raw_nn = raw(sm_idx);                    % N×K neighbor values
            valid = isfinite(raw_nn) & (sm_dist <= 3*SMOOTH_SIGMA);
            w = sm_weights .* valid;
            w_sum = sum(w, 2);
            smoothed = sum(w .* raw_nn, 2) ./ w_sum;
            smoothed(w_sum == 0 | ~isfinite(raw)) = raw(w_sum == 0 | ~isfinite(raw));
            DispSmoothed{comp} = smoothed;
        end
    else
        DispSmoothed = FinalResult.Displacement(f, :);
    end

    % Compute deformation gradient coefficients in SPECIMEN coordinate system
    % PlaneFit3_Quadtree with 'Specific Coordinate' uses DICpara.RforStrainCal
    % to rotate reference coords and displacements before plane fitting
    % coefficients{i} = [u_x, v_x, w_x; u_y, v_y, w_y; 0, 0, 0]
    tic;
    [~, coefficients, voidIndex] = PlaneFit3_Quadtree(STRAIN_SIZE, DICpara, ...
        coords2D_s, ...
        FinalResult.Coordinates(1, :), ...   % Reference 3D coords (camera system)
        FinalResult.Coordinates(f, :), ...   % Current 3D coords (unused inside function)
        DispSmoothed, ...                    % Smoothed 3D displacement (camera system)
        '', f, 'Specific Coordinate');
    t_strain = toc;

    % Compute Green-Lagrange strain from deformation gradient
    nPts = length(coefficients);
    strain_exx = NaN(nPts, 1);
    strain_eyy = NaN(nPts, 1);
    strain_exy = NaN(nPts, 1);
    strain_e1  = NaN(nPts, 1);
    strain_e2  = NaN(nPts, 1);
    strain_maxshear  = NaN(nPts, 1);
    strain_vonMises  = NaN(nPts, 1);

    for i = 1:nPts
        if voidIndex(i), continue; end
        c = coefficients{i, 1};
        if any(isnan(c(:))), continue; end

        u_x = c(1,1); u_y = c(2,1);
        v_x = c(1,2); v_y = c(2,2);
        w_x = c(1,3); w_y = c(2,3);

        % Deformation gradient F = I + grad(u)
        F_mat = [1+u_x, u_y, 0; ...
                 v_x, 1+v_y, 0; ...
                 w_x, w_y,   1];

        % Green-Lagrange strain: E = 0.5*(F'*F - I)
        E = 0.5 * (F_mat' * F_mat - eye(3));

        strain_exx(i) = E(1,1);
        strain_eyy(i) = E(2,2);
        strain_exy(i) = E(1,2);

        % Principal strains (2D in-plane)
        avg = 0.5*(E(1,1) + E(2,2));
        R_principal = sqrt((0.5*(E(1,1) - E(2,2)))^2 + E(1,2)^2);
        strain_e1(i)  = avg + R_principal;
        strain_e2(i)  = avg - R_principal;
        strain_maxshear(i) = R_principal;

        % Von Mises strain
        strain_vonMises(i) = sqrt(E(1,1)^2 + E(2,2)^2 - E(1,1)*E(2,2) + 3*E(1,2)^2);
    end

    % Store results
    FinalResult.ResultStrainWorld{f} = struct( ...
        'strain_exx', strain_exx, 'strain_eyy', strain_eyy, 'strain_exy', strain_exy, ...
        'strain_principal_max', strain_e1, 'strain_principal_min', strain_e2, ...
        'strain_maxshear', strain_maxshear, 'strain_vonMises', strain_vonMises);

    % --- Statistics ---
    validIdx = onSpecimen_s & isfinite(strain_exx);
    nValid = nnz(validIdx);
    fprintf('  PlaneFit: %.1f s | Valid strain: %d/%d on-specimen pts\n', t_strain, nValid, nnz(onSpecimen_s));
    if nValid > 0
        fprintf('  exx: [%.4f, %.4f], eyy: [%.4f, %.4f], exy: [%.4f, %.4f]\n', ...
            min(strain_exx(validIdx)), max(strain_exx(validIdx)), ...
            min(strain_eyy(validIdx)), max(strain_eyy(validIdx)), ...
            min(strain_exy(validIdx)), max(strain_exy(validIdx)));
        fprintf('  Max shear: [%.4f, %.4f], Von Mises: [%.4f, %.4f]\n', ...
            min(strain_maxshear(validIdx)), max(strain_maxshear(validIdx)), ...
            min(strain_vonMises(validIdx)), max(strain_vonMises(validIdx)));
    end

    % === Common setup: reference config in specimen coords ===
    Xref = FinalResult.CoordinatesNew{1,1};
    Yref = FinalResult.CoordinatesNew{1,2};
    Zref = FinalResult.CoordinatesNew{1,3};
    vertices3D = [Xref, Yref, Zref];

    % Displacement in specimen coords
    Usp = FinalResult.DisplacementNew{f,1};
    Vsp = FinalResult.DisplacementNew{f,2};
    Wsp = FinalResult.DisplacementNew{f,3};

    % Valid node masks
    nodeValidDisp   = onSpecimen_s & isfinite(Usp) & isfinite(Xref);
    nodeValidStrain = onSpecimen_s & isfinite(strain_exx) & isfinite(Xref);

    % Filter elements: keep only those with ALL 4 corner nodes valid
    elemFEM = RD_L.ResultFEMeshEachFrame{1,1}.elementsFEM;
    faces = elemFEM(:, 1:4);

    % --- Combined figure: Displacement + Strain (2x5 layout) ---
    % Row 1: Ux, Uy, Uz, e1, e2
    % Row 2: |U|, exx, eyy, exy, max shear
    elemValidD = all(nodeValidDisp(faces), 2);
    faces_disp = faces(elemValidD, :);
    elemValidS = all(nodeValidStrain(faces), 2);
    faces_plot = faces(elemValidS, :);
    fprintf('  Disp mesh: %d/%d, Strain mesh: %d/%d elements valid\n', ...
        nnz(elemValidD), size(faces,1), nnz(elemValidS), size(faces,1));

    if size(faces_disp, 1) < 5 || size(faces_plot, 1) < 5
        fprintf('  WARNING: too few valid elements, skipping combined plot\n');
        continue;
    end

    % Prepare displacement data (zero out invalid nodes for patch interp)
    magSp = sqrt(Usp.^2 + Vsp.^2 + Wsp.^2);
    Usp_p = Usp; Usp_p(~nodeValidDisp) = 0;
    Vsp_p = Vsp; Vsp_p(~nodeValidDisp) = 0;
    Wsp_p = Wsp; Wsp_p(~nodeValidDisp) = 0;
    magSp_p = magSp; magSp_p(~nodeValidDisp) = 0;

    % Prepare strain data
    strain_exx_p  = strain_exx;  strain_exx_p(~nodeValidStrain) = 0;
    strain_eyy_p  = strain_eyy;  strain_eyy_p(~nodeValidStrain) = 0;
    strain_exy_p  = strain_exy;  strain_exy_p(~nodeValidStrain) = 0;
    strain_e1_p   = strain_e1;   strain_e1_p(~nodeValidStrain)  = 0;
    strain_e2_p   = strain_e2;   strain_e2_p(~nodeValidStrain)  = 0;
    strain_ms_p   = strain_maxshear; strain_ms_p(~nodeValidStrain) = 0;

    % All 10 fields: {data, title, faces, clim_mode}
    % clim_mode: 'auto'=auto, [lo hi]=fixed range
    allFields = { ...
        Usp_p,        'U_x',                       faces_disp, 'auto'; ...
        Vsp_p,        'U_y',                       faces_disp, 'auto'; ...
        Wsp_p,        'U_z',                       faces_disp, 'auto'; ...
        strain_e1_p,  '\epsilon_1',                faces_plot, [-2.1, 2.1]; ...
        strain_e2_p,  '\epsilon_2',                faces_plot, [-0.45, 0.45]; ...
        magSp_p,      '|U|',                       faces_disp, 'auto'; ...
        strain_exx_p, '\epsilon_{xx}',             faces_plot, [-1.3, 1.3]; ...
        strain_eyy_p, '\epsilon_{yy}',             faces_plot, [-1.3, 1.3]; ...
        strain_exy_p, '\epsilon_{xy}',             faces_plot, [-0.7, 0.7]; ...
        strain_ms_p,  'Max Shear',                 faces_plot, 'auto'; ...
    };

    % Custom coolwarm colormap (blue → white → red, 256 levels)
    n_cw = 128;
    cool_side = [interp1([0 1], [0.2298 1.0], linspace(0,1,n_cw))', ...
                 interp1([0 1], [0.2987 1.0], linspace(0,1,n_cw))', ...
                 interp1([0 1], [0.7537 1.0], linspace(0,1,n_cw))'];
    warm_side = [interp1([0 1], [1.0 0.7057], linspace(0,1,n_cw))', ...
                 interp1([0 1], [1.0 0.0156], linspace(0,1,n_cw))', ...
                 interp1([0 1], [1.0 0.1502], linspace(0,1,n_cw))'];
    cmap_coolwarm = [cool_side; warm_side];

    figure('Name', sprintf('Frame %d: Disp+Strain', f), ...
        'Position', [20, 50, 2200, 700], 'Color', 'w');
    for s = 1:10
        subplot(2, 5, s);
        patch('Faces', allFields{s,3}, 'Vertices', vertices3D, ...
            'FaceVertexCData', allFields{s,1}, 'FaceColor', 'interp', 'EdgeColor', 'none');
        colorbar; colormap(gca, cmap_coolwarm);
        title(allFields{s,2}, 'Interpreter', 'tex', 'FontSize', 11);
        axis equal; axis tight; view(2);
        set(gca, 'FontSize', 9, 'XTickLabel', [], 'YTickLabel', []);

        % Apply color limits
        cm = allFields{s,4};
        if isnumeric(cm)
            clim(cm);
        end
    end
    sgtitle(sprintf('Frame %d/%d — Displacement & Strain (Specimen Coords, VSG=%dx%d px)', ...
        f, nFrames, VSG_length, VSG_length), 'FontSize', 13, 'FontWeight', 'bold');
    drawnow;

    fprintf('  Disp (specimen): Ux=[%.2f,%.2f], Uy=[%.2f,%.2f], Uz=[%.2f,%.2f]\n', ...
        min(Usp(nodeValidDisp)), max(Usp(nodeValidDisp)), ...
        min(Vsp(nodeValidDisp)), max(Vsp(nodeValidDisp)), ...
        min(Wsp(nodeValidDisp)), max(Wsp(nodeValidDisp)));

    fprintf('  e1: [%.4f, %.4f], e2: [%.4f, %.4f]\n', ...
        min(strain_e1(validIdx)), max(strain_e1(validIdx)), ...
        min(strain_e2(validIdx)), max(strain_e2(validIdx)));
end

% --- Save with strain ---
save_name_strain = sprintf('results/noise/coarse_test_noise_ws%d_ss%d_strain.mat', WINSIZE, WINSTEPSIZE);
save(save_name_strain, 'DICpara', 'FinalResult', 'RD_L', 'RD_R', 'StereoInfo', ...
    'reprojectionErrors', 'matchedPairs', 'STRAIN_SIZE', '-v7.3');
fprintf('\nStrain results saved to: %s\n', save_name_strain);
fprintf('============================================\n');
fprintf('  Strain computation complete!\n');
fprintf('============================================\n');
