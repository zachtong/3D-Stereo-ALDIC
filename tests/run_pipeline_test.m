%% run_pipeline_test.m
%
% Non-interactive regression test for the stereo-ALDIC pipeline.
%
% Behavior:
%   - First run: saves outputs to tests/baseline/baseline.mat
%   - Subsequent runs: diffs current outputs vs baseline; fails on regression.
%
% Dataset: examples/Stereo_DIC_Challenge_1.0_S3/ (3 frames, DICe calib).
% Runtime: ~3-6 minutes per run depending on hardware.
%
% Usage:
%   >> cd tests
%   >> run_pipeline_test
%
% Env options (use setenv before run):
%   RECAPTURE_BASELINE=1  force overwrite of baseline.mat
%
% Tolerances:
%   disp    <= 1e-8 absolute
%   strain  <= 1e-10 absolute (coefficients only — strain itself not in this test)
%   timing  regression alert if >20% slower than baseline

close all; clear; clc;

%% ========================================
%  0. Paths + env
%  ========================================
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
baseline_path = fullfile(script_dir, 'baseline', 'baseline.mat');

cd(project_root);
addpath('./examples', './func', './func_quadtree', ...
        './func_quadtree/rbfinterp', './func_quadtree/refinement', ...
        './plotFiles', './plotFiles/export_fig-d966721');

RECAPTURE = ~isempty(getenv('RECAPTURE_BASELINE')) && ...
            strcmp(getenv('RECAPTURE_BASELINE'), '1');
BASELINE_EXISTS = exist(baseline_path, 'file') == 2;

if BASELINE_EXISTS && ~RECAPTURE
    mode = 'COMPARE';
else
    mode = 'CAPTURE';
end
fprintf('\n========== Test mode: %s ==========\n', mode);

overall_tic = tic;

%% ========================================
%  1. Ensure mex is compiled
%  ========================================
mex_src = 'ba_interp2_spline.cpp';
mex_bin = ['ba_interp2_spline.', mexext];
if ~exist(mex_bin, 'file')
    fprintf('Compiling %s...\n', mex_src);
    try, mex('-O', mex_src); catch ME
        error('mex compile failed: %s', ME.message);
    end
end

%% ========================================
%  2. Fixed preset DICpara
%  ========================================
% Preset values chosen to match stored reference jpg (WS64_ST64) but
% reduced winsize for quicker test runs.

DICpara = struct();
DICpara.winsize         = 32;      % smaller than reference 64 for speed
DICpara.winstepsize     = 32;
DICpara.winsizeMin      = 8;
DICpara.DICIncOrNot     = 1;       % incremental (example has per-frame masks)
DICpara.Subpb2FDOrFEM   = 1;       % FEM
DICpara.ClusterNo       = 1;       % serial (reproducible)
DICpara.NewFFTSearch    = 1;
DICpara.ImgSeqIncUnit   = 1;
DICpara.ImgSeqIncROIUpdateOrNot = 0;
DICpara.LoadImgMethod   = 1;       % direct load
DICpara.InitFFTSearchMethod = 1;
DICpara.transformDisp   = 0;       % no coord transform in test
DICpara.showImgOrNot    = 0;       % headless

% FFT search zone: preset to bypass the interactive prompt in
% IntegerSearchQuadtree. 60 px is generous for stereo disparity on a
% 1920x1200 image and covers temporal motion between 3 near-identical frames.
DICpara.NewFFTSearchDistance   = [60, 60];
DICpara.fixSearchDistanceOrNot = 0;   % reuse preset for all subsequent frames

% Tuning fields (match current TemporalMatch defaults)
DICpara.DispFilterSize  = 0;
DICpara.DispFilterStd   = 0;
DICpara.StrainFilterSize = 0;
DICpara.StrainFilterStd  = 0;
DICpara.DispSmoothness   = 0;
DICpara.StrainSmoothness = 0;
DICpara.GaussPtOrder     = 2;

%% ========================================
%  3. Load images + masks
%  ========================================
example_dir = fullfile(project_root, 'examples', 'Stereo_DIC_Challenge_1.0_S3');
img_left_dir   = fullfile(example_dir, 'Images_Stereo_Sample3_images',    'Left');
img_right_dir  = fullfile(example_dir, 'Images_Stereo_Sample3_images',    'Right');
mask_left_dir  = fullfile(example_dir, 'Images_Stereo_Sample3_maskfiles', 'Left');
mask_right_dir = fullfile(example_dir, 'Images_Stereo_Sample3_maskfiles', 'Right');
calib_file     = fullfile(example_dir, 'calibration_DICe.xml');

frame_ids = {'0000', '0001', '0002'};
nFrames = length(frame_ids);

imageLeft  = cell(1, nFrames);
imageRight = cell(1, nFrames);
maskLeft   = cell(1, nFrames);
maskRight  = cell(1, nFrames);
fileNameLeft  = cell(2, nFrames);
fileNameRight = cell(2, nFrames);

for i = 1:nFrames
    fn_L = sprintf('%s_0.tif', frame_ids{i});
    fn_R = sprintf('%s_1.tif', frame_ids{i});
    imageLeft{i}  = double(imread(fullfile(img_left_dir,  fn_L))');
    imageRight{i} = double(imread(fullfile(img_right_dir, fn_R))');
    maskLeft{i}   = (imread(fullfile(mask_left_dir,  fn_L))' > 0);
    maskRight{i}  = (imread(fullfile(mask_right_dir, fn_R))' > 0);
    fileNameLeft{1, i}  = fn_L;  fileNameLeft{2, i}  = img_left_dir;
    fileNameRight{1, i} = fn_R;  fileNameRight{2, i} = img_right_dir;
end

DICpara.ImgRefMask = double(maskLeft{1});
DICpara.ImgSize    = size(imageLeft{1});

% ROI range from first left mask
[rows, cols] = find(maskLeft{1});
DICpara.gridxyROIRange.gridx = [min(rows), max(rows)];
DICpara.gridxyROIRange.gridy = [min(cols), max(cols)];

fprintf('Loaded %d frames, size=%dx%d\n', nFrames, DICpara.ImgSize);

%% ========================================
%  4. Normalize images
%  ========================================
tic;
[imgNormalized_L, DICpara.gridxyROIRange] = funNormalizeImg(imageLeft,  DICpara.gridxyROIRange);
[imgNormalized_R, DICpara.gridxyROIRange] = funNormalizeImg(imageRight, DICpara.gridxyROIRange);
timing.normalize = toc;

% Storage
numFrames = nFrames - 1;
RD_L = struct();
RD_L.ResultDisp = cell(numFrames, 1);
RD_L.ResultDefGrad = cell(numFrames, 1);
RD_L.ResultFEMeshEachFrame = cell(numFrames, 1);
RD_L.ResultFEMesh = cell(numFrames, 1);
RD_L.ResultDisp_inc = cell(numFrames, 1);

RD_R = RD_L;

%% ========================================
%  5. Stereo calibration
%  ========================================
tic;
stereoInfo = struct();
stereoInfo.cameraParams = cameraParamsFormatConvertFromDICe( ...
    fullfile(example_dir), 'calibration_DICe.xml');
timing.calib = toc;
fprintf('Calibration loaded: %.2fs\n', timing.calib);

%% ========================================
%  6. Stereo match (frame 1)
%  ========================================
tic;
[stereoInfo, RD_L, RD_R] = StereoMatch_STAQ( ...
    RD_L, RD_R, imgNormalized_L{1}, imgNormalized_R{1}, ...
    fileNameLeft, maskLeft{1}, maskRight{1}, DICpara, stereoInfo, 1);
timing.stereo_match = toc;
fprintf('Stereo match: %.2fs\n', timing.stereo_match);

%% ========================================
%  7. Temporal match (left + right)
%  ========================================
tic;
RD_L = TemporalMatch_quadtree_ST1(DICpara, fileNameLeft, maskLeft, ...
    imgNormalized_L, RD_L, stereoInfo, 'camera0', 1);
timing.temporal_L = toc;
fprintf('Temporal match (L): %.2fs\n', timing.temporal_L);

tic;
RD_R = TemporalMatch_quadtree_ST1(DICpara, fileNameRight, maskRight, ...
    imgNormalized_R, RD_R, stereoInfo, 'notCamera0', 1);
timing.temporal_R = toc;
fprintf('Temporal match (R): %.2fs\n', timing.temporal_R);

%% ========================================
%  8. Organize pairs + 3D reconstruction
%  ========================================
tic;
matchedPairs = organizeMatchedPairds_quadtree_ST1(RD_L, RD_R, DICpara);
[FinalResult, reprojectionErrors] = stereoReconstruction_quadtree( ...
    matchedPairs, stereoInfo.cameraParams);
timing.recon3d = toc;
fprintf('3D reconstruction: %.2fs\n', timing.recon3d);

timing.total = toc(overall_tic);
fprintf('\n========== Pipeline complete: %.2fs ==========\n', timing.total);

%% ========================================
%  9. Baseline capture or compare
%  ========================================

% Build snapshot of key outputs (small arrays only, not full meshes)
snapshot = struct();
snapshot.timing = timing;
snapshot.nFrames = nFrames;
snapshot.DICpara_winsize = DICpara.winsize;

% FinalResult: coords + disp (per-frame, per-axis)
snapshot.Coordinates = cell(nFrames, 3);
snapshot.Displacement = cell(nFrames, 3);
for i = 1:nFrames
    for axis = 1:3
        snapshot.Coordinates{i, axis}  = FinalResult.Coordinates{i, axis};
        snapshot.Displacement{i, axis} = FinalResult.Displacement{i, axis};
    end
end

% RD_L sizes + hashes of ResultDisp
snapshot.RD_L_nMeshPts = zeros(numFrames, 1);
snapshot.RD_L_U_hash = cell(numFrames, 1);
for k = 1:numFrames
    if isstruct(RD_L.ResultDisp{k}) && isfield(RD_L.ResultDisp{k}, 'U')
        U = RD_L.ResultDisp{k}.U;
        snapshot.RD_L_nMeshPts(k) = numel(U);
        snapshot.RD_L_U_hash{k} = sum(U(isfinite(U)));  % cheap content hash
    end
end

snapshot.reprojErr_mean = cellfun(@(e) mean(e, 'omitnan'), reprojectionErrors);

if strcmp(mode, 'CAPTURE')
    % ===== Capture mode: save baseline =====
    if ~exist(fileparts(baseline_path), 'dir'), mkdir(fileparts(baseline_path)); end
    save(baseline_path, '-struct', 'snapshot', '-v7.3');
    fprintf('\n✓ Baseline saved to %s\n', baseline_path);
    fprintf('  Run again (without RECAPTURE_BASELINE=1) to test for regression.\n');
    return;
end

% ===== Compare mode =====
fprintf('\n========== Regression diff ==========\n');
baseline = load(baseline_path);

DISP_TOL = 1e-8;
COORD_TOL = 1e-8;
TIME_REGRESS_PCT = 20;

failed = false;

% Diff Coordinates + Displacement
max_coord_diff = 0;
max_disp_diff = 0;
for i = 1:nFrames
    for axis = 1:3
        c_now = snapshot.Coordinates{i, axis};
        c_ref = baseline.Coordinates{i, axis};
        d_now = snapshot.Displacement{i, axis};
        d_ref = baseline.Displacement{i, axis};
        if ~isequal(size(c_now), size(c_ref))
            fprintf('  FAIL  frame %d axis %d coord SIZE mismatch: %s vs %s\n', ...
                i, axis, mat2str(size(c_now)), mat2str(size(c_ref)));
            failed = true;
            continue;
        end
        cd = max(abs(c_now(:) - c_ref(:)), [], 'omitnan');
        dd = max(abs(d_now(:) - d_ref(:)), [], 'omitnan');
        max_coord_diff = max(max_coord_diff, cd);
        max_disp_diff  = max(max_disp_diff,  dd);
    end
end
fprintf('  max |coord diff| = %.3e   (tol %.1e)\n', max_coord_diff, COORD_TOL);
fprintf('  max |disp  diff| = %.3e   (tol %.1e)\n', max_disp_diff,  DISP_TOL);
if max_coord_diff > COORD_TOL, fprintf('  FAIL: coord regression\n'); failed = true; end
if max_disp_diff  > DISP_TOL,  fprintf('  FAIL: disp regression\n');  failed = true; end

% Diff ResultDisp hashes (catches any inner-loop change)
for k = 1:numFrames
    if ~isempty(snapshot.RD_L_U_hash{k}) && ~isempty(baseline.RD_L_U_hash{k})
        h_now = snapshot.RD_L_U_hash{k};
        h_ref = baseline.RD_L_U_hash{k};
        if abs(h_now - h_ref) > 1e-6
            fprintf('  FAIL  RD_L frame %d U hash diff: %.3e vs %.3e\n', ...
                k, h_now, h_ref);
            failed = true;
        end
    end
end

% Diff reprojection error
reproj_diff = max(abs(snapshot.reprojErr_mean - baseline.reprojErr_mean));
fprintf('  max |reproj err diff| = %.3e\n', reproj_diff);

% Timing regression
t_stages = {'normalize', 'calib', 'stereo_match', 'temporal_L', 'temporal_R', 'recon3d', 'total'};
fprintf('\n  Timing (current / baseline / delta):\n');
for i = 1:length(t_stages)
    s = t_stages{i};
    t_now = snapshot.timing.(s);
    t_ref = baseline.timing.(s);
    pct = (t_now - t_ref) / t_ref * 100;
    marker = '';
    if pct > TIME_REGRESS_PCT, marker = ' <-- REGRESSION'; end
    fprintf('    %-14s %7.2fs  %7.2fs  %+6.1f%%%s\n', s, t_now, t_ref, pct, marker);
end

if failed
    fprintf('\n✗ TEST FAILED\n');
    error('Regression detected — see diffs above');
else
    fprintf('\n✓ ALL CHECKS PASSED\n');
end
