# Main Script Refactor + Git Hygiene — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Clean up `main_3D_ALDIC.m` (404 → ~250 lines) and prevent private research files from entering git.

**Architecture:** Two independent workstreams. **E** (git hygiene) is a 3-task sequential flow that happens first and is isolated. **C** (main-script refactor) has 14 independent sub-tasks grouped as C1 (cleanup, no behavior change) then C2 (user convenience). Each task touches only `main_3D_ALDIC.m` and at most one new helper file under `func/`. No modifications to the existing `func/` or `func_quadtree/` core algorithm files.

**Tech Stack:** MATLAB R2024+; git; `checkcode` for static linting.

**Design reference:** `docs/plans/2026-04-22-main-script-refactor-design.md`

**Validation strategy:** MATLAB research codebase — no pytest-style TDD. Per-task verification uses:
1. `checkcode` on the modified file (must not introduce new warnings).
2. MATLAB parse (`matlab -batch "p = mtree('main_3D_ALDIC.m'); assert(~isempty(p))"`) as a cheap syntax check.
3. **Before starting**: user runs `main_3D_ALDIC.m` on `examples/Stereo_DIC_Challenge_1.0_S3/` once, saves `FinalResult` as `baseline_before.mat` (optional but recommended).
4. **After finishing**: user re-runs same flow, saves as `baseline_after.mat`, diffs key fields.

---

## Pre-flight: Create baseline (optional, recommended)

**Why:** Most class-1 changes are pure text cleanup and should be bit-identical. If you want proof, capture output before any edit.

**Steps:**
1. In MATLAB, run `main_3D_ALDIC` once on `examples/Stereo_DIC_Challenge_1.0_S3/`. Accept defaults.
2. After completion: `save('baseline_before.mat', 'FinalResult', 'DICpara', 'RD_L', '-v7.3')`
3. After all tasks complete, re-run and `save('baseline_after.mat', 'FinalResult', 'DICpara', 'RD_L', '-v7.3')`
4. Diff:
```matlab
a = load('baseline_before.mat'); b = load('baseline_after.mat');
assert(isequal(size(a.FinalResult.Coordinates), size(b.FinalResult.Coordinates)), 'coord size mismatch');
for i = 1:numel(a.FinalResult.Displacement)
    d = max(abs(a.FinalResult.Displacement{i}(:) - b.FinalResult.Displacement{i}(:)), [], 'omitnan');
    assert(d < 1e-10, sprintf('frame %d displacement diff: %.2e', i, d));
end
fprintf('PASS: outputs bit-identical\n');
```

If you skip this, proceed — `checkcode` + manual visual review will catch the important issues.

---

## Workstream E: Git Hygiene

### Task E1: Extend `.gitignore`

**Files:**
- Modify: `.gitignore`

**Step 1: Read current `.gitignore`**

Run: `cat .gitignore`

**Step 2: Append the private-file patterns**

Append to `.gitignore`:
```gitignore

# Private research / personal work
Stereo_DIC_Challenge_2.1_Bespoke/
challenge_scripts/
run_coarse_test*.m
run_export_*.m
validate_csv_output.m
view_csv_gui.py
benchmark_single_frame.m
coarse_test_*.mat
results/
docs/plans/2026-04-16-csv-export-pipeline.md
CLAUDE.md

# Python / editor artifacts
__pycache__/
*.pyc
.claude/
```

**Step 3: Verify ignore patterns work**

Run: `git status --short | wc -l` before and after. Expected: significantly fewer untracked files shown.

Then: `git check-ignore -v run_export_csv.m CLAUDE.md Stereo_DIC_Challenge_2.1_Bespoke` — should report match for each.

**Step 4: Commit**

```bash
git add .gitignore
git -c commit.gpgsign=false commit -m "chore: gitignore private research and personal files"
```

---

### Task E2: Un-track already-committed private files

**Files:**
- Remove from git index (keep on disk):
  - `run_coarse_test_noise.m`
  - `run_export_csv.m`
  - `validate_csv_output.m`
  - `challenge_scripts/export_csv_pipeline.m`
  - `challenge_scripts/interpolate_to_dense_grid_from_mask.m`
  - `challenge_scripts/compute_strain_dense_grid_lowmem.m`
  - `challenge_scripts/export_challenge_csv.m`
  - `results/csv_output/README.md`

**Step 1: Verify files exist on disk**

Run: `ls run_coarse_test_noise.m run_export_csv.m validate_csv_output.m challenge_scripts/export_csv_pipeline.m challenge_scripts/interpolate_to_dense_grid_from_mask.m challenge_scripts/compute_strain_dense_grid_lowmem.m challenge_scripts/export_challenge_csv.m results/csv_output/README.md`

Expected: all 8 files listed with no errors.

**Step 2: Remove from git index**

```bash
git rm --cached \
  run_coarse_test_noise.m \
  run_export_csv.m \
  validate_csv_output.m \
  challenge_scripts/export_csv_pipeline.m \
  challenge_scripts/interpolate_to_dense_grid_from_mask.m \
  challenge_scripts/compute_strain_dense_grid_lowmem.m \
  challenge_scripts/export_challenge_csv.m \
  results/csv_output/README.md
```

Expected output: 8 `rm` lines.

**Step 3: Verify files still on disk**

Run same `ls` as step 1. Expected: all files still present.

**Step 4: Verify git now ignores them**

Run: `git status --short | grep -E "(run_|challenge_scripts/|csv_output)"`. Expected: empty (all caught by gitignore).

**Step 5: Commit**

```bash
git -c commit.gpgsign=false commit -m "chore: untrack private research files (kept on disk)"
```

---

### Task E3: Verify clean state

**Step 1: Check git status**

Run: `git status --short`

Expected: Only modifications to `main_3D_ALDIC.m`, `func/*.m`, etc. — no `run_export_*`, `challenge_scripts/*`, or `results/*` lines.

**Step 2: Confirm working tree builds mentally**

Spot-check: `git ls-files | grep -E "(run_export|challenge_scripts|coarse_test_|validate_csv)"` should return empty.

---

## Workstream C1: `main_3D_ALDIC.m` Cleanup (no behavior change)

All tasks modify `main_3D_ALDIC.m` only unless noted. Line numbers below refer to the file **at its current state (404 lines)**. After each task commits, subsequent line numbers will drift — use the code-context strings in each task to locate edits, not line numbers.

### Task C1.1: Remove commented-out dead code

**Files:**
- Modify: `main_3D_ALDIC.m`

**Step 1: Remove each of the following blocks**

Find and delete these commented-out lines (keep the active hardcoded defaults that follow them):

1. In Section 3.1 (calib_method):
```matlab
% calib_method = 0;
```

2. In Section 6 (after `DICpara.outputFilePath`):
```matlab
% DICpara.um2px = funParaInput('ConvertUnit');
```
(keep `DICpara.um2px = 1;`)

3. The smooth prompt:
```matlab
%DICpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot');
```
(keep `DICpara.DoYouWantToSmoothOnceMore = 0;`)

4. The strain method prompt:
```matlab
% ------ Choose strain computation method ------
%DICpara.MethodToComputeStrain = funParaInput('StrainMethodOp');
```

5. The strain type prompt:
```matlab
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
%DICpara.StrainType = funParaInput('StrainType');
```
(keep `DICpara.StrainType = 0;` and its comment)

6. Save fig format:
```matlab
%DICpara.MethodToSaveFig = 1;
```

7. Transparency block inside `if DICpara.MethodToSaveFig == 1`:
```matlab
% DICpara.OrigDICImgTransparency = funParaInput('OrigDICImgTransparency');
```

8. The smooth prompt block at the top of the for-loop (lines ~300-303):
```matlab
% prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
% DoYouWantToSmoothOnceMore = input(prompt);
% ------ Smooth displacements ------
```

9. The pre-display disp block (lines ~329-330):
```matlab
% close all; Plotdisp_show_3D(FinalResult.Displacement(ImgSeqNum,:),RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM,...
%     RD_L.ResultFEMeshEachFrame{1,1}.elementsFEM,DICpara,'NoEdgeColor');
```

10. The PlotWarpage block (lines ~355-357):
```matlab
% Zach Optional：
% PlotWarpage(FinalResult.warpage(ImgSeqNum,:),RD_L.ResultDisp{ImgSeqNum-1,1}.U,...
%     RD_L.ResultFEMeshEachFrame{1,1}, FullImageName_first,FullImageName_current,DICpara,voidIndex);
```

11. The final save block (lines ~399-401):
```matlab
% selectedFolder = uigetdir('', 'Select save folder');
% results_name = [selectedFolder,'\','results_',ImageName,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
% save(results_name, 'FullImageName_current','DICpara','FinalResult');
```

Also remove the inner `if DICpara.transformDisp == 0 ... else ... end` with all-commented branches around line 307-312 — just keep the active else branch body as the single call.

**Step 2: Lint**

Run: `matlab -batch "checkcode('main_3D_ALDIC.m', '-cyc', '-nostringlist')"`

Expected: no NEW warnings vs. before this task (pre-existing warnings OK).

**Step 3: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "refactor: remove commented-out dead code from main_3D_ALDIC"
```

---

### Task C1.2: Remove unused typo variable `shapeFuncOalized_rder`

**Files:**
- Modify: `main_3D_ALDIC.m`

**Step 1: Find and delete the line**

Find:
```matlab
shapeFuncOalized_rder = 1; % Curre1ntly, we only support 1st shape function
```

This line (in Section 4, right camera branch) is never referenced. Delete it entirely. The call below (`TemporalMatch_quadtree_ST1(..., shapeFuncOrder)`) uses `shapeFuncOrder` defined earlier in the left-camera branch, which is still in scope.

**Step 2: Verify the actual usage**

Run: `grep -n "shapeFuncO" main_3D_ALDIC.m`

Expected: exactly 3 matches (line assignment in Section 4 Left, call in Section 4 Left, call in Section 4 Right) — all using `shapeFuncOrder`, not the typo.

**Step 3: Lint**

`matlab -batch "checkcode('main_3D_ALDIC.m')"` — no new warnings.

**Step 4: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "refactor: remove unused typo variable shapeFuncOalized_rder"
```

---

### Task C1.3: Renumber Section 6 banners

**Files:**
- Modify: `main_3D_ALDIC.m`

**Step 1: Find the mis-labeled banners**

Currently in Section 6 (strain computation):
```matlab
fprintf('------------ Section 8 Start ------------ \n')
```
and near the end:
```matlab
fprintf('------------ Section 8 Done ------------ \n \n')
```

**Step 2: Replace with correct numbering**

Change `Section 8` → `Section 6` in both `fprintf` calls.

**Step 3: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "refactor: fix section numbering in strain computation banner"
```

---

### Task C1.4: Remove Challenge 2.1 Bespoke hardcoded base points

**Files:**
- Modify: `main_3D_ALDIC.m`

**Step 1: Find the leak**

In Section 6, inside `if DICpara.transformDisp == 1`:
```matlab
Base_Points2D = getBasePoints(imageLeft{1,1}',maskLeft{1,1}');

% Zach 20260127: Table values are 0-based, MATLAB is 1-based → +1
Base_Points2D = [1024+1, 1224+1; 1509.69+1, 1219.38+1; 1020.10+1, 470.48+1];
```

The second assignment **overwrites** the interactive `getBasePoints` result with values hardcoded for Challenge 2.1 Bespoke.

**Step 2: Delete the hardcoded override**

Remove these two lines:
```matlab

% Zach 20260127: Table values are 0-based, MATLAB is 1-based → +1
Base_Points2D = [1024+1, 1224+1; 1509.69+1, 1219.38+1; 1020.10+1, 470.48+1];
```

Keep the `Base_Points2D = getBasePoints(...)` call alone.

**Step 3: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "refactor: remove Challenge 2.1 hardcoded base points override"
```

---

### Task C1.5: Remove personal annotations

**Files:**
- Modify: `main_3D_ALDIC.m`

**Step 1: Find the annotations**

1. In Section 3.1 (calib_method case 0):
```matlab
% ----------  Zach comments  --------------------
% Do not use this line, please install CV toolbox and use stereo 
% camera calibrator app to calibrate your cameras.
% StereoInfo = StereoCameraCalibration; 
% ---------------------------------------
```

2. Already handled in C1.4: `% Zach 20260127: ...`

**Step 2: Rewrite the Section 3.1 comment**

Replace the "Zach comments" block with a professional comment:
```matlab
% NOTE: Install Computer Vision Toolbox and use the Stereo Camera
% Calibrator app to calibrate your cameras; then export via
% cameraParamsFormatConvertFromMatlabCV below.
```

(Also remove the `% StereoInfo = StereoCameraCalibration;` line — it's dead code.)

**Step 3: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "refactor: replace personal annotations with neutral comments"
```

---

### Task C1.6: Extract `initResultStorage` helper

**Files:**
- Create: `func/initResultStorage.m`
- Modify: `main_3D_ALDIC.m`

**Step 1: Create helper**

Write `func/initResultStorage.m`:
```matlab
function RD = initResultStorage(numImagePairs, isIncremental)
% INITRESULTSTORAGE Pre-allocate result containers for one camera.
%
%   RD = initResultStorage(numImagePairs, isIncremental)
%
% Inputs:
%   numImagePairs - length(imgNormalized_*)
%   isIncremental - true for DICpara.DICIncOrNot == 1
%
% Output:
%   RD - struct with ResultDisp, ResultDefGrad, ResultFEMeshEachFrame,
%        ResultFEMesh, and (if incremental) ResultDisp_inc.

numFrames = numImagePairs - 1;
RD.ResultDisp     = cell(numFrames, 1);
RD.ResultDefGrad  = cell(numFrames, 1);

if isIncremental
    RD.ResultFEMeshEachFrame = cell(numFrames, 1);
    RD.ResultFEMesh          = cell(numFrames, 1);
    RD.ResultDisp_inc        = cell(numFrames, 1);
else
    RD.ResultFEMeshEachFrame = cell(1, 1);
    RD.ResultFEMesh          = cell(1, 1);
end
end
```

**Step 2: Replace duplicated init in `main_3D_ALDIC.m`**

Find the block starting `RD_L.ResultDisp = cell(length(imgNormalized_L)-1,1);` through `RD_R.ResultDisp_inc = cell(length(imgNormalized_L)-1,1); end` (~20 lines, including `RD_L`, `RD_R`, and the two if/else blocks).

Replace with:
```matlab
% ====== Initialize variable storage ======
isIncremental = (DICpara.DICIncOrNot == 1);
RD_L = initResultStorage(length(imgNormalized_L), isIncremental);
RD_R = initResultStorage(length(imgNormalized_R), isIncremental);
```

**Step 3: Lint both files**

```bash
matlab -batch "checkcode('main_3D_ALDIC.m'); checkcode('func/initResultStorage.m')"
```

Expected: no errors.

**Step 4: Commit**

```bash
git add main_3D_ALDIC.m func/initResultStorage.m
git -c commit.gpgsign=false commit -m "refactor: extract initResultStorage helper from main script"
```

---

### Task C1.7: Replace hardcoded `'\'` with `fullfile()`

**Files:**
- Modify: `main_3D_ALDIC.m`

**Step 1: Find the hardcoded joins**

In the for-loop:
```matlab
FullImageName_current = [fileNameLeft{2,ImgSeqNum} ,'\', fileNameLeft{1,ImgSeqNum}];
FullImageName_first = [fileNameLeft{2,1} ,'\', fileNameLeft{1,1}];
```

**Step 2: Replace with `fullfile`**

```matlab
FullImageName_current = fullfile(fileNameLeft{2,ImgSeqNum}, fileNameLeft{1,ImgSeqNum});
FullImageName_first   = fullfile(fileNameLeft{2,1}, fileNameLeft{1,1});
```

**Step 3: Search for any other `'\'` joins**

Run: `grep -n "'\\\\'" main_3D_ALDIC.m` and `grep -nF "'\\'" main_3D_ALDIC.m`

If matches appear, apply the same `fullfile` fix.

**Step 4: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "refactor: use fullfile() instead of hardcoded path separator"
```

---

### Task C1.8: Extract 3D reconstruction verification block

**Files:**
- Create: `func/verify_stereo_calibration.m`
- Modify: `main_3D_ALDIC.m`

**Step 1: Create the helper**

Write `func/verify_stereo_calibration.m`:
```matlab
function [reconstructedPoints, reprojectionErrors] = verify_stereo_calibration(StereoInfo, showPlot)
% VERIFY_STEREO_CALIBRATION Reconstruct 3D points from stereo match
% (first frame only) and optionally display a scatter plot + mean
% reprojection error. Use after Section 3.2 as a sanity check.
%
%   [pts, errs] = verify_stereo_calibration(StereoInfo)           % plots
%   [pts, errs] = verify_stereo_calibration(StereoInfo, false)    % no plot

if nargin < 2, showPlot = true; end

RD0_L_Pts = StereoInfo.ResultFEMeshEachFrame.coordinatesFEM;
RD0_R_Pts = StereoInfo.ResultFEMesh_corr;
matchedPairs = { [RD0_L_Pts, RD0_R_Pts] };

cameraParams = StereoInfo.cameraParams;
K_left  = cameraParams.cameraParamsLeft.K;
K_right = cameraParams.cameraParamsRight.K;
R_left  = eye(3);
T_left  = zeros(3, 1);
R_right = cameraParams.rotationMatrix;
T_right = cameraParams.translationVector';

matchedPairs_undistort = funUndistortPoints(matchedPairs, cameraParams);

P_left  = K_left  * [R_left,  T_left];
P_right = K_right * [R_right, T_right];

reconstructedPoints = cell(size(matchedPairs_undistort, 1), 1);
reprojectionErrors  = cell(size(matchedPairs_undistort, 1), 1);

for i = 1:size(matchedPairs_undistort, 1)
    [reconstructedPoints{i,1}, reprojectionErrors{i}] = triangulate( ...
        matchedPairs_undistort{i,1}(:, 1:2), ...
        matchedPairs_undistort{i,1}(:, 3:4), ...
        P_left, P_right);
end

if showPlot
    figure; plot(RD0_L_Pts(:,1), RD0_L_Pts(:,2), 'o');
    hold on; plot(RD0_R_Pts(:,1), RD0_R_Pts(:,2), 'o');
    figure; scatter3(reconstructedPoints{1,1}(:,1), ...
                     reconstructedPoints{1,1}(:,2), ...
                     reconstructedPoints{1,1}(:,3));
    fprintf('mean reprojection error = %.4f\n', mean(reprojectionErrors{1}));
end
end
```

**Step 2: Replace the embedded block in `main_3D_ALDIC.m`**

Find the block between `%%%%%%%%%%%%%%%%%%%%% Test 3D construction: START %%%%%%%%%%%%%%%%%%%%%` and `%%%%%%%%%%%%%%%%%%%%% Test end %%%%%%%%%%%%%%%%%%%%%` (inclusive of both fence comments).

Replace with:
```matlab
% Optional: sanity-check stereo calibration on frame 1
if isfield(DICpara, 'verifyStereoReconstruction') && DICpara.verifyStereoReconstruction
    verify_stereo_calibration(StereoInfo);
end
```

**Step 3: Lint**

```bash
matlab -batch "checkcode('main_3D_ALDIC.m'); checkcode('func/verify_stereo_calibration.m')"
```

**Step 4: Commit**

```bash
git add main_3D_ALDIC.m func/verify_stereo_calibration.m
git -c commit.gpgsign=false commit -m "refactor: extract stereo verification block into helper function"
```

---

## Workstream C2: User Convenience

These tasks change behavior in small, intentional ways. Each is backward-compatible (defaults preserve old behavior).

### Task C2.1: Smart mex recompile

**Files:**
- Modify: `main_3D_ALDIC.m` (Section 1)

**Step 1: Replace the unconditional `mex` call**

Find:
```matlab
try
    mex -O ba_interp2_spline.cpp;
    warning('off');
    fprintf('Mex compilation successful.\n');
catch ME
    fprintf('Mex compilation failed: %s\n', ME.message);
    fprintf('Please check complier installation and path.\n');
end
```

Replace with:
```matlab
mex_src = 'ba_interp2_spline.cpp';
mex_bin = ['ba_interp2_spline.', mexext];
needRebuild = ~exist(mex_bin, 'file') || ...
              (exist(mex_bin, 'file') && ...
               dir(mex_src).datenum > dir(mex_bin).datenum);
if isfield(DICpara, 'forceMexRebuild') && DICpara.forceMexRebuild
    needRebuild = true;
end
if needRebuild
    try
        mex('-O', mex_src);
        warning('off');
        fprintf('Mex compilation successful.\n');
    catch ME
        fprintf('Mex compilation failed: %s\n', ME.message);
        fprintf('Please check compiler installation and path.\n');
    end
else
    fprintf('Mex binary up-to-date; skipping rebuild.\n');
end
```

**Note:** `DICpara` doesn't exist yet at Section 1 (loaded in Section 2). So `isfield(DICpara, 'forceMexRebuild')` will error. Replace with:
```matlab
if exist('DICpara', 'var') && isfield(DICpara, 'forceMexRebuild') && DICpara.forceMexRebuild
    needRebuild = true;
end
```

**Step 2: Test**

Delete `ba_interp2_spline.mexw64`. Run `main_3D_ALDIC` up through Section 1. Expect "Mex compilation successful". Run again. Expect "Mex binary up-to-date; skipping rebuild."

**Step 3: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "feat: skip mex rebuild if binary is up-to-date"
```

---

### Task C2.2: Make `strain_size` a `DICpara` field with prompt fallback

**Files:**
- Modify: `main_3D_ALDIC.m` (Section 6)

**Step 1: Replace the unconditional `input` prompt**

Find:
```matlab
% ------ Start main part ------
prompt = 'What is your strain size? e.g. 3,5,7...\nInput: ';
strain_size = input(prompt);
strain_length = (strain_size-1) * DICpara.winstepsize + 1;
VSG_length = strain_length + DICpara.winsize;
fprintf('Your strain size is %d * %d (Unit: Calculated Points) \nYour VSG size is %d * %d (Unit: Pixel) \n', strain_size,strain_size,VSG_length,VSG_length);
```

Replace with:
```matlab
% ------ Strain size: from DICpara if set, else prompt ------
if isfield(DICpara, 'strain_size') && ~isempty(DICpara.strain_size)
    strain_size = DICpara.strain_size;
else
    strain_size = input('What is your strain size? e.g. 3,5,7...\nInput: ');
    DICpara.strain_size = strain_size;
end
strain_length = (strain_size - 1) * DICpara.winstepsize + 1;
VSG_length    = strain_length + DICpara.winsize;
fprintf('strain_size = %d x %d pts  |  VSG = %d x %d px\n', ...
    strain_size, strain_size, VSG_length, VSG_length);
```

**Step 2: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "feat: allow preset DICpara.strain_size to skip interactive prompt"
```

---

### Task C2.3: Conditional `MW_MINGW64_LOC` setenv

**Files:**
- Modify: `main_3D_ALDIC.m` (Section 1)

**Step 1: Replace the unconditional `setenv`**

Find:
```matlab
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
```

Replace with:
```matlab
% Set MW_MINGW64_LOC only if not already set and the default path exists
if isempty(getenv('MW_MINGW64_LOC'))
    default_mingw = 'C:\TDM-GCC-64';
    if exist(default_mingw, 'dir')
        setenv('MW_MINGW64_LOC', default_mingw);
    else
        warning('MW_MINGW64_LOC not set and %s does not exist. If mex compilation fails, set this env var to your MinGW install path.', default_mingw);
    end
end
```

**Step 2: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "feat: conditionally set MW_MINGW64_LOC with fallback warning"
```

---

### Task C2.4: Default `outputFilePath` with timestamp

**Files:**
- Modify: `main_3D_ALDIC.m` (Section 6)

**Step 1: Replace the empty default**

Find:
```matlab
% Image save path
DICpara.outputFilePath = [];
```

Replace with:
```matlab
% Image save path: default to ./results/<timestamp>/ unless preset
if ~isfield(DICpara, 'outputFilePath') || isempty(DICpara.outputFilePath)
    DICpara.outputFilePath = fullfile('.', 'results', datestr(now, 'yyyymmdd_HHMMSS'));
end
if ~exist(DICpara.outputFilePath, 'dir')
    mkdir(DICpara.outputFilePath);
end
fprintf('Output folder: %s\n', DICpara.outputFilePath);
```

**Step 2: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "feat: default outputFilePath to timestamped results folder"
```

---

### Task C2.5: Add `DICpara.PlotMode` for headless support

**Files:**
- Modify: `main_3D_ALDIC.m` (Section 6, inside the for-loop)

**Step 1: Set default at top of Section 6**

After the other `DICpara.*` defaults (e.g., after the `DICpara.Image2PlotResults` assignment), add:
```matlab
% ------ Plot mode: 'show' (default) | 'save' | 'none' (headless) ------
if ~isfield(DICpara, 'PlotMode') || isempty(DICpara.PlotMode)
    DICpara.PlotMode = 'show';
end
```

**Step 2: Guard the Plot* calls inside the for-loop**

Wrap the `if DICpara.DICIncOrNot == 0 ... elseif == 1 ... end` block (the one containing `PlotdispQuadtreeMasks3D_*_ST1` and `PlotstrainQuadtreeMasks3D_*_ST1` calls) with:

```matlab
if ~strcmp(DICpara.PlotMode, 'none')
    if DICpara.DICIncOrNot == 0
        % ... existing acc branch (PlotdispQuadtreeMasks3D_acc_ST1 +
        %     PlotstrainQuadtreeMasks3D_acc_ST1)
    elseif DICpara.DICIncOrNot == 1
        % ... existing inc branch
    end
else
    % Still need to compute strain even in headless mode.
    % The PlotstrainQuadtreeMasks3D_* functions do strain computation
    % inline; when headless, we need an alternative strain-only call.
    % For v1 of this refactor: strain_* are NOT computed in 'none' mode.
    % This is a known limitation documented below.
    strain_exx = []; strain_eyy = []; strain_exy = [];
    strain_principal_max = []; strain_principal_min = [];
    strain_maxshear = []; strain_vonMises = [];
end
```

**Known limitation:** `PlotstrainQuadtreeMasks3D_*_ST1` couples strain computation with plotting. Truly headless strain requires a separate refactor (outside this task). Current v1 just skips both when `PlotMode = 'none'`; `FinalResult.ResultStrainWorld{ImgSeqNum,1}` will contain empty fields. Add a warning when the user sets `'none'`:

```matlab
if strcmp(DICpara.PlotMode, 'none')
    warning('PlotMode=''none'' also skips strain computation in v1. Set to ''save'' for strain + saved figures without display.');
end
```

Place this warning once, right after the PlotMode default assignment.

**Step 3: Also guard `close all` inside loop**

Find `close all;` inside the for-loop (lines ~291 and 391). Wrap:
```matlab
if ~strcmp(DICpara.PlotMode, 'none')
    close all;
end
```

**Step 4: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "feat: add DICpara.PlotMode for save/headless modes"
```

---

### Task C2.6: Uniform section banners with elapsed time

**Files:**
- Modify: `main_3D_ALDIC.m`

**Step 1: At the very top of Section 1**

Add:
```matlab
session_start = tic;
numSections = 6;
```

**Step 2: Define a local helper at top of file (before Section 1)**

Right after the big header comment block (line ~39), add:
```matlab
sectionBanner = @(n, name) fprintf('\n==== Section %d/%d: %s  |  elapsed: %.1fs ====\n', ...
    n, numSections, name, toc(session_start));
```

**Note:** This creates an anonymous function that captures `session_start` and `numSections`. Requires both to be defined first — so put the two declarations (`session_start = tic; numSections = 6;`) BEFORE the anonymous function.

**Step 3: Replace each `fprintf('------------ Section X Start ------------ \n')` and `Done` variant**

Section 1 Start → `sectionBanner(1, 'Environment & mex setup');`
Section 2 Start → `sectionBanner(2, 'Load images & DIC parameters');`
Section 3 Start (combine 3.1 + 3.2 conceptually; keep sub-banners) → `sectionBanner(3, 'Stereo calibration & matching');`
Section 4 Start → `sectionBanner(4, 'Temporal matching');`
Section 5 Start → `sectionBanner(5, '3D reconstruction');`
Section 6 Start → `sectionBanner(6, 'Strain & visualization');`

Remove the old `------------ Section X Start ------------ \n` and `Done` lines; the new single banner per section is enough.

At the very end of the script, add:
```matlab
fprintf('\n==== ALL DONE  |  total elapsed: %.1fs ====\n', toc(session_start));
```

**Step 4: Commit**

```bash
git add main_3D_ALDIC.m
git -c commit.gpgsign=false commit -m "feat: add uniform section banners with elapsed time"
```

---

## Final: Post-refactor validation

### Task V1: Static check + line count

**Step 1: Lint full file**

```bash
matlab -batch "msgs = checkcode('main_3D_ALDIC.m'); arrayfun(@(m) fprintf('%d: %s\n', m.line, m.message), msgs)"
```

Expected: warnings allowed only if they pre-existed. No `err` level messages.

**Step 2: Line count check**

```bash
wc -l main_3D_ALDIC.m
```

Expected: ~230-270 lines (was 404). If still > 300, spot-check what's left.

### Task V2: Runtime verification (user-driven)

**Step 1:** Run `main_3D_ALDIC` on `examples/Stereo_DIC_Challenge_1.0_S3/` with same inputs as `baseline_before.mat` capture.

**Step 2:** Save `FinalResult` as `baseline_after.mat`.

**Step 3:** Compare against `baseline_before.mat` using the diff script in the Pre-flight section.

Expected: max displacement diff < 1e-10.

**Step 4: Final commit (if nothing to commit, skip)**

After confirming validation, you're done.

---

## Summary

| Workstream | Tasks | Est. time |
|---|---|---|
| E (git hygiene) | 3 | 15 min |
| C1 (cleanup) | 8 | 60-90 min |
| C2 (user convenience) | 6 | 60-90 min |
| Validation | 2 | 30 min |
| **Total** | **19** | **~4 hours** |

Each task is independently committable and reversible via `git revert`.
