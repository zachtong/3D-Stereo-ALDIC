# Stereo-ALDIC Refactor — Implementation Plan (v2)

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans.

**Goal:** Execute the 8 workstreams from `docs/plans/2026-04-22-main-script-refactor-design.md`.

**Supersedes:** v1 plan at `2026-04-22-main-script-refactor-plan.md` (superseded by this v2).

**Architecture:** Each workstream is a sequence of small, independently-committable tasks. After each workstream, user runs `tests/run_pipeline_test.m` in MATLAB; if it passes, proceed to next workstream.

**Validation gate after every workstream:** user runs test; if baseline diff > tolerance, roll back and debug.

---

## Execution order

1. **T** Test (already written, user captures baseline first)
2. **E** Git hygiene (3 tasks)
3. **C6** Defaults + validation (2 tasks) — needed by C2/C4
4. **C1** Cleanup (cluster into 5 commits)
5. **C5** Smoothing params (2 tasks)
6. **C2** User convenience (5 tasks)
7. **C4** Light GUI guard (1 task)
8. **C3** Architectural (5 tasks)

---

## T — Test script

**Status:** ✅ Written as `tests/run_pipeline_test.m` (commit `f740bb2`).

**User action required NOW before proceeding:**
```
cd tests
run_pipeline_test
```
On first run, this captures `tests/baseline/baseline.mat`. Report success or any errors back.

---

## E — Git hygiene (3 commits)

### E1: Extend `.gitignore`
Append the 13 private-file patterns (see design doc).
Commit: `chore: gitignore private research and personal files`

### E2: Un-track committed private files
```
git rm --cached run_coarse_test_noise.m run_export_csv.m validate_csv_output.m \
  challenge_scripts/export_csv_pipeline.m \
  challenge_scripts/interpolate_to_dense_grid_from_mask.m \
  challenge_scripts/compute_strain_dense_grid_lowmem.m \
  challenge_scripts/export_challenge_csv.m \
  results/csv_output/README.md
```
Commit: `chore: untrack private research files`

### E3: Verify
`git status` should show no `run_export_*`, `challenge_scripts/*`, or `results/csv_output/*` entries.

---

## C6 — Defaults + validation (2 commits)

### C6.1: Create `func/setDICparaDefaults.m`
```matlab
function DICpara = setDICparaDefaults()
% Returns a DICpara struct with all fields initialized to defaults.
% Mirrors pyALDIC's al_dic.core.config.dicpara_default().

DICpara = struct();

% --- Algorithm ---
DICpara.winsize              = [];      % user must set
DICpara.winstepsize          = [];      % user must set
DICpara.winsizeMin           = 8;
DICpara.DICIncOrNot          = 0;       % 0=acc, 1=inc
DICpara.NewFFTSearch         = 1;
DICpara.UseGlobal            = true;
DICpara.Subpb2FDOrFEM        = 1;
DICpara.GaussPtOrder         = 2;
DICpara.ClusterNo            = 1;

% --- Image ---
DICpara.ImgRefMask           = [];      % set by caller from mask
DICpara.ImgSize              = [];
DICpara.gridxyROIRange       = struct('gridx', [], 'gridy', []);
DICpara.InitFFTSearchMethod  = 1;
DICpara.LoadImgMethod        = 1;

% --- Strain / output ---
DICpara.strain_size          = [];      % prompts if empty
DICpara.transformDisp        = 0;
DICpara.um2px                = 1;
DICpara.StrainType           = 0;       % 0 = Green-Lagrange (only supported)

% --- Smoothing (independent 2D/3D counts) ---
DICpara.Smooth2DTimes        = 0;       % inside ADMM, 0 = disabled
DICpara.Smooth3DTimes        = 3;       % post-ALDIC, current default
DICpara.DispSmoothness       = 0;
DICpara.StrainSmoothness     = 0;
DICpara.DispFilterSize       = 0;
DICpara.DispFilterStd        = 0;
DICpara.StrainFilterSize     = 0;
DICpara.StrainFilterStd      = 0;

% --- Plotting / output ---
DICpara.PlotMode             = 'show';  % 'none' | 'save' | 'show'
DICpara.showImgOrNot         = 0;       % debug figures inside algorithms
DICpara.plots_disp_to_generate   = {'u', 'v', 'w', 'magnitude'};
DICpara.plots_strain_to_generate = {'exx', 'eyy', 'exy', 'e1', 'e2'};
DICpara.Image2PlotResults    = 1;
DICpara.MethodToSaveFig      = 'png';
DICpara.OrigDICImgTransparency = 0.8;
DICpara.outputFilePath       = '';      % default fills with timestamp at runtime

% --- MEX ---
DICpara.forceMexRebuild      = false;

% --- Optional callback ---
DICpara.progressCallback     = [];      % @(frac, msg) printf()
end
```
Commit: `feat: add setDICparaDefaults.m (pyALDIC dicpara_default mirror)`

### C6.2: Create `func/validateDICpara.m`
```matlab
function validateDICpara(DICpara)
% Raises error on common misconfigurations.
% Mirrors pyALDIC's al_dic.core.config.validate_dicpara().

if ~isempty(DICpara.winsize) && (mod(DICpara.winsize, 2) ~= 0 || DICpara.winsize <= 0)
    error('winsize must be positive even integer, got %g', DICpara.winsize);
end
if ~isempty(DICpara.winstepsize)
    isPow2 = @(v) v > 0 && bitand(v, v-1) == 0;
    if ~isPow2(DICpara.winstepsize)
        error('winstepsize must be a positive power of 2, got %g', DICpara.winstepsize);
    end
    if DICpara.winsizeMin > DICpara.winstepsize
        error('winsizeMin (%d) must be <= winstepsize (%d)', ...
            DICpara.winsizeMin, DICpara.winstepsize);
    end
end
if ~ismember(DICpara.DICIncOrNot, [0 1])
    error('DICIncOrNot must be 0 (acc) or 1 (inc), got %g', DICpara.DICIncOrNot);
end
if ~ismember(DICpara.PlotMode, {'none', 'save', 'show'})
    error('PlotMode must be ''none'', ''save'', or ''show''');
end
if ~ismember(DICpara.StrainType, 0:3)
    error('StrainType must be 0-3');
end
if DICpara.Smooth2DTimes < 0 || DICpara.Smooth3DTimes < 0
    error('SmoothNTimes must be non-negative');
end
end
```
Commit: `feat: add validateDICpara.m for config validation`

---

## C1 — Cleanup (5 commits, grouped)

### C1a: Remove dead commented code from `main_3D_ALDIC.m`
See design doc C1 section 1.1-1.5, 1.7. Remove ~11 blocks.
Commit: `refactor(main): remove dead commented code, typos, hardcoded leaks`

### C1b: Remove dead code from `TemporalMatch_quadtree_ST1.m`
- L30-42 + L201-213: DataDriven mode (Q1 = don't need)
- L132-153: commented RBF initial guess
- L370-382: adaptive winsize commented code
- L349, 404, 422, 444, 585: commented `save(...)` debug dumps
- L63-77: replace "Zach Attention" hack with proper note (real fix comes in C3.2)
Commit: `refactor(temporal): remove data-driven mode, dead RBF/winsize blocks`

### C1c: Remove dead code from `StereoMatch_STAQ.m`, `stereoReconstruction_quadtree.m`, `setDICParas_STAQ.m`, `setDICParas_IncOrNot.m`, `funParaInput.m`
- StereoMatch L53-55: commented ResultFEMesh
- stereoReconstruction: drop `computeFundamentalEssentialMatrix` + `skewSymmetric` (unused)
- setDICParas_STAQ: drop L10-27 (old ROI UI) + L63-64 double-assign + L78-127 (migrated)
- setDICParas_IncOrNot: drop try/catch around switch (unreachable)
- funParaInput: fix L27 variable name, drop L223-230 `MaterialModel` case, fix "Eluerian" typo
- Remove 2nd-order shape function branches in `main` + `StereoMatch_STAQ`
Commit: `refactor: remove dead code from calib/params/stereo helpers`

### C1d: Remove dead code from `PlaneFit3_Quadtree.m`, plot functions
- PlaneFit3_Quadtree L21-68: 48-line commented manual-selection block
- PlaneFit3_Quadtree L93-126: Chinese comments → English
- PlaneFit3_Quadtree L116-117: `max(..)-1` → `find(..,1,'first')-1`
- PlotdispQuadtreeMasks3D_acc_ST1 L29-57: `JY!!!Mask` block (user confirmed remove)
- PlotdispQuadtreeMasks3D_inc_ST1: remove analogous mask-out block if present
Commit: `refactor(strain/plot): remove dead code, translate Chinese comments`

### C1e: Remove all remaining Chinese comments across codebase
Grep for `[\u4e00-\u9fff]` and translate to English.
Commit: `refactor: translate all remaining Chinese comments to English`

---

## C5 — Smoothing parameterization (2 commits)

### C5.1: Replace hardcoded `SmoothTimes < 3` in `main_3D_ALDIC.m` L307
```matlab
while DICpara.Smooth3DTimes > 0 && SmoothTimes < DICpara.Smooth3DTimes
    FinalResult.Displacement_smooth(ImgSeqNum,:) = funSmoothDisp_Quadtree(...);
    SmoothTimes = SmoothTimes + 1;
end
```
Also update `funSmoothDisp_Quadtree.m` to respect `DICpara.Smooth3DTimes` internal loop (L101 `if SmoothTimes > 2 break` → parameterized).
Commit: `feat(smooth): parameterize 3D smoothing count via DICpara.Smooth3DTimes`

### C5.2: Parameterize 2D smoothing in `TemporalMatch_quadtree_ST1.m`
Currently `DICpara.DispSmoothness=0; StrainSmoothness=0` is hardcoded L286-287. Remove those lines and rely on the defaults from C6.1 (which keep them 0). Update `funSmoothDispQuadtree.m` L86 to respect `DICpara.Smooth2DTimes`.
Commit: `feat(smooth): parameterize 2D smoothing count via DICpara.Smooth2DTimes`

---

## C2 — User convenience (5 commits)

### C2.1: Smart mex recompile in `main_3D_ALDIC.m`
Check `dir(mex_bin).datenum > dir(mex_src).datenum`, skip if up-to-date. Respect `DICpara.forceMexRebuild`.
Commit: `feat(main): skip mex rebuild if binary up-to-date`

### C2.2: Conditional MinGW setenv
Replace unconditional `setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')` with `if isempty(getenv(...)) && exist(path, 'dir'), setenv; end`. Emit warning if neither.
Commit: `feat(main): conditional MW_MINGW64_LOC setenv with fallback`

### C2.3: Default `outputFilePath` with timestamp
When `DICpara.outputFilePath` is empty at Section 6, fill with `./results/<yyyymmdd_HHMMSS>/` and `mkdir`.
Commit: `feat(main): default outputFilePath to timestamped results folder`

### C2.4: Section banners with elapsed time
Define anonymous function `sectionBanner = @(n, name) fprintf(...)` at top of script; replace all `------ Section X Start ------` strings. Add final "ALL DONE" with total elapsed.
Commit: `feat(main): uniform section banners with elapsed time`

### C2.5: Promote hardcoded overrides to `DICpara` fields
- `TemporalMatch_quadtree_ST1.m` L9: `UseGlobal = DICpara.UseGlobal;` (with default true from C6.1)
- `TemporalMatch_quadtree_ST1.m` L10: remove `DICpara.showImgOrNot = 0;` (respect caller)
- `TemporalMatch_quadtree_ST1.m` L103: remove `DICpara.NewFFTSearch = 1;` (respect setDICparaDefaults value)
- `setDICParas_IncOrNot.m` L38: remove `NewFFTSearch = 1; %tbd` (same field, different place)
Commit: `feat: promote UseGlobal/showImgOrNot/NewFFTSearch to DICpara fields`

---

## C4 — Light GUI guard (1 commit)

### C4.1: Each prompt guarded by `isempty(DICpara.X)`
- `main_3D_ALDIC.m` L277 (strain_size): `if isempty(DICpara.strain_size), prompt, DICpara.strain_size = ..., end`
- Section 3.1 calib method: `if isempty(DICpara.calibrationMethod), funParaInput, end`
- Section 6: `Image2PlotResults`, `SaveFigFormat`, `TransformDispOrNot` — all guarded
- `setDICParas_STAQ.m` L53 winsize, L58 winstepsize, `setDICParas_IncOrNot` L9 incOrNot — all guarded

This enables batch mode: pre-populate DICpara → main runs end-to-end without prompts.
Commit: `feat: skip interactive prompts when DICpara fields are preset`

---

## C3 — Architectural (5 commits)

### C3.1: Extract `func_quadtree/computeStrain3D.m` (B1)
```matlab
function [exx, eyy, exy, e1, e2, maxShear, vonMises, dwdx, dwdy] = ...
    computeStrain3D(coefficients)
% Computes Green-Lagrange strain components from displacement
% gradient coefficients. See PlaneFit3_Quadtree.m for format.

n = size(coefficients, 1);
[exx, eyy, exy, dwdx, dwdy] = deal(zeros(n, 1));
for i = 1:n
    u_x = coefficients{i,1}(1,1); u_y = coefficients{i,1}(2,1); u_z = coefficients{i,1}(3,1);
    v_x = coefficients{i,1}(1,2); v_y = coefficients{i,1}(2,2); v_z = coefficients{i,1}(3,2);
    w_x = coefficients{i,1}(1,3); w_y = coefficients{i,1}(2,3); w_z = coefficients{i,1}(3,3);
    F = [1+u_x, u_y, u_z; v_x, 1+v_y, v_z; w_x, w_y, 1+w_z];  % unified (A9.3)
    E = 0.5 * (F' * F - eye(3));
    exx(i) = E(1,1); eyy(i) = E(2,2); exy(i) = E(1,2);
    dwdx(i) = w_x;   dwdy(i) = w_y;
end
maxShear = sqrt((0.5*(exx-eyy)).^2 + exy.^2);
e1 = 0.5*(exx+eyy) + maxShear;
e2 = 0.5*(exx+eyy) - maxShear;
vonMises = sqrt(e1.^2 + e2.^2 - e1.*e2 + 3*maxShear.^2);
end
```
Update `PlotstrainQuadtreeMasks3D_acc_ST1.m` and `_inc_ST1.m` to call `computeStrain3D` (drop their inline strain loops).
Commit: `refactor(strain): extract computeStrain3D, unify acc/inc F matrix`

### C3.2: Right-camera DICpara isolation (Q5=B)
Change `TemporalMatch_quadtree_ST1.m` signature from:
```
function RD = TemporalMatch_quadtree_ST1(DICpara, fileName, ImgMask, ImgNormalized, RD, stereoInfo, camera0OrNot, shapeFuncOrder)
```
to:
```
function RD = TemporalMatch_quadtree_ST1(DICpara, fileName, ImgMask, ImgNormalized, RD, stereoInfo, camera0OrNot, shapeFuncOrder, refMask, roiRange)
```
Caller (main + test) passes per-camera `refMask` and `roiRange`. Remove the "fake update" hack at L63-77; inside function, use the passed `refMask`/`roiRange` directly.
Commit: `refactor(temporal): isolate right-camera mask/ROI via explicit args`

### C3.3: RBF → `scatteredInterpolant` in strain plot inc version (B5)
In `PlotstrainQuadtreeMasks3D_inc_ST1.m` L43-56, replace `rbfcreate`/`rbfinterp` with `scatteredInterpolant('natural', 'none')`.
Commit: `perf: replace RBF with scatteredInterpolant in strain plot (inc)`

### C3.4: Fix inc plot position (Q6=b)
In `PlotstrainQuadtreeMasks3D_inc_ST1.m` L38, change:
```matlab
coordinatesFEMWorldDef = [CurrentFEM.coordinatesFEM(:,1)+U_2D_L_inc(1:2:end), ...]
```
to:
```matlab
coordinatesFEMWorldDef = CurrentFEM.coordinatesFEM;
```
Commit: `fix(plot): use CurrentFEM position directly for inc strain plot`

### C3.5: Main script extractions
- `func/initResultStorage.m` (from v1 design C1.6)
- `func/verify_stereo_calibration.m` (from v1 design C1.8)
- Update `main_3D_ALDIC.m` to call them
Commit: `refactor(main): extract initResultStorage + verify_stereo_calibration helpers`

---

## Final validation

After all workstreams:
1. `tests/run_pipeline_test.m` passes with < 20% timing regression and zero disp/coord regression
2. `matlab -batch "checkcode('main_3D_ALDIC.m')"` shows no new warnings
3. `wc -l main_3D_ALDIC.m` shows ~230-270 (down from 404)
4. User manually runs `main_3D_ALDIC` interactively — all prompts still work

---

## Rollback protocol

If test fails after a commit:
```
git reset --hard HEAD~1   # roll back last commit (preserves later commits on a branch)
```
Debug before re-applying. All commits are small, so rollback is cheap.
