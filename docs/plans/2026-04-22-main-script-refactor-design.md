# Design: 3D-Stereo-ALDIC Refactor (Final)

**Date:** 2026-04-22 (updated 2026-04-23)
**Status:** Approved for implementation
**Supersedes:** previous design at `2026-04-22-main-script-refactor-design.md` (v1)

## Goal

Make the current MATLAB implementation cleaner, faster to use, and less painful to run — while keeping the eventual Python port (mirror of `pyALDIC`) in mind. No algorithm changes. No feature additions. No file reorganization yet.

## Non-goals

- Porting any logic to Python.
- Restructuring folders into `src/al_dic_stereo/` (deferred to future Python-port work).
- Adding new features, file formats, or output paths.
- Performance parallelization (parfor for L/R temporal match — skipped).
- 2nd-order shape functions (not implemented; remove placeholders).
- Data-driven initial-guess mode (not needed; remove placeholders).
- Adaptive winsize (not needed; remove commented placeholders).

---

## Guiding principles (for Python-port compatibility)

`pyALDIC` uses:
- **Frozen dataclasses** for all data — no mutation.
- **Pure functions** with explicit args/returns — no hidden state.
- `dicpara_default()` + `validate_dicpara()` as separate concerns.
- **Progress callbacks** as optional kwargs.
- Clean `core / solver / mesh / strain / io / export / gui` separation (deferred for us).

MATLAB-specific adaptations for 1:1 future porting:
1. Functions must be **pure** — no mutation of `DICpara` inside callees (MATLAB's pass-by-value hides this bug today but Python will surface it).
2. **Explicit arguments** — camera-specific data (mask, ROI) passed as params, not smuggled through `DICpara`.
3. Default values centralized in `setDICparaDefaults.m`; validation in `validateDICpara.m`.
4. Optional progress callback via `DICpara.progressCallback = @(frac, msg) ...`.
5. Keep field names **string-based where extension is likely**: e.g., `DICpara.RefMode = 'acc' | 'inc'` instead of boolean, so future `FrameSchedule`-style generalization is non-breaking.

---

## Workstreams (8 total, ~12.5h)

### T — Test / CI script (FIRST)

Single MATLAB script `tests/run_pipeline_test.m` that:
- Runs the full pipeline on `examples/Stereo_DIC_Challenge_1.0_S3/` with a fixed preset `DICpara` (no interactive prompts).
- Captures `FinalResult`, `RD_L`, `RD_R`, and per-section timings.
- First run saves `tests/baseline/baseline.mat`.
- Subsequent runs diff against baseline: `max(abs(...))` on disp + strain, relative timing deltas.
- Pass criterion: disp/strain diff < 1e-8, timing regression < 20%.

### E — Git hygiene (~15 min)

Same as v1 design doc:
- Extend `.gitignore` (13 entries for private files)
- `git rm --cached` on 8 already-tracked private files

### C1 — Cleanup (expanded, ~2h)

All of v1's C1 class 1 (8 items) PLUS:
- Remove `DataDrivenOrNot` branching in `TemporalMatch_quadtree_ST1` L30-42 + L201-213 (Q1)
- Remove 2nd-order shape function path: in `main_3D_ALDIC.m` + `StereoMatch_STAQ.m`, drop the `if stereoMatchShapeOrder == 1 ... else ...` and always use order 1 (Q4)
- Remove all adaptive winsize commented code in `TemporalMatch_quadtree_ST1` L370-382 (A2.8)
- Remove the `JY!!!Mask` commented block in `PlotdispQuadtreeMasks3D_acc_ST1` L29-57 (A9.1, user confirmed don't need)
- Remove inc-version mask-out block if any (user confirmed don't need)
- Delete dead code in `funParaInput.m`, `stereoReconstruction_quadtree.m`, `setDICParas_STAQ.m`, `PlaneFit3_Quadtree.m`, etc. (A3-A8)
- **All Chinese comments converted to English**

### C2 — User convenience (adjusted, ~2.5h)

- C2.1 Smart mex recompile (check mexw64 timestamp)
- C2.2 Conditional MinGW setenv with clear fallback warning
- C2.3 Default `outputFilePath = ./results/<timestamp>/`
- C2.4 Section banners with elapsed time
- C2.5 Promote hardcoded overrides into `DICpara` fields:
  - `DICpara.showImgOrNot` (default 0 for batch; was hardcoded 0 in TemporalMatch)
  - `DICpara.UseGlobal` (default true; was hardcoded 1)
  - `DICpara.strain_size` (was interactive prompt)
  - `DICpara.PlotMode ∈ {'none', 'save', 'show'}` (default 'show')
  - `DICpara.forceMexRebuild` (default false)
  - `DICpara.NewFFTSearch` (single source of truth; currently set in BOTH `setDICParas_IncOrNot` and `TemporalMatch_quadtree_ST1`)

### C3 — Architectural (Python-port friendly, ~3h)

- C3.1 **Split plot/compute** (B1). Extract `func_quadtree/computeStrain3D.m`:
  ```matlab
  function [exx, eyy, exy, e1, e2, maxShear, vonMises, dwdx, dwdy] = ...
      computeStrain3D(coefficients)
  ```
  Called by `PlotstrainQuadtreeMasks3D_acc_ST1` and `_inc_ST1`. Enables headless strain computation.

- C3.2 **Right-camera DICpara isolation** (Q5=B). Change `TemporalMatch_quadtree_ST1` signature:
  ```matlab
  function RD = TemporalMatch_quadtree_ST1(DICpara, ...
      fileName, maskSeries, imgSeries, RD, stereoInfo, ...
      cameraLabel, shapeFuncOrder, refMask, roiRange)
  ```
  where `refMask` and `roiRange` are per-camera and passed explicitly. Remove the "fake update" hack at L63-77.

- C3.3 **RBF → scatteredInterpolant everywhere** (B5): `PlotstrainQuadtreeMasks3D_inc_ST1` L43-56.

- C3.4 **Unify F matrix** (A9.3): both acc and inc versions use `F = [1+u_x, u_y, u_z; v_x, 1+v_y, v_z; w_x, w_y, 1+w_z]`. Since u_z/v_z/w_z are 0 under 'Local' coord (which main always uses), output is identical; eliminates latent 'Camera0' bug.

- C3.5 **Fix inc plot position** (Q6=b): change `CurrentFEM.coordinatesFEM + U_2D_L_inc` → `CurrentFEM.coordinatesFEM` in `PlotstrainQuadtreeMasks3D_inc_ST1` L38.

### C4 — Light GUI guard (~1h)

Each existing interactive prompt wrapped in `if isempty(DICpara.X) ... prompt ... end`, so a preset `DICpara` skips all prompts. No centralized dialog system.

### C5 — Independent smoothing controls (~1h)

- Add `DICpara.Smooth2DTimes` (default 0 — matches current disabled state inside ADMM)
- Add `DICpara.Smooth3DTimes` (default 3 — matches current `SmoothTimes < 3` loop in main)
- Replace `while SmoothTimes < 3` in main L307 with `while SmoothTimes < DICpara.Smooth3DTimes`
- Replace `while (SmoothTimes > 2, break)` in `funSmoothDispQuadtree` L86 with `DICpara.Smooth2DTimes` when called from ADMM

### C6 — Defaults + validation separation (~1h)

Two new helpers:
- `func/setDICparaDefaults.m`: returns a DICpara struct with all fields initialized. Called at top of `main_3D_ALDIC.m`.
- `func/validateDICpara.m`: asserts types, ranges, and mutual consistency (e.g., `winsize_min <= winstepsize`, both powers of 2).

Mirror pyALDIC's `dicpara_default()` + `validate_dicpara()`.

---

## Invariants (pre/post refactor)

1. `FinalResult.Coordinates`, `FinalResult.Displacement`, `FinalResult.ResultStrainWorld` bit-identical on the example dataset with default settings.
2. Total wall-clock time ≤ 120% of pre-refactor baseline (no major perf regression).
3. `main_3D_ALDIC.m` with zero preset `DICpara` still runs interactively (all prompts still fire).
4. `main_3D_ALDIC.m` with a fully preset `DICpara` runs batch (no prompts).

## Deliverables

- `tests/run_pipeline_test.m` (new, baseline-capture + diff)
- `tests/baseline/baseline.mat` (user-captured on first run; not checked in)
- `func/setDICparaDefaults.m` (new)
- `func/validateDICpara.m` (new)
- `func/initResultStorage.m` (new, from v1)
- `func/verify_stereo_calibration.m` (new, from v1)
- `func_quadtree/computeStrain3D.m` (new)
- Modified: `main_3D_ALDIC.m` (~250 lines, down from 404)
- Modified: `TemporalMatch_quadtree_ST1.m` (signature change + ~150 line reduction from dead code removal)
- Modified: `PlotstrainQuadtreeMasks3D_{acc,inc}_ST1.m` (unified F matrix, no inline compute, no RBF)
- Modified: `PlotdispQuadtreeMasks3D_{acc,inc}_ST1.m` (remove JY!!!Mask blocks)
- Modified: `PlaneFit3_Quadtree.m` (~100 line reduction, English-only comments)
- Modified: `StereoMatch_STAQ.m`, `stereoReconstruction_quadtree.m`, `setDICParas_STAQ.m`, `setDICParas_IncOrNot.m`, `funParaInput.m`, `funSmoothDispQuadtree.m`, `funSmoothDisp_Quadtree.m` (cleanup)
- Updated `.gitignore`
- Git history: ~25 small commits

## Validation protocol

After each workstream commit:
1. User runs `tests/run_pipeline_test.m` in MATLAB.
2. Script auto-diffs against baseline. Any regression is reported.
3. If pass, proceed to next workstream. If fail, debug before continuing.
