# Design: Main Script Refactor + Git Hygiene

**Date:** 2026-04-22
**Author:** Zach Tong
**Status:** Approved for implementation

## Goal

Clean up and refactor `main_3D_ALDIC.m` for improved readability, user convenience, and separation from private/research work. No new features, no behavior changes in the core algorithm.

## Scope

Two workstreams, executed in order:

1. **E — Git hygiene**: Prevent private research files from being committed.
2. **C — Main script refactor**: Cleanup (class 1) + user convenience (class 2) of `main_3D_ALDIC.m`.

Explicitly **out of scope**: performance optimization (parfor, checkpointing), config file parameterization, stage modularization, user-facing documentation.

---

## Workstream E: Git Hygiene

### E.1 Extend `.gitignore`

Add the following entries:

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

### E.2 Un-track already-committed private files

Files currently tracked by git but considered private:
- `run_coarse_test_noise.m`
- `run_export_csv.m`
- `validate_csv_output.m`
- `challenge_scripts/export_csv_pipeline.m`
- `challenge_scripts/interpolate_to_dense_grid_from_mask.m`
- `challenge_scripts/compute_strain_dense_grid_lowmem.m`
- `challenge_scripts/export_challenge_csv.m`
- `results/csv_output/README.md`

Action: `git rm --cached <file>` on each, then commit the removal. Files remain on disk.

---

## Workstream C: `main_3D_ALDIC.m` Refactor

### C.1 Cleanup (no behavior change)

| # | Location | Issue | Fix |
|---|----------|-------|-----|
| 1.1 | L124, L230, L237, L241, L244, L252, L257, L300-303, L329-330, L355-357, L399-401 | Commented-out dead code | Remove all. Keep the hardcoded default that is active. |
| 1.2 | L211 | Typo variable `shapeFuncOalized_rder` (unused) | Delete line; use the existing `shapeFuncOrder`. |
| 1.3 | L225, L396 | Section labeled "Section 8" but is actually Section 6 | Renumber to match header workflow. |
| 1.4 | L267-268 | Hardcoded Challenge 2.1 Bespoke base points `[1024+1, 1224+1; 1509.69+1, 1219.38+1; 1020.10+1, 470.48+1]` | Remove. Rely on `getBasePoints` interactive selection only. |
| 1.5 | L128, L267 | "Zach comments" / "Zach 20260127" personal annotations | Remove. |
| 1.6 | L91-111 | ~20 lines of duplicated init for `RD_L` and `RD_R` | Extract helper: `RD = initResultStorage(numFrames, isIncremental)`. Call twice. |
| 1.7 | L293-294 | Hardcoded path separator `'\'` | Replace with `fullfile(fileNameLeft{2,ImgSeqNum}, fileNameLeft{1,ImgSeqNum})`. |
| 1.8 | L168-203 | 3D reconstruction "Test" block embedded in production flow | Move to separate `verify_stereo_calibration.m` helper, gated by `DICpara.verifyStereoReconstruction` flag (default `false`). |

### C.2 User convenience

| # | Current behavior | New behavior |
|---|------------------|--------------|
| 2.1 | `mex -O ba_interp2_spline.cpp` always re-runs (5-15s) | Check if `ba_interp2_spline.mexw64` exists and is newer than `ba_interp2_spline.cpp`; skip if so. Still allow force rebuild via `DICpara.forceMexRebuild`. |
| 2.2 | `strain_size = input(prompt)` blocks execution mid-run | Read from `DICpara.strain_size` if set; only prompt if missing. Set at top of script or via config struct. |
| 2.3 | `setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')` hardcoded | Only setenv if the env var is unset AND the hardcoded path exists. Otherwise emit a clear warning. |
| 2.4 | `DICpara.outputFilePath = []` | Default to `./results/<yyyymmdd_HHMMSS>/`. Create if missing. User can override by setting before this line. |
| 2.5 | Plotting always runs every frame (close all + create figures) | Add `DICpara.PlotMode ∈ {'none', 'save', 'show'}`. Default `'show'`. Skip all `Plot*` calls when `'none'`; skip display when `'save'`. |
| 2.6 | No global progress/timing | Each section starts with a uniform banner: `fprintf('==== Section X/N: <name> | elapsed: Y.Ys ====\n')`. Use `tic/toc` around a `session_start` tic to track. |

### C.3 New helper functions

Two small helpers to be added under `func/`:

- `func/initResultStorage.m` — consolidates the RD_L/RD_R init duplication (for 1.6).
- `func/verify_stereo_calibration.m` — optional debug plot for stereo reconstruction (for 1.8).

---

## Invariants (must hold before and after)

1. **Algorithm outputs unchanged**: For the same inputs, `FinalResult.Coordinates`, `FinalResult.Displacement`, and `FinalResult.ResultStrainWorld` must be bit-identical (up to floating-point determinism).
2. **Function signatures unchanged**: No modifications to `func/*.m` or `func_quadtree/*.m` beyond the new `initResultStorage.m` helper.
3. **Interactive prompts remain**: The user can still run `main_3D_ALDIC` with no preset config and hit every GUI dialog exactly as before (2.2 only **enables skipping** by presetting).

## Validation

A small benchmark: run the refactored script against the existing `examples/Stereo_DIC_Challenge_1.0_S3/` example, compare `FinalResult.Displacement` and `FinalResult.ResultStrainWorld` against a reference `.mat` saved from the pre-refactor code.

## Deliverables

- Modified `main_3D_ALDIC.m` (expected: ~250 lines, down from 404)
- New `func/initResultStorage.m`
- New `func/verify_stereo_calibration.m`
- Updated `.gitignore`
- Git history showing `git rm --cached` for the 8 private files
- One validation run on `Stereo_DIC_Challenge_1.0_S3` confirming algorithm output is unchanged
