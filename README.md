# 3D-Stereo-ALDIC

![Demo video](assets/images/stereo-ALDIC_demo.gif)

## Overview

3D-Stereo-ALDIC (3D Stereo Adaptive Mesh Augmented Lagrangian Digital Image Correlation) is a MATLAB-based implementation for full-field 3D deformation measurement using stereo vision and adaptive quadtree mesh refinement. It provides accurate and efficient 3D displacement and strain measurements for experimental mechanics applications.

## Features

- **Two entry points**: a point-and-click GUI (`run_stereo_aldic_gui`) and a traditional script (`main_3D_ALDIC`).
- Augmented Lagrangian DIC with kinematic compatibility enforcement.
- Adaptive quadtree mesh refinement for sharp features and complex geometries.
- Both **accumulative** (small-deformation) and **incremental** (large-deformation) tracking.
- Calibration format support: MATLAB Stereo Calibrator, MatchID, MCC, DICe XML, OpenCorr.
- Built-in parallel ICGN using the Parallel Computing Toolbox.
- Post-pipeline visualization with frame slider, variable dropdown, color-range sliders, and CSV export.

## Requirements

- MATLAB R2024a or newer (tested on R2024b).
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Computer Vision Toolbox
- Parallel Computing Toolbox (optional, for multi-core ICGN)
- A C/C++ compiler for MEX compilation. On Windows we've used MinGW-w64 via [TDM-GCC](https://jmeubank.github.io/tdm-gcc/); set `MW_MINGW64_LOC` if the compiler lives at a non-default path.

---

## Quick start

### Option A — GUI (recommended for first-time users)

```matlab
run_stereo_aldic_gui
```

The GUI opens a tabbed window with three tabs:

1. **Configuration** — pick input folders and parameters.
2. **Log** — pipeline output (stage-by-stage progress).
3. **Visualize** — post-run result explorer (enabled after the pipeline finishes).

**Configuration tab, section by section:**

- **1. Data**
  - `Left images folder` / `Right images folder` — two folders, one per camera. Images within each folder must be named in sorted frame order.
  - `Left mask folder` / `Right mask folder` — binary masks (white = ROI, black = ignore). **Accumulative mode** only needs the first mask pair; **incremental mode** requires one mask pair per frame.
  - `Calibration file` + `Calibration format` — pick the stereo calibration file exported from your calibration tool; select the matching format.
  - `DIC mode` — Accumulative (reference = frame 0) or Incremental (reference = previous frame).
  - `View calib` — inspect the parsed intrinsics / extrinsics / distortion.
  - `Draw ROI` — optional polygon ROI drawer (accumulative mode only).

- **2. DIC parameters** — subset size, step size, quadtree min size, FFT search radii (independent for stereo matching and temporal matching), parallel workers, ALDIC global-step toggle, and temporal initial-guess strategy (`Per-frame FFT` or `Reuse previous frame` — the latter only applies in accumulative mode).

- **3. Strain & smoothing** — strain window size (in points), 2D/3D smoothing passes, optional transform to specimen coordinate frame.

- **4. Output** — plot mode (default `none` = skip per-frame legacy plots for speed — use the Visualize tab instead), figure file format, output folder.

Click **Run ▶**. The window switches to the Log tab while running, then jumps to the Visualize tab when the pipeline completes. Everything is saved to `FinalResult.mat` inside the output folder you chose.

**Visualize tab:**

- Variable dropdown: `U`, `V`, `W`, `|U|`, `epsilon_xx`, `epsilon_yy`, `epsilon_xy`, `epsilon_1`, `epsilon_2`, `maxshear`, `vonMises`.
- Frame slider (jumps between loaded frames).
- Color-range sliders (min / max), plus `Auto color` to reset.
- Background choice: current frame or first frame.
- `Export CSV…` — write the currently-displayed variable × frame as a 3-column CSV (`X_image, Y_image, value`).

### Option B — Script mode (batch / scripted runs)

```matlab
main_3D_ALDIC
```

The script walks through the same pipeline but collects parameters via the original sequence of command-line prompts and `uigetfile` dialogs. Useful for reproducible scripted workflows where you preset a `DICpara` struct in your own driver (each prompt is guarded by `isempty(DICpara.fieldname)`, so a fully-populated struct runs end-to-end non-interactively).

---

## Example dataset

The repo ships with a compact demo dataset under
`examples/Stereo_DIC_Challenge_1.0_S3/` (3 frame pairs, DICe calibration, masks). Point the GUI inputs at:

```
Left images folder:  examples/Stereo_DIC_Challenge_1.0_S3/Images_Stereo_Sample3_images/Left
Right images folder: examples/Stereo_DIC_Challenge_1.0_S3/Images_Stereo_Sample3_images/Right
Left mask folder:    examples/Stereo_DIC_Challenge_1.0_S3/Images_Stereo_Sample3_maskfiles/Left
Right mask folder:   examples/Stereo_DIC_Challenge_1.0_S3/Images_Stereo_Sample3_maskfiles/Right
Calibration file:    examples/Stereo_DIC_Challenge_1.0_S3/calibration_DICe.xml
Calibration format:  DICe XML (3)
```

Reference result figures from a known-good run are included under `Images_Stereo_Sample3_results/` for sanity checking.

---

## Input data format

- **Images**: any common 8-bit or 16-bit format (tif, png, jpg, bmp, jp2). Left and right go into separate folders.
- **Masks**: same formats, non-zero = inside ROI. Left and right go into separate folders.
- **Calibration**: one of MATLAB CV (`.mat`), MatchID (`.caldat`), MCC (`.mat`), DICe (`.xml`), or OpenCorr (`.csv`).

## Output

`FinalResult.mat` is saved to the configured output folder (default `./results/<timestamp>/`) and contains:

- `FinalResult` — 3D coordinates, displacements, and strain per frame.
- `DICpara` — full configuration used for the run.
- `RD_L`, `RD_R` — per-camera meshes + 2D displacements.
- `StereoInfo` — stereo match + calibration.
- `coefficientsPerFrame`, `strainPerFrame` — raw plane-fit outputs + derived strain components.
- `fileNameLeft`, `fileNameRight` — image path tables for Visualize-tab background reloading.

---

## Citation

If you use this code in your research, please cite:

1. Tong, Z. et al. *3D Stereo Adaptive Mesh Augmented Lagrangian Digital Image Correlation.* Exp. Mech. 2025. [https://doi.org/10.1007/s11340-025-01225-7](https://doi.org/10.1007/s11340-025-01225-7)
2. Tong, Z. et al. *Machine Learning-Aided Spatial Adaptation for Improved Digital Image Correlation Analysis of Complex Geometries.* [https://doi.org/10.21203/rs.3.rs-5566473/v1](https://doi.org/10.21203/rs.3.rs-5566473/v1)

**Previous related work:**

3. Yang, J. and Bhattacharya, K. *Augmented Lagrangian Digital Image Correlation.* Exp. Mech. 59: 187, 2018. [https://doi.org/10.1007/s11340-018-00457-0](https://doi.org/10.1007/s11340-018-00457-0)
4. Yang, J. *2D_ALDIC.* [https://github.com/jyang526843/2D_ALDIC](https://github.com/jyang526843/2D_ALDIC)
5. Yang, J. *FEM-based Global DIC.* [https://github.com/jyang526843/2D_FE_Global_DIC](https://github.com/jyang526843/2D_FE_Global_DIC)
6. Yang, J. *ALDVC code.* [https://github.com/FranckLab/ALDVC](https://github.com/FranckLab/ALDVC)

---

## Contributing

Pull requests and issue reports are welcome. The repo includes a regression test (`tests/run_pipeline_test.m`) that captures a baseline on the demo dataset and diffs future runs against it — please ensure it still passes before sending a PR.

## License

MIT License.

## Authors

- Zixiang (Zach) Tong (<zachtong@utexas.edu>) @ UT-Austin
- Jin Yang (<jin.yang@austin.utexas.edu>) @ UT-Austin
