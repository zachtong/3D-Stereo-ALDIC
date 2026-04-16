# Stereo DIC Challenge 2.1 Bespoke — Submission README

## Contents
- `VSG_23/`, `VSG_43/`, `VSG_63/` — 18 CSVs each (54 total)
- Frame list: `0000` (reference), `0001–0005` (noise baseline), `0200, 0600, 1000, 1400, 1800, 2200, 2600` (loading), `2620, 2640, 2660, 2680, 2700` (failure)

## Processing parameters

| Parameter | Value |
|-----------|-------|
| Subset (winsize) | 24 px |
| DIC step (coarse) | 8 px |
| Shape function | Quadratic |
| Correlation | ZNSSD + Inverse-Compositional Gauss-Newton |
| Global regularization | ALDIC ADMM |
| Dense grid step | 1 px (via scatteredInterpolant linear) |
| Strain tensor | Green-Lagrangian |
| Strain units | mm/mm |

## VSG labels vs. actual

| Label | Actual VSG | strain_size |
|-------|-----------:|------------:|
| 23    | **26 px**  | 2           |
| 43    | 43 px      | 19          |
| 63    | 63 px      | 39          |

**Note:** `VSG_23` is labeled 23 per challenge requirement, but the actual virtual strain gauge size is 26 px because `VSG ≥ winsize = 24` is a mathematical lower bound for this subset size. Strain at the labeled-23 size is the minimum achievable given winsize=24.

## Coordinate system

Specimen coordinate system established from 3 base points (see `PROJECT_REQUIREMENTS.md` §3). All X, Y, Z, U, V, W, and strain values are in this specimen frame.
