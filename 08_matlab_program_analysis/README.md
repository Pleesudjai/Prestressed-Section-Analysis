# MATLAB Program Analysis

All MATLAB scripts for prestressed beam stress analysis and feasibility design charts.

---

## Program Index

| Script | Purpose | Output Folder |
|--------|---------|---------------|
| `main_PrestressedBeamAnalysis.m` | Top-level driver — runs full Project 1 analysis | `Project1/` |
| `inputPrestressedBeam_Project1.m` | Input: full service loads (SW + SDL + LL) | — |
| `inputPrestressedBeam_Project1_DeadOnly.m` | Input: dead load only (SW + SDL) | — |
| `feasibilityDesignChart.m` | Standalone Magnel diagram + feasibility zone | `Project1/` |
| `analyzePrestressedBeam.m` | Core engine: stress + moment along beam | — |
| `plotSection.m` | Cross-section geometry plotter | — |
| `plotPrestressedBeamResults.m` | Results plotter (stress diagrams) | — |
| `plotSectionStressStrain.m` | Stress/strain distribution plotter | — |
| `createPolygonalSection.m` | Shoelace formula section property utility | — |
| `inputPrestressedBeam_HW1.m` | Moved to `../04_comparison_tables/HW_1/` | — |

## Output Subfolder

- `Project1/` — figures saved from Project 1 runs (`.fig` and `.png`)

## How to Run

1. Open MATLAB, set working directory to `08_matlab_program_analysis/`
2. Run `main_PrestressedBeamAnalysis.m` for full stress analysis
3. Run `feasibilityDesignChart.m` for Magnel diagram and feasibility zone
4. Figures auto-save to `Project1/`

## Sign Convention

- **Compression = positive**, tension = negative (throughout all scripts)
- **y = 0** at bottom of stems, positive upward
- **e = yc − y_tendon** → positive when tendon is below centroid
