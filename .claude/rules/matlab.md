# Domain Rules — MATLAB Scripts

Load this file when editing or creating .m files in `08_matlab_program_analysis/`.

## Working Directory
All scripts assume `cd` to `08_matlab_program_analysis/`. Relative paths to `Project1/` will fail from root.

## Sign Convention (ABSOLUTE — never violate)
- **Compression = POSITIVE (+)**
- **Tension = NEGATIVE (−)**
- This applies to: stresses, allowable limits, prestress force effects, Magnel line formulations

## Coordinate System
- Origin: `y = 0` at **bottom of stems**
- y-axis: positive **upward**
- Eccentricity: `e = yc − y_tendon` → positive when tendon is **below** centroid

## File Conventions
- Input files: `inputPrestressedBeam_[ProjectName].m` — returns 6 structs
- Analysis: `analyzePrestressedBeam.m` — core engine, do not modify lightly
- Plotting: `plotPrestressedBeamResults.m`, `plotSection.m`, `plotSectionStressStrain.m`
- Standalone scripts (like `feasibilityDesignChart.m`): hardcode inputs, no external deps
- Figures saved to project subfolders: `Project1/`, `project_Project2/`, etc.

## Function Signature Pattern
```matlab
function [beam, section, materials, prestress, reinforcement, loads] = inputPrestressedBeam_XXX()
function [results] = analyzePrestressedBeam(beam, section, materials, prestress, reinforcement, loads, stage)
function stages = defineDesignStages(materials, prestress)
```

## Section Properties
Always use **shoelace formula** via `computeSectionProps(vertices)` or `createPolygonalSection(vertices)`.
Vertices must be **counter-clockwise**.

## Prestress Rules
- At **transfer**: use full `F` (η = 1.0), moment = self-weight only
- At **service**: use `η·F` (η = 0.85 for 15% losses), moment = SW + SDL [+ LL]
- Tendon force: `Pe = Aps × fpi` per strand
- Never apply SDL at transfer for precast members

## Magnel Diagram Pitfall
Check denominator signs for all 4 lines:
- `dI = fci_allow − Mi/St`  → if negative, Line I flips from lower to upper bound
- `dII = Mi/Sb − fti_allow` → should be positive for reasonable spans
- `dIII = MT/Sb − fts_allow`
- `dIV = fcs_allow − MT/St`

## Code Style
- Helper functions go at bottom of file as local functions (MATLAB convention)
- Use `fprintf` for console output, always include units
- Figure names: set `'Name'` property for `matlab-save-figures` skill compatibility
