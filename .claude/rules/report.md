# Domain Rules — Reports & Deliverables

Load this file when writing content for `06_draft_text/` or preparing homework submissions.

## Prerequisites
- `03_concept_notes/` must be complete for the topic before writing in `06_draft_text/`
- Figures must exist in `05_figures/` or `08_matlab_program_analysis/Project*/` before referencing

## Report Naming
`CEE 530 Prestressed Concrete Spring 2026 - [Assignment] Report`

## Structure (typical project report)
1. Problem statement & given data
2. Section properties (table: Ac, Ic, yc, yt, yb, St, Sb)
3. Material properties & allowable stresses (table)
4. Load calculations (w_sw, w_SDL, w_LL, moments at critical sections)
5. Tendon layout description + profile figure
6. Stress analysis results (table: f_top, f_bot at transfer & service, pass/fail)
7. Magnel diagram figure + interpretation
8. Feasibility zone figure + interpretation
9. Conclusions

## Figure Standards
- All figures must have labeled axes with units
- Title format: descriptive, no MATLAB default titles
- Font size ≥ 11 pt for readability
- Include legend when multiple curves are shown
- Export as both .png (for report) and .fig (for editing)

## Sign Convention in Text
When writing prose: state stresses with sign and label.
Example: "The bottom fiber stress at transfer is −0.15 ksi (tension), which is within the allowable tension limit of −0.208 ksi."

## Units
- Stress: ksi
- Force: kip
- Length: in (convert to ft only for span description)
- Moment: kip-in (or kip-ft with explicit conversion)
- Area: in²
- Moment of inertia: in⁴
