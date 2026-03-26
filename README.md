# CEE 530 Prestressed Concrete — Section Analysis

This folder is the working brain for Project 1 and homework assignments in CEE 530 (ASU Spring 2026), focused on prestressed concrete Double-T beam analysis, Magnel diagrams, and feasibility design charts.

## Purpose

Use this workspace to:

- analyze prestressed Double-T section stresses at transfer and service
- develop the Magnel diagram and feasibility zone for tendon design
- verify stress limits per ACI 318 at critical sections along the span
- generate figures and reports for Project 1 deliverables

## Folder Map

- `01_sources`: ACI codes, lecture slides, reference PDFs (drop zone — never edit files here)
- `02_literature_review`: class notes, theory summaries, derivation handouts
- `03_concept_notes`: Maple worksheets and concept-level derivations for Magnel diagram, feasibility charts
- `04_comparison_tables`: homework assignments by number (HW_1, HW_2, etc.) with input files and output figures
- `05_figures`: final polished figures for reports and submissions
- `06_draft_text`: draft and final written reports
- `07_tasks`: open questions, to-do lists, next actions
- `08_matlab_program_analysis`: MATLAB scripts for beam analysis and feasibility design chart

## Recommended Workflow

1. Drop reference material into `01_sources`.
2. Take class notes in `02_class_notes`.
3. Work out derivations and concept checks in `03_concept_notes`.
4. Run MATLAB analyses from `08_matlab_program_analysis` — output figures go in the `Project1/` subfolder.
5. Move polished figures to `05_figures` for reports.
6. Write up deliverables in `06_reports`.

## Expected Outputs

- Stress diagrams at critical sections (transfer and service)
- Magnel diagram with 4 boundary lines and feasible zone shaded
- Feasibility zone envelope (e_min / e_max vs. x/L) with tendon profile overlay
- Console verification: pass/fail at each section for ACI 318 stress limits

## Current Focus (Spring 2026)

- Project 1: Double-T beam (120" wide, 28" deep), 4 strands (Aps = 0.153 in² each, Pe = 25 kip each)
- Span: 64 ft simply supported
- Check stresses at transfer (self-weight only) and service (SW + SDL + LL)
- Design report due end of semester

## Draft Naming

Use this subject-specific name for report artifacts in this folder:

- `CEE 530 Prestressed Concrete Spring 2026 - Project 1 Report`
