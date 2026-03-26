# Primary Sources — CEE 530 Prestressed Concrete

## NotebookLM

**URL:** https://notebooklm.google.com/notebook/6c3d0565-e74f-4483-9ace-080c39fe89cf
**Library ID:** `cee-530-prestressed-concrete`
**Sources confirmed in notebook:** ACI 318-19, Naaman 2nd ed, Wight RC 8th ed

---

## Sources Confirmed in NotebookLM

### 1. ACI 318-19 — Building Code Requirements for Structural Concrete [IN NOTEBOOKLM]
**Full Citation:** ACI Committee 318. *Building Code Requirements for Structural Concrete (ACI 318-19) and Commentary (ACI 318R-19)*. American Concrete Institute, 2019.
**Why it matters:** Governs all allowable stress limits used in the analysis:
- Transfer: f_ci_allow = 0.60 f'ci (compression), −3√f'ci ksi (tension)
- Service: f_cs_allow = 0.45 f'c (compression), −12√f'c psi (tension, Class C)
**Local copy:** Not yet in `01_sources/` — query via NotebookLM

---

### 2. Naaman — Prestressed Concrete Analysis and Design: Fundamentals (2nd ed) [IN NOTEBOOKLM]
**Full Citation:** Naaman, A.E. *Prestressed Concrete Analysis and Design: Fundamentals* (2nd ed). Techno Press 3000, 2004.
**Why it matters:** Primary course textbook. Covers working stress design, Magnel diagram derivation, tendon layout, stress analysis at transfer and service, prestress losses, feasibility zone.
**Key sections:** Flexural analysis, tendon design, Magnel diagram, stress limits
**Local copy:** `01_sources/Naaman_prestressedconcrete.pdf`

---

### 3. Wight — Reinforced Concrete: Mechanics and Design (8th ed) [IN NOTEBOOKLM]
**Full Citation:** Wight, J.K. *Reinforced Concrete: Mechanics and Design* (8th ed). Pearson, 2016.
**Why it matters:** Supplementary reference for ACI 318-19 code interpretation, shear design, and general RC design context.
**Local copy:** Not in `01_sources/` — query via NotebookLM

---

## Course Slides

### CEE 530 Spring 2024 Lecture Slides [LOCAL PDF]
**File:** `01_sources/Pretressed Concrete (CEE 530) spring 2024.pdf`
**Why it matters:** ASU course-specific content, problem setup, section properties for Double-T beam, Project 1 background.
**When to use:** For course-specific notation or problem parameters that may differ from textbook.

---

## Tracking Format

When a new source is obtained, add:
```
**Full Citation:** Author(s), Year. "Title." Publisher.
**Local copy:** 01_sources/filename.pdf
**Key equations used:** list relevant equations
**Notes:** any important clarifications
```
