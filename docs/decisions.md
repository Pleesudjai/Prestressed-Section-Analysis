# Decision Log

## 2026-03-16 — Project Setup & Folder Structure
**What:** Created 8-subfolder structure mirroring CEE 790 pattern. Moved all MATLAB scripts to `08_matlab_program_analysis/`. Renamed `Projec1/` → `Project1/`.
**Why:** Predictable navigation across all CEE courses. Keeps root clean.
**Files:** CLAUDE.md, README.md, all folders `01_sources/` through `08_matlab_program_analysis/`
**Next:** Complete Project 1 stress analysis and Magnel diagram.

## 2026-03-16 — Sign Convention Decision
**What:** Adopted compression = positive throughout all scripts and documentation.
**Why:** Consistent with ACI 318 prestress notation (Naaman textbook convention).
**Files:** All `.m` files in `08_matlab_program_analysis/`
**Next:** Never flip — all future scripts must follow this convention.

## 2026-03-16 — Coordinate Origin Decision
**What:** Set y = 0 at bottom of stems (not top of flange).
**Why:** Makes shoelace formula consistent, eccentricity `e = yc − y_tendon` is always positive-below-centroid.
**Files:** All input files, `analyzePrestressedBeam.m`, `feasibilityDesignChart.m`
**Next:** Ignore old comments in input files that reference y from top.

## 2026-03-23 — Feasibility Design Chart Created
**What:** Built standalone `feasibilityDesignChart.m` with Magnel diagram (Figure 1) and feasibility zone envelope (Figure 2). Implements all 4 main Magnel lines + 2 additional checks. Dead-load-only case (SW + SDL, no LL).
**Why:** Project 1 deliverable — verify that design point (F=100 kip, e=13.44 in) falls within feasible zone.
**Files:** `feasibilityDesignChart.m` (root, to be moved to `08_matlab_program_analysis/`)
**Next:** Run in MATLAB, verify plots, move to `08_matlab_program_analysis/`, save figures to `Project1/`.

## 2026-03-23 — WISC Framework Setup
**What:** Installed WISC context engineering scaffold: 5 slash commands (`/prime`, `/plan-feature`, `/execute`, `/handoff`, `/commit`), 3 domain rules (`matlab.md`, `prestress-theory.md`, `report.md`), decision log, handoff doc.
**Why:** Prevent context rot in long Claude Code sessions. Keep sessions productive.
**Files:** `.claude/commands/`, `.claude/rules/`, `docs/decisions.md`, `docs/handoff.md`
**Next:** Use `/prime` to start future sessions, `/handoff` to end them.

## 2026-03-25 — Vci Two-Tier Bounds (Floor + Cap)
**What:** Changed the 5.0 coefficient from a lower bound (floor) to an **upper bound (cap)** on Vci. Added 1.7 as the ACI 318-19 hard floor. MATLAB logic: `Vci = max(Vci, 1.7*...); Vci = min(Vci, 5.0*...);`
**Why:** User corrected the interpretation: ACI 318-19 Eq. 22.5.8.3.1a uses 1.7 as the minimum; the 5.0 (simplified-method upper cap from Sec. 22.5.8.1) is used by CEE 530 as a ceiling to prevent inflated Vci near supports.
**Files:** `shearDesign.m` (lines 270-276), `create_shear_report_v2.js` (Sec 2.1, 2.2, 3.6, Part 5)
**Impact:** Sections 2-4 changed from Vcw-governed to Vci-governed. Section 2 stirrup spacing tightened from 21" to 13".
**Next:** Review v3 report formatting; ultimate design report.

## 2026-03-25 — Detailed Calculations for All 4 Sections
**What:** Added Parts 3B, 3C, 3D to the shear report with full hand-calc detail for Sections 2 (x=10 ft), 3 (x=20 ft), and 4 (x=32 ft, midspan). Previously only Section 1 (x=1.5 ft) had detailed calcs.
**Why:** User requested full transparency for all summary table values.
**Files:** `create_shear_report_v2.js` (lines 540-737), `ShearDesign_Report_v3.docx`
**Next:** Verify formatting in Word.

## 2026-03-26 — Comprehensive Git Commit
**What:** First full commit of organized project: folder reorganization, core analysis code updates, Project 2 shear design (with Vci two-tier bounds), Project 3 feasibility, defineDesignStages.m, feasibilityDesignChart.m, CLAUDE.md, docs, and .claude rules.
**Why:** All prior work was staged but uncommitted. Consolidating into one clean commit.
**Files:** All `08_matlab_program_analysis/` scripts, `project_Project2/`, `project_Project3 Feasibility/`, `04_comparison_tables/`, `.claude/`, `docs/`, `CLAUDE.md`
**Next:** Ultimate design report (`ultimateDesign.m`).
