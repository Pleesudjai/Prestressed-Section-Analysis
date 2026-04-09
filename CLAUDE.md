# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## ▶ RESUME HERE

- **Current focus:** End block design (`endBlockDesign.m`) implemented and tested
- **Status:** Deep beam method (0.8h, Prof. Fafitis) active; Gergely-Sozen (3h/4, Naaman) available as toggle. 6 figures generated. Folder structure flattened — all .m files now in `08_matlab_program_analysis/`. Uncommitted changes in endBlockDesign.m (force plot, stress plot cleanup, deep beam switch).
- **Next action:** Commit latest endBlockDesign.m changes; finish HW3 end block submission; run full `main_PrestressedBeamAnalysis` to verify pipeline after restructure
- **Status updated:** 2026-04-08

---

## ⛔ DO NOT

- Do NOT update CLAUDE.md manually with the Write tool — always invoke the `make-claude-md` skill first
  - **Trigger phrases that MUST invoke the skill:** "update claude.md", "save progress", "save session", "save memory", "before I go", "don't forget", "remember this"
- Do NOT write to `06_draft_text/` before `03_concept_notes/` is complete for that topic
- Do NOT place new MATLAB scripts in the root — all scripts go in `08_matlab_program_analysis/`
- Do NOT use the old `Projec1/` folder name — it has been renamed to `08_matlab_program_analysis/Project1/`
- Do NOT mix sign convention: **compression = positive, tension = negative** throughout all scripts
- Do NOT apply prestress losses η at transfer — use full F; η applies only at service
- Do NOT use y measured from top of section — origin is at **bottom of stems, y = 0**

---

## Context Management (WISC Framework)

This project uses the **WISC** framework (Write, Isolate, Select, Compress) to prevent context rot.

| Strategy | How it works here |
|----------|-------------------|
| **Write** | Externalize memory → `docs/decisions.md`, `docs/handoff.md`, git commits |
| **Isolate** | Sub-agents for research only — return 3-5 bullet summaries, never raw dumps |
| **Select** | Load only what's needed: CLAUDE.md always, `.claude/rules/` only when touching that domain |
| **Compress** | Run `/handoff` + start fresh session when context exceeds ~50 messages |

### Slash Commands
| Command | When to use |
|---------|-------------|
| `/prime` | Start of every session — reads CLAUDE.md + decisions + handoff + git state |
| `/plan-feature` | Before any task >1 hour — creates spec in `specs/` |
| `/execute` | Implementation mode — load spec only, no planning context |
| `/handoff` | End of session — writes `docs/handoff.md`, updates RESUME HERE |
| `/commit` | After completing work — logs to `docs/decisions.md` + git commit |

### Domain Rules (load selectively, never all at once)
| Rule file | Load when... |
|-----------|-------------|
| `.claude/rules/matlab.md` | Editing or creating .m files |
| `.claude/rules/prestress-theory.md` | Working on stress analysis, Magnel, ACI checks |
| `.claude/rules/report.md` | Writing content for `06_draft_text/` |

### Key WISC Files
| File | Purpose |
|------|---------|
| `docs/decisions.md` | Running log of what was built, why, and what's next |
| `docs/handoff.md` | Latest session summary for next Claude instance |
| `specs/` | Feature specs from `/plan-feature` |

---

## Project Context

CEE 530 Prestressed Concrete, ASU PhD coursework, Spring 2026.
Student: Chidchanok Pleesudjai (Fen). Focus: Double-T precast prestressed beam analysis.

---

## Code Architecture

### Call Graph (main analysis pipeline)

```
main_PrestressedBeamAnalysis.m  (ROOT DRIVER — script)
│
├─→ inputPrestressedBeam_Project1()
│   └─ Returns: [beam, section, materials, prestress, reinforcement, loads]
│   └─ Embeds: calculateSectionProperties(), processTendonProfiles(), displayInputSummary()
│
├─→ defineDesignStages(materials, prestress)
│   └─ Returns: cell array of stage structs {Transfer, Service_Sustained, Service_Total}
│   └─ Each stage has: .eta, .include_sw/.include_sdl/.include_ll, .f_allow_compr/.f_allow_tens
│
└─→ FOR EACH ACTIVE STAGE:
    ├─→ analyzePrestressedBeam(beam, section, materials, prestress, reinforcement, loads, stage)
    │   ├─ calculatePrestressEffects()  → P, e_eff, M_prestress, V_prestress
    │   ├─ calculateInternalForces()    → V, M, reactions, w
    │   ├─ calculateStresses()          → f_top, f_bot, allowables
    │   └─ calculateCapacity()          → capacity struct (non-Transfer only)
    │
    ├─→ plotPrestressedBeamResults(results, options)   — 4×2 subplot: beam, N, V, M, stresses
    ├─→ plotSection(section, prestress, reinforcement, materials, opts)   — cross-section figure
    ├─→ plotSectionStressStrain(results, x_location, opts)   — vertical stress distribution at x
    └─→ saveStageFigures(data, stage_name)   — .png + .fig to Project1/
```

### Standalone Scripts (all in `08_matlab_program_analysis/`, run from there)

| Script | Purpose | Dependencies |
|--------|---------|-------------|
| `shearDesign.m` | ACI 318-19 shear design (Vci/Vcw), output to `project_Project2/output/` | `inputData()` |
| `ultimateDesign.m` | Nominal moment capacity (ACI 318-19), output to `project_Project2/output/` | `inputData()` |
| `runEndBlockDesign.m` | End block anchorage zone design (Deep Beam / Gergely-Sozen), output to `project_Project2/output/EndBlock/` | `inputData()`, `endBlockDesign()`, `plotSection()`, `plotSectionStressStrain()` |
| `feasibilityDesignChart.m` | Magnel diagram (1/F vs e) + feasibility zone (e vs x/L). Hardcoded inputs. | None (self-contained) |
| `example_sections.m` | Demo section definitions | `createPolygonalSection()` |
| `quickInputPE.m` | One-liner input: `quickInputPE(P, e, L, vertices)` | Same pipeline as full input |

### Data Flow

All analysis functions pass a **`results` struct** with fields:
`.x`, `.L`, `.N`, `.V`, `.M`, `.P`, `.e`, `.f_top_total`, `.f_bot_total`, `.stresses`, `.reactions`, `.stage_name`

Input functions return **6 structs**: `beam`, `section`, `materials`, `prestress`, `reinforcement`, `loads`.

### Project Subfolders

| Folder | Status |
|--------|--------|
| `project_Project2/` | Output + reports only (all .m files moved to `08_matlab_program_analysis/`) |
| `project_Project3 Feasibility/` | Multi-section feasibility (`feasibility4Sections.m`, `runAllCases.m`) |
| `project_Project4 Ultimate/` | Placeholder (empty) |

---

## MATLAB Rules

- All scripts must be run with working directory set to `08_matlab_program_analysis/`
- Each analysis saves figures to its own output subfolder (e.g., `Project1/`)
- Sign convention: **compression = positive, tension = negative** — never flip
- At transfer: use prestress force `F` (full); at service: use `η·F` (with losses)
- Eccentricity: `e = yc − y_tendon`, positive when tendon is **below** centroid
- Moment (simply supported): `M(x) = w·x·(L−x)/2`
- Section properties computed via **shoelace formula** (`computeSectionProps()` or `createPolygonalSection()`)

---

## Feasibility Design Chart — Algorithm

### Fiber Stress Equations (compression = +)

```
AT TRANSFER (force F, moment Mi = self-weight only):
  f_top = F/Ac  −  F·e/St  +  Mi/St
  f_bot = F/Ac  +  F·e/Sb  −  Mi/Sb

AT SERVICE (force η·F, moment MT = SW + SDL [+ LL]):
  f_top = η·F/Ac  −  η·F·e/St  +  MT/St
  f_bot = η·F/Ac  +  η·F·e/Sb  −  MT/Sb
```

### Four Magnel Lines (1/F as function of e — straight lines)

| Line | Condition | Formula | Bound |
|------|-----------|---------|-------|
| **I**   | Top @ transfer ≤ fci | `1/F = (1/Ac − e/St) / (fci − Mi/St)` | lower |
| **II**  | Bot @ transfer ≥ −fti | `1/F = (1/Ac + e/Sb) / (Mi/Sb − fti)` | upper |
| **III** | Bot @ service ≥ −fts | `1/F = η·(1/Ac + e/Sb) / (MT/Sb − fts)` | upper |
| **IV**  | Top @ service ≤ fcs | `1/F = η·(1/Ac − e/St) / (fcs − MT/St)` | lower |

Feasible zone: `max(I, IV) ≤ 1/F ≤ min(II, III)`, with `e ≤ yb − cover`.

### Eccentricity Bounds (for fixed F, varying x along beam)

| Condition | Bound | Formula |
|-----------|-------|---------|
| Top compr @ transfer | e ≥ | `St/Ac + (Mi − fci·St)/F` |
| Top tens @ transfer | e ≤ | `St/Ac + (Mi + fti·St)/F` |
| Bot tens @ transfer | e ≥ | `(Mi − fti·Sb)/F − Sb/Ac` |
| Bot compr @ transfer | e ≤ | `(Mi + fci·Sb)/F − Sb/Ac` |
| Top compr @ service | e ≥ | `St/Ac + (MT − fcs·St)/(η·F)` |
| Top tens @ service | e ≤ | `St/Ac + (MT + fts·St)/(η·F)` |
| **Bot tens @ service** | **e ≥** | `(MT − fts·Sb)/(η·F) − Sb/Ac` ← usually governs |
| Bot compr @ service | e ≤ | `(MT + fcs·Sb)/(η·F) − Sb/Ac` |

`e_min(x) = max(all lower bounds)`, `e_max(x) = min(all upper bounds, e_geo_max)`

---

## Coordinate System

- **Origin**: bottom of stems, `y = 0`
- **y-axis**: positive **upward** (bottom = 0 in, flange soffit = 26 in, top = 28 in)
- **Eccentricity**: `e = yc − y_tendon` → positive when tendon is **below** centroid

### Double-T Section (Project 1)

```
 _______________________________________________   y = 28 in (top)
|               TOP FLANGE  120" wide          |
|_______________________________________________|   y = 26 in (flange soffit)
       |  left stem  |       | right stem |
       |_____________|       |____________|        y = 0 in (stem bottom)
```

Vertices (counterclockwise):
```
(-60,26), (-33,26), (-32,0), (-28.25,0), (-27.25,26), (27.25,26),
(28.25,0), (32,0), (33,26), (60,26), (60,28), (-60,28)
```

---

## Project 1 Parameters

### Material Properties

| Property | Value | | Property | Value |
|----------|-------|-|----------|-------|
| f'c  | 6.0 ksi | | fpu | 270 ksi |
| f'ci | 4.8 ksi | | fpy | 243 ksi |
| Ec   | 4700 ksi | | Eps | 28500 ksi |

### Allowable Stresses (ACI 318)

| Condition | Compression | Tension |
|-----------|-------------|---------|
| Transfer  | +0.60 f'ci = 2.88 ksi | −3√f'ci = −0.208 ksi |
| Service   | +0.45 f'c = 2.70 ksi  | −12√f'c = −0.929 ksi (Class C) |

### Loads

| Load | Value | Moment at transfer | Moment at service |
|------|-------|--------------------|-------------------|
| Self-weight | `Ac × (150/1728)` kip/in | Mi = w_sw·L²/8 | included |
| SDL (2" topping) | `2 × 120 × (150/1728)` kip/in | not applied at transfer | included |
| Span | 768 in (64 ft) | | |
| η (losses) | 0.85 (15% total) | η=1.0 at transfer | η=0.85 at service |

### Tendon Layout (4 strands, Aps = 0.153 in² each, Pe = 25 kip each)

| Tendon | Profile | y at support | y at midspan |
|--------|---------|-------------|-------------|
| 1, 3 (straight) | constant | 6 in | 6 in |
| 2, 4 (harped) | trapezoidal drape at x = L/2 ± 144 in | 20.36 in | 6 in |

**Total F = 100 kip**, **e at midspan = yc − 6 ≈ 13.44 in**

---

## Common Pitfalls

1. **Denominator sign flip** in Magnel lines: if `Mi/St > fci_allow`, Line I denominator goes negative — inequality flips and Line I becomes an upper bound. Always check dI, dII, dIII, dIV signs.
2. **Coordinate origin mismatch**: eccentricity uses `y` from **bottom**. Some older input file comments reference y from top — ignore those comments, trust the vertex data.
3. **Losses at transfer**: η = 1.0 at transfer (full F). The `defineDesignStages.m` handles this automatically.
4. **Moment at transfer** = self-weight only. SDL is applied after transfer for precast members.
5. **Working directory**: all scripts assume `cd` to `08_matlab_program_analysis/`. Relative paths to `Project1/` will fail from root.

---

## Working Configuration

- **MATLAB working directory:** `08_matlab_program_analysis/`
- **Magnel diagram (standalone):** `feasibilityDesignChart.m`
- **Figure output folder:** `08_matlab_program_analysis/Project1/`
- **Concept derivations (Maple):** `03_concept_notes/FeasibilityChart.mw`, `FeasibilityChart2.mw`

### NotebookLM — Primary Reference

- **Notebook:** CEE 530 Prestressed Concrete
- **Library ID:** `cee-530-prestressed-concrete`
- **Sources loaded:** ACI 318-19, Naaman (2nd ed), Wight RC 8th ed
- **Query command:**
  ```bash
  PYTHONUTF8=1 python scripts/run.py ask_question.py --question "..." --notebook-id cee-530-prestressed-concrete
  ```
  Run from: `C:\Users\chidc\.claude\skills\notebooklm\`

### PDF Sources (direct read fallback)

- `01_sources/Naaman_prestressedconcrete.pdf` — Naaman textbook
- `01_sources/Pretressed Concrete (CEE 530) spring 2024.pdf` — Course lecture slides

**Workflow:** Query NotebookLM first. If answer is incomplete, read the PDF directly.

---

## Decisions Made

| Decision | Reason | Date |
|----------|--------|------|
| Mirror 8-subfolder structure from CEE 790 | Same cognitive pipeline every project | 2026-03-16 |
| All MATLAB scripts in `08_matlab_program_analysis/` | Keeps root clean | 2026-03-16 |
| Rename `Projec1/` → `Project1/` | Fix typo | 2026-03-16 |
| HW files in `04_comparison_tables/` | Matches CEE 790 reference | 2026-03-16 |
| Maple worksheets in `03_concept_notes/` | Concept-level, not analysis outputs | 2026-03-16 |
| Compression = positive sign convention | Consistent with ACI 318 prestress notation | — |
| y = 0 at bottom of stems | Consistent shoelace formula and eccentricity | — |
| Vci_min coeff = 5.0 (CEE 530), not ACI's 1.7 | Prof. Fafitis course convention; extracted as `Vci_min_coeff` variable | 2026-03-24 |
| Signed Vu in shear design plot | Match main analysis SFD convention (+left, −right) | 2026-03-24 |
| FEN v2 report is the base for shear submission | User's hand-polished version; do NOT rebuild from scratch | 2026-03-24 |
| End block: deep beam method (0.8h) as default | Prof. Fafitis CEE 530 course method; Gergely-Sozen (3h/4) available as toggle | 2026-04-08 |
| Flatten folder structure | All .m files in `08_matlab_program_analysis/`, output stays in `project_Project2/output/` | 2026-04-08 |
| Prof. Fafitis teaches CEE 530 | Not Mobasher — corrected attribution in all end block code and docs | 2026-04-08 |

---

## Task Checklist

### Project 1 — Double-T Beam
- [x] Section geometry and properties (shoelace formula)
- [x] Tendon layout (4 strands: 2 straight + 2 harped)
- [x] Stress analysis at transfer and service along beam length
- [x] Magnel diagram (4 lines, feasibility zone)
- [x] Feasibility zone envelope (e_min / e_max vs. x/L)
- [ ] Verify all critical sections pass ACI 318 stress limits
- [ ] Written report in `06_draft_text/`

### Project 2 — Shear, Ultimate, End Block
- [x] Shear design (Vci/Vcw, shearDesign.m)
- [x] End block design (endBlockDesign.m, runEndBlockDesign.m)
- [ ] Ultimate design report (ultimateDesign.m)

### Homework
- [x] HW 1 — input file + section figure (`04_comparison_tables/HW_1/`)
- [x] HW 2 — figures (`04_comparison_tables/HW_2/`, `HW_2_new/`, `HW_2_test/`)
- [ ] HW 3 — End block design (in progress, using runEndBlockDesign.m output)

---

## How to Set Up a New Assignment or Project Folder

```
ProjectName/
├── 01_sources/         ← PDFs, codes
├── 02_literature_review/
├── 03_concept_notes/   ← Maple worksheets, derivations
├── 04_comparison_tables/
├── 05_figures/
├── 06_draft_text/      ← stays empty until concept notes done
├── 07_tasks/inbox.md   ← write THIS FIRST
└── 08_matlab_program_analysis/
```

**Order:** `07_tasks` → `01_sources` → `03_concept_notes` → `04_comparison_tables` → `05_figures` → `06_draft_text`
