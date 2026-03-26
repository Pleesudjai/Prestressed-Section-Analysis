# Tasks — CEE 530 Prestressed Concrete

## Project 1 — Open Items

### Stress Analysis
- [ ] Run `main_PrestressedBeamAnalysis.m` with full service loads (SW + SDL + LL)
- [ ] Verify all critical sections pass ACI 318 stress limits at transfer and service
- [ ] Print console summary: pass/fail for each section

### Magnel Diagram
- [ ] Run `feasibilityDesignChart.m` — confirm 4 Magnel lines plot correctly
- [ ] Check feasibility zone is non-empty (design is feasible)
- [ ] Verify design point (e, 1/F) falls inside feasible zone

### Feasibility Zone
- [ ] Plot e_min and e_max envelope vs. x/L along full span
- [ ] Overlay tendon profile (straight + harped) on same plot
- [ ] Confirm tendon CG stays inside envelope at all x

### Report
- [ ] Write Project 1 report in `06_draft_text/`
- [ ] Include: section properties table, stress diagrams, Magnel diagram, feasibility zone plot

---

## Homework — Open Items

- [ ] HW 2 finalized? Check `04_comparison_tables/HW_2_new/HW2.pdf`

---

---

## MATLAB Program — Open Items (from 08_matlab_program_analysis/)

### Shear Design (not yet implemented)
- [ ] Implement shear design in `analyzePrestressedBeam.m` — compute Vci, Vcw, Vc along span
- [ ] Add simplified Vc method: `[0.6λ√f'c + 700·Vu·dp/Mu]·bw·dp`, bounded [2λ, 5λ]·√f'c·bw·dp
- [ ] Compute Vp(x) from tendon slope: `Vp = Pe · |dy/dx|` (harped: constant per segment; parabolic: derivative)
- [ ] Design stirrup spacing s(x) and enforce max spacing (0.75h/24 in, 0.375h/12 in)
- [ ] Add Av,min check (general + prestressed alternative)
- [ ] Add shear plots to `plotPrestressedBeamResults.m` — Vu, φVc, φVn diagrams

### Lsection Project — Residual Warning
- [ ] Investigate `M(L) = −144 kip-in ≠ 0` boundary condition warning in `project_Lsection`
  - Likely cause: partial-span live load (x = 48 to 432) reaction equilibrium edge case
  - Fix: verify how `analyzePrestressedBeam` handles partial distributed loads at endpoints

### Concept Notes — Still to Write
- [ ] `07_deflection.md` — short-term and long-term deflection, ACI limits, PCI multiplier method
- [ ] `08_losses.md` — prestress loss summary: elastic shortening, creep, shrinkage, relaxation (ACI + Naaman)

---

## Open Questions

1. Does the current tendon layout (4 strands) satisfy all ACI 318 service stress limits?
2. Is the feasible zone wide enough to allow tendon adjustment if needed?
3. What is the governing stress condition at midspan — bottom tension at service?

---

## Professor / Assignment Requirements

- Add items here as they come up from lectures or assignment sheets
