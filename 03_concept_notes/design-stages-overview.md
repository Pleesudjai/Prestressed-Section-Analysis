# Design Stages Overview — Prestressed Concrete Beam

**Sources:** ACI 318-19 (via NotebookLM), Naaman 2nd ed, CEE 530 lecture slides
**Sign convention:** Compression = positive (+), Tension = negative (−)

---

## The 5 Design Stages

A prestressed concrete beam must be checked at **5 distinct stages**. Each stage has different forces, moments, and code limits.

| # | Stage | Prestress Force | Concrete Strength | Moment |
|---|-------|----------------|-------------------|--------|
| 1 | **Transfer** | Fi (maximum — no losses yet) | f'ci (minimum) | M_sw only |
| 2 | **Handling / Transport** | Fi (before losses) | f'ci | M from lift points |
| 3 | **Service** | F = η·Fi (after all losses) | f'c (28-day) | M_sw + M_SDL + M_LL |
| 4 | **Ultimate** | fps (at steel yielding) | f'c | Mu = factored loads |
| 5 | **Fatigue** | fpe (effective) | f'c | Cyclic ΔM |

Detailed notes for each stage: see `01_transfer.md`, `02_handling-transport.md`, `03_service.md`, `04_ultimate.md`, `05_fatigue.md`

---

## Key Insight: Why Multiple Stages?

- At **transfer**: prestress is maximum but concrete is weakest → governs compression at top, tension at bottom near midspan
- At **service**: concrete is strongest but moment is maximum → governs bottom tension (Class C) or top compression
- At **ultimate**: prestress effect disappears as steel yields → beam behaves like an RC beam
- The **Magnel diagram** finds the (e, F) combination that satisfies transfer + service simultaneously

---

## Prestress Loss Factor η

All service-stage calculations use the **effective prestress force** after losses:
```
F = η · Fi       where η = 1 − total losses fraction
```
For Project 1: η = 1 − 0.15 = **0.85** (15% total losses assumed)

Long-term losses come from: elastic shortening, creep, shrinkage, steel relaxation.
