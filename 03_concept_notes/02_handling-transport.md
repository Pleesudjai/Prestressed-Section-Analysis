# Stage 2: Handling and Transportation

**Source:** Naaman Ch. 3, ACI 318-19 §26.10 — confirmed via NotebookLM
**Sign convention:** Compression = positive (+), Tension = negative (−)

---

## When Does This Occur?

Between fabrication and final installation:
- **Form removal / stripping** — beam is lifted off the form bed
- **Storage** — beam rests on dunnage (temporary supports)
- **Transportation** — beam is loaded on a truck or rail car
- **Erection** — beam is lifted and set in final position

Prestress force = **Fi** (before long-term losses)
Concrete strength = **f'ci** (same as transfer stage)

---

## What Makes This Stage Different from Transfer?

The **support locations change**. In the final structure the beam is simply supported at the ends. During handling, the beam may be:
- Lifted from two pickup points near the ends
- Stored on dunnage blocks at intermediate points
- Cantilevered over the truck during transport

The **bending moment diagram is completely different** from the in-service configuration. Lifting points near the ends create negative moment at midspan. Supports near midspan create positive moment at the ends — possibly causing tension at the top.

---

## Loads Considered

1. **Self-weight** (w_sw) — always present
2. **Dynamic impact factor** — additional load for acceleration/deceleration
   - Typical: 1.0 to 1.33× self-weight (no specific ACI value; use engineering judgment or project spec)
   - AASHTO provides IM = 33% for vehicular live loads on bridges
3. **Temporary loads** during erection (rigging weight, equipment)

---

## ACI 318-19 Stress Limits During Handling

Same temporary limits as transfer (f'ci governs, losses have not occurred yet):

| Type | Limit | Project 1 Value (f'ci = 4.8 ksi) |
|------|-------|----------------------------------|
| Compression | ≤ +0.60 f'ci | +2.88 ksi |
| Tension (no bonded rebar) | ≥ −0.0948√f'c (ksi) | ~−0.208 ksi |
| Tension (piles, AASHTO) | ≥ −0.158√f'c (ksi) | — |

---

## Effect of Lifting Point Location

**Optimal lift points** minimize the maximum moment magnitude:
- For a simply supported beam of length L, symmetric lift points at `x = 0.207L` from each end minimize the peak moment
- At those points: M_max = w·(0.207L)²/2 (negative cantilever) = w·(0.586L)²/8 (positive midspan) → both equal

> **Source note (NotebookLM verified):** The 0.207L rule and the 1.0–1.33× dynamic impact factor are from the **PCI Design Handbook**, not from Naaman's book or ACI 318-19. Naaman's sources only cite IM = 33% for vehicular live loads on bridges (AASHTO). For building precast, use PCI recommendations.

**For Project 1:** This stage is not the primary analysis focus, but should be considered if the beam is a precast Double-T transported to the site.

---

## Practical Notes

- Precast Double-T beams are particularly sensitive to lateral stability during lifting
- Tipping and rollover must be checked for long, slender members
- Lifting loops are typically at 0.10L to 0.20L from each end

---

## In MATLAB

Not explicitly modeled in `main_PrestressedBeamAnalysis.m` for Project 1.
If needed, modify the support locations and recompute M(x) with the new boundary conditions.
