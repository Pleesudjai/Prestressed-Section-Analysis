# Stage 5: Fatigue Limit State

**Source:** ACI 318-19, ACI Committee 215, Naaman — confirmed via NotebookLM
**Applies to:** Bridges, parking structures, industrial floors — any member under repetitive loading

---

## When Does This Matter?

Prestressed concrete members under **repeated cyclic loading** (vehicular traffic, machinery, wave loads) must be checked for fatigue. Prestressing steels do NOT have a true endurance limit — they can fail by fatigue at stress ranges below their static strength.

**Design fatigue life:** ≥ 2 million cycles (standard per ACI Committee 215)

---

## What Is Checked?

The **stress range** (Δf) — the difference between maximum and minimum stress in one load cycle — must stay below the allowable fatigue limit.

```
Δf = f_max − f_min    (during one load cycle)
```

This is NOT the absolute stress, but the fluctuation amplitude.

---

## Allowable Stress Ranges (ACI Committee 215)

| Steel Type | Limit | Formula |
|-----------|-------|---------|
| Prestressing strands | Δfps ≤ 0.10 fpu | For f'pu = 270 ksi → Δfps ≤ 27 ksi |
| Prestressing wires | Δfps ≤ 0.12 fpu | For fpu = 250 ksi → Δfps ≤ 30 ksi |
| Mild steel (non-prestressed) | Δfs ≤ 20–24 ksi | Depending on bar deformation pattern |

**For cracked (Class C) or partially prestressed members:**
- Limit steel stress range to **Δfps ≤ 16 ksi (112 MPa)**
- Also control crack width range to prevent fretting between wires and concrete

---

## Concrete Fatigue

> ⚠️ **Correction (NotebookLM/Naaman):** `f_max/f'c > 0.60` is the **static service** compression limit (ACI §24.5), NOT a fatigue limit.

For **repetitive fatigue loading**, the correct limits are:

| Authority | Limit |
|-----------|-------|
| ACI (bridges) | Maximum compression under service load ≤ **0.50 f'c** |
| ACI Committee 215 | Compressive stress range: `fcr ≤ 0.40·f'c + 0.47·f_min` |

Where `f_min` = minimum compressive stress in the cycle (positive = compression).

The static service limit of 0.45 f'c provides an indirect margin, but is not itself a fatigue criterion.

---

## Relationship to the Other Stages

| Stage | Force | Moment |
|-------|-------|--------|
| Service (max) | η·Fi | MT_max (all loads) |
| **Fatigue check** | η·Fi | ΔM = M_LL only (live load cycles) |
| Service (min) | η·Fi | MT_min (dead loads only, no LL) |

The stress **range** due to live load cycling is what matters — dead load and prestress are constant (no range contribution).

---

## Practical Implication for Project 1

Project 1 is a **floor/roof Double-T** in a building — not a highway bridge. Fatigue is generally not the governing design stage for building structures unless there is cyclic machinery loading. For this project, **fatigue does not govern** but should be noted as a required check in real design.

---

## In MATLAB (not implemented for Project 1)

```matlab
% Stress range due to live load cycling
f_bot_max = F_eff/Ac + F_eff*e/Sb - MT_max/Sb;   % max moment (all loads)
f_bot_min = F_eff/Ac + F_eff*e/Sb - M_dead/Sb;   % min moment (dead only)

delta_f = f_bot_max - f_bot_min;   % stress range at bottom (where tendons are)

% Allowable
delta_f_allow = 0.10 * fpu;   % for strands (compression positive convention)
```
