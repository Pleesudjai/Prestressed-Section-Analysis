# Stage 1: Transfer (Initial Prestress)

**Source:** ACI 318-19 §24.5, Naaman Ch. 3 — confirmed via NotebookLM
**Sign convention:** Compression = positive (+), Tension = negative (−)

---

## When Does This Occur?

Immediately after the prestressing force is released (transferred) to the concrete. Most critical condition because:
- Prestress force = **Fi (maximum)** — no losses have occurred yet
- Concrete strength = **f'ci (minimum)** — young concrete, typically 70–80% of f'c
- External moment = **self-weight only** (Mi = M_sw) — no SDL or LL yet

---

## Fiber Stress Equations

At any section x along the span:

```
f_top_i = Fi/Ac  −  Fi·e/St  +  Mi/St      (top fiber)
f_bot_i = Fi/Ac  +  Fi·e/Sb  −  Mi/Sb      (bottom fiber)
```

Where:
- `Fi` = initial prestress force (before losses)
- `e` = eccentricity = yc − y_tendon (positive when tendon is below centroid)
- `St = Ic/yt`, `Sb = Ic/yb` (section moduli, top and bottom)
- `Mi` = M_sw = self-weight moment = w_sw · x · (L−x) / 2

**Sign of each term:**
- `Fi/Ac` → compression → positive
- `−Fi·e/St` → if e > 0 (tendon below CG), this term is negative → tension at top
- `+Mi/St` → hogging moment creates compression at top → positive
- `+Fi·e/Sb` → compression at bottom → positive
- `−Mi/Sb` → tension at bottom from moment → negative

---

## Allowable Stress Limits at Transfer — Code Edition Comparison

> **Two editions in use.** The CEE 530 course uses an older ACI standard; ACI 318-19 differs in one limit.
> `inputData.m` has a `materials.code_edition` switch — set `'ACI-318-19'` or `'CEE530'`.

| Location | Type | ACI-318-19 | CEE530 (course edition) |
|----------|------|------------|-------------------------|
| General | Compression | ≤ **+0.60 f'ci** | ≤ **+0.60 f'ci** (same) |
| End regions | Compression | ≤ **+0.70 f'ci** | ≤ **+0.60 f'ci** (no boost) |
| General | Tension | ≥ **−3√f'ci** (psi) | ≥ **−3√f'ci** (psi) (same) |
| End regions | Tension | ≥ **−6√f'ci** (psi) | ≥ **−6√f'ci** (psi) (same) |

**Project 1 / Project 2:** f'ci = 4.8 ksi
- Compression (general): 0.60 × 4.8 = **+2.88 ksi** ← both editions
- Compression (ends, ACI-318-19): 0.70 × 4.8 = **+3.36 ksi**
- Compression (ends, CEE530): 0.60 × 4.8 = **+2.88 ksi** (same as general)
- Tension (general): −3√4800 = **−0.208 ksi** ← both editions
- Tension (ends): −6√4800 = **−0.416 ksi** ← both editions

> If tension exceeds −6√f'ci, bonded mild steel reinforcement must be added to carry all tension.

---

## What Governs?

- At **midspan**: Bottom fiber tension from moment is relieved by prestress → may be OK
- At **end (support)**: No moment → full prestress force acts → top fiber tension can be critical
- The tension limit at transfer is usually the **most restrictive** condition

---

## In MATLAB (`analyzePrestressedBeam.m`)

```matlab
f_top_i = Fi/Ac - Fi*e/St + Mi/St;
f_bot_i = Fi/Ac + Fi*e/Sb - Mi/Sb;

% Check limits
f_ci_allow = 0.60 * fci;    % compression limit
f_ti_allow = -3*sqrt(fci_psi)/1000;  % tension limit (ksi)
```
