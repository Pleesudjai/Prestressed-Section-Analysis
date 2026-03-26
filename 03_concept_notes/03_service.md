# Stage 3: Service (Serviceability Limit State)

**Source:** ACI 318-19 §24.5, Naaman Ch. 3–4 — confirmed via NotebookLM
**Sign convention:** Compression = positive (+), Tension = negative (−)

---

## When Does This Occur?

Under normal in-service loading conditions after all prestress losses have occurred:
- Prestress force = **F = η · Fi** (effective, after all losses)
- Concrete strength = **f'c** (28-day strength, full)
- Moment = **MT = M_sw + M_SDL + M_LL** (all loads)

---

## Fiber Stress Equations

```
f_top_s = η·Fi/Ac  −  η·Fi·e/St  +  MT/St      (top fiber)
f_bot_s = η·Fi/Ac  +  η·Fi·e/Sb  −  MT/Sb      (bottom fiber)
```

Where:
- `η` = loss factor (η = 1 − total losses fraction); **Project 1: η = 0.85**
- `MT = M_sw + M_SDL + M_LL` = total service moment
- `St = Ic/yt`, `Sb = Ic/yb`

---

## ACI 318-19 Allowable Stress Limits at Service

### Compression Limits
| Condition | Limit | Project 1 Value (f'c = 6.0 ksi) |
|-----------|-------|----------------------------------|
| Prestress + sustained loads (SW + SDL) | ≤ +0.45 f'c | +2.70 ksi |
| Prestress + total loads (SW + SDL + LL) | ≤ +0.60 f'c | +3.60 ksi |

> The 0.45 f'c limit prevents excessive creep under sustained loads.

### Tension Limits — Beam Classification

> **Two editions in use.** Class U boundary differs between the CEE 530 course (older ACI) and ACI 318-19.
> `inputData.m` has a `materials.code_edition` switch — set `'ACI-318-19'` or `'CEE530'`.

| Class | ACI-318-19 | CEE530 (course edition) | Behavior | Analysis |
|-------|------------|-------------------------|----------|----------|
| **U** (Uncracked) | ft ≤ **7.5√f'c** (psi) | ft ≤ **6√f'c** (psi) | No cracking | Gross section |
| **T** (Transition) | 7.5 < ft ≤ 12√f'c | 6 < ft ≤ 12√f'c | Transition | Cracked section for deflection |
| **C** (Cracked) | ft > **12√f'c** (psi) | ft > **12√f'c** (psi) (same) | Cracked | Cracked section + crack control |

**Project 1 / Project 2:** f'c = 6.0 ksi → √6000 = 77.46 psi^0.5
- Class U limit (ACI-318-19): −7.5 × 77.46 = **−0.581 ksi**
- Class U limit (CEE530): −6 × 77.46 = **−0.465 ksi**
- Class C limit (both): −12 × 77.46 = **−0.929 ksi**

> Two-way prestressed slabs must always be Class U (both editions).

---

## Dead-Only Service Case

When checking **SDL only** (before LL is applied):
```
MT_dead = M_sw + M_SDL
f_top_s_dead = η·Fi/Ac  −  η·Fi·e/St  +  MT_dead/St
f_bot_s_dead = η·Fi/Ac  +  η·Fi·e/Sb  −  MT_dead/Sb
```
This is handled by `inputPrestressedBeam_Project1_DeadOnly.m`.

---

## Magnel Diagram — 4 Conditions at Midspan

> **NotebookLM verified (Naaman):** Naaman's canonical line numbering is given below. The CLAUDE.md/MATLAB convention uses a different ordering — see note at end of this section.

Rearranged as **1/Fi or 1/F vs. e** using kern points: `kb = r²/yb`, `kt = r²/yt`.

### Naaman's Convention (governing conditions for typical beams)

**Line I** — Top fiber @ Transfer, **tension** limit (fti = −3√f'ci psi):

Physical reason: Near supports where Mi = 0, the eccentricity causes top tension.
```
fi_top = Fi/Ac − Fi·e/St  ≥  −fti_allow
→  1/Fi ≤ (1/Ac − e/St) / (−fti_allow)   = (Mmin/Zt − σ̄ti) / (e − kb)   [upper bound]
```

**Line II** — Bottom fiber @ Transfer, **compression** limit (fci = +0.60 f'ci):

Physical reason: Near supports, full prestress concentrates compression at the bottom.
```
fi_bot = Fi/Ac + Fi·e/Sb  ≤  fci_allow
→  1/Fi ≥ (1/Ac + e/Sb) / fci_allow   = (Mmin/Zb + σ̄ci) / (e − kt)   [lower bound]
```

**Line III** — Top fiber @ Service, **compression** limit (fcs = +0.45 f'c):
```
fs_top = F/Ac − F·e/St + MT/St  ≤  fcs_allow
→  1/F ≥ η·(1/Ac − e/St) / (fcs_allow − MT/St)   [lower bound]
```

**Line IV** — Bottom fiber @ Service, **tension** limit (fts = −12√f'c psi, Class C):
```
fs_bot = F/Ac + F·e/Sb − MT/Sb  ≥  fts_allow
→  1/F ≤ η·(1/Ac + e/Sb) / (MT/Sb − fts_allow)   [upper bound]
```

**Feasible zone:** `max(Line II, Line III) ≤ 1/F ≤ min(Line I, Line IV)`

> **Pitfall:** If a denominator changes sign as e varies, the inequality **flips**. Always check signs of all denominators.

---

### Why These 4? — Physical Intuition

| Stage | Fiber | Dominant Effect | Governing Limit |
|-------|-------|-----------------|-----------------|
| Transfer (x ≈ support) | Top | Prestress eccentricity → **tension** | Tension limit |
| Transfer (x ≈ support) | Bottom | Prestress → large **compression** | Compression limit |
| Service (x ≈ midspan) | Top | Large moment → **compression** | Compression limit |
| Service (x ≈ midspan) | Bottom | Large moment → **tension** | Tension limit |

---

### CLAUDE.md / MATLAB Convention

The `feasibilityDesignChart.m` uses a different labeling that checks conditions at midspan:

| Line | Fiber | Stage | Limit type | Bound |
|------|-------|-------|------------|-------|
| I | Top | Transfer | **Compression** | lower |
| II | Bottom | Transfer | **Tension** | upper |
| III | Bottom | Service | Tension | upper |
| IV | Top | Service | Compression | lower |

These are valid constraints (they capture midspan behavior) but Lines I and II are **less restrictive** than Naaman's for typical beams (the governing transfer conditions are at the support, not midspan). For a conservative design, **Naaman's conditions should govern**. The MATLAB code produces the correct Magnel diagram shape but may underestimate how tight the transfer bounds are near the support.

---

## Eccentricity of the Tendon Group — e_cg(x)

### Why We Need a Single e

The Magnel diagram has **two axes**: `1/F` (vertical) and `e` (horizontal).
`F` is the **total** prestress force in the beam, and `e` is the **eccentricity of the resultant of that total force**.

When there is only one tendon, `e` is simply `yc − y_tendon`.
When there are **multiple tendons at different heights and with different forces**, we must replace all of them with a single equivalent force `F` acting at a single equivalent eccentricity `e_cg`. That is the centroid of the tendon group.

---

### Physical Concept — Centre of a Force System

Think of each tendon `i` as carrying a force `Pi` at height `y_i` from the bottom.
The resultant `F = ΣPi` acts at the **force-weighted centroid**:

```
y_cg = Σ( Pi × y_i ) / Σ( Pi )
```

The eccentricity measured from the section centroid is:

```
e_cg = yc − y_cg
     = yc − [ Σ( Pi × y_i ) / Σ( Pi ) ]
     = Σ( Pi × (yc − y_i) ) / Σ( Pi )
     = Σ( Pi × e_i ) / Σ( Pi )
```

So equivalently:

```
            Σ( Pi_i × e_i(x) )
e_cg(x) = ─────────────────────
                 Σ( Pi_i )

where  Pi_i  = Aps_i × fpi_i     (initial force of tendon i, kip)
       e_i(x) = yc − y_i(x)      (eccentricity of tendon i at position x)
```

> **Why force-weighted, not area-weighted?**
> The Magnel diagram treats `F` as the total force and `e` as its line of action.
> A larger strand carries more force and therefore pulls the resultant centroid more.
> Using area alone would be wrong if strands have different `fpi` (e.g. different wire sizes or different initial stress levels).
> The correct weight is always the **force** `Pi = Aps_i × fpi_i`.

---

### How e_cg Varies Along the Beam

`e_cg(x)` is not constant — it follows the tendon profiles.

**Straight tendon:** `y_i(x) = constant` → `e_i(x) = constant`

**Harped (trapezoidal) tendon:**
```
y_h(x) = y_sup  +  (y_drape − y_sup) × x / x_drape      [0 ≤ x ≤ x_drape]
y_h(x) = y_drape                                          [x_drape < x ≤ L/2]
```
So `e_h(x) = yc − y_h(x)` varies linearly from support to drape point, then is constant.

For Project 1 (2 straight + 2 harped, all equal Pi):
```
e_cg(x) = [ Pi × e_straight  +  Pi × e_harped(x) ] / (2 × Pi)
         = [ e_straight + e_harped(x) ] / 2          (Pi cancels when equal)
```
Because Pi cancels only when all strands have the same force. In the general case it does NOT cancel.

---

### Step-by-Step: Project 1 Numerical Example

**Given:**
```
yc        = 20.362 in  (centroid from bottom of stems)
y_straight = 6.0 in    (strands 1, 3 — constant)
y_sup_h    = 20.362 in  (strands 2, 4 at support — same as yc)
y_drape_h  = 6.0 in    (strands 2, 4 at and beyond x = 240 in)
x_drape    = 240 in    (20 ft from support)
Aps_i      = 0.153 in² (all four the same)
fpi_i      = 163.4 ksi (all four the same)
Pi         = 0.153 × 163.4 = 25.0 kip  (all four equal)
```

**At x = 0 ft (support):**
```
y_h(0)   = 20.362 in   →   e_h = 20.362 − 20.362 = 0.000 in
e_straight = 20.362 − 6.0 = 14.362 in

e_cg = [ 2×25×14.362 + 2×25×0.000 ] / (4×25)
     = [ 718.1 + 0 ] / 100
     = 7.181 in
```

**At x = 10 ft (120 in):**
```
y_h(120) = 20.362 + (6.0 − 20.362)×120/240 = 13.181 in
e_h      = 20.362 − 13.181 = 7.181 in

e_cg = [ 2×25×14.362 + 2×25×7.181 ] / 100
     = [ 718.1 + 359.05 ] / 100
     = 10.772 in
```

**At x = 20 ft (240 in) — drape point:**
```
y_h(240) = 6.0 in   →   e_h = 14.362 in  (same as straight)

e_cg = [ 2×25×14.362 + 2×25×14.362 ] / 100
     = 14.362 in     (all tendons at same depth → e_cg = e_individual)
```

**At x = 32 ft (384 in, midspan):**
```
y_h = 6.0 in (constant past drape point)   →   e_h = 14.362 in

e_cg = 14.362 in   (same as drape point)
```

Summary table:

| x (ft) | y_h (in) | e_straight (in) | e_harped (in) | **e_cg (in)** |
|--------|----------|-----------------|---------------|---------------|
| 0  | 20.362 | 14.362 | 0.000  | **7.181**  |
| 10 | 13.181 | 14.362 | 7.181  | **10.772** |
| 20 |  6.000 | 14.362 | 14.362 | **14.362** |
| 32 |  6.000 | 14.362 | 14.362 | **14.362** |

---

### Why e_cg < e_individual at the Support

At the support, the straight tendons have their maximum eccentricity (14.36 in below centroid) but the harped tendons are **at the centroid** (e = 0). The combined group centroid is pulled upward toward `yc` by the harped strands, giving `e_cg = 7.18 in` — much less than the 14.36 in of the straight tendons alone.

This is the physical reason the Magnel feasibility zone at x = 0 is **much less constrained** than at midspan — a smaller e_cg means the prestress moment `F × e_cg` is smaller, so it is easier to satisfy the stress limits.

---

### What Happens with Unequal Strand Forces

If you have two groups of strands with **different Pi** (e.g. 0.5" strands at 25 kip mixed with 0.6" strands at 36 kip):

```
e_cg = (n_A × Pi_A × e_A  +  n_B × Pi_B × e_B)
       ──────────────────────────────────────────
             n_A × Pi_A  +  n_B × Pi_B
```

The group with higher force pulls `e_cg` more toward its own eccentricity.
**Do not use simple area averaging** — it will give the wrong resultant location.

---

### In MATLAB (feasibility4Sections.m)

```matlab
% Force-weighted centroid of tendon group at each analysis section
% Pi_i = Aps_i * fpi_i  (initial prestress force, may differ between tendons)

e_cg_sec = zeros(1, n_sec);
for k = 1:n_sec
    xi     = x_sec(k);
    Pe_sum = 0;    % sum of Pi * e_i
    P_sum  = 0;    % sum of Pi
    for i = 1:n_tendons
        t   = prestress.tendons{i};
        Pi  = t.Aps * t.fpi;                                   % force this tendon
        e_i = interp1(beam.x, t.e, xi, 'linear', 'extrap');   % e at this x
        Pe_sum = Pe_sum + Pi * e_i;
        P_sum  = P_sum  + Pi;
    end
    e_cg_sec(k) = Pe_sum / P_sum;
end
```

`t.e` is the eccentricity profile already computed by `processTendonProfiles` inside `inputData.m`, using `e_i(x) = yc − y_i(x)` for each tendon at every x in `beam.x`.
`interp1` evaluates it at the specific section location `xi`.

---

## Deflection and Camber (Serviceability Check)

| Effect | Direction | Cause |
|--------|-----------|-------|
| Camber (upward) | ↑ | Prestress eccentricity |
| Elastic deflection (downward) | ↓ | Applied loads |
| Long-term camber growth | ↑ | Creep under sustained prestress |
| Long-term deflection | ↓ | Creep under sustained loads |

**ACI limits (confirmed by Naaman):**

| Condition | Deflection limit |
|-----------|-----------------|
| Flat roofs, elements NOT likely to be damaged | L/180 (immediate live load) |
| Floors, elements NOT likely to be damaged | L/360 (immediate live load) |
| Roofs/floors supporting elements likely to be damaged | L/480 (immediate LL + long-term additional) |
| Roofs/floors NOT supporting damageable elements | L/240 (immediate LL + long-term additional) |

---

## In MATLAB

```matlab
F_eff = eta * Fi;          % effective prestress after losses
MT = M_sw + M_SDL + M_LL;  % total service moment

f_top_s = F_eff/Ac - F_eff*e/St + MT/St;
f_bot_s = F_eff/Ac + F_eff*e/Sb - MT/Sb;

% Check limits
f_cs_allow = 0.45 * fc;          % compression
f_ts_allow = -12*sqrt(fc_psi)/1000;  % tension Class C (ksi)
```
