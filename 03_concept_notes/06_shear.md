# Stage 6: Shear Design (Strength Limit State)

**Source:** ACI 318-19 §22.5 + Naaman Ch. 6 — cross-verified via NotebookLM
**Sign convention:** All shear/force values are positive magnitudes unless noted

---

## Why Shear in Prestressed Members is Different

In a prestressed beam:
1. **Inclined tendons** carry a vertical component **Vp** that directly resists shear
2. **Prestress compression (fpc)** suppresses diagonal tension → raises shear capacity
3. **Two distinct failure modes** → must check both:
   - **Flexure-shear cracking (Vci)** — inclined crack propagating from a flexural crack
   - **Web-shear cracking (Vcw)** — inclined crack initiating in the web before any flexural cracking

```
Vc = min(Vci, Vcw)
```

---

## Design Requirement

```
φVn ≥ Vu

Vn = Vc + Vs
φ = 0.75  (shear)
```

Where:
- `Vc` = concrete shear contribution (detailed or simplified method)
- `Vs` = transverse reinforcement contribution
- `Vp` = vertical component of prestress (already included inside Vcw)

---

## Naaman Notation vs ACI Notation

Naaman presents equations in **stress units** (divides by bw·dp); ACI uses **force units**. Same physics — different presentation.

| ACI Symbol | Naaman Symbol | Meaning |
|------------|---------------|---------|
| `Vc`, `Vci`, `Vcw` (kip) | `vc`, `vci`, `vcw` (psi) | Concrete shear capacity |
| `fpc` | `σg` | Effective prestress at centroid = F/Ac |
| `Mcre` | `ΔMcr` | Flexural cracking moment (above dead load) |
| `Sb = I/yb` | `Zb` | Bottom section modulus |
| `Vd`, `Vi`, `Mmax` | `VG`, `ΔVu`, `ΔMu` | Dead load shear / incremental factored loads |

**Naaman's load separation** (pedagogical advantage):
- `VG` = shear from self-weight only (unfactored)
- `ΔVu` = factored shear from superimposed dead + live load only
- Since SDL and LL are proportional loads, `ΔVu/ΔMu` stays constant as loads scale → simplifies Vci calculation

ACI lumps both into `Vd` (unfactored dead) and `Vi = Vu − Vd`.

---

## Method 1 — Detailed (§22.5.8.2 & §22.5.8.3)

### Vci — Flexure-Shear Cracking Capacity

Governs in **flexure-dominated regions** (high M, moderate V — typically midspan zone):

**ACI force form:**
```
Vci = 0.6λ√f'c · bw·dp  +  Vd  +  Vi · Mcre / Mmax

Vci ≥ 1.7λ√f'c · bw·dp    (minimum)
```

**Naaman stress form** (divide through by bw·dp):
```
vci = 0.6λ√f'c  +  VG/(bw·dp)  +  (ΔVu/ΔMu) · ΔMcr/(bw·dp)

vci ≥ 1.7λ√f'c    (minimum)
```

Variables (all in **psi and in** units):

| Symbol | ACI | Naaman | Meaning |
|--------|-----|--------|---------|
| `λ` | λ | λ | 1.0 normal-weight, 0.85 sand-LW, 0.75 all-LW |
| `bw` | bw | bw | web width (in) — use sum of all web widths for double-T |
| `dp` | dp | dp | compression fiber to PS steel centroid ≥ **0.80h** |
| `Vd` / `VG` | shear from unfactored DL | shear from self-weight only | (kip) |
| `Vi` / `ΔVu` | `Vu − Vd` | factored shear from SDL+LL | (kip) |
| `Mmax` / `ΔMu` | `Mu − Md` | factored moment from SDL+LL | (kip-in) |
| `Mcre` / `ΔMcr` | cracking moment | cracking moment above MG | (kip-in) |

---

#### Cracking Moment Mcre (ACI) = ΔMcr (Naaman)

**ACI form:**
```
Mcre = Sb · (6λ√f'c + fpe − fd)
```

> `Sb = I/yb` = bottom section modulus (ACI uses `I/yt` where `yt` means distance to extreme **tension** fiber; for a simply supported beam, extreme tension is at **bottom**, so `yt = yb` and `I/yt = Sb`).

| Symbol | Meaning |
|--------|---------|
| `Sb = I/yb` | section modulus at the **tension** (bottom) fiber (in³) |
| `6λ√f'c` | modulus of rupture for shear (psi) — note: **6**, not 7.5 used in flexure |
| `fpe` | effective prestress at **bottom** fiber (psi, compression = positive): `fpe = Pe/Ac + Pe·e/Sb` |
| `fd` | stress at bottom fiber from **unfactored dead load only** (psi, tension = positive): `fd = Md/Sb` |

**Naaman explicit form:**
```
ΔMcr = Sb · (6λ√f'c + fpe_bot) − MG
```

where `fpe_bot = F/Ac + F·e0/Sb` (effective prestress at bottom fiber, psi)
and `MG` = moment from self-weight only.

> Both forms are equivalent. Naaman's separates MG explicitly so the ratio `ΔVu/ΔMu` can be used directly.

---

### Vcw — Web-Shear Cracking Capacity

Governs in **shear-dominated regions** (near supports, low M):

**ACI force form:**
```
Vcw = (3.5λ√f'c + 0.3·fpc) · bw·dp  +  Vp
```

**Naaman stress form:**
```
vcw = 3.5λ√f'c + 0.3·σg + Vp/(bw·dp)
```

| Symbol | ACI | Naaman | Meaning |
|--------|-----|--------|---------|
| `fpc` / `σg` | fpc | σg | effective prestress at centroid (psi) = `Pe/Ac` after all losses |
| `Vp` | Vp | Vp | vertical component of effective prestress (kip) |

> If centroid is in the **flange** of a T or I-section: use fpc at the **junction of flange and web** instead of at the centroid.

---

#### Vertical Component of Prestress Vp

**General:** `Vp = Pe · sin α ≈ Pe · tan α = Pe · |dy/dx|` (small angle)

**For harped (trapezoidal) profile** — constant slope within each segment:
```
Vp = Pe · |Δe / Δx| = Pe · |y_high − y_low| / x_drape_length
```

**For parabolic profile** `y(x) = ax² + bx + c` — slope varies:
```
dy/dx = 2a·x + b
Vp(x) = Pe · |2a·x + b|
```
Parabolic equivalent uniform load (upward): `wp = 8·Pe·sag / L²`

> Naaman: sinα ≈ tanα = Δe/Δx for harped; use derivative dy/dx for parabolic.

> **Vp is maximum at the supports** (steepest slope) and **zero at midspan** (flat portion).

---

## Method 2 — Simplified (§22.5.8.1)

**Validity condition** (ACI §22.5.8.1):
```
Aps · fse  ≥  0.40 · (Aps · fpu + As · fy)
```
Naaman presents this simply as `fpe ≥ 0.40·fpu` (for members with prestress only, no mild steel).
If condition is not met → **must use detailed Vci/Vcw method**.

**Formula:**
```
Vc = [0.6λ√f'c  +  700 · (Vu · dp / Mu)] · bw · dp
```

**Cap on ratio:** `Vu · dp / Mu ≤ 1.0`
(this caps the bracket at `0.6λ√f'c + 700`, same as the upper bound expression)

**Bounds on Vc:**
```
2λ√f'c · bw · dp   ≤   Vc   ≤   5λ√f'c · bw · dp
```

> **Important:** The `1.7λ√f'c · bw · dp` bound applies to **Vci minimum** (detailed method), NOT the simplified method. Simplified minimum = **2λ√f'c · bw · dp**.

> Use `dp ≥ 0.80·h` in this formula.
> Simplified method is conservative near supports — prefer detailed method there.

---

## Transverse Reinforcement (Vs)

**Force equation** (ACI):
```
Vs = Av · fyt · dp / s
```

**Stress equation** (Naaman — equivalent):
```
vs = Av · fyt / (bw · s)    →    Av = (vu/φ − vc) · bw · s / fyt
```

**Required stirrup spacing:**
```
s = Av · fyt · dp / (Vu/φ − Vc)
```

**Maximum spacing** (§9.7.6.2.2):
```
s ≤ 0.75h  (≤ 24 in)           if  Vs ≤ 4√f'c · bw · dp
s ≤ 0.375h (≤ 12 in)           if  Vs > 4√f'c · bw · dp
```

Naaman equivalents: `0.75h` and `3h/8 = 0.375h` — same as ACI.

**Maximum shear steel can carry** (§22.5.1.2):
```
Vs,max = 8√f'c · bw · dp      (if exceeded → enlarge section)
```

---

## Minimum Transverse Reinforcement (§9.6.3)

Required when `Vu > 0.5·φ·Vc`:

**General minimum:**
```
Av,min = 0.75λ√f'c · bw · s / fyt    ≥    50 · bw · s / fyt      (psi, in)
```

**Prestressed member alternative** (when `fpe ≥ 0.40·fpu`, per Naaman):
```
Av,min = Aps · fpu · s / (80 · fyt · dp) · √(bw / dp)
```

ACI allows using the **lesser** of the two alternatives for prestressed members.

---

## Critical Section Location

| Code / Reference | Critical Section |
|-----------------|-----------------|
| **ACI §9.4.3.2** | `h/2` from face of support (for prestressed) |
| **Naaman** | `h/2` from face of support |

> ACI uses `d` from face for ordinary RC but **`h/2`** for prestressed members.
> Shear reinforcement designed at the `h/2` section is then applied conservatively to all sections closer to the support.

---

## Where to Check Along the Span

| Region | Governing Mode |
|--------|---------------|
| Near support (x ≤ ~L/5) | **Vcw** governs (high V, low M, steep tendon slope) |
| Midspan zone | **Vci** governs (high M, lower V, flat tendon) |
| Transition | Compute both; use `Vc = min(Vci, Vcw)` |

---

## Shear Flow — Composite Sections (precast + cast-in-place topping)

```
vu = Vu / (φ · Acv)             (horizontal shear stress at interface)
```

```
Vnh = μ · (Avf · fy + Pc)       (shear-friction, ACI §16.4)
```

| `μ` | Interface condition |
|-----|-------------------|
| 1.0 | Intentionally roughened (amplitude ≥ 1/4 in) |
| 0.6 | Not intentionally roughened |

---

## Summary: Shear Design Steps

```
1.  Compute Vu along span  (Vu = 1.2·VD + 1.6·VL)
2.  Identify critical section: h/2 from face of support
3.  Check need for stirrups: if Vu ≤ 0.5·φ·Vc → minimum Av or none
4.  Compute Vc:
       Simplified: Vc = [0.6λ√f'c + 700·min(Vu·dp/Mu, 1.0)]·bw·dp
                   bounded by [2λ√f'c·bw·dp, 5λ√f'c·bw·dp]
       Detailed:   Vc = min(Vci, Vcw) — compute at each section
5.  Check capacity: Vu ≤ φ·(Vc + 8√f'c·bw·dp) → else enlarge section
6.  Required Vs = Vu/φ − Vc
7.  Stirrup spacing: s = Av·fyt·dp / Vs
8.  Apply max spacing (0.75h or 0.375h) and min Av
9.  Report design spacing along span
```

---

## Project 1 Values (Double-T, f'c = 6 ksi, Pe = 85 kips)

| Property | Value |
|----------|-------|
| `bw` | 2 × 4.75 ≈ **9.5 in** (two stems) |
| `dp` | h − y_tendon = 28 − 6 = **22 in** |
| `dp check` | 0.80 × 28 = 22.4 in → use **22.4 in** (governs) |
| `fpc = σg` | 85 / Ac ksi → × 1000 for psi |
| `Vp` at support (harped) | Pe × (20.36 − 6)/240 = 85 × 0.0598 ≈ **5.1 kips** |
| `Vp` at midspan | ≈ 0 kips (flat portion of harped tendon) |
| `φ` (shear) | **0.75** |
| `λ` | **1.0** (normal-weight) |

---

## MATLAB Sketch

```matlab
%% Shear design — simplified method (ACI §22.5.8.1)
lambda  = 1.0;
fc_psi  = fc * 1000;      % ksi → psi
bw      = 2 * 4.75;       % in, two stems
h       = 28;             % in
dp      = max(22, 0.80*h);% in, dp >= 0.80h

% Ratio cap
ratio = min(Vu .* dp ./ max(abs(Mu), 1e-3), 1.0);   % cap at 1.0

% Simplified Vc (psi formula × bw·dp / 1000 = kip)
Vc_raw = (0.6*lambda*sqrt(fc_psi) + 700*ratio) .* bw .* dp / 1000;
Vc_lo  = 2.0*lambda*sqrt(fc_psi)*bw*dp / 1000;    % min = 2λ√f'c·bw·dp
Vc_hi  = 5.0*lambda*sqrt(fc_psi)*bw*dp / 1000;    % max = 5λ√f'c·bw·dp
Vc     = min(max(Vc_raw, Vc_lo), Vc_hi);

% Required Vs
phi    = 0.75;
Vs_req = max(0, Vu./phi - Vc);

% Stirrup spacing (Av = 0.40 in^2 for #4 @ 2 legs per stem = 4 legs total)
Av  = 0.40 * 2;   % in^2, 2 stems × 2 legs
fyt = 60;         % ksi
s_req = Av * fyt * dp ./ max(Vs_req, 1e-6);   % in

% Max spacing limits
Vs_thresh = 4*lambda*sqrt(fc_psi)*bw*dp/1000;   % threshold for reduced spacing
s_max = zeros(size(Vs_req));
s_max(Vs_req <= Vs_thresh) = min(0.75*h, 24);
s_max(Vs_req >  Vs_thresh) = min(0.375*h, 12);

s_design = min(s_req, s_max);   % governing spacing

fprintf('At midspan:\n');
fprintf('  Vc = %.1f kips\n', Vc(mid_idx));
fprintf('  Vs req = %.1f kips\n', Vs_req(mid_idx));
fprintf('  s design = %.1f in\n', s_design(mid_idx));
```

---

## Common Pitfalls

1. **f'c units**: Vci/Vcw formulas use **psi** for f'c and all stresses. Convert `fpe`, `fpc`/`σg` from ksi → psi before substituting.
2. **yt vs yb in Mcre**: ACI writes `I/yt` where `yt` = distance to extreme **tension** fiber. For a simply supported beam, tension is at the **bottom** → `yt = yb` → `I/yt = Sb`. Do not confuse with `St = I/yt_top`.
3. **Simplified method lower bound = 2, not 1.7**: `1.7λ√f'c·bw·dp` is the minimum for **Vci** (detailed). The simplified method minimum is `2λ√f'c·bw·dp`.
4. **Vu·dp/Mu ratio**: Must be **capped at 1.0** in the simplified formula. Without the cap, the formula can overestimate Vc near supports.
5. **dp ≥ 0.80h**: Both methods require dp used in the formula to be no less than 0.80h.
6. **Mmax ≠ Mu**: In detailed Vci, `Mmax = Mu − Md` (factored minus dead-load moment). Naaman calls this `ΔMu`.
7. **Vp sign**: Vp acts upward at supports (resists shear), zero at midspan flat zone. Do not double-count Vp — it is already embedded in Vcw.
8. **Two stems**: Double-T has `bw = 2 × stem width` and typically `Av = 4 legs` (2 per stem).
9. **Simplified validity**: Check `fpe ≥ 0.40·fpu` (or full ACI condition `Aps·fse ≥ 0.40·(Aps·fpu + As·fy)`) before using simplified method.
