# Shear Design Calculation Template
# ACI 318-19 §22.5 — Detailed Method (Vci / Vcw)

**Source:** ACI 318-19 §22.5.8.3 + Naaman Ch. 6
**Sign convention:** All shear/force values are positive magnitudes; compression = +, tension = −
**Units:** kip, in. throughout (convert f'c to ksi when using √f'c in ksi form)

---

## GIVEN INFORMATION

### Section Properties
| Symbol | Description | Value |
|--------|-------------|-------|
| `h`  | Total depth | ___ in. |
| `bw` | Min. web width (all stems combined) | ___ in. |
| `Ac` | Gross area | ___ in² |
| `Ix` | Moment of inertia (centroidal) | ___ in⁴ |
| `yc` | Centroid from bottom | ___ in. |
| `yb` | Dist. centroid → bottom fiber | ___ in. |
| `yt` | Dist. centroid → top fiber | ___ in. |
| `Sb = Ix/yb` | Bottom section modulus | ___ in³ |

### Material Properties
| Symbol | Value |
|--------|-------|
| `f'c` | ___ ksi  →  √f'c = ___ ksi |
| `f'ci` (at transfer) | ___ ksi |
| `λ` | 1.0 (normal-weight) |
| `φ` (shear) | 0.75 |
| `fy` (stirrups) | ___ ksi |
| `fpu` | ___ ksi |

### Prestress
| Symbol | Value |
|--------|-------|
| `Pi` total | ___ kips |
| `η = 1 − losses` | ___ |
| `Pe = η · Pi` | ___ kips |
| `fpc = Pe / Ac` | ___ ksi  (average prestress at centroid) |

### Tendon Profile
| Tendon | Type | y at support | y at midspan | Harp point |
|--------|------|-------------|-------------|-----------|
| ___ | Straight/Harped | ___ in. | ___ in. | x = ___ in. |

### Applied Loads
| Load | kip/in. | kip/ft |
|------|---------|--------|
| `wsw` (self-weight) | ___ | ___ |
| `wsdl` (superimposed DL) | ___ | ___ |
| `wll` (live load) | ___ | ___ |
| `wDL = wsw + wsdl` | ___ | ___ |
| `wu = 1.2wDL + 1.6wll` | ___ | ___ |
| `wext = 1.2wsdl + 1.6wll` | ___ | ___ |

---

## STEP 1 — Factored Shear Envelope

```
Vu(x)  =  wu × (L/2 − x)
        =  [wu] × (L/2 − x)

Vu,max  =  wu × L/2  =  [wu] × [L/2]  =  ___ kips   (at support face)
```

---

## STEP 2 — Critical Section Location

```
dp,min  =  0.8 h  =  0.8 × [h]  =  ___ in.
```

> Compute actual dp from tendon centroid at each section.
> Use dp = max(dp,calc, 0.8h) everywhere.

```
xcr  =  dp,min  =  ___ in.  =  ___ ft from support face
Vu at xcr  =  wu × (L/2 − xcr)  =  ___ kips
```

---

## STEP 3 — Vertical Prestress Component, Vp

For **harped tendons** (linear profile, constant slope within each segment):

```
θ  =  arctan[(y_support − y_mid) / x_harp]
   =  arctan[(___ − ___) / ___]
   =  ___ °   →   sin θ = ___

Pe per harped strand  =  Pi,strand × η  =  ___ × ___  =  ___ kips

Vp  (sloped zone, x ≤ x_harp)
   =  (no. harped strands) × Pe,strand × sin θ
   =  ___ × ___ × ___
   =  ___ kips   (upward)

Vp  (flat zone, x_harp < x ≤ L − x_harp)  =  0
```

For **parabolic tendons** (if used):
```
wp  =  8 · Pe · sag / L²         (equivalent uniform upward load)
Vp(x)  =  wp × (L/2 − x)
```

---

## STEP 4 — Modulus of Rupture and Key Constants

```
fr  =  6 λ √f'c    [f'c in ksi, result in ksi]
    =  6 × [λ] × √[f'c]
    =  ___ ksi

Vci,min  =  1.7 λ √f'c · bw · dp
          =  1.7 × [λ] × [√f'c] × [bw] × [dp]
          =  ___ kips   (constant lower bound for Vci)

Vcw base  =  (3.5 λ √f'c  +  0.3 fpc) · bw · dp
           =  (3.5 × [λ] × [√f'c]  +  0.3 × [fpc]) × [bw] × [dp]
           =  ___ kips   (add Vp for sloped zone)
```

---

## STEP 5 — Tendon Position at Each Section x

For a harped tendon sloping from y_s (at support) to y_m (at midspan) with harp point at x_h:

```
yharped(x)  =  y_s  +  (y_m − y_s) × (x / x_h)        for x ≤ x_h
            =  y_m                                       for x_h < x ≤ L/2
```

Centroid of all PS strands at section x:

```
yps(x)  =  [ Σ y_i(x) ] / n_strands
```

Eccentricity and effective depth:

```
e(x)    =  yc − yps(x)
dp(x)   =  max(h − yps(x),  0.8h)
```

---

## STEP 6 — Mcr: Cracking Moment at Each Section

```
fce(x)  =  Pe/Ac  +  Pe · e(x) · yb / Ix
         =  [Pe]/[Ac]  +  [Pe] × e(x) × [yb] / [Ix]

Md(x)   =  wDL × x × (L − x) / 2
fd(x)   =  Md(x) × yb / Ix

(fr + fce − fd) check:
   If  fr + fce − fd > 0:  Mcr = (Ix/yb) × (fr + fce − fd)  =  Sb × (fr + fce − fd)
   If  fr + fce − fd ≤ 0:  Mcr = 0  (section pre-cracked under DL alone → Vci = Vci,min)
```

---

## STEP 7 — Vci: Flexural-Shear Cracking Strength

**ACI 318-19 Eq. 22.5.8.3.1a:**

```
Vci  =  0.6λ√f'c · bw · dp   +   Vd   +   (Vi / Mmax) · Mcr

where:
  term1  =  0.6λ√f'c · bw · dp   =  ___ kips   (constant)
  Vd     =  wDL × (L/2 − x)      =  unfactored DL shear
  Vi     =  wext × (L/2 − x)     =  factored shear from loads applied after decompression
  Mmax   =  wext × x × (L − x)/2 =  factored moment from same loading as Vi

Vci  ≥  Vci,min  =  1.7λ√f'c · bw · dp  =  ___ kips
```

> **Mcr = 0** → Vci = term1 + Vd (no flexural-shear cracking reserve). Vci,min still applies.

---

## STEP 8 — Vcw: Web-Shear Cracking Strength

**ACI 318-19 Eq. 22.5.8.3.2:**

```
Vcw  =  (3.5λ√f'c  +  0.3 fpc) · bw · dp   +   Vp
      =  Vcw,base  +  Vp

Sloped zone (x ≤ x_harp):   Vcw  =  Vcw,base  +  Vp,harped
Flat zone   (x > x_harp):   Vcw  =  Vcw,base  +  0
```

---

## STEP 9 — Concrete Shear Capacity

```
Vc  =  min(Vci,  Vcw)
φVc  =  0.75 × Vc
```

Check: `Vu vs. φVc`

| Result | Action |
|--------|--------|
| Vu ≤ 0.5 φVc | No stirrups required |
| 0.5 φVc < Vu ≤ φVc | Minimum stirrups only |
| Vu > φVc | Design stirrups for Vs,req |

---

## STEP 10 — Stirrup Design

```
Vs,req  =  Vu / φ  −  Vc        (if Vu > φVc, else Vs,req = 0)

Av/s (demand)  =  Vs,req / (fy · dp)
```

**Minimum Av/s — largest of three (ACI §22.5.10.5) [f'c and fy in ksi]:**

```
(Av/s)₁  =  0.75√f'c · bw / fy     [f'c in ksi units inside √; fy in ksi]
          =  0.75 × [√f'c(ksi)] × [bw] / [fy]  ... but formula requires psi:
          =  0.75 × √(f'c×1000) × bw / (fy×1000)  in²/in

(Av/s)₂  =  50 · bw / (fy × 1000)   [f'c in psi convention]

(Av/s)₃  =  Aps · fpu / (80 · fy · dp) × √(dp / bw)

(Av/s)min  =  max of above three
```

> Quick values for f'c = 6 ksi, bw = ___ in., dp = ___ in., fy = 60 ksi, Aps = ___ in²:
> - (Av/s)₁ = 0.75 × √6000 × bw / 60000 = ___ in²/in
> - (Av/s)₂ = 50 × bw / 60000 = ___ in²/in
> - (Av/s)₃ = ___ × 270 / (80 × 60 × dp) × √(dp/bw) = ___ in²/in

**Stirrup spacing:**

```
s_demand  =  Av / (Av/s demand)
s_min     =  Av / (Av/s)min

USE  s  =  min(s_demand, s_min, s_max)
```

**Spacing limits (ACI §22.7.6.2):**

```
Vs,threshold  =  4λ√f'c · bw · dp   [kip]

If Vs,req ≤ Vs,threshold:  s_max = min(3h/4,  24 in.)  =  min(___, 24) = ___ in.
If Vs,req > Vs,threshold:  s_max = min(3h/8,  12 in.)  =  min(___, 12) = ___ in.
```

**Section size check (ACI §22.5.10.1):**

```
Vs,max  =  8λ√f'c · bw · dp
         =  8 × [λ] × [√f'c] × [bw] × [dp]
         =  ___ kips

If Vs,req > Vs,max  →  Enlarge section.
```

---

## RESULTS TABLE (fill one row per analysis section)

| x (ft) | Vu (k) | yps (in) | e (in) | fce (ksi) | fd (ksi) | Mcr (k-in) | term1 | Vd | Vi | Mmax | Vci (k) | Vcw (k) | Vc (k) | φVc (k) | Vs,req (k) | Av/s (in²/in) | s_use (in) |
|--------|--------|----------|--------|-----------|----------|------------|-------|----|----|------|---------|---------|--------|---------|-----------|--------------|-----------|
| xcr = ___ | | | | | | | 7.81 | | | | | | | | | | |
| | | | | | | | | | | | | | | | | | |
| L/2 | 0 | | | | | | | 0 | 0 | — | Vci,min | | | | 0 | (Av/s)min | smax |

---

## RECOMMENDED STIRRUP LAYOUT

| Zone | x from support | Spacing | Note |
|------|---------------|---------|------|
| Support — xcr | 0 to ___ in. | s = ___ in. | Vcw governs |
| High shear | ___ to ___ in. | s = ___ in. | Vci governs; demand |
| Transition | ___ to ___ in. | s = ___ in. | Vci,min; demand or smax |
| Mid-beam | ___ in. to L/2 | s = smax = ___ in. | Min Av/s |

Beam is symmetric — mirror layout for right half.

---

## COMMON PITFALLS

1. **f'c units in formulas**: the coefficient 0.6, 1.7, 3.5 assume f'c in **psi** in ACI original.
   When computing in ksi (√f'c in ksi), the numerical coefficient changes:
   - 0.6√f'c_psi = 0.6 × √(f'c_ksi × 1000) = 0.6 × 31.623 × √f'c_ksi = **18.97 √f'c_ksi**
   - Alternatively: use f'c_ksi directly and keep all forces in kips — √6000 psi = 77.46; √6 ksi = 2.449 → multiply by 1000/1000 = same.
   - **Safest**: keep f'c in ksi, use √f'c (ksi), multiply coefficient by √1000 = 31.623 internally, or just use the formula with √f'c_ksi and confirm units cancel to kips.

2. **dp minimum**: dp = max(h − yps, 0.8h). Must apply **at each section** — dp varies if yps varies.

3. **Mmax ≠ Mu (total)**: In Vci, Mmax is the moment from the **incremental factored loads** (loads applied after decompression = SDL + LL factored = wext). It is NOT the total factored Mu.

4. **Vi definition**: Vi = Vu,total − Vd (ACI) = factored shear from wext only. Since the beam is simply supported with uniform load: Vi = wext × (L/2 − x).

5. **Mcr = 0 zone**: Occurs near midspan when large DL moment exceeds prestress + rupture. Vci = Vci,min there. Vcw still governs at those sections if Vcw > Vci,min.

6. **Vp already inside Vcw**: Do NOT add Vp separately to Vc after taking min(Vci, Vcw). It is already embedded in Vcw.

7. **Two-stem sections (Double-T)**: bw = sum of both stem widths at the critical level. Av = 2 legs per stem × 2 stems = 4 legs (if using one stirrup per stem). Av = 2 legs per stem × 2 stems = 2 × 2 × Abar.

8. **Mcr sign check**: fr + fce − fd < 0 means bottom fiber is in net tension under DL alone even with prestress → section has zero cracking reserve → Mcr = 0, Vci = Vci,min. This is physically meaningful; the section is likely inadequate in flexure at service too.

---

## EQUATIONS REFERENCE CARD

```
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  SHEAR DESIGN — ACI 318-19 DETAILED METHOD
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  fr    = 6λ√f'c                            [ksi, f'c in ksi]

  fce   = Pe/Ac + Pe·e·yb / Ix             [ksi, compression +]

  fd    = Md · yb / Ix                      [ksi, tension magnitude]

  Mcr   = (Ix/yb) · (fr + fce − fd)        [kip-in, ≥ 0]
        = Sb · (fr + fce − fd)

  Vci   = 0.6λ√f'c·bw·dp + Vd + Vi·Mcr/Mmax    [kip]
        ≥ 1.7λ√f'c·bw·dp  (minimum)

  Vcw   = (3.5λ√f'c + 0.3·fpc)·bw·dp + Vp      [kip]

  Vc    = min(Vci, Vcw)                     [kip]

  Vs    = Vu/φ − Vc        (when Vu > φVc)  [kip]

  Av/s  = max(Vs/(fy·dp),  0.75√f'c·bw/fy,  50·bw/fy,
              Aps·fpu/(80·fy·dp)·√(dp/bw))  [in²/in]

  s     = Av / (Av/s)                       [in]
        ≤ min(3h/4, 24 in.)  if Vs ≤ 4λ√f'c·bw·dp
        ≤ min(3h/8, 12 in.)  if Vs > 4λ√f'c·bw·dp

  Vs,max = 8λ√f'c·bw·dp                    [kip, else enlarge section]
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

---

## PROJECT 2 QUICK CONSTANTS (Double-T, f'c = 6 ksi, bw = 7.5 in., dp = 22.40 in.)

| Constant | Formula | Value |
|----------|---------|-------|
| `√f'c` | √6.0 | 0.07746 ksi^0.5 → √6000 = 77.46 psi^0.5 |
| `fr` | 6×1.0×0.07746 | 0.4648 ksi |
| `term1 (Vci)` | 0.6×1.0×0.07746×7.5×22.40 | 7.808 kips |
| `Vci,min` | 1.7×1.0×0.07746×7.5×22.40 | 22.12 kips |
| `Vcw,base` | (3.5×0.07746 + 0.3×0.1745)×7.5×22.40 | 54.34 kips |
| `Vcw (sloped, Vp=2.538k)` | 54.34 + 2.538 | 56.88 kips |
| `Vcw (flat,   Vp=0)` | 54.34 + 0 | 54.34 kips |
| `Vs,max` | 8×1.0×0.07746×7.5×22.40 | 104.2 kips |
| `Vs,threshold` | 4×1.0×0.07746×7.5×22.40 | 52.17 kips |
| `smax (basic)` | min(3×28/4, 24) | 21.0 in. |
| `(Av/s)min` | max of three criteria | 0.00726 in²/in. |
| `s from (Av/s)min, #3` | 0.22/0.00726 | 30.3 in. → use 21 in. |
