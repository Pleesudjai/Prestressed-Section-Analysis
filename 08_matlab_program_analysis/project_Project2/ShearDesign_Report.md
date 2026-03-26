# Project 1 — Precast Prestressed Beam
## Assignment 4: Shear Design

**CEE 530 — Prestressed Concrete**
Due 3/23/2026

*Group 3*
- Hailey Crosson
- Reyaneh Kazemian
- Leohuns Robert Kambu
- Prince Aminu
- Chidchanok Pleesudjai (Fen)
- Haribhushan Sirigineedi

---

## 1. Given Information

### 1.1 Section Geometry

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Total section depth | h | 28.0 in. |
| Flange width | b_f | 120.0 in. |
| Flange thickness | h_f | 2.0 in. |
| Min. web width (2 stems) | b_w | 7.50 in. (2 × 3.75 in.) |
| Span length | L | 768 in. (64.0 ft) |
| Support condition | — | Simply supported |

### 1.2 Section Properties (Shoelace Formula)

| Property | Symbol | Value |
|----------|--------|-------|
| Gross area | A_c | 487.0 in² |
| Moment of inertia (centroidal) | I_x | 34,638.8 in⁴ |
| Centroid from bottom | y_c | 20.362 in. |
| Dist. centroid to bottom fiber | y_b | 20.362 in. |
| Dist. centroid to top fiber | y_t | 7.638 in. |
| Bottom section modulus S_b = I_x/y_b | S_b | 1,701.5 in³ |
| Top section modulus S_t = I_x/y_t | S_t | 4,532.6 in³ |

### 1.3 Material Properties

| Material | Property | Value |
|----------|----------|-------|
| Concrete | f'c | 6.0 ksi (6,000 psi) |
| | E_c | 4,700 ksi |
| | λ | 1.0 (normal-weight concrete) |
| Prestressing steel | f_pu | 270 ksi |
| | E_ps | 28,500 ksi |
| Mild steel (stirrups) | f_y | 60 ksi |

### 1.4 Allowable Stresses

| Stage | f'c or f'ci | Compression (ksi) | Tension (ksi) |
|-------|-------------|-------------------|---------------|
| Transfer | f'ci = 4.8 ksi | +0.60 f'ci = +2.880 | −3√f'ci = −0.208 |
| Service **(use CEE530 factor)** | f'c = 6.0 ksi | +0.45 f'c = +2.700 | −12√f'c = −0.929 |

### 1.5 Prestressing Tendons

| Tendon | x-location | Profile | y (support) | y (midspan) | P_e/strand |
|--------|-----------|---------|-------------|-------------|------------|
| 1 (left) | −30.125 in. | Straight | 6.00 in. | 6.00 in. | 21.25 kips |
| 2 (left) | −30.125 in. | Harped | 20.36 in. | 6.00 in. | 21.25 kips |
| 3 (right) | +30.125 in. | Straight | 6.00 in. | 6.00 in. | 21.25 kips |
| 4 (right) | +30.125 in. | Harped | 20.36 in. | 6.00 in. | 21.25 kips |

4 strands × A_ps = 0.153 in²/strand → A_ps,total = 0.612 in².
Harp point at x = 240 in. (20 ft) from each support.
P_i = 100.0 kips; η = 0.85; P_e = 85.0 kips.

### 1.6 Applied Loads

| Load | Symbol | kip/in. | kip/ft |
|------|--------|---------|--------|
| Self-weight SW | w_sw | 0.04227 | 0.5073 |
| Superimposed dead load SDL | w_sdl | 0.02083 | 0.2500 |
| Live load LL | w_ll | 0.03500 | 0.4200 |
| Dead load total DL = SW+SDL | w_DL | 0.06310 | 0.7573 |
| Factored w_u = 1.2DL + 1.6LL | w_u | 0.13173 | 1.5807 |
| Fact. var. w_ext = 1.2SDL + 1.6LL | w_ext | 0.08100 | 0.9720 |

---

## 2. Shear Design Procedure

For prestressed members, ACI 318-19 §22.5.8.3 defines two failure modes:

- **V_ci** (Flexural-shear cracking) — inclined crack starts from a flexural crack.
- **V_cw** (Web-shear cracking) — inclined crack forms in web before any flexure crack.

Concrete shear capacity: V_c = min(V_ci, V_cw). Stirrups required when V_u > φV_c, carrying V_s = V_u/φ − V_c.
And when V_u ≤ φ(V_c + V_s).

### 2.1 Critical Section Location

> **Comment:** Considering as at 1.5 ft. You need to draw the actual shear force diagram and track down the moment and shear that happen at the design section.

The critical section for shear is located at d_p from the face of support, where d_p shall not be taken less than 0.8h.

d_p,min = 0.8 h = 0.8 × 28.0 in. = 22.40 in. ← governs throughout (actual d_p < 22.40 in. near supports)

All calculations use d_p = 22.40 in. Critical section location:

x_cr = d_p,min = 22.40 in. = 1.87 ft from support face (I use 1.5 ft)

### 2.2 Factored Shear Envelope

For a simply supported beam with uniform factored load, the shear envelope is:

V_u(x) = w_u × (L/2 − x) = 0.13173 × (384 − x) [kips, x in inches]

V_u,max (at support face, x = 0) = w_u × L/2 = 0.13173 × 384 = 50.58 kips

V_u (at x_cr = 22.40 in.) = 0.13173 × (384 − 22.40) = 0.13173 × 361.60 = 47.63 kips

### 2.3 Vertical Component of Effective Prestress, V_p

Tendons 2 and 4 are harped; tendons 1 and 3 are straight (V_p = 0). Harped tendons slope from y = 20.36 in. at supports to y = 6.00 in. at harp point x = 240 in.

θ_harp = arctan[(y_end − y_mid) / x_harp] = arctan[(20.36 − 6.00) / 240] = arctan(14.36 / 240) = arctan(0.05983) = 3.424° → sin θ = 0.05973

P_e per harped strand = P_i,strand × η = 25.0 × 0.85 = 21.25 kips

V_p (for x ≤ 240 in., sloped zone) = 2 strands × P_e,strand × sin θ = 2 × 21.25 × 0.05973 = 2.538 kips (upward — reduces net V_u on web)

For 240 in. < x ≤ 528 in. (flat zone): V_p = 0.

### 2.4 Modulus of Rupture and Cracking Moment

All cracking moment equations reference ACI 318-19 §22.5.8.3.3 and §22.5.8.3.1b.

f_r = 6 λ √f'c [ksi, f'c in ksi] = 6 × 1.0 × √6.0 = 0.4648 ksi

Effective prestress compressive stress at extreme tension fiber (bottom):

f_ce = P_e/A_c + P_e·e·y_b / I_x = 85.0/487.0 + (85.0 × e × 20.362) / 34,638.8 ← evaluated at each section (e varies with tendon profile)

Unfactored DL moment and flexural tensile stress at bottom fiber:

M_d(x) = w_DL × x × (L − x) / 2 = 0.06310 × x × (768 − x) / 2 [kip-in., x in inches]

f_d(x) = M_d × y_b / I_x = M_d × 20.362 / 34,638.8

Cracking moment (ACI 318-19 Eq. 22.5.8.3.1b):

M_cr = (I_x / y_b) × (f_r + f_ce − f_d) = (34,638.8 / 20.362) × (f_r + f_ce − f_d) = 1,701.5 × (f_r + f_ce − f_d) [kip-in.]

Note: If f_r + f_ce − f_d < 0, section is already cracked under DL alone → M_cr = 0.

### 2.5 Flexural-Shear Cracking Strength, V_ci — ACI 318-19 Eq. 22.5.8.3.1a

V_ci = 0.6 λ √f'c · b_w · d_p + V_d + (V_i / M_max) · M_cr

where:

| Symbol | Definition | Formula |
|--------|-----------|---------|
| term1 | 0.6λ√f'c·b_w·d_p (concrete contribution) | 0.6×1.0×0.07746×7.5×22.40 = 7.808 kips (constant) |
| V_d | Unfactored DL shear | w_DL×(L/2−x) = 0.06310×(384−x) |
| V_i | Factored shear from loads applied after decompression | w_ext×(L/2−x) = 0.08100×(384−x) |
| M_max | Factored moment at x from same loading as V_i | w_ext×x×(L−x)/2 = 0.08100×x×(768−x)/2 |

Minimum value (always applies as lower bound):

> **⚠️ CHANGE:** V_ci,min = **5** λ √f'c · b_w · d_p (changed from 1.7 coefficient)
>
> = 5 × 1.0 × 0.07746 × 7.5 × 22.40

### 2.6 Web-Shear Cracking Strength, V_cw — ACI 318-19 Eq. 22.5.8.3.2

V_cw = (3.5 λ √f'c + 0.3 f_pc) · b_w · d_p + V_p

where f_pc = average prestress at centroid:

f_pc = P_e / A_c = 85.0 / 487.0 = 0.1745 ksi (constant throughout beam)

Concrete + prestress term (constant, since b_w, d_p, f_pc are constant):

(3.5λ√f'c + 0.3f_pc)·b_w·d_p = (3.5×1.0×0.07746 + 0.3×0.1745) × 7.5 × 22.40
= (0.2711 + 0.05235) × 168.0 = 0.32345 × 168.0 = 54.34 kips

V_cw (sloped zone, x ≤ 240 in., V_p = 2.538 kips) = 54.34 + 2.538 = 56.88 kips

V_cw (flat zone, 240 < x ≤ 528 in., V_p = 0) = 54.34 + 0 = 54.34 kips

### 2.7 Minimum Transverse Reinforcement — ACI 318-19 §22.5.10.5

Three criteria; the largest governs:

(A_v/s)₁ = 0.75 √f'c · b_w / f_y = 0.75 × √6,000 × 7.5 / 60,000 = 0.00726 in²/in. ← GOVERNS

(A_v/s)₂ = 50 · b_w / f_y = 50 × 7.5 / 60,000 = 0.00625 in²/in.

(A_v/s)₃ = A_ps·f_pu / (80·f_y·d_p) × √(d_p/b_w) = 0.612×270 / (80×60×22.40) × √(22.40/7.50) = 0.00266 in²/in.

Stirrup selected: #3 closed U-stirrups → A_v = 2 × 0.11 = 0.22 in².

s (from min A_v/s) = A_v / (A_v/s)_min = 0.22 / 0.00726 = 30.3 in. → limited by s_max

### 2.8 Spacing Limits — ACI 318-19 §22.7.6.2

s_max (basic) = min(3h/4, 24 in.) = min(3×28/4, 24) = min(21.0, 24.0) = 21.0 in.

V_s,threshold (for tight zone) = 4 λ √f'c · b_w · d_p = 4 × 1.0 × 0.07746 × 7.5 × 22.40 = 52.17 kips

When V_s > 52.17 kips: s_max,tight = min(3h/8, 12 in.) = 10.5 in. (Does NOT apply here — max V_s,req = 22.77 kips < 52.17 kips throughout.)

---

## 3. Detailed Calculations at Selected Sections

Results at five sections are shown. The beam is symmetric; the right half mirrors the left.

### 3.1 Critical Section: x = 22.40 in. (1.87 ft from support)

**Step A — Tendon Position and Effective Depth**

y_harped(x) = 20.36 + (6.00 − 20.36) × (22.40 / 240) = 20.36 − 1.340 = 19.02 in.

y_ps = (6.00 + 19.02 + 6.00 + 19.02) / 4 = 12.51 in.

e = y_c − y_ps = 20.362 − 12.51 = 7.852 in. (tendon below centroid → positive)

d_p,calc = h − y_ps = 28.0 − 12.51 = 15.49 in. < 0.8h = 22.40 in. → use d_p = 22.40 in.

**Step B — f_ce: Effective Prestress at Bottom Fiber**

f_ce = P_e/A_c + P_e·e·y_b / I_x = 85.0/487.0 + (85.0 × 7.852 × 20.362) / 34,638.8 = 0.1745 + 0.3924 = 0.5669 ksi (compression +)

**Step C — f_d: DL Flexural Stress at Bottom Fiber**

M_d = w_DL × x × (L − x) / 2 = 0.06310 × 22.40 × (768 − 22.40) / 2 = 0.06310 × 22.40 × 372.80 = 527.0 kip-in.

f_d = M_d × y_b / I_x = 527.0 × 20.362 / 34,638.8 = 0.3098 ksi (tension −)

**Step D — M_cr: Cracking Moment**

f_r + f_ce − f_d = 0.4648 + 0.5669 − 0.3098 = 0.7219 ksi > 0 → M_cr > 0

M_cr = (I_x / y_b) × (f_r + f_ce − f_d) = (34,638.8 / 20.362) × 0.7219 = 1,701.5 × 0.7219 = 1,228 kip-in.

**Step E — V_ci: Flexural-Shear Cracking Strength**

V_d = w_DL × (L/2 − x) = 0.06310 × (384 − 22.40) = 0.06310 × 361.60 = 22.82 kips

V_i = w_ext × (L/2 − x) = 0.08100 × (384 − 22.40) = 0.08100 × 361.60 = 29.29 kips

M_max = w_ext × x × (L − x) / 2 = 0.08100 × 22.40 × (768 − 22.40) / 2 = 0.08100 × 22.40 × 372.80 = 676.4 kip-in.

term1 = 0.6λ√f'c·b_w·d_p = 0.6 × 1.0 × 0.07746 × 7.5 × 22.40 = 7.808 kips

term2 = V_d = 22.82 kips

term3 = (V_i/M_max)·M_cr = (29.29 / 676.4) × 1,228 = 0.04331 × 1,228 = 53.18 kips

V_ci = term1 + term2 + term3 = 7.808 + 22.82 + 53.18 = 83.81 kips

V_ci,min = 1.7λ√f'c·b_w·d_p = 1.7 × 1.0 × 0.07746 × 7.5 × 22.40 = 22.12 kips < 83.81 kips → V_ci = 83.81 kips

**Step F — V_cw: Web-Shear Cracking Strength**

x = 22.40 in. ≤ 240 in. → sloped tendon zone → V_p = 2.538 kips.

V_cw = (3.5λ√f'c + 0.3f_pc)·b_w·d_p + V_p = 54.34 + 2.538 = 56.88 kips

**Step G — V_c, φV_c, and Stirrup Design**

V_c = min(V_ci, V_cw) = min(83.81, 56.88) = 56.88 kips (V_cw governs)

φV_c = 0.75 × 56.88 = 42.66 kips

V_u = 47.63 kips > φV_c = 42.66 kips → Stirrups required.

V_s,req = V_u/φ − V_c = 47.63/0.75 − 56.88 = 63.51 − 56.88 = 6.63 kips

A_v/s (demand) = V_s,req / (f_y · d_p) = 6.63 / (60 × 22.40) = 0.00493 in²/in. < A_v/s,min = 0.00726 in²/in.

s = A_v / (A_v/s,min) = 0.22 / 0.00726 = 30.3 in. → limited by s_max → **USE s = 21 in.**

### 3.2 Section at x = 63.6 in. (5.30 ft)

**Step A — Tendon Position**

y_harped = 20.36 + (6.00 − 20.36) × (63.6 / 240) = 20.36 − 3.805 = 16.56 in.

y_ps = (6.00 + 16.56 + 6.00 + 16.56) / 4 = 11.28 in.

e = y_c − y_ps = 20.362 − 11.28 = 9.082 in.

**Step B — f_ce**

f_ce = 85.0/487.0 + (85.0 × 9.082 × 20.362) / 34,638.8 = 0.1745 + 0.4537 = 0.6282 ksi

**Step C — f_d**

M_d = 0.06310 × 63.6 × (768 − 63.6) / 2 = 0.06310 × 63.6 × 352.2 = 1,413.5 kip-in.

f_d = 1,413.5 × 20.362 / 34,638.8 = 0.8309 ksi

**Step D — M_cr**

f_r + f_ce − f_d = 0.4648 + 0.6282 − 0.8309 = 0.2621 ksi > 0 → M_cr > 0

M_cr = 1,701.5 × 0.2621 = 446.2 kip-in.

**Step E — V_ci**

V_d = 0.06310 × (384 − 63.6) = 0.06310 × 320.4 = 20.22 kips

V_i = 0.08100 × 320.4 = 25.95 kips

M_max = 0.08100 × 63.6 × (768 − 63.6) / 2 = 0.08100 × 63.6 × 352.2 = 1,814.3 kip-in.

term1 = 7.808 kips (constant)

term2 = V_d = 20.22 kips

term3 = (V_i/M_max)·M_cr = (25.95 / 1,814.3) × 446.2 = 0.01431 × 446.2 = 6.384 kips

V_ci = 7.808 + 20.22 + 6.384 = 34.41 kips > V_ci,min = 22.12 kips → V_ci = 34.41 kips

**Step F — V_cw**

x = 63.6 in. ≤ 240 in. → V_p = 2.538 kips.

V_cw = 54.34 + 2.538 = 56.88 kips

**Step G — V_c and Stirrup Design**

V_c = min(34.41, 56.88) = 34.41 kips ← V_ci governs

φV_c = 0.75 × 34.41 = 25.81 kips

V_u = 42.21 kips > φV_c = 25.81 kips → Stirrups required.

V_s,req = 42.21/0.75 − 34.41 = 56.28 − 34.41 = 21.87 kips

A_v/s (demand) = 21.87 / (60 × 22.40) = 0.01627 in²/in. > A_v/s,min = 0.00726 → demand governs

s = 0.22 / 0.01627 = 13.5 in. → **USE s = 13 in.** (≤ s_max = 21 in. ✓)

### 3.3 Section at x = 128.4 in. (10.70 ft)

**Step A — Tendon Position**

y_harped = 20.36 + (6.00 − 20.36) × (128.4 / 240) = 20.36 − 7.683 = 12.68 in.

y_ps = (6.00 + 12.68 + 6.00 + 12.68) / 4 = 9.34 in.

e = 20.362 − 9.34 = 11.02 in.

**Step B — f_ce**

f_ce = 85.0/487.0 + (85.0 × 11.02 × 20.362) / 34,638.8 = 0.1745 + 0.5508 = 0.7253 ksi

**Step C — f_d**

M_d = 0.06310 × 128.4 × (768 − 128.4) / 2 = 0.06310 × 128.4 × 319.8 = 2,591.2 kip-in.

f_d = 2,591.2 × 20.362 / 34,638.8 = 1.5232 ksi

**Step D — M_cr (Section cracked under DL alone)**

f_r + f_ce − f_d = 0.4648 + 0.7253 − 1.5232 = −0.333 ksi < 0 → M_cr = 0

Interpretation: The effective prestress (f_ce = 0.7253 ksi) plus modulus of rupture (f_r = 0.4648 ksi) together cannot overcome the DL tensile stress (f_d = 1.5232 ksi) at the bottom fiber. The (V_i/M_max)·M_cr term vanishes. V_ci reduces to its minimum.

V_ci = 1.7λ√f'c·b_w·d_p = 22.12 kips (minimum governs)

**Step E — V_cw**

x = 128.4 in. ≤ 240 in. → V_p = 2.538 kips.

V_cw = 54.34 + 2.538 = 56.88 kips

**Step F — V_c and Stirrup Design**

V_c = min(22.12, 56.88) = 22.12 kips ← V_ci,min governs

φV_c = 0.75 × 22.12 = 16.59 kips

V_u = 33.67 kips > φV_c = 16.59 kips → Stirrups required.

V_s,req = 33.67/0.75 − 22.12 = 44.89 − 22.12 = 22.77 kips (MATLAB exact at this x: 20.95 kips)

A_v/s (demand) = 20.95 / (60 × 22.40) = 0.01559 in²/in. > A_v/s,min → demand governs

s = 0.22 / 0.01559 = 14.1 in. → **USE s = 14 in.** (≤ s_max = 21 in. ✓)

### 3.4 Section at x = 192.0 in. (16.0 ft)

**Step A — Tendon Position**

y_harped = 20.36 + (6.00 − 20.36) × (192.0 / 240) = 20.36 − 11.486 = 8.874 in.

y_ps = (6.00 + 8.874 + 6.00 + 8.874) / 4 = 7.437 in.

e = 20.362 − 7.437 = 12.925 in.

**Step B — f_ce**

f_ce = 85.0/487.0 + (85.0 × 12.925 × 20.362) / 34,638.8 = 0.1745 + 0.6459 = 0.8204 ksi

**Step C — f_d**

M_d = 0.06310 × 192.0 × (768 − 192.0) / 2 = 0.06310 × 192.0 × 288.0 = 3,489.4 kip-in.

f_d = 3,489.4 × 20.362 / 34,638.8 = 2.051 ksi

**Step D — M_cr**

f_r + f_ce − f_d = 0.4648 + 0.8204 − 2.051 = −0.766 ksi < 0 → M_cr = 0 → V_ci = V_ci,min = 22.12 kips

**Step E — V_cw**

x = 192.0 in. ≤ 240 in. → V_p = 2.538 kips. → V_cw = 56.88 kips.

**Step F — V_c and Stirrup Design**

V_c = min(22.12, 56.88) = 22.12 kips (V_ci,min governs)

φV_c = 0.75 × 22.12 = 16.59 kips

V_u = 25.29 kips > φV_c = 16.59 kips → Stirrups required.

V_s,req = 25.29/0.75 − 22.12 = 33.72 − 22.12 = 11.60 kips

A_v/s (demand) = 11.60 / (60 × 22.40) = 0.00863 in²/in. > A_v/s,min = 0.00726 → demand governs

s = 0.22 / 0.00863 = 25.5 in. → limited by s_max → **USE s = 21 in.**

### 3.5 Midspan: x = 384.0 in. (32.0 ft)

At midspan all four tendons are at y = 6.0 in. V_u = 0 (by symmetry). V_p = 0 (flat zone).

**Step A — Tendon Position**

y_ps (midspan) = 6.0 in. (all strands)

e = y_c − y_ps = 20.362 − 6.0 = 14.362 in. ← maximum eccentricity

**Step B — f_ce**

f_ce = 85.0/487.0 + (85.0 × 14.362 × 20.362) / 34,638.8 = 0.1745 + 0.7178 = 0.8923 ksi

**Step C — f_d**

M_d = w_DL × L² / 8 = 0.06310 × 768² / 8 = 0.06310 × 73,728 / 8 = 4,652.6 kip-in. (maximum DL moment)

f_d = 4,652.6 × 20.362 / 34,638.8 = 2.735 ksi

**Step D — M_cr**

f_r + f_ce − f_d = 0.4648 + 0.8923 − 2.735 = −1.378 ksi < 0 → M_cr = 0

The maximum DL tensile stress at the bottom (2.735 ksi) far exceeds the combined prestress and rupture capacity (0.892 + 0.465 = 1.357 ksi). V_ci defaults to its minimum value.

V_ci = 22.12 kips (minimum governs)

**Step E — V_cw (flat zone)**

V_cw (V_p = 0) = (3.5λ√f'c + 0.3f_pc)·b_w·d_p + 0 = 0.32345 × 168.0 + 0 = 54.34 kips

**Step F — V_c and Stirrup Check**

V_c = min(22.12, 54.34) = 22.12 kips (V_ci,min governs)

φV_c = 0.75 × 22.12 = 16.59 kips

V_u = 0 kips < φV_c = 16.59 kips → No stirrups required by strength; use minimum.

s (minimum) = A_v / (A_v/s,min) = 0.22 / 0.00726 = 30.3 in. → **USE s = 21 in.** (limited by s_max)

---

## 4. Shear Design Summary

### 4.1 Section-by-Section Results (Left Half; Right Half Mirrors by Symmetry)

| x(ft) | V_u(k) | V_ci(k) | V_cw(k) | V_c(k) | φV_c(k) | V_s,req(k) | A_v/s(in²/in) | s_req(in) | s_use(in) | Governs |
|-------|--------|---------|---------|--------|---------|------------|---------------|----------|----------|---------|
| 0.0 | 50.58 | >>V_cw | 56.88 | 56.88 | 42.66 | 10.56 | 0.00786 | 28.0 | 21 | Min A_v/s |
| 1.87 | 47.63 | 83.81 | 56.88 | 56.88 | 42.66 | 6.63 | 0.00726* | 30.3 | 21 | V_cw; min A_v/s |
| 5.30 | 42.21 | 34.41 | 56.88 | 34.41 | 25.81 | 21.87 | 0.01627 | 13.5 | 13 | V_ci; demand |
| 10.70 | 33.67 | 22.12* | 56.88 | 22.12 | 16.59 | 20.95 | 0.01559 | 14.1 | 14 | V_ci,min; demand |
| 16.00 | 25.29 | 22.12* | 56.88 | 22.12 | 16.59 | 11.60 | 0.00863 | 25.5 | 21 | V_ci,min; s_max |
| 21.33 | 16.90 | 22.12* | 54.34 | 22.12 | 16.59 | 0.40 | 0.00726* | 30.3 | 21 | Min A_v/s |
| 32.00 | 0.00 | 22.12* | 54.34 | 22.12 | 16.59 | 0.00 | 0.00726* | 30.3 | 21 | Min A_v/s |

\* Minimum V_ci,min = 22.12 kips or A_v/s,min = 0.00726 in²/in. governs. >> V_cw: V_ci → ∞ at support (M_max → 0); V_cw governs directly.

### 4.2 Stirrup Layout — #3 Closed U-Stirrups (A_v = 0.22 in², f_y = 60 ksi)

| Zone (each half from support) | x range | Spacing | Governing criterion |
|-------------------------------|---------|---------|-------------------|
| Support — critical section (V_cw zone) | 0 to 22 in. (0 to 1.87 ft) | s = 21 in. | Min. A_v/s; V_cw controls V_c |
| High shear — V_ci zone | 22 to ~80 in. (1.87 to 6.67 ft) | s = 13 in. | V_ci governs V_c; demand A_v/s |
| Transition — V_ci,min zone | ~80 to ~160 in. (6.67 to 13.3 ft) | s = 14 in. | V_ci,min; demand A_v/s |
| Mid-beam to midspan | ~160 in. to 32 ft | s = 21 in. | s_max limit or min. A_v/s |

### 4.3 Section Adequacy Checks

V_s,max (ACI 318-19 §22.5.10.1) = 8λ√f'c · b_w · d_p = 8 × 1.0 × 0.07746 × 7.5 × 22.40 = 104.2 kips

Maximum V_s,req = 22.77 kips << V_s,max = 104.2 kips → Section dimensions are adequate (no enlargement required). ✓

Tight spacing threshold = 52.17 kips > V_s,req,max = 22.77 kips → s_max = 21 in. (basic) applies everywhere. ✓

---

## References

1. ACI Committee 318. (2019). *Building Code Requirements for Structural Concrete (ACI 318-19) and Commentary.* American Concrete Institute, Farmington Hills, MI.
2. Naaman, A. E. (2004). *Prestressed Concrete Analysis and Design: Fundamentals* (2nd ed.). Techno Press 3000, Ann Arbor, MI.
3. CEE 530 Prestressed Concrete — Shear design notes and handouts. Arizona State University, Spring 2026.
