# End Block Design — Specification

**Created:** 2026-04-08
**Course:** CEE 530 Prestressed Concrete, ASU Spring 2026
**References:** ACI 318-19 Sec. 25.9; Naaman Ch. 4 (Sec. 4.17, Ex. 4.17.3); Wight Ch. 17 (STM); Prof. Mobasher lecture notes (S 2020)
**Fact-check scope:** Checked against the files in `01_sources/End Block` only. Material cited to ACI 318-19, Naaman, Wight, or AASHTO but not visible in that folder is kept as a study note, but was not independently verified here.

---

## 1. Background and Terminology

The **end block** (anchorage zone) is the D-region at the end of a prestressed member where the concentrated prestressing force spreads to a uniform stress distribution. Saint-Venant's principle: the general zone extends approximately one member depth **h** from the loaded face.

### 1.1 Local Zone vs. General Zone

| Zone | Definition | Design responsibility |
|------|-----------|----------------------|
| **Local zone** | Rectangular prism immediately surrounding the anchorage device and confining reinforcement | Bearing check; often tested by anchorage supplier per ACI 423.7 / AASHTO LRFD |
| **General zone** | Includes local zone + the broader region where force fans out to full section; extends ~h from end face | Engineer designs bursting + spalling reinforcement |

### 1.2 Types of Tensile Forces in the General Zone

| Force | Location | Cause |
|-------|----------|-------|
| **Bursting (T_burst)** | Ahead of anchor, perpendicular to tendon axis; peak at ~0.3h from end | Transverse tension as compressive trajectories spread from narrow anchor to full section depth |
| **Spalling** | Near loaded face, away from anchor | Compatibility of deformations near member edge |

---

## 2. ACI 318-19 Provisions (Section 25.9)

*Fact-check note:* The provided class reference folder does not include the ACI 318-19 text, so Section 2 should be confirmed directly against the code before it is used as a design basis.

### 2.1 Factored Prestress Force

The factored prestressing force for anchorage zone design:

```
P_pu = 1.2 * P_jack

but not exceeding:
  P_pu <= 1.2 * 0.94 * f_py * A_ps
  P_pu <= 1.2 * 0.80 * f_pu * A_ps
```

### 2.2 Approximate Method — Bursting Force (ACI 318-19 Sec. 25.9.4.4)

When the general zone contains no severely curved tendons or massive discontinuities:

```
T_burst = 0.25 * SUM(P_pu) * (1 - h_anc / h)

d_burst = 0.5 * (h - e_anc)
```

**Variables:**
| Symbol | Definition |
|--------|-----------|
| `SUM(P_pu)` | Sum of factored prestressing forces from individual tendons |
| `h_anc` | Depth of anchorage device (or group of closely spaced devices) in the direction considered |
| `h` | Total depth of the cross-section |
| `e_anc` | Eccentricity of anchorage device (or group) from section centroid; **always positive** |
| `d_burst` | Distance from anchor face to centroid of bursting force |

**Grouping rule:** Anchors with center-to-center spacing <= 1.5 * (width of anchorage device) are treated as a single group.

### 2.3 Bursting Reinforcement Design

```
A_s = T_burst / (phi * f_y)
```

where `phi = 0.75` (strut-and-tie provisions) and `f_y` is the yield strength of nonprestressed reinforcement.

- Concrete tensile strength is **neglected** (ACI 318 Sec. 25.9.4.1)
- Reinforcement distributed within the general zone, centered around `d_burst`
- Nominal compressive stress in general zone concrete limited to `0.7 * lambda * f'ci`

### 2.4 Spalling Reinforcement

When the tendon centroid lies within the kern of the section:

```
T_spalling >= 0.02 * SUM(P_pu)
```

Provide reinforcement near the end face to resist this force.

### 2.5 Local Zone Bearing Check

Without special confinement:

```
f_bi = 0.5 * f'ci * sqrt(A2 / A1)  <=  1.0 * f'ci
```

| Symbol | Definition |
|--------|-----------|
| `A1` | Gross bearing area of anchor plate |
| `A2` | Maximum area of the supporting surface geometrically similar to and concentric with A1 (frustum with 1V:2H side slopes) |

With confinement: must satisfy ACI 423.7 or AASHTO LRFD acceptance testing.

### 2.6 Pre-tensioned Members (AASHTO Provision)

ACI 318 does not provide explicit formulas for pre-tensioned end zones. AASHTO specifies:

```
Vertical stirrups at f_s = 20 ksi to resist >= 0.04 * F_i
placed within d_p/4 from beam end
```

where `F_i` = total initial prestressing force, `d_p` = depth to prestressing steel.

---

## 3. Gergely-Sozen Method (Linear Elastic / Free-Body Approach)

This is the method taught in CEE 530 (Prof. Mobasher). It uses equilibrium of a free-body diagram of the end block.

### 3.1 Procedure

1. **Cut the end block** at a distance `h` from the loaded face
2. On the cut face: apply the **elastic bending stresses** due to prestress (P/A +/- Pe*y/I)
3. On the end face: apply the **concentrated prestress force** at its eccentricity
4. **Neglect** the vertical component of the prestressing force
5. Compute **net moment on horizontal planes** at multiple y-locations by summing:
   - Moment from concrete stress resultants above the cut plane
   - Moment from prestress force (if the cut plane is below the tendon)
6. The **maximum net moment** identifies the location of the potential horizontal splitting crack (typically near the CGC)
7. Equate the maximum net moment to the **internal couple**:
   - `T` (tension in stirrups) acts at `h/4` from the end face
   - `C` (compression in concrete) acts at `h` from the end face
   - Lever arm = `h - h/4 = 3h/4`

```
M_max = T * (3h/4)

Therefore:  T = M_max / (3h/4) = 4 * M_max / (3h)
```

8. Size stirrups using the adopted steel stress level. In Naaman's working-stress example, `A_s = T / f_s` with `f_s = 20 ksi`.

### 3.2 Worked Example (Naaman Example 4.17.3)

*Fact-check note:* Verified against `Naaman_prestressedconcrete.pdf`, Chapter 4, Sec. 4.17.3 (book pp. 202-204; PDF pp. 235-237).

**Given:**
- Beam depth: h = 40 in
- CGC location: y_c = 27.1 in from bottom
- Prestress force: F_i = 276.5 kip
- Tendon location at the support: y_p = 19.2 in from bottom
- Eccentricity: e_0 = 7.9 in below CGC
- End-zone length taken as h = 40 in
- Prestress force assumed to have constant eccentricity over the end zone; vertical component neglected
- Allowable stirrup stress used in the example: f_s = 20 ksi

**Net moments at horizontal planes (y from bottom):**

| y (in) | M_net (kip-in) |
|--------|----------------|
| 4 | 76.07 |
| 8 | 295.17 |
| 12 | 643.78 |
| 16 | 1107.97 |
| 19.2 | **1553.69** (max) |
| 24 | 1002.14 |
| 27.1 (CGC) | 705.31 |
| 32 | 311.04 |
| 36.5 | 48.17 |
| 40 | ~0 |

At y = 27.1 in (CGC):
- Moment from concrete stresses = 2889.66 kip-in
- Moment from prestress force = -F_i * (27.1 - 19.2) = -276.5 * 7.9 = -2184.35 kip-in
- **Net moment = 2889.66 - 2184.35 = 705.31 kip-in**

The **governing** net moment is larger and occurs at `y = 19.2 in`:
```text
M_max = 1553.69 kip-in
```

**Bursting force:**
```
T = M_max / (3h/4) = 1553.69 / (3*40/4) = 1553.69 / 30 = 51.79 kip
```

**Stirrup design:**
```
A_s = T / f_s = 51.79 / 20 = 2.59 in^2
```

**Selected:** 7 #4 closed stirrups at 6 in spacing (covers ~42 in, approximately h)

### 3.3 CEE 530 Lecture Notes — Deep Beam Analogy

From Prof. Mobasher's End Block Design notes (S 2020), the end block is also analyzed as a **deep beam** loaded by the prestress force:

**Given (lecture example):**
- I-section: I = 28,410 in^4, A = 315 in^2
- S_t = 1627 in^3, S_b = 2265 in^3
- e = 8.043 in
- P = 225.8 kip
- h = 30 in

**Stress distribution at end:**
- Top flange: +400 psi -> +4800 lb/in
- Other values written on the note:
  - y = 27 in: +208.2 psi -> +2498 lb/in
  - y = 24 in: +16.4 psi -> +98 lb/in
  - y = 12 in: -750 psi -> -4505 lb/in
  - y = 6 in: -1134 psi -> -20419 lb/in
  - y = 0 in: -1518 psi -> -27324 lb/in

**Bending moment diagram of end block:**
- Computed net moments on horizontal planes
- Moment magnitude called out near the zero-shear location = 165 kip-in
- Moment magnitude called out at 30 in from anchorage = 250 kip-in

**Deep beam design:**
```
z = 0.2h    (distance used on the sketch, so h - z = 0.8h)
f_s = 20,000 psi    (steel stress, working stress design)
T = M / (h - z)    (tensile force in steel)
```

**Steel near the anchorage:**
```
T = 165 / (0.8 * 30) = 6.875 kip
n = 6875 / (2 * 0.11 * 20000) = 1.56  -->  Use 2 or 3 #3 U-stirrups
```

**Steel at 30 in from anchorage:**
```
T = 250 / (0.8 * 30) = 10.4 kip
n = 10400 / (2 * 0.11 * 20000) = 2.36  -->  Use 3 or 4 #3 U-stirrups
```

---

## 4. Strut-and-Tie Model (STM) Approach

*Fact-check note:* The provided class reference folder does not include Wight Chapter 17 or ACI Chapter 23 text, so Section 4 was not independently verified here.

### 4.1 Concept

The STM treats the end block D-region as a pin-jointed truss:
- **Struts:** Concrete compressive stress fields radiating from anchorage
- **Ties:** Transverse reinforcement resisting bursting/spalling tension
- **Nodal zones:** Confined regions where struts and ties meet

### 4.2 Node Classifications (ACI 318 Ch. 23)

| Node type | Description | Strength coefficient beta_n |
|-----------|-------------|---------------------------|
| C-C-C | All compression | 1.0 |
| C-C-T | Two compression + one tension | 0.80 |
| C-T-T | One compression + two tension | 0.60 |

**Effective compressive strength of node:**
```
f_ce = 0.85 * beta_n * f'c
```

**Effective compressive strength of strut** (with crack-control reinforcement):
```
f_ce = 0.85 * beta_s * f'c     where beta_s = 0.75
```

### 4.3 Design Procedure

1. Determine the factored load P_pu = 1.2 * P_jack
2. Sketch the STM geometry: anchorage -> diagonal struts -> transverse tie
3. Compute tie force from equilibrium (geometry of the truss)
4. Size reinforcement: `A_s = T_tie / (phi * f_y)`
5. Check strut capacities: `phi * f_ce * A_strut >= F_strut`
6. Check nodal zone dimensions: ensure bearing area sufficient
7. Verify anchorage of tie reinforcement (development length for hooks)

---

## 5. HW Assignment — I-Beam End Block

From the homework assignment image:

**Problem:** Design the end block using the linear elastic method.

**Given directly on the image:**
- Top and bottom flange width = 18 in
- Top vertical dimensions shown = 4.5 in and 3.5 in
- Web region is marked with three 6 in horizontal segments
- Bottom vertical dimensions shown = 6 in and 4 in
- Prestressing force: P = 376,110 lb
- Force location: point P is 6.35 in above the bottom face

**Required:**
1. Draw the bending moment diagram (of the end block free body)
2. Indicate the amount and position of #3 stirrups

---

## 6. Project #2 Context

From the PROJ.png image - the sketch shows one symmetric half of a post-tensioned member with a parabolic tendon profile:

**Section properties written on the image:** A = 5904 in^2, I = 3,704,918 in^4
- y_b = 36.561 in, y_t = -29.439 in
- S_b = 101,335 in^3, S_t = -125,850 in^3
- k_b = 21.36 in, k_t = -17.16 in

**Half-model geometry shown:**
- Cross-section half-width to the symmetry line: 48 + 9 + 52 + 9 + 26 = 144 in
- Total depth shown: 6 + 3 + 48 + 9 = 66 in
- Cable-profile segment lengths: AB = 48 ft, BC = 57 ft, CD = 15 ft
- Total length from A to the symmetry line = 120 ft
- Cover callouts shown on the sketch: 4 in near B and 4 in near D

**Tasks:**
1. Draw cross-section and longitudinal section
2. Compute section parameters (A, I, y_c, y_b, k_b, S_b, S_t)
3. Give equations of parabolic segments AB, BC, CD
4. For F = 100 kip, compute equivalent loads for segments AB, BC, CD
5. Compute and draw bending moment diagram due to equivalent loads
6. For f_i = 189 ksi, compute stresses at sections A, B, C, D, C', B', A' and draw stress diagrams (mu = 0.15, k = 7.5E-4 rad/ft)
7. Use exact exponential formula for frictional losses

---

## 7. Summary of Design Methods

| Method | When to use | Key formula | Reference |
|--------|-------------|-------------|-----------|
| **ACI Approximate** | Post-tensioned, simple geometry | T_burst = 0.25*P_pu*(1 - h_anc/h) | ACI 318 Sec. 25.9.4.4 |
| **Gergely-Sozen** | Any member; hand-calc approach | T = M_max / (3h/4) from free-body | Naaman Sec. 4.17 |
| **Deep Beam Analogy** | CEE 530 course method | T = M / (0.8h), f_s = 20 ksi (WSD) | Prof. Mobasher notes |
| **Strut-and-Tie** | Complex geometry, multiple anchors | Truss equilibrium | ACI 318 Ch. 23, Wight Ch. 17 |
| **Pre-tensioned (AASHTO)** | Pre-tensioned beams | Stirrups for >= 0.04*F_i within d_p/4 | AASHTO LRFD |

---

## 8. Key Formulas Quick Reference

```
--- FACTORED LOAD ---
P_pu = 1.2 * P_jack

--- ACI APPROXIMATE ---
T_burst = 0.25 * SUM(P_pu) * (1 - h_anc/h)
d_burst = 0.5 * (h - e_anc)
A_s     = T_burst / (phi * f_y)

--- GERGELY-SOZEN ---
T = 4 * M_max / (3 * h)
  where M_max = max net moment on horizontal planes through end block

--- DEEP BEAM (CEE 530) ---
T = M / (h - z),   z = 0.2h
A_s = T / f_s,   f_s = 20 ksi (WSD)

--- BEARING ---
f_bi = 0.5 * f'ci * sqrt(A2/A1) <= 1.0 * f'ci

--- SPALLING ---
T_spalling >= 0.02 * SUM(P_pu)   (when tendon centroid within kern)

--- PRE-TENSIONED (AASHTO) ---
Stirrups for >= 0.04 * F_i within d_p/4 from end, at f_s = 20 ksi
```
