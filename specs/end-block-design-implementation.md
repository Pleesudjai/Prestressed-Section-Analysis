# End Block Design — Implementation Spec

**Created:** 2026-04-08
**Updated:** 2026-04-08 (verified against Naaman Ex. 4.17.3 PDF pp. 234–237 and NotebookLM)
**Theory reference:** `specs/end-block-design.md`
**Target folder:** `08_matlab_program_analysis/project_Project2/`

---

## What & Why

Create a standalone MATLAB function `endBlockDesign.m` that designs the end-block (anchorage zone) reinforcement for the Double-T prestressed beam using the **Gergely-Sozen linear elastic method** (CEE 530 course method).

The function reads the same 6-struct input from `inputData()` and produces:
- Net moment diagram on horizontal planes through the end block
- Bursting force T and required stirrup reinforcement
- Spalling reinforcement check
- Console output with full calculation detail
- Plots: (1) stress distribution at end face, (2) net moment diagram vs. y, (3) end block reinforcement sketch

---

## Inputs

All from existing `inputData()`:

| Struct | Fields used |
|--------|-----------|
| `beam` | `L` |
| `section` | `vertices`, `A`, `Ix`, `yc`, `yb`, `yt` |
| `materials` | `fci`, `fc`, `fy` (stirrup steel) |
| `prestress` | `tendons{i}.Aps`, `tendons{i}.fpi`, `tendons{i}.y` (y at x=0), `losses` |
| `reinforcement` | (not used directly) |
| `loads` | (not used — end block is a local D-region analysis) |

**Derived quantities:**
- `F_i` = total initial prestress force = SUM(Aps_i * fpi_i)
- `y_ps` = resultant tendon y-location at the support (x=0) = weighted average of tendon y-positions at x=0, weighted by each tendon's force (Aps_i * fpi_i)
- `e_0` = eccentricity at support = `yc - y_ps` (positive when tendon is below centroid)
- `h` = total section depth = max(vertices_y) - min(vertices_y)
- `St` = Ix / yt (top section modulus)
- `Sb` = Ix / yb (bottom section modulus)

---

## Outputs

### Console Output
```
END BLOCK DESIGN — GERGELY-SOZEN METHOD
========================================
Section Properties:
  h = XX.XX in, A = XXX.X in^2, Ix = XXXXX.X in^4
  yc = XX.XX in, St = XXXX.X in^3, Sb = XXXX.X in^3

Prestress at Support (x = 0):
  F_i = XXX.X kip (total initial, before losses)
  y_ps = XX.XX in (resultant tendon location)
  e_0 = XX.XX in (eccentricity below centroid)

Stress Distribution at End Face (P/A +/- Pe/S):
  f_top = +X.XXX ksi (compression)
  f_bot = -X.XXX ksi (tension)

Net Moment on Horizontal Planes:
  y (in)   |  M_concrete (kip-in)  |  M_prestress (kip-in)  |  M_net (kip-in)
  ---------|----------------------|------------------------|----------------
  ...      |  ...                 |  ...                   |  ...

Maximum Net Moment:
  M_max = XXXX.XX kip-in at y = XX.XX in

Bursting Force (Gergely-Sozen):
  T = M_max / (3h/4) = XXXX.XX / XX.XX = XX.XX kip

Stirrup Design:
  f_s = 20 ksi (WSD)
  A_s_req = T / f_s = XX.XX / 20 = X.XX in^2
  Using #X U-stirrups (A_v = X.XX in^2 per stirrup):
  n = ceil(A_s_req / A_v) = X stirrups
  Spacing = h / n = XX.X in
  Zone = 0 to h = XX in from end face

Spalling Reinforcement (ACI 318 Sec. 25.9.4.4.5):
  T_spalling = 0.02 * 1.2 * F_i = X.XX kip
  A_s_spalling = T_spalling / f_y = X.XX in^2
```

### Figures

**Figure 1 — End Block Analysis**
- 3-panel subplot:
  - (a) End face stress distribution: stress (ksi) vs. y (in) — vertical axis is y
  - (b) Net moment diagram: M_net (kip-in) vs. y (in) — shows peak location
  - (c) Reinforcement layout: schematic showing stirrup positions within h from end

### Return Struct

```matlab
eb.F_i          % Total initial prestress force (kip)
eb.y_ps         % Resultant tendon y at support (in)
eb.e_0          % Eccentricity at support (in)
eb.h            % Section depth (in)
eb.f_top        % Top fiber stress at end (ksi)
eb.f_bot        % Bottom fiber stress at end (ksi)
eb.y_net        % y-coordinates for net moment evaluation (in)
eb.M_net        % Net moment at each y (kip-in)
eb.M_max        % Maximum net moment (kip-in)
eb.y_Mmax       % y-location of M_max (in)
eb.T_burst      % Bursting force (kip)
eb.f_s          % Steel design stress (ksi, default 20 for WSD)
eb.As_req       % Required stirrup area (in^2)
eb.bar_size     % Selected bar size string (e.g., '#3')
eb.Av_bar       % Area per stirrup (in^2)
eb.n_stirrups   % Number of stirrups
eb.spacing      % Stirrup spacing (in)
eb.T_spalling   % Spalling force (kip)
eb.As_spalling  % Required spalling reinforcement area (in^2)
```

---

## Equations / Theory Reference

All from `specs/end-block-design.md`, Section 3 (Gergely-Sozen Method).
Verified against Naaman *Prestressed Concrete*, 2nd ed., Sec. 4.17, Example 4.17.3 (book pp. 201–204, PDF pp. 234–237).

### Naaman Reference Example (Ex. 4.17.3) — Verification Benchmark

Before describing the algorithm, the reference example is summarized so our implementation can be validated against it.

**Section:** Inverted T-beam (Naaman Fig. 4.9a)
- h = 40 in, A = 711 in², I = 102,916 in⁴
- Bottom flange: 18 in wide × 9 in tall
- Web: 6 in wide (some parts noted as b = 8 in due to adjacent transitions)
- Top flange: 60 in wide × 6 in tall (Note: b(y) at y = 36.5 in is wider than 8 in — Naaman Table 4.6 footnote)
- CGC: y_c = 27.1 in from bottom
- Z_t = 7868 in³, Z_b = 4490 in³

**Prestress:** F_i = 276.5 kip, e_0 = 7.9 in below CGC → y_ps = 19.2 in from bottom

**Elastic stresses at end face (linear bending):**
- f_bot (y=0): F_i/A + F_i*e_0*yc/Ix = 276.5/711 + 276.5*7.9*27.1/102916 = 0.389 + 0.575 = **+1.224 ksi** (compression)
- f(y=19.2): interpolated = **+0.713 ksi**
- f_top (y=40): F_i/A - F_i*e_0*yt/Ix = **+0.159 ksi** (small compression)

**Net moment table (Naaman Table 4.6):**

| y (in) | M_concrete (kip-in) | M_prestress (kip-in) | M_net (kip-in) |
|--------|--------------------|-----------------------|----------------|
| 4      | 76.06              | —                     | 76.07          |
| 8      | 295.17             | —                     | 295.17         |
| 12     | 643.78             | —                     | 643.78         |
| 16     | 1107.97            | —                     | 1107.97        |
| **19.2** | **1553.69**      | **0** (tendon on cut) | **1553.69**    |
| 24     | 2329.34            | -1327.2               | 1002.14        |
| 27.1   | 2889.66            | -2184.35              | 705.31         |
| 32     | 3850.24            | -3539.20              | 311.04         |
| 36.5   | 4831.62            | -4783.45              | 48.17          |
| 40     | 5714.5             | -5751.2               | ~0             |

**Key insight:** M_max = 1553.69 kip-in occurs at y = 19.2 in (the tendon level), NOT at the CGC. This is because below y_ps the prestress force does not contribute a moment (lever arm = 0), so M_net = M_concrete, which peaks just before the tendon.

**How M_concrete is computed at y = 19.2 in:**
The elastic stress at y = 19.2 is σ = 0.713 ksi, at y = 0 is σ = 1.224 ksi, and web width = 8 in.
The trapezoidal stress block is decomposed into a rectangle (σ_top * b) + triangle ((σ_bot − σ_top) * b / 2):

```
M_conc = [b/2 * (σ_bot - σ_top) * 2/3 + b * σ_top / 1] * y²
       = [8/2 * (1.224 - 0.713) * 2/3 + 4 * 0.713] * (19.2)²
       = [1.3627 + 2.852] * 368.64
       = 1553.69 kip-in
```

This closed-form works for a constant-width strip. For our Double-T (width changes at y=26), we use numerical integration instead.

**Bursting force and stirrups:**
```
T = M_max / (3h/4) = 1553.69 / 30 = 51.79 kip
A_s = T / f_s = 51.79 / 20 = 2.59 in²
→ 7 #4 closed stirrups at 6 in spacing (Fig. 4.32c)
```

---

### Step-by-step Algorithm (for arbitrary polygonal section)

**Step 1: Compute elastic stresses at the end face due to prestress only**

The stress varies linearly over the full depth (standard bending formula):
```
f(y) = F_i/A + F_i * e_0 * (y - yc) / Ix
```

Compression = positive (consistent with project sign convention).
- At top (y = h): `f_top = F_i/A - F_i*e_0*yt/Ix`
- At bot (y = 0): `f_bot = F_i/A + F_i*e_0*yb/Ix`

Note: e_0 is positive when tendon is below centroid. With e_0 > 0, bottom fiber gets more compression, top fiber gets less — physically correct for a low tendon.

**Step 2: Compute section width b(y) at each height**

For each horizontal strip at height y, intersect a horizontal line with the polygon edges to find the total width `b(y)`.

For the Double-T section:
- y = 0 to 26 in: b(y) ≈ 7.5 in (two stems, each ~3.75 in, with slight taper)
- y = 26 to 28 in: b(y) = 120 in (full flange width)

The `computeWidthAtY(vertices, y)` helper handles arbitrary polygons.

**Important:** Use fine y-spacing near the flange transition (y = 26 in) to capture the abrupt width change accurately.

**Step 3: Compute net moment on horizontal planes**

For a horizontal cut at height `y_cut`, consider the free body of the sub-region below the cut (Naaman Fig. 4.31, Sec. 4.17):

- **Concrete stress moment:** `M_concrete(y_cut) = integral from 0 to y_cut of f(y) * b(y) * (y_cut - y) dy`
  This is the moment of the stress resultant on the vertical cut face (at distance h from end), taken about the horizontal cut line.

- **Prestress moment:** `M_prestress(y_cut) = -F_i * max(0, y_cut - y_ps)`
  If the horizontal cut is above the tendon (y_cut > y_ps), the concentrated force F_i on the end face produces a moment about the cut line with lever arm = (y_cut - y_ps). The sign is negative because the prestress (compression, acting horizontally into the end face) creates a counterclockwise moment about the cut when below it.

- **Net moment:** `M_net(y_cut) = M_concrete(y_cut) + M_prestress(y_cut)`

Evaluate at a fine grid of y values from 0 to h (e.g., 200+ points, with refinement near y_ps and flange transitions).

**Numerical integration:** Use cumulative trapezoidal quadrature. For each y_cut, the integral `M_concrete` can be computed efficiently:
```matlab
% Fine y-grid
y_grid = linspace(0, h, N);
f_grid = F_i/A + F_i*e_0*(y_grid - yc)/Ix;   % stress at each y
b_grid = arrayfun(@(y) computeWidthAtY(vertices, y), y_grid);  % width at each y
q_grid = f_grid .* b_grid;  % force per unit height (kip/in)

% For each y_cut, M_concrete = integral_0^y_cut of q(y)*(y_cut - y) dy
% This can be computed as: y_cut * integral_0^y_cut(q dy) - integral_0^y_cut(q*y dy)
% Using cumtrapz for efficiency.
```

**Step 4: Find maximum net moment**

```
[M_max, idx] = max(abs(M_net));
y_Mmax = y_grid(idx);
```

This identifies the most likely location of a horizontal splitting crack. Typically occurs at or near y_ps (tendon level).

**Step 5: Compute bursting force**

```
T = M_max / (3*h/4)
```

(Gergely-Sozen: tension T in stirrups acts at h/4 from end face, compression C in concrete acts at h from end face, lever arm = 3h/4)

**Step 6: Size stirrups**

```
A_s = T / f_s       where f_s = 20 ksi (WSD, per Naaman/Mobasher)
n = ceil(A_s / A_v)  where A_v = 2 * A_bar (two legs of U-stirrup)
spacing = h / n
```

Default: #3 U-stirrups → A_bar = 0.11 in², A_v = 0.22 in² per stirrup.

**Step 7: Spalling check (ACI 318 Sec. 25.9.4.4.5)**

```
P_pu = 1.2 * F_i                   (factored for end-zone design)
T_spalling = 0.02 * P_pu           (when tendon centroid is within the kern)
A_s_spalling = T_spalling / f_y    (strength design, f_y = 60 ksi)
```

**Also print ACI strength-design comparison:**
```
T_burst_ACI = 0.25 * P_pu * (1 - h_anc/h)     (ACI 318 Sec. 25.9.4.4)
A_s_ACI = T_burst_ACI / (phi * f_y)            (phi = 0.75, f_y = 60 ksi)
```
where h_anc is estimated from the tendon spread at the end face.

---

## Implementation Steps

- [ ] 1. Create `endBlockDesign.m` in `project_Project2/` with function signature:
      `function [eb] = endBlockDesign(beam, section, materials, prestress)`
- [ ] 2. Implement helper `computeWidthAtY(vertices, y)` — polygon intersection to get total width at height y
      - Must handle the Double-T correctly: two stems at y < 26 in (total ~7.5 in), full flange at y >= 26 in (120 in)
      - Edge case: at exact vertex y-coordinates, use a small epsilon offset
- [ ] 3. Compute F_i, y_ps (force-weighted centroid of tendons at x=0), e_0
      - For custom-profile tendons, use `tendon.y(1)` (the y-value at the first x-position)
      - For linear tendons, use `tendon.y_position`
- [ ] 4. Compute elastic stress: `f(y) = F_i/A + F_i*e_0*(y - yc)/Ix` on a fine y-grid
      - Use 500+ points with refinement near y_ps, y=26 (flange transition), and y_c
- [ ] 5. Compute b(y) at each grid point via `computeWidthAtY`
- [ ] 6. Compute M_concrete(y_cut) via cumulative trapezoidal integration:
      ```
      Q(y_cut) = cumtrapz(y, q)        % cumulative force
      S(y_cut) = cumtrapz(y, q .* y)   % cumulative first moment
      M_concrete(y_cut) = y_cut * Q(y_cut) - S(y_cut)
      ```
- [ ] 7. Compute M_prestress(y_cut) = -F_i * max(0, y_cut - y_ps)
- [ ] 8. Compute M_net = M_concrete + M_prestress; find M_max and y_Mmax
- [ ] 9. Compute T_burst = M_max / (3*h/4)
- [ ] 10. Size stirrups: A_s = T/f_s, n = ceil(A_s/A_v), spacing = h/n
- [ ] 11. Compute spalling reinforcement per ACI 318
- [ ] 12. Print ACI approximate method comparison (T_burst_ACI = 0.25*P_pu*(1-h_anc/h))
- [ ] 13. Print full console output with table of M_net at key y-locations
- [ ] 14. Generate 3-panel figure:
      - (a) End face stress f(y) vs. y, with section outline overlay showing b(y)
      - (b) Net moment diagram M_net vs. y, with M_max marked
      - (c) Reinforcement layout: cross-section with stirrup positions within h from end
- [ ] 15. Add call to `endBlockDesign` in `main_PrestressedBeamAnalysis.m` (after stage loop, before report)

---

## Verification Plan

After implementation, verify against Naaman Ex. 4.17.3:
1. Create a temporary test with Naaman's inverted-T section (A=711, I=102916, yc=27.1, h=40)
2. Set F_i=276.5 kip, y_ps=19.2 in
3. Confirm M_net at y=19.2 = 1553.69 kip-in (within 1%)
4. Confirm T = 51.79 kip
5. Confirm A_s = 2.59 in²

---

## Open Questions — Resolved

1. **Bar size for stirrups:** Default #3 U-stirrups (matching HW assignment). Should this be configurable?
   → **Decision:** Default #3, but accept optional parameter to override.

2. **Steel stress f_s:** Naaman and Mobasher both use 20 ksi (WSD). ACI 318 uses phi*f_y = 0.75*60 = 45 ksi (strength design). Which to use?
   → **Decision:** Use f_s = 20 ksi (WSD) as primary (matches CEE 530 course). Also print the ACI strength-design result for comparison.

3. **Which prestress force?** F_i (initial, at transfer) is correct for end block design since the end block stresses are highest at jacking/transfer.
   → **Decision:** Use F_i = SUM(Aps * fpi), no losses applied.

4. **Where does M_max occur?** (Resolved by reading Naaman Ex. 4.17.3)
   → M_max occurs at or near y_ps (tendon height), NOT at the CGC. Below the tendon, M_prestress = 0, so M_net = M_concrete grows unchecked. Above the tendon, M_prestress subtracts and M_net drops. The peak is at y_ps where the prestress moment first kicks in with zero lever arm.

5. **How does width b(y) affect the integration?** (Resolved by Naaman Table 4.6 footnote)
   → The stress is linear (standard bending formula), but forces and moments depend on b(y). Naaman's footnote at y=36.5 in says "beam width at this level is larger than eight inches" — confirming the method accounts for actual section shape. For our Double-T, the flange transition at y=26 will cause a large jump in force intensity.

---

## Out of Scope

- Strut-and-tie model (STM) — different methodology, not required for this assignment
- Local zone bearing check — requires anchorage plate dimensions not defined in inputData
- Post-tensioned friction losses along the end block — this is a pre-tensioned beam
- 3D effects (biaxial splitting in double-T stems)
