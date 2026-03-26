# Stage 4: Ultimate Strength (Strength Limit State)

**Source:** ACI 318-19 §22.3–22.5, Naaman Ch. 5 — confirmed via NotebookLM + CEE530 course notes (03 CEE530_UltimateDesign.pdf)
**Sign convention:** strength design, all values positive magnitudes

---

## When Does This Occur?

Under factored (extreme) load combinations. The member must not **collapse**.
At ultimate, the initial compression from prestress essentially vanishes as the steel yields — the beam behaves similarly to a conventionally reinforced concrete beam.

---

## ACI 318-19 Factored Load Combinations (U)

The required strength U must satisfy the most critical of:

| Combination | Formula |
|-------------|---------|
| 1 | U = 1.4D |
| **2 (usually governs)** | **U = 1.2D + 1.6L + 0.5(Lr or S or R)** |
| 3 | U = 1.2D + 1.6(Lr or S or R) + (1.0L or 0.5W) |
| 4 | U = 1.2D + 1.0W + 1.0L + 0.5(Lr or S or R) |
| 5 | U = 1.2D + 1.0E + 1.0L + 0.2S |
| 6 | U = 0.9D + 1.0W |
| 7 | U = 0.9D + 1.0E |

**For Project 1 (gravity only):** `Mu = 1.2·M_D + 1.6·M_L`

---

## Stress in Prestressing Steel at Nominal Strength (fps)

ACI 318-19 approximate formula for **bonded tendons** (valid when fse ≥ 0.5·fpu):

```
fps = fpu · [1 − (γp/β1) · (ρp·fpu/f'c  +  (ds/dp)·(ωs − ω'))]
```

> ⚠️ **Bug fix (cross-checked vs CEE530 course notes):** The mild steel depth ratio is **(ds/dp)** — mild steel depth divided by strand depth — NOT (dp/ds). The two are usually close to 1.0 but the direction matters.

Where:
- `fpu` = ultimate tensile strength of prestressing steel
- `γp` = factor depending on strand type (see table below)
- `β1` = 0.65–0.85 depending on f'c (= 0.75 for f'c = 6 ksi)
- `ρp = Aps/(b·dp)` = prestressing steel ratio (use **web width bw** for T-sections)
- `ds` = effective depth of mild (non-prestressed) tension steel
- `dp` = effective depth of prestressing steel (strand centroid from compression face)
- `ωs = As·fy/(b·dp·f'c)`, `ω' = As'·fy'/(b·dp·f'c)` (mild steel terms)

### γp Table (ACI 318-19 Table 22.3.2.1)

| fpy/fpu | γp | Strand type |
|---------|-----|-------------|
| ≥ 0.90  | **0.28** | Low-relaxation (Grade 270 LR) ← **Project 1/2** |
| ≥ 0.85  | **0.40** | Stress-relieved (older strand, used in CEE530 course example) |
| ≥ 0.80  | 0.55 | Wire |

**Project 1/2:** fpu = 270 ksi, fpy = 243 ksi → fpy/fpu = 0.90 → **γp = 0.28**

**Simplified (prestress only, no mild steel):**
```
fps = fpu · [1 − (γp/β1) · ρp·fpu/f'c]
```

---

## Nominal Moment Capacity (Mn) — Two-Step Procedure

### Step 1: Try rectangular (assume a ≤ hf, compression in flange only)

Force equilibrium assuming full flange width b:
```
a_trial = (Aps·fps + As·fy − As'·fy') / (0.85·f'c·b)
```

**If a_trial ≤ hf → rectangular section applies:**
```
Mn = Aps·fps·(dp − a/2)  +  As·fy·(ds − a/2)
```

### Step 2: If a_trial > hf → T-section required (neutral axis in web)

> This is common for Double-T and I-sections. Recompute fps using **web-only** ρpw = Aps/(bw·dp), then:

```
a = (Aps·fps + As·fy − As'·fy' − 0.85·f'c·(b−bw)·hf) / (0.85·f'c·bw)
```

**Moment about centroid of compression zone** (NOT a/2 — must locate T-block centroid):
```
Mn = Aps·fps·(dp − ȳ)  +  As·fy·(ds − ȳ)
```
where ȳ = centroid of T-shaped compression block from top fiber:
```
ȳ = [(b−bw)·hf·(hf/2) + bw·a·(a/2)] / [(b−bw)·hf + bw·a]
```

> **CEE530 example (f'c = 5 ksi):** a_trial = 7.69 in > hf = 6.25 in → T-section. After T-section recompute: fps dropped from 242 ksi to 189.8 ksi, ȳ = 3.16 in, Mn = 16,901 kip-in. The T-section correction is significant — do not skip.

---

## ACI Design Check

```
φMn ≥ Mu
```

Where φ depends on the net tensile strain at the extreme tension steel (εt).

### Ductility Check (Net Tensile Strain):

From strain compatibility at ultimate (εcu = 0.003):
```
c  = a / β1
εt = (dt − c) / c · 0.003
```

| εt range | Section type | φ |
|----------|-------------|---|
| εt ≥ 0.005 | Tension-controlled | **φ = 0.90** |
| 0.004 ≤ εt < 0.005 | Transition (flexure OK, φ interpolated) | 0.65 < φ < 0.90 |
| εt < 0.004 | **Not permitted** for flexural members | — |
| εt ≤ 0.002 | Compression-controlled | φ = 0.65 |

> **From CEE530 course notes:** The **minimum allowable εt = 0.004** for flexural members (ACI 318-19 §9.3.3). The φ = 0.90 threshold is εt ≥ 0.005. Between 0.004 and 0.005 you still have a valid design but φ < 0.90.

**CEE530 example check:**
```
c = 6.95/0.80 = 8.06 in   (β1 = 0.80 for f'c = 5 ksi)
εt = (37.6 − 8.06)/8.06 × 0.003 = 0.0110 > 0.005  →  φ = 0.90  ✓
```

---

## Minimum Flexural Strength (Anti-Brittle Failure)

To prevent sudden brittle failure at first cracking:
```
φMn ≥ 1.2·Mcr
```

**Cracking moment** — two equivalent forms:

**ACI compact form:**
```
Mcr = (fr + fce) · Sb
```
Where:
- `fr = 7.5λ√f'c` (psi) = modulus of rupture (λ = 1.0 normal weight)
- `fce = Fe/Ac + Fe·e/Sb` = effective **concrete** compressive stress at bottom fiber (compression = positive, ksi → psi)
- `Sb = Ic/yb` = bottom section modulus

**Derivation form (CEE530 course notes):**

Set bottom fiber stress = rupture at moment Mcr (tension positive in this derivation):
```
fr = −fe_bot + Mcr/Sb
Mcr = fr·Sb + fe_bot·Sb = (fr + fce)·Sb      ← same result
```

**Naaman form:**
```
Mcr = Fe·(e0 + kt) + fr·Zb
```
where `kt = r²/yt` (upper kern), `Zb = Sb`. All three forms are equivalent.

> ⚠️ **Notation trap:** In Naaman, `fpe` usually means **steel** effective prestress (~150 ksi). ACI 318-19 uses `fpe` for the **concrete** stress at bottom fiber. Always use `fce` (concrete) to avoid confusion.

**Project 1/2 (f'c = 6000 psi):** fr = 7.5×√6000 = **+0.581 ksi**
**CEE530 example (f'c = 5000 psi):** fr = 7.5×√5000 = **+0.530 ksi**

---

## Summary Table

| Check | Requirement | φ factor |
|-------|-------------|----------|
| Flexural strength | φMn ≥ Mu | 0.90 (if εt ≥ 0.005) |
| Minimum strength | φMn ≥ 1.2·Mcr | 0.90 |
| Ductility floor | εt ≥ 0.004 (must satisfy) | — |
| Shear strength | φVn ≥ Vu | 0.75 |

---

## CEE530 vs ACI-318-19 Parameter Comparison (Exam Reference)

> This is the most common source of differences between textbook examples and ACI code calculations. Know which γp to use!

### β1 — Compression Block Factor (SAME in both)

| f'c (ksi) | β1 | Notes |
|-----------|-----|-------|
| 4 | 0.85 | — |
| 5 | 0.80 | CEE530 course example |
| **6** | **0.75** | **Project 1/2** |
| 7 | 0.70 | — |
| ≥ 8 | 0.65 (min) | lower bound |

Formula: `β1 = 0.85 − 0.05·(f'c − 4 ksi)`, clamped to [0.65, 0.85]

---

### γp — Prestress Factor (KEY DIFFERENCE)

| Parameter | CEE530 | ACI-318-19 |
|-----------|--------|------------|
| **fpy/fpu assumed** | **≥ 0.85** (stress-relieved) | **≥ 0.90** (low-relaxation) |
| **γp** | **0.40** | **0.28** |
| Strand type | Older stress-relieved | Modern low-relaxation Grade 270 |
| fpy for Project 1 (fpu=270) | 0.85×270 = **229.5 ksi** | 0.90×270 = **243 ksi** |

> **Why it matters:** Larger γp → smaller fps → smaller Mn. CEE530 uses γp=0.40 (conservative), ACI-318-19 uses γp=0.28 (higher fps, higher Mn for same strands).

**fps formula (prestress only, no mild steel):**
```
fps = fpu · [1 − (γp/β1) · ρp · fpu/f'c]
```

**Project 1 result (13 strands, Aps=1.989 in², dp≈22 in, f'c=6 ksi, b=120 in):**

| Quantity | CEE530 (γp=0.40) | ACI-318-19 (γp=0.28) |
|----------|-----------------|---------------------|
| ρp = Aps/(b·dp) | 0.000755 | 0.000755 |
| fps | ~261 ksi | ~265 ksi |
| a_trial = Aps·fps/(0.85·f'c·b) | ~0.845 in | ~0.862 in |
| Section type | **Rectangular** (a < hf=2 in) | **Rectangular** |
| Mn (kip-in) | ~11,225 | ~11,374 |
| φMn (kip-in) | ~10,103 | ~10,237 |
| Mu (kip-in) | ~9,712 | ~9,712 |
| D/C | ~0.961 | ~0.949 |
| **All checks** | **PASS** | **PASS** |

---

### Allowable Stresses — Service Stage

| Condition | CEE530 | ACI-318-19 |
|-----------|--------|------------|
| Compression (sustained, SW+SDL) | +0.45·f'c | +0.45·f'c (**SAME**) |
| Compression (total, SW+SDL+LL) | +0.60·f'c | +0.60·f'c (**SAME**) |
| Tension Class U boundary | −**6.0**√f'c (psi)/1000 | −**7.5**√f'c (psi)/1000 |
| Tension Class C | −12√f'c (psi)/1000 | −12√f'c (psi)/1000 (**SAME**) |

**f'c = 6 ksi:**
- Tension Class U: CEE530 = −0.465 ksi, ACI = −0.581 ksi (ACI is less conservative)
- Tension Class C: −0.930 ksi (same)

---

### Allowable Stresses — Transfer Stage

| Condition | CEE530 | ACI-318-19 |
|-----------|--------|------------|
| Compression (general) | +0.60·f'ci | +0.60·f'ci (**SAME**) |
| Compression (end zones) | +0.60·f'ci (same as general) | +**0.70**·f'ci |
| Tension (general) | −3√f'ci (psi)/1000 | −3√f'ci (psi)/1000 (**SAME**) |
| Tension (end zones) | −6√f'ci (psi)/1000 | −6√f'ci (psi)/1000 (**SAME**) |

**f'ci = 4.8 ksi:**
- Transfer compression (ends): CEE530 = +2.88 ksi, ACI = +3.36 ksi

---

### ductility / φ factor (SAME in both codes — ACI 318-19)

| εt range | Classification | φ |
|----------|---------------|---|
| εt ≥ 0.005 | Tension-controlled | **0.90** |
| 0.004 ≤ εt < 0.005 | Transition | 0.65 + (εt−0.002)/0.003 × 0.25 |
| εt < 0.004 | **NOT PERMITTED** (flexure) | — |
| εt ≤ 0.002 | Compression-controlled | 0.65 |

Minimum εt for flexural members: **0.004** (ACI §9.3.3)

---

### Cracking Moment Mcr (SAME formula)

```
Mcr = (fr + fce) · Sb
fr  = 7.5·λ·√f'c_psi / 1000   [ksi]   λ=1.0 normal weight
fce = Fe/Ac + Fe·e/Sb          [ksi]   (concrete compression at bottom, + = compression)
```

**f'c = 6 ksi:** fr = 7.5×√6000/1000 = **+0.581 ksi**

> ⚠️ `fce` uses the **effective** prestress force Fe = η·F (after losses), not the initial force.

---

## In MATLAB

```matlab
% --- fps: bonded tendons, prestress only (no mild steel) ---
gamma_p = 0.28;          % fpy/fpu >= 0.90 (low-relaxation, Project 1/2)
% gamma_p = 0.40;        % fpy/fpu >= 0.85 (stress-relieved, older strand)
beta1   = 0.85 - 0.05*(fc - 4);
beta1   = max(0.65, min(0.85, beta1));   % ACI limits

rho_p = Aps / (b * dp);
fps   = fpu * (1 - (gamma_p/beta1) * rho_p * fpu / fc);

% --- Step 1: try rectangular ---
a_trial = Aps * fps / (0.85 * fc * b);

if a_trial <= hf
    % Rectangular (a in flange)
    a  = a_trial;
    Mn = Aps * fps * (dp - a/2);
else
    % T-section: recompute with web-only rho_p
    rho_pw = Aps / (bw * dp);
    % also add mild steel terms if present
    fps = fpu * (1 - (gamma_p/beta1) * (rho_pw * fpu/fc));
    a   = (Aps*fps - 0.85*fc*(b-bw)*hf) / (0.85*fc*bw);
    % centroid of T compression block
    A_flange = (b-bw) * hf;
    A_web    = bw * a;
    ybar     = (A_flange*hf/2 + A_web*a/2) / (A_flange + A_web);
    Mn       = Aps * fps * (dp - ybar);
end

% --- Factored moment ---
Mu = 1.2 * M_dead + 1.6 * M_live;

% --- Ductility check ---
c     = a / beta1;
eps_t = (dp - c) / c * 0.003;    % use dp as extreme tension steel depth
if eps_t >= 0.005
    phi = 0.90;
elseif eps_t >= 0.004
    phi = 0.65 + (eps_t - 0.002) / (0.005 - 0.002) * (0.90 - 0.65);
else
    error('eps_t = %.4f < 0.004: section is over-reinforced (not permitted)', eps_t);
end

% --- Design checks ---
phi_Mn = phi * Mn;
fprintf('phi*Mn = %.1f kip-in  |  Mu = %.1f kip-in  |  D/C = %.2f\n', ...
    phi_Mn, Mu, Mu/phi_Mn);
assert(phi_Mn >= Mu, 'FAIL: phi*Mn < Mu');

% --- Minimum reinforcement ---
fr  = 7.5 * sqrt(fc * 1000) / 1000;   % ksi
fce = Fe/Ac + Fe*e/Sb;                 % ksi (compression at bottom fiber)
Mcr = (fr + fce) * Sb;                % kip-in
assert(phi_Mn >= 1.2*Mcr, 'FAIL: phi*Mn < 1.2*Mcr (minimum reinforcement)');
fprintf('1.2*Mcr = %.1f kip-in  -->  %s\n', 1.2*Mcr, ...
    ternary(phi_Mn >= 1.2*Mcr, 'OK', 'FAIL'));
```
