# Stress Limit Comparison — ACI 318 Transfer vs. Service

## Sign Convention: Compression = Positive, Tension = Negative

---

## Allowable Stress Table

| Condition | Location | Type | Allowable | Formula | Value (ksi) |
|-----------|----------|------|-----------|---------|-------------|
| Transfer | Top fiber | Compression | f_ci_allow | 0.60 × f'ci | +2.88 ksi |
| Transfer | Bottom fiber | Tension | -f_ti_allow | -3√f'ci (ksi) | -0.208 ksi |
| Service | Top fiber | Compression | f_cs_allow | 0.45 × f'c | +2.70 ksi |
| Service | Bottom fiber | Tension | -f_ts_allow | -12√f'c (psi) × (1/1000) | -0.928 ksi |

Notes:
- f'ci = 4.8 ksi (strength at transfer)
- f'c = 6.0 ksi (28-day strength)
- Class C member (tension allowed at service)
- √f'ci in ksi units for tension formula

---

## Section Properties (computed from shoelace formula)

| Property | Symbol | Value | Units |
|----------|--------|-------|-------|
| Cross-sectional area | Ac | (computed) | in² |
| Centroid from bottom | yb | (computed) | in |
| Centroid from top | yt | (computed) | in |
| Moment of inertia | Ic | (computed) | in⁴ |
| Section modulus (top) | St = Ic/yt | (computed) | in³ |
| Section modulus (bottom) | Sb = Ic/yb | (computed) | in³ |
| Kern top | k_t = r²/yt | (computed) | in |
| Kern bottom | k_b = r²/yb | (computed) | in |

Fill in from MATLAB console output when `main_PrestressedBeamAnalysis.m` is run.

---

## Stress Check at Midspan (x = 384 in)

| Condition | Fiber | Computed Stress | Allowable | Status |
|-----------|-------|----------------|-----------|--------|
| Transfer | Top | | +2.88 ksi | |
| Transfer | Bottom | | -0.208 ksi | |
| Service (full) | Top | | +2.70 ksi | |
| Service (full) | Bottom | | -0.928 ksi | |

Fill in from MATLAB console output.
