# Unit Conversions — SI and US Customary for Structural & Civil Engineering

**Fact-checked against:** NIST SP 1038, ASCE Unit Conversion Appendix, ACI 318-19, Naaman
**Date verified:** 2026-03-16
**Sign convention note:** All conversions are magnitude-only (positive values)

---

## Quick Reference Card

| Quantity | US → SI | SI → US |
|----------|---------|---------|
| 1 in | = 25.4 mm | 1 mm = 0.03937 in |
| 1 ft | = 0.3048 m | 1 m = 3.2808 ft |
| 1 lbf | = 4.4482 N | 1 kN = 224.81 lbf |
| 1 kip | = 4.4482 kN | 1 kN = 0.22481 kip |
| 1 psi | = 0.006895 MPa | 1 MPa = 145.04 psi |
| **1 ksi** | **= 6.8948 MPa** | **1 MPa = 0.14504 ksi** |
| 1 kip·ft | = 1.3558 kN·m | 1 kN·m = 0.73756 kip·ft |
| 1 kip·in | = 0.11298 kN·m | 1 kN·m = 8.8507 kip·in |
| 150 pcf | = 23.56 kN/m³ | — |

---

## 1. Length

| From | Exact Factor | To | Notes |
|------|-----------|----|-------|
| 1 in | × 25.4 | mm | **Exact** by international definition (1959) |
| 1 in | × 2.54 | cm | Exact |
| 1 in | × 0.0254 | m | Exact |
| 1 ft | × 0.3048 | m | Exact (= 12 × 0.0254) |
| 1 yd | × 0.9144 | m | Exact |
| 1 mile | × 1.60934 | km | |
| 1 mm | × 0.03937 | in | = 1/25.4 |
| 1 m | × 3.28084 | ft | = 1/0.3048 |
| 1 m | × 39.3701 | in | |
| 1 km | × 0.62137 | mile | |

**Structural use:** Beam spans, member dimensions, rebar spacing, cover depth

---

## 2. Force

| From | Exact Factor | To | Notes |
|------|-----------|----|-------|
| 1 lbf | × 4.44822 | N | = 0.45359 kg × 9.80665 m/s² |
| 1 lbf | × 0.00444822 | kN | |
| 1 kip (= 1000 lbf) | × 4.44822 | kN | **Most used in structural** |
| 1 kip | × 0.00444822 | MN | |
| 1 N | × 0.22481 | lbf | |
| 1 kN | × 0.22481 | kip | |
| 1 MN | × 224.81 | kip | |
| 1 ton (short, 2000 lb) | × 8.8964 | kN | |
| 1 tonne (metric, 1000 kg) | × 9.8067 | kN | |

**Structural use:** Reactions, tendon force, shear, axial load

---

## 3. Stress and Pressure

| From | Exact Factor | To | Notes |
|------|-----------|----|-------|
| 1 psi (lb/in²) | × 0.0068948 | MPa | = 4.44822 N / (25.4 mm)² |
| 1 psi | × 6.8948 | kPa | |
| **1 ksi (kip/in²)** | **× 6.8948** | **MPa** | **Key conversion for structural** |
| 1 ksi | × 0.0068948 | GPa | |
| 1 MPa | × 145.038 | psi | |
| 1 MPa | × 0.14504 | ksi | |
| 1 GPa | × 145.038 | ksi | |
| 1 kPa | × 0.14504 | psi | |
| 1 bar | × 14.5038 | psi | ≈ 0.1 MPa |
| 1 atm | × 14.696 | psi | = 101.325 kPa |

**Structural use:** Concrete strength f'c, steel yield fy, allowable stresses, soil pressure

### Common Material Properties — Both Units

| Property | US | SI |
|----------|----|----|
| Normal concrete f'c | 4–8 ksi | 28–55 MPa |
| Steel yield fy (mild) | 60 ksi | 414 MPa |
| Prestress strand fpu | 270 ksi | 1862 MPa |
| Concrete Ec | 4700√f'c (ksi) | 4700√f'c (MPa) ← same coeff! |
| Steel Es | 29,000 ksi | 200,000 MPa |
| Eps (strand) | 28,500 ksi | 196,500 MPa |

---

## 4. Moment and Torque

| From | Factor | To | Notes |
|------|--------|----|-------|
| 1 kip·in | × 0.11298 | kN·m | = 4.44822 kN × 0.0254 m |
| 1 kip·ft | × 1.35582 | kN·m | = 4.44822 kN × 0.3048 m |
| 1 lb·ft | × 0.0013558 | kN·m | |
| 1 lb·in | × 0.00011298 | kN·m | |
| 1 kN·m | × 8.8507 | kip·in | |
| 1 kN·m | × 0.73756 | kip·ft | |
| 1 N·m | × 0.73756 | lb·ft | |

**Structural use:** Bending moment, torsion, applied moment from loads

---

## 5. Distributed Load (Linear — per unit length)

| From | Factor | To | Notes |
|------|--------|----|-------|
| 1 lb/ft (plf) | × 0.014594 | kN/m | |
| 1 kip/ft (klf) | × 14.5939 | kN/m | **Most common in beam design** |
| 1 kip/in | × 175.127 | kN/m | |
| 1 lb/in | × 0.17513 | kN/m | |
| 1 kN/m | × 0.068522 | kip/ft | |
| 1 kN/m | × 0.005710 | kip/in | |
| 1 N/m | × 0.068522 | lb/ft | |

**Structural use:** Self-weight w (kip/in or kN/m), uniform load on beams

---

## 6. Area Load (Pressure — per unit area)

| From | Factor | To | Notes |
|------|--------|----|-------|
| 1 psf (lb/ft²) | × 0.047880 | kPa | |
| 1 ksf (kip/ft²) | × 47.880 | kPa | |
| 1 psi (lb/in²) | × 6.8948 | kPa | = same as stress |
| 1 kPa | × 20.885 | psf | |
| 1 kPa | × 0.020885 | ksf | |

**Structural use:** Floor live loads, snow loads, wind pressure, soil pressure

---

## 7. Area (Cross-Section)

| From | Factor | To | Notes |
|------|--------|----|-------|
| 1 in² | × 645.16 | mm² | = 25.4² — **Exact** |
| 1 in² | × 6.4516 | cm² | |
| 1 in² | × 6.4516×10⁻⁴ | m² | |
| 1 ft² | × 0.092903 | m² | = 0.3048² — Exact |
| 1 mm² | × 0.0015500 | in² | |
| 1 cm² | × 0.15500 | in² | |
| 1 m² | × 10.7639 | ft² | |
| 1 m² | × 1550.0 | in² | |

**Structural use:** Cross-sectional area Ac, steel area Aps, As, rebar area

---

## 8. Second Moment of Area (Moment of Inertia — I)

| From | Factor | To | Notes |
|------|--------|----|-------|
| 1 in⁴ | × 416,231 | mm⁴ | = 25.4⁴ |
| 1 in⁴ | × 41.6231 | cm⁴ | = 2.54⁴ |
| 1 cm⁴ | × 10,000 | mm⁴ | Exact |
| 1 mm⁴ | × 2.4025×10⁻⁶ | in⁴ | |

**Structural use:** Ic, Ig for deflection and stress calculations

---

## 9. Section Modulus (S or Z)

| From | Factor | To | Notes |
|------|--------|----|-------|
| 1 in³ | × 16,387.1 | mm³ | = 25.4³ |
| 1 in³ | × 16.3871 | cm³ | = 2.54³ |
| 1 cm³ | × 1,000 | mm³ | Exact |
| 1 mm³ | × 6.1024×10⁻⁵ | in³ | |

**Structural use:** St, Sb (section moduli for top/bottom fiber stress)

---

## 10. Unit Weight and Density

| From | Factor | To | Notes |
|------|--------|----|-------|
| 1 pcf (lb/ft³) | × 0.15709 | kN/m³ | |
| 1 pcf | × 16.0185 | kg/m³ | |
| **150 pcf** | **= 23.56 kN/m³** | | Normal weight concrete |
| 145 pcf | = 22.78 kN/m³ | | Lightweight concrete typical |
| 1 lb/in³ | × 271.45 | kN/m³ | |
| 1 kN/m³ | × 6.3658 | pcf | |
| 1 kg/m³ | × 0.062428 | pcf | |

**Structural use:** Self-weight calculation: w = γ × Ac

---

## 11. Temperature

| Conversion | Formula |
|-----------|---------|
| °F → °C | °C = (°F − 32) × 5/9 |
| °C → °F | °F = °C × 9/5 + 32 |
| °F → K | K = (°F + 459.67) × 5/9 |
| °C → K | K = °C + 273.15 |

| Reference Point | °F | °C |
|----------------|----|----|
| Concrete placement (min) | 50°F | 10°C |
| Standard curing | 68°F | 20°C |
| Hot weather concreting | 90°F | 32°C |
| Boiling water | 212°F | 100°C |

---

## 12. ACI 318 Code — √f'c Conversion (Critical!)

ACI 318-19 writes allowable stress formulas differently in US and SI:

| Limit | US Customary | SI Metric |
|-------|-------------|-----------|
| Modulus of rupture | fr = 7.5 λ√f'c (psi) | fr = 0.625 λ√f'c (MPa) |
| Transfer tension (general) | −3√f'ci (psi) | −0.25√f'ci (MPa) |
| Transfer tension (ends) | −6√f'ci (psi) | −0.50√f'ci (MPa) |
| Service tension Class C | −12√f'c (psi) | −1.0√f'c (MPa) |
| Service compression | 0.45 f'c (ksi) | 0.45 f'c (MPa) |

**Why the coefficients differ:** The factor between US and SI coefficients = √(1/0.006895) = **12.01 ≈ 12**

```
fr_US (psi) = 7.5 × √f'c_psi
fr_SI (MPa) = 0.625 × √f'c_MPa

Check: 7.5 / 0.625 = 12.0 ✓   (matches the √(psi/MPa) conversion factor)
```

---

## 13. Velocity and Flow (Civil Engineering)

| From | Factor | To |
|------|--------|----|
| 1 ft/s | × 0.3048 | m/s |
| 1 mph | × 1.60934 | km/h |
| 1 mph | × 0.44704 | m/s |
| 1 ft³/s (cfs) | × 0.028317 | m³/s |
| 1 gallon/min (gpm) | × 6.3090×10⁻⁵ | m³/s |
| 1 MGD (million gal/day) | × 0.043813 | m³/s |
| 1 m³/s | × 35.3147 | cfs |

**Civil use:** Hydraulics, drainage design, stormwater

---

## 14. Soil and Geotechnical

| Quantity | US | SI |
|----------|----|----|
| Soil unit weight (typical) | 110–130 pcf | 17.3–20.4 kN/m³ |
| Saturated soil | ~125 pcf | ~19.6 kN/m³ |
| Water unit weight | 62.4 pcf | 9.81 kN/m³ |
| Bearing capacity | ksf or psi | kPa or MPa |
| 1 ton/ft² (tsf) | × 95.76 | kPa | |
| 1 kPa | × 0.01044 | tsf | |

---

## 15. MATLAB Unit Conversion Helper

```matlab
% ============================================================
% UNIT CONVERSION CONSTANTS — fact-checked 2026-03-16
% ============================================================

% Length
in_to_mm  = 25.4;          % exact
ft_to_m   = 0.3048;        % exact
mm_to_in  = 1/25.4;
m_to_ft   = 1/0.3048;

% Force
lbf_to_N  = 4.44822;
kip_to_kN = 4.44822;
kN_to_kip = 1/4.44822;

% Stress
psi_to_MPa = 0.0068948;
ksi_to_MPa = 6.8948;       % = 1000 × psi_to_MPa
MPa_to_ksi = 1/6.8948;     % = 0.14504

% Moment
kipin_to_kNm = 0.11298;
kipft_to_kNm = 1.35582;
kNm_to_kipft = 1/1.35582;

% Area
in2_to_mm2 = 645.16;       % = 25.4^2, exact
mm2_to_in2 = 1/645.16;

% Inertia
in4_to_mm4 = 25.4^4;       % = 416231.4
mm4_to_in4 = 1/25.4^4;

% Section modulus
in3_to_mm3 = 25.4^3;       % = 16387.1
mm3_to_in3 = 1/25.4^3;

% Unit weight
pcf_to_kNm3 = 0.15709;
kNm3_to_pcf = 1/0.15709;

% ACI sqrt(f'c) — note: use psi input for US, MPa input for SI
% fr_US = 7.5 * lambda * sqrt(fc_psi)   [result in psi]
% fr_SI = 0.625 * lambda * sqrt(fc_MPa) [result in MPa]
% Ratio: 7.5 / 0.625 = 12.0 (= sqrt(psi/MPa conversion))
```

---

---

## 16. Load Types — Physical Meaning and Units

Civil and structural engineers describe loads by how they are **distributed in space**. The type determines the unit and how it is used in equilibrium equations.

---

### Point Load (Concentrated Load) — P

**Physical meaning:** Force applied at a single point. In reality all loads have some contact area, but when that area is small compared to the structure, it is idealized as a point.

```
Dimension: [Force]
US units:  kip (k),  lb
SI units:  kN,  N
```

**Examples:**
- Column reaction on a footing
- Wheel load from a truck axle on a bridge girder
- Prestressing tendon force at anchorage
- Reaction at a beam support

**In equilibrium:**
```
ΣF = 0    (sum of point forces)
ΣM = 0    (moment = P × distance)
```

**In MATLAB:**
```matlab
P = 100;   % kip   (point load at midspan)
M_mid = P * L / 4;  % kip·in  (simply supported, midspan moment)
```

---

### Line Load (Distributed Load per Unit Length) — w

**Physical meaning:** Force spread uniformly (or varying) along the length of a member. Arises from weight of slabs on beams, soil pressure on walls, wind on building facades.

```
Dimension: [Force / Length]
US units:  kip/ft (klf),  kip/in,  lb/ft (plf)
SI units:  kN/m,  N/m
```

**To convert to equivalent point load:**
```
P_equiv = w × L       (total resultant, acts at centroid of distribution)
```

**For uniform w on simply supported beam of span L:**
```
R = w·L/2             (reaction at each support)
M_max = w·L²/8        (midspan moment)
V_max = w·L/2         (shear at support)
```

**Examples:**
- Self-weight: w_sw = γ_c × Ac  [kip/in = (kip/in³) × in²]
  - Project 1: w_sw = (150/1728) × Ac  (ρ in lb/in³ × area in in²)
- SDL from topping slab: w_SDL = γ_c × t_slab × b_flange
- Live load on floor: LL (psf) × tributary width = kip/ft

**Varying (trapezoidal or triangular):**
```
w(x) = w1 + (w2 - w1) × x/L    (linearly varying from w1 to w2)
```

**In MATLAB:**
```matlab
w = 0.05;          % kip/in  (uniform distributed load)
x = linspace(0, L, 1000);
M = w .* x .* (L - x) / 2;    % moment at each x (simply supported)
```

---

### Area Load (Surface Load) — q

**Physical meaning:** Force applied over a surface area. Used for floor loads, soil pressure on slabs, wind/snow on roofs.

```
Dimension: [Force / Area]
US units:  psf (lb/ft²),  ksf (kip/ft²),  psi (lb/in²)
SI units:  kPa (kN/m²),  MPa (N/mm²)
```

**To convert area load → line load on a beam:**
```
w (kip/ft) = q (ksf) × tributary width (ft)
```

**Examples:**
- Floor live load: 40–100 psf depending on occupancy (ASCE 7)
- Snow load: 20–40 psf depending on ground snow
- Soil pressure on basement wall: q = γ_soil × H (increases with depth)
- Prestress bearing stress: q = P / (bearing area)

**Common ACI/ASCE values:**
| Load Type | Typical US | Typical SI |
|-----------|-----------|-----------|
| Office live load | 50 psf | 2.4 kPa |
| Residential live load | 40 psf | 1.9 kPa |
| Parking garage | 40 psf | 1.9 kPa |
| Roof live load | 20 psf | 1.0 kPa |
| 2" concrete topping | 25 psf | 1.2 kPa |

---

### Body Load / Volume Load — γ (Unit Weight)

**Physical meaning:** Force generated by gravity acting through the entire volume of a material. This is how self-weight is computed.

```
Dimension: [Force / Volume]
US units:  lb/ft³ (pcf),  lb/in³,  kip/ft³
SI units:  kN/m³,  N/m³
```

**Self-weight calculation:**
```
w (line load) = γ (unit weight) × A_cross-section

w [kip/in] = γ [kip/in³] × Ac [in²]
w [kN/m]   = γ [kN/m³]  × Ac [m²]
```

**Common material unit weights (fact-checked):**

| Material | US (pcf) | SI (kN/m³) |
|----------|----------|-----------|
| Normal weight concrete | **150** | **23.56** |
| Lightweight concrete | 90–115 | 14.1–18.1 |
| Sand / gravel | 100–120 | 15.7–18.9 |
| Steel | 490 | 77.0 |
| Water | 62.4 | 9.81 |
| Wood (Douglas fir) | 34 | 5.3 |
| Brick masonry | 120 | 18.9 |
| Soil (dry) | 85–100 | 13.3–15.7 |
| Soil (saturated) | 110–130 | 17.3–20.4 |

**Project 1 (Double-T self-weight):**
```matlab
gamma_c = 150/1728;   % lb/in³ → kip/in³ = 0.08681 kip/in³
w_sw = gamma_c * Ac;  % kip/in  (Ac in in²)
```

---

### Moment Load — M

**Physical meaning:** A couple applied at a point — pure rotation without net force. Also arises from eccentric point loads.

```
Dimension: [Force × Length]
US units:  kip·in,  kip·ft,  lb·ft
SI units:  kN·m,  N·m
```

**Examples:**
- Prestress moment: M_ps = F × e  (force × eccentricity)
- Eccentric column load on footing
- Applied end moment in a frame analysis

---

### How Load Types Relate — Hierarchy

```
Body load  (γ, kN/m³)
    × cross-section area (m²)
    ─────────────────────────
    = Line load  (w, kN/m)
         × tributary width
    ─────────────────────────
    equivalent from Area load  (q, kN/m²)
         × length or area
    ─────────────────────────
    = Point load  (P, kN)  ← resultant of any distributed load
         × moment arm
    ─────────────────────────
    = Moment  (M, kN·m)
```

**Practical conversion example (Project 1 SDL):**
```
SDL = 2 in topping × 120 in flange width × (150 lb/ft³ ÷ 1728 in³/ft³)
    = 2 × 120 × 0.08681 kip/in³ × in²
    = 20.83 lb/in = 0.02083 kip/in

Check in SI:
q_SDL = 2 in × (1/12 ft/in) × 150 pcf = 25 psf = 1.197 kPa
b_flange = 120 in = 10 ft
w_SDL = 1.197 kPa × (10 ft × 0.3048) = 1.197 × 3.048 = 3.65 kN/m
Convert: 0.02083 kip/in × 175.127 = 3.65 kN/m  ✓
```

---

## Sources (Fact-Checked)

| Source | Reference |
|--------|-----------|
| NIST SP 1038 | SI Units and Conversion Factors |
| ASCE | Appendix — Conversion Factors, Customary to Metric |
| ACI 318-19 | Appendix on SI conversions |
| Naaman 2nd ed | Unit conventions throughout |
| Engineering Toolbox | Area moment of inertia converter |
