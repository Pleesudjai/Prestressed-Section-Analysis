# Domain Rules — Prestressed Concrete Theory

Load this file when working on stress analysis, Magnel diagrams, feasibility charts, or ACI 318 code checks.

## Reference Sources (priority order)
1. NotebookLM: `cee-530-prestressed-concrete` (ACI 318-19, Naaman 2nd ed, Wight 8th ed)
2. `01_sources/Naaman_prestressedconcrete.pdf` — Naaman textbook (direct read)
3. `01_sources/Pretressed Concrete (CEE 530) spring 2024.pdf` — Lecture slides

## Fiber Stress Equations (compression = +)

```
AT TRANSFER (force F, moment Mi = self-weight only):
  f_top = F/Ac  −  F·e/St  +  Mi/St
  f_bot = F/Ac  +  F·e/Sb  −  Mi/Sb

AT SERVICE (force η·F, moment MT):
  f_top = η·F/Ac  −  η·F·e/St  +  MT/St
  f_bot = η·F/Ac  +  η·F·e/Sb  −  MT/Sb
```

## ACI 318 Allowable Stresses

| Stage | Compression (+) | Tension (−) | ACI Section |
|-------|-----------------|-------------|-------------|
| Transfer | 0.60·f'ci | −3√f'ci ksi | §24.5.3 |
| Service (sustained) | 0.45·f'c | −12√f'c ksi (Class C) | §24.5.4 |
| Service (total) | 0.60·f'c | −7.5√f'c ksi (Class U) | §24.5.4 |

## Magnel Lines (1/F vs e at critical section)

| Line | Condition | Formula | Type |
|------|-----------|---------|------|
| I   | Top@transfer ≤ fci | `(1/Ac − e/St) / (fci − Mi/St)` | lower bound |
| II  | Bot@transfer ≥ −fti | `(1/Ac + e/Sb) / (Mi/Sb − fti)` | upper bound |
| III | Bot@service ≥ −fts | `η·(1/Ac + e/Sb) / (MT/Sb − fts)` | upper bound |
| IV  | Top@service ≤ fcs | `η·(1/Ac − e/St) / (fcs − MT/St)` | lower bound |

Feasible: `max(I, IV) ≤ 1/F ≤ min(II, III)`

## Eccentricity Bounds (for fixed F, varying x)

| Condition | Type | Formula |
|-----------|------|---------|
| Top compr @ transfer | e ≥ | `St/Ac + (Mi − fci·St)/F` |
| Top tens @ transfer | e ≤ | `St/Ac + (Mi + fti·St)/F` |
| Bot tens @ transfer | e ≥ | `(Mi − fti·Sb)/F − Sb/Ac` |
| Bot compr @ transfer | e ≤ | `(Mi + fci·Sb)/F − Sb/Ac` |
| Top compr @ service | e ≥ | `St/Ac + (MT − fcs·St)/(η·F)` |
| Top tens @ service | e ≤ | `St/Ac + (MT + fts·St)/(η·F)` |
| **Bot tens @ service** | **e ≥** | `(MT − fts·Sb)/(η·F) − Sb/Ac` ← usually governs |
| Bot compr @ service | e ≤ | `(MT + fcs·Sb)/(η·F) − Sb/Ac` |

## Design Stages (defineDesignStages.m)

| Stage | η | Loads | Compression limit | Tension limit |
|-------|---|-------|-------------------|---------------|
| Transfer | 1.0 | SW only | 0.60·f'ci | −3√f'ci |
| Service_Sustained | 1−losses | SW + SDL | 0.45·f'c | −12√f'c |
| Service_Total | 1−losses | SW + SDL + LL | 0.60·f'c | −7.5√f'c |

## Common Mistakes to Catch
1. Applying SDL moment at transfer (wrong — SDL applied after transfer for precast)
2. Using η < 1.0 at transfer (wrong — full F at transfer)
3. Mixing y-origin (some comments say top, code uses bottom)
4. Forgetting that sagging moment → compression at TOP (+Mi/St) and tension at BOTTOM (−Mi/Sb)
5. Denominator sign flip in Magnel lines when stress from moment exceeds allowable
