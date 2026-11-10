# Discrete-time three-species commensalism–amensalism model: verification code

This repository contains the Wolfram Mathematica and Python code behind every computation,
numerical example and figure in the manuscript

> **F. Ahmad, M. K. Jamil** — *Stability, bifurcation and chaos control of a discrete-time
> three-species commensalism–amensalism model with Beddington–DeAngelis functional response*
> (manuscript in preparation; working title).

The code is arranged in the same order as the paper, so each section can be checked on its own.

## The model

The continuous model of Shofiyya et al. (2025) is discretized with the forward Euler scheme, step size `h`:

```
x_{k+1} = x_k + h [ x_k (a1 - b1 x_k) + beta z_k ]
y_{k+1} = y_k + h y_k [ a2 - b2 y_k + c x_k / (m x_k + n y_k + 1) - alpha z_k ]
z_{k+1} = z_k + h z_k [ a3 / x_k - b3 z_k ]
```

Parameter order used everywhere in the code: `(a1, b1, a2, b2, a3, b3, beta, alpha, c, m, n)`.

| Example (Section 8) | Parameters | Bifurcation of the coexistence point E4 |
|---|---|---|
| 8.1 | `(1.2, 0.6, 0.5, 0.4, 0.9, 0.7, 0.3, 0.6, 2.0, 0.5, 0.4)` | flip in the y-direction at h = 1.0658 |
| 8.2 | `(3.0, 1.5, 1.0, 0.5, 1.5, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2)` | flip in the (x,z)-block at h = 0.5644 |
| 8.3 | `(2.0, 3.0, 0.5, 0.2, 1.5, 3.0, 0.5, 0.5, 3.0, 1.0, 2.0)` | Neimark–Sacker at h = 0.7806 |

## Repository layout

```
mathematica/
  00_model.wl              shared definitions (model, fixed points, Jacobian, forms B and C)
  01_fixed_points.wl       Sections 3-4: symbolic checks of Theorem 4.1, fixed points of the examples
  02_stability.wl          Section 5:   factorization (5.2), trace/determinant (5.6), thresholds h*, 2/rho4, h_E4
  03_bifurcation.wl        Section 6:   flip coefficients (6.7)-(6.8), Lyapunov coefficient (6.10),
                                        validation on a known map, all conditions of Theorems 6.1-6.3
  04_chaos_control.wl      Section 7:   Lemma 7.1 and Theorem 7.1 for the three examples
  run_all.wls              runs 01-04 in order
  output/                  reference output of each script
python/
  model.py                 shared definitions (vectorized map, Jacobian, bifurcation and MLE routines)
  report.py                prints results and writes output/*.txt with PASS/FAIL checks
  01_fixed_points_and_stability.py   Sections 4-5, compared with the values printed in the paper
  02_bifurcation_checks.py           Section 6 and Section 8: coefficients, period-two orbits,
                                     invariant curve, cascades, period-5 locking
  03_chaos_control_checks.py         Section 7 and Section 8: MLE, B*, control on both sides of B*
  04_make_figures.py                 regenerates all eight figures of Section 8
  run_all.py                         runs 01-04 in order
  output/                  reference output of each script
figures/                   the eight figures of the paper (PDF, with PNG previews)
```

## Where each result of the paper is checked

| Paper | Mathematica | Python |
|---|---|---|
| Proposition 3.1, Theorem 4.1, Remark 4 | `01_fixed_points.wl` | `01_fixed_points_and_stability.py` |
| Eqs. (5.1)–(5.9), Lemma 5.2, Theorems 5.1–5.4 | `02_stability.wl` | `01_fixed_points_and_stability.py` |
| Theorem 6.1 and coefficients (6.7) | `03_bifurcation.wl` (a), (c) | `02_bifurcation_checks.py` |
| Theorem 6.2 and coefficient (6.8) | `03_bifurcation.wl` (c) | `02_bifurcation_checks.py` |
| Theorem 6.3 and coefficient (6.10) | `03_bifurcation.wl` (b), (c) | `02_bifurcation_checks.py` |
| Remark 9 (subcritical case) | `03_bifurcation.wl` (c) | `02_bifurcation_checks.py` |
| Lemma 7.1, Theorem 7.1 | `04_chaos_control.wl` | `03_chaos_control_checks.py` |
| Examples 8.1–8.3, Table 1 | all four scripts | all four scripts |
| Figures 1–8 | – | `04_make_figures.py` |

The Mathematica and Python implementations of the normal-form coefficients were written independently
and give the same values. The Neimark–Sacker formula is also validated on the delayed logistic map,
whose bifurcation is known to be supercritical.

## How to run

**Python** (3.9 or later; tested with Python 3.13, NumPy 2, Matplotlib 3.11):

```bash
cd python
pip install -r requirements.txt
python run_all.py          # or run 01_... to 04_... one by one
```

Each verification script prints its results with `[PASS]`/`[FAIL]` tags and saves them in `python/output/`.
The whole run takes under a minute. Currently: 39 + 21 + 15 checks, all passing.

**Mathematica** (version 13 or later; tested with Wolfram Engine 15.0):

```bash
cd mathematica
wolframscript -file run_all.wls      # or: wolframscript -file 01_fixed_points.wl, ...
```

The scripts can also be opened and evaluated in a Mathematica notebook front end. Each script loads
`00_model.wl` from its own folder. Reference output is in `mathematica/output/`.

## Numerical settings

Bifurcation diagrams use 1400 parameter values, start from E4 + (0.01, 0.01, 0.01), discard 3000–4000
transient iterates and plot the next 250–300. Maximal Lyapunov exponents are computed by the
tangent-vector method, averaged over 5000–6000 iterates. Orbits that leave the positive octant are
discarded, and each diagram is cut where this first happens.

## Licence

The code is released under the MIT licence (see `LICENSE`). If you use it, please cite the paper above.
