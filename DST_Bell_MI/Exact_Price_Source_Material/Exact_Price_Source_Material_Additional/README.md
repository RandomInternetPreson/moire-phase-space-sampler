# Exact Price Source Material

Companion code and data for:

**The exact price of local realism in CHSH experiments: a measurement-dependence–detection
trade-off surface, a moiré phase-locking mechanism that saturates it, and an unmeasured
fringe in the fourfold coincidence sum** — Aaron Alai (2026).

Every numerical claim in the paper is reproducible from the scripts below. Requirements:
Python 3 with numpy and scipy (HiGHS via `scipy.optimize.linprog`). Exact verification
steps use only integer arithmetic and Python's `fractions.Fraction` (plus symbolic
Q(sqrt 2) reconstruction); floating point is used solely to propose candidates. The
Section 7 pipeline is deterministic end to end (seedless T = 0 dynamics); independent
re-executions reproduce every table entry bit-for-bit. Total wall time for the full
Section 7 pipeline is under two hours on a 56-core workstation; the remaining scripts
run in seconds to minutes on commodity hardware.

## Map: paper section → script → claim

### Section 3 — The price surface
- `exact_frontier_lp.py` — the S = 2·sqrt(2) section over the complete 81-strategy
  space (Fig. 1): exact LP values vs. the closed form and the falsified additive
  hypothesis.
- `price_surface_lp.py` — the full surface M(S, eta_eff) on the 32-point grid
  (S in [2.1, 3.4], eta_eff in [0.85, 1]); verifies Eq. (1) to machine precision.
- `md_floor_lp.py` — Hall's floor at eta_eff = 1 for the deterministic sign model
  (the M = (S − 2)/6 edge).

### Section 4 — Exact certification
- `certify_frontier.py` — Theorem 1. Float LP proposes the optimal support; every
  nonzero primal weight and dual multiplier is reconstructed in Q(sqrt 2) and verified
  symbolically. Certifies M(2·sqrt(2), 9/10) = (27·sqrt(2) − 33)/100 with matching
  primal and dual objectives.
- `certify_ladder.py` — Corollary 1. The same pipeline in plain `Fraction` arithmetic
  at the six rational (S, eta_eff) points: 1/60, 1/12, 209/4800, 369/2000, 2489/24000,
  1003/12000.

The certificates are reconstructed and verified in exact arithmetic at runtime: running
these two scripts *is* the certificate check.

### Section 5 — The exact detection sector
- `solve_h_spectral.py` — solves the master equation Cf = cos·Cg in arc form
  (Chebyshev series for h, Gauss–Jacobi(1/2, 1/2) quadrature); residual 1.5e−12,
  correlation match to −cos 2Δ at 4e−11, h(0)/h(1) = 1.17450854538894.
- `detection_profile_solve.py` — the direct solve for the smooth profile
  D(m) = sqrt(m)·h(m) with independent Monte Carlo verification.
- `sqrt_law_test.py` — numerical verification of the sqrt(m) edge law (Lemma 1).

### Section 6 — Squares, spheres, impossibility
- `tsirelson_pinning.py` — the adversarial check: 300 random positive harmonic
  mixtures with CHSH angles optimized against the bound; max never exceeds
  2·sqrt(2). Also demonstrates the positive-definite ⇒ Gram-function route.
- `pearle_band_solve.py` — the hard edge-band model reaching S = 4 with 48%
  coincidence-rate modulation (super-quantum postselection is purchased with
  angle-dependent rates).

### Section 7 — Moiré phase locking (the mechanism)
- `breed_t0.py` — greedy selection of the interaction weights (23 active beat terms)
  under the faithfulness objective. Emits `greedy_results.json`.
- `flight_time.py` — the uniform-register price table. `NANG=41` and the tolerance
  sweep `EPS` (env vars) reproduce the table
  (eps = 0.10 / 0.05 / 0.02 / 0.01 → 1.000 / 1.015 / 1.235 / 1.375). With the
  constraint grid set to [2°, 178°] it is the second, independently coded solver
  cross-validating 1.375 at eps = 0.01.
- `extract_pattern.py` — re-solves at eps = 0.01 and emits the residual-pattern
  fingerprint of Fig. 3 (`pattern.json`).
- `orchestra.py` — the detection-composition scan D_alpha = (1 − alpha) + alpha·D*:
  register price at alpha = 0 (independently recovering the eps = 0.02 price,
  M = 0.17048), M = 0 at eta_eff = 0.797 at alpha = 1.
- `flight_time_dt.py`, `extract_pattern_dt.py` — timestep-robustness variants
  (`DT` env var): DT = 0.004 → 1.375, DT = 0.002 → 1.345, first-order
  extrapolation ≈ 1.32, as quoted in the paper.
- `slice_check.py` — the adversarial covariance test: evaluates the eps = 0.01
  solution (from `pattern.json`) on the unconstrained a = 45° slice; residuals up
  to 0.065, as stated in Sec. 7.
- `greedy_results.json` — output of `breed_t0.py`, consumed downstream.
- `pattern.json` — reference copy of the Fig. 3 fingerprint emitted by
  `extract_pattern.py`.

### Section 8 — The guilt trade-off
- `guilt_tradeoff.py` — minimal measurement dependence vs. fourfold-sum anisotropy
  along the interpolating detection family (Fig. 4); the delta_f < 2% / 1% / 0.5%
  → 45% / 68% / 80%-of-floor bounds.

### Section 10 — The factorization-vetoed protocol
- `protocol_sim.py` — the settings-torus mock experiment (Fig. 5): factorizable
  apparatus ripple confined to the log-Fourier axes, physical signal on the
  anti-diagonal pixels, null recovery (0.001 ± 0.015)% with Poisson noise at
  1e6 counts per point.

## data/ — reference outputs

Saved outputs from the author's runs, so results can be diffed without re-running:
- `exact_frontier.npz`, `price_surface.npz` — Sec. 3 LP values.
- `detection_profile.npz`, `h_spectral.npz`, `sqrt_law_h.npz` — Sec. 5 solutions.
- `pearle_band_solution.npz` — Sec. 6 edge-band model.

## Note on provenance

The detection structure investigated here was originally suggested by a beat-pattern
(moiré) heuristic; none of the results depend on it. See Sec. 11 of the paper.
