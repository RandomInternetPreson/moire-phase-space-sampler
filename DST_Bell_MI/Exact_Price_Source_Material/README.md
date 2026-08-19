# Exact Price Source Material

Companion code for:

**The exact price of local realism in CHSH experiments: a measurement-dependence–detection
trade-off surface, a moiré phase-locking mechanism that saturates it, and an unmeasured
fringe in the fourfold coincidence sum** — Aaron Alai (2026).

This directory contains the Section 7 (moiré phase locking) pipeline named in the
paper's Code and data availability section. The pipeline is deterministic end to end
(seedless T = 0 dynamics); independent re-executions reproduce every table entry
bit-for-bit. Requirements: Python 3 with numpy and scipy. Total wall time is under
two hours on a 56-core workstation.

## The Section 7 pipeline

- `breed_t0.py` — selects the interaction weights (23 active beat terms) under the
  faithfulness objective.
- `flight_time.py` — produces the uniform-register price table. Run with `NANG=41`
  and the tolerance sweep `EPS` (environment variables) to reproduce
  eps = 0.10 / 0.05 / 0.02 / 0.01 → M/M_min = 1.000 / 1.015 / 1.235 / 1.375.
  With its constraint grid set to [2°, 178°] it serves as the second,
  independently coded solver cross-validating 1.375 at eps = 0.01.
- `extract_pattern.py` — re-solves at eps = 0.01 and emits the residual-pattern
  fingerprint of Fig. 3 (`pattern.json`).
- `orchestra.py` — the detection-composition scan D_alpha = (1 − alpha) + alpha·D*;
  independently recovers the eps = 0.02 price at alpha = 0 (M = 0.17048) and
  returns M = 0 at eta_eff = 0.797 at alpha = 1.
- `pattern.json` — reference copy of the Fig. 3 fingerprint emitted by
  `extract_pattern.py`, so results can be diffed without re-running.

## Everything else

The full section-by-section source material — the Section 3 price-surface linear
programs, the Section 4 exact-certification pipeline (Theorem 1 in Q(sqrt 2) and the
Corollary 1 rational ladder), the Section 5 spectral solvers, the Section 6
adversarial checks, the Section 8 guilt trade-off, the Section 10 protocol
simulation, and reference output data — is in:

**[`Exact_Price_Source_Material_Additional/`](./Exact_Price_Source_Material_Additional/)**

See the README there for a script-by-script map from each paper section to the
claims it verifies.

`Visualizer_WIP/` is an exploratory visualization not referenced by the paper.
