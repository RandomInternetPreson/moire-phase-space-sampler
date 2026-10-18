# Settings-torus brief and protocol code

Experimental brief and protocol scripts for the fourfold-coincidence-sum test proposed in
arXiv:2608.18886 (Secs. 9–10).

## Current files (v4 of the paper, October 2026)

- `settings_torus_brief_v2_2.pdf` / `.tex` — the two-page brief: the observable and what
  existing data say about it, Design A (the 2×2 cross-ratio at 0° and 45°, about 10^8
  fourfold coincidences for 5σ at 0.1%) and Design B (the settings torus, about
  2.3×10^8), the two-channel efficiency-imbalance term and its singles-based
  calibration, controls, requirements, and the three outcomes.
- `protocol_sensitivity_v2_2.py` — the protocol script named in the paper's Code and
  data availability section: outcome-dependent channel efficiencies, Poisson singles,
  the uncalibrated and singles-calibrated (2,−2) log-Fourier estimators, the cross-ratio,
  the imbalance nulls (constant, angle-dependent, and the equal-phase case), and the
  coincidence-marginal check. Reference figures: with 3% factorizable ripple plus 20%
  constant and 3% angle-dependent imbalance per side at 10^6 coincidences per point, the
  uncalibrated null is (−0.041 ± 0.017)% and the calibrated null (0.000 ± 0.006)%.
- `hance_checks.py` — numpy port of the independent reviewer's Mathematica checks
  (sections 3–5 of their script), reproducing each printed figure: the imbalance pixel
  values, the multi-pair identity, the Larsson / symmetrised Gisin–Gisin efficiency
  4/(π+2) = 0.778, and the 2×2 count estimate.

The exact-profile harmonics used by the brief and the script
(a_1 = 0.0577, a_2 = 0.0098, a_3 = 0.0031, …) are produced by
`../Exact_Price_Source_Material_Additional/solve_h_spectral.py` and stored in
`../Exact_Price_Source_Material_Additional/data/fringe_harmonics.npy`.

## Record

- `settings_torus_brief.pdf` / `.tex`, `protocol_sensitivity.py` — the August 2026
  originals (paper v1–v3), kept unchanged. The v2.2 brief corrects them: the fringe
  takes equal values (not zero) at the CHSH relative angles; Aspect et al. did check the
  sum's constancy (±1.5%) and Brida et al. (2007) do bound its single-channel analogue
  (Δ_exp = 0.0074), so existing data exclude the full-amplitude profile; flatness does
  not force measurement dependence (Larsson 1999); the non-factorizable imbalance term
  and its calibration are new; the time estimate is 1.3 h at 5×10^4 coincidences/s
  (the "7 hours" of the original was a carry-over); and Design A is new.
