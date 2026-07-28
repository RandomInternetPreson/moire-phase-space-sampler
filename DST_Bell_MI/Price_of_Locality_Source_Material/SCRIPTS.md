# Price of Locality — script inventory

## Certification (self-certifying assertions; these establish the paper's exact values)
- `exact_certificates.py` — integer-arithmetic squeeze certifying every Mermin floor n=3–13.
- `delta_R_check.py` — Mermin product identity / frustration-gap identity, exhaustive n=3–9.
- `ultrametric_check.py` — landscape geometry: no bit-hierarchy ultrametric grades the
  violation shells (0 of 120 bit orderings at n=4, exhaustive, definition stated in the
  docstring); Hamming-proximity failure (maximal violation at distance one).
- `relaxation_check.py` — exact zero-temperature non-ergodicity of single-site dynamics
  on constraint sectors (n=3,4) plus seeded contextual Metropolis with collective moves
  relaxing to uniform contextual ground shells and regenerating F = 1/3.
- `inheritance_check.py` — coherence-support inheritance (stabilizer support; emission model).
- `eta_check.py`, `finite_beta.py`, `vise_toy.py` — finite-β analysis of the companion program.

## LP lineage (discovery; double-precision, superseded for certification by the squeeze)
- `quantumtestv1.py` — full strategy-level faithful LP, all marginal constraints explicit
  (n=2,3,4). Independent numerical confirmation of the reduction theorem at small n.
- `quantumtestv2.py` — sparse strategy-level LP (n through ~5–6).
- `quantumtestv3.py` — fiber-reduced class-level LP, 2^(n+1) classes, marginals free (n through 7–8).
- `quantumtestv4.py` — class-level LP, production version.
- `quantumtestv5.py` — S_n symmetry-reduced exact LP; n=9 in ~200 density variables; recomputes
  the staircase n=3–13 in under a minute.
