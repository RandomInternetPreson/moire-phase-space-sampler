# Price of Locality — script inventory

Companion code for arXiv:2608.00124, *Exact minimum measurement dependence for faithful
local deterministic models of multipartite GHZ–Mermin correlations* (A. Alai, 2026; v2,
October 2026). Requirements: Python 3 with numpy and scipy (HiGHS via scipy ≥ 1.9).
Every script asserts the numbers the paper quotes from it; an `AssertionError` is a
finding. The scripts import one another by filename — keep them in one directory.

## Exact certification of the Mermin values (Theorem 8; Table I)
- `exact_certificates.py` — integer-arithmetic certificate of every Mermin value n = 3–13
  (tabulates the satisfaction of every setting at every class).
- `exact_certificates_fwht_v2.py` — the same certificate for n ≤ 22 (violation counts from
  an integer Walsh–Hadamard transform, memory O(2^n)); asserts the Lemma 10–11 counts
  (|ground| = 2^n or 2^{n−1}, g = R(R+1), minimum pair overlap R(R+2)/2, equidistance,
  complementary pairs at even n). `--upto 22` takes about three minutes.

## Cluster-state cores (Theorem 9; Table II; Sec. X.B)
- `core_census_v1.py` — complete census of the minimal contradictions (cores) of the C4 chain
  and C5 ring at all sizes, by GF(2) kernel enumeration; asserts 8/8/8/16 and
  60/180/572/1980/5220/9640/8912.
- `cluster_cores_v4.py` — cores of C4, C5, C6 up to `--max-size` (default 5); for each core the
  Theorem 3 lower bound and the ground-state construction with every marginal checked in
  integer arithmetic (certificate i). Asserts the quoted counts when run without `--sample`
  (`--no-lp` recommended, ~25 min). `--lp-only START STRIDE` runs the collapsed-space
  faithful LP on the C6 cores where the construction is not faithful (20 parallel
  processes: `i 20`).
- `c6_sweep_tally_v1.py <dir>` — tallies the 20 `c6_lp_*.log` files; asserts 80 four-cores at
  1/3, 192 five-cores at 1/4, 112 five-cores at 1/3; writes `c6_above_cores.txt`.
- `core_dual_certificate_v2.py "<graph>" <k> <index|--core ...> [--expect P/Q]` — exact lower
  bound on the faithful minimum of one core from a rationalized LP dual, verified in
  rational arithmetic (certificate ii, lower half).
- `core_tight_bound_v3.py "<graph>" <sizes> <sample> [P/Q]` — the tightened construction
  (exact faithful model; certificate ii, upper half) on cores where the plain construction
  is not faithful; the fourth argument asserts its value on every core run.
- `c6_above_certify_v3.py c6_above_cores.txt START STRIDE` — both halves of certificate ii on
  each of the 112 C6 five-cores with LP value 1/3; a full run asserts all 112.
- `c6_exact_primal_v2.py <dir> START STRIDE` — exact rational faithful models at 1/(k−1) on
  the 160 C6 cores (28 four-cores, 132 five-cores) whose LP value equals the Theorem 9
  bound (certificate iii); a full run asserts 28 + 132.
- `c6_above_cores.txt` — the 112 five-cores, one per line.

## Scope computations (Sec. X.B–X.C)
- `pol_scope_checks_v3.py` — Mermin-5 minus one constraint (LP 2/5; exact B = 4/11 over all
  2^15 subsets; 180 cores = 60 of size 4 + 120 of size 6; construction feasible at 1/2 on
  12 ground classes) and the C4 four-core at strategy level (LP 1/3; construction faithful
  and optimal); asserts all of it. `--padded` additionally LPs all 55 padded scenarios
  (~5 min each).
- `c4_padding_construction_v2.py` — the C4 four-core padded with two further stabilizers, all
  55 choices: construction and marginals in integer arithmetic; asserts 43 faithful at 1/3,
  8 not faithful, 4 faithful at 1/2 or 1, s = 5/6 throughout.
- `c4_bad_pairs_lp_v2.py <start> <stride>` — the faithful LP on the 12 padded scenarios where
  the construction is not optimal; asserts 1/3 on each.

## Landscape (Sec. X.D) and LP lineage (discovery; superseded for certification)
- `ultrametric_check.py` — no bit hierarchy grades the violation shells (0 of 120 orderings at
  n = 4); Hamming-proximity failure.
- `relaxation_check.py`, `delta_R_check.py`, `inheritance_check.py`, `eta_check.py`,
  `finite_beta.py`, `vise_toy.py` — earlier landscape and companion-program checks, not
  cited by v2 of the paper.
- `quantumtestv1.py` … `quantumtestv5.py` — the strategy-level, class-level and
  S_n-symmetry-reduced LPs used for discovery (v5 recomputes n = 3–13 in under a minute).

## Run logs
`run_logs/aaron_machine/` and `run_logs/container/` hold the outputs of every computation
quoted in the paper on two independent machines.
