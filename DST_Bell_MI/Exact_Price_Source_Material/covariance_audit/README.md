# covariance_audit — v3 supporting runs (August 26–27, 2026)

Version 2 of arXiv:2608.18886 attributed a 0.065 correlation residual to broken
rotational covariance of the frozen offsets. That attribution was wrong: the
mechanism is exactly covariant (a global rotation shifts every phase phi_mn+- by
2(m+-n)s, matching the beat argument's 2ks in both branches), verified here on
four slices (a = 0, 22.5, 45, 67.5 deg) which return identical correlation
tables to four decimals. The 0.065 is a between-pin register excursion, present
on the constrained slice itself. Version 3 corrects the text and adds:

- Density curve, 16x16 lattice, eps = 0.01: MD/floor = 1.375 / 1.511 / 1.635
  at 41 / 81 / 161 pinned angles (slice_solve2.py; pattern_rho_original.json,
  pattern_rho_dense81.json, pattern_rho_dense161.json; MD values recorded in
  run transcripts).
- Two-slice register (a = 0 and a = 45 pinned): identical solution and price to
  single-slice at 16x16 — direct consequence of exact covariance
  (pattern_rho_twoslice.json).
- 1024 random offsets, 161 pins, eps = 0.01: MD = 0.1380711874577208 =
  (sqrt(2)-1)/3 to twelve significant digits; support 656/15360 atoms, top-12
  carry 10% (dither_sweep_rho.py; rho_nsh32_random0_eps0.01_nang161.json).
- Between-pin probes of that floor solution: sampled residuals <= 0.026 in the
  historically worst window (140-148 deg, 0.2 deg steps) and <= 0.017 over the
  full range (slice_check_mp.py; slice_check_a0_d140-148.json,
  slice_check_a0_d10-170.json — the "pat" field records which solution each
  probe evaluated).

File notes. pattern_rho_original.json is the 16x16, 41-angle solution behind
the paper's 1.375 price (schema: rho/nsh/nw/eps; distinct from ../pattern.json,
which is extract_pattern.py's fingerprint file in a different schema —
extract_pattern_rho.py is the producer of the rho schema).
pattern_rho_dense81.json and pattern_rho_dense161.json are the 81- and
161-angle lattice solutions (MD = 0.208598 and 0.225798, stored in-file).

Covariance verification: slice_check_orig_a{0,22.5,45,67.5}_d10-170.json are
four probes of pattern_rho_original (the "pat" field records this), returning
identical tables, worst residual 0.0651 on every slice. The historical April
2026 probe of the same solution on a = 45 (slice_check.json +
slice_partial.json, worst 0.06514, produced by the superseded single-threaded
slice_check.py in the parent folder) is retained as provenance; its file
names carry no window suffix because they predate the windowed prober.
slice_check_a22.5.json (two_slice = 1) probes the two-slice solution and
returns the same table — the a = 45 constraint rows duplicate the a = 0 rows
under exact covariance and the LP returns the same vertex.
slice_check_a0.json probes the dense-161 lattice solution (worst 0.0569).

Between-pin probes of the 1024-offset floor solution (see each file's "pat"
field): slice_check_floor1024_a0_d140-148.json (worst 0.0258, the
historically worst window at 0.2 deg steps) and
slice_check_floor1024_a0_d10-170.json (worst 0.0173, full range) — the
residuals cited in the paper.

The "ratio" field inside rho_nsh32_random0_eps0.01_nang161.json was computed
with a floor constant rounded to six digits (0.138071) and reads
1.0000013577; the "md" field is the twelve-digit value, md =
0.1380711874577208 vs (sqrt(2)-1)/3 = 0.1380711874576983...
dither_sweep_rho.py now uses the exact floor, so subsequent dumps carry a
consistent ratio.

Dependencies shipped for standalone reproduction: greedy_results.json (the 23
beat terms, weights, and v) and flight_time.py (the module imported as ft:
banks, correlator signs, S_QM, min_md). dither_dense161.json and the
eps = 0.02 / 0.10 rho files complete the record when the corresponding LPs
finish (both values are forced to the floor: Hall's bound lower-bounds any
eps, and monotonicity upper-bounds them by the eps = 0.01 value).


# v5 supporting runs — census structure (October 15–21, 2026)

Version 5 of arXiv:2608.18886 adds to Sec. 7 what the floor-attaining census
looks like, that it is a function of the terrain while the price is not, and
that it cannot be built from the two stations separately. The scripts, solutions
and transcripts below are also deposited, with a short note, at
https://doi.org/10.5281/zenodo.23537611 (version 1, 21 Oct 2026; concept DOI
10.5281/zenodo.23537610). Every quoted number is produced by the named script
from the named input; each script asserts the stored floor price of its input
before reporting.

## Corrections and additions (21 Oct 2026, after two cold reads; deposit version 2)

- Result A, statement (2): the sum term is cos(4(xi_A + xi_B) - phi), phase outside
  the factor of 4; the preferred value of xi_A + xi_B is phi/4 ~ 40 deg (mod 90 deg)
  in the pinned frame. Result A table: active atoms (weight > 0) random1 675,
  strat0 684, strat1 669 (656 in the paper's Sec. 7 is the count at weight > 1e-6).
  Result C control: MD - floor = 2e-14.
- Result D: the NB = 32 run used 14 starts per class (uniform, rho*-marginals,
  random0-11), not 20; "max slack" in the transcript is the excess over the eps band,
  so the correlators miss their targets by eps + slack = 3-13% at eps = 0.01; at
  eps = 0.05 the station-local classes pay 2.20 / 3.07 times the FLOOR (the binned
  joint control hit its time limit there). The NB = 16 run (product_census_v2.py,
  6 starts) is recorded in the author's notes only; its output was not retained.
- Per-pixel null for the (2,-2) term (pixel_null_v1.py, same measure, 400 draws,
  seed 20261015, independent re-pairing): |F(2,-2)| = 0.127 / 0.152 / 0.143 / 0.124
  on random0 / random1 / strat0 / strat1 against null 95th percentiles
  0.117 / 0.115 / 0.099 / 0.106 (p = 0.038, 0.005, 0.008, 0.023), phases
  +14 / +45 / +50 / +48 deg. A weaker relative-phase term cos(4(xi_A - xi_B) - psi)
  is present at about a third of the sum term. Transcript:
  pixel_null_v1_run_2026-10-21.txt. Two-machine (container and author's machine,
  identical to all printed digits).
- Result F, sum/difference coordinates (product_census_v4.py, COORDS=sumdiff,
  random0 banks, NB = 32, 14 starts per class, 17,655 s): at eps = 0.05 semi 1.965
  (14/14 feasible), product 2.749 (12/14); at eps = 0.01 joint 1.107, no semi or
  product census feasible from any start (targets missed by 1.7-3.0% and 2.2-3.9%).
  Closer to separable than in station coordinates, but not separable: the
  dominant (2,2)_AB pixel is the pure-S harmonic, the other four skeleton pixels
  are cross terms in (S, D). A census depending on xi_A - xi_B alone is the
  f(S) = const sub-case of the class that fails. Transcript and result file:
  product_census_v4_sumdiff_random0_10_21_26.txt, product_census_v4_sumdiff_rho_nsh32_random0_eps0.01_nang161.json.
- Deposit: version 2, https://doi.org/10.5281/zenodo.23540603 (concept DOI
  10.5281/zenodo.23537610, cited by the paper, resolves to the latest version).
  Note: census_structure_note_zenodo_v2.pdf replaces v1 in this folder.

## Result A — the floor-attaining census has a universal low-order skeleton

Four independent 1024-offset solutions at the floor (dither_sweep_rho.py,
NANG = 161, NSH = 32, eps = 0.01, SAVE_RHO = 1): random seeds 0 and 1 and
stratified seeds 0 and 1 — rho_nsh32_{random0,random1,strat0,strat1}_eps0.01_nang161.json,
each with MD/floor = 1.000 to twelve digits, sharing no offsets.
census_structure_v2.py (v1 + phases of the top pixels) forms the weighted
empirical measure over (xi_A, xi_B), its circular Fourier plane |F(k_a,k_b)|
for |k| <= 6, and a 400-draw seeded permutation null (weights fixed, offsets
re-paired). On all four: top pixel (2,2) at |F| = 0.373 / 0.352 / 0.349 / 0.384
against a null 95th percentile of 0.17–0.18, with phase +156.5 / +161.1 /
+157.1 / +158.0 deg (spread 4.5 deg); second pixel (2,0) at 0.24–0.31; the
anti-diagonal (k,-k) share at the null (z = +0.9 / -0.7 / -0.7 / -0.9), the
diagonal (k,k) share at z = +5.8 / +6.2 / +5.8 / +7.1. So rho* is not a
function of xi_A - xi_B alone; its dominant structure is cos 4(xi_A + xi_B - phi)
with phi fixed to within 5 deg, plus a nonuniform xi_A marginal cos 4 xi_A.
Transcripts: census_structure_v2_{random0,random1,strat0,strat1}_run_2026-10-21.txt
(regenerated 21 Oct from the stored solutions; identical to the October 15–16
runs). The (2,0) > (0,2) asymmetry is a property of the solve (a pinned at 0,
b swept), not of the stations.

## Result B — the terrain is not a product of two screens

terrain_factor_v1.py on greedy_results.json (23 terms, v = 15.2339): a product
of two real cosine screens predicts w_mn+ = w_mn- for every (m,n); observed
ratios run 0.05 to 183; difference beats (k = m - n) carry 72% of the weight.
A statement about the greedy-found terrain only. Transcript:
terrain_factor_v1_run_2026-10-21.txt.

## Result C — terrain perturbation: the census follows the terrain, the price does not

flight_time_perturb_v1.py is flight_time.py with a PERTURB environment
variable (identity when unset; module bit-identical to the original in that
case); dither_sweep_perturb_v1.py is dither_sweep_rho.py driven by it, with
tagged outputs. Two re-solves on the random0 offsets (NANG = 161, NSH = 32,
eps = 0.01), matched control random0:
- Run A, (2,2,+) term phase shifted +45 deg
  (rho_perturb_2_2_1_shift_45_nsh32_random0_eps0.01_nang161.json):
  MD - floor = 5e-8; the (2,2) census pixel moves by -38 deg (the sign the
  convention sin(2k lambda - ph) dictates) and grows to 0.42.
- Run B, (2,2,+) term removed
  (rho_perturb_2_2_1_weight_0_nsh32_random0_eps0.01_nang161.json):
  MD - floor = 3e-14; a (2,2) pixel persists at 0.31 at a different phase,
  from harmonic mixing of the remaining terms.
The floor is attained on all three terrains; the census skeleton is
phase-locked to the terrain.

## Result D — no station-local census reaches the floor

product_census_v3.py tests, on the random0 banks
(banks_nsh32_random0_nang161_nw2048.npz — the August banks, the matched pair
of rho_nsh32_random0; self-check reproduces the stored price to twelve digits),
whether a census of the form q(tau,s) f(xi_A) g(xi_B) ("product"), or one with
xi_B independent of everything else ("semi"), reproduces the singlet
correlators at tolerance eps, against a binned joint census as control.
NB = 32 bins, 12 starts per class, eps = 0.05 and 0.01
(product_census_v3_random0_NB32_oct_17_16.txt,
product_census_v3_rho_nsh32_random0_eps0.01_nang161.json): joint 1.119x the
floor at eps = 0.01 (the excess is the binning); semi 2.203 and product 3.069
at eps = 0.05; no semi or product census feasible at eps = 0.01 (final
correlator violations 2–15%). The NB = 16 run (joint 1.59 at eps = 0.01;
station-local classes infeasible at eps <= 0.02, 2.2–3.0x the joint price where
feasible) was made with product_census_v2.py and is recorded in the paper and
the deposit note; v1 and v2 of the script are withdrawn (v1 had a normalization
error; v3 = v2 + a time limit on the joint control LP, a SKIP_JOINT switch, and
zero-mass columns dropped) and are not shipped. product_census_v4.py is v3 +
COORDS=sumdiff (the sum/difference-coordinate test); its results are not part
of v5.

## Result E — a terrain bred to favour station-local censuses does not change this

breed_product_v1.py breeds a terrain from greedy_results.json scored by the
product-class price (breeding scale NSH = 8, NANG = 9, eps = 0.10);
pipeline_v1.py / RUN_PIPELINE.bat run smoke -> breed -> validate unattended
(pipeline_status.txt and pipeline_FINAL.txt record every transition). Run of
18 Oct 2026 (14.1 h): breeding stopped after 5 rounds with no improving move,
joint 1.631x, product 1.927x, ratio 1.18 ("promising" under the pre-registered
rule, ratio <= 1.3); full-scale validation of the champion (NANG = 161, 1024
random offsets, NB = 32, 12 starts): eps = 0.05 joint 1.261 / semi 2.492 /
product 2.864; eps = 0.01 joint 1.591, semi and product infeasible. The bred
terrain, its banks, its solution and its census test live in
A_breed_greedy/ and VA_validate_greedy/ (the files there carry the same names
as the paper's files but are the bred-terrain versions — compare by hash:
the paper's greedy_results.json is sha256 e2d53682..., the bred champion
e3c735a4...). Stage B (minimal start) was not run; the exclusion reads
"from the greedy start".

## Standing conclusion (as stated in v5, Sec. 7)

On the paper's terrain, on four independent censuses, and on a terrain bred
against the conclusion, no product census reproduces the singlet at
eps <= 0.02, and where feasible it pays 2–3x the joint price: the frozen
offsets of the two stations must be correlated at the source. Nothing here
derives rho; the skeleton is what any derivation must produce.

## Files added

Scripts: census_structure_v1.py, census_structure_v2.py, terrain_factor_v1.py,
flight_time_perturb_v1.py, dither_sweep_perturb_v1.py, product_census_v3.py,
product_census_v4.py, breed_product_v1.py, pipeline_v1.py, RUN_PIPELINE.bat,
check_folders_v1.py (folder/hash audit). Solutions: the three new rho files
(random1, strat0, strat1), the two perturbation rho files, and their banks.
Results: product_census_v3 outputs, A_breed_greedy/, VA_validate_greedy/,
breedrun/logs/, pipeline_FINAL.txt, pipeline_status.txt. Transcripts as named
above. Note: census_structure_note_zenodo_v1.pdf (the deposit note).


