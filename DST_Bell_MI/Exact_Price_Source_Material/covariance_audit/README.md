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

# Deposit version 3 — the bandwidth of the census (October 21–23, 2026)

Zenodo version 3 of the census record (concept DOI 10.5281/zenodo.23537610, which
resolves to the latest version). Everything below was computed on the author's
workstation from the deposited random0 / random1 / strat0 / strat1 banks; every log
begins with the self-check that reproduces md = 0.138071187458. Note:
census_bandwidth_note_zenodo_v3.pdf (one page, all tables). Nothing here is in
arXiv v5; the sentences for v6 are held until the arXiv cool-off.

## Results I–L — band-limited censuses (census_bandlimit_v1–v4.py)

Question: for a census built only from harmonics exp(2i(k_A xi_A + k_B xi_B)) with
|k_A|, |k_B| <= K, how closely can it reproduce the singlet correlators, and what
does it pay above the floor when it can? Per (tau, s) cell the census over the 1024
offsets is a degree-K trigonometric polynomial; non-negativity at every atom,
normalization, correlator band (eps), CHSH and MD rows exactly as in
dither_sweep_rho.py. Two LPs per K: the price LP (min MD subject to the band) and
the miss LP (min band width any K-band census needs). v4 re-verifies both from the
reconstructed census and certifies infeasibility by value; v1–v3's "infeasible"
labels were read from a solver status and are superseded by the v4 certificates.
Launchers: run_parallel_v1/v2/v3.py, RUN_PARALLEL*.bat (batches 1–3).

- Minimum miss of the best K-band census, random0: K = 3 0.1052, K = 4 0.0174,
  K = 5 0.0120; K = 5 on random1 / strat0 / strat1: 0.0120 / 0.0131 / 0.0110.
  At K = 2 the miss LP itself is infeasible (no non-negative census satisfies the
  CHSH row). So reproducing the singlet to 1% requires K >= 6 on this terrain.
- Price (MD/floor) of the smoothest attaining census at eps = 0.01, random0:
  K = 6 2.274, 8 1.861, 10 1.598, 12 1.356, 14 1.121; K = 6 on random1 / strat0 /
  strat1: 2.094 / 2.019 / 2.081. The excess falls roughly linearly in K and is
  sampling-limited beyond K ~ 14 (the unrestricted census is representable at
  K ~ 16 by construction).
- (eps, K) surface, random0: K = 6: 2.274 / 1.913 / 1.590 / 1.451; K = 8: 1.861 /
  1.645 / 1.389 / 1.275; K = 10: 1.598 / 1.428 / 1.230 / 1.132 at
  eps = 0.01 / 0.02 / 0.05 / 0.10. Monotone in both variables.
- The skeleton (leading pixels (2,2), (2,0), (2,4), (1,1), (0,2), (2,-2), (1,-1),
  (3,-1) and their phases) is the same at every feasible K.
- Logs: census_bandlimit_v1_random0_10_21_26.txt, census_bandlimit_v2_random0_10_22_26.txt,
  parallel_console.txt, parallel2_console.txt, parallel3_console.txt, bl_*.txt.
  Coefficient files: census_bandlimit_v{1,2,3,4}_*rho_*.json (the v4 min-miss files
  carry the certificate values).
- Name collisions, stated for the record: bl_r0_e10_k6.txt and bl_r0_e10_k8.txt were
  written by batch 2 (v3, status misread as infeasible) and overwritten by batch 3
  (v4, feasible); bl_r1_k5.txt, bl_s0_k5.txt, bl_s1_k5.txt were likewise overwritten
  by batch 3. The batch-1/2 runs survive as the SUMMARY excerpts in
  parallel_console.txt / parallel2_console.txt and are superseded by the v4
  certificates. Files are deposited as they now stand plus the three console files.

## Result G — rank-k mixtures in sum/difference coordinates (product_census_rank_v2.py)

Mixtures of k separable censuses in (S, D) = (xi_A + xi_B, xi_A - xi_B), k = 1, 2, 3,
8 starts each, random0 banks: none feasible at eps = 0.01; best misses 2.65% ->
2.46% -> 2.37%. No short mixture reaches the floor (8 starts; not a certificate).
Log: product_census_rank_sumdiff_10_21_26.txt; result file:
product_census_rank_v2_sumdiff_rho_nsh32_random0_eps0.01_nang161.json
(product_census_rank_v1.py is the first version, kept for the record).

## Result H — pooled (S, D) spectrum (census_sumdiff_v1.py)

The four censuses pooled: eleven cross pixels above the permutation null; S preferred
at 39.5 deg mod 90 deg; D tracks S inside the lobes. Log:
census_sumdiff_v1_four_10_21_26.txt.

## Result M — the solve frame (dither_sweep_swap_v1.py)

Re-solving random0's offsets with Bob's setting pinned at 0 and Alice's swept (the
paper's solves pin Alice) attains the same floor (md = 0.138071187458) and the same
census: (2,2) at phi = 157.4 deg against ~157 deg in the pinned frame, preferred sum
39.35 deg against 39.5 deg, the same (S, D) cross pixels. The station asymmetry of
the census ((2,0) = 0.26 vs (0,2) = 0.13) stays on the same side, so it is a property
of the terrain, not of the solve frame; whether the terrain's own asymmetry is
physical or inherited from the breeding run (fitness evaluated with Alice pinned)
remains open. Files: banks_swap_nsh32_random0_nang161_nw2048.npz,
rho_swap_nsh32_random0_eps0.01_nang161.json, dither_swap_random0_10_22_26.txt
(console), RUN_SWAP_windows.txt (launch), census_swap_random0_10_23_26.txt and
sumdiff_swap_random0_10_23_26.txt (analyses), census_structure_v2_out_swap.json and
census_sumdiff_v1_out_swap.json (the analysis outputs, renamed from the 21 Oct
file names they overwrote; the 21 Oct content is in the 21 Oct transcripts).

## Scope

One terrain (the paper's greedy terrain), four censuses plus one roles-swapped
solve, 1024-offset banks. Feasibility and infeasibility of band-limited censuses are
certified by LP optima re-verified from the census; the rank-k statement is "no
feasible mixture from 8 starts". Nothing here derives rho; the certificates bound
the smoothness of any formation process from below.

### Version 4 (27 Oct 2026) — sampling correction, lattice census, dial dependence

**Correction.** The band-limited prices of version 3 (Results K–L) were computed on 1024 random offsets and include a quadrature cost: a band-limited census on N random atoms evaluates its correlators as an N-point random quadrature, and the LP pays MD to cancel the ~N^−1/2 noise. On a 32×32 lattice of offsets (`pipeline_4096_v2.py` with the grid option, job 400; `RUN_GRID32.bat`) the effect is absent. Lattice prices at ε = 0.01: 1.700, 1.429, 1.258, 1.112, 1.000 at K = 6, 8, 10, 12, 14 (`bl_n32_k*.txt`). 4096 random offsets (job 399; `RUN_4096_v2.bat`, `RUN_4096_MISS.bat`): 1.828, 1.428, 1.301, 1.219, 1.139, 1.071, 1.012 at K = 6–18 (`bl_n64_k*.txt`). The K ≤ 5 certificates hold on every sampling: best K = 5 miss 1.20 / 1.03 / 1.055% (1024 random / 4096 random / lattice), K = 4: 1.74 / 1.43 / 1.17%, K = 3: 10.5 / 12.6 / 8.5% (`bl_n64_k*_miss.txt`, `bl_n32_k2–5.txt`). The version-3 values stand for the 1024-random census; the lattice values are the model's smoothness price. The (ε, K) surface of version 3 is a 1024-random surface; its shape stands, its values do not transfer. Rows K = 4, 5 of that surface (`bl_r0_e*_k4.txt`, `bl_r0_e*_k5.txt`) are deposited for completeness.

**Result N — bandwidth, three samplings.** Floor attained at 1.000000 on 1024 random (×4), 4096 random (762 active of 61440), lattice (699 active of 15360). All three K-curves monotone. The floor census needs harmonics to K ≈ 14–18: angular features of 5–6°.

**Result O — the lattice census is the same census.** `census_structure_v2.py` and `census_sumdiff_v1.py` on `rho_nsh32_grid0`: (2,2) |F| = 0.314 at φ = 159.4°, preferred sum 39.84° mod 90°, (2,0) 0.271 vs (0,2) 0.162 (same-side station asymmetry), diagonal z = +4.2, anti-diagonal z = −0.6, same (S, D) cross pixels and both push families (`census_structure_v2_nsh32_grid0.txt`, `census_sumdiff_v1_nsh32_grid0.txt`, `census_structure_v2_out_nsh32.json`, `census_sumdiff_v1_out_nsh32.json`).

**Result P — the census depends on the dial** (`census_cells_v1.py`). 15 × 1024 matrix R[cell, offset] of the smooth K = 6 census: rank-1 shape share 0.38 (random0) / 0.51 (lattice); cell cosine similarity 0.33 / 0.50, and on the lattice 0.70 within a scale vs 0.43 within a τ. Best rank-1 approximation misses the singlet by 5.0% / 5.7%; rank 3 by 3.5% / 4.0% (`census_cells_random0_10_23_26.txt`, `census_cells_grid0_10_27_26.txt`, `census_cells_v1_out_grid0.json`). The lock axis rising with scale is suggestive, not established.

**Result Q — price of a dial-blind birth** (`census_cellsep_v1.py`, free shapes; alternating LP, 6 starts; upper bounds). One or two offset densities shared by all cells and reweighted per cell: lattice rank 1 miss 0.38%, price 1.473; rank 2 miss 0.08%, price 1.196; random0 1.78–1.84 / 1.55–1.67 with the quadrature cost (`census_cellsep_random0_10_25_26.txt`, `cs_free_s2.txt`, `cs_g_free_s1.txt`, `census_cellsep_v1_rho_*.json`). The sum-locked / difference-locked split of the random0 rank-2 solution is one local optimum, not replicated.

**Result R — smooth and dial-blind** (`census_cellsep_v2.py`, BAND = K). Rank 1 not found under 1.2% miss at K = 6 or 8 on either sampling; rank 2 at K = 6: 0.992% miss, 2.372×; rank 1 at K = 10: 0.97%, 2.269× (`cs_k*_r*.txt`, `cs_g_k*_r*.txt`, `census_cellsep_v2_*.json`). Best found, not certificates.

**Result S — parametric birth model** (`birth_model_v1.py` cold, `birth_model_v2.py` warm-started from the K = 6 census). ~30-parameter von Mises push model, cell-separable by construction: 2.4–2.6% miss from cold starts and 2.42–2.50% warm; fits the K = 6 census with 0.70 relative error; finds the axis (38.1°), station slots (38.4° mod 90°) and (2,2) phase (156°) but over-sharpens and lacks the high-D harmonics (`bm_*.txt`, `bm2_*.txt`, `birth_model_v*_*.json`). A floor-attaining birth model must carry harmonics to K ≈ 14–18 and vary with the terrain depth.

**Batches.** `run_parallel_v4.py`–`v8.py` with `RUN_BATCH4/5/6/8.bat`; consoles `parallel4–8_console.txt`. Batch 7 (`parallel7_console.txt`) is the 4096 K-curve and certificate batch; the 4096 pipeline launcher is `RUN_4096_v2.bat` / `RUN_4096_MISS.bat`.

**Housekeeping.**
- `pipeline_4096_v1.py` had no `__main__` guard and spawn-looped for 46 h under Windows multiprocessing (job 398, no output). Withdrawn; not deposited. `pipeline_4096_v2.py` is the tested replacement (tested under spawn; banks bit-identical).
- `census_cells_v1.py` writes `census_cells_v1_out.json` beside its first input; the 27 Oct lattice run overwrote the 23 Oct random0 file. The deposited `census_cells_v1_out_grid0.json` is the lattice run; the random0 numbers are in `census_cells_random0_10_23_26.txt`.
- Results Q, R, S are "best found"; only the K ≤ 5 misses and the LP prices are certificates or optima.
- K = 20 on 4096 (job 399) and the 4096 structure outputs (`census_structure_v2_nsh64_random0.txt`, `census_sumdiff_v1_nsh64_random0.txt`, `*_out_nsh64.json`) are added when the job completes — version 4 if before publishing, otherwise version 5.

Deposit: version 4 of 10.5281/zenodo.23537610. Note: `census_bandwidth_note_zenodo_v4.pdf`.

### Version 5 (29 Oct 2026) — the fine structure of the census: lattice surface, parametric wall, pleat, grain, crest, and a negative replication

**Result T — the lattice (ε, K) surface** (`run_parallel_v9.py`, `bl_g_e*_k*.txt`, `census_bandlimit_v4_g_*.json`). MD/floor at ε = 0.01 / 0.02 / 0.05 / 0.10: K = 5: infeasible (miss 1.055%) / 1.790 / 1.369 / 1.237; K = 6: 1.700 / 1.482 / 1.230 / 1.127; K = 8: 1.429 / 1.268 / 1.102 / 1.009; K = 10: 1.258 / 1.135 / 1.009 / ≤ 1.009 (the ε = 0.10 cell was stopped bracketed at 1.000–1.009; `bl_g_e10_k10.txt`, `bl_g_e10_k10_ds.txt` as they stand). Monotone in both variables. Replaces the 1024-random surface of version 3.

**Result U — a few-parameter law does not reach the floor** (`birth_model_v3.py`, `v4.py`; batches 10 and 12; `bm3_*.txt`, `bm4_*.txt`). Von Mises push model with a preferred axis, push harmonics 1–4, every parameter depending on terrain depth and lock duration, free cell weights and a station factor (57–77 parameters), warm-started from the K = 6 census: miss 2.05–2.2% (2.2–2.6% with lobes allowed down to 1.5°, v4). The fit to the K = 6 census improves with depth dependence (relative error 0.80 → 0.46); the miss does not. Best-found.

**Result V — not a ripple, not a crease** (`census_pleat_v1.py`, batch 11, `bl_pl_N*_lines.txt`, `bl_pl_N*_blocks.txt`). K = 6 census plus the four harmonics of an N-fold periodic ripple, N = 7–20, each cell free in amplitude and phase: 1.665 at N = 8 and 1.694–1.699 otherwise, against 1.700. Every sideband of ripple × smooth: contiguous bands to 14 or 18 reach 1.000; a band covering 10–22 with a hole at 7–9 pays 1.134. The fine structure is broadband and contiguous from ≈ 7 to ≈ 14.

**Result W — which edge belongs to the mesh** (`pipeline_terrain_v1.py`, `RUN_TERRAIN.bat`; `terrain_ord3/`, `terrain_ord2/`). Terrain cut to station indices ≤ 3 (15 of 23 terms) and ≤ 2 (6 terms); both attain the floor. First feasible K = 6 / 5 / 4 for heaviest index 5 / 3 / 2 (full / cut 3 / cut 2); the top of the band moves only from ≈ 14 to ≈ 12.5 (K = 12: 1.112 / 1.136 / 1.034). The coarsest harmonic the census needs is set by the finest fence it meets; the finest structure it carries is not. The pre-registered prediction that the top would track the cut failed and is recorded as such.

**Result X — localized features** (`census_bumps_v1.py`, `v2.py`; batches 13–15; `bumps_cr_*.txt`, `bumps_ctl_*.txt`, `bumps_r0ctl_*.txt`, `bumps_n64ctl_*.txt`; `beads_inspect_v1.py`, `beads_control_v1–v3.py`). K = 6 plus a dictionary of von Mises bumps (half-width 12°, 9°, 6°, 3°; 16×16 or 32×32 centres; weights ≥ 0): 12° bumps unused; from 6° down the LP moves 67–100% of the mass into the bumps and the price falls with narrowness to the atom scale (1.700 → 1.502 → 1.300 → 1.048). Same on random0 (2.274 → 1.855) and 4096 random (1.828 → 1.384). The active bumps have no mean station phase in any τ (no traveling feature); their count does not follow the mesh (10.6 / 12.9 / 11.5 per cell on full / cut 3 / cut 2); they sit on the sum-lock ridge (S mod 90° at 22–45°, top three 2.8° bins carrying 0.36 / 0.74 / 0.63 of the mass on lattice / random0 / 4096); off-lattice their differences cluster at D ≈ 14° and ≈ 104°. Every bump has a twin at (ξ_A + 90°, ξ_B + 90°), an exact symmetry of the model. Vertices stable under objective perturbation (same-cell overlap 0.86–0.94 on all three samplings).

**Retired — Results Y and Z of note v19 (partition of difference values across cells; difference phase sorts the depth).** Found on the lattice with the vertex-perturbation control; does not replicate off-lattice:

| | lattice | random0 | 4096 random |
|---|---|---|---|
| between-cell D-overlap, z vs shuffle null | −5.2 | +0.1 | +1.3 |
| cells differing in scale only | −4.1 | −1.0 | −0.2 |
| cells differing in τ only | −2.0 | +0.3 | +0.4 |
| sorting map (D-bin → scale) | structured | unstructured | unstructured |

The v19 "sorting strength" statistic (0.58) is reproduced by unsorted censuses (0.54, 0.58) and had no null. Lattice-specific (bump centres coincide with atoms and the lattice symmetries tie the cells' columns); recorded with its numbers, not a result. Transcripts `beads_control_v3_random0_10_29_26.txt`, `beads_control_v3_nsh64_10_29_26.txt`, `beads_inspect_offlattice_10_29_26.txt`.

**Method note.** A lattice result on the census's fine structure must be replicated on random offsets before it is called a result, as a random-offset result on the smooth structure had to be replicated on the lattice (version 4). The two samplings fail in opposite ways: quadrature noise on random atoms, symmetry ties on the lattice. A statistic needs a null before it is read.

**Also.** `cellsep_inspect_v1.py` (`cellsep_inspect_10_27_26.txt`) inspected the rank-2 dial-blind decompositions of version 4 and retired the sum-locked / difference-locked reading.

**Still pending.** K = 20 on the 4096 census (job 399) and its structure transcripts → version 6.

Deposit: version 5 of 10.5281/zenodo.23537610. Note: `census_fine_structure_note_zenodo_v5.pdf`.
