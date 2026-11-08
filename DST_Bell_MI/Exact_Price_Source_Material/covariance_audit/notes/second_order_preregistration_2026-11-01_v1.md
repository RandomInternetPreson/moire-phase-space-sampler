# Pre-registration — the mesa's second-order geometry and the disc's edge (Nov 1, 2026, v1; written before the runs)

Follows `bulge_note_2026-11-01_v3.md`. Aaron's hypotheses of 1 Nov, each with its slot in the model, the code that tests it, and a prediction to be marked kept / failed.

## A. Sweep 2 — the mesa's second-order geometry (`run_bulge2_v1.py` 3159b1b6 → `census_profile_v3.py` 5a2eceeb; `RUN_BULGE2.bat` 1582212b; on the job-403 histogram store; 12 readers; full + K = 6, 8, 10)

| idea (Aaron) | reader | port rule | prediction |
|---|---|---|---|
| the disc tilted against the mesa's top | `tilt:c` — cap height falls linearly from the rear edge to the front edge by c = 0.05, 0.1 | kept | **infeasible**, like the dome: a tilted cap is a cap that is not flat; miss ≈ 0.4c (dome scaling) |
| a front and a back: the rear wall mirrored at a greater angle | `awall:w1:w2` — front wall 3.6° (≈ tanh:8), rear 7.2° / 14.3°; `awall:1.8:7.2` | **broken** when w1 ≠ w2 | the price is set by the *wider* wall: `awall:3.6:14.3` ≈ tanh:6 (1.46) or worse; marginals become physical (⟨A⟩_sym ≠ 0, correlators move under the +90° symmetrization) |
| the twin born upside down | `mawall:w1:w2` — Bob's reader is Alice's mirrored, D_B(θ) = D_A(−θ) | broken | same prices as `awall` to within a few percent; no gain — the correlator sees the product of the two walls and the mirror only swaps which edge meets which |
| banded walls (slope changes partway down) | `bwall:wt:wb` — the +1 half of the front wall has width wt, the −1 half wb (1.8/7.2 and 7.2/1.8) | kept | price between tanh:8 and the price of the softer half taken alone; the two kink orientations differ by < 10% |
| is the mesa a square or a circle? | `shift:d` — the wall's position follows the photon's approach phase, step(θ − d sin 2φ), d = 5°, 10°, 20° | kept | the *marginal* price falls below 1.000 (the LP exploits the extra freedom) — **that number is not the price**: the reader reads ξ, so the joint-(λ, atom) MD is the measure, and it is far above the floor (≳ 5×). Recorded as a method result: a ξ-reading detector moves the MD accounting, it does not buy the singlet cheaply |
| the grain question, all rows | K6/full and the absolute K = 6 | — | no reader brings K = 6 below 1.70; every shaped reader raises it (as tanh:8 did, 2.56) |

## B. The disc's edge — a non-uniform birth of λ₀ (`pipeline_fuzz_v1.py` 8a36d412 → `census_fuzz_v1.py` 4cc93176; `RUN_FUZZ.bat`; new particle store; step reader)
Slot: the paper starts every pair with λ₀ uniform and independent of (ξ, τ, s). Aaron's thickness / knife-edge / blunt / multimodal rim is the birth fuzz p(λ₀); the richest version is a joint birth over (λ₀, ξ_A, ξ_B, τ, s). The LP carries the joint density as w(atom, λ₀) = Σ_l w_l(atom) f_l(λ₀) with L harmonics in λ₀, band limit K in ξ, positivity at 64L λ₀-points (dips between points printed; ≤ 0.5% of the mean density in the smoke test), mass, CHSH exact, band ε, MD on the λ(τ) marginal (the reader reads λ(τ) only, so the paper's measure stands). Scopes: `cell` (rim depends on the dial only), `atom` (on the approach phases too). This **enlarges the model class**: if it helps, the paper must say the birth carries λ₀ structure.
Predictions: (1) L = 0 reproduces 1.700 / 1.429 / 1.258 (control). (2) `cell` scope changes the K = 6 price by < 2%: a rim that knows only the dial cannot stand in for structure in ξ. (3) `atom` scope lowers the K = 6 price measurably (the LP is given a genuinely new axis to put structure on) but **not** to 1.000 — my guess 1.4–1.6 at L = 2 — because the floor's fine structure was shown to be in the *difference* of the approach phases and λ₀ is one angle shared by the pair; it can shift the lock, not split the two families. (4) The fitted rim is not uniform: pooled harmonic amplitudes > 0.1 at L ≥ 1, with the h₁ phase tracking the sum-lock axis. If (3) fails downward — the atom-scope fuzz reaches 1.000 at K = 6 — the grain was the disc's edge, and October's "the grain is the birth's" is right for the wrong row: the birth's λ₀, not its ξ.

## C. The axis the beads sit on (no compute — two existing files)
`terrain_ord3\census_structure_v2_nsh32_grid0.txt` and `terrain_ord2\census_structure_v2_nsh32_grid0.txt` carry the (2,2) pixel's phase φ for the cut terrains. Prediction: φ moves with the cut (the axis is the terrain's, set by the phases of the heavy difference beats), by more than the 4.5° spread seen across samplings. If it does not move, the axis belongs to the lock or to v, and "vary v" becomes the next computation.

## D. The breed's angles at fixed order (not coded today)
The terrain's 23 beats are oblique plane waves in (λ, ξ_A, ξ_B), 35–75° off the travel axis in two families (difference beats steep, sum beats shallow). Changing the obliquities at fixed order means re-breeding a terrain with a different (m, n) term list — the breeding pipeline, not the census pipeline. Queued behind A–C; pre-register before breeding.

## Run order
`RUN_BULGE2.bat` now (60 GiB, histograms only, 2–5 h). `RUN_FUZZ.bat` now or after (100 GiB; 7.5 h build + LPs); both can run beside job 399. Send: `bulge2_console.txt`'s table, `pipeline_fuzz_console.txt`'s table, and the two cut-terrain structure transcripts for C.
