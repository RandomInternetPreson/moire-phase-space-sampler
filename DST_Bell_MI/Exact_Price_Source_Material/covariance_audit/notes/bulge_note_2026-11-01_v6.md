# The bulge — dated note v6 (Nov 1, 2026, evening): sweep 2, the axis's owner, the siding, and the crack bracketed

Supersedes v5 (same day; adds §12b from `bulge_stair18_2.txt`, `bulge_stair24_2.txt`, `bulge_stair36_2.txt`, `bulge_stair12_2_K6.txt`, and corrects §12a's reading). v5 added §12a; v4 §10–13. Adds §10 (sweep 2, job 404, `bulge2_console.txt`; pre-registration `second_order_preregistration_2026-11-01_v1.md` §A), §11 (the axis, §C of the same), §12 (siding: `bulge_stair12_4.txt`, `bulge_shingle12_4_6.txt`), and the scorecard §13. v3 added §3a, the wall threshold. v1 was Oct 31. Job 403 complete (`pipeline_bulge_console.txt`; 14 profiles; the slowest, tanh:4, 1052 min). Companion: `bulge_test_preregistration_2026-10-31_v1.md` (written before the run); scorecard in §6. Tier: computed, Aaron's machine, lattice sampling (32×32), ε = 0.01, full terrain. Transcripts: `bulge_<tag>.txt` × 14, `bulge_sign_eps0.txt`, `bulge_dome0p2.txt`, `bulge_dome0p05.txt`, `bulge_tanh5_full.txt`, `bulge_tanh6_full.txt`, `bulge_tanh7_full.txt`; jsons `census_profile_v1_<tag>.json` × 14, `census_profile_v2_sign.json`, `census_profile_v2_dome0p2.json`, `census_profile_v2_dome0p05.json`, `census_profile_v2_tanh5.json`, `census_profile_v2_tanh6.json`, `census_profile_v2_tanh7.json`.

## 1. The table (MD/floor; `m` = singlet unreachable at ε = 0.01, minimum miss shown)

| profile | wall width (\|D\| < 0.9) | S_max/2√2 | relaxed miss | full | K=6 | K=8 | K=10 | K=12 | K=14 |
|---|---|---|---|---|---|---|---|---|---|
| **sign** (the step) | 0° | 1.4142 | 0.0000 | **1.0000** | 1.6996 | 1.4290 | 1.2579 | 1.1122 | **1.0000** |
| tanh:8 | 10.6° | 1.4142 | 0.0001 | **1.0476** | 2.5615 | 2.0004 | 1.7401 | 1.5451 | 1.3038 |
| tanh:4 | 21.5° | 1.4105 | 0.0131 | m0.0131 | m0.047 | m0.033 | m0.026 | m0.022 | m0.017 |
| tanh:2 | 41° | 1.3516 | 0.1339 | m0.167 | m0.338 | m0.252 | m0.226 | m0.211 | m0.194 |
| tanh:1 | 57° | 1.2306 | 0.2579 | m0.401† | m0.288 | m0.277 | m0.571† | m0.508† | m0.470† |
| smooth (cos 2θ) | 64° | 1.1263 | 0.3154 | m0.502† | m0.347 | m0.336 | m0.329 | m0.325 | m0.320 |
| dome:0.05 (v1) | 0° (cap sags 2.8% at 28°) | 1.3982 | 0.0206 | m0.0206 | m0.034 | m0.030 | m0.028 | m0.026 | m0.024 |
| dome:0.2 (v1) | 0° (cap sags 11% at 28°) | 1.3511 | 0.0817 | m0.082 | m0.141 | m0.123 | m0.111 | m0.102 | m0.091 |
| notch:3:0.25 | 0° (notched cap) | 1.2802 | 0.1397 | m0.141 | m0.217 | m0.263† | m0.227 | m0.200 | m0.174 |
| notch:5:0.25 | | 1.2793 | 0.1890 | m0.190 | m0.214 | m0.219 | m0.202 | m0.195 | m0.192 |
| notch:7:0.25 | | 1.2813 | 0.2131 | m0.213 | m0.230 | m0.231 | m0.222 | m0.217 | m0.215 |
| notch:5:0.1 | | 1.3599 | 0.0506 | m0.051 | m0.088 | m0.076 | m0.070 | m0.063 | m0.057 |
| anotch:3:0.25 | | 1.2798 | 0.2347 | m0.235 | m0.242 | m0.242 | m0.241 | m0.240 | m0.238 |
| anotch:5:0.25 | | 1.3070 | 0.2345 | m0.235 | m0.242 | m0.242 | m0.241 | m0.240 | m0.237 |
| anotch:7:0.25 | | 1.2912 | 0.2345 | m0.235 | m0.241 | m0.241 | m0.239 | m0.238 | m0.237 |
| anotch:5:0.1 | | 1.3709 | 0.0975 | m0.098 | m0.099 | m0.099 | m0.099 | m0.098 | m0.098 |

† Non-monotone in K (a larger class cannot miss by more): the min-miss LP with the CHSH equality is numerically unreliable far from feasibility; the *relaxed miss* column (no CHSH equality, status optimal, verified from the reconstructed census) is the number to read for those rows. All † rows are ≥ 14% off the singlet either way. Where the ladder is monotone (tanh:4, tanh:2, dome, notch:5, anotch) the certificates stand.

## 2. The control — kept exactly
The rebuilt histograms contracted with the step reproduce the deposited banks to 2.2e-15; full 1.0000; K-curve 1.6996 / 1.4290 / 1.2579 / 1.1122 / 1.0000 against the deposited 1.700 / 1.429 / 1.258 / 1.112 / 1.000; same skeleton ((2,0), (2,2) leading, phases 162–177°, (0,4) third). **And the floor is attained at ε = 0** (v1 §2; `bulge_sign_eps0.txt`): the price LP with a zero-width band is feasible at MD/floor 1.0000. The ε-surface's corner is the floor.

## 3. The wall must be sharp — and a softer wall asks the birth for *more* grain, not less
The only shaped reader that reaches the singlet at all is tanh:8 (wall 10.6° wide): its floor price is 1.0476, and its K-curve lies **entirely above the step's** — 2.56 at K = 6 (step 1.70), 1.30 at K = 14 (step 1.00); it has not reached its own full price by K = 14. A 21.5° wall (tanh:4) already misses the singlet by 1.3%; wider walls by 13–32%. So: (i) the wall is sharper than ~20°, and ~11° already costs 5% of floor; (ii) blurring the reader does not let the census be smoother — the census must be *finer* to compensate. The bulge does not absorb the grain; a shaped bulge demands more of it. **Prediction 4's premise ("blunting lowers the smoothness premium") failed — the premium rose from 1.70 to 2.45; its conclusion ("no profile brings the absolute K = 6 price toward 1.000") is kept, at 2.56.**

### 3a. The wall, bracketed (full census only; ε = 0.01)
| wall width (\|D\| < 0.9) | reader | relaxed miss | floor price | excess over floor |
|---|---|---|---|---|
| 0° | step | 0.0000 | 1.0000 | 0 |
| 10.6° | tanh:8 | 0.0001 | 1.0476 | 0.048 |
| 12.1° | tanh:7 | 0.0003 | 1.1686 | 0.169 |
| 14.2° | tanh:6 | 0.0012 | 1.4568 | 0.457 |
| 17.1° | tanh:5 | 0.0039 | 2.3514 | 1.351 |
| 21.5° | tanh:4 | 0.0131 | — (singlet lost) | ∞ |
The excess over the floor roughly triples for every ~2° of wall width (0.048 → 0.169 → 0.457 → 1.351): quasi-exponential. The singlet is lost between 17° and 21.5°. A wall of ≈ 10° is the widest the model tolerates at a 5% price; the step is the only reader at the floor. Symmetrization held in all three (marginals ≤ 3e-16 at identical price).

## 4. The cap must be flat (domes, v1) and may carry no ripples (notches)
Every cap feature tried kills the singlet outright, census or no census, at any bandwidth: domes sagging 2.8% at 28° miss by 2.1% (linear in the sag); notches of depth 0.25 miss by 14–21%, of depth 0.1 by 5%; the K-ladders only walk down toward the full-census limit. **The anti-pair is dead by an identity:** at the aligned setting the two readers see the same θ, and with Bob's notches between Alice's the product of the two caps is g_A g_B = 1 − E + E²(1 − cos²4Jθ)/4 ≤ (1 − E/2)², so |E| ≤ (1 − E/2)² and the miss at perfect anticorrelation is at least 1 − (1 − E/2)² = 0.2344 (E = 0.25) / 0.0975 (E = 0.1) — the computed misses, 0.2345–0.2347 and 0.0975, to four digits, for every J. Twins carrying complementary features can never both be at full height for the same λ. Prediction 5 overshot in the "worse" direction: not pricier, infeasible; the anti-pair is systematically worse (an identity), not "not systematically better".

## 5. The smooth bulge — prediction 2 failed as stated, kept in substance
A cos 2θ reader *does* reach the CHSH value (S_max = 1.126 × 2√2; the lock is sharp enough for four settings) — the pre-registered "S_max < 2√2" is **failed, recorded**. It cannot reach the singlet's shape over the register: relaxed miss 0.315, the worst in the table.

## 6. Scorecard (pre-registration of Oct 31)
1. Control reproduces the deposited numbers — **kept** (all K, to 4 digits; banks to 1e-15).
2. Smooth bulge cannot reach 2√2 — **failed as stated** (S_max 1.126 × 2√2); singlet unreachable by 31%.
3. Blunting raises the full price monotonically — **kept, stronger**: 1.000 → 1.048 (10.6° wall) → infeasible beyond.
4. Premium falls, absolute K = 6 stays ≥ 1.5 — **premise failed** (premium 1.70 → 2.45), **conclusion kept** (2.56).
5. Notches raise both prices; anti not systematically better — **overshot**: infeasible; anti systematically worse by an identity.
6. Symmetrized census has zero marginals at the same price — **kept** (17 of 17 runs; raw vertex marginals up to 0.78, symmetrized ≤ 3e-16, correlators unchanged to ≤ 2e-14).
Added after pre-registration: the dome family (direction predicted from the toy run; magnitude — infeasible — not).

## 7. Reading (for the casebook; tier: reading on computed)
Aaron's question was whether the bulge could have features — a blunted-spiky bead, the twin carrying the anti-version — that carry the 7–14 band instead of the birth. The answer is no in every direction tried, and the geometry of the no is specific. **The reader is a mesa**: a wall sharper than ~17° (10° costs 5% of floor, 14° costs 46%, 17° costs 135%; §3a) and a tabletop flat to the precision of the singlet (1% for ε = 0.01; 0.1% at polarimetric precision), with no ripples, and with the twin's tabletop identical, not complementary. Where λ locked inside the pass band is invisible to the click; only which side of the wall counts. The casebook's "orientation and nothing else" (3.17, for the static D(m)) is now a computed property of the dynamical model with the lock included — the lock does not rescue a shaped reader. And the one shaped reader that survives (the 10.6° wall) *raises* the census's bandwidth demand: the grain is the birth's, and a blurred bulge needs more of it. In drop language: not a bead, a mesa; and the sharper the mesa's edge, the less the birth has to know.

## 10. Sweep 2 — the mesa's second-order geometry (job 404; `census_profile_v3.py` 5a2eceeb, `run_bulge2_v1.py` 3159b1b6; full + K = 6, 8, 10; ε = 0.01; lattice)

| idea | reader | port rule | S_max/2√2 | relaxed miss | full | K=6 | K=8 | K=10 | joint MD | ⟨A⟩ sym | dE sym |
|---|---|---|---|---|---|---|---|---|---|---|---|
| cap tilted 5% | tilt:0.05 | kept | 1.3758 | 0.0305 | m0.031 | m0.045 | m0.040 | m0.037 | — | 3e-17 | 2e-16 |
| cap tilted 10% | tilt:0.1 | kept | 1.3381 | 0.0612 | m0.061 | m0.089 | m0.079 | m0.074 | — | 2e-17 | 3e-16 |
| front 3.6° / rear 7.2° | awall:3.6:7.2 | **broken** | 1.4142 | 0.0019 | **2.083** | m0.023 | m0.015 | m0.011 | 7.22 | 0.023 | 0.012 |
| front 3.6° / rear 14.3° | awall:3.6:14.3 | broken | 1.4127 | 0.0350 | m0.035 | m0.093 | m0.071 | m0.059 | — | 0.11 | 0.13 |
| front 1.8° / rear 7.2° | awall:1.8:7.2 | broken | 1.4142 | 0.0011 | **1.791** | m0.020 | m0.013 | 3.430 | 7.20 | 0.031 | 0.019 |
| twin upside down, 3.6/7.2 | mawall:3.6:7.2 | broken | 1.4141 | 0.0027 | 2.128 | m0.024 | m0.015 | m0.011 | 7.14 | 0.021 | 0.012 |
| twin upside down, 3.6/14.3 | mawall:3.6:14.3 | broken | 1.4012 | 0.0356 | m0.036 | m0.101 | m0.079 | m0.066 | — | 0.11 | 0.13 |
| kinked wall, top 1.8° / bottom 7.2° | bwall:1.8:7.2 | kept | 1.4142 | 0.0010 | **1.450** | m0.015 | 3.129 | 2.559 | 6.97 | 2e-16 | 2e-16 |
| kinked wall, top 7.2° / bottom 1.8° | bwall:7.2:1.8 | kept | 1.4142 | 0.0003 | **1.333** | m0.014 | 2.980 | 2.466 | 7.06 | 2e-16 | 3e-16 |
| wall follows approach phase, 5° | shift:5 | kept | 1.4142 | 0.0000 | (0.542)† | 2.073 | 1.569 | 1.283 | 6.85 | 0.064 | 0.12 |
| 10° | shift:10 | kept | 1.4142 | 0.0000 | (0.586)† | m0.019 | 2.576 | 1.841 | 6.81 | 0.11 | 0.18 |
| 20° | shift:20 | kept | 1.4142 | 0.0000 | (1.937)† | m0.134 | m0.064 | m0.032 | 6.37 | 0.28 | 0.27 |

† Not a price. The `shift` reader reads the approach phase ξ as well as λ, so the Bell hidden variable is (λ, ξ) and the measurement dependence must be measured on the joint: 6.4–6.9× the floor. The λ-marginal "price" below 1.000 is the LP exploiting a variable the MD rows do not see; Hall's floor is not touched. Recorded as a method result: a ξ-reading detector moves the accounting, it does not buy the singlet.

Readings. (i) **A tilted cap is a cap that is not flat**: infeasible, miss ≈ 0.61 c (dome: linear in the sag as well). (ii) **A front and a back**: priced by the wider wall and worse — front 3.6°/rear 7.2° costs 2.08 where both-3.6° (tanh:8) cost 1.05 and both-7.2° (tanh:4) is infeasible; sharpening the front to 1.8° with the same soft rear brings 1.79, so the front edge is not idle. Breaking the port rule has the predicted physical consequence: the +90° symmetrization is no longer free (correlators move by 0.01–0.13) and the singles are polarized (⟨A⟩ 0.02–0.11) — a front-backed mesa predicts polarized singles, which the lab does not see. (iii) **The twin born upside down** changes the prices by 2–3%: nothing. (iv) **Kinked walls** (port rule kept) price between the sharp wall and the soft half alone (1.33–1.45); the two kink orientations differ by 8%. (v) **Every shaped reader again raises the band-limited prices**: K = 6 infeasible for all of them, K = 8–10 at 2.5–3.4 where the step pays 1.43–1.26. The grain is not in the mesa's second-order geometry either.

## 11. The axis the beads sit on — owner found in existing files (§C of the pre-registration)
`terrain_ord3\census_structure_v2_nsh32_grid0.txt` (SHA-256 df086304…) and `terrain_ord2\…` (72c8bb84…), both lattice, both at the floor:

| | full terrain (6 censuses) | ord3 (m, n ≤ 3; 15 terms) | ord2 (m, n ≤ 2; 6 terms) |
|---|---|---|---|
| top pixel | (2,2) 0.33–0.38 | (2,2) 0.328 | **(4,4) 0.324** — (2,2) not in the top five |
| sum-lock axis ξ_A + ξ_B | 39.5° mod 90° (φ 157–159°) | 36.6° mod 90° (φ +146°) | 33.7° **mod 45°** (φ₄₄ −22.5°): period halved |
| Alice marginal (2,0) phase | +148…+154° | +153° | **−19°** (near sign reversal) |
| anti-diagonal (relative phase) share | at the null (z ≈ 0) | z = +2.4; (2,−2) 0.220 | **z = +3.7; (2,−2) 0.295** |
| weight at s = 1.0 | ≈ 0.3 | 0.29 | **0.000** (six terms cannot lock the shallow mesh) |

Prediction C ("φ moves with the cut by more than the 4.5° sampling spread"): **ord3 — not kept** (3°, inside the spread); **ord2 — kept, and more**: the skeleton reorganizes rather than drifting. Together with the Oct 16 perturbation run (the (2,2) pixel tracks the (2,2,+) beat's phase; removing that beat leaves a weaker (2,2) at another phase): **the axis belongs to the terrain's order-3 structure** — the (2,2,+) sum beat sets its phase, and the heavy difference beats (1,3,−), (2,3,−), (3,2,−) — the steep oblique waves, 58–75° off the travel axis — keep it on top; strip them and the census falls back on the shallow sum beats (45° period, hence (4,4)) and begins to use relative phase. First pairing of the mesh's two obliquity families with the census's two structural families; testable directly by the breed-angles job (§D of the pre-registration).

## 12. Siding — the first shaped reader at the floor (`census_profile_v4.py` a3413fc5; full census only)

| reader | port | relaxed miss | full | joint MD | ⟨A⟩ sym | dE sym |
|---|---|---|---|---|---|---|
| **stair:12:4** — 12° wall as a staircase of 4 risers (5 click levels) | kept | **0.0000** | **1.0048** | 6.98 | 2e-17 | 2e-16 |
| shingle:12:4:6 — the same, rows staggered 6° with sin 2φ | kept | 0.0000 | (1.0045)† | 6.77 | 0.051 | 0.059 |
| tanh:7 — 12.1° smooth wall (§3a) | kept | 0.0003 | 1.1686 | — | — | — |

**Pre-registered prediction failed, recorded:** "a 12° staircase prices like a 12° smooth wall (≈ 1.17)". It prices at **1.0048**, four times closer to the floor than the 10.6° smooth wall (1.048), and it reaches the singlet *exactly* (relaxed miss 0, like the step; every smooth wall misses by ≥ 0.0001). Reading (after the data, marked so): what the floor demands of the mesa is **sharp risers, not a narrow wall** — the staircase keeps the step's high harmonics in λ, the smooth wall throws them away. Width tolerance and sharpness tolerance are different specs, and §3a measured the wrong one as "the" wall tolerance. Whether a staircase reader also changes the *band-limited* prices — the grain question — is not yet run (K = 6, 8, 10 for stair:12:4) and is the next command. The shingle's stagger buys nothing beyond the stair (1.0045 vs 1.0048 on the marginal) and reads ξ († as in §10).

### 12a. Siding, second round — width vs risers vs the grain (full census unless marked; ε = 0.01)

| reader | wall | risers / levels | relaxed miss | full | K=6 | K=8 | K=10 |
|---|---|---|---|---|---|---|---|
| step (control) | 0° | — | 0.0000 | 1.0000 | 1.6996 | 1.4290 | 1.2579 |
| **stair:12:2** | 12° | 2 / {+1, 0, −1}: a 6°-wide coin-flip zone centred on each edge | 0.0000 | **1.0000** | — | — | — |
| stair:12:4 | 12° | 4 / 5 levels | 0.0000 | 1.0048 | 2.524 | 1.920 | 1.679 |
| stair:12:8 | 12° | 8 / 9 levels | 0.0000 | 1.0046 | — | — | — |
| stair:24:4 | 24° | 4 / 5 levels (6° plateaus) | 0.0005 | **2.657** | — | — | — |
| tanh:7 (smooth, §3a) | 12.1° | — | 0.0003 | 1.169 | — | — | — |
| tanh:8 (smooth, §3a) | 10.6° | — | 0.0001 | 1.048 | 2.562 | 2.000 | 1.740 |

- **The grain is not in the siding.** stair:12:4's band-limited prices (2.52 / 1.92 / 1.68) lie above the step's (1.70 / 1.43 / 1.26) and within 2% of tanh:8's (2.56 / 2.00 / 1.74): a shaped wall of either kind asks the birth for more structure. Pre-registered (§12): kept.
- **Riser count is irrelevant at fixed width** (4 vs 8 risers: 1.0048 vs 1.0046); **width is not**: 12° → 24° at 4 risers costs 1.00 → 2.66, the same steep rise the smooth walls showed (§3a).
- **A coin-flip crack at the edge is free.** The three-level reader — definite clicks except within ±3° of the 45° boundary, where the outcome is a fair coin — attains the floor **exactly** (1.0000, relaxed miss 0) with marginals that symmetrize to 2e-17. The 5- and 9-level readers of the same width cost 0.5%; the smooth wall of the same width costs 17%. So the floor's demand on the wall is not "sharp" in the sense of narrow, and not "no intermediate levels": it is **no graded slope**. A plateau at zero is a loss of information at the edge and is tolerated; a graded value is wrong-shaped information and is paid for.
- Reading (after the data, marked so): the mesa's edge may be a crack — a band around the boundary where the reader does not decide — and the crack's width is a new tolerance with its own number, not yet bracketed: free at 6°; the 24°/4-riser wall (which contains a 6° zero plateau plus two 6° half-plateaus) is not free, but that is a different reader. Next: `stair:18:2`, `stair:24:2`, `stair:36:2` full-only (the crack width at which the floor is lost), and `stair:12:2` at K = 6 (does a crack change the band-limited price?). Pre-registered: the crack is free to somewhere between 12° and 24° and then rises steeply; stair:12:2 at K = 6 ≥ 1.70.

### 12b. The crack bracketed, and §12a's reading corrected

| three-level reader | crack (coin-flip zone per edge) | relaxed miss | full | K=6 |
|---|---|---|---|---|
| stair:12:2 | 6° | 0.0000 | **1.0000** | **2.2105** |
| stair:18:2 | 9° | 0.0000 | 1.1636 | — |
| stair:24:2 | 12° | 0.0000 | 1.9360 | — |
| stair:36:2 | 18° | 0.0275 | infeasible (min miss 0.028) | — |

- **The crack is free at 6° and not beyond:** 9° costs 16%, 12° costs 94%, 18° loses the singlet. Pre-registered "free to somewhere between 12° and 24°": **failed** — the free range ends between 6° and 9°.
- **The crack does not lower the band-limited price:** stair:12:2 at K = 6 is 2.21 (step 1.70). Pre-registered ≥ 1.70: kept. The grain is the birth's, under every reader tried (now 30).
- **§12a's reading was too strong and is withdrawn in that form.** It read "a zero plateau is tolerated, a graded value is charged". Comparing readers by the *area* they lose against the step per edge (∫|D − step| dθ, in degrees of full loss): the three staircases of lost area 6° price 1.000 / 1.005 / 1.005 (2 / 4 / 8 risers); the smooth walls of lost area 5.0° and 5.7° (tanh:8, tanh:7) price 1.048 and 1.169; the 9° crack prices 1.16 where the smooth wall of area 8° (tanh:5) prices 2.35; the 12° crack 1.94 where the smooth wall of area 10° (tanh:4) has lost the singlet. So the correct statement is: **at equal lost area, a sharp-edged loss (any number of risers, including a flat zero) is far cheaper than a smooth graded loss; and about 6° of sharp-edged loss per edge is tolerated at ≤ 0.5%, after which the price rises steeply whatever the shape.** Separately, loss *on the cap* is in a different class: tilt (area 2.3°), dome (3°) and notch (4.5°) are infeasible with less lost area than the free staircases — where the loss sits matters more than how much.
- The spec, as the machinist would write it: tabletop flat (no feature of any size); wall vertical faces, no bevel; kerf at the edge ≤ 6° of lost area, clean-edged; both edges and both twins alike; no dependence of any of it on the approach phase.

## 13. Scorecard for the second-order pre-registration (§A, §C)
A-tilt **kept** (infeasible; 0.61c). A-awall **kept** (priced by the wider wall, worse; polarized singles). A-mawall **kept** (no gain). A-bwall **kept** (between; orientations differ 8%). A-shift **kept** (marginal "price" below the floor is not a price; joint MD 6.4–6.9). A-grain **kept** (no reader brings K = 6 below 1.70; all raise it). C **half**: ord3 held (3°), ord2 reorganized. Siding (added after the pre-registration, predicted ≈ 1.17): **failed** — 1.0048.

## 8. Open after today
- ~~Threshold of the wall~~ — done, §3a.
- ~~Siding~~ — done, §12a–b; the mesa's spec is written. **Readers tried: 30**; none lowers the band-limited price.
- **The swirl (Aaron, 1 Nov evening):** the disc's own rotation is v (λ̇ = v on a flat register; v = 15.23); a time-varying swirl v(t) is a new dynamics slot (modified `flight_time`, one rebuild per variant). Pre-register before building; the first discriminant is still constant v varied (does the 6° grain follow v?).
- **Breed angles at fixed order** (§11 makes it concrete): a terrain with the order-3 difference beats re-aimed.
- The *same* readers on the 4096 random sampling (the smooth-structure rule: a lattice result on a *smooth* quantity is replicated on random offsets — here the prices are smooth-structure numbers; the infeasibilities are robust by margin). Needs a 4096 histogram store (4× the lattice build).
- Malus for one photon through two polarizers in sequence still needs a single-photon station (v1 §"what the test cannot do"); unchanged.
- For the paper: the ε = 0 floor (one sentence, next revision). For Zenodo version 6: this sweep + K = 20 (job 399, running) + the three command-line runs.

## 9. Files for deposit v6 (bulge part)
Scripts `pipeline_bulge_v1.py` (d166dac8), `census_profile_v1.py` (b1ed2db3), `census_profile_v2.py` (bc2f617b), `RUN_BULGE.bat` (d26f3131). Transcripts: `pipeline_bulge_console.txt`; `bulge_sign.txt`, `bulge_smooth.txt`, `bulge_tanh1.txt`, `bulge_tanh2.txt`, `bulge_tanh4.txt`, `bulge_tanh8.txt`, `bulge_notch3_e25.txt`, `bulge_notch5_e25.txt`, `bulge_notch7_e25.txt`, `bulge_anotch3_e25.txt`, `bulge_anotch5_e25.txt`, `bulge_anotch7_e25.txt`, `bulge_notch5_e10.txt`, `bulge_anotch5_e10.txt`; `bulge_sign_eps0.txt`, `bulge_dome0p2.txt`, `bulge_dome0p05.txt`, `bulge_tanh5_full.txt`, `bulge_tanh6_full.txt`, `bulge_tanh7_full.txt`. Data: the 20 jsons named above. **Not** the histogram store (`hist_nsh32_grid0_nang161_nw2048.npz`, size as on disk) unless Zenodo's limit allows it; it is reproducible from `pipeline_bulge_v1.py` and the deposited `greedy_results.json`, and the bit-for-bit check against the deposited banks is in `bulge_sign.txt`. Cosmetic slip in all transcripts: the profile header says "step 15"; the 13 printed angles are 13.85° apart; values correct.
