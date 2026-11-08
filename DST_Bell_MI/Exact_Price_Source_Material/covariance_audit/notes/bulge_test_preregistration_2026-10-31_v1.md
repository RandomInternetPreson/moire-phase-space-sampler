# The bulge test — pre-registration (Oct 31, 2026, v1; written before the run)

**Origin.** Aaron (31 Oct): is there a transition we are not seeing, bulge → disc? Could the bulge have features — a blunted-spiky bulge, with the twin photon carrying the anti-version — that carry the 7–14 band instead of the birth? **Translation into the model.** The disc is λ under Eq. (4); the bulge is what reads the locked λ: the detection rule A = sign(cos 2(λ − a)), B = sign(cos 2(λ + π/2 − b)). In the model as written the bulge is a *step*: every odd harmonic of λ, no structure beyond the edge. Nothing in September–October varied it, because the banks store the 161 correlator angles already contracted with the sign. This is the first computation in the program that varies the detection rule.

**Code.** `pipeline_bulge_v1.py` (d166dac8) rebuilds the 32×32 lattice trajectories once and stores every λ-histogram as counts (15 cells × 165 settings × 1024 offsets × 256 bins, uint16; `hist_nsh32_grid0_nang161_nw2048.npz`), then runs `census_profile_v1.py` (b1ed2db3) on 14 profiles in parallel; `RUN_BULGE.bat` (d26f3131), one launcher job, 80 GiB, ~4 h build + 4–8 h sweep. Smoke-tested in the container at nsh 8 random and nsh 4 grid (nang 13, nw 256): the step reproduces the stored banks to 9e-16; the +90°/+90° symmetrization leaves every correlator and the MD unchanged to 2e-16 and zeroes the marginals.

**The class.** D(θ) on the station-relative angle θ, |D| ≤ 1, D(θ + 90°) = −D(θ) (the other port). A value strictly inside (−1, 1) is a stochastic local response (Bell's A(λ,a) ∈ [−1, 1]); the correlator kernel is D_A(λ − a) · D_B(λ + π/2 − b) and the MD rows are unchanged (they see only the λ-distribution). Profiles:
- `sign` — the paper's step (control).
- `smooth` — cos 2θ (the fully smooth bulge).
- `tanh:B`, B = 1, 2, 4, 8 — blunted steps tanh(B cos 2θ)/tanh B (edge widths ≈ 30°, 15°, 7°, 3.6°).
- `notch:J:E`, J = 3, 5, 7, E = 0.25 (and J = 5, E = 0.1) — the step with periodic notches of depth E at λ-frequency 4J (features at λ-harmonics ≈ 2J ± 1: 5–7, 9–11, 13–15), the same on both stations.
- `anotch:J:E` — the twin carries the anti-version: Bob's notches sit where Alice's do not.

**Per profile:** S_max (largest CHSH value any census reaches with this bulge; below 2√2 the singlet is unreachable); the relaxed minimum miss over 165 correlators (no CHSH equality); the full-census price (the floor price of the bulge itself); the band-limited prices at K = 6, 8, 10, 12, 14 (positivity at atoms, mass, CHSH exact, ε = 0.01 band; infeasible cells certified by value, as in `census_bandlimit_v4`); the marginals ⟨A⟩(t), ⟨B⟩(t) at all 161 angles for the raw vertex and for the census symmetrized with its (+90°, +90°) twin.

**Pre-registered predictions (to be marked kept / failed when the run lands).**
1. `sign` reproduces the deposited lattice numbers exactly: full 1.000000; K = 6, 8, 10, 12, 14 → 1.700 / 1.429 / 1.258 / 1.112 / 1.000.
2. `smooth` cannot reach the singlet: S_max < 2√2 (the lock is not sharp enough for a cos² reader to give unit correlation).
3. Blunting raises the full price monotonically (tanh:8 < tanh:4 < tanh:2 < tanh:1 in price, all > 1.000): a stochastic reader needs more measurement dependence for the same S.
4. Blunting lowers the *smoothness premium* K6/full (a blunter reader averages the lock position over ξ, so the census's fine structure buys less) — but the *absolute* K = 6 price does not fall below 1.5 for any profile. If it does fall toward 1.000 for some profile, the grain lives in the bulge and the October reading of the fine structure as the birth's is wrong.
5. Notches (same and anti) raise both the full price and the K = 6 price; the anti-pair is not systematically better than the same-pair.
6. Marginals: the raw vertex may carry nonzero ⟨A⟩, ⟨B⟩; the symmetrized census has ⟨A⟩ = ⟨B⟩ = 0 at the same price and the same correlators (exact symmetry of the model). This is a verification, not a prediction with risk.

**What the test cannot do.** Malus's law from the mechanism for *one photon through two sequential polarizers* (casebook open item 3) needs a single-photon station the paper does not define (the terrain is the pair's: its phases carry both settings). The heralded single-photon Malus law *is* the 161-angle correlator band and is enforced at ε. The sequential-polarizer computation stays open and needs a modelling decision (which terrain a lone photon meets) before it is coded.

**Reading rule.** The grain is in the bulge iff some profile's absolute K = 6 price approaches 1.000 while the full price stays near 1.000. A fall in the premium K6/full accompanied by a rise in the full price is *not* the grain moving: it is the reader blurring the terrain and paying for it in measurement dependence.
