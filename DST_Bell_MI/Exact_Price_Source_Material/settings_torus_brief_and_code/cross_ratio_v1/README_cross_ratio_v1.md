# cross_ratio_v1.py — what to export, and what it computes

Companion to arXiv:2608.18886 (Sec. 10, Design A) and settings_torus_brief_v2_4.pdf.
Python 3 + numpy. `python cross_ratio_v1.py --selftest` runs it on synthetic data with
planted fringes and planted imbalance; `--example table` / `--example peaks` print
input templates.

## The quantity

For a two-channel, passive-basis apparatus, R4(a,b) is the sum of the four coincidence
channels in one basis combination. The cross-ratio

    X = R4(HV,HV) · R4(DA,DA) / [ R4(HV,DA) · R4(DA,HV) ]

cancels every factorizable apparatus response — splitter ratios, path and fibre losses,
detector efficiencies — without knowing any of them. Quantum mechanics with a
factorizable response: X = 1. The Sec. 5 detection profile: X = 1.278. A mixture with
detection share s: X(s), inverted by the script to report s with its error and a 95%
upper bound. 10^8 fourfold coincidences in all give 5σ at 0.1% peak-to-peak.

Unknown analyser frame rotations of δ (e.g. from fibre neutralisation) scale the
signal by cos 4δ — 94% at 5°, 77% at 10° — and do not bias the null.

## What to export

**Four detectors per side (one basis per detector pair):** per run, the 2×2 coincidence
table for each of the four basis combinations (A outcome × B outcome, 16 numbers), and
the singles on each of the eight detectors. Optionally the run duration and coincidence
window, for an accidentals subtraction. No timestamps are needed. Two runs at different
pump powers are the one useful extra: multi-pair events add a term of the same sign as
the signal, of relative size ~Rτ/8.

**Time-multiplexed basis (two detectors per user, Bristol-style):** per user pair, the
coincidence counts in the central peak (same basis) and in each of the two side peaks
(cross basis), each split by detector pair where possible (2×2), plus the singles on each
user's two detectors. Label the side peaks by which user took the long path. With three
or more users the per-user long/short path efficiency ratios x_u are solved from the
network itself — side-peak ratios give x_b/x_a, central peaks overdetermine the rest —
with a χ² consistency check. With two users only, the script reports an upper bound on
X (exact if either path ratio is 1), or a measurement if the two ratios are supplied.

## What it does with the one non-factorizable term

Outcome-dependent channel efficiency gives R4 = ¼[S_A S_B + D_A D_B E] (paper Eq. 6),
which does not cancel in X. The script takes δ = D/S for each side and basis from the
channel singles (the local states of a maximally entangled pair are maximally mixed, so
singles ratios are efficiency ratios), takes E from the efficiency-reweighted
coincidence table, divides each R4 by (1 + δ_A δ_B E), and reports X before and after.
As a cross-check it also reads δ off the coincidence marginals of the cross-basis
combinations, where E = 0 makes that estimate exact; the two δ estimates should agree.

## Self-test results (synthetic, 2.5×10^7 coincidences per setting)

| planted share | raw X | calibrated X | recovered share |
|---|---|---|---|
| 0 | 1.0795 (20% imbalance bias) | 0.9996 ± 0.0008 | 0.000 ± 0.003 |
| 0.01 | 1.0822 | 1.0021 ± 0.0008 | 0.008 ± 0.003 |
| 0.10 | 1.1066 | 1.0247 ± 0.0008 | 0.098 ± 0.003 |
| 1.00 | 1.3795 | 1.2774 ± 0.0008 | 0.998 ± 0.003 |

Network mode, four users with planted path ratios 0.75 / 0.90 / 0.80 / 0.85: recovered
0.747 / 0.896 / 0.797 / 0.846 (±0.004–0.005), X = 0.998 ± 0.002 for a null, share
0.05 → 0.041 ± 0.009, share 1 → 0.991 ± 0.009.
