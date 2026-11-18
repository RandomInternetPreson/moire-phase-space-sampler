# From raw time-tags to the cross-ratio: two scripts

Companion to arXiv:2608.18886 (Sec. 10, Design A) and settings_torus_brief_v2_5.pdf.
Python 3 + numpy. Each script runs its own self-test on synthetic data.

    time-tag file(s)  --timetag_to_table_v1.py-->  table.json  --cross_ratio_v2.py-->  X, share, bound

## 1. timetag_to_table_v1.py — counting only

Reads a time-tag stream (a time and a channel number per detection), subtracts a fixed
delay per channel, and counts coincidences inside a window. It interprets nothing.

- **Formats:** text (two numeric columns, any order, any delimiter, any time unit),
  binary fixed-length records (described by a numpy dtype), or .npz arrays. Files are
  read in chunks, so a time-ordered file larger than memory is fine. Several files are
  consecutive segments of a run by default, or simultaneous streams on one clock with
  `"merge": "parallel"`. Two taggers on separate clocks must be synchronised first, as
  the lab's own sifting software already does.
- **Mode "table"** (four detectors per side): the 2x2 coincidences for each of the four
  basis combinations, the eight detector singles, run duration and window.
- **Mode "peaks"** (time-multiplexed basis, two detectors per user): per user pair, the
  2x2 coincidences in the central peak and in the two side peaks at plus and minus the
  long-path delay, and each user's two detector singles. Convention: tau = t_b - t_a,
  so "side_ab" (A short, B long) is the peak at +delay.
- **`--scan config.json`** histograms t_B - t_A for every channel pair on the first
  chunk and prints where the peak sits, to set `delays_s`. A channel's delay is
  subtracted from its times, so delay_B - delay_A = peak position.
- **`--example-config table|peaks`** prints a config template.

Self-test (`--selftest`, about 15 s): 7 million synthetic tags with planted per-channel
delays, timing jitter, unpaired singles and dark counts. Checks that a run split into
19 chunks gives the identical table to a single pass; that a text file in seconds with
swapped columns gives the identical table to the binary file; that every extracted cell
matches the true pair list times the window capture fraction plus accidentals (worst
cell 0.3 sigma); that the delay scan recovers the planted delays to one bin; and, if
cross_ratio_v2.py is in the same folder, that the whole chain recovers a planted
detection share of 0.10 as 0.102 +/- 0.019 and planted path ratios 0.8 / 0.9 / 0.7 as
0.802 / 0.902 / 0.702.

## 2. cross_ratio_v2.py — what changed from v1

v1 removed the one non-factorizable apparatus term (paper Eq. 6, outcome-dependent
channel efficiency) using delta = (eta+ - eta-)/(eta+ + eta-) taken from the channel
**singles**. That is correct only if a detector's singles are proportional to its
efficiency. Dark counts, or any background common to both detectors of a basis, pull
the singles ratio toward 1 and under-correct. With 20% imbalance per side:

| background / signal singles | planted share | v1 share | v2 share |
|---|---|---|---|
| 0 | 0 | 0.000 | 0.000 |
| 0.01 | 0 | 0.004 | 0.000 |
| 0.1 | 0 | 0.052 | 0.000 |
| 1 | 0 | 0.232 | 0.000 |
| 10 | 0 | 0.307 | 0.000 |
| 10 | 0.10 | 0.408 | 0.098 |

The size of the effect scales with the product of the two imbalances, and its sign
with their relative sign, so it can fake a signal or hide one. v1 did print the
coincidence-marginal estimate next to the singles estimate, so the disagreement was
visible in its output, but it used the singles value.

v2 calibrates from the **coincidence tables alone**. A 2x2 table
n_ij = K (1 + s_i dA)(1 + s_j dB)(1 + s_i s_j E) determines dA, dB and E in closed form
for any E, and its sum is 4K (1 + dA dB E). The deltas come from the two cross-basis
tables, where all four cells are well populated; E comes from each table itself; each
fourfold sum is divided by (1 + dA dB E). Background singles do not enter. The singles
estimate is kept as a printed cross-check, with a flag when the two disagree by more
than 0.01 or five standard errors. Singles are now optional: they are needed only for
that cross-check and for the accidentals subtraction.

Everything else is as in v1: the raw cross-ratio, the Sec. 5 share inversion and 95%
upper bound, the peaks mode with the network self-calibration of per-user long/short
path ratios, the two-user upper bound. On clean inputs v1 and v2 agree to the fourth
decimal in X.

## What to export (unchanged)

Four detectors per side: the 16-entry coincidence table and, if convenient, the eight
detector singles, run duration and window. Time-multiplexed: per user pair, central
and side-peak counts split by detector pair. Or the raw time-tags, and the first
script does the counting.
