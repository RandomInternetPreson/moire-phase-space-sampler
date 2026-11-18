# Cross-ratio analysis for passive-basis links

Scripts for Design A of arXiv:2608.18886 (Sec. 10): the two-setting cross-ratio X of the
fourfold coincidence sum, computed from archived or new two-channel data.

    time-tag file(s)  --timetag_to_table_v1.py-->  table.json  --cross_ratio_v2.py-->  X, share, bound

The folder keeps the name `cross_ratio_v1` so that links already sent out still resolve.
**The current analysis script is `cross_ratio_v2.py`.**

## Current files

| File | What it is |
|---|---|
| `cross_ratio_v2.py` | The analysis. Takes the coincidence table (or the central and side-peak counts of a time-multiplexed network) and returns X, its error, the equivalent Sec. 5 detection share and a 95% upper bound. |
| `timetag_to_table_v1.py` | The counting. Turns raw time-tags (text, binary records or .npz, read in chunks) into the table the analysis reads. `--scan` locates per-channel delays. |
| `README_timetag_and_cross_ratio_v2.md` | What to export, how each script works, and the self-test results. |
| `ex_table_v2.json`, `ex_peaks_v2.json` | Example inputs for the two modes. |

Python 3 and numpy. Each script checks itself on synthetic data:

    python cross_ratio_v2.py --selftest
    python timetag_to_table_v1.py --selftest

## Superseded, kept for the record

`cross_ratio_v1.py`, `README_cross_ratio_v1.md`, `ex_table.json`, `ex_peaks.json` (2 November 2026).

**If you ran v1, use the v2 number.** v1 took the detector-imbalance calibration from the
channel singles. That under-corrects when the singles carry dark counts or any other
background common to both detectors of a basis, and the residual can fake a signal or
hide one. With 20% imbalance per side and a true null, v1 reports a detection share of
0.004, 0.052 and 0.232 for background at 0.01, 0.1 and 1 times the signal singles; v2
reports 0.000 in each case. v2 calibrates from the coincidence tables alone, which
background singles do not touch, and prints the singles estimate beside it with a flag
when the two disagree. On inputs without background the two versions agree to the
fourth decimal in X.

The paper (Sec. 10, through v7) and the brief before v2.6 describe the calibration from
channel singles with the coincidence marginals as a cross-check. `cross_ratio_v2.py`
and brief v2.6 use the coincidence estimate as the calibration and the singles as the
cross-check; `../protocol_sensitivity_v2_3.py` tests the same reversal on the full torus.
