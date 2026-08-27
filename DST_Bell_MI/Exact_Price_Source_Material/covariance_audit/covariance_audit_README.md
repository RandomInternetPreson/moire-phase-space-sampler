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
