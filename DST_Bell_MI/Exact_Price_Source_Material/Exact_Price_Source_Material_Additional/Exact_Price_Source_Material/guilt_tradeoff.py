#!/usr/bin/env python3
"""
guilt_tradeoff.py - The two forensic signatures traded against each other.

Hybrid family: D_alpha(m) = (1-alpha) + alpha*D*(m) (the exact profile at
full strength alpha=1). For each alpha:
  - rate anisotropy eps(alpha): modulation of TOTAL coincidence rate vs
    relative angle Delta, N(Delta) ~ C_g(2 Delta)  [QM: 0]
  - MD(alpha): minimal measurement dependence to restore exact QM
    correlations (from the rate-constrained frontier LP, already computed
    in frontier_data_v2.npz for exactly this family)
Deliverable: MD floor as a function of the experimentally bounded
anisotropy eps -> a Bell bound driven by a LOGGED observable.
Also verify: (a) the anisotropy is a cos(4 Delta) fringe; (b) it vanishes
identically at the four CHSH angles (the hiding mechanism).
"""
import numpy as np

NX = 8192
x = np.linspace(0, 2*np.pi, NX, endpoint=False)
dat = np.load("sqrt_law_h.npz")
Dstar = lambda m: np.sqrt(m)*np.interp(m, dat["m"], dat["h"]/dat["h"][-1])
mvals = np.abs(np.cos(x))

def rate_stats(alpha):
    g = (1-alpha) + alpha*Dstar(mvals)
    Gf = np.fft.rfft(g)
    Cg = np.fft.irfft(np.abs(Gf)**2, n=NX)/NX
    th = x                                  # theta = 2 Delta
    N = Cg
    eps = (N.max()-N.min())/N.mean()
    # harmonic content of N in Delta: coefficients of cos(4k Delta) = cos(2k theta)
    h1 = 2*np.mean(N*np.cos(2*th))/N.mean()     # cos(4 Delta) fringe amplitude
    h2 = 2*np.mean(N*np.cos(4*th))/N.mean()     # cos(8 Delta)
    # values at CHSH angles Delta = 22.5, 67.5  (theta = 45, 135 deg)
    i45 = NX//8; i135 = 3*NX//8
    chsh_gap = abs(N[i45]-N[i135])/N.mean()
    return eps, h1, h2, chsh_gap

fr = np.load("frontier_data_v2.npz")["interp"]   # rows: (eta, MD) for alpha 0..1 (21)
alphas = np.linspace(0, 1, 21)

print(f"{'alpha':>6} {'eta':>7} {'MD_floor':>9} {'eps(rate)':>10} "
      f"{'cos4D amp':>10} {'cos8D amp':>10}")
rows = []
for a, (eta, md) in zip(alphas, fr):
    eps, h1, h2, gap = rate_stats(a)
    rows.append((a, eta, md, eps, h1, h2))
    print(f"{a:6.2f} {eta:7.4f} {md:9.5f} {eps*100:9.2f}% {h1*100:9.2f}% {h2*100:9.2f}%")

eps_full, h1f, h2f, gap = rate_stats(1.0)
print(f"\nsignature checks at full detection strength:")
print(f"  cos(4D) fringe = {h1f*100:.2f}%  vs cos(8D) = {h2f*100:.2f}% "
      f"-> leading harmonic dominates {abs(h1f/h2f):.1f}:1")
print(f"  rate difference between the two CHSH angle classes "
      f"(D=22.5 vs 67.5): {gap*100:.4f}%  -> invisible to four-point CHSH runs")

rows = np.array(rows)
# the combined bound: MD floor as a function of enforced anisotropy bound
print(f"\nTHE COMBINED BOUND (this family): if experiments bound total-rate")
print(f"modulation below eps, measurement dependence must be at least:")
for eps_b in [0.10, 0.05, 0.02, 0.01, 0.005]:
    idx = np.where(rows[:,3] <= eps_b)[0]
    md_min = rows[idx, 2].max() if len(idx) else rows[0,2]
    # smallest MD consistent with eps <= eps_b = MD at the LARGEST allowed alpha
    ok = rows[rows[:,3] <= eps_b]
    md_req = ok[:,2].min() if len(ok) else rows[0,2]
    print(f"  eps < {eps_b*100:4.1f}%  ->  MD >= {md_req:.4f}   "
          f"({md_req/0.138071*100:.0f}% of the full Hall floor)")

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7.5, 5))
ax.plot(rows[:,3]*100, rows[:,2], "o-", color="#1D9E75", lw=2, ms=4)
ax.axhline(0.138071, color="#888", ls=":", lw=1, label="Hall floor (pure MD)")
ax.axvline(eps_full*100, color="#7F77DD", ls=":", lw=1, label="pure detection (12.4%)")
ax.set_xlabel("total coincidence-rate anisotropy eps (%) — a logged observable")
ax.set_ylabel("minimal measurement dependence")
ax.set_title("The guilt trade-off: rate-flatness bought with measurement dependence")
ax.legend(fontsize=9); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig("/home/claude/guilt_tradeoff.png", dpi=150)
print("saved guilt_tradeoff.png")
