#!/usr/bin/env python3
"""
beads_inspect_v1.py -- the beads on the sum-lock ridge: are their DIFFERENCE values a fixed set, or scattered afresh per cell?

Input: census_bumps_v1 jsons (the 32x32 runs, whose 5.6-degree grid does not pre-quantize D as coarsely as 16x16).
For each json: the active bumps (position, mass, cell); S = xi_A + xi_B and D = xi_A - xi_B mod 180.
  (1) S histogram (mod 90) -- the ridge; (2) D histogram, pooled, 5.625-degree bins -- a fixed set shows as repeated peaks;
  (3) D histogram per scale and per tau; (4) the pairwise agreement of the D-sets between cells (do the same D values recur?):
      for each pair of cells, the mass-weighted overlap of their D histograms (1 = identical sets, ~bin-filling fraction = random);
  (5) the D values of the top beads per cell, so the sets can be read by eye.
Usage (covariance_audit): python beads_inspect_v1.py census_bumps_v1_cr_k6_w6_nc32_rho_nsh32_grid0_eps0.01_nang161.json census_bumps_v1_cr_k6_w3_nc32_rho_nsh32_grid0_eps0.01_nang161.json > beads_inspect_10_29_26.txt
"""
import sys, json, os
import numpy as np
TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]; SCALES = [1.0, 1.8, 2.6]; NB = 32; EDGES = np.linspace(0, 180, NB + 1)

def hist(vals, w):
    h, _ = np.histogram(np.asarray(vals) % 180, bins=EDGES, weights=w); return h / max(h.sum(), 1e-300)
def bar(h, width=48):
    m = h.max() if h.max() > 0 else 1; return "".join("#" if x > 0.5 * m else ("+" if x > 0.2 * m else ("." if x > 0 else " ")) for x in h)

for path in sys.argv[1:]:
    d = json.load(open(path)); rows = d["bumps"]   # (tau, scale, xiA_deg, xiB_deg, mass)
    tau = np.array([r[0] for r in rows]); sc = np.array([r[1] for r in rows]); a = np.array([r[2] for r in rows]); b = np.array([r[3] for r in rows]); w = np.array([r[4] for r in rows])
    S = (a + b) % 180; D = (a - b) % 180
    print(f"\n===== {os.path.basename(path)}: K={d['K']} + bumps w={d['bump_w']} nc={d['bump_nc']}; price {d['price']:.4f}; {len(rows)} active beads, mass {w.sum():.3f} =====")
    hS = hist(S % 90 * 2, w)   # S mod 90 stretched to [0,180) for the same bins: bin width 2.8 deg in S
    print("(1) S mod 90 (bins of 2.8 deg, 0..90):     " + bar(hS) + f"   peak at S = {EDGES[np.argmax(hS)]/2:.1f} deg, top-3 bins carry {np.sort(hS)[-3:].sum():.2f} of the mass")
    hD = hist(D, w)
    print("(2) D pooled (bins of 5.6 deg, 0..180):     " + bar(hD) + f"   top-5 bins carry {np.sort(hD)[-5:].sum():.2f} of the mass; bins occupied {int((hD>0).sum())}/{NB}")
    print("    D bin centres with > 5% of the mass: " + ", ".join(f"{(EDGES[i]+EDGES[i+1])/2:.0f}deg {hD[i]:.2f}" for i in np.argsort(hD)[::-1] if hD[i] > 0.05))
    print("(3) D per scale and per tau:")
    for s in SCALES:
        m = sc == s
        if w[m].sum() > 0: print(f"    s={s:<4}  " + bar(hist(D[m], w[m])) + f"   mass {w[m].sum():.3f}")
    for t in TAUS:
        m = tau == t
        if w[m].sum() > 0: print(f"    tau={t:<4}" + bar(hist(D[m], w[m])) + f"   mass {w[m].sum():.3f}")
    # (4) overlap of D-sets between cells
    cells = sorted(set(zip(tau, sc))); H = {}
    for c in cells:
        m = (tau == c[0]) & (sc == c[1])
        if w[m].sum() > 1e-4: H[c] = hist(D[m], w[m])
    ov = []; fill = []
    keys = list(H)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            ov.append(np.minimum(H[keys[i]], H[keys[j]]).sum())
    for k in keys: fill.append((H[k] > 0).mean())
    # null: shuffle D within the pooled set, re-assign to cells, recompute the mean overlap (200 draws)
    rng = np.random.default_rng(1); null = []
    for _ in range(200):
        Dp = rng.permutation(D); Hp = {}
        for c in keys:
            m = (tau == c[0]) & (sc == c[1]); Hp[c] = hist(Dp[m], w[m])
        null.append(np.mean([np.minimum(Hp[keys[i]], Hp[keys[j]]).sum() for i in range(len(keys)) for j in range(i + 1, len(keys))]))
    print(f"(4) D-set overlap between cells (histogram intersection; 1 = same set): mean {np.mean(ov):.2f}, max {np.max(ov):.2f};   "
          f"null (D values shuffled across cells): mean {np.mean(null):.2f} +- {np.std(null):.2f}   -> z = {(np.mean(ov)-np.mean(null))/max(np.std(null),1e-9):+.1f}")
    print("(5) top beads per cell, as (S, D) in deg with mass:")
    for c in cells:
        m = (tau == c[0]) & (sc == c[1]); idx = np.argsort(w[m])[::-1][:6]
        print(f"    tau={c[0]:<4} s={c[1]:<4} n={int(m.sum()):3d}  " + "  ".join(f"(S{S[m][i]:5.1f}, D{D[m][i]:6.1f}) {w[m][i]:.4f}" for i in idx))
