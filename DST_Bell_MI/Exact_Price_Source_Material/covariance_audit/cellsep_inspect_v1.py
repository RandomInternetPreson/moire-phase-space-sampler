#!/usr/bin/env python3
"""
cellsep_inspect_v1.py -- print the shapes of the dial-blind decompositions already on disk (census_cellsep_v1/v2 jsons).

Question: the rank-2 free-shape decomposition on random0 (seed 1) split into a 76% SUM-locked shape and a 24%
DIFFERENCE-locked shape (shift and tilt).  Does the split replicate on seed 2 and on the lattice?
For each json given: per rank, per shape -- total weight, leading harmonics in (k_A, k_B) with (S,D) labels,
and the share of the shape's low-order power on the diagonal (sum-locked) vs the anti-diagonal (difference-locked).

Usage (covariance_audit):
  python cellsep_inspect_v1.py census_cellsep_v1_rho_nsh32_random0_eps0.01_nang161.json census_cellsep_v2_free_s2_rho_nsh32_random0_eps0.01_nang161.json census_cellsep_v2_g_free_s1_rho_nsh32_grid0_eps0.01_nang161.json > cellsep_inspect_10_27_26.txt
"""
import os, sys, json
import numpy as np

def spectrum(w, offs, K=4):
    xa, xb = offs[:, 0], offs[:, 1]; pix = []
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            F = np.sum(w * np.exp(2j * (ka * xa + kb * xb))) / w.sum(); pix.append((abs(F), ka, kb, np.degrees(np.angle(F))))
    pix.sort(reverse=True); return pix

for path in sys.argv[1:]:
    d = json.load(open(path)); rho = json.load(open(os.path.join(os.path.dirname(path) or ".", d["file"])))
    offs = np.asarray(rho["offsets"], float)
    print(f"\n===== {os.path.basename(path)}  (census {d['file']}; band {d.get('band', 0) or 'free'}) =====")
    for R, rec in sorted(d["ranks"].items(), key=lambda kv: int(kv[0])):
        if "F" not in rec: print(f"  rank {R}: {rec.get('status')}"); continue
        W = np.asarray(rec["W"]); F = np.asarray(rec["F"])
        print(f"  rank {R}: min miss {rec['min_miss']:.5f}" + (f", price {rec['price']:.4f}" if "price" in rec else ""))
        for r in range(F.shape[0]):
            pix = spectrum(F[r], offs)
            tot = sum(a * a for a, ka, kb, ph in pix); diag = sum(a * a for a, ka, kb, ph in pix if ka == kb); anti = sum(a * a for a, ka, kb, ph in pix if ka == -kb)
            print(f"    shape {r}: weight {W[r].sum():.3f}   power share |k|<=4: diagonal (sum-locked) {diag/tot:.2f}, anti-diagonal (difference-locked) {anti/tot:.2f}")
            print("      " + "  ".join(f"({ka:+d},{kb:+d}) {a:.2f}@{ph:+.0f}" for a, ka, kb, ph in pix[:6]))
            print("      cell weights by tau (rows) x scale (cols): " + " | ".join(" ".join(f"{W[r][ti*3+si]:.3f}" for si in range(3)) for ti in range(W.shape[1] // 3)))
