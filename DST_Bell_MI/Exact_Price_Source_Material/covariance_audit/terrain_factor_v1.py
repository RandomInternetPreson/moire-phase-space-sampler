#!/usr/bin/env python3
"""terrain_factor_v1.py -- Item 1: do the 23 beat weights w_{mn+-} factor as a product of two screens?
Usage: python terrain_factor_v1.py path/to/greedy_results.json"""
import json, sys
import numpy as np
g = json.load(open(sys.argv[1]))
terms = [tuple(t) for t in g["terms"]]; w = np.array(g["w"], float)
print(f"{len(terms)} terms, v = {g['v']:.6f}, total weight {w.sum():.4f}")
# 1. real-screen test: product of two REAL cosine screens gives w_{mn+} = w_{mn-} for every (m,n)
d = {(m, n, s): x for (m, n, s), x in zip(terms, w)}
pairs = sorted({(m, n) for m, n, s in terms if (m, n, 1) in d and (m, n, -1) in d})
print("\n1. +/- pairs present (product of two real screens predicts w+ = w-):")
print("   (m,n)   w(m+n)    w(m-n)    ratio w-/w+")
for m, n in pairs:
    print(f"   ({m},{n})  {d[(m,n,1)]:.4f}   {d[(m,n,-1)]:.4f}   {d[(m,n,-1)]/d[(m,n,1)]:8.2f}")
wp = sum(x for (m, n, s), x in zip(terms, w) if s == 1); wm = sum(x for (m, n, s), x in zip(terms, w) if s == -1)
print(f"   total weight on sum beats (k=m+n): {wp:.4f}   on difference beats (k=m-n): {wm:.4f}   share difference {wm/(wp+wm):.3f}")
# 2. rank-1 (log-linear) fit within each sign class: log w_{mn,s} = alpha_m + beta_n  (complex-screen product)
for s, name in ((1, "sum beats (+)"), (-1, "difference beats (-)")):
    idx = [i for i, t in enumerate(terms) if t[2] == s]
    ms = sorted({terms[i][0] for i in idx}); ns = sorted({terms[i][1] for i in idx})
    A = np.zeros((len(idx), len(ms) + len(ns)))
    y = np.log(w[idx])
    for r, i in enumerate(idx):
        A[r, ms.index(terms[i][0])] = 1; A[r, len(ms) + ns.index(terms[i][1])] = 1
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef; resid = y - pred
    print(f"\n2. {name}: {len(idx)} terms, m in {ms}, n in {ns}, free params {len(ms)+len(ns)-1}")
    print(f"   log-space rms residual {np.sqrt(np.mean(resid**2)):.3f}  (factor {np.exp(np.sqrt(np.mean(resid**2))):.2f}x); max |resid| {np.max(np.abs(resid)):.3f} (factor {np.exp(np.max(np.abs(resid))):.2f}x)")
    for r, i in enumerate(idx):
        print(f"      {terms[i]}  w={w[i]:.4f}  rank1={np.exp(pred[r]):.4f}  ratio={w[i]/np.exp(pred[r]):.2f}")
# 3. the ten heaviest terms
order = np.argsort(w)[::-1]
print("\n3. heaviest terms:", ", ".join(f"{terms[i]}:{w[i]:.3f}" for i in order[:10]))
