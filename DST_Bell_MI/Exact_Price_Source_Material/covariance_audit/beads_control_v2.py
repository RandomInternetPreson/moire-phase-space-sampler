#!/usr/bin/env python3
"""
beads_control_v2.py -- v1 + the SPLIT of the partition: between cells of different tau at the same scale, vs different scale at the same tau.

--- v1 docstring follows ---
beads_control_v1.py -- is the partition of difference values across cells physical, or the solver's thrift?

Input: several census_bumps_v2 jsons from the SAME class with different PSEED (different optimal vertices).
  (1) per run: the between-cell D-overlap and its shuffle-null z (as beads_inspect_v1)
  (2) stability: for each cell, the overlap of its D-histogram between runs (same cell, different vertices), averaged;
      compared with the between-cell overlap within a run.  Physical partition -> same-cell/different-run overlap HIGH,
      between-cell LOW.  Thrift -> same-cell/different-run overlap as low as between-cell (the assignment reshuffles).
  (3) the averaged census (mean of the runs' bead masses per (cell, position)): its between-cell D-overlap and z.
      Physical -> still below the null.  Thrift -> at the null.
Usage: python beads_control_v1.py census_bumps_v2_ctl_p1_*.json census_bumps_v2_ctl_p2_*.json ... > beads_control_<date>.txt
"""
import sys, json, os
import numpy as np
NB = 32; EDGES = np.linspace(0, 180, NB + 1)
def hist(vals, w):
    h, _ = np.histogram(np.asarray(vals) % 180, bins=EDGES, weights=w); return h / max(h.sum(), 1e-300)
runs = []
for path in sys.argv[1:]:
    d = json.load(open(path)); rows = d["bumps"]
    runs.append(dict(name=os.path.basename(path), price=d["price"], tau=np.array([r[0] for r in rows]), sc=np.array([r[1] for r in rows]),
                     a=np.array([r[2] for r in rows]), b=np.array([r[3] for r in rows]), w=np.array([r[4] for r in rows])))
cells = sorted(set((t, s) for r in runs for t, s in zip(r["tau"], r["sc"])))
def cell_hists(r):
    H = {}
    for c in cells:
        m = (r["tau"] == c[0]) & (r["sc"] == c[1])
        if r["w"][m].sum() > 1e-4: H[c] = hist(r["a"][m] - r["b"][m], r["w"][m])
    return H
def between(H):
    k = list(H); return np.mean([np.minimum(H[k[i]], H[k[j]]).sum() for i in range(len(k)) for j in range(i + 1, len(k))]) if len(k) > 1 else np.nan
def null_z(r, obs, ndraw=200, seed=1):
    rng = np.random.default_rng(seed); D = (r["a"] - r["b"]); vals = []
    for _ in range(ndraw):
        Dp = rng.permutation(D); H = {}
        for c in cells:
            m = (r["tau"] == c[0]) & (r["sc"] == c[1])
            if r["w"][m].sum() > 1e-4: H[c] = hist(Dp[m], r["w"][m])
        vals.append(between(H))
    return np.mean(vals), np.std(vals), (obs - np.mean(vals)) / max(np.std(vals), 1e-9)
print(f"{len(runs)} runs; {len(cells)} cells")
Hs = []
print("(1) per run: price, between-cell D-overlap, shuffle null, z")
for r in runs:
    H = cell_hists(r); Hs.append(H); ob = between(H); m, s, z = null_z(r, ob)
    print(f"    {r['name'][:40]:40s} price {r['price']:.4f}   overlap {ob:.3f}   null {m:.3f}+-{s:.3f}   z = {z:+.1f}")
print("(2) stability of each cell's D-set across runs (same cell, different vertices): overlap averaged over run pairs")
same = []
for c in cells:
    hs = [H[c] for H in Hs if c in H]
    if len(hs) > 1:
        v = np.mean([np.minimum(hs[i], hs[j]).sum() for i in range(len(hs)) for j in range(i + 1, len(hs))]); same.append(v)
        print(f"    tau={c[0]:<4} s={c[1]:<4} same-cell overlap across runs {v:.3f}")
print(f"    MEAN same-cell/different-run overlap {np.mean(same):.3f}   vs   mean between-cell/same-run overlap {np.mean([between(H) for H in Hs]):.3f}")
print("    reading: same-cell >> between-cell -> the cells keep their D-sets (physical); same-cell ~ between-cell -> the assignment reshuffles (thrift)")
# (3) averaged census
print("(3) averaged census over runs (bead masses averaged per (cell, position)):")
acc = {}
for r in runs:
    for t, s, a, b, w in zip(r["tau"], r["sc"], r["a"], r["b"], r["w"]):
        key = (t, s, round(a, 2), round(b, 2)); acc[key] = acc.get(key, 0.0) + w / len(runs)
avg = dict(tau=np.array([k[0] for k in acc]), sc=np.array([k[1] for k in acc]), a=np.array([k[2] for k in acc]), b=np.array([k[3] for k in acc]), w=np.array(list(acc.values())))
H = cell_hists(avg); ob = between(H); m, s, z = null_z(avg, ob)
print(f"    beads in the average {len(acc)}; between-cell D-overlap {ob:.3f}   null {m:.3f}+-{s:.3f}   z = {z:+.1f}")
# (4) split of the between-cell overlap in the averaged census: same scale / different tau  vs  same tau / different scale
def split_overlap(H):
    k = list(H); st, ss, oth = [], [], []
    for i in range(len(k)):
        for j in range(i + 1, len(k)):
            v = np.minimum(H[k[i]], H[k[j]]).sum()
            if k[i][1] == k[j][1]: st.append(v)        # same scale, different tau
            elif k[i][0] == k[j][0]: ss.append(v)      # same tau, different scale
            else: oth.append(v)
    return (np.mean(st) if st else np.nan, len(st)), (np.mean(ss) if ss else np.nan, len(ss)), (np.mean(oth) if oth else np.nan, len(oth))
obs = split_overlap(H)
rng = np.random.default_rng(2); D = avg["a"] - avg["b"]; nulls = {0: [], 1: [], 2: []}
for _ in range(200):
    Dp = rng.permutation(D); Hp = {}
    for c in cells:
        mm = (avg["tau"] == c[0]) & (avg["sc"] == c[1])
        if avg["w"][mm].sum() > 1e-4: Hp[c] = hist(Dp[mm], avg["w"][mm])
    sp = split_overlap(Hp)
    for q in range(3): nulls[q].append(sp[q][0])
print("(4) the split (averaged census): which pairs of cells avoid each other's differences?")
for q, name in enumerate(("same scale, different tau   (the difference phase sorts the lock DURATION)", "same tau, different scale   (the difference phase sorts the DEPTH)", "different tau and scale")):
    o, n = obs[q]; mu, sd = np.nanmean(nulls[q]), np.nanstd(nulls[q])
    print(f"    {name:70s} pairs {n:3d}  overlap {o:.3f}   null {mu:.3f}+-{sd:.3f}   z = {(o-mu)/max(sd,1e-9):+.1f}")
print("    reading: z still strongly negative -> partition survives averaging (physical); z ~ 0 -> it was the vertex (thrift)")
