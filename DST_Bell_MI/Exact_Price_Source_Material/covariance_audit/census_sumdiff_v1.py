#!/usr/bin/env python3
"""
census_sumdiff_v1.py -- the floor-attaining census written as a density in (S, D) = (xi_A + xi_B, xi_A - xi_B),
and the conditional structure of D given S (the "anisotropy of the push").

For each rho_*.json given (dither_sweep_rho.py schema), pooled over tau and scale:
  (1) the truncated Fourier series  rho(S, D) ~ 1 + sum_k 2|F_k| cos(k_S S + k_D D - phi_k)   for |k_S|,|k_D| <= KMAX
      (both S and D have period pi on the torus; harmonics are e^{i(2 k_S S + 2 k_D D)} in the paper's convention,
      written here with k_S, k_D counting multiples of 2),
      listed by magnitude, with the permutation null (offsets re-paired, weights fixed) for scale;
  (2) the S marginal: preferred S (mod 90 deg for the k=2 harmonic) and concentration;
  (3) conditional slices: S binned into NS bins; in each, the weighted D-distribution's circular mean and
      concentration at harmonics cos 2D and cos 4D  ->  mu(S) and kappa(S), the push's size and spread vs the axis;
  (4) a product-vs-joint diagnostic: |F(k_S,k_D)| vs |F(k_S,0)||F(0,k_D)| for every cross pixel.
If several files are given, (1)-(3) are also reported for the pooled census (all files, equal weight).

Usage:  python census_sumdiff_v1.py rho_nsh32_random0_eps0.01_nang161.json [more rho files...]
Env:    KMAX (default 4), NS (default 8), NPERM (default 400), SEED
Writes: census_sumdiff_v1_out.json next to the first file.  Seconds to run; numpy only.
"""
import os, sys, json
import numpy as np

KMAX = int(os.environ.get("KMAX", 4)); NS = int(os.environ.get("NS", 8))
NPERM = int(os.environ.get("NPERM", 400)); SEED = int(os.environ.get("SEED", 20261021))
TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]; SCALES = [1.0, 1.8, 2.6]

def load(path):
    d = json.load(open(path))
    rho = np.asarray(d["rho"], float); offs = np.asarray(d["offsets"], float); NV = len(offs)
    nflat = rho.size // NV
    wo = rho.reshape(nflat, NV).sum(axis=0); wo = wo / wo.sum()
    S = (offs[:, 0] + offs[:, 1]) % np.pi; D = (offs[:, 0] - offs[:, 1]) % np.pi
    return d, wo, S, D

def plane(S, D, w, kmax):
    ks = np.arange(-kmax, kmax + 1)
    F = np.zeros((len(ks), len(ks)), complex)
    for i, a in enumerate(ks):
        for j, b in enumerate(ks):
            F[i, j] = np.sum(w * np.exp(2j * (a * S + b * D)))
    return ks, F

def circ(x, w, k):
    """weighted circular mean and resultant length of k-th harmonic of x (period pi): returns (mean_deg, R)"""
    z = np.sum(w * np.exp(2j * k * x)) / w.sum()
    return np.degrees(np.angle(z) / (2 * k)) % (180.0 / k), abs(z)

def report(name, S, D, w, rng):
    print(f"\n===== {name} =====")
    neff = 1 / np.sum(w ** 2)
    print(f"effective number of offsets {neff:.1f}")
    ks, F = plane(S, D, w, KMAX)
    # null
    null = []
    for _ in range(NPERM):
        p = rng.permutation(len(w))
        _, Fn = plane(S, D, w[p], KMAX)   # weights re-assigned to other offsets: breaks weight-position structure
        Fn[KMAX, KMAX] = 0; null.append(np.abs(Fn).max())
    n95 = np.percentile(null, 95)
    print(f"permutation null for any single pixel |F|: 95th pct {n95:.3f}")
    # (1) series
    items = []
    for i, a in enumerate(ks):
        for j, b in enumerate(ks):
            if a < 0 or (a == 0 and b <= 0): continue   # one of each conjugate pair
            items.append((abs(F[i, j]), a, b, np.degrees(np.angle(F[i, j]))))
    items.sort(reverse=True)
    print("\n(1) rho(S,D) ~ 1 + sum 2|F| cos(2 k_S S + 2 k_D D - phi)   [pixels above the null]")
    print("     k_S  k_D    |F|     phi[deg]   type")
    out_terms = []
    for amp, a, b, ph in items:
        if amp < n95: continue
        typ = "pure S" if b == 0 else ("pure D" if a == 0 else "cross (S-D coupling)")
        print(f"    {a:+3d}  {b:+3d}   {amp:.3f}   {ph:+8.2f}   {typ}")
        out_terms.append(dict(kS=int(a), kD=int(b), amp=float(amp), phi_deg=float(ph), type=typ))
    # (2) S marginal
    print("\n(2) S marginal (axis of the birth, pinned frame):")
    for k in (1, 2, 3):
        m, R = circ(S, w, k)
        print(f"    harmonic cos {2*k}S: preferred S = {m:6.2f} deg (mod {180/k:.0f}), concentration {R:.3f}")
    # (3) conditional slices
    print(f"\n(3) D given S  ({NS} bins of S):")
    print("     S-bin [deg]      mass    <D> at cos2D [deg]  R2     <D> at cos4D [deg]  R4")
    edges = np.linspace(0, np.pi, NS + 1); slices = []
    for s0, s1 in zip(edges[:-1], edges[1:]):
        m = (S >= s0) & (S < s1); mass = w[m].sum()
        if mass < 1e-6:
            print(f"    [{np.degrees(s0):5.1f},{np.degrees(s1):5.1f})   {mass:.3f}   (empty)"); continue
        mu2, R2 = circ(D[m], w[m], 1); mu4, R4 = circ(D[m], w[m], 2)
        print(f"    [{np.degrees(s0):5.1f},{np.degrees(s1):5.1f})   {mass:.3f}      {mu2:6.1f}        {R2:.3f}      {mu4:6.1f}        {R4:.3f}")
        slices.append(dict(S_lo=float(np.degrees(s0)), S_hi=float(np.degrees(s1)), mass=float(mass), muD2=float(mu2), R2=float(R2), muD4=float(mu4), R4=float(R4)))
    # (4) product diagnostic
    print("\n(4) cross pixels vs product of marginals  |F(kS,kD)|  vs  |F(kS,0)||F(0,kD)|:")
    for amp, a, b, ph in items:
        if a == 0 or b == 0 or amp < n95: continue
        pa = abs(F[list(ks).index(a), KMAX]) * abs(F[KMAX, list(ks).index(b)])
        print(f"    ({a:+d},{b:+d}): {amp:.3f} vs {pa:.3f}   ratio {amp/max(pa,1e-9):.1f}")
    return dict(name=name, n_eff=float(neff), null95=float(n95), terms=out_terms, slices=slices)

def main(paths):
    rng = np.random.default_rng(SEED)
    results = []; pooled = []
    for p in paths:
        d, w, S, D = load(p)
        print(f"\nfile {os.path.basename(p)}: md = {d.get('md')}, nsh {d.get('nsh')}, mode {d.get('mode')}{d.get('seed')}")
        results.append(report(os.path.basename(p), S, D, w, rng)); pooled.append((w, S, D))
    if len(paths) > 1:
        w = np.concatenate([x[0] for x in pooled]) / len(pooled)
        S = np.concatenate([x[1] for x in pooled]); D = np.concatenate([x[2] for x in pooled])
        results.append(report(f"POOLED ({len(paths)} censuses, equal weight)", S, D, w, rng))
    json.dump(results, open(os.path.join(os.path.dirname(paths[0]) or ".", "census_sumdiff_v1_out.json"), "w"), indent=1)

if __name__ == "__main__":
    main(sys.argv[1:])
