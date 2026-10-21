#!/usr/bin/env python3
"""
census_structure_v1.py  --  Item 2: what does the LP-optimal census rho* depend on?

Reads a rho_*.json produced by dither_sweep_rho.py (schema: rho, offsets, ordering=tau_major,
nsh, eps, nang, md) and asks, for the weighted empirical measure over the quenched variables
(xi_A, xi_B, tau, scale):

  (a) marginals in xi_A and xi_B      -- uniform?            (U(1) at each station)
  (b) joint structure in (xi_A, xi_B) -- Fourier plane of the weighted measure:
        power on the anti-diagonal (k,-k)  <-> rho = f(xi_A - xi_B)   [relative phase only]
        power on the      diagonal (k, k)  <-> rho depends on xi_A + xi_B [absolute phase]
        power on the axes (k,0),(0,k)      <-> product / marginal structure
        everything else                    <-> general joint dependence
  (c) weight by tau and by coupling scale
  (d) a permutation null: same weights, offsets re-paired at random (seeded), 400 draws,
      so every Fourier statistic is reported against its noise floor.

Deterministic. Prints a report; writes census_structure_v1_out.json next to the input.
Usage:  python census_structure_v1.py path/to/rho_nsh32_random0_eps0.01_nang161.json
"""
import json, sys, os
import numpy as np

TAUS   = [0.5, 1.0, 2.0, 4.0, 8.0]      # ft.TAUS   (flight_time.py)
SCALES = [1.0, 1.8, 2.6]                # ft.SCALES
KMAX   = 6
NPERM  = 400
SEED   = 20261015

def fourier_plane(xa, xb, w, kmax):
    ks = np.arange(-kmax, kmax + 1)
    F = np.zeros((len(ks), len(ks)), complex)
    for i, ka in enumerate(ks):
        for j, kb in enumerate(ks):
            F[i, j] = np.sum(w * np.exp(2j * (ka * xa + kb * xb)))
    return ks, F

def classify(ks, P):
    """split |F|^2 (excluding (0,0)) into anti-diagonal, diagonal, axes, other"""
    out = dict(anti=0.0, diag=0.0, axes=0.0, other=0.0)
    for i, ka in enumerate(ks):
        for j, kb in enumerate(ks):
            if ka == 0 and kb == 0: continue
            p = P[i, j]
            if ka == -kb:   out["anti"] += p
            elif ka == kb:  out["diag"] += p
            elif ka == 0 or kb == 0: out["axes"] += p
            else:           out["other"] += p
    return out

def main(path):
    d = json.load(open(path))
    rho = np.asarray(d["rho"], float)
    offs = np.asarray(d["offsets"], float)            # radians in [0, pi)
    nv, ntau, nsc = len(offs), len(TAUS), len(SCALES)
    assert d.get("ordering", "tau_major") == "tau_major"
    assert rho.size == nv * ntau * nsc, (rho.size, nv, ntau, nsc)
    rho = rho / rho.sum()
    R = rho.reshape(ntau, nsc, nv)                     # [tau, scale, offset]

    print(f"file: {os.path.basename(path)}")
    print(f"md = {d.get('md')}, eps = {d.get('eps')}, nang = {d.get('nang')}, nsh = {d.get('nsh')}, mode = {d.get('mode')}")
    act = rho > 1e-9
    neff = 1.0 / np.sum(rho ** 2)
    print(f"atoms total {rho.size}, active {act.sum()}, effective number 1/sum(w^2) = {neff:.1f}")

    # (c) tau and scale
    wt = R.sum(axis=(1, 2)); ws = R.sum(axis=(0, 2))
    print("\n(c) weight by tau:   " + "  ".join(f"tau={t}: {x:.4f}" for t, x in zip(TAUS, wt)))
    print("    weight by scale: " + "  ".join(f"s={s}: {x:.4f}" for s, x in zip(SCALES, ws)))
    # tau x scale table
    print("    tau x scale table (rows tau, cols scale):")
    for ti, t in enumerate(TAUS):
        print(f"      tau={t:<4}" + " ".join(f"{R[ti, si].sum():8.4f}" for si in range(nsc)))

    # weight per offset (pooled over tau, scale)
    wo = R.sum(axis=(0, 1))
    xa, xb = offs[:, 0], offs[:, 1]
    print(f"\n    offsets with any weight: {(wo > 1e-9).sum()} of {nv}; effective number over offsets {1/np.sum(wo**2):.1f}")

    # (a) marginals: circular Fourier coefficients of each marginal, period pi -> e^{2ik xi}
    print("\n(a) marginals  |<e^{2ik xi}>|  (uniform => 0; noise floor ~ 1/sqrt(n_eff))")
    print("     k   xi_A      xi_B      floor")
    for k in range(1, KMAX + 1):
        ca = abs(np.sum(wo * np.exp(2j * k * xa))); cb = abs(np.sum(wo * np.exp(2j * k * xb)))
        print(f"     {k}   {ca:.4f}    {cb:.4f}    {1/np.sqrt(1/np.sum(wo**2)):.4f}")

    # (b) Fourier plane of the joint weighted measure over offsets
    ks, F = fourier_plane(xa, xb, wo, KMAX)
    P = np.abs(F) ** 2
    cls = classify(ks, P)
    tot = sum(cls.values())
    print("\n(b) joint Fourier plane |F(k_a,k_b)|^2, |k|<=%d, (0,0) excluded" % KMAX)
    print("    share: anti-diagonal (k,-k) %.3f | diagonal (k,k) %.3f | axes %.3f | other %.3f"
          % (cls["anti"]/tot, cls["diag"]/tot, cls["axes"]/tot, cls["other"]/tot))
    # print the plane (magnitudes), rows k_a descending
    print("    |F| plane (rows k_a = +%d..-%d top to bottom, cols k_b = -%d..+%d):" % (KMAX, KMAX, KMAX, KMAX))
    for i in range(len(ks))[::-1]:
        print("      k_a=%+d " % ks[i] + " ".join(f"{abs(F[i,j]):.3f}" for j in range(len(ks))))
    # strongest pixels
    flat = [(abs(F[i, j]), ks[i], ks[j]) for i in range(len(ks)) for j in range(len(ks)) if not (ks[i] == 0 and ks[j] == 0)]
    flat.sort(reverse=True)
    print("    top 10 pixels (|F|, k_a, k_b): " + ", ".join(f"({a:.3f},{ka:+d},{kb:+d})" for a, ka, kb in flat[:10]))
    # relative-offset histogram
    rel = np.mod(xa - xb, np.pi); summ = np.mod(xa + xb, np.pi)
    nb = 12
    hr, _ = np.histogram(rel, bins=nb, range=(0, np.pi), weights=wo)
    hs, _ = np.histogram(summ, bins=nb, range=(0, np.pi), weights=wo)
    print("    weighted histogram of xi_A - xi_B (12 bins over [0,pi)): " + " ".join(f"{x:.3f}" for x in hr))
    print("    weighted histogram of xi_A + xi_B (12 bins over [0,pi)): " + " ".join(f"{x:.3f}" for x in hs))

    # (d) permutation null: same weights wo, offsets re-paired at random among all nv offsets
    rng = np.random.default_rng(SEED)
    null = {k: [] for k in ("anti", "diag", "axes", "other")}
    null_marg = []
    null_top = []
    for _ in range(NPERM):
        perm = rng.permutation(nv)
        xa_p, xb_p = xa[perm], xb[rng.permutation(nv)]   # independent re-pairing of each coordinate
        ksn, Fn = fourier_plane(xa_p, xb_p, wo, KMAX)
        Pn = np.abs(Fn) ** 2
        c = classify(ksn, Pn); t = sum(c.values())
        for kk in null: null[kk].append(c[kk] / t)
        null_marg.append(max(abs(np.sum(wo * np.exp(2j * k * xa_p))) for k in range(1, KMAX + 1)))
        null_top.append(np.sqrt(np.max(np.delete(Pn.ravel(), KMAX * (2 * KMAX + 1) + KMAX))))
    print("\n(d) permutation null (%d draws, weights fixed, offsets re-paired):" % NPERM)
    for kk in null:
        a = np.array(null[kk])
        print(f"    share {kk:<5}: observed {cls[kk]/tot:.3f}   null mean {a.mean():.3f}  sd {a.std():.3f}  "
              f"z = {(cls[kk]/tot - a.mean())/a.std():+.2f}")
    nm = np.array(null_marg); nt = np.array(null_top)
    obs_marg = max(max(abs(np.sum(wo*np.exp(2j*k*xa))), abs(np.sum(wo*np.exp(2j*k*xb)))) for k in range(1, KMAX+1))
    print(f"    max marginal |coef|: observed {obs_marg:.4f}   null 95th pct {np.percentile(nm,95):.4f}")
    print(f"    max joint pixel |F|: observed {flat[0][0]:.4f}   null 95th pct {np.percentile(nt,95):.4f}")

    # per-tau Fourier plane shares (does the structure differ with flight time?)
    print("\n    per-tau anti/diag/axes/other shares (weights within each tau renormalized):")
    for ti, t in enumerate(TAUS):
        w_t = R[ti].sum(axis=0)
        if w_t.sum() < 1e-9: print(f"      tau={t}: no weight"); continue
        w_t = w_t / w_t.sum()
        kst, Ft = fourier_plane(xa, xb, w_t, KMAX); ct = classify(kst, np.abs(Ft)**2); tt = sum(ct.values())
        print(f"      tau={t:<4} anti {ct['anti']/tt:.3f} diag {ct['diag']/tt:.3f} axes {ct['axes']/tt:.3f} other {ct['other']/tt:.3f}")

    out = dict(file=os.path.basename(path), md=d.get("md"), n_active=int(act.sum()), n_eff=float(neff),
               weight_by_tau=dict(zip(map(str, TAUS), map(float, wt))),
               weight_by_scale=dict(zip(map(str, SCALES), map(float, ws))),
               shares=dict((k, float(v / tot)) for k, v in cls.items()),
               null_shares=dict((k, [float(np.mean(v)), float(np.std(v))]) for k, v in null.items()),
               top_pixels=[(float(a), int(ka), int(kb)) for a, ka, kb in flat[:10]],
               hist_rel=hr.tolist(), hist_sum=hs.tolist())
    json.dump(out, open(os.path.join(os.path.dirname(path) or ".", "census_structure_v1_out.json"), "w"), indent=1)

if __name__ == "__main__":
    main(sys.argv[1])
