#!/usr/bin/env python3
"""
rim_skeleton_v1.py -- the Fourier skeleton of the RIM: the lean field of a joint-birth solution saved by census_fuzz_v2/v3.

The solution npz holds u[l, atom]: the joint density's lambda_0 harmonics per atom, w(atom, lambda_0) = u_0 + u_1 cos 2 lambda_0 +
u_2 sin 2 lambda_0 (+ u_3 cos 4 lambda_0 + u_4 sin 4 lambda_0 ...).  The lean field is z = u_1 + i u_2 = |z| e^{2 i psi}: psi(xi_A, xi_B)
is the preferred birth phase at that atom, |z| its strength.  This script prints: the modulation depth; the fraction of atoms at the
positivity limit (a half-wave rim); the Fourier plane |F(k_A, k_B)| of z pooled over the cells (|k| <= 8), with a within-cell shuffle
null; the parity rule (the +90 twin symmetry puts the census w_0 on k_A + k_B EVEN pixels and the lean on ODD ones); the top pixel per
cell; and, for L >= 2, the same for the second harmonic z_2 = u_3 + i u_4 (even pixels).
Usage:  python rim_skeleton_v1.py census_fuzz_v2_L1_K6_atom_solution.npz > rim_skeleton_L1_K6_atom.txt
"""
import sys, numpy as np
f = sys.argv[1]; z = np.load(f); u = z["u"]; offs = z["offsets"]; L = int(z["L"]); K = int(z["K"]); scope = str(z["scope"])
NV = len(offs); nflat = u.shape[1] // NV; NT = nflat * NV; KK = 8
print(f"rim_skeleton_v1 on {f}: L={L} K={K} scope={scope}; {nflat} cells x {NV} offsets")
w0 = u[0]; act = w0 > 1e-7; xa = np.tile(offs[:, 0], nflat); xb = np.tile(offs[:, 1], nflat)
taus = [0.5, 1, 2, 4, 8]; scs = [1.0, 1.8, 2.6]
print(f"atoms with mass {act.sum()}/{NT}; mass {w0.sum():.6f}")

def plane(fld, norm):
    P = np.zeros((2 * KK + 1, 2 * KK + 1), complex)
    for i, ka in enumerate(range(-KK, KK + 1)):
        for j, kb in enumerate(range(-KK, KK + 1)):
            P[i, j] = np.sum(fld * np.exp(-2j * (ka * xa + kb * xb))) / norm
    return P

def report(fld, name, parity_expected):
    amp = np.abs(fld); r = amp / np.maximum(w0, 1e-15); is_rim = not np.allclose(fld.imag, 0)
    print(f"\n=== {name} ===")
    if is_rim:
        print(f"modulation depth (mass-weighted sum|z| / sum w0): {amp.sum()/w0.sum():.3f}")
        print(f"relative modulation |z|/w0 at active atoms: percentiles 10/50/90 = {np.percentile(r[act], 10):.3f} / {np.percentile(r[act], 50):.3f} / {np.percentile(r[act], 90):.3f};  fraction > 0.95: {(r[act] > 0.95).mean():.3f}")
    P = plane(fld, amp.sum()); A = np.abs(P)
    par = np.add.outer(np.arange(-KK, KK + 1), np.arange(-KK, KK + 1)) % 2
    share_odd = (A[par == 1] ** 2).sum() / (A ** 2).sum()
    print(f"share of the spectrum on odd-parity pixels (k_A + k_B odd): {share_odd:.3f}   (expected {'odd' if parity_expected == 1 else 'even'} by the +90 twin symmetry)")
    order = np.argsort(-A.ravel()); shown = 0
    print("top pixels  |F|  (k_A, k_B)  phase[deg]:")
    for idx in order:
        i, j = divmod(idx, 2 * KK + 1); ka, kb = i - KK, j - KK
        if ka == 0 and kb == 0: continue
        print(f"  {A[i, j]:.3f}  ({ka:+d},{kb:+d})  {np.degrees(np.angle(P[i, j])):+.0f}"); shown += 1
        if shown >= 10: break
    rng = np.random.default_rng(1); nulls = []
    for _ in range(200):
        fs = fld.copy()
        for ts in range(nflat):
            sl = slice(ts * NV, (ts + 1) * NV); fs[sl] = rng.permutation(fld[sl])
        An = np.abs(plane(fs, amp.sum())); An[KK, KK] = 0; nulls.append(An.max())
    print(f"shuffle null (values permuted among the offsets within each cell, 200 draws): max pixel 95th pct {np.percentile(nulls, 95):.3f}, mean {np.mean(nulls):.3f}")
    print(f"|F| on the {'odd' if parity_expected == 1 else 'even'} pixels, |k| <= 4 (rows k_A = +4..-4, cols k_B = -4..+4; '.' = other parity):")
    for ka in range(4, -5, -1):
        print(f"  k_A={ka:+d} " + " ".join(("  .  " if (ka + kb) % 2 != parity_expected else f"{A[ka + KK, kb + KK]:.3f}") for kb in range(-4, 5)))
    print("per cell (tau, s): mass, depth, top pixel:")
    for ts in range(nflat):
        sl = slice(ts * NV, (ts + 1) * NV); m = w0[sl].sum()
        if m < 1e-6: continue
        best = None
        for ka in range(-KK, KK + 1):
            for kb in range(-KK, KK + 1):
                if ka == 0 and kb == 0: continue
                v = abs(np.sum(fld[sl] * np.exp(-2j * (ka * offs[:, 0] + kb * offs[:, 1])))) / max(amp[sl].sum(), 1e-300)
                if best is None or v > best[0]: best = (v, ka, kb)
        print(f"  tau={taus[ts // 3]:<4} s={scs[ts % 3]:<4} mass {m:.4f}" + (f" depth {amp[sl].sum()/m:.3f}" if is_rim else "") + f"  top |F| {best[0]:.3f} at ({best[1]:+d},{best[2]:+d})")

report(w0.astype(complex), "census w_0 (the K-band-limited census; even pixels)", 0)
report(u[1] + 1j * u[2], "lean field z_1 = u_1 + i u_2 (first harmonic of lambda_0; odd pixels)", 1)
if L >= 2: report(u[3] + 1j * u[4], "second harmonic z_2 = u_3 + i u_4 (cos/sin 4 lambda_0; even pixels)", 0)
