#!/usr/bin/env python3
"""
census_bumps_v1.py -- the CREST test: smooth census + a dictionary of localized bumps.  How many crests, and do they travel?

Aaron's idea (29 Oct): the fine structure of the birth is a single localized event -- a crest -- that travels around the
disc's hem (or splits and travels both ways).  A localized bump of angular width w has a FLAT spectrum up to k ~ 180/w
and then falls: a contiguous band, which is what the floor census needs (7..14, no holes).  A traveling crest would appear
in the census as one or two bumps per (tau, scale) cell whose positions advance with tau.

The LP: rho(cell, xi) = [trig polynomial of degree K, free coefficients, itself >= 0 at every atom] + sum_i b_i B_i(xi),  b_i >= 0,
with B_i a von Mises bump of half-width BUMP_W degrees centred on a BUMP_NC x BUMP_NC grid of (xi_A, xi_B);
positivity at every atom; the paper's correlator / CHSH / MD rows (min_md_full).  One LP; exact optimum.
Reported: price; the number of bumps carrying weight per cell and in all; the share of mass in bumps vs the smooth part;
the top bumps per cell (position, weight); and the DRIFT: per tau, the weighted circular mean of the active bumps'
positions in xi_A, xi_B and S -- a crest that travels shows a mean that advances with tau.
Caveat: LP vertices are sparse by construction, so the active count is an UPPER bound on the number of crests needed; what
is decisive is whether it is ~1-2 per cell or dozens, and whether the positions march.

Env: K (default 6; 0 = bumps only), BUMP_W (deg, default 6), BUMP_NC (default 16), EPS (0.01), METHOD (highs-ipm), OUT_TAG.
Usage (covariance_audit):  python census_bumps_v1.py rho_nsh32_grid0_eps0.01_nang161.json > bumps_<tag>.txt
Compare (lattice): K=6 alone 1.700; K=6 + one ripple 1.665-1.700; K=14 1.000; the unrestricted census has 699 active atoms (47/cell).
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
import census_bandlimit_v4 as v4            # load, md_of, reduce, spectrum, verify helpers (env-free parts)
from census_bandlimit_v4 import FLOOR, S_QM, RED, _deg

K = int(os.environ.get("K", 6)); BUMP_W = float(os.environ.get("BUMP_W", 6)); BUMP_NC = int(os.environ.get("BUMP_NC", 16))
EPS = float(os.environ.get("EPS", 0.01)); METHOD = os.environ.get("METHOD", "highs-ipm"); OUT_TAG = os.environ.get("OUT_TAG", f"k{K}_w{BUMP_W:g}_nc{BUMP_NC}")
TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]; SCALES = [1.0, 1.8, 2.6]

def trig_basis(offs, K):
    xa, xb = offs[:, 0], offs[:, 1]; cols = [np.ones(len(offs))]
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            arg = 2 * (ka * xa + kb * xb); cols.append(np.cos(arg)); cols.append(np.sin(arg))
    return np.column_stack(cols)

def bump_basis(offs, w_deg, nc):
    kappa = np.log(2) / (1 - np.cos(2 * w_deg * _deg))            # half-maximum at +-w in each coordinate
    cs = (np.arange(nc) + 0.5) * np.pi / nc; centers = [(ca, cb) for ca in cs for cb in cs]
    xa, xb = offs[:, 0], offs[:, 1]
    B = np.column_stack([np.exp(kappa * ((np.cos(2 * (xa - ca)) - 1) + (np.cos(2 * (xb - cb)) - 1))) for ca, cb in centers])
    B /= B.max(axis=0, keepdims=True)                              # peak 1 (unit-mass scaling is badly conditioned against the O(1) trig columns)
    return B, np.array(centers), kappa, B.sum(axis=0)              # colmass: mass of each bump at the atoms

def price_lp(Qr, Sr, Er, Mr, P, angs, eps, nflat, nb, ntrig, method):
    """v4's price LP with mixed bounds: trig coefficients free, bump coefficients >= 0; positivity rows on rho."""
    nc = nflat * nb; nang = len(angs)
    cmb = [(i, j) for i in range(4) for j in range(i + 1, 4)]; n_s = len(cmb) * RED; nv = nc + 1 + n_s
    Qc = [Q.reshape(nc, RED, -1).mean(axis=2) for Q in Qr]; binw = np.pi / RED; blocks, bub = [], []
    for ci, (i, j) in enumerate(cmb):
        Dif = (Qc[i] - Qc[j]) * binw
        for sign in (1, -1):
            A = sp.lil_matrix((RED, nv)); A[:, :nc] = sign * Dif.T
            for l in range(RED): A[l, nc + 1 + ci * RED + l] = -1
            blocks.append(A.tocsr()); bub += [0.0] * RED
    T = sp.lil_matrix((len(cmb), nv))
    for ci in range(len(cmb)): T[ci, nc + 1 + ci * RED:nc + 1 + (ci + 1) * RED] = 0.5; T[ci, nc] = -1
    blocks.append(T.tocsr()); bub += [0.0] * len(cmb)
    tgt = -np.cos(2 * angs * _deg)
    Ep = sp.hstack([sp.csr_matrix(Er), sp.csr_matrix((nang, 1 + n_s))]).tocsr(); blocks += [Ep, -Ep]; bub += list(tgt + eps) + list(-tgt + eps)
    blocks.append(sp.hstack([P, sp.csr_matrix((P.shape[0], 1 + n_s))]).tocsr()); bub += [0.0] * P.shape[0]
    Aub = sp.vstack(blocks).tocsr()
    Aeq = sp.csr_matrix(np.vstack([np.concatenate([Mr, [0], np.zeros(n_s)]), np.concatenate([Sr, [0], np.zeros(n_s)])]))
    cobj = np.zeros(nv); cobj[nc] = 1
    bounds = []
    for ts in range(nflat): bounds += [(None, None)] * ntrig + [(0, None)] * (nb - ntrig)
    bounds += [(0, None)] * (1 + n_s)
    t0 = time.time(); r = linprog(cobj, A_ub=Aub, b_ub=np.array(bub), A_eq=Aeq, b_eq=[1.0, -S_QM], bounds=bounds, method=method)
    print(f"    [{method}] status {r.status}; {r.message.strip()}   [{time.time()-t0:.0f}s]", flush=True)
    return r

def main(path):
    t0 = time.time()
    d, offs, nflat, Qs, Svec, Ecoef, angs = v4.load(path); NV = len(offs); NT = nflat * NV
    T = trig_basis(offs, K) if K > 0 else np.ones((NV, 1)); B, centers, kappa, colmass = bump_basis(offs, BUMP_W, BUMP_NC)
    Phi = np.hstack([T, B]); ntrig = T.shape[1]; nb = Phi.shape[1]
    print(f"file: {os.path.basename(path)}   atoms {NT} = {nflat} x {NV}")
    print(f"class: smooth K={K} ({ntrig} functions/cell, free) + {B.shape[1]} bumps/cell of half-width {BUMP_W:g} deg (kappa {kappa:.1f}) on a {BUMP_NC}x{BUMP_NC} grid, coefficients >= 0;  eps {EPS}", flush=True)
    Qr, Sr, Er, Mr, _ = v4.reduce(Qs, Svec, Ecoef, Phi, nflat, NV)
    # positivity on the SMOOTH PART ALONE (-T c <= 0 at every atom); bumps are >= 0 by their bounds, so rho >= 0 follows.
    # (positivity on the sum lets free trig coefficients and bumps cancel each other -- a degenerate LP that stalls both solvers)
    P = sp.block_diag([sp.hstack([sp.csr_matrix(-T), sp.csr_matrix((NV, nb - ntrig))])] * nflat).tocsr()
    print("  price LP (smooth part >= 0, crests >= 0):", flush=True)
    r = price_lp(Qr, Sr, Er, Mr, P, angs, EPS, nflat, nb, ntrig, METHOD)
    if r.status != 0: sys.exit(f"LP not optimal (status {r.status})")
    c = r.x[:nflat * nb].reshape(nflat, nb); rho, v = v4.verify(c, Phi, nflat, Qs, Er, angs, EPS)
    price = v["price"]
    bw = c[:, ntrig:] * colmass[None, :]                            # bump masses per cell
    mass_bumps = bw.sum(); mass_smooth = 1 - mass_bumps
    thr = 1e-6 * max(bw.max(), 1e-30); active = (bw > thr)
    print(f"\nRESULT  K={K} + bumps w={BUMP_W:g} deg: MD/floor = {price:.4f}   (K=6 alone 1.700; K=14 1.000)")
    print(f"  mass in bumps {mass_bumps:.3f}, in the smooth part {mass_smooth:.3f};  active bumps: {active.sum()} in all, per cell " + " ".join(str(int(a)) for a in active.sum(axis=1)) + f"   (atoms: 699 active = 47/cell)")
    print("\n  top bumps per cell (xi_A, xi_B in deg; weight):")
    rows = []
    for cidx in range(nflat):
        idx = np.argsort(bw[cidx])[::-1][:4]; tau, sc = TAUS[cidx // 3], SCALES[cidx % 3]
        s = "  ".join(f"({np.degrees(centers[i][0]):5.1f},{np.degrees(centers[i][1]):5.1f}) {bw[cidx, i]:.4f}" for i in idx if bw[cidx, i] > thr)
        print(f"    tau={tau:<4} s={sc:<4} n={int(active[cidx].sum()):3d}  {s}")
        for i in np.where(active[cidx])[0]: rows.append((tau, sc, float(np.degrees(centers[i][0])), float(np.degrees(centers[i][1])), float(bw[cidx, i])))
    # drift with tau: weighted circular mean of the active bumps' positions (period 180 deg), pooled over scale
    print("\n  drift: weighted circular mean of active bump positions vs tau (pooled over scale) -- a traveling crest advances with tau:")
    print("     tau     <xi_A>   R_A    <xi_B>   R_B    <S=xi_A+xi_B>  R_S    mass")
    drift = []
    for tau in TAUS:
        sel = [(a, b, w) for (t, sc, a, b, w) in rows if t == tau]
        if not sel: print(f"    {tau:<6} (no bumps)"); continue
        a = np.radians([x[0] for x in sel]); b = np.radians([x[1] for x in sel]); w = np.array([x[2] for x in sel])
        def cm(x): z = np.sum(w * np.exp(2j * x)) / w.sum(); return (np.degrees(np.angle(z)) / 2) % 180, abs(z)
        ma, ra = cm(a); mb, rb = cm(b); ms, rs = cm(a + b)
        print(f"    {tau:<6} {ma:7.1f}  {ra:.2f}  {mb:7.1f}  {rb:.2f}    {ms:7.1f}     {rs:.2f}   {w.sum():.3f}")
        drift.append(dict(tau=tau, xiA=ma, RA=ra, xiB=mb, RB=rb, S=ms, RS=rs, mass=float(w.sum())))
    pix = v4.spectrum(rho, offs, nflat, NV, 14)
    print("\n  leading harmonics of the pooled census (|F|, k_A, k_B, phi):")
    for a, ka, kb, ph in pix[:10]: print(f"     {a:.3f}  ({ka:+d},{kb:+d})  {ph:+8.2f}")
    out = dict(file=os.path.basename(path), K=K, bump_w=BUMP_W, bump_nc=BUMP_NC, eps=EPS, price=float(price), mass_bumps=float(mass_bumps),
               active_total=int(active.sum()), active_per_cell=[int(a) for a in active.sum(axis=1)], bumps=rows, drift=drift,
               leading=[(float(a), int(ka), int(kb), float(ph)) for a, ka, kb, ph in pix[:20]], rho=rho.tolist(), offsets=d["offsets"], md=float(price * FLOOR),
               nsh=d["nsh"], mode=d["mode"], seed=d["seed"], nang=d["nang"], nw=d["nw"], ordering="tau_major")
    jpath = os.path.join(os.path.dirname(path) or ".", f"census_bumps_v1_{OUT_TAG}_{os.path.splitext(os.path.basename(path))[0]}.json"); json.dump(out, open(jpath, "w"))
    print(f"\nSUMMARY  K={K} + bumps w={BUMP_W:g} nc={BUMP_NC}: price {price:.4f}  bump mass {mass_bumps:.3f}  active bumps {active.sum()} ({active.sum()/nflat:.1f}/cell)   [{time.time()-t0:.0f}s]")
    print(f"wrote {os.path.basename(jpath)}")

if __name__ == "__main__":
    main(sys.argv[1])
