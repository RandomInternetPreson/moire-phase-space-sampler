#!/usr/bin/env python3
"""
census_bandlimit_v1.py -- the smoothest census that attains the floor: minimal Fourier bandwidth of the birth.

Question.  The floor-attaining census is not a product in station coordinates, nor in sum/difference
coordinates, nor a short mixture of either (Results D-G).  Its (S,D) spectrum has power along at least
eleven distinct directions in harmonic space, so no linear change of coordinates makes it separable.
The right question is therefore not "which coordinates" but "how much structure": what is the smallest
Fourier bandwidth K such that a census built only from harmonics e^{2i(k_A xi_A + k_B xi_B)} with
|k_A|, |k_B| <= K reproduces the singlet at eps and pays Hall's floor?

Method.  For each (tau, scale) cell the census over the 1024 offsets is written as a trigonometric
polynomial  rho_ts(xi_A, xi_B) = sum_j c_{ts,j} phi_j(xi_A, xi_B)  over the real Fourier basis up to
degree K (phi_j = 1, cos/sin 2(k_A xi_A + k_B xi_B)).  The coefficients are FREE (sign-unconstrained);
non-negativity is imposed at every atom (rho >= 0 at the 1024 offsets), normalization is sum rho = 1,
and the correlator / CHSH / measurement-dependence constraints are exactly those of the paper's LP
(min_md_full in dither_sweep_rho.py).  This is ONE convex LP per K -- no starts, no alternation, no
local optima: infeasible means no band-limited census of that K exists on this terrain at this eps.

Outputs, per K: MD/floor or "infeasible", the number of active harmonics, and the leading harmonics of
the pooled census (magnitude, phase, direction in (S,D) harmonic space).  The minimal K and its census
are the formation-model target: the smoothest density a birth has to produce.

Usage (covariance_audit folder; rho file next to its banks_*.npz):
  set K_LIST=2,3,4,5,6
  set EPS_LIST=0.01
  set TIME_LIMIT=7200        (seconds per LP; 0 = none)
  python census_bandlimit_v1.py rho_nsh32_random0_eps0.01_nang161.json > census_bandlimit_v1_<date>.txt
K=6 has 169 basis functions per (tau,scale) cell (2535 coefficients in all) against 15360 atoms.
Writes census_bandlimit_v1_<stem>.json with the coefficients of each feasible census.
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

_deg = np.pi / 180; G = 256; RED = 256; DL = np.pi / G
lam_c = (np.arange(G) + 0.5) * np.pi / G
settings = [(0.0, 22.5), (0.0, 67.5), (45.0, 22.5), (45.0, 67.5)]; signs = [+1, -1, +1, +1]
S_QM = 2 * np.sqrt(2); FLOOR = (np.sqrt(2.0) - 1.0) / 3.0
REGISTER = os.environ.get("REGISTER", "paper")
K_LIST = [int(x) for x in os.environ.get("K_LIST", "2,3,4,5,6").split(",")]
EPS_LIST = [float(x) for x in os.environ.get("EPS_LIST", "0.01").split(",")]
TIME_LIMIT = float(os.environ.get("TIME_LIMIT", 7200))

def load(path):
    d = json.load(open(path)); offs = np.asarray(d["offsets"], float); NV = len(offs)
    cachef = os.path.join(os.path.dirname(path) or ".", f"banks_nsh{d['nsh']}_{d['mode']}{d['seed']}_nang{d['nang']}_nw{d['nw']}.npz")
    if not os.path.exists(cachef): sys.exit(f"banks checkpoint not found: {cachef}")
    z = np.load(cachef); nflat = int(z["nflat"])
    Qs = [np.vstack([z[f"q{k}_{p}"] for k in range(nflat)]) for p in range(4)]
    Svec = np.concatenate([z[f"s{k}"] for k in range(nflat)]); Ecoef = np.hstack([z[f"e{k}"] for k in range(nflat)])
    angs = np.linspace(2, 178, int(d["nang"])) if REGISTER == "paper" else np.linspace(6, 174, int(d["nang"]))
    return d, offs, nflat, Qs, Svec, Ecoef, angs

def md_of(Qs, rho):
    qm = [Qs[p].T @ rho for p in range(4)]
    return max(0.5 * np.sum(np.abs(qm[a] - qm[b])) * DL for a in range(4) for b in range(a + 1, 4))

def basis(offs, K):
    """real Fourier basis on the torus (period pi in each offset): columns = [1, cos/sin 2(kA xi_A + kB xi_B)]
       over one representative of each +-k pair; returns Phi (NV, nb) and the list of (kA, kB, 'c'/'s')."""
    xa, xb = offs[:, 0], offs[:, 1]; cols = [np.ones(len(offs))]; labels = [(0, 0, "1")]
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            arg = 2 * (ka * xa + kb * xb)
            cols.append(np.cos(arg)); labels.append((ka, kb, "c"))
            cols.append(np.sin(arg)); labels.append((ka, kb, "s"))
    return np.column_stack(cols), labels

def solve_bandlimit(Qs, Svec, Ecoef, angs, eps, Phi, nflat, NV):
    """one LP: variables c (nflat*nb, free), t, slacks s (6*RED >= 0).  rho = blockdiag(Phi) c >= 0."""
    nb = Phi.shape[1]; nc = nflat * nb; nang = len(angs)
    cmb = [(i, j) for i in range(4) for j in range(i + 1, 4)]; n_s = len(cmb) * RED
    nv = nc + 1 + n_s
    # reduced banks: for cell ts, Q_red[ts*nb + j, :] = sum_i Phi[i, j] Q[ts*NV + i, :]
    Qr = [np.vstack([Phi.T @ Q[ts * NV:(ts + 1) * NV] for ts in range(nflat)]) for Q in Qs]        # (nc, G)
    Sr = np.concatenate([Phi.T @ Svec[ts * NV:(ts + 1) * NV] for ts in range(nflat)])                 # (nc,)
    Er = np.hstack([Ecoef[:, ts * NV:(ts + 1) * NV] @ Phi for ts in range(nflat)])                    # (nang, nc)
    Mr = np.concatenate([Phi.sum(axis=0) for _ in range(nflat)])                                      # (nc,)
    Qc = [Q.reshape(nc, RED, -1).mean(axis=2) for Q in Qr]
    rows, bub = [], []; binw = np.pi / RED
    # MD rows (dense but few)
    blocks = []
    for ci, (i, j) in enumerate(cmb):
        Dif = (Qc[i] - Qc[j]) * binw                                                                   # (nc, RED)
        A1 = sp.hstack([sp.csr_matrix(Dif.T), sp.csr_matrix((RED, 1)), sp.csr_matrix((RED, n_s))]).tolil()
        A2 = sp.hstack([sp.csr_matrix(-Dif.T), sp.csr_matrix((RED, 1)), sp.csr_matrix((RED, n_s))]).tolil()
        for l in range(RED):
            A1[l, nc + 1 + ci * RED + l] = -1; A2[l, nc + 1 + ci * RED + l] = -1
        blocks += [A1.tocsr(), A2.tocsr()]; bub += [0.0] * (2 * RED)
    T = sp.lil_matrix((len(cmb), nv))
    for ci in range(len(cmb)):
        T[ci, nc + 1 + ci * RED:nc + 1 + (ci + 1) * RED] = 0.5; T[ci, nc] = -1
    blocks.append(T.tocsr()); bub += [0.0] * len(cmb)
    tgt = -np.cos(2 * angs * _deg)
    Ep = sp.hstack([sp.csr_matrix(Er), sp.csr_matrix((nang, 1 + n_s))]).tocsr()
    blocks += [Ep, -Ep]; bub += list(tgt + eps) + list(-tgt + eps)
    # positivity at every atom: -rho <= 0, rho = blockdiag(Phi) c
    P = sp.block_diag([sp.csr_matrix(-Phi)] * nflat).tocsr()
    blocks.append(sp.hstack([P, sp.csr_matrix((nflat * NV, 1 + n_s))]).tocsr()); bub += [0.0] * (nflat * NV)
    Aub = sp.vstack(blocks).tocsr(); bv = np.array(bub)
    Aeq = sp.csr_matrix(np.vstack([np.concatenate([Mr, [0], np.zeros(n_s)]), np.concatenate([Sr, [0], np.zeros(n_s)])]))
    beq = [1.0, -S_QM]
    cobj = np.zeros(nv); cobj[nc] = 1
    bounds = [(None, None)] * nc + [(0, None)] * (1 + n_s)
    opts = {"time_limit": TIME_LIMIT} if TIME_LIMIT > 0 else None
    res = linprog(cobj, A_ub=Aub, b_ub=bv, A_eq=Aeq, b_eq=beq, bounds=bounds, method="highs", options=opts)
    if res.status != 0:
        return ("time limit" if res.status == 1 else "infeasible"), None, None
    c = res.x[:nc].reshape(nflat, nb)
    rho = np.concatenate([Phi @ c[ts] for ts in range(nflat)])
    return res.x[nc], c, rho

def main(path):
    t0 = time.time()
    d, offs, nflat, Qs, Svec, Ecoef, angs = load(path); NV = len(offs); NT = nflat * NV
    rho_star = np.asarray(d["rho"], float); rho_star /= rho_star.sum(); md_star = md_of(Qs, rho_star)
    print(f"file: {os.path.basename(path)}   atoms {NT} = {nflat} x {NV}   K list {K_LIST}   eps {EPS_LIST}")
    print(f"check: MD of the stored rho* on these banks = {md_star:.12f}   (file says {d['md']:.12f}; floor {FLOOR:.12f})")
    out = {"file": os.path.basename(path), "md_star": md_star, "results": []}
    for eps in EPS_LIST:
        for K in K_LIST:
            Phi, labels = basis(offs, K); nb = Phi.shape[1]
            print(f"\n=== eps={eps}  K={K}   ({nb} basis functions per (tau,scale) cell, {nflat*nb} coefficients) ===", flush=True)
            tk = time.time()
            t, c, rho = solve_bandlimit(Qs, Svec, Ecoef, angs, eps, Phi, nflat, NV)
            if c is None:
                print(f"  {t}   [{time.time()-tk:.0f}s]", flush=True)
                out["results"].append({"eps": eps, "K": K, "status": t}); continue
            assert rho.min() > -1e-7 and abs(rho.sum() - 1) < 1e-7
            md = md_of(Qs, np.clip(rho, 0, None))
            assert md >= FLOOR - 1e-7
            # pooled spectrum of the band-limited census over offsets
            wo = rho.reshape(nflat, NV).sum(axis=0); wo = wo / wo.sum()
            xa, xb = offs[:, 0], offs[:, 1]
            pix = []
            for ka in range(-K, K + 1):
                for kb in range(-K, K + 1):
                    if ka < 0 or (ka == 0 and kb <= 0): continue
                    F = np.sum(wo * np.exp(2j * (ka * xa + kb * xb)))
                    pix.append((abs(F), ka, kb, np.degrees(np.angle(F))))
            pix.sort(reverse=True)
            nact = int(np.sum(np.abs(c) > 1e-6 * np.abs(c).max()))
            print(f"  MD/floor = {md/FLOOR:.4f}   feasible   active coefficients {nact}/{nflat*nb}   [{time.time()-tk:.0f}s]")
            print("  leading harmonics of the pooled census  (|F|, k_A, k_B, phi[deg])  -- (k_S,k_D) = ((k_A+k_B)/2, (k_A-k_B)/2):")
            for a, ka, kb, ph in pix[:10]:
                print(f"     {a:.3f}  ({ka:+d},{kb:+d})  {ph:+8.2f}    (S,D)-harmonic ({(ka+kb)/2:+.1f},{(ka-kb)/2:+.1f})")
            out["results"].append({"eps": eps, "K": K, "status": "feasible", "price": md / FLOOR, "active": nact,
                                   "leading": [(float(a), int(ka), int(kb), float(ph)) for a, ka, kb, ph in pix[:20]],
                                   "coeff_labels": [list(map(str, l)) for l in labels], "coeff": c.tolist()})
            json.dump(out, open(os.path.join(os.path.dirname(path) or ".", f"census_bandlimit_v1_{os.path.splitext(os.path.basename(path))[0]}.json"), "w"), default=float)
    print("\nSUMMARY  (MD/floor of the smoothest census of bandwidth K; '-' = infeasible)")
    print(f"{'eps':>6} {'K':>3} {'price':>9}")
    for r in out["results"]:
        print(f"{r['eps']:>6} {r['K']:>3} " + (f"{r['price']:9.4f}" if r.get("price") is not None else f"   {r['status']}"))
    print(f"\nreference: rho* (unrestricted, 1024 atoms per cell) MD/floor = {md_star/FLOOR:.6f}   total {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main(sys.argv[1])
