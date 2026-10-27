#!/usr/bin/env python3
"""
census_cellsep_v2.py -- the price of a birth that is dial-blind AND smooth.

v2 adds BAND=K: the shared shapes f_r are band-limited trigonometric polynomials of degree K (free coefficients,
positivity at every atom), so the class is "same offset density at every dial setting, smooth to K harmonics".
BAND=0 is v1 (free shapes).  v1 result on random0 (25 Oct): free rank-1 misses 0.22%, price 1.836; rank-2 price 1.548.
Band-limited (not separable) prices: K=6 2.274, K=8 1.861, K=10 1.598.  The question: do the two prices compound?

--- v1 docstring follows ---
census_cellsep_v1.py -- can a birth that stamps the SAME offsets at every dial setting attain the floor?

The class.  A cell-separable census of rank R:
    rho(cell, xi) = sum_{r=1..R} W_r[cell] * f_r(xi),     W_r >= 0, f_r >= 0,  sum_v f_r = 1,  sum_{r,c} W_r[c] = 1
with the shapes f_r FREE over the 1024 offsets (no parametric form) and the weights W_r free over the 15 (tau, scale)
cells.  R = 1 is "one offset density, reweighted per cell" -- the assumption of every birth model so far (v1, v2).
R = 2, 3 are mixtures of two or three dial-independent shapes.  census_cells_v1 showed the K = 6 LP census is far from
this class (best rank-1 approximation misses 5%), but that was an approximation of ONE census; this script optimizes
over the whole class.

Method.  Bilinear in (W, f), so alternating LPs from random starts (as product_census_v4 did for station-locality):
  miss step   W fixed -> LP in (f, sigma): min sigma  s.t. |E - target| <= eps + sigma, CHSH exact, masses;  then f fixed -> LP in (W, sigma)
  price step  (only once feasible) same alternation with the paper's MD objective (min_md_full rows) at the band eps
Local optima, so the result is an UPPER bound on the class's minimum miss; NSTART starts and the spread across them are
reported.  Compare: best K=5 band-limited census misses 0.0120 (a certificate); birth models v1/v2 (separable by form) 0.024.

Usage (covariance_audit):  python census_cellsep_v1.py rho_nsh32_random0_eps0.01_nang161.json > census_cellsep_random0_<date>.txt
Env: RANKS (default 1,2,3), NSTART (6), ITERS (40), EPS (0.01), SEED (1), REGISTER (paper).  scipy/HiGHS; minutes.
Writes census_cellsep_v1_<stem>.json with, per rank, the best W (cells) and f (offsets) and their misses/prices.
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

_deg = np.pi / 180; G = 256; RED = 256; DL = np.pi / G
S_QM = 2 * np.sqrt(2); FLOOR = (np.sqrt(2.0) - 1.0) / 3.0
RANKS = [int(x) for x in os.environ.get("RANKS", "1,2").split(",")]; BAND = int(os.environ.get("BAND", 0)); NSTART = int(os.environ.get("NSTART", 6))
ITERS = int(os.environ.get("ITERS", 40)); EPS = float(os.environ.get("EPS", 0.01)); SEED = int(os.environ.get("SEED", 1))
REGISTER = os.environ.get("REGISTER", "paper"); OUT_TAG = os.environ.get("OUT_TAG", ""); TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]; SCALES = [1.0, 1.8, 2.6]

def load(path):
    d = json.load(open(path)); offs = np.asarray(d["offsets"], float); NV = len(offs)
    cachef = os.path.join(os.path.dirname(path) or ".", f"banks_nsh{d['nsh']}_{d['mode']}{d['seed']}_nang{d['nang']}_nw{d['nw']}.npz")
    if not os.path.exists(cachef): sys.exit(f"banks checkpoint not found: {cachef}")
    z = np.load(cachef); nflat = int(z["nflat"])
    Qs = [np.vstack([z[f"q{k}_{p}"] for k in range(nflat)]) for p in range(4)]
    Svec = np.concatenate([z[f"s{k}"] for k in range(nflat)]); Ecoef = np.hstack([z[f"e{k}"] for k in range(nflat)])
    angs = np.linspace(2, 178, int(d["nang"])) if REGISTER == "paper" else np.linspace(6, 174, int(d["nang"]))
    return d, offs, nflat, Qs, Svec, Ecoef, -np.cos(2 * angs * _deg)

def md_of(Qs, rho):
    qm = [Qs[p].T @ rho for p in range(4)]
    return max(0.5 * np.sum(np.abs(qm[a] - qm[b])) * DL for a in range(4) for b in range(a + 1, 4))

def basis(offs, K):
    xa, xb = offs[:, 0], offs[:, 1]; cols = [np.ones(len(offs))]
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            arg = 2 * (ka * xa + kb * xb); cols.append(np.cos(arg)); cols.append(np.sin(arg))
    return np.column_stack(cols)

def band_maps(Er, Sr, Qr, Meq, meq, Phi, R):
    """f-step with f_r = Phi c_r: re-express the maps in the coefficients c (R*nb, free), add positivity rows -Phi c <= 0."""
    BD = sp.block_diag([sp.csr_matrix(Phi)] * R).tocsr()
    return Er @ BD, BD.T @ Sr, [BD.T @ Q for Q in Qr], Meq @ BD, meq, (-BD).tocsr()

# ---- the linear maps from the step's variables to the atoms ----
def maps_f(W, nflat, NV, Qs, Svec, Ecoef):
    """W fixed (R x nflat).  Variables f (R*NV).  rho[c*NV+v] = sum_r W[r,c] f[r*NV+v]."""
    R = W.shape[0]
    Er = np.hstack([sum(W[r, c] * Ecoef[:, c * NV:(c + 1) * NV] for c in range(nflat)) for r in range(R)])        # (nang, R*NV)
    Sr = np.concatenate([sum(W[r, c] * Svec[c * NV:(c + 1) * NV] for c in range(nflat)) for r in range(R)])
    Qr = [np.vstack([sum(W[r, c] * Q[c * NV:(c + 1) * NV] for c in range(nflat)) for r in range(R)]) for Q in Qs]  # (R*NV, G)
    Meq = np.zeros((R, R * NV))
    for r in range(R): Meq[r, r * NV:(r + 1) * NV] = 1.0
    return Er, Sr, Qr, Meq, np.ones(R)
def maps_W(F, nflat, NV, Qs, Svec, Ecoef):
    """f fixed (R x NV, each summing to 1).  Variables W (R*nflat).  rho[c*NV+v] = sum_r W[r*nflat+c] f[r,v]."""
    R = F.shape[0]
    Er = np.column_stack([Ecoef[:, c * NV:(c + 1) * NV] @ F[r] for r in range(R) for c in range(nflat)])             # (nang, R*nflat)
    Sr = np.array([Svec[c * NV:(c + 1) * NV] @ F[r] for r in range(R) for c in range(nflat)])
    Qr = [np.vstack([Q[c * NV:(c + 1) * NV].T @ F[r] for r in range(R) for c in range(nflat)]) for Q in Qs]        # (R*nflat, G)
    return Er, Sr, Qr, np.ones((1, R * nflat)), np.ones(1)

def miss_lp(Er, Sr, Meq, meq, tgt, eps, Apos=None):
    """Apos (sparse, rows x n): extra rows Apos x <= 0 (positivity of a band-limited shape); then x is free, else x >= 0."""
    n = Er.shape[1]; nang = Er.shape[0]
    Aub = sp.vstack([sp.hstack([sp.csr_matrix(Er), sp.csr_matrix(-np.ones((nang, 1)))]), sp.hstack([sp.csr_matrix(-Er), sp.csr_matrix(-np.ones((nang, 1)))])]).tocsr()
    bub = np.concatenate([tgt + eps, -tgt + eps])
    if Apos is not None:
        Aub = sp.vstack([Aub, sp.hstack([Apos, sp.csr_matrix((Apos.shape[0], 1))])]).tocsr(); bub = np.concatenate([bub, np.zeros(Apos.shape[0])])
    Aeq = np.vstack([np.hstack([Meq, np.zeros((Meq.shape[0], 1))]), np.concatenate([Sr, [0.0]])[None, :]])
    beq = np.concatenate([meq, [-S_QM]])
    c = np.zeros(n + 1); c[n] = 1
    bounds = ([(None, None)] if Apos is not None else [(0, None)]) * n + [(0, None)]
    r = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds, method="highs")
    return (r.x[:n], r.x[n]) if r.status == 0 else (None, np.inf)

def price_lp(Er, Sr, Qr, Meq, meq, tgt, eps, Apos=None):
    """paper's objective: min t with MD rows (reduced banks), band at eps, CHSH exact, masses; variables >= 0 (or free with positivity rows Apos)."""
    n = Er.shape[1]; nang = Er.shape[0]
    cmb = [(i, j) for i in range(4) for j in range(i + 1, 4)]; n_s = len(cmb) * RED; nv = n + 1 + n_s
    Qc = [Q.reshape(n, RED, -1).mean(axis=2) for Q in Qr]; binw = np.pi / RED
    blocks, bub = [], []
    for ci, (i, j) in enumerate(cmb):
        Dif = (Qc[i] - Qc[j]) * binw
        for sign in (1, -1):
            A = sp.lil_matrix((RED, nv)); A[:, :n] = sign * Dif.T
            for l in range(RED): A[l, n + 1 + ci * RED + l] = -1
            blocks.append(A.tocsr()); bub += [0.0] * RED
    T = sp.lil_matrix((len(cmb), nv))
    for ci in range(len(cmb)): T[ci, n + 1 + ci * RED:n + 1 + (ci + 1) * RED] = 0.5; T[ci, n] = -1
    blocks.append(T.tocsr()); bub += [0.0] * len(cmb)
    Ep = sp.hstack([sp.csr_matrix(Er), sp.csr_matrix((nang, 1 + n_s))]).tocsr(); blocks += [Ep, -Ep]; bub += list(tgt + eps) + list(-tgt + eps)
    if Apos is not None:
        blocks.append(sp.hstack([Apos, sp.csr_matrix((Apos.shape[0], 1 + n_s))]).tocsr()); bub += [0.0] * Apos.shape[0]
    Aub = sp.vstack(blocks).tocsr()
    Aeq = sp.csr_matrix(np.vstack([np.hstack([Meq, np.zeros((Meq.shape[0], 1 + n_s))]), np.concatenate([Sr, np.zeros(1 + n_s)])[None, :]]))
    c = np.zeros(nv); c[n] = 1
    bounds = ([(None, None)] if Apos is not None else [(0, None)]) * n + [(0, None)] * (1 + n_s)
    r = linprog(c, A_ub=Aub, b_ub=np.array(bub), A_eq=Aeq, b_eq=np.concatenate([meq, [-S_QM]]), bounds=bounds, method="highs")
    return (r.x[:n], r.x[n]) if r.status == 0 else (None, np.inf)

def rho_of(W, F, nflat, NV):
    return np.concatenate([sum(W[r, c] * F[r] for r in range(W.shape[0])) for c in range(nflat)])

def misses(rho, Ecoef, Svec, tgt):
    return float(np.abs(Ecoef @ rho - tgt).max()), float(abs(Svec @ rho + S_QM))

def spectrum(w, offs, K=4):
    xa, xb = offs[:, 0], offs[:, 1]; pix = []
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            Fv = np.sum(w * np.exp(2j * (ka * xa + kb * xb))) / w.sum(); pix.append((abs(Fv), ka, kb, np.degrees(np.angle(Fv))))
    pix.sort(reverse=True); return pix

def main(path):
    t0 = time.time(); rng = np.random.default_rng(SEED)
    d, offs, nflat, Qs, Svec, Ecoef, tgt = load(path); NV = len(offs)
    rho_star = np.asarray(d["rho"], float); rho_star /= rho_star.sum()
    print(f"file: {os.path.basename(path)}   atoms {nflat*NV} = {nflat} x {NV};  rho* misses: band {misses(rho_star, Ecoef, Svec, tgt)[0]:.5f}, price {md_of(Qs, rho_star)/FLOOR:.4f}")
    Phi = basis(offs, BAND) if BAND > 0 else None
    print(f"class: cell-separable censuses of rank R, shapes {'band-limited to K=' + str(BAND) + ' (' + str(Phi.shape[1]) + ' harmonics)' if BAND else 'free'} over {NV} offsets;  ranks {RANKS}, {NSTART} starts, {ITERS} alternations, eps {EPS}")
    out = {"file": os.path.basename(path), "eps": EPS, "band": BAND, "ranks": {}}
    def fstep(W, price=False):
        Er, Sr, Qr, Meq, meq = maps_f(W, nflat, NV, Qs, Svec, Ecoef); R = W.shape[0]; Apos = None
        if Phi is not None: Er, Sr, Qr, Meq, meq, Apos = band_maps(Er, Sr, Qr, Meq, meq, Phi, R)
        x, val = price_lp(Er, Sr, Qr, Meq, meq, tgt, EPS, Apos) if price else miss_lp(Er, Sr, Meq, meq, tgt, 0.0, Apos)
        if x is None: return None, val
        F = np.vstack([Phi @ x[r * Phi.shape[1]:(r + 1) * Phi.shape[1]] for r in range(R)]) if Phi is not None else x.reshape(R, NV)
        return np.clip(F, 0, None), val
    for R in RANKS:
        print(f"\n=== rank {R} ===", flush=True); best = (np.inf, None, None)
        for s in range(NSTART):
            W = rng.dirichlet(np.ones(nflat), size=R) * rng.dirichlet(np.ones(R))[:, None]   # sum_{r,c} W = 1
            F = None; hist = []; sig = np.inf
            for it in range(ITERS):
                F2, sig_f = fstep(W)                                   # eps=0: sigma IS the miss
                if F2 is None: hist.append("f-step infeasible (CHSH row unattainable for these weights)"); break
                F = F2
                Er, Sr, Qr, Meq, meq = maps_W(F, nflat, NV, Qs, Svec, Ecoef)
                w, sig_w = miss_lp(Er, Sr, Meq, meq, tgt, 0.0)
                if w is None: hist.append("W-step infeasible"); break
                W = w.reshape(R, nflat); hist.append(sig_w)
                if sig - sig_w < 1e-7 and it > 2: sig = sig_w; break
                sig = sig_w
            if F is None or not np.isfinite(sig): print(f"  start {s}: {hist[-1]}"); continue
            rho = rho_of(W, F, nflat, NV); mE, mS = misses(rho, Ecoef, Svec, tgt)
            print(f"  start {s}: sigma* {sig:.6f} -> min miss {mE:.5f} (CHSH {mS:.1e}) after {len(hist)} alternations   [{time.time()-t0:.0f}s]", flush=True)
            if mE < best[0]: best = (mE, W.copy(), F.copy())
        mE, W, F = best
        if W is None: print(f"  rank {R}: no start produced a feasible alternation"); out["ranks"][R] = {"status": "no feasible start"}; continue
        rec = {"min_miss": float(mE), "W": W.tolist(), "F": F.tolist()}
        print(f"\n  rank {R}: best min miss {mE:.5f}   {'<= eps: FEASIBLE at eps' if mE <= EPS + 1e-9 else '> eps: not feasible at eps (upper bound on the class minimum)'}")
        for r in range(R):
            Wr = W[r]; print(f"    shape {r}: total weight {Wr.sum():.3f}; cell weights (rows tau, cols scale):")
            if nflat == 15:
                for ti, tau in enumerate(TAUS): print(f"      tau={tau:<4}" + "".join(f"  {Wr[ti*3+si]:.3f}" for si in range(3)))
            else: print("      " + " ".join(f"{x:.3f}" for x in Wr))
            pix = spectrum(F[r], offs); print("      leading harmonics (|F|, k_A, k_B, phi):  " + "  ".join(f"({ka:+d},{kb:+d}) {a:.2f}@{ph:+.0f}" for a, ka, kb, ph in pix[:5]))
        if mE <= EPS + 1e-9:
            price = np.inf
            for it in range(6):
                F2, t = fstep(W, price=True)
                if F2 is None: break
                F = F2
                Er, Sr, Qr, Meq, meq = maps_W(F, nflat, NV, Qs, Svec, Ecoef); w, t2 = price_lp(Er, Sr, Qr, Meq, meq, tgt, EPS)
                if w is None: break
                W = w.reshape(R, nflat)
                if price - t2 < 1e-6: price = t2; break
                price = t2
            rho = rho_of(W, F, nflat, NV); pr = md_of(Qs, rho) / FLOOR; mE2, mS2 = misses(rho, Ecoef, Svec, tgt)
            print(f"    price at eps={EPS}: MD/floor = {pr:.4f} (re-measured; band miss {mE2:.5f})   [{time.time()-t0:.0f}s]")
            rec.update(price=float(pr), W=W.tolist(), F=F.tolist())
        out["ranks"][R] = rec
    print(f"\nSUMMARY  cell-separable censuses ({'K=' + str(BAND) if BAND else 'free shapes'}), min miss by rank:  " + "  ".join(f"R={R}: {out['ranks'][R].get('min_miss', float('nan')):.5f}" + (f" (price {out['ranks'][R]['price']:.4f})" if 'price' in out['ranks'][R] else "") for R in RANKS))
    print(f"compare: K=5 band-limited (certified) 0.01196; birth models v1/v2 (separable by form) 0.0242; rank-1 approx of the K=6 census 0.0499   [{time.time()-t0:.0f}s total]")
    json.dump(out, open(os.path.join(os.path.dirname(path) or ".", f"census_cellsep_v2_{OUT_TAG + '_' if OUT_TAG else ''}{os.path.splitext(os.path.basename(path))[0]}.json"), "w"))

if __name__ == "__main__":
    main(sys.argv[1])
