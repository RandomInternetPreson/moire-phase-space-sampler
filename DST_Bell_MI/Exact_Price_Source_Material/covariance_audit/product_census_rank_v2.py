#!/usr/bin/env python3
"""
product_census_rank_v2.py -- how many independent modes does the birth need?

Rank-R census classes on the saved banks, in either coordinate system:
    w(k, i) = sum_{r=1..R}  q_r[k] * f_r[C1(i)] * g_r[C2(i)]        (all factors >= 0)
with (C1, C2) = (xi_A, xi_B) for COORDS=ab or (S, D) = (xi_A + xi_B, xi_A - xi_B) for COORDS=sumdiff,
binned NB x NB; k indexes (tau, scale).  R = 1 is the product class of product_census_v3/v4.  Rank R
is the class "the census is a mixture of R separable births".  The smallest R that attains the floor
at eps = 0.01 is the number of independent modes a formation model must contain.

Solved by alternating penalty LPs (blocks q, f, g, each jointly over all R modes; all starts random
except the first, which is uniform), feasible = slacks < FEAS_TOL; Hall's bound asserted on every
reported census.  v2 (same math as v1): prints the correlator MISS = eps + slack next to the slack, since the
slack is the excess over the eps band and was misreported as the miss in earlier prose.  Same machinery as product_census_v2-v4 (normalization via group mass), extended to
sums of products.

Usage (covariance_audit folder, rho file next to its banks_*.npz):
  set COORDS=sumdiff          (or ab)
  set RANK_LIST=1,2,3
  set NB_LIST=32
  set EPS_LIST=0.05,0.01
  set NSTART=8
  set MAXIT=40
  set SEED=11
  python product_census_rank_v1.py rho_nsh32_random0_eps0.01_nang161.json > product_census_rank_<coords>_<date>.txt
No joint control (known: 1.12 / 1.11 at eps=0.01 on random0, NB=32).  Writes product_census_rank_v1_<coords>_<stem>.json.
"""
import os, sys, json, time
import numpy as np
from scipy.optimize import linprog

_deg = np.pi / 180; G = 256; RED = 256; DL = np.pi / G
lam_c = (np.arange(G) + 0.5) * np.pi / G
settings = [(0.0, 22.5), (0.0, 67.5), (45.0, 22.5), (45.0, 67.5)]; signs = [+1, -1, +1, +1]
S_QM = 2 * np.sqrt(2); FLOOR = (np.sqrt(2.0) - 1.0) / 3.0
REGISTER = os.environ.get("REGISTER", "paper")
COORDS = os.environ.get("COORDS", "sumdiff"); assert COORDS in ("ab", "sumdiff")
RANK_LIST = [int(x) for x in os.environ.get("RANK_LIST", "1,2,3").split(",")]
NB_LIST = [int(x) for x in os.environ.get("NB_LIST", "32").split(",")]
EPS_LIST = [float(x) for x in os.environ.get("EPS_LIST", "0.05,0.01").split(",")]
NSTART = int(os.environ.get("NSTART", 8)); MAXIT = int(os.environ.get("MAXIT", 40))
PEN = float(os.environ.get("PEN", 100.0)); FEAS_TOL = float(os.environ.get("FEAS_TOL", 1e-6))
SEED = int(os.environ.get("SEED", 11))

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

def reduce_multi(Qs, Svec, Ecoef, grp, cs, nth):
    """w_v = sum_r theta[r*nth + grp[v]] * cs[r][v]  ->  banks over the R*nth theta entries"""
    R = len(cs); N = R * nth
    Qr = [np.zeros((N, G)) for _ in range(4)]; Sr = np.zeros(N); Er = np.zeros((Ecoef.shape[0], N)); Mr = np.zeros(N)
    for r, c in enumerate(cs):
        g = r * nth + grp
        for p in range(4): np.add.at(Qr[p], g, Qs[p] * c[:, None])
        np.add.at(Sr, g, Svec * c); np.add.at(Er.T, g, (Ecoef * c[None, :]).T); np.add.at(Mr, g, c)
    return Qr, Sr, Er, Mr

def solve_lp(Qs, Svec, Ecoef, Mr, angs, eps):
    keep = np.where(Mr > 0)[0]; full_n = len(Mr)
    Qs = [Q[keep] for Q in Qs]; Svec = Svec[keep]; Ecoef = Ecoef[:, keep]; Mr = Mr[keep]
    NV = len(Svec); nang = len(angs)
    Qc = [Q.reshape(NV, RED, -1).mean(axis=2) for Q in Qs]
    cmb = [(i, j) for i in range(4) for j in range(i + 1, 4)]; n_s = len(cmb) * RED; nu = nang + 1
    nv = NV + 1 + n_s + nu; rows, bub = [], []; binw = np.pi / RED
    for ci, (i, j) in enumerate(cmb):
        Dif = (Qc[i] - Qc[j]) * binw
        for l in range(RED):
            r = np.zeros(nv); r[:NV] = Dif[:, l]; r[NV + 1 + ci * RED + l] = -1; rows.append(r); bub.append(0.0)
            r = np.zeros(nv); r[:NV] = -Dif[:, l]; r[NV + 1 + ci * RED + l] = -1; rows.append(r); bub.append(0.0)
    for ci in range(len(cmb)):
        r = np.zeros(nv); r[NV + 1 + ci * RED:NV + 1 + (ci + 1) * RED] = 0.5; r[NV] = -1; rows.append(r); bub.append(0.0)
    tgt = -np.cos(2 * angs * _deg); u0 = NV + 1 + n_s
    for ai in range(nang):
        r = np.zeros(nv); r[:NV] = Ecoef[ai]; r[u0 + ai] = -1; rows.append(r); bub.append(tgt[ai] + eps)
        r = np.zeros(nv); r[:NV] = -Ecoef[ai]; r[u0 + ai] = -1; rows.append(r); bub.append(-tgt[ai] + eps)
    r = np.zeros(nv); r[:NV] = Svec; r[u0 + nang] = -1; rows.append(r); bub.append(-S_QM)
    r = np.zeros(nv); r[:NV] = -Svec; r[u0 + nang] = -1; rows.append(r); bub.append(S_QM)
    Aeq = np.zeros((1, nv)); Aeq[0, :NV] = Mr
    cobj = np.zeros(nv); cobj[NV] = 1; cobj[u0:] = PEN
    res = linprog(cobj, A_ub=np.array(rows), b_ub=np.array(bub), A_eq=Aeq, b_eq=[1.0], bounds=[(0, None)] * nv, method="highs")
    if res.status != 0: return None, None, None
    theta = np.zeros(full_n); theta[keep] = res.x[:NV]
    return res.x[NV], theta, float(res.x[u0:].max())

def alternate_rank(Qs, Svec, Ecoef, angs, eps, blocks, R, init, maxit):
    """blocks: [(name, grp, nth)], theta[name] shape (R, nth)"""
    th = {n: init[n].copy() for n, _, _ in blocks}; NT = len(blocks[0][1]); last = None
    for it in range(maxit):
        for name, grp, nth in blocks:
            cs = []
            for r in range(R):
                c = np.ones(NT)
                for n2, g2, _ in blocks:
                    if n2 != name: c *= th[n2][r][g2]
                cs.append(c)
            Qr, Sr, Er, Mr = reduce_multi(Qs, Svec, Ecoef, grp, cs, nth)
            t, theta, viol = solve_lp(Qr, Sr, Er, Mr, angs, eps)
            if t is None: return None
            th[name] = theta.reshape(R, nth)
        w = np.zeros(NT)
        for r in range(R):
            c = np.ones(NT)
            for n2, g2, _ in blocks: c *= th[n2][r][g2]
            w += c
        assert abs(w.sum() - 1) < 1e-6, w.sum()
        md = md_of(Qs, w)
        if viol < FEAS_TOL: assert md >= FLOOR - 1e-7, f"Hall bound violated: {md}"
        if last is not None and abs(md - last[0]) < 1e-7 and abs(viol - last[1]) < 1e-9: break
        last = (md, viol)
    return md, w, viol, it + 1, th

def main(path):
    t0 = time.time()
    d, offs, nflat, Qs, Svec, Ecoef, angs = load(path); NV = len(offs); NT = nflat * NV
    rho_star = np.asarray(d["rho"], float); rho_star /= rho_star.sum()
    md_star = md_of(Qs, rho_star)
    print(f"file: {os.path.basename(path)}   atoms {NT}   coords {COORDS}   ranks {RANK_LIST}")
    print(f"check: MD of the stored rho* on these banks = {md_star:.12f}   (file says {d['md']:.12f}; floor {FLOOR:.12f})")
    k_of = np.repeat(np.arange(nflat), NV); rng = np.random.default_rng(SEED)
    if COORDS == "ab": c1, c2 = offs[:, 0] % np.pi, offs[:, 1] % np.pi
    else: c1, c2 = (offs[:, 0] + offs[:, 1]) % np.pi, (offs[:, 0] - offs[:, 1]) % np.pi
    out = {"file": os.path.basename(path), "coords": COORDS, "md_star": md_star, "results": []}
    for NB in NB_LIST:
        A_of = np.tile((c1 / np.pi * NB).astype(int) % NB, nflat); B_of = np.tile((c2 / np.pi * NB).astype(int) % NB, nflat)
        blocks = [("q", k_of, nflat), ("f", A_of, NB), ("g", B_of, NB)]
        for eps in EPS_LIST:
            for R in RANK_LIST:
                print(f"\n=== NB={NB}  eps={eps}  rank={R} ===", flush=True)
                best = None; rows = []
                for s in range(NSTART):
                    if s == 0: init = dict(q=np.ones((R, nflat)), f=np.ones((R, NB)), g=np.ones((R, NB)))
                    else: init = dict(q=rng.random((R, nflat)), f=rng.random((R, NB)), g=rng.random((R, NB)))
                    ta = time.time(); r = alternate_rank(Qs, Svec, Ecoef, angs, eps, blocks, R, init, MAXIT)
                    if r is None: print(f"  start {s}: LP failed"); rows.append((s, None)); continue
                    md, w, viol, iters, th = r; feas = viol < FEAS_TOL
                    rows.append((s, md / FLOOR, viol, iters, feas))
                    print(f"  start {s:2d}: MD/floor = {md/FLOOR:.4f}   max slack {viol:.2e} (correlator miss = eps + slack = {eps+viol:.4f}) {'feasible' if feas else 'NOT feasible'}   {iters} rounds  [{time.time()-ta:.0f}s]", flush=True)
                    if feas and (best is None or md < best[0]): best = (md, w, s, th)
                rec = {"NB": NB, "eps": eps, "rank": R, "starts": rows}
                if best is None:
                    minviol = min((x[2] for x in rows if x[1] is not None), default=None)
                    print(f"  RESULT rank {R}: no feasible census at eps={eps}  (smallest slack {minviol:.2e}; smallest correlator miss {eps+minviol:.4f})" if minviol is not None else f"  RESULT rank {R}: LP failures")
                    rec["price"] = None; rec["min_slack"] = minviol
                else:
                    md, w, s, th = best
                    # mode masses: share of the census carried by each mode
                    masses = []
                    for r in range(R):
                        c = np.ones(NT)
                        for n2, g2, _ in blocks: c *= th[n2][r][g2]
                        masses.append(float(c.sum()))
                    print(f"  RESULT rank {R}: best feasible MD/floor = {md/FLOOR:.4f} (start {s}); mode masses {np.round(masses, 3).tolist()}")
                    rec["price"] = md / FLOOR; rec["mode_masses"] = masses
                out["results"].append(rec)
                json.dump(out, open(os.path.join(os.path.dirname(path) or ".", f"product_census_rank_v2_{COORDS}_{os.path.splitext(os.path.basename(path))[0]}.json"), "w"), indent=1, default=float)
    print(f"\nSUMMARY (MD/floor; '-' = no feasible census)   coords = {COORDS}")
    print(f"{'NB':>4} {'eps':>6} {'rank':>5} {'price':>9}  min slack if infeasible")
    for r in out["results"]:
        print(f"{r['NB']:>4} {r['eps']:>6} {r['rank']:>5} " + (f"{r['price']:9.4f}" if r["price"] is not None else "    -    ") + (f"  {r['min_slack']:.2e}" if r.get("min_slack") is not None else ""))
    print(f"\nreference: rho* on these banks MD/floor = {md_star/FLOOR:.6f}   total {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main(sys.argv[1])
