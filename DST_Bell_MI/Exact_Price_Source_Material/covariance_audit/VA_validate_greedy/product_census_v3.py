#!/usr/bin/env python3
"""
product_census_v3.py -- Can a STATION-LOCAL (product) census attain Hall's floor?

Question.  The LP-optimal census rho* of arXiv:2608.18886 is far from a product over the two
stations (|F(2,2)| = 0.37 vs the product prediction 0.02).  Any station-local origin of the
frozen offsets (lattice, collection cone, disorder) can only supply a census of the form
rho(xi_A, xi_B) = f(xi_A) g(xi_B).  This script asks what such a census costs, on the SAME
trajectory banks the floor solution used, so the comparison is like for like.

Three classes of census on the same 1024 offsets x 15 (tau, scale) atoms, with the offsets
binned into NB x NB cells of the (xi_A, xi_B) torus (period pi):

  joint    w(k, i) = W[k, A(i), B(i)]            free over 15*NB*NB cells      -- ONE LP (exact)
  semi     w(k, i) = H[k, A(i)] * g[B(i)]        station B independent of all  -- alternating LPs
  product  w(k, i) = q[k] * f[A(i)] * g[B(i)]    fully factorized              -- alternating LPs

k indexes (tau, scale) in the banks' tau-major order; A(i), B(i) are the bins of offset i.
'joint' is the control: it says what the binning alone costs.  'semi' and 'product' are
non-convex (rank-1) classes solved by alternation from several starts; a price of 1.000 is an
existence proof, a high price is "no product census found", not a theorem.  The alternation
uses a penalty form (slack on the correlator and CHSH constraints, weight PEN) so that it is
always feasible; a result counts as feasible at eps only if the final slacks are < FEAS_TOL.

Reads:  the rho_*.json written by dither_sweep_rho.py (for offsets and md) and the matching
        banks_nsh*_nang*_nw*.npz checkpoint (same directory).
Writes: product_census_v2_<rho-stem>.json

v2 (Oct 17): fixes v1, whose normalization row in the reduced LPs read sum(theta) = 1 instead of
sum(w) = 1 (w_v = theta[grp(v)] * c_v); v1 solutions had sum(w) = 1 + eps and prices below Hall's
floor, which is impossible.  v2 carries the group mass m_j = sum_{v in j} c_v into the equality
row and asserts sum(w) = 1 and the Hall bound on every reported census.
v3 (Oct 17): the joint control LP at NB=32 is as large as the paper's full LP and can run for hours;
v3 gives it a time limit (JOINT_TIME_LIMIT seconds, default 1800; 0 = no limit) and a SKIP_JOINT=1
switch, drops zero-mass columns from every LP, and otherwise is v2 unchanged. and a printed report.

Usage (from the covariance_audit folder):
  set REGISTER=paper         (the [2,178] register the banks were built with; default)
  set NB_LIST=16             (bin counts to try, comma separated; default 16)
  set EPS_LIST=0.1,0.05,0.02,0.01
  set NSTART=6               (random starts per class, plus the two deterministic starts)
  set MAXIT=40
  set JOINT_TIME_LIMIT=1800   (seconds for the joint control LP; 0 = unlimited)
  set SKIP_JOINT=0
  python product_census_v3.py rho_nsh32_random0_eps0.01_nang161.json
"""
import os, sys, json, time
import numpy as np
from scipy.optimize import linprog

_deg = np.pi / 180
G = 256; RED = 256
DL = np.pi / G
lam_c = (np.arange(G) + 0.5) * np.pi / G
settings = [(0.0, 22.5), (0.0, 67.5), (45.0, 22.5), (45.0, 67.5)]
signs = [+1, -1, +1, +1]
S_QM = 2 * np.sqrt(2)
FLOOR = (np.sqrt(2.0) - 1.0) / 3.0
SGA = [np.sign(np.cos(2 * (lam_c - a * _deg))) for a, b in settings]
SGB = [np.sign(np.cos(2 * (lam_c + np.pi / 2 - b * _deg))) for a, b in settings]

REGISTER = os.environ.get("REGISTER", "paper")
NB_LIST = [int(x) for x in os.environ.get("NB_LIST", "16").split(",")]
EPS_LIST = [float(x) for x in os.environ.get("EPS_LIST", "0.1,0.05,0.02,0.01").split(",")]
NSTART = int(os.environ.get("NSTART", 6))
MAXIT = int(os.environ.get("MAXIT", 40))
PEN = float(os.environ.get("PEN", 100.0))
FEAS_TOL = float(os.environ.get("FEAS_TOL", 1e-6))
SEED = int(os.environ.get("SEED", 20261017))
JOINT_TIME_LIMIT = float(os.environ.get("JOINT_TIME_LIMIT", 1800))
SKIP_JOINT = int(os.environ.get("SKIP_JOINT", 0))


# ----------------------------------------------------------------------------------------------
def load(path):
    d = json.load(open(path))
    offs = np.asarray(d["offsets"], float); NV = len(offs)
    nang = int(d["nang"]); nw = int(d["nw"]); nsh = int(d["nsh"]); mode = d["mode"]; seed = int(d["seed"])
    cachef = os.path.join(os.path.dirname(path) or ".", f"banks_nsh{nsh}_{mode}{seed}_nang{nang}_nw{nw}.npz")
    if not os.path.exists(cachef):
        sys.exit(f"banks checkpoint not found: {cachef}")
    z = np.load(cachef)
    nflat = int(z["nflat"])
    Qs = [np.vstack([z[f"q{k}_{p}"] for k in range(nflat)]) for p in range(4)]   # (nflat*NV, G)
    Svec = np.concatenate([z[f"s{k}"] for k in range(nflat)])                      # (nflat*NV,)
    Ecoef = np.hstack([z[f"e{k}"] for k in range(nflat)])                          # (nang, nflat*NV)
    assert Qs[0].shape[0] == nflat * NV
    if REGISTER == "paper":
        angs = np.linspace(2, 178, nang)
    else:
        angs = np.linspace(6, 174, nang)
    return d, offs, nflat, Qs, Svec, Ecoef, angs


def md_of(Qs, rho):
    qm = [Qs[p].T @ rho for p in range(4)]
    return max(0.5 * np.sum(np.abs(qm[a] - qm[b])) * DL for a in range(4) for b in range(a + 1, 4))


def reduce_banks(Qs, Svec, Ecoef, grp, c, nth):
    """w_v = theta[grp[v]] * c[v]  ->  banks over theta."""
    Qr = []
    for p in range(4):
        R = np.zeros((nth, G)); np.add.at(R, grp, Qs[p] * c[:, None]); Qr.append(R)
    Sr = np.zeros(nth); np.add.at(Sr, grp, Svec * c)
    Er = np.zeros((Ecoef.shape[0], nth)); np.add.at(Er.T, grp, (Ecoef * c[None, :]).T)
    Mr = np.bincount(grp, weights=c, minlength=nth)      # mass of each theta: sum_v w_v = Mr . theta
    return Qr, Sr, Er, Mr


def solve_lp(Qs, Svec, Ecoef, Mr, angs, eps, penalty, time_limit=0.0):
    """min t  (+ PEN * slacks if penalty)  over theta >= 0.
       hard form (penalty=False): exactly min_md_full of dither_sweep_rho.py.
       Zero-mass columns (theta entries no atom maps to) are dropped and restored as 0."""
    keep = np.where(Mr > 0)[0]
    full_n = len(Mr)
    Qs = [Q[keep] for Q in Qs]; Svec = Svec[keep]; Ecoef = Ecoef[:, keep]; Mr = Mr[keep]
    NV = len(Svec); nang = len(angs)
    Qc = [Q.reshape(NV, RED, -1).mean(axis=2) for Q in Qs]
    cmb = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    n_s = len(cmb) * RED
    nu = (nang + 1) if penalty else 0
    nv = NV + 1 + n_s + nu
    rows, bub = [], []
    binw = np.pi / RED
    for ci, (i, j) in enumerate(cmb):
        Dif = (Qc[i] - Qc[j]) * binw
        for l in range(RED):
            r = np.zeros(nv); r[:NV] = Dif[:, l]; r[NV + 1 + ci * RED + l] = -1; rows.append(r); bub.append(0.0)
            r = np.zeros(nv); r[:NV] = -Dif[:, l]; r[NV + 1 + ci * RED + l] = -1; rows.append(r); bub.append(0.0)
    for ci in range(len(cmb)):
        r = np.zeros(nv); r[NV + 1 + ci * RED:NV + 1 + (ci + 1) * RED] = 0.5; r[NV] = -1; rows.append(r); bub.append(0.0)
    tgt = -np.cos(2 * angs * _deg)
    u0 = NV + 1 + n_s
    for ai in range(nang):
        r = np.zeros(nv); r[:NV] = Ecoef[ai]
        if penalty: r[u0 + ai] = -1
        rows.append(r); bub.append(tgt[ai] + eps)
        r = np.zeros(nv); r[:NV] = -Ecoef[ai]
        if penalty: r[u0 + ai] = -1
        rows.append(r); bub.append(-tgt[ai] + eps)
    Aeq = [np.zeros(nv)]; Aeq[0][:NV] = Mr; beq = [1.0]          # sum of the CENSUS weights = 1
    if penalty:
        r = np.zeros(nv); r[:NV] = Svec; r[u0 + nang] = -1; rows.append(r); bub.append(-S_QM)
        r = np.zeros(nv); r[:NV] = -Svec; r[u0 + nang] = -1; rows.append(r); bub.append(S_QM)
    else:
        r = np.zeros(nv); r[:NV] = Svec; Aeq.append(r); beq.append(-S_QM)
    Aub = np.array(rows); bv = np.array(bub)
    cobj = np.zeros(nv); cobj[NV] = 1
    if penalty: cobj[u0:] = PEN
    opts = {"time_limit": float(time_limit)} if time_limit and time_limit > 0 else None
    res = linprog(cobj, A_ub=Aub, b_ub=bv, A_eq=np.array(Aeq), b_eq=beq,
                  bounds=[(0, None)] * nv, method="highs", options=opts)
    if res.status != 0:
        return (None, None, "time limit") if res.status == 1 else (None, None, None)
    theta = np.zeros(full_n); theta[keep] = res.x[:NV]
    viol = float(res.x[u0:].max()) if penalty else 0.0
    return res.x[NV], theta, viol


# ----------------------------------------------------------------------------------------------
def run_joint(Qs, Svec, Ecoef, angs, eps, k_of, A_of, B_of, nflat, NB):
    grp = (k_of * NB + A_of) * NB + B_of
    nth = nflat * NB * NB
    Qr, Sr, Er, Mr = reduce_banks(Qs, Svec, Ecoef, grp, np.ones(len(grp)), nth)
    t, theta, flag = solve_lp(Qr, Sr, Er, Mr, angs, eps, penalty=False, time_limit=JOINT_TIME_LIMIT)
    if t is None: return ("time limit" if flag == "time limit" else None), None
    w = theta[grp]
    assert abs(w.sum() - 1) < 1e-8, w.sum()
    md = md_of(Qs, w)
    assert md >= FLOOR - 1e-9, f"Hall bound violated: {md}"
    return md, w


def alternate(Qs, Svec, Ecoef, angs, eps, blocks, init, maxit):
    """blocks: list of (name, grp_array, nth); w_v = prod_b theta_b[grp_b[v]].
       init: dict name -> theta.  Returns (md, w, viol, iters, history)."""
    th = {n: init[n].copy() for n, _, _ in blocks}
    hist = []; last = None
    for it in range(maxit):
        for name, grp, nth in blocks:
            c = np.ones(len(grp))
            for n2, g2, _ in blocks:
                if n2 != name: c *= th[n2][g2]
            Qr, Sr, Er, Mr = reduce_banks(Qs, Svec, Ecoef, grp, c, nth)
            t, theta, viol = solve_lp(Qr, Sr, Er, Mr, angs, eps, penalty=True)
            if t is None:
                return None
            th[name] = theta
        w = np.ones(len(blocks[0][1]))
        for n2, g2, _ in blocks: w *= th[n2][g2]
        assert abs(w.sum() - 1) < 1e-6, w.sum()
        md = md_of(Qs, w)
        if viol < FEAS_TOL: assert md >= FLOOR - 1e-7, f"Hall bound violated by a feasible census: {md}"
        hist.append((md, viol))
        if last is not None and abs(md - last[0]) < 1e-7 and abs(viol - last[1]) < 1e-9:
            break
        last = (md, viol)
    return md, w, viol, it + 1, hist


def main(path):
    t0 = time.time()
    d, offs, nflat, Qs, Svec, Ecoef, angs = load(path)
    NV = len(offs); NT = nflat * NV
    rho_star = np.asarray(d["rho"], float); rho_star = rho_star / rho_star.sum()
    print(f"file: {os.path.basename(path)}   atoms {NT} = {nflat} (tau,scale) x {NV} offsets   nang {len(angs)}  register [{angs[0]:.0f},{angs[-1]:.0f}]")
    md_star = md_of(Qs, rho_star)
    print(f"check: MD of the stored rho* on these banks = {md_star:.12f}   (file says {d['md']:.12f}; floor {FLOOR:.12f})")
    k_of = np.repeat(np.arange(nflat), NV)
    out = {"file": os.path.basename(path), "md_star": md_star, "floor": FLOOR, "results": []}
    rng = np.random.default_rng(SEED)
    for NB in NB_LIST:
        A_of = np.tile((offs[:, 0] / np.pi * NB).astype(int) % NB, nflat)
        B_of = np.tile((offs[:, 1] / np.pi * NB).astype(int) % NB, nflat)
        # starts: uniform; rho* marginals; random
        fA = np.bincount(A_of, weights=rho_star, minlength=NB); fB = np.bincount(B_of, weights=rho_star, minlength=NB)
        qk = np.bincount(k_of, weights=rho_star, minlength=nflat)
        Hk = np.zeros(nflat * NB); np.add.at(Hk, k_of * NB + A_of, rho_star)
        starts = [("uniform", dict(q=np.ones(nflat), f=np.ones(NB), g=np.ones(NB), h=np.ones(nflat * NB))),
                  ("rho*-marginals", dict(q=qk + 1e-9, f=fA + 1e-9, g=fB + 1e-9, h=Hk + 1e-9))]
        for s in range(NSTART):
            starts.append((f"random{s}", dict(q=rng.random(nflat), f=rng.random(NB), g=rng.random(NB), h=rng.random(nflat * NB))))
        for eps in EPS_LIST:
            print(f"\n=== NB={NB}  eps={eps} ===", flush=True)
            # joint (control)
            tj = time.time()
            if SKIP_JOINT:
                mdj = "skipped"
            else:
                mdj, wj = run_joint(Qs, Svec, Ecoef, angs, eps, k_of, A_of, B_of, nflat, NB)
            rec = {"NB": NB, "eps": eps, "joint": (mdj / FLOOR if isinstance(mdj, float) else mdj)}
            print(f"  joint   (binned, free)   : " + (f"MD/floor = {mdj/FLOOR:.4f}" if isinstance(mdj, float) else ("infeasible" if mdj is None else mdj)) + f"   [{time.time()-tj:.0f}s]", flush=True)
            for cls, blocks in (("semi", [("h", k_of * NB + A_of, nflat * NB), ("g", B_of, NB)]),
                                ("product", [("q", k_of, nflat), ("f", A_of, NB), ("g", B_of, NB)])):
                best = None; rows = []
                for sname, init in starts:
                    ta = time.time()
                    r = alternate(Qs, Svec, Ecoef, angs, eps, blocks, init, MAXIT)
                    if r is None:
                        rows.append((sname, None)); print(f"  {cls:8s} start {sname:15s}: LP failed"); continue
                    md, w, viol, iters, hist = r
                    feas = viol < FEAS_TOL
                    rows.append((sname, md / FLOOR, viol, iters, feas))
                    print(f"  {cls:8s} start {sname:15s}: MD/floor = {md/FLOOR:.4f}   max slack {viol:.2e} {'feasible' if feas else 'NOT feasible'}   {iters} rounds  [{time.time()-ta:.0f}s]", flush=True)
                    if feas and (best is None or md < best[0]): best = (md, w, sname)
                if best is None:
                    print(f"  {cls:8s} RESULT: no feasible {cls} census found at eps={eps}")
                    rec[cls] = None
                else:
                    md, w, sname = best
                    # Fourier check of the best census: |F(2,2)| vs product prediction
                    wo = np.zeros(NV); np.add.at(wo, np.tile(np.arange(NV), nflat), w); wo /= wo.sum()
                    xa, xb = offs[:, 0], offs[:, 1]
                    F22 = abs(np.sum(wo * np.exp(4j * (xa + xb)))); F20 = abs(np.sum(wo * np.exp(4j * xa))); F02 = abs(np.sum(wo * np.exp(4j * xb)))
                    print(f"  {cls:8s} RESULT: best feasible MD/floor = {md/FLOOR:.4f} (start {sname}); |F(2,2)| = {F22:.3f} vs product |F(2,0)||F(0,2)| = {F20*F02:.3f}")
                    rec[cls] = md / FLOOR; rec[cls + "_F22"] = F22; rec[cls + "_F20F02"] = F20 * F02
                rec[cls + "_starts"] = rows
            out["results"].append(rec)
            json.dump(out, open(os.path.join(os.path.dirname(path) or ".", f"product_census_v3_{os.path.splitext(os.path.basename(path))[0]}.json"), "w"), indent=1, default=float)
    print("\nSUMMARY (MD/floor; '-' = no feasible census in class)")
    print(f"{'NB':>4} {'eps':>6} {'joint':>8} {'semi':>8} {'product':>8}")
    for r in out["results"]:
        f = lambda x: "   -    " if x is None else (f"{x:>8s}" if isinstance(x, str) else f"{x:8.4f}")
        print(f"{r['NB']:>4} {r['eps']:>6} {f(r['joint'])} {f(r['semi'])} {f(r['product'])}")
    print(f"\nreference: rho* on these banks MD/floor = {md_star/FLOOR:.6f}    total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main(sys.argv[1])
