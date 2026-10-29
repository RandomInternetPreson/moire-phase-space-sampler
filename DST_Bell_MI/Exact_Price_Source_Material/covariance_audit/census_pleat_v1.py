#!/usr/bin/env python3
"""
census_pleat_v1.py -- a smooth census plus a PLEAT: does one high harmonic (a periodic ripple on the disc's hem) buy the floor?

Aaron's idea (28 Oct): the disc is not flat; its edge is pleated.  A pleat with N pleats around the disc modulates the
birth's offset as rho(xi) = smooth(xi) * (1 + a cos(2 N xi - phi)): three parameters, one HIGH harmonic -- exactly the
property every few-parameter form lacked (the floor needs structure at K ~ 14-18, i.e. 5-6 degrees).
This script is the certifiable LP version: the band-limited basis of degree K plus the harmonics a pleat at N contributes.
  PLEAT_MODE=lines   adds the pixels (N,0), (0,N), (N,N), (N,-N)          -- the pleat lines themselves (8 functions per cell)
  PLEAT_MODE=blocks  adds the whole K-block shifted by (N,0), (0,N), (N,N), (N,-N)  -- every sideband of pleat x smooth (a linear
                     superset of the product form, so its price is a lower bound on what a pleated smooth census pays)
PLEAT_N=0 reproduces census_bandlimit_v4 exactly.  Price vs N is a measurement of the pleat count on this terrain.
Compare (lattice, eps=0.01): K=6 alone 1.700; K=8 1.429; K=10 1.258; K=14 1.000.

--- census_bandlimit_v4 docstring follows ---
census_bandlimit_v4.py -- the smoothest census that attains the floor: minimal Fourier bandwidth of the birth.

v4 (22 Oct 2026) fixes one thing and adds two.
  FIX.  v1-v3 reported every non-optimal solver status except "time limit" as "infeasible".  HiGHS also
        returns a status for numerical failure, and the eps = 0.10 batch showed it: K = 6 and K = 8 came back
        "infeasible" at eps = 0.10 although the eps = 0.01 solutions are feasible points of the eps = 0.10 LP
        (eps enters only as the width of the correlator band; widening it can only enlarge the feasible set).
        v4 prints the solver's status and message for every LP and never calls anything "infeasible" on a
        status code alone.
  ADD 1. Every "infeasible" is now certified BY VALUE: a second LP with the same positivity, normalization
        and CHSH rows minimizes the extra width sigma that the correlator band would need; its optimum is
        the smallest miss any census of that bandwidth can achieve, and the miss is re-measured from the
        reconstructed census, independently of the solver.  sigma > 0 is the certificate; sigma = 0 means
        the price LP failed numerically, in which case the price LP is retried with the other HiGHS method.
        With MISS_ONLY=1 only this LP runs (the miss curve below K_min).
  ADD 2. Every feasible price is re-verified from the reconstructed census: rho >= 0, sum rho = 1,
        MD >= floor, and max |E(theta) - (-cos 2 theta)| <= eps, all printed.

Question / method are those of v1-v3: for each (tau, scale) cell the census over the offsets is a
trigonometric polynomial of degree K with free coefficients; non-negativity at every atom; normalization,
correlator band, CHSH and measurement-dependence rows exactly as in min_md_full (dither_sweep_rho.py).

Usage (covariance_audit folder; rho file next to its banks_*.npz):
  set K_LIST=2,3,4,5,6
  set EPS_LIST=0.01
  set METHOD=highs-ipm        primary solver;  METHOD2=highs-ds is the retry after a numerical failure
  set MISS_ONLY=0             1 = only the min-miss LP (certificates / miss curve), no price LP
  set TIME_LIMIT=0            seconds per LP, 0 = none;  TIME_LIMIT2 for the retry (default 7200)
  set OUT_TAG=<tag>           output json name tag (parallel runs on one rho file must not collide)
  python census_bandlimit_v4.py rho_nsh32_random0_eps0.01_nang161.json > bl_<tag>.txt
Writes census_bandlimit_v4_<tag>_<stem>.json with the coefficients of each feasible census.
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
TIME_LIMIT = float(os.environ.get("TIME_LIMIT", 0)); TIME_LIMIT2 = float(os.environ.get("TIME_LIMIT2", 7200))
METHOD = os.environ.get("METHOD", "highs-ipm"); METHOD2 = os.environ.get("METHOD2", "highs-ds")
MISS_ONLY = os.environ.get("MISS_ONLY", "0") == "1"
OUT_TAG = os.environ.get("OUT_TAG", "")
PLEAT_N = int(os.environ.get("PLEAT_N", 0)); PLEAT_MODE = os.environ.get("PLEAT_MODE", "lines")
STATUS = {0: "optimal", 1: "iteration/time limit", 2: "infeasible (solver status)", 3: "unbounded (solver status)", 4: "numerical difficulties"}

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

def pixel_set(K):
    """canonical (k_A, k_B) pairs: the K-block, plus the pleat pixels / shifted blocks for PLEAT_N > 0."""
    def canon(ka, kb):
        if ka < 0 or (ka == 0 and kb < 0): ka, kb = -ka, -kb
        return (ka, kb)
    pix = set(canon(ka, kb) for ka in range(-K, K + 1) for kb in range(-K, K + 1))
    if PLEAT_N > 0:
        N = PLEAT_N; shifts = [(N, 0), (0, N), (N, N), (N, -N)]
        if PLEAT_MODE == "lines":
            for sa, sb in shifts: pix.add(canon(sa, sb))
        else:
            for sa, sb in shifts:
                for ka in range(-K, K + 1):
                    for kb in range(-K, K + 1):
                        pix.add(canon(ka + sa, kb + sb))
    pix.discard((0, 0)); return sorted(pix)

def basis(offs, K):
    """real Fourier basis on the torus over the pixel set: columns = [1, cos/sin 2(kA xi_A + kB xi_B)]; returns Phi and labels."""
    xa, xb = offs[:, 0], offs[:, 1]; cols = [np.ones(len(offs))]; labels = [(0, 0, "1")]
    for ka, kb in pixel_set(K):
        arg = 2 * (ka * xa + kb * xb)
        cols.append(np.cos(arg)); labels.append((ka, kb, "c"))
        cols.append(np.sin(arg)); labels.append((ka, kb, "s"))
    return np.column_stack(cols), labels

def reduce(Qs, Svec, Ecoef, Phi, nflat, NV):
    """band-limited banks: everything in terms of the coefficients c (nflat*nb)."""
    nb = Phi.shape[1]
    Qr = [np.vstack([Phi.T @ Q[ts * NV:(ts + 1) * NV] for ts in range(nflat)]) for Q in Qs]        # (nc, G)
    Sr = np.concatenate([Phi.T @ Svec[ts * NV:(ts + 1) * NV] for ts in range(nflat)])                 # (nc,)
    Er = np.hstack([Ecoef[:, ts * NV:(ts + 1) * NV] @ Phi for ts in range(nflat)])                    # (nang, nc)
    Mr = np.concatenate([Phi.sum(axis=0) for _ in range(nflat)])                                      # (nc,)
    P = sp.block_diag([sp.csr_matrix(-Phi)] * nflat).tocsr()                                          # -rho <= 0
    return Qr, Sr, Er, Mr, P

def run_lp(cobj, Aub, bub, Aeq, beq, bounds, method, tlim):
    opts = {"time_limit": tlim} if tlim > 0 else None
    t = time.time()
    res = linprog(cobj, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds, method=method, options=opts)
    print(f"    [{method}] status {res.status} = {STATUS.get(res.status, '?')}; {res.message.strip()}   [{time.time()-t:.0f}s]", flush=True)
    return res

def price_lp(Qr, Sr, Er, Mr, P, angs, eps, nflat, nb, method, tlim):
    """min t  s.t.  MD rows, band, positivity, mass, CHSH.  Returns (res, nc)."""
    nc = nflat * nb; nang = len(angs)
    cmb = [(i, j) for i in range(4) for j in range(i + 1, 4)]; n_s = len(cmb) * RED
    nv = nc + 1 + n_s
    Qc = [Q.reshape(nc, RED, -1).mean(axis=2) for Q in Qr]
    binw = np.pi / RED; blocks = []; bub = []
    for ci, (i, j) in enumerate(cmb):
        Dif = (Qc[i] - Qc[j]) * binw
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
    blocks.append(sp.hstack([P, sp.csr_matrix((P.shape[0], 1 + n_s))]).tocsr()); bub += [0.0] * P.shape[0]
    Aub = sp.vstack(blocks).tocsr()
    Aeq = sp.csr_matrix(np.vstack([np.concatenate([Mr, [0], np.zeros(n_s)]), np.concatenate([Sr, [0], np.zeros(n_s)])]))
    cobj = np.zeros(nv); cobj[nc] = 1
    bounds = [(None, None)] * nc + [(0, None)] * (1 + n_s)
    return run_lp(cobj, Aub, np.array(bub), Aeq, [1.0, -S_QM], bounds, method, tlim), nc

def miss_lp(Sr, Er, Mr, P, angs, eps, nc, method, tlim):
    """min sigma  s.t.  |E - target| <= eps + sigma, positivity, mass, CHSH.  Always feasible (sigma large).
       sigma* > 0  <=>  no census of this bandwidth reproduces the singlet to eps: a certificate by value."""
    nang = len(angs); tgt = -np.cos(2 * angs * _deg)
    Ep = sp.hstack([sp.csr_matrix(Er), sp.csr_matrix(-np.ones((nang, 1)))]).tocsr()
    En = sp.hstack([sp.csr_matrix(-Er), sp.csr_matrix(-np.ones((nang, 1)))]).tocsr()
    Pp = sp.hstack([P, sp.csr_matrix((P.shape[0], 1))]).tocsr()
    Aub = sp.vstack([Ep, En, Pp]).tocsr(); bub = np.concatenate([tgt + eps, -tgt + eps, np.zeros(P.shape[0])])
    Aeq = sp.csr_matrix(np.vstack([np.concatenate([Mr, [0]]), np.concatenate([Sr, [0]])]))
    cobj = np.zeros(nc + 1); cobj[nc] = 1
    return run_lp(cobj, Aub, bub, Aeq, [1.0, -S_QM], [(None, None)] * nc + [(0, None)], method, tlim)

def verify(c, Phi, nflat, Qs, Er, angs, eps, band_required=True):
    """everything re-measured from the reconstructed census, independent of the solver's claims.
       band_required=False (min-miss census): report max |E - target| and judge only positivity, mass, floor."""
    rho = np.concatenate([Phi @ c[ts] for ts in range(nflat)])
    tgt = -np.cos(2 * angs * _deg); E = Er @ c.reshape(-1)
    v = dict(rho_min=float(rho.min()), mass=float(rho.sum()), band_miss=float(np.abs(E - tgt).max() - eps))
    rho_c = np.clip(rho, 0, None); rho_c /= rho_c.sum()
    v["md"] = float(md_of(Qs, rho_c)); v["price"] = v["md"] / FLOOR
    v["ok"] = (v["rho_min"] > -1e-7) and abs(v["mass"] - 1) < 1e-7 and v["md"] >= FLOOR - 1e-7 and (v["band_miss"] <= 1e-6 or not band_required)
    bandtxt = f"band miss over eps {v['band_miss']:+.2e}" if band_required else f"max |E - target| {v['band_miss']+eps:.6f}"
    print(f"    verify: rho_min {v['rho_min']:+.2e}  mass {v['mass']:.9f}  {bandtxt}  MD/floor {v['price']:.4f}   {'OK' if v['ok'] else '** VIOLATED **'}", flush=True)
    return rho_c, v

def spectrum(rho, offs, nflat, NV, K):
    wo = rho.reshape(nflat, NV).sum(axis=0); wo = wo / wo.sum(); xa, xb = offs[:, 0], offs[:, 1]; pix = []
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            F = np.sum(wo * np.exp(2j * (ka * xa + kb * xb))); pix.append((abs(F), ka, kb, np.degrees(np.angle(F))))
    pix.sort(reverse=True); return pix

def main(path):
    t0 = time.time()
    d, offs, nflat, Qs, Svec, Ecoef, angs = load(path); NV = len(offs); NT = nflat * NV
    rho_star = np.asarray(d["rho"], float); rho_star /= rho_star.sum(); md_star = md_of(Qs, rho_star)
    print(f"file: {os.path.basename(path)}   atoms {NT} = {nflat} x {NV}   K list {K_LIST}   eps {EPS_LIST}   solver {METHOD} (retry {METHOD2})   miss-only {MISS_ONLY}   PLEAT N={PLEAT_N} mode={PLEAT_MODE}")
    print(f"check: MD of the stored rho* on these banks = {md_star:.12f}   (file says {d['md']:.12f}; floor {FLOOR:.12f})")
    out = {"file": os.path.basename(path), "md_star": md_star, "results": []}
    jpath = os.path.join(os.path.dirname(path) or ".", f"census_pleat_v1_{OUT_TAG + '_' if OUT_TAG else ''}{os.path.splitext(os.path.basename(path))[0]}.json")
    for eps in EPS_LIST:
        for K in K_LIST:
            Phi, labels = basis(offs, K); nb = Phi.shape[1]; nc = nflat * nb
            print(f"\n=== eps={eps}  K={K}   ({nb} basis functions per (tau,scale) cell, {nc} coefficients) ===", flush=True)
            Qr, Sr, Er, Mr, P = reduce(Qs, Svec, Ecoef, Phi, nflat, NV)
            rec = {"eps": eps, "K": K}
            res = None
            if not MISS_ONLY:
                print("  price LP:", flush=True)
                res, _ = price_lp(Qr, Sr, Er, Mr, P, angs, eps, nflat, nb, METHOD, TIME_LIMIT)
            if res is not None and res.status == 0:
                c = res.x[:nc].reshape(nflat, nb); rho, v = verify(c, Phi, nflat, Qs, Er, angs, eps)
                rec.update(status="feasible" if v["ok"] else "feasible-unverified", price=v["price"], verify=v, solver=METHOD)
            else:
                print("  min-miss LP:", flush=True)
                rm = miss_lp(Sr, Er, Mr, P, angs, eps, nc, METHOD, TIME_LIMIT)
                if rm.status != 0:
                    rm = miss_lp(Sr, Er, Mr, P, angs, eps, nc, METHOD2, TIME_LIMIT2)
                if rm.status != 0:
                    rec.update(status="undetermined (both solvers failed on the min-miss LP)"); print(f"  {rec['status']}", flush=True)
                else:
                    c = rm.x[:nc].reshape(nflat, nb); rho, v = verify(c, Phi, nflat, Qs, Er, angs, 0.0, band_required=False)   # eps=0: band_miss = max|E-target|
                    miss = v["band_miss"]; sigma = float(rm.x[nc])
                    print(f"    min-miss: sigma* = {sigma:.6f}  ->  best K={K} census misses the singlet by {miss:.6f} (= eps {eps} + {miss-eps:+.6f})", flush=True)
                    rec.update(sigma=sigma, min_miss=miss, miss_verify=v)
                    if miss > eps + 1e-6:
                        rec.update(status="infeasible (certified by value)")
                        print(f"  INFEASIBLE at eps={eps}: certified, min miss {miss:.5f} > eps", flush=True)
                    else:
                        print(f"  the band IS attainable at K={K} (min miss {miss:.2e} <= eps)" + ("" if MISS_ONLY else f"; the price LP failed numerically -> retry with {METHOD2}"), flush=True)
                        rec.update(upper_bound_price=v["price"])   # the min-miss census is feasible: its MD is an upper bound on the price
                        if MISS_ONLY:
                            rec.update(status="feasible (miss-only run; price not computed)")
                        else:
                            res2, _ = price_lp(Qr, Sr, Er, Mr, P, angs, eps, nflat, nb, METHOD2, TIME_LIMIT2)
                            if res2.status == 0:
                                c = res2.x[:nc].reshape(nflat, nb); rho, v = verify(c, Phi, nflat, Qs, Er, angs, eps)
                                rec.update(status="feasible" if v["ok"] else "feasible-unverified", price=v["price"], verify=v, solver=METHOD2)
                            else:
                                rec.update(status=f"feasible, price undetermined (both solvers failed); upper bound {v['price']:.4f} from the min-miss census")
                                print(f"  {rec['status']}", flush=True)
            if "price" in rec:
                pix = spectrum(rho, offs, nflat, NV, max(K, PLEAT_N + K)); nact = int(np.sum(np.abs(c) > 1e-6 * np.abs(c).max()))
                print(f"  MD/floor = {rec['price']:.4f}   {rec['status']}   active coefficients {nact}/{nc}   [{time.time()-t0:.0f}s total]")
                print("  leading harmonics of the pooled census  (|F|, k_A, k_B, phi[deg])  -- (k_S,k_D) = ((k_A+k_B)/2, (k_A-k_B)/2):")
                for a, ka, kb, ph in pix[:10]:
                    print(f"     {a:.3f}  ({ka:+d},{kb:+d})  {ph:+8.2f}    (S,D)-harmonic ({(ka+kb)/2:+.1f},{(ka-kb)/2:+.1f})")
                rec.update(active=nact, leading=[(float(a), int(ka), int(kb), float(ph)) for a, ka, kb, ph in pix[:20]],
                           coeff_labels=[list(map(str, l)) for l in labels], coeff=c.tolist())
            out["results"].append(rec)
            json.dump(out, open(jpath, "w"), default=float)
    print("\nSUMMARY  (MD/floor of the smoothest census of bandwidth K; infeasible entries carry the certified minimum miss)")
    print(f"{'eps':>6} {'K':>3} {'price':>9}   note")
    for r in out["results"]:
        if "price" in r: print(f"{r['eps']:>6} {r['K']:>3} {r['price']:9.4f}   {r['status']} [{r.get('solver','')}]")
        elif "min_miss" in r: print(f"{r['eps']:>6} {r['K']:>3} {'-':>9}   {r['status']}; min miss {r['min_miss']:.5f}" + (f"; price <= {r['upper_bound_price']:.4f}" if 'upper_bound_price' in r else ""))
        else: print(f"{r['eps']:>6} {r['K']:>3} {'-':>9}   {r['status']}")
    print(f"\nreference: rho* (unrestricted, {NV} atoms per cell) MD/floor = {md_star/FLOOR:.6f}   total {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main(sys.argv[1])
