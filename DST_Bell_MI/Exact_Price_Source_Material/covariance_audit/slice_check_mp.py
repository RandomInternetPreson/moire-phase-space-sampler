#!/usr/bin/env python3
"""slice_check_mp.py — parallel covariance probe (math of slice_check.py,
byte-identical column ordering; reads pattern_rho.json written by
slice_solve.py to avoid the pattern.json schema collision).
Evolves NEW trajectories on slice a=A_SLICE at NTEST relative angles and
evaluates E under the frozen rho.  Env: A_SLICE=45 NTEST=7 CORES=n.
"""
import os, json
import numpy as np
from multiprocessing import Pool

PAT = os.environ.get("PAT", "pattern_rho.json")
pat = json.load(open(PAT))
rho = np.array(pat["rho"]); NSH = pat["nsh"]; NW = pat["nw"]; EPS = pat.get("eps", 0.01)
_deg = np.pi/180; G = 256
TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]; NT = len(TAUS)
SCALES = [1.0, 1.8, 2.6]
lam_c = (np.arange(G)+0.5)*np.pi/G; DL = np.pi/G
gd = json.load(open("greedy_results.json"))
TERMS = [tuple(t) for t in gd["terms"]]
W0 = np.array(gd["w"]); V0 = float(gd["v"])
if "offsets" in pat:
    OFFS = [tuple(o) for o in pat["offsets"]]
else:
    _xi = np.arange(NSH)*np.pi/NSH
    OFFS = [(xa, xb) for xa in _xi for xb in _xi]
NVo = len(OFFS)
ORDER = pat.get("ordering", "scale_major")   # scale_major: si*NT*NVo+ti*NVo+j ; tau_major: ti*NSC*NVo+si*NVo+j
A_SLICE = float(os.environ.get("A_SLICE", 45.0))
NTEST = int(os.environ.get("NTEST", 7))
CORES = int(os.environ.get("CORES", os.cpu_count() or 8))
DMIN = float(os.environ.get("DMIN", 10.0))
DMAX = float(os.environ.get("DMAX", 170.0))
test_d = np.linspace(DMIN, DMAX, NTEST)

def force(lam, w, xa, xb, a, b):
    dU = np.zeros_like(lam)
    for wt, (m, n, sg) in zip(w, TERMS):
        if wt < 1e-12: continue
        k = m+n if sg > 0 else m-n
        ph = 2*m*(a+xa) + sg*2*n*(b+xb) + sg*n*np.pi
        dU += wt*2*k*np.sin(2*k*lam - ph)
    return V0 - dU

def evolve(args):
    w, xa, xb, a, b = args
    lam = np.linspace(0, np.pi, NW, endpoint=False).copy()
    dt = 0.004; t = 0.0; out = []
    for tau in TAUS:
        for _ in range(int(round((tau-t)/dt))):
            lam += force(lam, w, xa, xb, a, b)*dt
        t = tau
        h = np.bincount(((lam % np.pi)/np.pi*G).astype(int) % G,
                        minlength=G).astype(float)
        out.append(h/(h.sum()*DL))
    return out

if __name__ == "__main__":
    a0 = A_SLICE*_deg
    sga = np.sign(np.cos(2*(lam_c - a0)))
    print(f"probing slice a={A_SLICE} deg at {NTEST} relative angles "
          f"(two_slice register: {pat.get('two_slice', 0)})")
    print(f"{'Delta':>7} {'E_model':>9} {'E_QM':>9} {'resid':>9}  (band +-{EPS})")
    worst = 0.0; results = {}
    for td in test_d:
        b0 = a0 + td*_deg
        ob = np.sign(np.cos(2*(lam_c + np.pi/2 - b0)))
        jobs = [(W0*sc, xa, xb, a0, b0)
                for sc in SCALES for (xa, xb) in OFFS]
        with Pool(CORES) as pool:
            outs = pool.map(evolve, jobs, chunksize=4)
        E = 0.0; idx = 0; NSC = len(SCALES)
        for si in range(NSC):
            for j in range(NVo):
                o = outs[idx]; idx += 1
                for ti in range(NT):
                    col = (si*NVo*NT + ti*NVo + j) if ORDER == "scale_major" \
                          else (ti*NSC*NVo + si*NVo + j)
                    E += rho[col]*np.sum(o[ti]*sga*ob)*DL
        eqm = -np.cos(2*td*_deg); r = E - eqm
        worst = max(worst, abs(r)); results[float(td)] = float(E)
        flag = "" if abs(r) <= EPS + 1e-6 else "  <- OUTSIDE BAND"
        print(f"{td:7.1f} {E:+9.4f} {eqm:+9.4f} {r:+9.4f}{flag}", flush=True)
    print(f"\nworst |residual| on a={A_SLICE} slice: {worst:.4f} vs band {EPS}")
    print("VERDICT:", "within band" if worst <= EPS + 1e-6
          else "OUTSIDE band on this slice")
    json.dump(dict(a_slice=A_SLICE, worst=float(worst), results=results,
                   two_slice=pat.get("two_slice", 0), pat=PAT),
              open(f"slice_check_a{A_SLICE:g}_d{DMIN:g}-{DMAX:g}.json", "w"), indent=1)
