#!/usr/bin/env python3
"""slice_check.py — adversarial covariance test (finding F6 probe).
The census (frozen offsets xa, xb) breaks rotational covariance, and the
melody constraints pin E only on the a=0 slice. This script loads the
solved rho from the reduced C2 pattern.json, evolves NEW trajectories on
the a=45deg slice at NTEST relative angles, and evaluates the model's
E(45, 45+Delta) under the SAME rho. If |E - (-cos 2 Delta)| blows past
EPS on this unconstrained slice, faithfulness holds only on the
constrained slice — a scope caveat for the paper's "entire correlation
curve" sentence. Column ordering replicates extract_pattern exactly.
"""
import os, json
import numpy as np

pat = json.load(open("pattern.json"))
rho = np.array(pat["rho"]); NSH = pat["nsh"]; NW = pat["nw"]; EPS = pat["eps"]
_deg = np.pi/180; G = 256
TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]; NT = len(TAUS)
SCALES = [1.0, 1.8, 2.6]
lam_c = (np.arange(G)+0.5)*np.pi/G; DL = np.pi/G
gd = json.load(open("greedy_results.json"))
TERMS = [tuple(t) for t in gd["terms"]]
W0 = np.array(gd["w"]); V0 = float(gd["v"])
xi = np.arange(NSH)*np.pi/NSH
NVo = NSH*NSH
A_SLICE = float(os.environ.get("A_SLICE", 45.0))
NTEST = int(os.environ.get("NTEST", 7))
test_d = np.linspace(10, 170, NTEST)      # relative angles to probe

def force(lam, w, xa, xb, a, b):
    dU = np.zeros_like(lam)
    for wt, (m, n, sg) in zip(w, TERMS):
        if wt < 1e-12: continue
        k = m+n if sg > 0 else m-n
        ph = 2*m*(a+xa) + sg*2*n*(b+xb) + sg*n*np.pi
        dU += wt*2*k*np.sin(2*k*lam - ph)
    return V0 - dU

def evolve(w, xa, xb, a, b):
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

a0 = A_SLICE*_deg
sga = np.sign(np.cos(2*(lam_c - a0)))
print(f"probing UNCONSTRAINED slice a={A_SLICE} deg at {NTEST} relative angles")
print(f"{'Delta':>7} {'E_model':>9} {'E_QM':>9} {'resid':>9}  (band +-{EPS})")
worst = 0.0
done = {}
if os.path.exists("slice_partial.json"):
    done = {float(k): v for k, v in json.load(open("slice_partial.json")).items()}
for td in test_d:
    if float(td) in done:
        E = done[float(td)]
        eqm = -np.cos(2*td*_deg); r = E - eqm; worst = max(worst, abs(r))
        print(f"{td:7.1f} {E:+9.4f} {eqm:+9.4f} {r:+9.4f}  (cached)", flush=True)
        continue
    b0 = a0 + td*_deg
    ob = np.sign(np.cos(2*(lam_c + np.pi/2 - b0)))
    E = 0.0
    for si, sc in enumerate(SCALES):
        w = W0*sc
        j = 0
        for xa in xi:
            for xb in xi:
                o = evolve(w, xa, xb, a0, b0)
                for ti in range(NT):
                    col = si*NVo*NT + ti*NVo + j
                    E += rho[col]*np.sum(o[ti]*sga*ob)*DL
                j += 1
    eqm = -np.cos(2*td*_deg)
    r = E - eqm
    worst = max(worst, abs(r))
    flag = "" if abs(r) <= EPS + 1e-6 else "  <- OUTSIDE BAND"
    print(f"{td:7.1f} {E:+9.4f} {eqm:+9.4f} {r:+9.4f}{flag}", flush=True)
    done[float(td)] = float(E)
    json.dump({str(k): v for k, v in done.items()}, open("slice_partial.json", "w"))
print(f"\nworst |residual| on a={A_SLICE} slice: {worst:.4f} vs band {EPS}")
print("VERDICT:", "within band — covariance survives off-slice" if worst <= EPS + 1e-6
      else "VIOLATED off-slice — melody holds only on the constrained slice")
json.dump(dict(a_slice=A_SLICE, deltas=test_d.tolist(), worst=float(worst)),
          open("slice_check.json", "w"), indent=1)
