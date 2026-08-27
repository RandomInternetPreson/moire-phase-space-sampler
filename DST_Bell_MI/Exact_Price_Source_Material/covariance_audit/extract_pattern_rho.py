#!/usr/bin/env python3
"""
extract_pattern.py — THE PREDICTION.
Re-solves the pooled-census LP at the empirical register (EPS, default
0.01), then evaluates the optimal mixture's correlation curve on a FINE
angle grid and reports the residual pattern E_model(Delta) - E_QM(Delta):
where the mechanism must press against quantum mechanics, with what sign,
at which angles. Output: pattern table + pattern.png + pattern.json.
This is simultaneously the paper's prediction section and the pure-MD
stake for R3: "if nature runs moire phase locking at this register,
deviations concentrate HERE, shaped like THIS, at the <=1% level."
Run in breedrun: python extract_pattern.py   (options: EPS=0.01 NSH=12
NANGF=41)
"""
import os, json, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.optimize import linprog

_deg = np.pi/180; G = 256; RED = 256; NANG = 13
NSH = int(os.environ.get("NSH", 12))
NW = int(os.environ.get("NW", 2048))
NANGF = int(os.environ.get("NANGF", 41))
EPS = float(os.environ.get("EPS", 0.01))
DENSE = int(os.environ.get("DENSE", 1))
TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]
settings = [(0.0,22.5),(0.0,67.5),(45.0,22.5),(45.0,67.5)]
signs = [+1,-1,+1,+1]
S_QM = 2*np.sqrt(2); FLOOR = 0.138071
lam_c = (np.arange(G)+0.5)*np.pi/G; DL = np.pi/G
SGA = [np.sign(np.cos(2*(lam_c-a*_deg))) for a,b in settings]
SGB = [np.sign(np.cos(2*(lam_c+np.pi/2-b*_deg))) for a,b in settings]

gd = json.load(open("greedy_results.json"))
TERMS = [tuple(t) for t in gd["terms"]]
W0 = np.array(gd["w"]); V0 = float(gd["v"])
SCALES = [1.0, 1.8, 2.6]
xi = np.arange(NSH)*np.pi/NSH
angC = np.linspace(6, 174, NANG)          # constraint angles
angF = np.linspace(2, 178, NANGF)         # fine evaluation grid

def force(lam, w, xa, xb, a, b):
    dU = np.zeros_like(lam)
    for wt,(m,n,sg) in zip(w, TERMS):
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
    cores = int(os.environ.get("CORES", os.cpu_count() or 8))
    NVo = NSH*NSH; NT = len(TAUS); NV = NVo*NT*len(SCALES)
    jobs, meta = [], []
    for si, sc in enumerate(SCALES):
        w = W0*sc
        for pi,(a,b) in enumerate(settings):
            for xa in xi:
                for xb in xi:
                    jobs.append((w, xa, xb, a*_deg, b*_deg))
                    meta.append(("p", si, pi))
        for ai, td in enumerate(angF):
            for xa in xi:
                for xb in xi:
                    jobs.append((w, xa, xb, 0.0, td*_deg))
                    meta.append(("f", si, ai))
    print(f"evolving {len(jobs)} trajectories (constraint pairs + "
          f"{NANGF} fine angles)...")
    t0 = time.time()
    with Pool(cores) as pool:
        outs = pool.map(evolve, jobs, chunksize=4)
    print(f"  done in {time.time()-t0:.0f}s")
    Qs = [np.zeros((NV, G)) for _ in range(4)]
    Ef = np.zeros((NANGF, NV))
    cp = {}; cf = {}
    for (kind, si, idx), o in zip(meta, outs):
        for ti in range(NT):
            col_base = si*NVo*NT + ti*NVo
            if kind == "p":
                j = cp.setdefault((si, idx, ti), 0)
                Qs[idx][col_base + j % NVo] = o[ti]
                cp[(si, idx, ti)] += 1
            else:
                j = cf.setdefault((si, idx, ti), 0)
                tr = angF[idx]*_deg
                ob = np.sign(np.cos(2*(lam_c+np.pi/2-tr)))
                Ef[idx, col_base + j % NVo] = np.sum(o[ti]*SGA[0]*ob)*DL
                cf[(si, idx, ti)] += 1
    Svec = np.zeros(NV)
    for i, sg in enumerate(signs):
        Svec += sg*(Qs[i]*(SGA[i]*SGB[i])[None,:]).sum(axis=1)*DL
    # constraints: DENSE=1 (default) pins EVERY fine angle - uniform
    # fidelity; DENSE=0 reproduces the old 13-sample register.
    if DENSE:
        idxC = list(range(NANGF))
    else:
        idxC = [int(np.argmin(np.abs(angF - a))) for a in angC]
    cmb = [(i,j) for i in range(4) for j in range(i+1,4)]
    nv = NV + 1 + len(cmb)*RED
    rows, bub = [], []
    binw = np.pi/RED
    for ci,(i,j) in enumerate(cmb):
        Dif = (Qs[i]-Qs[j])*binw
        for l in range(RED):
            r = np.zeros(nv); r[:NV] = Dif[:,l]; r[NV+1+ci*RED+l] = -1
            rows.append(r); bub.append(0.0)
            r = np.zeros(nv); r[:NV] = -Dif[:,l]; r[NV+1+ci*RED+l] = -1
            rows.append(r); bub.append(0.0)
    for ci in range(len(cmb)):
        r = np.zeros(nv); r[NV+1+ci*RED:NV+1+(ci+1)*RED] = 0.5
        r[NV] = -1; rows.append(r); bub.append(0.0)
    for k in idxC:
        tgt = -np.cos(2*angF[k]*_deg)
        r = np.zeros(nv); r[:NV] = Ef[k]; rows.append(r); bub.append(tgt+EPS)
        r = np.zeros(nv); r[:NV] = -Ef[k]; rows.append(r); bub.append(-tgt+EPS)
    Aub = np.array(rows); bv = np.array(bub)
    Aeq = np.zeros((2,nv)); Aeq[0,:NV] = 1; Aeq[1,:NV] = Svec
    cobj = np.zeros(nv); cobj[NV] = 1
    res = linprog(cobj, A_ub=Aub, b_ub=bv, A_eq=Aeq, b_eq=[1.0,-S_QM],
                  bounds=[(0,None)]*nv, method="highs")
    if res.status != 0:
        raise SystemExit("LP infeasible - raise EPS")
    rho = res.x[:NV]
    qm = [Qs[p].T@rho for p in range(4)]
    md = max(0.5*np.sum(np.abs(qm[a]-qm[b]))*DL
             for a in range(4) for b in range(a+1,4))
    resid = (Ef@rho) - (-np.cos(2*angF*_deg))
    print(f"\nsolution: MD = {md:.5f} ({md/FLOOR:.3f}x floor) at dev<={EPS} "
          f"({'UNIFORM ' + str(NANGF) + '-angle' if DENSE else '13-sample'} register)")
    print(f"\nTHE PREDICTED DEVIATION PATTERN (E_model - E_QM):")
    print(f"{'Delta':>7} {'residual':>10}  note")
    act = []
    for k, t in enumerate(angF):
        note = ""
        if abs(abs(resid[k]) - EPS) < 1e-4 and k in idxC:
            note = "<- ACTIVE (pressed against the bound)"
            act.append(round(float(t),1))
        if abs(resid[k]) > EPS*0.5 or note:
            print(f"{t:7.1f} {resid[k]:+10.4f}  {note}")
    print(f"\nactive angles: {act}")
    json.dump(dict(eps=EPS, md=float(md), angles=angF.tolist(),
                   residual=resid.tolist(), active=act,
                   rho=rho.tolist(), nsh=NSH, nw=NW),
              open("pattern.json","w"), indent=1)
    fig, ax = plt.subplots(figsize=(11, 4.6))
    ax.axhspan(-EPS, EPS, color="#1D9E75", alpha=0.12,
               label=f"allowed band ±{EPS}")
    ax.plot(angF, resid, color="#D85A30", lw=2)
    for a, b in settings:
        ax.axvline(abs(a-b), color="#888", ls=":", lw=0.8)
    ax.axhline(0, color="#444", lw=0.7)
    ax.set_xlabel("relative angle Δ (deg)")
    ax.set_ylabel("E_model − E_QM")
    ax.set_title(f"Moiré phase locking: predicted deviation pattern "
                 f"(MD = {md/FLOOR:.3f}× floor, register {EPS})", fontsize=11)
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig("pattern.png", dpi=150)
    print("saved pattern.png + pattern.json")
    print("\nTHE STAKE: if nature runs this mechanism at this register,")
    print("correlation deviations concentrate at the angles above, with")
    print("the printed signs, at the <=1% level - measurable by any")
    print("modern Bell apparatus willing to scan between the Bell angles.")
