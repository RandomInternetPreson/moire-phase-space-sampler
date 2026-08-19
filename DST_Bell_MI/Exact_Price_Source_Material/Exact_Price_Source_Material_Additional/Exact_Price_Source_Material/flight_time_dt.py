#!/usr/bin/env python3
"""
flight_time.py — Aaron's distance-traveled hypothesis, built properly.
Components are no longer infinite-time stationary laws: a pair detected
after finite flight time tau carries the PARTIALLY RELAXED distribution
— the deterministic (T=0) flow applied to uniform birth fuzz for time
tau. This interpolates uniform -> locked, a strictly richer family, and
tau is physical (source-analyzer geometry). Two questions:
  (1) per-tau: is there a sweet flight time (accumulation without
      over-relaxation)?
  (2) pooled: does a SPREAD of flight times (tau as one more quenched
      census coordinate — real labs have path-length distributions)
      beat any single tau?
Benchmarks: tau=infinity (residence law) = 1.19-1.22x at this grid;
melody register dev<=0.10 throughout.
Falsifier on the books: any tau-dependence must saturate fast
(Weihs/Innsbruck at 400 m saw pure QM).
Run: CORES=56 python flight_time.py     (~20-40 min)
Options: NSH=12 NW=2048 TAUS="0.5,1,2,4,8"
"""
import os, json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linprog

_deg = np.pi/180; G = 256; RED = 256
NANG = int(os.environ.get('NANG', 13))
NSH = int(os.environ.get("NSH", 12))
NW = int(os.environ.get("NW", 2048))
TAUS = [float(x) for x in os.environ.get("TAUS","0.5,1,2,4,8").split(",")]
EPS = float(os.environ.get("EPS", 0.10))
settings = [(0.0,22.5),(0.0,67.5),(45.0,22.5),(45.0,67.5)]
signs = [+1,-1,+1,+1]
S_QM = 2*np.sqrt(2); FLOOR = 0.138071
lam_c = (np.arange(G)+0.5)*np.pi/G; DL = np.pi/G
SGA = [np.sign(np.cos(2*(lam_c-a*_deg))) for a,b in settings]
SGB = [np.sign(np.cos(2*(lam_c+np.pi/2-b*_deg))) for a,b in settings]

gd = json.load(open("greedy_results.json"))
TERMS = [tuple(t) for t in gd["terms"]]
W0 = np.array(gd["w"]); V0 = float(gd["v"])
SCALES = [float(x) for x in os.environ.get("SCALES","1.0,1.8,2.6").split(",")]
xi = np.arange(NSH)*np.pi/NSH
angs = np.linspace(6, 174, NANG)

def force(lam, w, xa, xb, a, b):
    dU = np.zeros_like(lam)
    for wt,(m,n,sg) in zip(w, TERMS):
        if wt < 1e-12: continue
        k = m+n if sg > 0 else m-n
        ph = 2*m*(a+xa) + sg*2*n*(b+xb) + sg*n*np.pi
        dU += wt*2*k*np.sin(2*k*lam - ph)
    return V0 - dU

def evolve_component(args):
    """returns histograms at each tau snapshot for one (scale, xa, xb, pair)"""
    w, xa, xb, a, b = args
    lam = np.linspace(0, np.pi, NW, endpoint=False).copy()
    dt = float(os.environ.get("DT", 0.004))
    out = []
    t = 0.0
    for tau in TAUS:
        steps = int(round((tau - t)/dt))
        for _ in range(steps):
            lam += force(lam, w, xa, xb, a, b)*dt
        t = tau
        h = np.bincount(((lam % np.pi)/np.pi*G).astype(int) % G,
                        minlength=G).astype(float)
        out.append(h/(h.sum()*DL))
    return out

def build_banks(cores):
    """BANKS[tau_idx][scale_idx] = (Qs[4][NV,G], Svec, Ecoef)"""
    NV = NSH*NSH
    jobs, meta = [], []
    for si, sc in enumerate(SCALES):
        w = W0*sc
        for pi,(a,b) in enumerate(settings):
            for xa in xi:
                for xb in xi:
                    jobs.append((w, xa, xb, a*_deg, b*_deg))
                    meta.append(("pair", si, pi))
        for ai, tdeg in enumerate(angs):
            for xa in xi:
                for xb in xi:
                    jobs.append((w, xa, xb, 0.0, tdeg*_deg))
                    meta.append(("ang", si, ai))
    print(f"evolving {len(jobs)} component trajectories "
          f"({len(TAUS)} tau snapshots each)...")
    t0 = time.time()
    with Pool(cores) as pool:
        outs = pool.map(evolve_component, jobs, chunksize=4)
    print(f"  done in {time.time()-t0:.0f}s")
    BANKS = [[None]*len(SCALES) for _ in TAUS]
    for ti in range(len(TAUS)):
        for si in range(len(SCALES)):
            Qs = [np.zeros((NV, G)) for _ in range(4)]
            Ec = np.zeros((NANG, NV))
            cnt_p = [0]*4; cnt_a = [0]*NANG
            for (kind, s2, idx), o in zip(meta, outs):
                if s2 != si: continue
                if kind == "pair":
                    Qs[idx][cnt_p[idx]] = o[ti]; cnt_p[idx] += 1
                else:
                    tr = angs[idx]*_deg
                    ob = np.sign(np.cos(2*(lam_c+np.pi/2-tr)))
                    Ec[idx, cnt_a[idx]] = np.sum(o[ti]*SGA[0]*ob)*DL
                    cnt_a[idx] += 1
            Svec = np.zeros(NV)
            for i, sg in enumerate(signs):
                Svec += sg*(Qs[i]*(SGA[i]*SGB[i])[None,:]).sum(axis=1)*DL
            BANKS[ti][si] = (Qs, Svec, Ec)
    return BANKS

def min_md(banklist, eps=EPS):
    """banklist: list of (Qs, Svec, Ecoef) to POOL as one census."""
    Qs = [np.vstack([b[0][p] for b in banklist]) for p in range(4)]
    Svec = np.concatenate([b[1] for b in banklist])
    Ecoef = np.hstack([b[2] for b in banklist])
    NV = len(Svec)
    Qc = [Q.reshape(NV, RED, -1).mean(axis=2) for Q in Qs]
    cmb = [(i,j) for i in range(4) for j in range(i+1,4)]
    nv = NV + 1 + len(cmb)*RED
    rows, bub = [], []
    binw = np.pi/RED
    for ci,(i,j) in enumerate(cmb):
        Dif = (Qc[i]-Qc[j])*binw
        for l in range(RED):
            r = np.zeros(nv); r[:NV] = Dif[:,l]; r[NV+1+ci*RED+l] = -1
            rows.append(r); bub.append(0.0)
            r = np.zeros(nv); r[:NV] = -Dif[:,l]; r[NV+1+ci*RED+l] = -1
            rows.append(r); bub.append(0.0)
    for ci in range(len(cmb)):
        r = np.zeros(nv); r[NV+1+ci*RED:NV+1+(ci+1)*RED] = 0.5
        r[NV] = -1; rows.append(r); bub.append(0.0)
    tgt = -np.cos(2*angs*_deg)
    for ai in range(NANG):
        r = np.zeros(nv); r[:NV] = Ecoef[ai]; rows.append(r)
        bub.append(tgt[ai]+eps)
        r = np.zeros(nv); r[:NV] = -Ecoef[ai]; rows.append(r)
        bub.append(-tgt[ai]+eps)
    Aub = np.array(rows); bv = np.array(bub)
    Aeq = np.zeros((2,nv)); Aeq[0,:NV] = 1; Aeq[1,:NV] = Svec
    cobj = np.zeros(nv); cobj[NV] = 1
    res = linprog(cobj, A_ub=Aub, b_ub=bv, A_eq=Aeq, b_eq=[1.0,-S_QM],
                  bounds=[(0,None)]*nv, method="highs")
    if res.status != 0: return None
    rho = res.x[:NV]
    qm = [Qs[p].T@rho for p in range(4)]
    return max(0.5*np.sum(np.abs(qm[a]-qm[b]))*DL
               for a in range(4) for b in range(a+1,4))

if __name__ == "__main__":
    cores = int(os.environ.get("CORES", os.cpu_count() or 8))
    BANKS = build_banks(cores)
    print(f"\nFLIGHT-TIME MAP (melody dev<={EPS}; benchmark register: dev 0.10)")
    print(f"{'tau':>6} {'MD/floor':>9}")
    results = {}
    for ti, tau in enumerate(TAUS):
        md = min_md(BANKS[ti])
        results[f"tau{tau}"] = md
        print(f"{tau:6.1f} " + (f"{md/FLOOR:9.3f}" if md else "  infeas"))
    # pooled: tau as a census coordinate (path-length spread)
    md = min_md([BANKS[ti][si] for ti in range(len(TAUS))
                 for si in range(len(SCALES))][:0] or
                [b for row in BANKS for b in row])
    results["pooled_all_tau"] = md
    print(f"{'POOLED':>6} " + (f"{md/FLOOR:9.3f}" if md else "  infeas") +
          "   <- tau as a quenched census coordinate")
    json.dump(results, open("flight_time_results.json","w"), indent=1)
    print("\nREADING: a sweet single tau below the tau=inf benchmark =>")
    print("partial relaxation is a resource (accumulation without")
    print("over-relaxation). POOLED below every single tau => the path-")
    print("length SPREAD itself is doing work — distance traveled joins")
    print("the census. Saturation check: if tau=4 ~ tau=8 ~ tau=inf, the")
    print("Weihs 400 m falsifier is automatically satisfied.")
