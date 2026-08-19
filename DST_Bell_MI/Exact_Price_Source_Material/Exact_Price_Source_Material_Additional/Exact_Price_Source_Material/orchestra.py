#!/usr/bin/env python3
"""
orchestra.py — THE FULL ORCHESTRA.
Every instrument at once: the pooled census (frozen offsets x flight
times x scales) PLUS the solved sliding-detection law D(m) = sqrt(m)h(m)
at strength alpha (the exact-cosine channel). With detection, correlators
are ratios, so constraints are rate-normalized (linear in rho) plus the
physical flat-coincidence-rate condition. Scored against the CERTIFIED
SURFACE: at effective efficiency eta < 1 the floor drops below 0.13807,
so the scoreboard is MD versus M(2*sqrt2, eta_eff) — the certified
minimum for that efficiency.
Question: does detection, carrying melody natively, pull the uniform-1%
price toward the surface? Run: python orchestra.py
Options: NSH=12 EPS=0.01 ALPHAS="0.5,0.75,1.0" NANGF=41
"""
import os, json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linprog

_deg = np.pi/180; G = 256; RED = 256
NSH = int(os.environ.get("NSH", 12))
NW = int(os.environ.get("NW", 2048))
NANGF = int(os.environ.get("NANGF", 41))
EPS = float(os.environ.get("EPS", 0.01))
ALPHAS = [float(x) for x in os.environ.get("ALPHAS","0.5,0.75,1.0").split(",")]
TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]
settings = [(0.0,22.5),(0.0,67.5),(45.0,22.5),(45.0,67.5)]
signs = [+1,-1,+1,+1]
S_QM = 2*np.sqrt(2)
lam_c = (np.arange(G)+0.5)*np.pi/G; DL = np.pi/G

# solved detection profile D(m) = sqrt(m) h(m), h normalized to h(1)=1
H_M = np.array([0.0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
                0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0])
H_V = np.array([1.17451,1.17372,1.17106,1.16496,1.15789,1.15012,1.14180,
                1.13300,1.12377,1.11414,1.10414,1.09379,1.08311,1.07211,
                1.06081,1.04922,1.03735,1.02521,1.01281,1.00016,1.00008,1.0])
def Dstar(m):
    m = np.clip(m, 0, 1)
    return np.sqrt(m)*np.interp(m, H_M, H_V)

# certified surface floor at eta: closed form M(2sqrt2, eta)
def floor_at(eta):
    return max(0.0, eta*((1+np.sqrt(2))*eta - 2)/3)

gd = json.load(open("greedy_results.json"))
TERMS = [tuple(t) for t in gd["terms"]]
W0 = np.array(gd["w"]); V0 = float(gd["v"])
SCALES = [1.0, 1.8, 2.6]
xi = np.arange(NSH)*np.pi/NSH
angF = np.linspace(2, 178, NANGF)

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
                    jobs.append((w, xa, xb, a*_deg, b*_deg)); meta.append(("p",si,pi))
        for ai, td in enumerate(angF):
            for xa in xi:
                for xb in xi:
                    jobs.append((w, xa, xb, 0.0, td*_deg)); meta.append(("f",si,ai))
    print(f"evolving {len(jobs)} trajectories once (reused for all alphas)...")
    t0 = time.time()
    with Pool(cores) as pool:
        outs = pool.map(evolve, jobs, chunksize=4)
    print(f"  done in {time.time()-t0:.0f}s")

    # raw source distributions per column (detection applied per alpha below)
    Qs_raw = [np.zeros((NV, G)) for _ in range(4)]
    Qf_raw = [np.zeros((NV, G)) for _ in range(NANGF)]
    cp, cf = {}, {}
    for (kind, si, idx), o in zip(meta, outs):
        for ti in range(NT):
            col = si*NVo*NT + ti*NVo
            if kind == "p":
                j = cp.setdefault((si, idx, ti), 0)
                Qs_raw[idx][col + j % NVo] = o[ti]; cp[(si,idx,ti)] += 1
            else:
                j = cf.setdefault((si, idx, ti), 0)
                Qf_raw[idx][col + j % NVo] = o[ti]; cf[(si,idx,ti)] += 1

    print(f"\n{'alpha':>6} {'eta_eff':>8} {'MD':>9} {'floor(eta)':>10} "
          f"{'MD/floor':>9}  (uniform dev<={EPS})")
    for alpha in ALPHAS:
        def wgt(ca, cb):
            Da = (1-alpha) + alpha*Dstar(np.abs(ca))
            Db = (1-alpha) + alpha*Dstar(np.abs(cb))
            return Da*Db
        # per pair: N (signed, weighted), Z (weight) per column
        Np, Zp = [], []
        for i,(a,b) in enumerate(settings):
            ca = np.cos(2*(lam_c-a*_deg)); cb = np.cos(2*(lam_c+np.pi/2-b*_deg))
            W = wgt(ca, cb)
            Np.append((Qs_raw[i]*(W*np.sign(ca)*np.sign(cb))[None,:]).sum(axis=1)*DL)
            Zp.append((Qs_raw[i]*W[None,:]).sum(axis=1)*DL)
        Nf, Zf = [], []
        for k, td in enumerate(angF):
            ca = np.cos(2*lam_c); cb = np.cos(2*(lam_c+np.pi/2-td*_deg))
            W = wgt(ca, cb)
            Nf.append((Qf_raw[k]*(W*np.sign(ca)*np.sign(cb))[None,:]).sum(axis=1)*DL)
            Zf.append((Qf_raw[k]*W[None,:]).sum(axis=1)*DL)
        # LP: rho >= 0, sum=1; flat rates Z1=Z2=Z3=Z4; S = 2sqrt2 * Zcommon;
        # melody |Nf + cos * Zf| <= eps * Zf per angle; minimize source MD.
        cmb = [(i,j) for i in range(4) for j in range(i+1,4)]
        nv = NV + 1 + len(cmb)*RED
        rows, bub = [], []
        binw = np.pi/RED
        for ci,(i,j) in enumerate(cmb):
            Dif = (Qs_raw[i]-Qs_raw[j])*binw
            for l in range(RED):
                r = np.zeros(nv); r[:NV] = Dif[:,l]; r[NV+1+ci*RED+l] = -1
                rows.append(r); bub.append(0.0)
                r = np.zeros(nv); r[:NV] = -Dif[:,l]; r[NV+1+ci*RED+l] = -1
                rows.append(r); bub.append(0.0)
        for ci in range(len(cmb)):
            r = np.zeros(nv); r[NV+1+ci*RED:NV+1+(ci+1)*RED] = 0.5
            r[NV] = -1; rows.append(r); bub.append(0.0)
        for k, td in enumerate(angF):
            tgt = -np.cos(2*td*_deg)
            r = np.zeros(nv); r[:NV] = Nf[k] - (tgt+EPS)*Zf[k]
            rows.append(r); bub.append(0.0)
            r = np.zeros(nv); r[:NV] = -(Nf[k] - (tgt-EPS)*Zf[k])
            rows.append(r); bub.append(0.0)
        Aub = np.array(rows); bv = np.array(bub)
        S_row = np.zeros(NV)
        for i, sg in enumerate(signs): S_row += sg*Np[i]
        Aeq = [np.concatenate([np.ones(NV), np.zeros(nv-NV)])]
        beq = [1.0]
        for i in range(1, 4):
            Aeq.append(np.concatenate([Zp[0]-Zp[i], np.zeros(nv-NV)]))
            beq.append(0.0)
        Aeq.append(np.concatenate([S_row + S_QM*Zp[0], np.zeros(nv-NV)]))
        beq.append(0.0)   # natural branch: S = -2sqrt2 * Zcommon
        cobj = np.zeros(nv); cobj[NV] = 1
        res = linprog(cobj, A_ub=Aub, b_ub=bv, A_eq=np.array(Aeq), b_eq=beq,
                      bounds=[(0,None)]*nv, method="highs")
        if res.status != 0:
            print(f"{alpha:6.2f}   infeasible at this register")
            continue
        rho = res.x[:NV]
        qm = [Qs_raw[p].T@rho for p in range(4)]
        md = max(0.5*np.sum(np.abs(qm[a]-qm[b]))*DL
                 for a in range(4) for b in range(a+1,4))
        # effective single-arm efficiency: mean detection prob
        ca = np.cos(2*lam_c)
        Dbar = ((1-alpha) + alpha*Dstar(np.abs(ca)))
        eta_eff = float(np.sum((qm[0]/max(qm[0].sum()*DL,1e-12))*Dbar)*DL)
        fl = floor_at(eta_eff)
        print(f"{alpha:6.2f} {eta_eff:8.3f} {md:9.5f} {fl:10.5f} "
              f"{md/fl if fl>1e-6 else float('inf'):9.3f}")
    print("\nREADING: MD/floor near 1.0 at some alpha => the orchestra sits")
    print("ON the certified surface with the melody held uniformly - the")
    print("complete story. alpha=1 approaches the paper's exact-cosine")
    print("detection channel; lower alpha = weaker detection shaping.")
