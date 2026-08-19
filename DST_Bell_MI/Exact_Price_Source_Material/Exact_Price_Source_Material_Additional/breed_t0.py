#!/usr/bin/env python3
"""
breed_t0.py — breed the base FOR the family, at last.
Every number so far was achieved with a base bred for the WRONG
objective (kinetic single-distribution efficiency). This greedy
evolution optimizes the thing we actually care about: the pooled-census
mixture's minimum MD at S = 2*sqrt(2) with melody dev <= 0.10 (the
working register). Census per candidate: per-station offsets x flight-
time spread {1, 4, residence-infinity} — the winning structure from
flight_time.py.
Moves per round: add a new beat term (3 trial weights), rescale an
existing term (x0.5 / x1.5), or nudge the drift (+-20%). Best move is
kept; checkpointed every round to breed_results.json; the champion is
also written in greedy_results.json format to bred_champion.json so all
downstream tools can use it (rename to greedy_results.json to adopt).
Run: CORES=56 ROUNDS=6 python breed_t0.py    (~30-60 min/round)
"""
import os, json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linprog

_deg = np.pi/180; G = 256; RED = 128; NANG = 9
NSH = int(os.environ.get("NSH", 8))
NW = int(os.environ.get("NW", 1024))
TAUS = [1.0, 4.0]                     # plus residence-infinity, added below
settings = [(0.0,22.5),(0.0,67.5),(45.0,22.5),(45.0,67.5)]
signs = [+1,-1,+1,+1]
S_QM = 2*np.sqrt(2); FLOOR = 0.138071
lam_c = (np.arange(G)+0.5)*np.pi/G; DL = np.pi/G
SGA = [np.sign(np.cos(2*(lam_c-a*_deg))) for a,b in settings]
SGB = [np.sign(np.cos(2*(lam_c+np.pi/2-b*_deg))) for a,b in settings]
xi = np.arange(NSH)*np.pi/NSH
angs = np.linspace(8, 172, NANG)
MAXO = 6
ALL_TERMS = []
for m in range(1, MAXO+1):
    for n in range(1, MAXO+1):
        ALL_TERMS.append((m,n,+1))
        if m != n: ALL_TERMS.append((m,n,-1))

GLOB = {}
def init_worker(g): GLOB.update(g)

def _force(lam, terms, w, V, xa, xb, a, b):
    dU = np.zeros_like(lam)
    for wt,(m,n,sg) in zip(w, terms):
        if wt < 1e-12: continue
        k = m+n if sg > 0 else m-n
        ph = 2*m*(a+xa) + sg*2*n*(b+xb) + sg*n*np.pi
        dU += wt*2*k*np.sin(2*k*lam - ph)
    return V - dU

def component_job(args):
    """one (xa, xb, a, b): returns [hist@tau1, hist@tau4, residence-inf]"""
    xa, xb, a, b = args
    terms = GLOB["terms"]; w = GLOB["w"]; V = GLOB["V"]
    lam = np.linspace(0, np.pi, NW, endpoint=False).copy()
    dt = 0.005; t = 0.0; out = []
    for tau in TAUS:
        for _ in range(int(round((tau-t)/dt))):
            lam += _force(lam, terms, w, V, xa, xb, a, b)*dt
        t = tau
        h = np.bincount(((lam % np.pi)/np.pi*G).astype(int) % G,
                        minlength=G).astype(float)
        out.append(h/(h.sum()*DL))
    # residence law (tau = infinity), exact for running / delta for locked
    lg = np.linspace(0, np.pi, G, endpoint=False)
    vel = _force(lg, terms, w, V, xa, xb, a, b)
    if vel.min() > 1e-6:
        p = 1.0/vel
    else:
        z = np.where(np.diff(np.sign(vel)) < 0)[0]
        p = np.zeros(G)
        if len(z):
            p[z[0]] = 1.0
            p = np.convolve(p, np.exp(-0.5*(np.arange(-6,7)/2.0)**2), "same")
        else:
            p = np.full(G, 1.0)
    out.append(p/(p.sum()*DL))
    return out

def evaluate(terms, w, V, _unused=None):
    """objective: min MD at S=2sqrt2, dev<=0.10, over the pooled census."""
    g = dict(terms=terms, w=list(w), V=V)
    jobs, meta = [], []
    for pi,(a,b) in enumerate(settings):
        for xa in xi:
            for xb in xi:
                jobs.append((xa, xb, a*_deg, b*_deg)); meta.append(("p", pi))
    for ai, tdeg in enumerate(angs):
        for xa in xi:
            for xb in xi:
                jobs.append((xa, xb, 0.0, tdeg*_deg)); meta.append(("a", ai))
    with Pool(GLOB_CORES, initializer=init_worker, initargs=(g,)) as p2:
        outs = p2.map(component_job, jobs, chunksize=4)
    NT = len(TAUS)+1
    NVo = NSH*NSH
    NV = NVo*NT
    Qs = [np.zeros((NV, G)) for _ in range(4)]
    Ec = np.zeros((NANG, NV))
    cp = [0]*4; ca = [0]*NANG
    for (kind, idx), o in zip(meta, outs):
        if kind == "p":
            for ti in range(NT):
                Qs[idx][ti*NVo + cp[idx]] = o[ti]
            cp[idx] += 1
        else:
            tr = angs[idx]*_deg
            ob = np.sign(np.cos(2*(lam_c+np.pi/2-tr)))
            for ti in range(NT):
                Ec[idx, ti*NVo + ca[idx]] = np.sum(o[ti]*SGA[0]*ob)*DL
            ca[idx] += 1
    Svec = np.zeros(NV)
    for i, sg in enumerate(signs):
        Svec += sg*(Qs[i]*(SGA[i]*SGB[i])[None,:]).sum(axis=1)*DL
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
        r = np.zeros(nv); r[:NV] = Ec[ai]; rows.append(r); bub.append(tgt[ai]+0.10)
        r = np.zeros(nv); r[:NV] = -Ec[ai]; rows.append(r); bub.append(-tgt[ai]+0.10)
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

GLOB_CORES = 8

if __name__ == "__main__":
    GLOB_CORES = int(os.environ.get("CORES", os.cpu_count() or 8))
    rounds = int(os.environ.get("ROUNDS", 6))
    gd = json.load(open("greedy_results.json"))
    terms = [tuple(t) for t in gd["terms"]]
    w = list(gd["w"]); V = float(gd["v"])
    keep = [i for i, x in enumerate(w) if x > 1e-3]
    terms = [terms[i] for i in keep]; w = [w[i] for i in keep]
    dummy = None
    base_md = evaluate(terms, w, V)
    print(f"starting base ({len(terms)} terms): "
          + (f"{base_md/FLOOR:.3f}x floor" if base_md else "infeasible"))
    history = [dict(round=0, md=base_md, move="start")]
    for rd in range(1, rounds+1):
        t0 = time.time()
        moves = []
        for t in ALL_TERMS:
            if t in terms: continue
            for tw in (0.3, 0.7, 1.2):
                moves.append(("add", t, tw))
        for j in range(len(terms)):
            moves.append(("scale", j, 0.5))
            moves.append(("scale", j, 1.5))
        moves.append(("drift", None, 0.8))
        moves.append(("drift", None, 1.2))
        best = (base_md, None)
        for mv in moves:
            kind, arg, val = mv
            t2, w2, V2 = list(terms), list(w), V
            if kind == "add": t2 = terms + [arg]; w2 = w + [val]
            elif kind == "scale": w2 = list(w); w2[arg] *= val
            else: V2 = V*val
            md = evaluate(t2, w2, V2)
            if md and md < best[0]:
                best = (md, (t2, w2, V2, mv))
        if best[1] is None:
            print(f"round {rd}: no improving move ({time.time()-t0:.0f}s) "
                  f"- BRED CEILING at {base_md/FLOOR:.3f}x")
            break
        base_md, (terms, w, V, mv) = best
        print(f"round {rd}: {mv[0]} {mv[1]} x{mv[2]} -> "
              f"{base_md/FLOOR:.3f}x floor ({time.time()-t0:.0f}s)")
        history.append(dict(round=rd, md=base_md, move=str(mv)))
        json.dump(dict(history=history), open("breed_results.json","w"),
                  indent=1)
        json.dump(dict(history=[dict(eff=FLOOR/base_md)],
                       terms=[list(t) for t in terms], w=w, v=V, T=0.0),
                  open("bred_champion.json","w"), indent=1)
    print(f"\nfinal: {base_md/FLOOR:.3f}x floor at dev<=0.10 "
          f"(entering benchmark: 1.204 pooled / 1.431 banked)")
    print("champion saved to bred_champion.json (greedy_results.json format;")
    print("rename to adopt it in all downstream tools).")
