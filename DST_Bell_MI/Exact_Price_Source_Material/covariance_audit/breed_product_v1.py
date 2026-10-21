#!/usr/bin/env python3
"""
breed_product_v1.py -- breed a terrain FOR a station-local (product) census.

Question.  On the greedy terrain of arXiv:2608.18886 no product census rho = q(tau) f(xi_A) g(xi_B)
reproduces the singlet at eps <= 0.02, and where one does (eps >= 0.05) it pays 2-3x the joint price
(product_census_v2/v3, Oct 17).  That is a statement about ONE terrain.  This script asks whether a
terrain can be bred whose product-class price approaches Hall's floor -- the only route from
"not on this terrain" to "not on any terrain we can find".

Method.  breed_t0.py's greedy evolution (add a beat term at 3 trial weights / rescale a term x0.5, x1.5
/ nudge the drift +-20%; keep the best move; checkpoint every round), at breeding scale (NSH x NSH grid
of offsets, NANG register angles, flight times {1, 4, residence-infinity}, dev <= EPS), with ONE change:
the score of a candidate terrain is the best PRODUCT-class price found by alternating LPs from NSTART
starts, instead of the joint LP price.  The joint price is also computed each time, for reference, so
the gap product/joint is tracked round by round.  On the grid support the product class is exact
(w[t, i, j] = q[t] f[i] g[j]).

Scoring order: feasible (max correlator/CHSH slack < FEAS_TOL) beats infeasible; among feasible, lower
MD; among infeasible, lower slack.  So the breed first drives toward feasibility, then toward the floor.

Outputs: breed_product_results.json (history), bred_product_champion.json (greedy_results.json format,
usable by every downstream tool), and a printed round log.

Run (Windows cmd, from covariance_audit, with greedy_results.json present):
  set CORES=56
  set ROUNDS=6
  set NSH=8
  set NW=1024
  set NANG=9
  set EPS=0.10
  set NSTART=4
  set MAXIT=15
  set MAXO=6            (max beat order of terms offered as "add" moves)
  set START=greedy      (greedy = start from greedy_results.json; minimal = start from the single term (1,1,+) w=1)
  python breed_product_v1.py > breed_product_<date>.txt
Roughly 30 s per candidate, ~200 candidates per round at MAXO=6 -> ~2 h per round on 56 cores.
"""
import os, json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linprog

_deg = np.pi / 180; G = 256; RED = 128
NANG = int(os.environ.get("NANG", 9))
NSH = int(os.environ.get("NSH", 8))
NW = int(os.environ.get("NW", 1024))
EPS = float(os.environ.get("EPS", 0.10))
NSTART = int(os.environ.get("NSTART", 4))
MAXIT = int(os.environ.get("MAXIT", 15))
MAXO = int(os.environ.get("MAXO", 6))
START = os.environ.get("START", "greedy")
PEN = float(os.environ.get("PEN", 100.0))
FEAS_TOL = float(os.environ.get("FEAS_TOL", 1e-6))
SEED = int(os.environ.get("SEED", 20261018))
TAUS = [1.0, 4.0]
settings = [(0.0, 22.5), (0.0, 67.5), (45.0, 22.5), (45.0, 67.5)]
signs = [+1, -1, +1, +1]
S_QM = 2 * np.sqrt(2); FLOOR = (np.sqrt(2.0) - 1.0) / 3.0
lam_c = (np.arange(G) + 0.5) * np.pi / G; DL = np.pi / G
SGA = [np.sign(np.cos(2 * (lam_c - a * _deg))) for a, b in settings]
SGB = [np.sign(np.cos(2 * (lam_c + np.pi / 2 - b * _deg))) for a, b in settings]
xi = np.arange(NSH) * np.pi / NSH
angs = np.linspace(8, 172, NANG)
ALL_TERMS = []
for m in range(1, MAXO + 1):
    for n in range(1, MAXO + 1):
        ALL_TERMS.append((m, n, +1))
        if m != n: ALL_TERMS.append((m, n, -1))

GLOB = {}
def init_worker(g): GLOB.update(g)

def _force(lam, terms, w, V, xa, xb, a, b):
    dU = np.zeros_like(lam)
    for wt, (m, n, sg) in zip(w, terms):
        if wt < 1e-12: continue
        k = m + n if sg > 0 else m - n
        ph = 2 * m * (a + xa) + sg * 2 * n * (b + xb) + sg * n * np.pi
        dU += wt * 2 * k * np.sin(2 * k * lam - ph)
    return V - dU

def component_job(args):
    """one (xa, xb, a, b): returns [hist@tau1, hist@tau4, residence-inf]  (breed_t0.py verbatim)"""
    xa, xb, a, b = args
    terms = GLOB["terms"]; w = GLOB["w"]; V = GLOB["V"]
    lam = np.linspace(0, np.pi, NW, endpoint=False).copy()
    dt = 0.005; t = 0.0; out = []
    for tau in TAUS:
        for _ in range(int(round((tau - t) / dt))):
            lam += _force(lam, terms, w, V, xa, xb, a, b) * dt
        t = tau
        h = np.bincount(((lam % np.pi) / np.pi * G).astype(int) % G, minlength=G).astype(float)
        out.append(h / (h.sum() * DL))
    lg = np.linspace(0, np.pi, G, endpoint=False)
    vel = _force(lg, terms, w, V, xa, xb, a, b)
    if vel.min() > 1e-6:
        p = 1.0 / vel
    else:
        z = np.where(np.diff(np.sign(vel)) < 0)[0]
        p = np.zeros(G)
        if len(z):
            p[z[0]] = 1.0
            p = np.convolve(p, np.exp(-0.5 * (np.arange(-6, 7) / 2.0) ** 2), "same")
        else:
            p = np.full(G, 1.0)
    out.append(p / (p.sum() * DL))
    return out

# ----------------------------------------------------------------------------------------------
def build_banks(terms, w, V, cores):
    """returns Qs[4] (NV,G), Ec (NANG,NV), Svec (NV,); atom index = ti*NVo + ia*NSH + ib"""
    g = dict(terms=terms, w=list(w), V=V)
    jobs, meta = [], []
    for pi, (a, b) in enumerate(settings):
        for xa in xi:
            for xb in xi:
                jobs.append((xa, xb, a * _deg, b * _deg)); meta.append(("p", pi))
    for ai, tdeg in enumerate(angs):
        for xa in xi:
            for xb in xi:
                jobs.append((xa, xb, 0.0, tdeg * _deg)); meta.append(("a", ai))
    with Pool(cores, initializer=init_worker, initargs=(g,)) as p2:
        outs = p2.map(component_job, jobs, chunksize=4)
    NT = len(TAUS) + 1; NVo = NSH * NSH; NV = NVo * NT
    Qs = [np.zeros((NV, G)) for _ in range(4)]; Ec = np.zeros((NANG, NV))
    cp = [0] * 4; ca = [0] * NANG
    for (kind, idx), o in zip(meta, outs):
        if kind == "p":
            for ti in range(NT): Qs[idx][ti * NVo + cp[idx]] = o[ti]
            cp[idx] += 1
        else:
            tr = angs[idx] * _deg
            ob = np.sign(np.cos(2 * (lam_c + np.pi / 2 - tr)))
            for ti in range(NT): Ec[idx, ti * NVo + ca[idx]] = np.sum(o[ti] * SGA[0] * ob) * DL
            ca[idx] += 1
    Svec = np.zeros(NV)
    for i, sg in enumerate(signs):
        Svec += sg * (Qs[i] * (SGA[i] * SGB[i])[None, :]).sum(axis=1) * DL
    return Qs, Ec, Svec

def md_of(Qs, rho):
    qm = [Qs[p].T @ rho for p in range(4)]
    return max(0.5 * np.sum(np.abs(qm[a] - qm[b])) * DL for a in range(4) for b in range(a + 1, 4))

def reduce_banks(Qs, Svec, Ecoef, grp, c, nth):
    Qr = []
    for p in range(4):
        R = np.zeros((nth, G)); np.add.at(R, grp, Qs[p] * c[:, None]); Qr.append(R)
    Sr = np.zeros(nth); np.add.at(Sr, grp, Svec * c)
    Er = np.zeros((Ecoef.shape[0], nth)); np.add.at(Er.T, grp, (Ecoef * c[None, :]).T)
    Mr = np.bincount(grp, weights=c, minlength=nth)
    return Qr, Sr, Er, Mr

def solve_lp(Qs, Svec, Ecoef, Mr, eps, penalty):
    keep = np.where(Mr > 0)[0]; full_n = len(Mr)
    Qs = [Q[keep] for Q in Qs]; Svec = Svec[keep]; Ecoef = Ecoef[:, keep]; Mr = Mr[keep]
    NV = len(Svec); nang = len(angs)
    Qc = [Q.reshape(NV, RED, -1).mean(axis=2) for Q in Qs]
    cmb = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    n_s = len(cmb) * RED; nu = (nang + 1) if penalty else 0
    nv = NV + 1 + n_s + nu
    rows, bub = [], []; binw = np.pi / RED
    for ci, (i, j) in enumerate(cmb):
        Dif = (Qc[i] - Qc[j]) * binw
        for l in range(RED):
            r = np.zeros(nv); r[:NV] = Dif[:, l]; r[NV + 1 + ci * RED + l] = -1; rows.append(r); bub.append(0.0)
            r = np.zeros(nv); r[:NV] = -Dif[:, l]; r[NV + 1 + ci * RED + l] = -1; rows.append(r); bub.append(0.0)
    for ci in range(len(cmb)):
        r = np.zeros(nv); r[NV + 1 + ci * RED:NV + 1 + (ci + 1) * RED] = 0.5; r[NV] = -1; rows.append(r); bub.append(0.0)
    tgt = -np.cos(2 * angs * _deg); u0 = NV + 1 + n_s
    for ai in range(nang):
        r = np.zeros(nv); r[:NV] = Ecoef[ai]
        if penalty: r[u0 + ai] = -1
        rows.append(r); bub.append(tgt[ai] + eps)
        r = np.zeros(nv); r[:NV] = -Ecoef[ai]
        if penalty: r[u0 + ai] = -1
        rows.append(r); bub.append(-tgt[ai] + eps)
    Aeq = [np.zeros(nv)]; Aeq[0][:NV] = Mr; beq = [1.0]
    if penalty:
        r = np.zeros(nv); r[:NV] = Svec; r[u0 + nang] = -1; rows.append(r); bub.append(-S_QM)
        r = np.zeros(nv); r[:NV] = -Svec; r[u0 + nang] = -1; rows.append(r); bub.append(S_QM)
    else:
        r = np.zeros(nv); r[:NV] = Svec; Aeq.append(r); beq.append(-S_QM)
    cobj = np.zeros(nv); cobj[NV] = 1
    if penalty: cobj[u0:] = PEN
    res = linprog(cobj, A_ub=np.array(rows), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=beq,
                  bounds=[(0, None)] * nv, method="highs")
    if res.status != 0: return None, None, None
    theta = np.zeros(full_n); theta[keep] = res.x[:NV]
    return res.x[NV], theta, (float(res.x[u0:].max()) if penalty else 0.0)

def joint_price(Qs, Ec, Svec):
    NV = len(Svec)
    t, theta, _ = solve_lp(Qs, Svec, Ec, np.ones(NV), EPS, penalty=False)
    return None if t is None else md_of(Qs, theta)

def product_price(Qs, Ec, Svec, rng):
    """best product census q[t] f[i] g[j] by alternating penalty LPs; returns (feasible, md, viol, start)"""
    NT = len(TAUS) + 1; NVo = NSH * NSH; NV = NVo * NT
    t_of = np.repeat(np.arange(NT), NVo)
    ia = np.tile(np.repeat(np.arange(NSH), NSH), NT)
    ib = np.tile(np.tile(np.arange(NSH), NSH), NT)
    blocks = [("q", t_of, NT), ("f", ia, NSH), ("g", ib, NSH)]
    starts = [("uniform", dict(q=np.ones(NT), f=np.ones(NSH), g=np.ones(NSH)))]
    for s in range(NSTART - 1):
        starts.append((f"random{s}", dict(q=rng.random(NT), f=rng.random(NSH), g=rng.random(NSH))))
    best = None
    for sname, init in starts:
        th = {n: init[n].copy() for n, _, _ in blocks}; last = None; ok = True
        for it in range(MAXIT):
            for name, grp, nth in blocks:
                c = np.ones(NV)
                for n2, g2, _ in blocks:
                    if n2 != name: c *= th[n2][g2]
                Qr, Sr, Er, Mr = reduce_banks(Qs, Svec, Ec, grp, c, nth)
                t, theta, viol = solve_lp(Qr, Sr, Er, Mr, EPS, penalty=True)
                if t is None: ok = False; break
                th[name] = theta
            if not ok: break
            w = np.ones(NV)
            for n2, g2, _ in blocks: w *= th[n2][g2]
            assert abs(w.sum() - 1) < 1e-6
            md = md_of(Qs, w)
            if last is not None and abs(md - last[0]) < 1e-7 and abs(viol - last[1]) < 1e-9: break
            last = (md, viol)
        if not ok: continue
        feas = viol < FEAS_TOL
        if feas: assert md >= FLOOR - 1e-7
        cand = (0 if feas else 1, md if feas else viol, md, viol, sname)
        if best is None or cand[:2] < best[:2]: best = cand
    if best is None: return None
    return dict(feasible=(best[0] == 0), md=best[2], viol=best[3], start=best[4])

def evaluate(terms, w, V, cores, rng):
    Qs, Ec, Svec = build_banks(terms, w, V, cores)
    pj = joint_price(Qs, Ec, Svec)
    pp = product_price(Qs, Ec, Svec, rng)
    return pj, pp

def score(pp):
    """lexicographic: (infeasible?, md or viol); None = LP failure = worst"""
    if pp is None: return (2, np.inf)
    return (0, pp["md"]) if pp["feasible"] else (1, pp["viol"])

def fmt(pj, pp):
    sj = "infeas" if pj is None else f"{pj/FLOOR:.3f}"
    if pp is None: sp = "LPfail"
    elif pp["feasible"]: sp = f"{pp['md']/FLOOR:.3f}"
    else: sp = f"infeas(slack {pp['viol']:.3f})"
    return f"joint {sj} | product {sp}"

if __name__ == "__main__":
    cores = int(os.environ.get("CORES", os.cpu_count() or 8))
    rounds = int(os.environ.get("ROUNDS", 6))
    rng = np.random.default_rng(SEED)
    if START == "greedy":
        gd = json.load(open("greedy_results.json"))
        terms = [tuple(t) for t in gd["terms"]]; w = list(gd["w"]); V = float(gd["v"])
        keep = [i for i, x in enumerate(w) if x > 1e-3]
        terms = [terms[i] for i in keep]; w = [w[i] for i in keep]
    else:
        terms = [(1, 1, +1)]; w = [1.0]; V = 15.0
    print(f"breed_product_v1  NSH={NSH} NANG={NANG} NW={NW} EPS={EPS} NSTART={NSTART} MAXIT={MAXIT} MAXO={MAXO} START={START} cores={cores}", flush=True)
    t0 = time.time()
    pj, pp = evaluate(terms, w, V, cores, rng)
    print(f"start ({len(terms)} terms, v={V:.4f}): {fmt(pj, pp)}   [{time.time()-t0:.0f}s]", flush=True)
    base = score(pp)
    history = [dict(round=0, joint=pj, product=pp, move="start", terms=[list(t) for t in terms], w=w, v=V)]
    for rd in range(1, rounds + 1):
        t0 = time.time()
        moves = []
        for t in ALL_TERMS:
            if t in terms: continue
            for tw in (0.3, 0.7, 1.2): moves.append(("add", t, tw))
        for j in range(len(terms)):
            moves.append(("scale", j, 0.5)); moves.append(("scale", j, 1.5))
        moves.append(("drift", None, 0.8)); moves.append(("drift", None, 1.2))
        best = (base, None)
        for mi, mv in enumerate(moves):
            kind, arg, val = mv
            t2, w2, V2 = list(terms), list(w), V
            if kind == "add": t2 = terms + [arg]; w2 = w + [val]
            elif kind == "scale": w2 = list(w); w2[arg] *= val
            else: V2 = V * val
            pj2, pp2 = evaluate(t2, w2, V2, cores, rng)
            sc = score(pp2)
            if sc < best[0]:
                best = (sc, (t2, w2, V2, mv, pj2, pp2))
                print(f"  round {rd} move {mi+1}/{len(moves)} {kind} {arg} x{val}: {fmt(pj2, pp2)}  <- new best", flush=True)
            elif (mi + 1) % 25 == 0:
                print(f"  round {rd} move {mi+1}/{len(moves)}  ({time.time()-t0:.0f}s)", flush=True)
        if best[1] is None:
            print(f"round {rd}: no improving move ({time.time()-t0:.0f}s) - PRODUCT-BRED CEILING", flush=True)
            break
        base, (terms, w, V, mv, pj, pp) = best
        print(f"round {rd}: {mv[0]} {mv[1]} x{mv[2]} -> {fmt(pj, pp)}  ({len(terms)} terms, v={V:.4f})  ({time.time()-t0:.0f}s)", flush=True)
        history.append(dict(round=rd, joint=pj, product=pp, move=str(mv), terms=[list(t) for t in terms], w=w, v=V))
        json.dump(dict(history=history, meta=dict(NSH=NSH, NANG=NANG, NW=NW, EPS=EPS, NSTART=NSTART, MAXIT=MAXIT, MAXO=MAXO, START=START)),
                  open("breed_product_results.json", "w"), indent=1, default=float)
        json.dump(dict(history=[dict(product_md=pp["md"] if pp else None, joint_md=pj)], terms=[list(t) for t in terms], w=w, v=V, T=0.0),
                  open("bred_product_champion.json", "w"), indent=1, default=float)
    print(f"\nfinal: {fmt(pj, pp)}   champion in bred_product_champion.json (greedy_results.json format)")
    print("Validate the champion at scale: copy it to a fresh folder as greedy_results.json and run dither_sweep_rho.py")
    print("(NANG=161, NSH_LIST=32, MODES=random, SAVE_RHO=1) then product_census_v3.py on the rho file.")
