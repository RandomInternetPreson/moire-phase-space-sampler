#!/usr/bin/env python3
"""
birth_model_v3.py -- a birth model that knows the terrain's depth, with the harmonics the floor census uses.

What v1/v2 established (batches 4-5, note v17 Result S): a ~30-parameter von Mises push model, cell-separable by
construction, reaches a 2.4-2.6% miss cold or warm-started, finds the LP's axis / slots / (2,2) phase, and over-sharpens.
What the LP classes established: a FREE dial-blind shape reaches 0.4% (price 1.47 on the lattice), so the wall was the
form; the smooth K=6 census pays 1.70; the census's shape varies with the dial, scale more than tau; the free shape
uses (4,-4), (4,2), (2,-4) -- harmonics beyond cos 4D.

v3, per family j and per (tau, scale) cell:
    f_j(S, D | cell) = exp( logw + kS [cos(2S-2muS) - 1] + kS2 [cos(4S-4muS) - 1]                       (axis, with its 90-degree harmonic)
                          + sum_{n=1..4} kDn [cos(2n(D - d)) - 1] ),   d = muD + beta sin(2(S - muS))    (push, cos 2D .. cos 8D, axis-coupled)
    every parameter p = p0 + p_s * ds [+ p_t * dt],   ds = (scale - 1.8)/0.8,  dt = log2(tau)/3          (the birth knows the depth; optionally the duration)
    rho(cell, xi) = W[cell] * (1 + alpha cos(4 xi_A - phi)) * sum_j f_j(S, D | cell)
Parameters: 10 per family (+10 scale slopes, +10 tau slopes if TAU_DEP) + 15 cell weights + 2 station.  J=2, scale only: 57.

Stages (as v2): FIT the form to the K=6 LATTICE census by least squares (NSTART starts) -> report relative error and the
missing pixels -> differential evolution in MODE (miss / price) seeded around the fit -> Powell polish -> verify, spectrum,
json in the LP census schema.

Env: FAMILIES (2), SCALE_DEP (1), TAU_DEP (0), STATION (1), MODE (miss|price), EPS (0.01), SEED (1), TIME_LIMIT (10800),
     POP (10), NSTART (12), K6_JSON (census_bandlimit_v4_n32_k6_rho_nsh32_grid0_eps0.01_nang161.json), K6_K (6), OUT_TAG.
Usage (covariance_audit):  python birth_model_v3.py rho_nsh32_grid0_eps0.01_nang161.json > bm3_<tag>.txt
Compare: lattice K=6 census 1.700 (feasible); free dial-blind 1.473 (0.38% miss); K=5 census misses 1.055%; v1/v2 form 2.4%.
"""
import os, sys, json, time
import numpy as np
from scipy.optimize import least_squares, differential_evolution, minimize

_deg = np.pi / 180; G = 256; DL = np.pi / G
S_QM = 2 * np.sqrt(2); FLOOR = (np.sqrt(2.0) - 1.0) / 3.0
REGISTER = os.environ.get("REGISTER", "paper")
J = int(os.environ.get("FAMILIES", 2)); SCALE_DEP = os.environ.get("SCALE_DEP", "1") == "1"; TAU_DEP = os.environ.get("TAU_DEP", "0") == "1"
STATION = os.environ.get("STATION", "1") == "1"; MODE = os.environ.get("MODE", "miss"); EPS = float(os.environ.get("EPS", 0.01)); PEN = float(os.environ.get("PEN", 20))
SEED = int(os.environ.get("SEED", 1)); TIME_LIMIT = float(os.environ.get("TIME_LIMIT", 10800)); POP = int(os.environ.get("POP", 10)); NSTART = int(os.environ.get("NSTART", 12))
K6_JSON = os.environ.get("K6_JSON", "census_bandlimit_v4_n32_k6_rho_nsh32_grid0_eps0.01_nang161.json"); K6_K = int(os.environ.get("K6_K", 6))
OUT_TAG = os.environ.get("OUT_TAG", f"j{J}_s{int(SCALE_DEP)}t{int(TAU_DEP)}_{MODE}_s{SEED}")
TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]; SCALES = [1.0, 1.8, 2.6]
PNAMES = ["logw", "muS", "kS", "kS2", "muD", "kD1", "kD2", "kD3", "kD4", "beta"]
BASE_LO = [-3, 0, 0, -30, 0, -30, -30, -30, -30, -np.pi / 2]; BASE_HI = [3, np.pi, 40, 30, np.pi, 40, 30, 30, 30, np.pi / 2]
SLOPE_LO = [-2, -np.pi / 4, -15, -15, -np.pi / 4, -15, -15, -15, -15, -np.pi / 4]; SLOPE_HI = [2, np.pi / 4, 15, 15, np.pi / 4, 15, 15, 15, 15, np.pi / 4]

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

def basis(offs, K):
    xa, xb = offs[:, 0], offs[:, 1]; cols = [np.ones(len(offs))]
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            arg = 2 * (ka * xa + kb * xb); cols.append(np.cos(arg)); cols.append(np.sin(arg))
    return np.column_stack(cols)

class Model:
    def __init__(self, offs, nflat):
        self.xa, self.xb = offs[:, 0], offs[:, 1]; self.nflat = nflat; self.NV = len(offs)
        self.S2 = 2 * (self.xa + self.xb); self.D2 = 2 * (self.xa - self.xb); self.A4 = 4 * self.xa
        nsc = len(SCALES); self.ds = np.array([(SCALES[c % nsc] - 1.8) / 0.8 for c in range(nflat)])
        self.dt = np.array([np.log2(TAUS[c // nsc]) / 3 for c in range(nflat)])
        self.names, lo, hi = [], [], []
        for j in range(J):
            for n, a, b in zip(PNAMES, BASE_LO, BASE_HI): self.names.append(f"{n}{j}"); lo.append(a); hi.append(b)
            if SCALE_DEP:
                for n, a, b in zip(PNAMES, SLOPE_LO, SLOPE_HI): self.names.append(f"{n}{j}_s"); lo.append(a); hi.append(b)
            if TAU_DEP:
                for n, a, b in zip(PNAMES, SLOPE_LO, SLOPE_HI): self.names.append(f"{n}{j}_t"); lo.append(a); hi.append(b)
        for c in range(nflat): self.names.append(f"cell{c}"); lo.append(-4); hi.append(4)
        self.names += ["alpha_raw", "phi"]; lo += [-3, 0]; hi += [3, 2 * np.pi]
        self.lo = np.array(lo, float); self.hi = np.array(hi, float); self.P = len(self.names)
        self.nfam = 10 * (1 + int(SCALE_DEP) + int(TAU_DEP))
    def census(self, x):
        rho = np.zeros((self.nflat, self.NV)); k = 0; S2 = self.S2[None, :]; D2 = self.D2[None, :]
        for j in range(J):
            base = x[k:k + 10]; k += 10
            sl_s = x[k:k + 10] if SCALE_DEP else np.zeros(10); k += 10 if SCALE_DEP else 0
            sl_t = x[k:k + 10] if TAU_DEP else np.zeros(10); k += 10 if TAU_DEP else 0
            P = base[None, :] + sl_s[None, :] * self.ds[:, None] + sl_t[None, :] * self.dt[:, None]      # (nflat, 10): per-cell parameters
            logw, muS, kS, kS2, muD, kD1, kD2, kD3, kD4, beta = (P[:, i][:, None] for i in range(10))     # each (nflat, 1)
            d2 = 2 * (muD + beta * np.sin(S2 - 2 * muS)); dd = D2 - d2
            e = (logw + kS * (np.cos(S2 - 2 * muS) - 1) + kS2 * (np.cos(2 * S2 - 4 * muS) - 1)
                 + kD1 * (np.cos(dd) - 1) + kD2 * (np.cos(2 * dd) - 1) + kD3 * (np.cos(3 * dd) - 1) + kD4 * (np.cos(4 * dd) - 1))
            rho += np.exp(np.clip(e, -60, 60))
        cells = x[k:k + self.nflat]; k += self.nflat
        W = np.exp(cells - cells.max()); W /= W.sum()
        rho *= W[:, None]
        if STATION: rho *= (1 + np.tanh(x[k]) * np.cos(self.A4 - x[k + 1]))[None, :]
        rho = rho.reshape(-1); s = rho.sum()
        return rho / s if s > 0 else rho

class Objective:
    def __init__(self, m, Qs, Svec, Ecoef, tgt):
        self.m = m; self.Qs = Qs; self.Svec = Svec; self.Ecoef = Ecoef; self.tgt = tgt
        self.n = 0; self.t0 = time.time(); self.best = (np.inf, None); self.best_miss = (np.inf, None)
    def misses(self, rho):
        return float(np.abs(self.Ecoef @ rho - self.tgt).max()), float(abs(self.Svec @ rho + S_QM))
    def __call__(self, x):
        self.n += 1; rho = self.m.census(x); mE, mS = self.misses(rho); miss = max(mE, mS)
        if miss < self.best_miss[0]: self.best_miss = (miss, x.copy())
        f = miss if MODE == "miss" else md_of(self.Qs, rho) / FLOOR + PEN * (max(0.0, mE - EPS) + max(0.0, mS - EPS)) / EPS
        if f < self.best[0]: self.best = (f, x.copy())
        return f
    def stop(self, xk=None, convergence=None): return (time.time() - self.t0) > TIME_LIMIT

def spectrum(rho, offs, nflat, NV, K=6):
    wo = rho.reshape(nflat, NV).sum(axis=0); wo = wo / wo.sum(); xa, xb = offs[:, 0], offs[:, 1]; pix = []
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            F = np.sum(wo * np.exp(2j * (ka * xa + kb * xb))); pix.append((abs(F), ka, kb, np.degrees(np.angle(F))))
    pix.sort(reverse=True); return pix

def target_census(path, offs, nflat):
    d = json.load(open(K6_JSON)); r = next((r for r in d["results"] if r.get("K") == K6_K and "coeff" in r), None)
    if r is None: sys.exit(f"no feasible K={K6_K} entry with coefficients in {K6_JSON}")
    c = np.asarray(r["coeff"], float); Phi = basis(offs, K6_K)
    rho = np.clip(np.concatenate([Phi @ c[ts] for ts in range(nflat)]), 0, None); return rho / rho.sum(), r.get("price")

def main(path):
    t0 = time.time(); rng = np.random.default_rng(SEED)
    d, offs, nflat, Qs, Svec, Ecoef, angs = load(path); NV = len(offs); NT = nflat * NV; tgt = -np.cos(2 * angs * _deg)
    m = Model(offs, nflat); obj = Objective(m, Qs, Svec, Ecoef, tgt)
    rho_t, price_t = target_census(path, offs, nflat); tE, tS = obj.misses(rho_t)
    print(f"file: {os.path.basename(path)}   atoms {NT} = {nflat} x {NV}")
    print(f"target: K={K6_K} census from {K6_JSON}: price {price_t:.4f}, re-measured {md_of(Qs, rho_t)/FLOOR:.4f}; misses band {tE:.5f} CHSH {tS:.1e}")
    print(f"model v3: J={J}, scale-dependence={SCALE_DEP}, tau-dependence={TAU_DEP}, station={STATION}; {m.P} parameters; mode={MODE} eps={EPS} seed={SEED} pop={POP}x{m.P} time limit {TIME_LIMIT:.0f}s", flush=True)
    # stage 1: fit
    def resid(x): return (m.census(x) - rho_t) * NT
    best = (np.inf, None); t1 = time.time()
    print(f"\nstage 1: least-squares fit to the K={K6_K} census, {NSTART} starts", flush=True)
    for i in range(NSTART):
        x0 = m.lo + rng.random(m.P) * (m.hi - m.lo)
        if SCALE_DEP or TAU_DEP:   # start slopes near zero: the cell-independent form first
            for n, name in enumerate(m.names):
                if name.endswith("_s") or name.endswith("_t"): x0[n] *= 0.1
        try:
            r = least_squares(resid, x0, bounds=(m.lo, m.hi), max_nfev=1500, xtol=1e-8, ftol=1e-8)
            rel = np.linalg.norm(m.census(r.x) - rho_t) / np.linalg.norm(rho_t); fe, fs = obj.misses(m.census(r.x))
            print(f"  start {i}: rel. L2 error {rel:.4f}   fitted census misses band {fe:.5f} CHSH {fs:.5f}   [{time.time()-t1:.0f}s]", flush=True)
            if rel < best[0]: best = (rel, r.x.copy())
        except Exception as e:
            print(f"  start {i}: failed ({e})")
    rel, xfit = best; rho_fit = m.census(xfit); fe, fs = obj.misses(rho_fit); pfit = md_of(Qs, rho_fit) / FLOOR
    print(f"\n  best fit: rel. L2 error {rel:.4f}  (v2 form on random0: 0.70);  fitted census: band miss {fe:.5f}  CHSH miss {fs:.5f}  MD/floor {pfit:.4f}")
    pt = spectrum(rho_t, offs, nflat, NV); pf = spectrum(rho_fit, offs, nflat, NV); ft = {(ka, kb): (a, ph) for a, ka, kb, ph in pf}
    print("  leading harmonics, target vs fit  (|F| target, |F| fit, phase target, phase fit):")
    for a, ka, kb, ph in pt[:12]:
        af, phf = ft.get((ka, kb), (0.0, 0.0)); flag = "" if af > 0.5 * a else ("   <-- missing in the form" if af < 0.25 * a else "   <-- weak in the form")
        print(f"     ({ka:+d},{kb:+d})   {a:.3f}  {af:.3f}   {ph:+8.1f}  {phf:+8.1f}{flag}")
    # stage 2: DE seeded around the fit
    print(f"\nstage 2: differential evolution ({MODE}) seeded around the fit", flush=True)
    M = POP * m.P; init = np.empty((M, m.P)); init[0] = xfit; n_pert = int(0.6 * M)
    for i in range(1, n_pert + 1): init[i] = np.clip(xfit + rng.normal(0, 0.05, m.P) * (m.hi - m.lo), m.lo, m.hi)
    for i in range(n_pert + 1, M): init[i] = m.lo + rng.random(m.P) * (m.hi - m.lo)
    obj.t0 = time.time(); obj(xfit)
    differential_evolution(obj, bounds=list(zip(m.lo, m.hi)), seed=SEED, popsize=POP, maxiter=100000, tol=0, polish=False, init=init,
                           mutation=(0.5, 1.0), recombination=0.7, callback=obj.stop, updating="immediate")
    print(f"  DE: {obj.n} evaluations, best objective {obj.best[0]:.6f}  [{time.time()-t0:.0f}s]", flush=True)
    try: minimize(obj, obj.best[1], method="Powell", bounds=list(zip(m.lo, m.hi)), options={"maxiter": 20, "xtol": 1e-4, "ftol": 1e-7})
    except Exception as e: print(f"  polish skipped: {e}")
    x = obj.best[1]; rho = m.census(x); mE, mS = obj.misses(rho); price = md_of(Qs, rho) / FLOOR
    feasible = (mE <= EPS + 1e-9) and (mS <= EPS + 1e-9); bm, _ = obj.best_miss
    print(f"\nRESULT  v3 J={J} scale-dep={SCALE_DEP} tau-dep={TAU_DEP} mode={MODE} eps={EPS}")
    print(f"  model census: band miss {mE:.5f}   CHSH miss {mS:.5f}   -> {'FEASIBLE' if feasible else 'not feasible'} at eps={EPS}")
    print(f"  MD/floor = {price:.4f}   (lattice: K=6 census 1.700; free dial-blind 1.473; LP optimum 1.000)")
    print(f"  smallest max-miss seen: {bm:.5f}   (v1/v2 form: 0.024; free dial-blind shape: 0.0038; K=5 census: 0.01055)")
    print("\n  parameters (base; per-unit-scale slope; per-unit-tau slope):"); k = 0
    for j in range(J):
        base = x[k:k + 10]; k += 10; ss = x[k:k + 10] if SCALE_DEP else np.zeros(10); k += 10 if SCALE_DEP else 0; st = x[k:k + 10] if TAU_DEP else np.zeros(10); k += 10 if TAU_DEP else 0
        print(f"    family {j}: weight {np.exp(base[0]):.3f}; axis muS {np.degrees(base[1]):6.1f} deg (slope_s {np.degrees(ss[1]):+5.1f}, slope_t {np.degrees(st[1]):+5.1f}); kS {base[2]:5.1f} kS2 {base[3]:+5.1f}; push muD {np.degrees(base[4]):6.1f} deg (slope_s {np.degrees(ss[4]):+5.1f}); kD1-4 {base[5]:+5.1f} {base[6]:+5.1f} {base[7]:+5.1f} {base[8]:+5.1f}; beta {np.degrees(base[9]):+6.1f} deg")
        print(f"      axis by scale: " + "  ".join(f"s={s}: {np.degrees(base[1] + ss[1] * (s - 1.8) / 0.8) % 180:.1f}" for s in SCALES) + "   push by scale: " + "  ".join(f"s={s}: {np.degrees(base[4] + ss[4] * (s - 1.8) / 0.8) % 180:.1f}" for s in SCALES))
    cells = x[k:k + nflat]; W = np.exp(cells - cells.max()); W /= W.sum(); k += nflat
    if nflat == 15:
        print("    cell weights (rows tau, cols scale):")
        for ti, tau in enumerate(TAUS): print(f"      tau={tau:<4}" + "".join(f"  {W[ti*3+si]:.3f}" for si in range(3)))
    if STATION: print(f"    station factor: alpha {np.tanh(x[k]):+.3f}  phi {np.degrees(x[k+1]):6.1f} deg")
    pix = spectrum(rho, offs, nflat, NV)
    print("\n  leading harmonics of the final model census  (|F|, k_A, k_B, phi[deg]):")
    for a, ka, kb, ph in pix[:10]: print(f"     {a:.3f}  ({ka:+d},{kb:+d})  {ph:+8.2f}")
    out = dict(d); out.update(rho=rho.tolist(), md=float(price * FLOOR), model=dict(version=3, J=J, scale_dep=SCALE_DEP, tau_dep=TAU_DEP, station=STATION, mode=MODE, eps=EPS, seed=SEED,
               target=K6_JSON, fit_rel_error=float(rel), fit_band_miss=float(fe), fit_price=float(pfit), params=dict(zip(m.names, map(float, x))),
               band_miss=float(mE), chsh_miss=float(mS), price=float(price), feasible=bool(feasible), min_miss_seen=float(bm), evaluations=obj.n,
               leading=[(float(a), int(ka), int(kb), float(ph)) for a, ka, kb, ph in pix[:20]]))
    jpath = os.path.join(os.path.dirname(path) or ".", f"birth_model_v3_{OUT_TAG}_{os.path.splitext(os.path.basename(path))[0]}.json"); json.dump(out, open(jpath, "w"))
    print(f"\nSUMMARY  v3 J={J} sdep={int(SCALE_DEP)} tdep={int(TAU_DEP)} mode={MODE} eps={EPS}: fit rel.err {rel:.4f} | final band miss {mE:.5f}  CHSH miss {mS:.5f}  price {price:.4f}  {'FEASIBLE' if feasible else 'infeasible'}  min-miss-seen {bm:.5f}   [{time.time()-t0:.0f}s, {obj.n} evals]")
    print(f"wrote {os.path.basename(jpath)}")

if __name__ == "__main__":
    main(sys.argv[1])
