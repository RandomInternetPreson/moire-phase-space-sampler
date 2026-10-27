#!/usr/bin/env python3
"""
birth_model_v1.py -- the first CONSTRUCTIVE test: can a low-parameter birth model attain Hall's floor?

Everything since 15 Oct has been exclusion (not station-local, not (S,D)-separable, not a short mixture, not smooth
below K = 6).  This script asks the positive question.  A birth model here is a parametric density over the SAME
15360 atoms (5 tau x 3 scale cells x 1024 offsets) that the LP censuses live on, built from the geometry the
censuses showed:

    rho(xi_A, xi_B, cell) = W[cell] * g(xi_A) * sum_j w_j * L_S,j(S) * L_D,j(D ; S)

    S = xi_A + xi_B, D = xi_A - xi_B  (period pi each; harmonics in 2S, 2D)
    L_S,j(S)     = exp( kS_j cos(2S - 2 muS_j) )                              a preferred birth axis
    L_D,j(D; S)  = exp( kD_j cos(2D - 2 d_j(S)) + kD2_j cos(4D - 4 d_j(S)) )   a push, with first and second harmonic
    d_j(S)       = muD_j + beta_j sin(2(S - muS_j))                            push depends on the axis orientation (the S-D coupling)
    g(xi_A)      = 1 + alpha cos(4 xi_A - phi)                                 optional station factor (the 90-degree slots seen on Alice)
    W[cell]      = softmax over the 15 (tau, scale) cells                      how the birth populates the dial

J families (FAMILIES env): J = 1 with beta = 0 is a pure (S,D) product (already excluded at eps = 0.01 by the LP -- a
control); J = 2 is the aligned + crossed pair seen in the slices; J = 3 one spare.  7 parameters per family + 15 + 2.

Two objectives (MODE env):
    miss   minimize  max( max_theta |E(theta) + cos 2theta| ,  |S_CHSH - 2 sqrt 2| )     -- the smallest tolerance the model class reaches
    price  minimize  MD/floor + PEN * (excess of both misses over eps)/eps                 -- the price at eps, once feasible
Differential evolution (scipy, serial), then a Powell polish; time-capped (TIME_LIMIT s).  No LP: every evaluation is
a few matrix products against the banks (~10 ms), so a run is ~10^5 evaluations.  Prices here are a priori >= the
LP's (the LP is the exact optimum over ALL censuses); what is new is how close a model with ~30 parameters gets.

Constraint note: the LP imposes CHSH = 2 sqrt 2 exactly; a parametric model cannot hit an equality, so the model is
held to |S - 2 sqrt 2| <= eps as well as the band.  Stated in the output.

Usage (covariance_audit folder; rho file next to its banks):
  set FAMILIES=2      set COUPLING=1      set STATION=1      set MODE=price      set EPS=0.01
  set SEED=1          set TIME_LIMIT=3000 set OUT_TAG=bm_j2_price_s1
  python birth_model_v1.py rho_nsh32_random0_eps0.01_nang161.json > bm_<tag>.txt
Writes birth_model_v1_<tag>_<stem>.json: parameters, misses, price, the model census's leading harmonics, and
the model census itself (15360 weights) so census_structure_v2.py / census_sumdiff_v1.py can read it like an LP census
(same schema: "rho", "offsets", "md", ...).
"""
import os, sys, json, time
import numpy as np
from scipy.optimize import differential_evolution, minimize

_deg = np.pi / 180; G = 256; DL = np.pi / G
S_QM = 2 * np.sqrt(2); FLOOR = (np.sqrt(2.0) - 1.0) / 3.0
REGISTER = os.environ.get("REGISTER", "paper")
J = int(os.environ.get("FAMILIES", 2)); COUPLING = os.environ.get("COUPLING", "1") == "1"; STATION = os.environ.get("STATION", "1") == "1"
MODE = os.environ.get("MODE", "price"); EPS = float(os.environ.get("EPS", 0.01)); PEN = float(os.environ.get("PEN", 20))
SEED = int(os.environ.get("SEED", 1)); TIME_LIMIT = float(os.environ.get("TIME_LIMIT", 3000))
POP = int(os.environ.get("POP", 12)); OUT_TAG = os.environ.get("OUT_TAG", f"j{J}_{MODE}_s{SEED}")
TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]; SCALES = [1.0, 1.8, 2.6]

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

# ---------------- the model ----------------
class Model:
    def __init__(self, offs, nflat):
        self.xa, self.xb = offs[:, 0], offs[:, 1]; self.nflat = nflat; self.NV = len(offs)
        self.S2 = 2 * (self.xa + self.xb); self.D2 = 2 * (self.xa - self.xb); self.A4 = 4 * self.xa
        # parameter layout
        self.names = []; self.lo = []; self.hi = []
        for j in range(J):
            self.names += [f"logw{j}", f"muS{j}", f"kS{j}", f"muD{j}", f"kD{j}", f"kD2{j}", f"beta{j}"]
            self.lo += [-3, 0, 0, 0, 0, -10, -np.pi / 2]; self.hi += [3, np.pi, 50, np.pi, 50, 10, np.pi / 2]
        for c in range(nflat):
            self.names.append(f"cell{c}"); self.lo.append(-4); self.hi.append(4)
        self.names += ["alpha_raw", "phi"]; self.lo += [-3, 0]; self.hi += [3, 2 * np.pi]
        self.lo = np.array(self.lo, float); self.hi = np.array(self.hi, float); self.P = len(self.names)
    def census(self, x):
        f = np.zeros(self.NV); k = 0
        for j in range(J):
            logw, muS, kS, muD, kD, kD2, beta = x[k:k + 7]; k += 7
            if not COUPLING: beta = 0.0
            d = muD + beta * np.sin(self.S2 - 2 * muS)
            f += np.exp(logw + kS * (np.cos(self.S2 - 2 * muS) - 1) + kD * (np.cos(self.D2 - 2 * d) - 1) + kD2 * (np.cos(2 * self.D2 - 4 * d) - 1))
        cells = x[k:k + self.nflat]; k += self.nflat
        W = np.exp(cells - cells.max()); W /= W.sum()
        if STATION:
            alpha = np.tanh(x[k]); phi = x[k + 1]
            f = f * (1 + alpha * np.cos(self.A4 - phi))
        rho = (W[:, None] * f[None, :]).reshape(-1)
        s = rho.sum()
        return rho / s if s > 0 else rho

class Objective:
    def __init__(self, model, Qs, Svec, Ecoef, tgt):
        self.m = model; self.Qs = Qs; self.Svec = Svec; self.Ecoef = Ecoef; self.tgt = tgt
        self.n = 0; self.t0 = time.time(); self.best = (np.inf, None); self.best_miss = (np.inf, None)
    def misses(self, rho):
        E = self.Ecoef @ rho; missE = np.abs(E - self.tgt).max()
        missS = abs(self.Svec @ rho + S_QM)        # LP convention: Svec . rho = -S_QM
        return missE, missS
    def __call__(self, x):
        self.n += 1
        rho = self.m.census(x); missE, missS = self.misses(rho); miss = max(missE, missS)
        if miss < self.best_miss[0]: self.best_miss = (miss, x.copy())
        if MODE == "miss":
            f = miss
        else:
            price = md_of(self.Qs, rho) / FLOOR
            f = price + PEN * (max(0.0, missE - EPS) + max(0.0, missS - EPS)) / EPS
        if f < self.best[0]: self.best = (f, x.copy())
        return f
    def stop(self, xk=None, convergence=None):
        return (time.time() - self.t0) > TIME_LIMIT

def spectrum(rho, offs, nflat, NV, K=6):
    wo = rho.reshape(nflat, NV).sum(axis=0); wo = wo / wo.sum(); xa, xb = offs[:, 0], offs[:, 1]; pix = []
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            F = np.sum(wo * np.exp(2j * (ka * xa + kb * xb))); pix.append((abs(F), ka, kb, np.degrees(np.angle(F))))
    pix.sort(reverse=True); return pix

def main(path):
    t0 = time.time()
    d, offs, nflat, Qs, Svec, Ecoef, angs = load(path); NV = len(offs)
    rho_star = np.asarray(d["rho"], float); rho_star /= rho_star.sum(); md_star = md_of(Qs, rho_star)
    tgt = -np.cos(2 * angs * _deg)
    m = Model(offs, nflat); obj = Objective(m, Qs, Svec, Ecoef, tgt)
    print(f"file: {os.path.basename(path)}   atoms {nflat*NV} = {nflat} x {NV}")
    print(f"check: MD of the stored rho* = {md_star:.12f} (file {d['md']:.12f}; floor {FLOOR:.12f});  rho* misses: band {obj.misses(rho_star)[0]:.5f}  CHSH {obj.misses(rho_star)[1]:.2e}")
    print(f"model: J={J} families, coupling={COUPLING}, station factor={STATION}, {m.P} parameters;  mode={MODE} eps={EPS} seed={SEED} pop={POP}x{m.P} time limit {TIME_LIMIT:.0f}s", flush=True)
    rng = np.random.default_rng(SEED)
    res = differential_evolution(obj, bounds=list(zip(m.lo, m.hi)), seed=SEED, popsize=POP, maxiter=100000, tol=0, polish=False,
                                 init="latinhypercube", mutation=(0.5, 1.0), recombination=0.7, callback=obj.stop, updating="immediate")
    print(f"  DE: {obj.n} evaluations, best objective {obj.best[0]:.6f}  [{time.time()-t0:.0f}s]", flush=True)
    # polish (bounded Powell from the DE best, short)
    x0 = obj.best[1]; TL_left = max(60.0, TIME_LIMIT * 0.15); obj.t0 = time.time(); obj_time = TIME_LIMIT;
    try:
        r2 = minimize(obj, x0, method="Powell", bounds=list(zip(m.lo, m.hi)), options={"maxiter": 20, "xtol": 1e-4, "ftol": 1e-7})
    except Exception as e:
        print(f"  polish skipped: {e}")
    x = obj.best[1]; rho = m.census(x); missE, missS = obj.misses(rho); price = md_of(Qs, rho) / FLOOR
    feasible = (missE <= EPS + 1e-9) and (missS <= EPS + 1e-9)
    print(f"  polish: {obj.n} evaluations total, best objective {obj.best[0]:.6f}  [{time.time()-t0:.0f}s]")
    print(f"\nRESULT  mode={MODE} J={J} coupling={COUPLING} station={STATION} eps={EPS}")
    print(f"  model census: band miss {missE:.5f}   CHSH miss {missS:.5f}   -> {'FEASIBLE' if feasible else 'not feasible'} at eps={EPS}")
    print(f"  MD/floor = {price:.4f}   (LP optimum over all censuses = 1.000; smoothest K=6 census = 2.274 at eps=0.01 on random0)")
    bm, xbm = obj.best_miss; print(f"  smallest max-miss seen in this run: {bm:.5f}   (compare: best K=5 band-limited census misses 0.0120, K=4 0.0174)")
    print("\n  parameters:")
    k = 0
    for j in range(J):
        logw, muS, kS, muD, kD, kD2, beta = x[k:k + 7]; k += 7
        print(f"    family {j}: weight {np.exp(logw):.3f}  axis muS {np.degrees(muS):6.1f} deg (kappa_S {kS:5.1f})  push muD {np.degrees(muD):6.1f} deg (kappa_D {kD:5.1f}, kappa_D2 {kD2:+5.1f})  coupling beta {np.degrees(beta) if COUPLING else 0:+6.1f} deg")
    cells = x[k:k + nflat]; W = np.exp(cells - cells.max()); W /= W.sum(); k += nflat
    if nflat == len(TAUS) * len(SCALES):
        print("    cell weights (rows tau, cols scale):")
        for ti, tau in enumerate(TAUS): print(f"      tau={tau:<4}" + "".join(f"  {W[ti*len(SCALES)+si]:.3f}" for si in range(len(SCALES))))
    else:
        print("    cell weights: " + " ".join(f"{w:.3f}" for w in W))
    if STATION: print(f"    station factor: alpha {np.tanh(x[k]):+.3f}  phi {np.degrees(x[k+1]):6.1f} deg  (slots of cos 4 xi_A at {np.degrees(x[k+1])/4 % 90:.1f} deg mod 90)")
    pix = spectrum(rho, offs, nflat, NV)
    print("\n  leading harmonics of the model census  (|F|, k_A, k_B, phi[deg])  -- LP skeleton: (2,2)~0.37@157, (2,0), (2,4), (1,1), (0,2):")
    for a, ka, kb, ph in pix[:10]: print(f"     {a:.3f}  ({ka:+d},{kb:+d})  {ph:+8.2f}")
    # json in the LP census schema so the census tools read it
    out = dict(d); out.update(rho=rho.tolist(), md=float(price * FLOOR), model=dict(J=J, coupling=COUPLING, station=STATION, mode=MODE, eps=EPS, seed=SEED,
               params=dict(zip(m.names, map(float, x))), band_miss=float(missE), chsh_miss=float(missS), price=float(price), feasible=bool(feasible),
               min_miss_seen=float(bm), evaluations=obj.n, leading=[(float(a), int(ka), int(kb), float(ph)) for a, ka, kb, ph in pix[:20]]))
    jpath = os.path.join(os.path.dirname(path) or ".", f"birth_model_v1_{OUT_TAG}_{os.path.splitext(os.path.basename(path))[0]}.json")
    json.dump(out, open(jpath, "w"))
    print(f"\nSUMMARY  mode={MODE} J={J} c={int(COUPLING)} st={int(STATION)} eps={EPS}: band miss {missE:.5f}  CHSH miss {missS:.5f}  price {price:.4f}  {'FEASIBLE' if feasible else 'infeasible'}  min-miss-seen {bm:.5f}   [{time.time()-t0:.0f}s, {obj.n} evals]")
    print(f"wrote {os.path.basename(jpath)}")

if __name__ == "__main__":
    main(sys.argv[1])
