#!/usr/bin/env python3
"""
birth_model_v2.py -- the birth model warm-started from the smoothest floor-attaining census (the K = 6 LP census).

v1 (batch 4) searched the model's parameter space from scratch and reached a 2.5% miss: it found the LP's axis
(38 deg), its station slots (38 deg mod 90) and its (2,2) phase (156 vs 157 deg), but got there by over-sharpening
(|F(1,1)| = 0.74, |F(2,2)| = 0.70 against the LP's 0.20 and 0.37) and collapsing the cell weights onto two cells.
Non-convex, so 2.5% is an upper bound on what the form can do, not a property of the form.

v2 separates the two questions:
  stage 1  FIT.  Can the parametric form REPRESENT the K = 6 LP census (smooth, reaches 1%, price 2.274)?
           Least squares of the model density to rho_K6 over all 15360 atoms, from NSTART random starts (+ the v1
           best, if its json is present).  Reports the relative L2 error of the best fit, the fitted census's own
           miss and price, and the leading harmonics of target vs fit -- the pixels the form cannot carry are listed.
  stage 2  REFINE.  Differential evolution in MODE (miss or price) with the population seeded around the best fit
           (1 exact copy, 60% perturbed copies, the rest uniform), then a Powell polish.  Reports as v1 does.
Reading: if the fit is good (rel. error small, fit miss ~1%), stage 2 should reach <= 1% -- the v1 result was the
optimizer.  If the fit is poor, the form is missing the listed pixels -- the v1 result was the form.

Usage (covariance_audit):  same env as v1 (FAMILIES, COUPLING, STATION, MODE, EPS, SEED, TIME_LIMIT, OUT_TAG) plus
  set K6_JSON=census_bandlimit_v4_r0_e01_k6_v4_rho_nsh32_random0_eps0.01_nang161.json     (the target census file)
  set K6_K=6            which K entry of that json to use
  set NSTART=12         random least-squares starts
  set WARM=birth_model_v1_j2_miss_s1_rho_nsh32_random0_eps0.01_nang161.json   (optional: a v1 result to add as a start)
  python birth_model_v2.py rho_nsh32_random0_eps0.01_nang161.json > bm2_<tag>.txt
Writes birth_model_v2_<tag>_<stem>.json in the LP census schema (readable by census_structure_v2.py / census_sumdiff_v1.py).
"""
import os, sys, json, time
import numpy as np
from scipy.optimize import least_squares, differential_evolution, minimize
import birth_model_v1 as v1                      # Model, Objective, load, md_of, spectrum; reads the same env at import
from birth_model_v1 import J, COUPLING, STATION, MODE, EPS, SEED, TIME_LIMIT, POP, FLOOR, S_QM, _deg, TAUS, SCALES

K6_JSON = os.environ.get("K6_JSON", "census_bandlimit_v4_r0_e01_k6_v4_rho_nsh32_random0_eps0.01_nang161.json")
K6_K = int(os.environ.get("K6_K", 6)); NSTART = int(os.environ.get("NSTART", 12)); WARM = os.environ.get("WARM", "")
OUT_TAG = os.environ.get("OUT_TAG", f"j{J}_{MODE}_s{SEED}_warm")

def basis(offs, K):
    xa, xb = offs[:, 0], offs[:, 1]; cols = [np.ones(len(offs))]
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            arg = 2 * (ka * xa + kb * xb); cols.append(np.cos(arg)); cols.append(np.sin(arg))
    return np.column_stack(cols)

def target_census(offs, nflat):
    d = json.load(open(K6_JSON))
    r = next((r for r in d["results"] if r.get("K") == K6_K and "coeff" in r), None)
    if r is None: sys.exit(f"no feasible K={K6_K} entry with coefficients in {K6_JSON}")
    c = np.asarray(r["coeff"], float); Phi = basis(offs, K6_K)
    rho = np.concatenate([Phi @ c[ts] for ts in range(nflat)]); rho = np.clip(rho, 0, None); rho /= rho.sum()
    return rho, r.get("price")

def main(path):
    t0 = time.time()
    d, offs, nflat, Qs, Svec, Ecoef, angs = v1.load(path); NV = len(offs); NT = nflat * NV
    tgt = -np.cos(2 * angs * _deg)
    m = v1.Model(offs, nflat); obj = v1.Objective(m, Qs, Svec, Ecoef, tgt)
    rho_t, price_t = target_census(offs, nflat)
    mt_E, mt_S = obj.misses(rho_t)
    print(f"file: {os.path.basename(path)}   atoms {NT} = {nflat} x {NV}")
    print(f"target: K={K6_K} census from {K6_JSON}: price {price_t:.4f} (file), re-measured {v1.md_of(Qs, rho_t)/FLOOR:.4f}; misses band {mt_E:.5f} CHSH {mt_S:.2e}")
    print(f"model: J={J} families, coupling={COUPLING}, station factor={STATION}, {m.P} parameters;  mode={MODE} eps={EPS} seed={SEED} time limit {TIME_LIMIT:.0f}s", flush=True)
    rng = np.random.default_rng(SEED)
    # ---------------- stage 1: fit ----------------
    scale = NT                                       # densities are ~1/NT; residuals ~O(1)
    def resid(x): return (m.census(x) - rho_t) * scale
    starts = [m.lo + rng.random(m.P) * (m.hi - m.lo) for _ in range(NSTART)]
    if WARM and os.path.exists(WARM):
        p = json.load(open(WARM)).get("model", {}).get("params", {})
        if p and all(n in p for n in m.names): starts.insert(0, np.array([p[n] for n in m.names])); print(f"  warm start added from {WARM}")
    best = (np.inf, None); t1 = time.time()
    print(f"\nstage 1: least-squares fit of the form to the K={K6_K} census, {len(starts)} starts", flush=True)
    for i, x0 in enumerate(starts):
        try:
            r = least_squares(resid, x0, bounds=(m.lo, m.hi), max_nfev=1500, xtol=1e-8, ftol=1e-8)
            rel = np.linalg.norm(m.census(r.x) - rho_t) / np.linalg.norm(rho_t)
            fe, fs = obj.misses(m.census(r.x))
            print(f"  start {i}: rel. L2 error {rel:.4f}   fitted census misses band {fe:.5f} CHSH {fs:.5f}   [{time.time()-t1:.0f}s]", flush=True)
            if rel < best[0]: best = (rel, r.x.copy())
        except Exception as e:
            print(f"  start {i}: failed ({e})")
    rel, xfit = best; rho_fit = m.census(xfit); fe, fs = obj.misses(rho_fit); pfit = v1.md_of(Qs, rho_fit) / FLOOR
    print(f"\n  best fit: rel. L2 error {rel:.4f};  fitted census: band miss {fe:.5f}  CHSH miss {fs:.5f}  MD/floor {pfit:.4f}")
    pt = v1.spectrum(rho_t, offs, nflat, NV); pf = v1.spectrum(rho_fit, offs, nflat, NV)
    ft = {(ka, kb): (a, ph) for a, ka, kb, ph in pf}
    print("  leading harmonics, target vs fit  (|F| target, |F| fit, phase target, phase fit):")
    for a, ka, kb, ph in pt[:12]:
        af, phf = ft.get((ka, kb), (0.0, 0.0))
        flag = "" if af > 0.5 * a else "   <-- missing in the form" if af < 0.25 * a else "   <-- weak in the form"
        print(f"     ({ka:+d},{kb:+d})   {a:.3f}  {af:.3f}   {ph:+8.1f}  {phf:+8.1f}{flag}")
    # ---------------- stage 2: refine ----------------
    print(f"\nstage 2: differential evolution ({MODE}) seeded around the fit", flush=True)
    M = POP * m.P; init = np.empty((M, m.P)); init[0] = xfit
    n_pert = int(0.6 * M)
    for i in range(1, n_pert + 1):
        init[i] = np.clip(xfit + rng.normal(0, 0.05, m.P) * (m.hi - m.lo), m.lo, m.hi)
    for i in range(n_pert + 1, M):
        init[i] = m.lo + rng.random(m.P) * (m.hi - m.lo)
    obj.t0 = time.time(); obj.n = 0; obj.best = (np.inf, None); obj.best_miss = (np.inf, None)
    obj(xfit)                                        # make sure the fit itself is on record
    differential_evolution(obj, bounds=list(zip(m.lo, m.hi)), seed=SEED, popsize=POP, maxiter=100000, tol=0, polish=False,
                           init=init, mutation=(0.5, 1.0), recombination=0.7, callback=obj.stop, updating="immediate")
    print(f"  DE: {obj.n} evaluations, best objective {obj.best[0]:.6f}  [{time.time()-t0:.0f}s]", flush=True)
    try:
        minimize(obj, obj.best[1], method="Powell", bounds=list(zip(m.lo, m.hi)), options={"maxiter": 20, "xtol": 1e-4, "ftol": 1e-7})
    except Exception as e:
        print(f"  polish skipped: {e}")
    x = obj.best[1]; rho = m.census(x); missE, missS = obj.misses(rho); price = v1.md_of(Qs, rho) / FLOOR
    feasible = (missE <= EPS + 1e-9) and (missS <= EPS + 1e-9); bm, _ = obj.best_miss
    print(f"\nRESULT  mode={MODE} J={J} coupling={COUPLING} station={STATION} eps={EPS}  (warm start from K={K6_K})")
    print(f"  model census: band miss {missE:.5f}   CHSH miss {missS:.5f}   -> {'FEASIBLE' if feasible else 'not feasible'} at eps={EPS}")
    print(f"  MD/floor = {price:.4f}   (target K={K6_K} census {pfit if False else price_t:.4f}; LP optimum 1.000)")
    print(f"  smallest max-miss seen: {bm:.5f}   (v1 cold start: 0.0256 at J=2; best K=5 census 0.0120)")
    print("\n  parameters:")
    k = 0
    for j in range(J):
        logw, muS, kS, muD, kD, kD2, beta = x[k:k + 7]; k += 7
        print(f"    family {j}: weight {np.exp(logw):.3f}  axis muS {np.degrees(muS):6.1f} deg (kappa_S {kS:5.1f})  push muD {np.degrees(muD):6.1f} deg (kappa_D {kD:5.1f}, kappa_D2 {kD2:+5.1f})  coupling beta {np.degrees(beta) if COUPLING else 0:+6.1f} deg")
    cells = x[k:k + nflat]; W = np.exp(cells - cells.max()); W /= W.sum(); k += nflat
    if nflat == len(TAUS) * len(SCALES):
        print("    cell weights (rows tau, cols scale):")
        for ti, tau in enumerate(TAUS): print(f"      tau={tau:<4}" + "".join(f"  {W[ti*len(SCALES)+si]:.3f}" for si in range(len(SCALES))))
    if STATION: print(f"    station factor: alpha {np.tanh(x[k]):+.3f}  phi {np.degrees(x[k+1]):6.1f} deg")
    pix = v1.spectrum(rho, offs, nflat, NV)
    print("\n  leading harmonics of the final model census  (|F|, k_A, k_B, phi[deg]):")
    for a, ka, kb, ph in pix[:10]: print(f"     {a:.3f}  ({ka:+d},{kb:+d})  {ph:+8.2f}")
    out = dict(d); out.update(rho=rho.tolist(), md=float(price * FLOOR), model=dict(version=2, J=J, coupling=COUPLING, station=STATION, mode=MODE, eps=EPS, seed=SEED,
               target=K6_JSON, target_K=K6_K, fit_rel_error=float(rel), fit_params=dict(zip(m.names, map(float, xfit))), fit_band_miss=float(fe), fit_chsh_miss=float(fs), fit_price=float(pfit),
               params=dict(zip(m.names, map(float, x))), band_miss=float(missE), chsh_miss=float(missS), price=float(price), feasible=bool(feasible),
               min_miss_seen=float(bm), evaluations=obj.n, leading=[(float(a), int(ka), int(kb), float(ph)) for a, ka, kb, ph in pix[:20]]))
    jpath = os.path.join(os.path.dirname(path) or ".", f"birth_model_v2_{OUT_TAG}_{os.path.splitext(os.path.basename(path))[0]}.json")
    json.dump(out, open(jpath, "w"))
    print(f"\nSUMMARY  v2 warm(K={K6_K}) mode={MODE} J={J} c={int(COUPLING)} st={int(STATION)} eps={EPS}: fit rel.err {rel:.4f} (fit miss {max(fe,fs):.5f}) | final band miss {missE:.5f}  CHSH miss {missS:.5f}  price {price:.4f}  {'FEASIBLE' if feasible else 'infeasible'}  min-miss-seen {bm:.5f}   [{time.time()-t0:.0f}s]")
    print(f"wrote {os.path.basename(jpath)}")

if __name__ == "__main__":
    main(sys.argv[1])
