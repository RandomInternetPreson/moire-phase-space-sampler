#!/usr/bin/env python3
"""
dither_sweep.py — Is the offset lattice load-bearing or an aliasing artifact?

The published price table uses quenched offsets on the regular lattice
xi_j = j*pi/NSH. With beat orders up to m=5, n=6, the phases 2m*xi are
commensurate with that lattice (m=4 -> only 4 distinct phases at NSH=16).
This is exactly the sampling-lattice lock seen in the garden. Test:

  grid    : the published lattice (control; must reproduce the table)
  strat   : stratified — one uniform random point per lattice cell
  random  : NSH*NSH offsets drawn i.i.d. uniform on [0,pi)^2
  refine  : grid at NSH_LIST (e.g. 16,24,32) to look for tongue edges

Same NW, dt, TAUS, SCALES, register and tolerance as flight_time.py.
For each run we record MD/floor at every EPS, the number of active
atoms in the optimal rho (weight > 1e-6), and the top atoms, so we can
see whether the optimizer concentrates on a few frozen phases or spreads.

Pre-registered reading (written before running):
  A. strat/random prices == grid prices within the dt uncertainty
     (~0.03 of floor at eps=0.01): lattice is scaffolding, table stands.
  B. off-lattice prices HIGHER: the table is partly an aliasing
     artifact; v2 must report off-lattice prices.
  C. prices jump at specific NSH in the refinement: tongue edges.

Run (Windows cmd):
  set NANG=41
  set NSH_LIST=16,24,32
  set MODES=grid,strat,random
  set SEEDS=0,1,2
  set EPS_LIST=0.10,0.05,0.02,0.01
  set CORES=56
  python dither_sweep.py
"""
import os, sys, json, time
import numpy as np
from scipy.optimize import linprog

# reuse the canonical machinery unchanged
import flight_time as ft

NSH_LIST = [int(x) for x in os.environ.get("NSH_LIST", "16").split(",")]
MODES    = os.environ.get("MODES", "grid,strat,random").split(",")
SEEDS    = [int(x) for x in os.environ.get("SEEDS", "0,1,2").split(",")]
EPS_LIST = [float(x) for x in os.environ.get("EPS_LIST", "0.10,0.05,0.02,0.01").split(",")]
REGISTER = os.environ.get("REGISTER", "paper")   # paper:[2,178]  companion:[6,174]
OUT      = os.environ.get("OUT", "dither_sweep_results.json")

_deg = ft._deg; G = ft.G; RED = ft.RED; DL = ft.DL
settings, signs = ft.settings, ft.signs
SGA, SGB, S_QM, FLOOR = ft.SGA, ft.SGB, ft.S_QM, ft.FLOOR

if REGISTER == "paper":
    ft.angs = np.linspace(2, 178, ft.NANG)
angs = ft.angs


def make_offsets(nsh, mode, seed):
    rng = np.random.default_rng(seed)
    if mode == "grid":
        g = np.arange(nsh) * np.pi / nsh
        return [(xa, xb) for xa in g for xb in g]
    if mode == "gridshift":          # same lattice, phase-shifted by half a cell
        g = (np.arange(nsh) + 0.5) * np.pi / nsh
        return [(xa, xb) for xa in g for xb in g]
    if mode == "strat":
        cell = np.pi / nsh
        return [((i + rng.random()) * cell, (j + rng.random()) * cell)
                for i in range(nsh) for j in range(nsh)]
    if mode == "random":
        pts = rng.random((nsh * nsh, 2)) * np.pi
        return [tuple(p) for p in pts]
    raise ValueError(mode)


def build_banks_offsets(offsets, cores):
    """Same as ft.build_banks but over an explicit offset list."""
    from multiprocessing import Pool
    NV = len(offsets)
    jobs, meta = [], []
    for si, sc in enumerate(ft.SCALES):
        w = ft.W0 * sc
        for pi, (a, b) in enumerate(settings):
            for (xa, xb) in offsets:
                jobs.append((w, xa, xb, a * _deg, b * _deg)); meta.append(("pair", si, pi))
        for ai, tdeg in enumerate(angs):
            for (xa, xb) in offsets:
                jobs.append((w, xa, xb, 0.0, tdeg * _deg)); meta.append(("ang", si, ai))
    print(f"  evolving {len(jobs)} trajectories ({len(ft.TAUS)} tau snapshots)...", flush=True)
    t0 = time.time()
    with Pool(cores) as pool:
        outs = pool.map(ft.evolve_component, jobs, chunksize=4)
    print(f"  done in {time.time()-t0:.0f}s", flush=True)
    BANKS = [[None] * len(ft.SCALES) for _ in ft.TAUS]
    for ti in range(len(ft.TAUS)):
        for si in range(len(ft.SCALES)):
            Qs = [np.zeros((NV, G)) for _ in range(4)]
            Ec = np.zeros((ft.NANG, NV))
            cnt_p = [0] * 4; cnt_a = [0] * ft.NANG
            for (kind, s2, idx), o in zip(meta, outs):
                if s2 != si: continue
                if kind == "pair":
                    Qs[idx][cnt_p[idx]] = o[ti]; cnt_p[idx] += 1
                else:
                    tr = angs[idx] * _deg
                    ob = np.sign(np.cos(2 * (ft.lam_c + np.pi / 2 - tr)))
                    Ec[idx, cnt_a[idx]] = np.sum(o[ti] * SGA[0] * ob) * DL
                    cnt_a[idx] += 1
            Svec = np.zeros(NV)
            for i, sg in enumerate(signs):
                Svec += sg * (Qs[i] * (SGA[i] * SGB[i])[None, :]).sum(axis=1) * DL
            BANKS[ti][si] = (Qs, Svec, Ec)
    return BANKS


def min_md_full(banklist, eps):
    """ft.min_md, but also returns the optimal rho so we can inspect support."""
    Qs = [np.vstack([b[0][p] for b in banklist]) for p in range(4)]
    Svec = np.concatenate([b[1] for b in banklist])
    Ecoef = np.hstack([b[2] for b in banklist])
    NV = len(Svec)
    Qc = [Q.reshape(NV, RED, -1).mean(axis=2) for Q in Qs]
    cmb = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    nv = NV + 1 + len(cmb) * RED
    rows, bub = [], []
    binw = np.pi / RED
    for ci, (i, j) in enumerate(cmb):
        Dif = (Qc[i] - Qc[j]) * binw
        for l in range(RED):
            r = np.zeros(nv); r[:NV] = Dif[:, l]; r[NV + 1 + ci * RED + l] = -1
            rows.append(r); bub.append(0.0)
            r = np.zeros(nv); r[:NV] = -Dif[:, l]; r[NV + 1 + ci * RED + l] = -1
            rows.append(r); bub.append(0.0)
    for ci in range(len(cmb)):
        r = np.zeros(nv); r[NV + 1 + ci * RED:NV + 1 + (ci + 1) * RED] = 0.5
        r[NV] = -1; rows.append(r); bub.append(0.0)
    tgt = -np.cos(2 * angs * _deg)
    for ai in range(ft.NANG):
        r = np.zeros(nv); r[:NV] = Ecoef[ai]; rows.append(r); bub.append(tgt[ai] + eps)
        r = np.zeros(nv); r[:NV] = -Ecoef[ai]; rows.append(r); bub.append(-tgt[ai] + eps)
    Aub = np.array(rows); bv = np.array(bub)
    Aeq = np.zeros((2, nv)); Aeq[0, :NV] = 1; Aeq[1, :NV] = Svec
    cobj = np.zeros(nv); cobj[NV] = 1
    res = linprog(cobj, A_ub=Aub, b_ub=bv, A_eq=Aeq, b_eq=[1.0, -S_QM],
                  bounds=[(0, None)] * nv, method="highs")
    if res.status != 0:
        return None, None
    rho = res.x[:NV]
    qm = [Qs[p].T @ rho for p in range(4)]
    md = max(0.5 * np.sum(np.abs(qm[a] - qm[b])) * DL
             for a in range(4) for b in range(a + 1, 4))
    min_md_full.last_residual = (Ecoef @ rho - tgt).tolist()   # E_model - E_QM per register angle
    return md, rho


def main():
    cores = int(os.environ.get("CORES", os.cpu_count() or 8))
    results = {"meta": {"NANG": ft.NANG, "NW": ft.NW, "TAUS": ft.TAUS,
                        "SCALES": ft.SCALES, "register": list(map(float, [angs[0], angs[-1]])),
                        "EPS_LIST": EPS_LIST, "floor": FLOOR}, "runs": []}
    if os.path.exists(OUT):
        results = json.load(open(OUT))
        print(f"resuming: {len(results['runs'])} runs already in {OUT}")
    bad = [r for r in results["runs"] if r.get("nang", ft.NANG) != ft.NANG]
    if bad:
        sys.exit(f"ERROR: {OUT} contains runs made at NANG={bad[0].get('nang')} but this "
                 f"session has NANG={ft.NANG}. Set NANG to match or use a different OUT.")
    done = {(r["nsh"], r["mode"], r["seed"]) for r in results["runs"]}

    for nsh in NSH_LIST:
        for mode in MODES:
            seeds = [0] if mode in ("grid", "gridshift") else SEEDS
            for seed in seeds:
                if (nsh, mode, seed) in done: continue
                print(f"\n=== NSH={nsh}  mode={mode}  seed={seed}  NANG={ft.NANG} ===", flush=True)
                offsets = make_offsets(nsh, mode, seed)
                banks = build_banks_offsets(offsets, cores)
                flat = [b for row in banks for b in row]
                nv_off = len(offsets); ntau = len(ft.TAUS); nsc = len(ft.SCALES)
                run = {"nsh": nsh, "mode": mode, "seed": seed, "nang": ft.NANG,
                       "price": {}, "support": {}}
                for eps in EPS_LIST:
                    md, rho = min_md_full(flat, eps)
                    if md is None:
                        run["price"][str(eps)] = None; print(f"  eps={eps}: infeasible"); continue
                    ratio = md / FLOOR
                    act = np.where(rho > 1e-6)[0]
                    # decode atom index -> (tau_idx, scale_idx, offset_idx); pooled order is tau-major
                    top = np.argsort(rho)[::-1][:12]
                    atoms = []
                    for k in top:
                        ti, rem = divmod(int(k), nsc * nv_off)
                        si, oi = divmod(rem, nv_off)
                        xa, xb = offsets[oi]
                        atoms.append({"w": float(rho[k]), "tau": ft.TAUS[ti],
                                      "scale": ft.SCALES[si],
                                      "xa_deg": float(xa / _deg), "xb_deg": float(xb / _deg)})
                    run["price"][str(eps)] = ratio
                    run.setdefault("residual", {})[str(eps)] = {"angles_deg": angs.tolist(),
                                                               "E_model_minus_E_QM": min_md_full.last_residual}
                    run["support"][str(eps)] = {"n_active": int(len(act)),
                                                "n_total": int(len(rho)),
                                                "top_weight_share": float(rho[top].sum()),
                                                "top": atoms}
                    print(f"  eps={eps:<5} MD/floor={ratio:7.3f}   active atoms {len(act)}/{len(rho)}"
                          f"   top-12 carry {rho[top].sum():.2f}", flush=True)
                results["runs"].append(run)
                json.dump(results, open(OUT, "w"), indent=1)

    # summary table
    print("\n\nSUMMARY  (MD/floor)")
    hdr = f"{'NSH':>4} {'mode':>7} {'seed':>4} " + " ".join(f"{'eps='+str(e):>10}" for e in EPS_LIST)
    print(hdr)
    for r in results["runs"]:
        row = f"{r['nsh']:>4} {r['mode']:>7} {r['seed']:>4} "
        row += " ".join(f"{(r['price'][str(e)] if r['price'][str(e)] is not None else float('nan')):>10.3f}"
                        for e in EPS_LIST)
        print(row)
    print(f"\npublished (grid, NSH=16, [2,178]): 1.000 1.015 1.235 1.375")
    print("A: off-lattice == grid (within ~0.03)  -> lattice is scaffolding")
    print("B: off-lattice > grid                   -> aliasing in the table; v2 needed")
    print("C: jumps across NSH                     -> tongue edges")


if __name__ == "__main__":
    main()
