#!/usr/bin/env python3
"""
pipeline_terrain_v1.py -- the grain test: does the census's bandwidth track the TERRAIN's finest beat, or is it the photon's?

Finding (28 Oct): the floor-attaining census needs harmonics 7..14 in the offsets, contiguous (one ripple buys < 2%; a hole at
7-9 costs 13%); the greedy terrain's station indices run to m = 5, n = 6.  Conjecture: the stamp resolves the finest fence it
meets, so the needed bandwidth is set by the terrain's highest order.  Test: cut the terrain to max(m, n) <= MMAX, rebuild the
lattice banks, solve the LP, run the K-curve.  Pre-registered: if the grain is the mesh's, K_needed falls with MMAX
(~7-8 at MMAX=3, ~5 at MMAX=2); if it is the photon's, it stays ~14.

Mechanics: a subfolder terrain_<tag>/ gets a FILTERED greedy_results.json and copies of flight_time.py, dither_sweep_rho.py,
census_bandlimit_v4.py, census_structure_v2.py, census_sumdiff_v1.py; the pipeline chdirs there before importing
flight_time (which reads greedy_results.json from the working directory), so the parent AND the spawned workers see the
filtered terrain.  The pool verifies the workers' term count before building.  Stages as pipeline_4096_v2 (skip-if-exists).
Outputs in terrain_<tag>/: banks_nsh32_grid0_*.npz, rho_nsh32_grid0_*.json, bl_n32_k<K>.txt, structure transcripts,
pipeline_terrain_<tag>_console.txt.

Usage (from covariance_audit, through the launcher -- see RUN_TERRAIN.bat):
  python pipeline_terrain_v1.py --mmax 3 --tag ord3 --nsh 32 --mode grid --klist 3,4,5,6,8,10,12,14 --cores 56
"""
import os, sys, json, time, argparse, subprocess, datetime
ap = argparse.ArgumentParser()
ap.add_argument("--nsh", type=int, default=64); ap.add_argument("--mode", default="random"); ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--nang", type=int, default=161); ap.add_argument("--nw", type=int, default=2048); ap.add_argument("--eps", type=float, default=0.01)
ap.add_argument("--cores", type=int, default=0); ap.add_argument("--klist", default="6,8,10,12,14,16,18,20")
ap.add_argument("--register", default="paper"); ap.add_argument("--chunk", type=int, default=16)
ap.add_argument("--mmax", type=int, required=True); ap.add_argument("--tag", required=True)
A = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
SUB = os.path.join(HERE, f"terrain_{A.tag}")
_full = json.load(open(os.path.join(HERE, "greedy_results.json")))
_keep = [i for i, (m, n, sg) in enumerate(_full["terms"]) if max(m, n) <= A.mmax]
_filt = dict(_full); _filt["terms"] = [_full["terms"][i] for i in _keep]; _filt["w"] = [_full["w"][i] for i in _keep]
_filt["filter"] = dict(rule=f"max(m,n) <= {A.mmax}", kept=len(_keep), of=len(_full["terms"]), weight_kept=float(sum(_filt["w"]) / sum(_full["w"])))
if __name__ == "__main__":                      # only the parent writes the folder; spawned workers just use it
    import shutil
    os.makedirs(SUB, exist_ok=True)
    json.dump(_filt, open(os.path.join(SUB, "greedy_results.json"), "w"), indent=1)
    for f in ("flight_time.py", "dither_sweep_rho.py", "census_bandlimit_v4.py", "census_structure_v2.py", "census_sumdiff_v1.py"):
        shutil.copy(os.path.join(HERE, f), os.path.join(SUB, f))
os.chdir(SUB); sys.path.insert(0, SUB)
# flight_time reads these from the environment at import time -- set them BEFORE importing
os.environ["NANG"] = str(A.nang); os.environ["NW"] = str(A.nw); os.environ["REGISTER"] = A.register
import numpy as np
import flight_time as ft
def _terms_count(): return len(ft.TERMS)          # run in a worker: must equal the parent's
import dither_sweep_rho as dsr            # make_offsets, min_md_full, settings, signs, SGA, SGB, angs (register applied)

def main():
    class Tee:                                   # everything printed also goes to pipeline_4096_console.txt (append: resubmits accumulate)
        def __init__(self, path): self.f = open(path, "a"); self.o = sys.__stdout__
        def write(self, s): self.o.write(s); self.f.write(s)
        def flush(self): self.o.flush(); self.f.flush()
    sys.stdout = Tee(f"pipeline_terrain_{A.tag}_console.txt")
    def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    def log(s): print(f"[{now()}] {s}", flush=True)
    cores = A.cores or int(os.environ.get("CORES", 0)) or os.cpu_count() or 8
    NV = A.nsh * A.nsh; ntau, nsc = len(ft.TAUS), len(ft.SCALES); G, DL, _deg = ft.G, ft.DL, ft._deg
    angs = dsr.angs; settings, SGA = dsr.settings, dsr.SGA
    bankf = f"banks_nsh{A.nsh}_{A.mode}{A.seed}_nang{A.nang}_nw{A.nw}.npz"
    rhof = f"rho_nsh{A.nsh}_{A.mode}{A.seed}_eps{A.eps}_nang{A.nang}.json"
    log(f"pipeline_terrain_v1 [{A.tag}]: terrain cut max(m,n) <= {A.mmax}: {_filt['filter']['kept']} of {_filt['filter']['of']} terms, {100*_filt['filter']['weight_kept']:.1f}% of the weight; parent sees {len(ft.TERMS)} terms")
    log(f"nsh={A.nsh} ({NV} offsets) mode={A.mode} seed={A.seed} nang={A.nang} nw={A.nw} eps={A.eps} cores={cores} register={A.register}  folder {SUB}")
    log(f"shared-compute job {os.environ.get('SHARED_COMPUTE_JOB_ID','(none)')}, reserved {os.environ.get('SHARED_COMPUTE_MEM_GIB','?')} GiB")

    # ---------- 1-2. offsets and streamed banks ----------
    offsets = dsr.make_offsets(A.nsh, A.mode, A.seed)
    if os.path.exists(bankf):
        log(f"banks exist: {bankf} (skipping build)")
    else:
        import multiprocessing
        ctx = multiprocessing.get_context("spawn")
        jobs, meta = [], []
        for si, sc in enumerate(ft.SCALES):
            w = ft.W0 * sc
            for pi, (a, b) in enumerate(settings):
                for (xa, xb) in offsets:
                    jobs.append((w, xa, xb, a * _deg, b * _deg)); meta.append(("pair", si, pi))
            for ai, tdeg in enumerate(angs):
                for (xa, xb) in offsets:
                    jobs.append((w, xa, xb, 0.0, tdeg * _deg)); meta.append(("ang", si, ai))
        # bank arrays, filled as results stream in (same placement as build_banks_offsets)
        Qs = np.zeros((ntau, nsc, 4, NV, G)); Ec = np.zeros((ntau, nsc, len(angs), NV))
        cnt_p = np.zeros((nsc, 4), int); cnt_a = np.zeros((nsc, len(angs)), int)
        OB = [np.sign(np.cos(2 * (ft.lam_c + np.pi / 2 - tdeg * _deg))) for tdeg in angs]
        log(f"evolving {len(jobs)} trajectories on {cores} cores ({ntau} tau snapshots each), streaming into the banks")
        t0 = time.time(); nxt = 0
        with ctx.Pool(cores) as pool:
            wc = pool.apply(_terms_count)
            log(f"worker check: workers see {wc} terms (parent {len(ft.TERMS)})")
            if wc != len(ft.TERMS): sys.exit("ABORT: workers are not using the filtered terrain")
            for k, o in enumerate(pool.imap(ft.evolve_component, jobs, chunksize=A.chunk)):
                kind, si, idx = meta[k]
                if kind == "pair":
                    c = cnt_p[si, idx]
                    for ti in range(ntau): Qs[ti, si, idx, c] = o[ti]
                    cnt_p[si, idx] += 1
                else:
                    c = cnt_a[si, idx]
                    for ti in range(ntau): Ec[ti, si, idx, c] = np.sum(o[ti] * SGA[0] * OB[idx]) * DL
                    cnt_a[si, idx] += 1
                if k + 1 >= nxt:
                    el = time.time() - t0; rate = (k + 1) / max(el, 1e-9)
                    log(f"  {k+1}/{len(jobs)} trajectories  ({el/3600:.2f} h elapsed, ~{(len(jobs)-k-1)/rate/3600:.2f} h left)"); nxt += len(jobs) // 40
        log(f"build done in {(time.time()-t0)/3600:.2f} h")
        assert cnt_p.min() == NV and cnt_a.min() == NV
        arrs = {"nflat": np.array(ntau * nsc)}; k = 0
        for ti in range(ntau):
            for si in range(nsc):
                Svec = np.zeros(NV)
                for p, sg in enumerate(dsr.signs):
                    Svec += sg * (Qs[ti, si, p] * (dsr.SGA[p] * dsr.SGB[p])[None, :]).sum(axis=1) * DL
                for p in range(4): arrs[f"q{k}_{p}"] = Qs[ti, si, p]
                arrs[f"s{k}"] = Svec; arrs[f"e{k}"] = Ec[ti, si]; k += 1
        np.savez_compressed(bankf, **arrs); log(f"banks checkpointed to {bankf}")
        del Qs, Ec, arrs

    # ---------- 3. the LP ----------
    if os.path.exists(rhof):
        log(f"census exists: {rhof} (skipping LP)"); d = json.load(open(rhof)); log(f"  price {d['ratio']:.6f}")
    else:
        z = np.load(bankf); nflat = int(z["nflat"])
        flat = [([z[f"q{k}_{p}"] for p in range(4)], z[f"s{k}"], z[f"e{k}"]) for k in range(nflat)]
        log(f"LP: min_md_full over {nflat * NV} atoms at eps={A.eps} (HiGHS)"); t0 = time.time()
        md, rho = dsr.min_md_full(flat, A.eps)
        if md is None: sys.exit("LP infeasible -- stop")
        ratio = md / dsr.FLOOR; act = int(np.sum(rho > 1e-6))
        json.dump(dict(rho=rho.tolist(), md=float(md), ratio=float(ratio), offsets=[[float(x), float(y)] for x, y in offsets],
                       ordering="tau_major", nsh=A.nsh, mode=A.mode, seed=A.seed, eps=A.eps, nang=A.nang, nw=A.nw), open(rhof, "w"))
        log(f"LP done in {(time.time()-t0)/60:.0f} min: MD/floor = {ratio:.6f}, active atoms {act}/{len(rho)}; saved {rhof}")

    # ---------- 4. the K-curve, in parallel ----------
    KL = [int(x) for x in A.klist.split(",")]
    procs = {}
    for K in KL:
        tag = f"n{A.nsh}_k{K}"; logf = f"bl_{tag}.txt"
        if os.path.exists(logf) and "SUMMARY" in open(logf).read():
            log(f"K={K}: {logf} already complete (skipping)"); continue
        env = dict(os.environ); env.update(K_LIST=str(K), EPS_LIST=str(A.eps), OUT_TAG=tag, TIME_LIMIT="0", METHOD="highs-ipm", METHOD2="highs-ds", MISS_ONLY="0", REGISTER=A.register)
        fh = open(logf, "w")
        procs[K] = (subprocess.Popen([sys.executable, "census_bandlimit_v4.py", rhof], env=env, stdout=fh, stderr=subprocess.STDOUT), fh, logf)
        log(f"started K={K} -> {logf}")
    t0 = time.time(); done = set()
    while len(done) < len(procs):
        time.sleep(30)
        for K, (p, fh, logf) in procs.items():
            if K not in done and p.poll() is not None:
                done.add(K); fh.close(); log(f"K={K} finished exit={p.returncode} after {(time.time()-t0)/60:.0f} min")

    # ---------- 5. structure ----------
    for script, out in (("census_structure_v2.py", f"census_structure_v2_nsh{A.nsh}_{A.mode}{A.seed}.txt"), ("census_sumdiff_v1.py", f"census_sumdiff_v1_nsh{A.nsh}_{A.mode}{A.seed}.txt")):
        if os.path.exists(script) and not os.path.exists(out):
            with open(out, "w") as fh: subprocess.run([sys.executable, script, rhof], stdout=fh, stderr=subprocess.STDOUT)
            log(f"{script} -> {out}")
            for src in ("census_structure_v2_out.json", "census_sumdiff_v1_out.json"):
                dst = src.replace("_out.json", f"_out_nsh{A.nsh}.json")
                if os.path.exists(src) and not os.path.exists(dst): os.replace(src, dst)

    # ---------- 6. summary ----------
    print("\nSUMMARY  K-curve of the band-limited census, MD/floor at eps=%s" % A.eps)
    print(f"{'K':>4} {'terrain '+A.tag+' (max(m,n)<='+str(A.mmax)+')':>30} {'full terrain, lattice':>22}   (prices are MD / Hall floor; this terrain's own unrestricted price is printed above)")
    ref = {3: "8.5% miss", 4: "1.17% miss", 5: "1.055% miss", 6: 1.700, 8: 1.429, 10: 1.258, 12: 1.112, 14: 1.000}
    for K in KL:
        lines = open(f"bl_n{A.nsh}_k{K}.txt").read().splitlines() if os.path.exists(f"bl_n{A.nsh}_k{K}.txt") else []
        i = next((j for j, l in enumerate(lines) if l.startswith("SUMMARY")), None)
        val = lines[i + 2].split()[2] if i is not None and i + 2 < len(lines) else "?"
        print(f"{K:>4} {val:>26} {ref.get(K, '-'):>22}")
    print("reading: the K at which this terrain's curve reaches its own unrestricted price is K_needed; if it falls with MMAX the grain is the mesh's, if it stays ~14 it is the photon's.")
    log(f"pipeline_terrain_v1 [{A.tag}] finished")


if __name__ == "__main__":
    main()
