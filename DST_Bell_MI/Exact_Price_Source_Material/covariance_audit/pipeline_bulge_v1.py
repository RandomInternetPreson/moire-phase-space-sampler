#!/usr/bin/env python3
"""
pipeline_bulge_v1.py -- the first test of the BULGE: does the shape of the detection rule matter?

In the model the bulge is the detection rule A = sign(cos 2(lambda - a)): a step, with every odd harmonic of lambda and no
structure beyond the step.  The banks store the 161 correlator angles already contracted with that sign, so no computation
so far could vary it.  This pipeline rebuilds the lattice trajectories once and stores the full lambda-HISTOGRAM of every
(tau, scale, angle, offset) -- 15 x 165 x NV x 256 counts, uint16 -- so that any detection profile D_A, D_B is a
re-contraction (seconds) followed by the usual LPs (minutes to an hour each).  Then it runs census_profile_v1.py on a list
of profiles, in parallel, and prints one summary table.

Profiles (census_profile_v1.py): sign (the paper's step; the control -- must reproduce the banks bit-for-bit and the known
K-curve 1.700/1.429/1.258/1.112/1.000) | smooth (D = cos 2theta: the fully smooth bulge) | tanh:B (blunted step,
tanh(B cos 2theta)/tanh B) | notch:J:E (step with notches of depth E at frequency 4J in theta: a blunted-SPIKY bulge, same on
both stations) | anotch:J:E (the twin carries the anti-version: Bob's notches sit where Alice's do not).

Stages (each skipped if its output exists):  1. offsets + streamed build -> hist_nsh<N>_<mode><seed>_nang<nang>_nw<nw>.npz
  2. for each profile: census_profile_v1.py -> bulge_<ptag>.txt + census_profile_v1_<ptag>.json   3. SUMMARY table.
Console: pipeline_bulge_console.txt.  Run through the launcher (RUN_BULGE.bat); ~4 h build, then the sweep.

Usage (from covariance_audit):
  python pipeline_bulge_v1.py --nsh 32 --mode grid --cores 56 --klist 6,8,10,12,14 --maxpar 14      (default profile list: 14 bulges)
"""
import os, sys, json, time, argparse, subprocess, datetime
ap = argparse.ArgumentParser()
ap.add_argument("--nsh", type=int, default=32); ap.add_argument("--mode", default="grid"); ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--nang", type=int, default=161); ap.add_argument("--nw", type=int, default=2048); ap.add_argument("--eps", type=float, default=0.01)
ap.add_argument("--cores", type=int, default=0); ap.add_argument("--klist", default="6,8,10,12,14"); ap.add_argument("--full", type=int, default=1)
ap.add_argument("--register", default="paper"); ap.add_argument("--chunk", type=int, default=16)
ap.add_argument("--profiles", default="sign,smooth,tanh:1,tanh:2,tanh:4,tanh:8,notch:3:0.25,notch:5:0.25,notch:7:0.25,anotch:3:0.25,anotch:5:0.25,anotch:7:0.25,notch:5:0.1,anotch:5:0.1")
ap.add_argument("--maxpar", type=int, default=14); ap.add_argument("--build_only", type=int, default=0)
A = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, HERE)
os.environ["NANG"] = str(A.nang); os.environ["NW"] = str(A.nw); os.environ["REGISTER"] = A.register    # before importing flight_time
import numpy as np
import flight_time as ft
import dither_sweep_rho as dsr
from census_profile_v1 import ptag as _ptag

def main():
    class Tee:
        def __init__(self, path): self.f = open(path, "a"); self.o = sys.__stdout__
        def write(self, s): self.o.write(s); self.f.write(s)
        def flush(self): self.o.flush(); self.f.flush()
    sys.stdout = Tee("pipeline_bulge_console.txt")
    def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    def log(s): print(f"[{now()}] {s}", flush=True)
    cores = A.cores or int(os.environ.get("CORES", 0)) or os.cpu_count() or 8
    NV = A.nsh * A.nsh; ntau, nsc = len(ft.TAUS), len(ft.SCALES); G, DL, _deg = ft.G, ft.DL, ft._deg
    angs = dsr.angs; settings = dsr.settings; nang = len(angs)
    histf = f"hist_nsh{A.nsh}_{A.mode}{A.seed}_nang{A.nang}_nw{A.nw}.npz"
    log(f"pipeline_bulge_v1: nsh={A.nsh} ({NV} offsets) mode={A.mode} seed={A.seed} nang={A.nang} nw={A.nw} eps={A.eps} cores={cores} register={A.register}")
    log(f"shared-compute job {os.environ.get('SHARED_COMPUTE_JOB_ID','(none)')}, reserved {os.environ.get('SHARED_COMPUTE_MEM_GIB','?')} GiB")
    offsets = dsr.make_offsets(A.nsh, A.mode, A.seed)

    # ---------- 1. the histogram store ----------
    if os.path.exists(histf):
        log(f"histograms exist: {histf} (skipping build)")
    else:
        import multiprocessing
        ctx = multiprocessing.get_context("spawn")
        jobs, meta = [], []
        for si, sc in enumerate(ft.SCALES):
            w = ft.W0 * sc
            for pi, (a, b) in enumerate(settings):
                for (xa, xb) in offsets: jobs.append((w, xa, xb, a * _deg, b * _deg)); meta.append(("pair", si, pi))
            for ai, tdeg in enumerate(angs):
                for (xa, xb) in offsets: jobs.append((w, xa, xb, 0.0, tdeg * _deg)); meta.append(("ang", si, ai))
        Hpair = np.zeros((ntau, nsc, 4, NV, G), np.uint16); Hang = np.zeros((ntau, nsc, nang, NV, G), np.uint16)
        cnt_p = np.zeros((nsc, 4), int); cnt_a = np.zeros((nsc, nang), int)
        log(f"evolving {len(jobs)} trajectories on {cores} cores ({ntau} tau snapshots each); storing counts "
            f"({Hpair.nbytes/1e9:.2f} + {Hang.nbytes/1e9:.2f} GB in memory as uint16)")
        t0 = time.time(); nxt = 0
        with ctx.Pool(cores) as pool:
            for k, o in enumerate(pool.imap(ft.evolve_component, jobs, chunksize=A.chunk)):
                kind, si, idx = meta[k]
                if kind == "pair":
                    c = cnt_p[si, idx]; cnt_p[si, idx] += 1
                    for ti in range(ntau): Hpair[ti, si, idx, c] = np.rint(o[ti] * A.nw * DL)
                else:
                    c = cnt_a[si, idx]; cnt_a[si, idx] += 1
                    for ti in range(ntau): Hang[ti, si, idx, c] = np.rint(o[ti] * A.nw * DL)
                if k + 1 >= nxt:
                    el = time.time() - t0; rate = (k + 1) / max(el, 1e-9)
                    log(f"  {k+1}/{len(jobs)} trajectories  ({el/3600:.2f} h elapsed, ~{(len(jobs)-k-1)/rate/3600:.2f} h left)"); nxt += len(jobs) // 40
        log(f"build done in {(time.time()-t0)/3600:.2f} h")
        assert cnt_p.min() == NV and cnt_a.min() == NV
        # every histogram must hold exactly NW counts (the density was h / (NW * DL))
        bad = int(np.sum(Hpair.reshape(-1, G).sum(axis=1, dtype=np.int64) != A.nw) + np.sum(Hang.reshape(-1, G).sum(axis=1, dtype=np.int64) != A.nw))
        log(f"count check: {bad} histograms with a total != {A.nw} (must be 0)")
        if bad: sys.exit("ABORT: histogram counts do not reconstruct the densities")
        t0 = time.time()
        np.savez_compressed(histf, Hpair=Hpair.reshape(ntau * nsc, 4, NV, G), Hang=Hang.reshape(ntau * nsc, nang, NV, G),
                            offsets=np.array(offsets, float), angs=np.array(angs, float), taus=np.array(ft.TAUS), scales=np.array(ft.SCALES),
                            nw=np.array(A.nw), G=np.array(G), ordering=np.array("tau_major"), register=np.array(A.register))
        log(f"histograms saved to {histf} ({os.path.getsize(histf)/1e9:.2f} GB on disk, {(time.time()-t0)/60:.1f} min)")
        del Hpair, Hang
    if A.build_only: log("build only -- done"); return

    # ---------- 2. the profile sweep ----------
    profs = [p for p in A.profiles.split(",") if p]
    pending = []; running = {}
    for p in profs:
        tag = _ptag(p); logf = f"bulge_{tag}.txt"
        if os.path.exists(logf) and "SUMMARY" in open(logf).read(): log(f"{p}: {logf} already complete (skipping)"); continue
        pending.append((p, tag, logf))
    t0 = time.time()
    while pending or running:
        while pending and len(running) < A.maxpar:
            p, tag, logf = pending.pop(0)
            env = dict(os.environ); env.update(REGISTER=A.register)
            fh = open(logf, "w")
            cmd = [sys.executable, "census_profile_v1.py", histf, p, "--klist", A.klist, "--eps", str(A.eps), "--full", str(A.full)]
            running[tag] = (subprocess.Popen(cmd, env=env, stdout=fh, stderr=subprocess.STDOUT), fh, logf, p)
            log(f"started {p} -> {logf}")
        time.sleep(30)
        for tag in list(running):
            pr, fh, logf, p = running[tag]
            if pr.poll() is not None:
                fh.close(); del running[tag]; log(f"{p} finished exit={pr.returncode} after {(time.time()-t0)/60:.0f} min")

    # ---------- 3. summary ----------
    KL = [int(x) for x in A.klist.split(",")]
    cols = (["full"] if A.full else []) + [f"K={K}" for K in KL]
    print(f"\nSUMMARY  bulge sweep on {histf}: MD/floor at eps={A.eps} for each detection profile (full census, then band-limited K)")
    print(f"{'profile':<16} {'S_max/2rt2':>10} {'miss':>8} " + " ".join(f"{c:>8}" for c in cols) + f" {'K'+str(KL[0])+'/full':>8}" + "   <A>max  <B>max  (raw vertex; symmetrized census)")
    for p in profs:
        tag = _ptag(p); jf = f"census_profile_v1_{tag}.json"
        if not os.path.exists(jf): print(f"{p:<16} (no json)"); continue
        d = json.load(open(jf)); r = {x["label"]: x for x in d["results"]}
        def cell(lbl):
            x = r.get(lbl)
            if x is None: return "?"
            if "price" in x: return f"{x['price']:.4f}"
            if "min_miss" in x: return f"m{x['min_miss']:.4f}"
            return "-"
        mg = d.get("marginals", {}); ms = d.get("marginals_sym")
        prem = "-"
        if "price" in r.get("full", {}) and "price" in r.get(f"K={KL[0]}", {}): prem = f"{r[f'K={KL[0]}']['price']/r['full']['price']:.3f}"
        print(f"{p:<16} {d['smax']/dsr.S_QM:10.6f} {d.get('miss_relaxed', float('nan')):8.4f} " + " ".join(f"{cell(c):>8}" for c in cols)
              + f" {prem:>8}" + f"   {mg.get('A_max', float('nan')):.4f}  {mg.get('B_max', float('nan')):.4f}" + (f"  sym {ms['A_max']:.0e} (MD/floor {ms['md_over_floor']:.4f})" if ms else ""))
    print("reading: 'sign' is the paper's rule and must give the known lattice K-curve (1.000 full; 1.700/1.429/1.258/1.112/1.000 at K=6..14).")
    print("         'full' = the floor price of the bulge itself (blunting is expected to raise it); K/full = the smoothness premium under that bulge.")
    print("         if a shaped bulge brings the ABSOLUTE K=6 price toward 1.000 the grain can live in the bulge; if not, the grain is the birth's.")
    print("         'm' = this bulge cannot reach the singlet to eps: minimum miss shown.  Pre-registered predictions: bulge_test_preregistration note (31 Oct).")
    log("pipeline_bulge_v1 finished")

if __name__ == "__main__":
    main()
