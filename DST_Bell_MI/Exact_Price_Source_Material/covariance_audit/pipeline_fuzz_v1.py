#!/usr/bin/env python3
"""
pipeline_fuzz_v1.py -- the disc's edge: build the PARTICLE store, then run the joint-birth LPs (census_fuzz_v1.py).

The histogram store (pipeline_bulge_v1) forgot which initial phase each particle started from.  To test a non-uniform birth
distribution of lambda_0 -- the disc's rim: knife-edge, blunt, multimodal, or correlated with the approach phases -- every
initial particle's final bin is kept: Bpair (cells, 4, NV, NW) and Bang (cells, nang, NV, NW), uint8 (G = 256 bins), 5.2 GB in
memory for the 32x32 lattice; any birth fuzz is then a reweighting, and the LP can carry the fuzz as free harmonics in lambda_0.
Same trajectories as job 403 (checked: the bins' histograms must equal the stored histograms bit for bit when hist_*.npz is present).

Stages (each skipped if its output exists):  1. particle store part_nsh<N>_<mode><seed>_nang<nang>_nw<nw>.npz
  2. contraction with Lmax harmonics (census_fuzz_v1.py, cached)   3. the LPs, (L, K, scope) grid, MAXPAR at a time -> fuzz_L<L>_K<K>_<scope>.txt
  4. SUMMARY table.   Console: pipeline_fuzz_console.txt.   Run through the launcher (RUN_FUZZ.bat): ~7.5 h build + LPs.
Usage (from covariance_audit):
  python pipeline_fuzz_v1.py --nsh 32 --mode grid --cores 56 --Llist 0,1,2 --klist 6,8,10 --scopes cell,atom --maxpar 12
"""
import os, sys, json, time, argparse, subprocess, datetime
ap = argparse.ArgumentParser()
ap.add_argument("--nsh", type=int, default=32); ap.add_argument("--mode", default="grid"); ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--nang", type=int, default=161); ap.add_argument("--nw", type=int, default=2048); ap.add_argument("--eps", type=float, default=0.01)
ap.add_argument("--cores", type=int, default=0); ap.add_argument("--register", default="paper"); ap.add_argument("--chunk", type=int, default=16)
ap.add_argument("--Llist", default="0,1,2"); ap.add_argument("--Lmax", type=int, default=3); ap.add_argument("--klist", default="6,8,10"); ap.add_argument("--scopes", default="cell,atom")
ap.add_argument("--maxpar", type=int, default=12); ap.add_argument("--build_only", type=int, default=0)
A = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, HERE)
os.environ["NANG"] = str(A.nang); os.environ["NW"] = str(A.nw); os.environ["REGISTER"] = A.register
import numpy as np
import flight_time as ft
import dither_sweep_rho as dsr

def evolve_bins(args):
    """ft.evolve_component's loop, returning each particle's lambda-bin at every tau snapshot (ntau, NW) uint8 instead of the histogram."""
    w, xa, xb, a, b = args
    lam = np.linspace(0, np.pi, ft.NW, endpoint=False).copy(); dt = 0.004; out = []; t = 0.0
    for tau in ft.TAUS:
        steps = int(round((tau - t) / dt))
        for _ in range(steps): lam += ft.force(lam, w, xa, xb, a, b) * dt
        t = tau
        out.append((((lam % np.pi) / np.pi * ft.G).astype(int) % ft.G).astype(np.uint8))
    return np.stack(out)

def main():
    class Tee:
        def __init__(self, path): self.f = open(path, "a"); self.o = sys.__stdout__
        def write(self, s): self.o.write(s); self.f.write(s)
        def flush(self): self.o.flush(); self.f.flush()
    sys.stdout = Tee("pipeline_fuzz_console.txt")
    def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    def log(s): print(f"[{now()}] {s}", flush=True)
    cores = A.cores or int(os.environ.get("CORES", 0)) or os.cpu_count() or 8
    NV = A.nsh * A.nsh; ntau, nsc = len(ft.TAUS), len(ft.SCALES); G, DL, _deg = ft.G, ft.DL, ft._deg
    angs = dsr.angs; settings = dsr.settings; nang = len(angs); NW = A.nw
    assert G <= 256, "bins must fit in uint8"
    partf = f"part_nsh{A.nsh}_{A.mode}{A.seed}_nang{A.nang}_nw{A.nw}.npz"
    histf = f"hist_nsh{A.nsh}_{A.mode}{A.seed}_nang{A.nang}_nw{A.nw}.npz"
    log(f"pipeline_fuzz_v1: nsh={A.nsh} ({NV} offsets) mode={A.mode} seed={A.seed} nang={A.nang} nw={A.nw} eps={A.eps} cores={cores} register={A.register}")
    log(f"shared-compute job {os.environ.get('SHARED_COMPUTE_JOB_ID','(none)')}, reserved {os.environ.get('SHARED_COMPUTE_MEM_GIB','?')} GiB")
    offsets = dsr.make_offsets(A.nsh, A.mode, A.seed)

    # ---------- 1. the particle store ----------
    if os.path.exists(partf):
        log(f"particle store exists: {partf} (skipping build)")
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
        Bpair = np.zeros((ntau, nsc, 4, NV, NW), np.uint8); Bang = np.zeros((ntau, nsc, nang, NV, NW), np.uint8)
        cnt_p = np.zeros((nsc, 4), int); cnt_a = np.zeros((nsc, nang), int)
        log(f"evolving {len(jobs)} trajectories on {cores} cores ({ntau} tau snapshots each); storing every particle's final bin "
            f"({(Bpair.nbytes + Bang.nbytes)/1e9:.2f} GB in memory as uint8)")
        t0 = time.time(); nxt = 0
        with ctx.Pool(cores) as pool:
            for k, o in enumerate(pool.imap(evolve_bins, jobs, chunksize=A.chunk)):
                kind, si, idx = meta[k]
                if kind == "pair": c = cnt_p[si, idx]; cnt_p[si, idx] += 1; Bpair[:, si, idx, c] = o
                else: c = cnt_a[si, idx]; cnt_a[si, idx] += 1; Bang[:, si, idx, c] = o
                if k + 1 >= nxt:
                    el = time.time() - t0; rate = (k + 1) / max(el, 1e-9)
                    log(f"  {k+1}/{len(jobs)} trajectories  ({el/3600:.2f} h elapsed, ~{(len(jobs)-k-1)/rate/3600:.2f} h left)"); nxt += len(jobs) // 40
        log(f"build done in {(time.time()-t0)/3600:.2f} h")
        assert cnt_p.min() == NV and cnt_a.min() == NV
        Bpair = Bpair.reshape(ntau * nsc, 4, NV, NW); Bang = Bang.reshape(ntau * nsc, nang, NV, NW)
        if os.path.exists(histf):          # the bins' histograms must be the stored histograms, bit for bit
            z = np.load(histf); bad = 0
            for ti in range(ntau * nsc):
                for p in range(4):
                    H = np.stack([np.bincount(Bpair[ti, p, v], minlength=G) for v in range(NV)]); bad += int(np.any(H != z["Hpair"][ti, p]))
                for ai in range(0, nang, max(1, nang // 8)):
                    H = np.stack([np.bincount(Bang[ti, ai, v], minlength=G) for v in range(NV)]); bad += int(np.any(H != z["Hang"][ti, ai]))
            log(f"check against {histf}: {bad} (cell, setting) blocks differ (must be 0; all 4 CHSH settings and every ~8th angle compared)")
            if bad: sys.exit("ABORT: particle bins do not reproduce the stored histograms")
        else: log(f"(no {histf} in the folder: histogram check skipped)")
        t0 = time.time()
        np.savez_compressed(partf, Bpair=Bpair, Bang=Bang, offsets=np.array(offsets, float), angs=np.array(angs, float), taus=np.array(ft.TAUS), scales=np.array(ft.SCALES),
                            nw=np.array(NW), G=np.array(G), ordering=np.array("tau_major"), register=np.array(A.register))
        log(f"particle store saved to {partf} ({os.path.getsize(partf)/1e9:.2f} GB on disk, {(time.time()-t0)/60:.1f} min)")
        del Bpair, Bang
    if A.build_only: log("build only -- done"); return

    # ---------- 2. contraction (once) ----------
    cf = f"fuzz_contracted_Lmax{A.Lmax}_{os.path.splitext(partf)[0]}.npz"
    if not os.path.exists(cf):
        log(f"contracting with Lmax={A.Lmax} -> {cf}")
        subprocess.run([sys.executable, "-c", f"import census_fuzz_v1 as c; c.contract('{partf}', {A.Lmax}, '{cf}')"], check=True)
    else: log(f"contraction exists: {cf}")

    # ---------- 3. the LPs ----------
    LL = [int(x) for x in A.Llist.split(",")]; KL = [int(x) for x in A.klist.split(",")]; SC = A.scopes.split(",")
    grid = [(L, K, s) for L in LL for K in KL for s in (SC if L > 0 else ["cell"])]      # L = 0 has no fuzz: one scope suffices
    pending = []; running = {}
    for L, K, s in grid:
        logf = f"fuzz_L{L}_K{K}_{s}.txt"
        if os.path.exists(logf) and "SUMMARY" in open(logf).read(): log(f"L={L} K={K} {s}: complete (skipping)"); continue
        pending.append((L, K, s, logf))
    t0 = time.time()
    while pending or running:
        while pending and len(running) < A.maxpar:
            L, K, s, logf = pending.pop(0); fh = open(logf, "w")
            running[logf] = (subprocess.Popen([sys.executable, "census_fuzz_v1.py", partf, "--L", str(L), "--K", str(K), "--scope", s, "--eps", str(A.eps), "--Lmax", str(A.Lmax)], stdout=fh, stderr=subprocess.STDOUT), fh)
            log(f"started L={L} K={K} {s} -> {logf}")
        time.sleep(30)
        for logf in list(running):
            pr, fh = running[logf]
            if pr.poll() is not None: fh.close(); del running[logf]; log(f"{logf} finished exit={pr.returncode} after {(time.time()-t0)/60:.0f} min")

    # ---------- 4. summary ----------
    print(f"\nSUMMARY  the disc's edge: MD/floor at eps={A.eps}, step reader, lattice; rows = birth fuzz (L harmonics in lambda_0, scope), columns = census bandwidth K")
    print(f"{'fuzz':<14} " + " ".join(f"{'K='+str(K):>9}" for K in KL))
    ref = {6: 1.700, 8: 1.429, 10: 1.258, 12: 1.112, 14: 1.000}
    print(f"{'paper (ref)':<14} " + " ".join(f"{ref.get(K, float('nan')):9.3f}" for K in KL))
    for L in LL:
        for s in (SC if L > 0 else ["cell"]):
            cells = []
            for K in KL:
                jf = f"census_fuzz_v1_L{L}_K{K}_{s}.json"
                if not os.path.exists(jf): cells.append("?"); continue
                d = json.load(open(jf)); cells.append(f"{d['price']:.4f}" if "price" in d else (f"m{d['min_miss']:.4f}" if "min_miss" in d else "-"))
            print(f"{('L='+str(L)+' '+s) if L > 0 else 'L=0 (uniform)':<14} " + " ".join(f"{c:>9}" for c in cells))
    print("reading: L=0 must reproduce the paper's K-curve.  If a shaped rim (L >= 1) lowers the K=6 price toward 1.000, the grain was partly the disc's edge;")
    print("         if the prices stay put, the birth fuzz is not where the last percent hides.  'cell' = rim depends on the dial only; 'atom' = on the approach phases too.")
    log("pipeline_fuzz_v1 finished")

if __name__ == "__main__":
    main()
