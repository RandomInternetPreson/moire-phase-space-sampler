#!/usr/bin/env python3
"""
run_bulge2_v1.py -- bulge sweep 2 (1 Nov): the mesa's second-order geometry, on the stored lattice histograms.
Runs census_profile_v3.py for each profile in PROFILES, MAXPAR at a time, into bulge2_<ptag>.txt, then prints one table.
No trajectories are evolved: needs hist_nsh32_grid0_nang161_nw2048.npz (job 403) in the folder.  Console: bulge2_console.txt.
Profiles (Aaron's 1 Nov list):  tilt (cap tilted against the disc) | awall (front and back walls of different width -- breaks the port
rule) | mawall (same, twin born upside down) | bwall (kinked wall, port rule kept) | shift (wall position follows the approach phase:
square vs circle -- reader reads xi, so the JOINT MD column is the one to read for it).
Usage (from covariance_audit, through the launcher -- RUN_BULGE2.bat):  python run_bulge2_v1.py
"""
import os, sys, json, time, subprocess, datetime
HIST = os.environ.get("HIST", "hist_nsh32_grid0_nang161_nw2048.npz")
KLIST = os.environ.get("KLIST", "6,8,10"); EPS = os.environ.get("EPS", "0.01"); MAXPAR = int(os.environ.get("MAXPAR", "12"))
PROFILES = os.environ.get("PROFILES", "tilt:0.05,tilt:0.1,awall:3.6:7.2,awall:3.6:14.3,awall:1.8:7.2,mawall:3.6:7.2,mawall:3.6:14.3,bwall:1.8:7.2,bwall:7.2:1.8,shift:5,shift:10,shift:20").split(",")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_profile_v3 import ptag

def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
def log(s): print(f"[{now()}] {s}", flush=True)

def main():
    class Tee:
        def __init__(self, path): self.f = open(path, "a"); self.o = sys.__stdout__
        def write(self, x): self.o.write(x); self.f.write(x)
        def flush(self): self.o.flush(); self.f.flush()
    sys.stdout = Tee("bulge2_console.txt")
    if not os.path.exists(HIST): sys.exit(f"histogram store not found: {HIST}")
    log(f"bulge sweep 2: {len(PROFILES)} profiles, K list {KLIST}, eps {EPS}, {MAXPAR} at a time, on {HIST}")
    pending = []; running = {}
    for p in PROFILES:
        tag = ptag(p); logf = f"bulge2_{tag}.txt"
        if os.path.exists(logf) and "SUMMARY" in open(logf).read(): log(f"{p}: {logf} complete (skipping)"); continue
        pending.append((p, tag, logf))
    t0 = time.time()
    while pending or running:
        while pending and len(running) < MAXPAR:
            p, tag, logf = pending.pop(0); fh = open(logf, "w")
            running[tag] = (subprocess.Popen([sys.executable, "census_profile_v3.py", HIST, p, "--klist", KLIST, "--eps", EPS, "--full", "1"], stdout=fh, stderr=subprocess.STDOUT), fh, logf, p)
            log(f"started {p} -> {logf}")
        time.sleep(30)
        for tag in list(running):
            pr, fh, logf, p = running[tag]
            if pr.poll() is not None: fh.close(); del running[tag]; log(f"{p} finished exit={pr.returncode} after {(time.time()-t0)/60:.0f} min")
    KL = [int(x) for x in KLIST.split(",") if x]; cols = ["full"] + [f"K={K}" for K in KL]
    print(f"\nSUMMARY  bulge sweep 2 on {HIST}: MD/floor at eps={EPS} (lambda-marginal MD, the paper's); 'm' = singlet unreachable, min miss shown")
    print(f"{'profile':<16} {'port':>6} {'S_max/2rt2':>10} {'miss':>8} " + " ".join(f"{c:>8}" for c in cols) + f" {'K'+str(KL[0])+'/full':>8} {'jointMD':>8}   <A>raw  <A>sym  dE_sym")
    for p in PROFILES:
        jf = f"census_profile_v3_{ptag(p)}.json"
        if not os.path.exists(jf): print(f"{p:<16} (no json)"); continue
        d = json.load(open(jf)); r = {x["label"]: x for x in d["results"]}
        def cell(l):
            x = r.get(l); return "?" if x is None else (f"{x['price']:.4f}" if "price" in x else (f"m{x['min_miss']:.4f}" if "min_miss" in x else "-"))
        prem = f"{r[f'K={KL[0]}']['price']/r['full']['price']:.3f}" if "price" in r.get("full", {}) and "price" in r.get(f"K={KL[0]}", {}) else "-"
        jm = f"{r['full']['md_joint_over_floor']:.3f}" if "md_joint_over_floor" in r.get("full", {}) else "-"
        mg = d.get("marginals", {}); ms = d.get("marginals_sym", {})
        port = "broken" if "PORT RULE BROKEN" in d.get("desc", "") else "kept"
        print(f"{p:<16} {port:>6} {d['smax']/(2*2**0.5):10.6f} {d.get('miss_relaxed', float('nan')):8.4f} " + " ".join(f"{cell(c):>8}" for c in cols)
              + f" {prem:>8} {jm:>8}   {mg.get('A_max', float('nan')):.4f}  {ms.get('A_max', float('nan')):.1e}  {ms.get('max_dE', float('nan')):.1e}")
    print("reference (sweep 1, step): full 1.0000, K=6 1.6996, K=8 1.4290, K=10 1.2579; tanh:8 (10.6 deg wall, both edges) 1.0476 / 2.5615 / 2.0004 / 1.7401.")
    print("reading: a port-broken reader has physical marginals (<A>sym != 0 and dE_sym != 0 mean the +90 symmetrization is no longer free);")
    print("         for 'shift' the reader reads xi, so its price must be read against jointMD, not the marginal MD.  The grain question is K6/full and the absolute K=6.")
    log("run_bulge2_v1 finished")

if __name__ == "__main__":
    main()
