#!/usr/bin/env python3
"""
run_parallel_v9.py -- batch 9: the (eps, K) surface on the LATTICE census (census_bandlimit_v4.py).

Why.  v3's (eps, K) surface was computed on 1024 random offsets and carries the quadrature cost; v6 of the paper dropped it.
This is the clean surface: K = 5, 6, 8, 10 at eps = 0.02, 0.05, 0.10 on the 32x32 lattice (eps = 0.01 column already in
bl_n32_k*.txt: K=5 infeasible (1.055%), 6: 1.700, 8: 1.429, 10: 1.258).  12 single-threaded processes, < 2 GB each, ~1 h.
Needs in covariance_audit: census_bandlimit_v4.py, rho_nsh32_grid0_eps0.01_nang161.json, banks_nsh32_grid0_nang161_nw2048.npz.

Usage (from covariance_audit):   python run_parallel_v9.py        (or RUN_BATCH9.bat)
Monitoring: bl_g_e<eps>_k<K>.txt per job, parallel9_status.txt, parallel9_console.txt.
"""
import os, sys, time, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
R0 = "rho_nsh32_random0_eps0.01_nang161.json"
G = "rho_nsh32_grid0_eps0.01_nang161.json"; BL = "census_bandlimit_v4.py"
def bl(K, eps):
    tag = f"g_e{eps.replace('0.', '')}_k{K}"
    return (tag, BL, G, dict(K_LIST=str(K), EPS_LIST=eps, OUT_TAG=tag, TIME_LIMIT="0", METHOD="highs-ipm", METHOD2="highs-ds", MISS_ONLY="0"), f"bl_{tag}.txt")
JOBS = [bl(K, e) for K in (5, 6, 8, 10) for e in ("0.02", "0.05", "0.10")]       # 12 jobs; K=5 is feasible at these eps (min miss 1.055%)
if len(sys.argv) > 1 and sys.argv[1] == "--tiny":
    T = "rho_nsh8_random0_eps0.3_nang13.json"
    JOBS = [(t, s, T, dict(e, EPS_LIST="0.3"), log) for (t, s, _, e, log) in JOBS[:2]]
STATUS, CONSOLE = "parallel9_status.txt", "parallel9_console.txt"

def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
missing = [f for f in {j[1] for j in JOBS} | {j[2] for j in JOBS} if not os.path.exists(f)]
if missing: sys.exit(f"missing: {missing}")
procs = {}; t0 = time.time(); done = {}
for tag, script, rho, envadd, log in JOBS:
    env = dict(os.environ); env.update(envadd)
    f = open(log, "w")
    procs[tag] = (subprocess.Popen([sys.executable, script, rho], env=env, stdout=f, stderr=subprocess.STDOUT), f, script, log)
    print(f"[{now()}] started {tag}: {script} {rho} {envadd}", flush=True)
def write_status():
    with open(STATUS, "w") as f:
        f.write(f"parallel run 9, elapsed {(time.time()-t0)/3600:.2f} h, {len(done)}/{len(JOBS)} finished\n")
        for tag, (p, fh, script, log) in procs.items():
            f.write(f"  {tag:18s} {script:24s}: " + (f"finished exit={done[tag][0]} after {done[tag][1]/60:.0f} min" if tag in done else "running") + f"   (log {log})\n")
write_status()
while len(done) < len(JOBS):
    time.sleep(20)
    for tag, (p, fh, script, log) in procs.items():
        if tag not in done and p.poll() is not None:
            done[tag] = (p.returncode, time.time() - t0); fh.close()
            print(f"[{now()}] finished {tag} exit={p.returncode} after {done[tag][1]/60:.0f} min", flush=True); write_status()
print(f"\nall {len(JOBS)} jobs finished in {(time.time()-t0)/3600:.2f} h\n")
print("SUMMARY lines from each log:")
for tag, (p, fh, script, log) in procs.items():
    lines = open(log).read().splitlines()
    i = next((k for k, l in enumerate(lines) if l.startswith("SUMMARY")), None)
    print(f"--- {tag}"); print("\n".join(lines[i:i+6]) if i is not None else "  (no SUMMARY -- see log)")
