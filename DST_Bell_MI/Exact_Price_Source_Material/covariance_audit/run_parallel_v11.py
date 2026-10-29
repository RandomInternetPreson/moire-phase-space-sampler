#!/usr/bin/env python3
"""
run_parallel_v11.py -- batch 11: the PLEAT scan (census_pleat_v1.py on the lattice census).

Question: does one high harmonic -- a periodic ripple with N pleats around the disc's hem -- buy the floor for a smooth
(K=6) census?  Price vs N is a measurement of the pleat count on this terrain.
Jobs (rho_nsh32_grid0, eps=0.01, K=6):
  pl_N<N>_lines   N = 7..20: the pleat pixels (N,0),(0,N),(N,N),(N,-N) added to the K=6 basis   -- 14 jobs, ~10 min each
  pl_N<N>_blocks  N = 8, 12, 16: every sideband of pleat x smooth (K=6 block shifted by the four pleat vectors) -- 3 jobs, ~4 h each
Reference: K=6 alone 1.700; K=8 1.429; K=10 1.258; K=12 1.112; K=14 1.000 (lattice).
Reading: a sharp minimum in price vs N (lines) names the pleat count; a flat curve says the fine structure is not one ripple.
Needs in covariance_audit: census_pleat_v1.py, rho_nsh32_grid0_eps0.01_nang161.json + banks_nsh32_grid0_nang161_nw2048.npz.
Usage (from covariance_audit):   python run_parallel_v11.py        (or RUN_BATCH11.bat)
Monitoring: bl_pl_<tag>.txt per job, parallel11_status.txt, parallel11_console.txt.
"""
import os, sys, time, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
R0 = "rho_nsh32_random0_eps0.01_nang161.json"
G = "rho_nsh32_grid0_eps0.01_nang161.json"; PL = "census_pleat_v1.py"
def pl(N, mode):
    tag = f"pl_N{N}_{mode}"
    return (tag, PL, G, dict(K_LIST="6", EPS_LIST="0.01", PLEAT_N=str(N), PLEAT_MODE=mode, OUT_TAG=tag, TIME_LIMIT="0", METHOD="highs-ipm", METHOD2="highs-ds", MISS_ONLY="0"), f"bl_{tag}.txt")
JOBS = [pl(N, "lines") for N in range(7, 21)] + [pl(N, "blocks") for N in (8, 12, 16)]
if len(sys.argv) > 1 and sys.argv[1] == "--tiny":
    T = "rho_nsh8_random0_eps0.3_nang13.json"
    JOBS = [(t, s, T, dict(e, EPS_LIST="0.3", K_LIST="2"), log) for (t, s, _, e, log) in JOBS[:2]]
STATUS, CONSOLE = "parallel11_status.txt", "parallel11_console.txt"

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
        f.write(f"parallel run 11, elapsed {(time.time()-t0)/3600:.2f} h, {len(done)}/{len(JOBS)} finished\n")
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
