#!/usr/bin/env python3
"""
run_parallel_v2.py -- launch independent LP jobs as parallel processes (one per core-group), wait, summarize.

Each job is one census_bandlimit_v3.py run with its own environment and its own output file. HiGHS is
essentially single-threaded, so N independent jobs on N cores finish in the time of the longest one.
Edit the JOBS list to add/remove runs. Memory: a K=10-12 job needs a few GB; 8 jobs are fine on 256 GB.

Usage (from covariance_audit, where the rho files and their banks are):   python run_parallel_v1.py
Monitoring: each job writes bl_<tag>.txt live; parallel_status.txt is rewritten when any job finishes.
"""
import os, sys, time, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
SCRIPT = "census_bandlimit_v3.py"
R0, R1, S0, S1 = ("rho_nsh32_random0_eps0.01_nang161.json", "rho_nsh32_random1_eps0.01_nang161.json",
                  "rho_nsh32_strat0_eps0.01_nang161.json", "rho_nsh32_strat1_eps0.01_nang161.json")
# tag, rho file, K_LIST, EPS_LIST
JOBS = [
    ("s0_k5",      S0, "5",  "0.01"),   # certificate on the two stratified censuses (3 min each)
    ("s1_k5",      S1, "5",  "0.01"),
    ("r0_e02_k8",  R0, "8",  "0.02"),   # fill the (eps, K) surface: eps = 0.02 row
    ("r0_e02_k10", R0, "10", "0.02"),
    ("r0_e10_k6",  R0, "6",  "0.10"),   # eps = 0.10 row
    ("r0_e10_k8",  R0, "8",  "0.10"),
    ("r0_e10_k10", R0, "10", "0.10"),
]
if len(sys.argv) > 1 and sys.argv[1] == "--tiny":        # self-test on toy banks
    JOBS = [("t1", "rho_nsh8_random0_eps0.3_nang13.json", "2,3", "0.3"), ("t2", "rho_nsh8_random0_eps0.3_nang13.json", "4", "0.3")]

def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
missing = [f for f in {SCRIPT} | {j[1] for j in JOBS} if not os.path.exists(f)]
if missing: sys.exit(f"missing: {missing}")
procs = {}; t0 = time.time(); done = {}
for tag, rho, ks, eps in JOBS:
    env = dict(os.environ); env.update(K_LIST=ks, EPS_LIST=eps, OUT_TAG=tag, TIME_LIMIT="0", METHOD="highs-ipm")
    log = open(f"bl_{tag}.txt", "w")
    procs[tag] = (subprocess.Popen([sys.executable, SCRIPT, rho], env=env, stdout=log, stderr=subprocess.STDOUT), log, rho, ks, eps)
    print(f"[{now()}] started {tag}: {rho} K={ks} eps={eps}", flush=True)
def write_status():
    with open("parallel2_status.txt", "w") as f:
        f.write(f"parallel run, elapsed {(time.time()-t0)/3600:.2f} h, {len(done)}/{len(JOBS)} finished\n")
        for tag, (p, log, rho, ks, eps) in procs.items():
            f.write(f"  {tag:12s} {rho} K={ks} eps={eps}: " + (f"finished exit={done[tag][0]} after {done[tag][1]/60:.0f} min" if tag in done else "running") + f"   (log bl_{tag}.txt)\n")
write_status()
while len(done) < len(JOBS):
    time.sleep(20)
    for tag, (p, log, rho, ks, eps) in procs.items():
        if tag not in done and p.poll() is not None:
            done[tag] = (p.returncode, time.time() - t0); log.close()
            print(f"[{now()}] finished {tag} exit={p.returncode} after {done[tag][1]/60:.0f} min", flush=True); write_status()
print(f"\nall {len(JOBS)} jobs finished in {(time.time()-t0)/3600:.2f} h\n")
print("SUMMARY lines from each log:")
for tag in procs:
    lines = open(f"bl_{tag}.txt").read().splitlines()
    i = next((k for k, l in enumerate(lines) if l.startswith("SUMMARY")), None)
    print(f"--- {tag}"); print("\n".join(lines[i:i+8]) if i is not None else "  (no SUMMARY -- see log)")
