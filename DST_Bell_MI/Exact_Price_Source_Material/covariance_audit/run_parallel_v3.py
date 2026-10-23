#!/usr/bin/env python3
"""
run_parallel_v3.py -- batch 3: repair the eps = 0.10 row and certify every "infeasible" by value (census_bandlimit_v4.py).

Why.  Batch 2 returned "infeasible" for eps = 0.10 at K = 6 and 8 while eps = 0.01 is feasible at both -- impossible,
since eps only widens the correlator band.  v1-v3 read every non-optimal solver status as "infeasible"; v4 prints the
status, certifies infeasibility by the minimum achievable miss (a second LP, re-measured from the census), and retries
a numerically failed price LP with dual simplex.

Jobs (all on the 1024-offset banks, OUT_TAG = job tag):
  r0_e10_k6, r0_e10_k8      price at eps = 0.10 (expected feasible; bracketed by 1.1317 at K = 10 and the eps = 0.05 row)
  r0_e01_k6_v4              re-run of the pass-2 K = 6 price (2.274) through v4's verifier -- a cross-check, not a new number
  r0_k2 .. r0_k5            min miss of the best K-band census on random0 (MISS_ONLY): the miss curve eps_min(K) below K = 6;
                            the miss is independent of eps, so one run gives the whole column of the (eps, K) surface
  r1_k5, s0_k5, s1_k5       the K = 5 certificate by value on the other three censuses
Memory: each job holds the banks (~0.5 GB) plus the LP; 10 jobs are fine on 256 GB -- but do NOT run alongside the
swap job (56 workers + 130 GB), which is what exhausted the page file on 22 Oct.  Start this after the swap finishes.

Usage (from covariance_audit):   python run_parallel_v3.py        (or RUN_PARALLEL3.bat)
Monitoring: bl_<tag>.txt per job (live), parallel3_status.txt rewritten when any job finishes, parallel3_console.txt.
"""
import os, sys, time, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
SCRIPT = "census_bandlimit_v4.py"
R0, R1, S0, S1 = ("rho_nsh32_random0_eps0.01_nang161.json", "rho_nsh32_random1_eps0.01_nang161.json",
                  "rho_nsh32_strat0_eps0.01_nang161.json", "rho_nsh32_strat1_eps0.01_nang161.json")
# tag, rho file, K_LIST, EPS_LIST, MISS_ONLY
JOBS = [
    ("r0_e10_k6",    R0, "6",  "0.10", "0"),
    ("r0_e10_k8",    R0, "8",  "0.10", "0"),
    ("r0_e01_k6_v4", R0, "6",  "0.01", "0"),
    ("r0_k2",        R0, "2",  "0.01", "1"),
    ("r0_k3",        R0, "3",  "0.01", "1"),
    ("r0_k4",        R0, "4",  "0.01", "1"),
    ("r0_k5",        R0, "5",  "0.01", "1"),
    ("r1_k5",        R1, "5",  "0.01", "1"),
    ("s0_k5",        S0, "5",  "0.01", "1"),
    ("s1_k5",        S1, "5",  "0.01", "1"),
]
if len(sys.argv) > 1 and sys.argv[1] == "--tiny":        # self-test on toy banks
    JOBS = [("t1", "rho_nsh8_random0_eps0.3_nang13.json", "2,3", "0.3", "0"), ("t2", "rho_nsh8_random0_eps0.3_nang13.json", "1,2", "0.001", "1")]
STATUS, CONSOLE = "parallel3_status.txt", "parallel3_console.txt"

def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
missing = [f for f in {SCRIPT} | {j[1] for j in JOBS} if not os.path.exists(f)]
if missing: sys.exit(f"missing: {missing}")
procs = {}; t0 = time.time(); done = {}
for tag, rho, ks, eps, mo in JOBS:
    env = dict(os.environ); env.update(K_LIST=ks, EPS_LIST=eps, OUT_TAG=tag, TIME_LIMIT="0", METHOD="highs-ipm", METHOD2="highs-ds", MISS_ONLY=mo)
    log = open(f"bl_{tag}.txt", "w")
    procs[tag] = (subprocess.Popen([sys.executable, SCRIPT, rho], env=env, stdout=log, stderr=subprocess.STDOUT), log, rho, ks, eps, mo)
    print(f"[{now()}] started {tag}: {rho} K={ks} eps={eps} miss_only={mo}", flush=True)
def write_status():
    with open(STATUS, "w") as f:
        f.write(f"parallel run 3, elapsed {(time.time()-t0)/3600:.2f} h, {len(done)}/{len(JOBS)} finished\n")
        for tag, (p, log, rho, ks, eps, mo) in procs.items():
            f.write(f"  {tag:14s} {rho} K={ks} eps={eps}{' miss-only' if mo == '1' else ''}: " + (f"finished exit={done[tag][0]} after {done[tag][1]/60:.0f} min" if tag in done else "running") + f"   (log bl_{tag}.txt)\n")
write_status()
while len(done) < len(JOBS):
    time.sleep(20)
    for tag, (p, log, rho, ks, eps, mo) in procs.items():
        if tag not in done and p.poll() is not None:
            done[tag] = (p.returncode, time.time() - t0); log.close()
            print(f"[{now()}] finished {tag} exit={p.returncode} after {done[tag][1]/60:.0f} min", flush=True); write_status()
print(f"\nall {len(JOBS)} jobs finished in {(time.time()-t0)/3600:.2f} h\n")
print("SUMMARY lines from each log:")
for tag in procs:
    lines = open(f"bl_{tag}.txt").read().splitlines()
    i = next((k for k, l in enumerate(lines) if l.startswith("SUMMARY")), None)
    print(f"--- {tag}"); print("\n".join(lines[i:i+8]) if i is not None else "  (no SUMMARY -- see log)")
