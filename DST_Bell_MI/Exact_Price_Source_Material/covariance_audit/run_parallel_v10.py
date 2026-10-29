#!/usr/bin/env python3
"""
run_parallel_v10.py -- batch 10: birth_model_v3.py on the LATTICE banks (overnight; 9 jobs x 3 h, light).

The model: per family, an axis with its 90-degree harmonic, a push with cos 2D..cos 8D, axis-coupled; every parameter linear
in the terrain scale (and optionally log tau); free cell weights; station factor.  Warm-started from the lattice K=6 census.
Jobs (rho_nsh32_grid0; MODE miss unless noted; eps 0.01; TIME_LIMIT 10800 s):
  a_j2_s0t0        J=2, NO dial dependence   -- control: the new harmonics alone (v1/v2 form reached 2.4%)
  b_j2_s1t0        J=2, scale dependence
  c_j2_s1t1        J=2, scale + tau dependence
  c2_j2_s1t1_s2    same, seed 2
  d_j1_s1t1        J=1, scale + tau
  e_j3_s1t1        J=3, scale + tau
  g_j2_s1t1_price  J=2, scale + tau, PRICE at eps=0.01
  h_j2_s1t1_p02    J=2, scale + tau, PRICE at eps=0.02
  i_j3_s1t1_price  J=3, scale + tau, PRICE at eps=0.01
Targets on the lattice: K=6 census 1.700 (feasible at 1%); free dial-blind 1.473 (miss 0.38%); K=5 census misses 1.055%.
Reading: a vs b vs c separates "harmonics" from "depth" from "duration"; the price jobs say what a model birth pays.
Needs in covariance_audit: birth_model_v3.py, rho_nsh32_grid0_eps0.01_nang161.json + banks, census_bandlimit_v4_n32_k6_rho_nsh32_grid0_eps0.01_nang161.json.
Usage (from covariance_audit):   python run_parallel_v10.py        (or RUN_BATCH10.bat)
Monitoring: bm3_<tag>.txt per job, parallel10_status.txt, parallel10_console.txt.
"""
import os, sys, time, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
R0 = "rho_nsh32_random0_eps0.01_nang161.json"
G = "rho_nsh32_grid0_eps0.01_nang161.json"; BM3 = "birth_model_v3.py"; K6 = "census_bandlimit_v4_n32_k6_rho_nsh32_grid0_eps0.01_nang161.json"
def bm(tag, J, sdep, tdep, mode, eps, seed):
    return (tag, BM3, G, dict(FAMILIES=str(J), SCALE_DEP=sdep, TAU_DEP=tdep, STATION="1", MODE=mode, EPS=eps, SEED=str(seed), TIME_LIMIT="10800",
                              NSTART="12", POP="10", K6_JSON=K6, K6_K="6", OUT_TAG=tag), f"bm3_{tag}.txt")
JOBS = [
    bm("a_j2_s0t0", 2, "0", "0", "miss", "0.01", 1),
    bm("b_j2_s1t0", 2, "1", "0", "miss", "0.01", 1),
    bm("c_j2_s1t1", 2, "1", "1", "miss", "0.01", 1),
    bm("c2_j2_s1t1_s2", 2, "1", "1", "miss", "0.01", 2),
    bm("d_j1_s1t1", 1, "1", "1", "miss", "0.01", 1),
    bm("e_j3_s1t1", 3, "1", "1", "miss", "0.01", 1),
    bm("g_j2_s1t1_price", 2, "1", "1", "price", "0.01", 1),
    bm("h_j2_s1t1_p02", 2, "1", "1", "price", "0.02", 1),
    bm("i_j3_s1t1_price", 3, "1", "1", "price", "0.01", 1),
]
if len(sys.argv) > 1 and sys.argv[1] == "--tiny":
    T = "rho_nsh8_random0_eps0.3_nang13.json"
    JOBS = [(t, s, T, dict(e, EPS="0.3", TIME_LIMIT="20", NSTART="2", K6_JSON="census_bandlimit_v4_t4_rho_nsh8_random0_eps0.3_nang13.json", K6_K="3"), log) for (t, s, _, e, log) in JOBS[:2]]
STATUS, CONSOLE = "parallel10_status.txt", "parallel10_console.txt"

def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
missing = [f for f in {j[1] for j in JOBS} | {j[2] for j in JOBS} | {j[3]["K6_JSON"] for j in JOBS} if not os.path.exists(f)]
if missing: sys.exit(f"missing: {missing}")
procs = {}; t0 = time.time(); done = {}
for tag, script, rho, envadd, log in JOBS:
    env = dict(os.environ); env.update(envadd)
    f = open(log, "w")
    procs[tag] = (subprocess.Popen([sys.executable, script, rho], env=env, stdout=f, stderr=subprocess.STDOUT), f, script, log)
    print(f"[{now()}] started {tag}: {script} {rho} {envadd}", flush=True)
def write_status():
    with open(STATUS, "w") as f:
        f.write(f"parallel run 10, elapsed {(time.time()-t0)/3600:.2f} h, {len(done)}/{len(JOBS)} finished\n")
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
