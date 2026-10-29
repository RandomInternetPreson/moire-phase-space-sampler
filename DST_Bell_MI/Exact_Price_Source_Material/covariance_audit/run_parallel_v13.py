#!/usr/bin/env python3
"""
run_parallel_v13.py -- batch 13: the CREST test (census_bumps_v1.py on the lattice census): smooth K=6 + a dictionary of
localized bumps; how many crests does the floor need, and do their positions travel with tau?

Jobs (rho_nsh32_grid0, eps=0.01):
  cr_k6_w6_nc16     K=6 + 256 bumps/cell of half-width 6 deg      (~30 min)
  cr_k6_w12_nc16    K=6 + 256 bumps/cell, 12 deg
  cr_k6_w9_nc16     K=6 + 256 bumps/cell, 9 deg
  cr_k6_w6_nc32     K=6 + 1024 bumps/cell, 6 deg (finer placement; ~K=16-sized LP, hours)
  cr_k6_w3_nc32     K=6 + 1024 bumps/cell, 3 deg
  cr_k0_w6_nc16     bumps only (no smooth part), 6 deg -- control
Reading: price near 1.0 with ~1-3 active bumps per cell whose mean position advances with tau -> a traveling crest;
dozens of bumps everywhere -> the fine structure is not a single localized event.  Active counts are upper bounds (LP vertices).
Needs in covariance_audit: census_bumps_v1.py, census_bandlimit_v4.py, rho_nsh32_grid0_eps0.01_nang161.json + banks.
Usage (from covariance_audit):   python run_parallel_v13.py        (or RUN_BATCH13.bat)
Monitoring: bumps_<tag>.txt per job, parallel13_status.txt, parallel13_console.txt.
"""
import os, sys, time, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
R0 = "rho_nsh32_random0_eps0.01_nang161.json"
G = "rho_nsh32_grid0_eps0.01_nang161.json"; CB = "census_bumps_v1.py"
def cr(K, w, nc):
    tag = f"cr_k{K}_w{w}_nc{nc}"
    return (tag, CB, G, dict(K=str(K), BUMP_W=str(w), BUMP_NC=str(nc), EPS="0.01", METHOD="highs-ipm", OUT_TAG=tag), f"bumps_{tag}.txt")
JOBS = [cr(6, 6, 16), cr(6, 12, 16), cr(6, 9, 16), cr(6, 6, 32), cr(6, 3, 32), cr(0, 6, 16)]
if len(sys.argv) > 1 and sys.argv[1] == "--tiny":
    T = "rho_nsh8_random0_eps0.3_nang13.json"
    JOBS = [(t, s, T, dict(e, EPS="0.3", K="2", BUMP_NC="8", BUMP_W="10"), log) for (t, s, _, e, log) in JOBS[:1]]
STATUS, CONSOLE = "parallel13_status.txt", "parallel13_console.txt"

def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
missing = [f for f in {j[1] for j in JOBS} | {j[2] for j in JOBS} | {"census_bandlimit_v4.py"} if not os.path.exists(f)]
if missing: sys.exit(f"missing: {missing}")
procs = {}; t0 = time.time(); done = {}
for tag, script, rho, envadd, log in JOBS:
    env = dict(os.environ); env.update(envadd)
    f = open(log, "w")
    procs[tag] = (subprocess.Popen([sys.executable, script, rho], env=env, stdout=f, stderr=subprocess.STDOUT), f, script, log)
    print(f"[{now()}] started {tag}: {script} {rho} {envadd}", flush=True)
def write_status():
    with open(STATUS, "w") as f:
        f.write(f"parallel run 13, elapsed {(time.time()-t0)/3600:.2f} h, {len(done)}/{len(JOBS)} finished\n")
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
