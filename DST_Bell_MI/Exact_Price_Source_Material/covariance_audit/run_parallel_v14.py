#!/usr/bin/env python3
"""
run_parallel_v14.py -- batch 14: the partition control -- census_bumps_v2.py (K=6 + 6-degree bumps, 16x16) with six different
objective perturbations (seeds 1-3 at 1e-7, seeds 4-6 at 1e-5), so six different optimal vertices of the same class.  ~20 min each, light.
Then:  python beads_control_v1.py census_bumps_v2_ctl_p*_rho_nsh32_grid0_eps0.01_nang161.json > beads_control_10_29_26.txt
Needs in covariance_audit: census_bumps_v2.py, census_bandlimit_v4.py, rho_nsh32_grid0_eps0.01_nang161.json + banks.
Usage (from covariance_audit):   python run_parallel_v14.py        (or RUN_BATCH14.bat)
Monitoring: bumps_ctl_p<seed>.txt per job, parallel14_status.txt, parallel14_console.txt.
"""
import os, sys, time, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
R0 = "rho_nsh32_random0_eps0.01_nang161.json"
G = "rho_nsh32_grid0_eps0.01_nang161.json"; CB = "census_bumps_v2.py"
def ctl(seed, pert):
    tag = f"ctl_p{seed}"
    return (tag, CB, G, dict(K="6", BUMP_W="6", BUMP_NC="16", EPS="0.01", METHOD="highs-ipm", PERTURB=pert, PSEED=str(seed), OUT_TAG=tag), f"bumps_{tag}.txt")
JOBS = [ctl(s, "1e-7") for s in (1, 2, 3)] + [ctl(s, "1e-5") for s in (4, 5, 6)]     # two perturbation sizes: seeds 1-3 tiny, 4-6 larger
if len(sys.argv) > 1 and sys.argv[1] == "--tiny":
    T = "rho_nsh8_random0_eps0.3_nang13.json"
    JOBS = [(t, s, T, dict(e, EPS="0.3", K="2", BUMP_NC="8", BUMP_W="10"), log) for (t, s, _, e, log) in JOBS[:2]]
STATUS, CONSOLE = "parallel14_status.txt", "parallel14_console.txt"

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
        f.write(f"parallel run 14, elapsed {(time.time()-t0)/3600:.2f} h, {len(done)}/{len(JOBS)} finished\n")
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
