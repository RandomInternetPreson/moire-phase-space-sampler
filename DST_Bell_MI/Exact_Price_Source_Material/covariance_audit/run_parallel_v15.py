#!/usr/bin/env python3
"""
run_parallel_v15.py -- batch 15: the sorting rule on the other two samplings (census_bumps_v2.py, K=6 + 6-degree bumps, 16x16,
three perturbed vertices each) -- random0 (1024 random offsets) and the 4096-offset census.  Six light jobs; the 4096 ones ~1 h.
Then:
  python beads_control_v3.py census_bumps_v2_r0ctl_p1_rho_nsh32_random0_eps0.01_nang161.json census_bumps_v2_r0ctl_p2_rho_nsh32_random0_eps0.01_nang161.json census_bumps_v2_r0ctl_p3_rho_nsh32_random0_eps0.01_nang161.json > beads_control_v3_random0_10_29_26.txt
  python beads_control_v3.py census_bumps_v2_n64ctl_p1_rho_nsh64_random0_eps0.01_nang161.json census_bumps_v2_n64ctl_p2_rho_nsh64_random0_eps0.01_nang161.json census_bumps_v2_n64ctl_p3_rho_nsh64_random0_eps0.01_nang161.json > beads_control_v3_nsh64_10_29_26.txt
Needs in covariance_audit: census_bumps_v2.py, census_bandlimit_v4.py, both rho files + their banks.
Usage (from covariance_audit):   python run_parallel_v15.py        (or RUN_BATCH15.bat)
"""
import os, sys, time, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
R0 = "rho_nsh32_random0_eps0.01_nang161.json"
R0 = "rho_nsh32_random0_eps0.01_nang161.json"; N64 = "rho_nsh64_random0_eps0.01_nang161.json"; CB = "census_bumps_v2.py"
def ctl(prefix, rho, seed, pert):
    tag = f"{prefix}_p{seed}"
    return (tag, CB, rho, dict(K="6", BUMP_W="6", BUMP_NC="16", EPS="0.01", METHOD="highs-ipm", PERTURB=pert, PSEED=str(seed), OUT_TAG=tag), f"bumps_{tag}.txt")
JOBS = [ctl("r0ctl", R0, 1, "1e-7"), ctl("r0ctl", R0, 2, "1e-7"), ctl("r0ctl", R0, 3, "1e-5"),
        ctl("n64ctl", N64, 1, "1e-7"), ctl("n64ctl", N64, 2, "1e-7"), ctl("n64ctl", N64, 3, "1e-5")]
if len(sys.argv) > 1 and sys.argv[1] == "--tiny":
    T = "rho_nsh8_random0_eps0.3_nang13.json"
    JOBS = [(t, s, T, dict(e, EPS="0.3", K="2", BUMP_NC="8", BUMP_W="10"), log) for (t, s, _, e, log) in JOBS[:1]]
STATUS, CONSOLE = "parallel15_status.txt", "parallel15_console.txt"

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
        f.write(f"parallel run 15, elapsed {(time.time()-t0)/3600:.2f} h, {len(done)}/{len(JOBS)} finished\n")
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
