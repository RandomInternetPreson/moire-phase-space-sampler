#!/usr/bin/env python3
"""
pipeline_v1.py -- self-driving terrain search for a station-local census (Oct 18, 2026).

Runs, in order, with the results of each stage deciding the next:

  0. SMOKE   : tiny-scale run of every script the pipeline will use (breed_product_v1.py,
               dither_sweep_rho.py, product_census_v3.py), in a sandbox folder, so a broken script
               fails in the first few minutes instead of hour four.  Any failure here STOPS the pipeline.
  A. BREED-G : breed_product_v1.py from the greedy terrain (START=greedy), ROUNDS_A rounds.
       -> if PROMISING (see rule)  : V-A. VALIDATE the champion at full scale, then FINISH.
       -> else                      : stage B.
  B. BREED-M : breed_product_v1.py from a one-term terrain (START=minimal), ROUNDS_B rounds.
       -> if PROMISING              : V-B. VALIDATE, then FINISH.
       -> else                      : FINISH with "greedy terrain search finds no station-local terrain".

  PROMISING rule (pre-registered): the breed's final product-class census is feasible at the breeding
  EPS and its MD is <= PROMISING_RATIO x the joint price on the same terrain (default 1.3).

  VALIDATE = copy the champion to a fresh folder as greedy_results.json, build full-scale banks
  (dither_sweep_rho.py: NANG=161, NSH=32 random, eps 0.01, SAVE_RHO), then product_census_v3.py
  (NB=32, eps 0.05 and 0.01, 12 starts).  The full-scale product row at eps=0.01 is the number that counts.

Monitoring while it runs:
  pipeline_status.txt   -- rewritten at every transition: what is running, since when, which log to open
  logs/<stage>.txt      -- the live output of the running stage (scripts flush every line)
  STOP.txt              -- create this file in the folder and the pipeline will not start another stage
  pipeline_FINAL.txt    -- written at the end: branch taken, numbers, verdict

Requirements in this folder: breed_product_v1.py, flight_time.py, dither_sweep_rho.py,
product_census_v3.py, greedy_results.json.  Environment: CORES (default 56).  Set PIPELINE_TINY=1 to
run the whole decision tree at toy scale in a few minutes (used to test the pipeline itself).
"""
import os, sys, json, time, shutil, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
CORES = os.environ.get("CORES", "56")
TINY = int(os.environ.get("PIPELINE_TINY", 0))
PROMISING_RATIO = float(os.environ.get("PROMISING_RATIO", 1.3))
ROUNDS_A = os.environ.get("ROUNDS_A", "8"); ROUNDS_B = os.environ.get("ROUNDS_B", "8")
FLOOR = 0.1380711874576983
REQUIRED = ["breed_product_v1.py", "flight_time.py", "dither_sweep_rho.py", "product_census_v3.py", "greedy_results.json"]
LOGS = os.path.join(HERE, "logs"); os.makedirs(LOGS, exist_ok=True)
STATUS = os.path.join(HERE, "pipeline_status.txt")
T0 = time.time()
TRANSITIONS = []

def now(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def status(stage, detail, log=None):
    line = f"[{now()}] {stage}: {detail}" + (f"   (log: {log})" if log else "")
    TRANSITIONS.append(line); print(line, flush=True)
    with open(STATUS, "w") as f:
        f.write("PIPELINE STATUS  (rewritten at every transition; elapsed %.1f h)\n" % ((time.time() - T0) / 3600))
        f.write("To watch the running stage, open the log named on the last line.\n")
        f.write("To stop after the current stage, create a file named STOP.txt in this folder.\n\n")
        f.write("\n".join(TRANSITIONS) + "\n")

def stop_requested():
    return os.path.exists(os.path.join(HERE, "STOP.txt"))

def run(stage, script, env_extra, cwd, logname, must_succeed=True):
    """run a script as a subprocess with extra env; stream output to logs/<logname>.txt; return exit code"""
    env = dict(os.environ); env.update({k: str(v) for k, v in env_extra.items()}); env["CORES"] = CORES
    log = os.path.join(LOGS, logname + ".txt")
    status(stage, f"started {script} in {os.path.relpath(cwd, HERE) or '.'}", log)
    t = time.time()
    with open(log, "w") as f:
        f.write(f"# {now()}  {script}  env: {env_extra}\n"); f.flush()
        p = subprocess.Popen([sys.executable, script], cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT)
        rc = p.wait()
    status(stage, f"finished {script} exit={rc} in {(time.time()-t)/3600:.2f} h", log)
    if rc != 0 and must_succeed:
        status(stage, "FAILED - pipeline stopped; read the log")
        finish("FAILED at stage " + stage + " (see " + log + ")")
        sys.exit(1)
    return rc

def finish(verdict, extra=""):
    with open(os.path.join(HERE, "pipeline_FINAL.txt"), "w") as f:
        f.write(f"PIPELINE FINAL  {now()}   elapsed {(time.time()-T0)/3600:.2f} h\n\n{verdict}\n\n{extra}\n\nTransitions:\n" + "\n".join(TRANSITIONS) + "\n")
    print("\n" + verdict, flush=True)

def fresh_dir(name, files):
    d = os.path.join(HERE, name)
    if os.path.exists(d): shutil.rmtree(d)
    os.makedirs(d)
    for src, dst in files: shutil.copy(os.path.join(HERE, src) if not os.path.isabs(src) else src, os.path.join(d, dst))
    return d

# --------------------------------------------------------------------------------------- parameters
if TINY:
    BREED = dict(NSH=4, NW=128, NANG=5, EPS=0.3, NSTART=2, MAXIT=3, MAXO=2)
    ROUNDS_A = ROUNDS_B = "1"
    SWEEP = dict(NANG=13, NSH_LIST=6, NW=128, TAUS="0.5,1", EPS_LIST=0.5, MODES="random", SEEDS=0, SAVE_RHO=1, REGISTER="paper")
    RHO_NAME = "rho_nsh6_random0_eps0.5_nang13.json"
    CENSUS = dict(NB_LIST=4, EPS_LIST="0.5", NSTART=1, MAXIT=3, SEED=7, JOINT_TIME_LIMIT=60)
    CENSUS_EPS_KEY = "0.5"
else:
    BREED = dict(NSH=8, NW=1024, NANG=9, EPS=0.10, NSTART=4, MAXIT=15, MAXO=6)
    SWEEP = dict(NANG=161, NSH_LIST=32, NW=2048, TAUS="0.5,1,2,4,8", EPS_LIST=0.01, MODES="random", SEEDS=0, SAVE_RHO=1, REGISTER="paper")
    RHO_NAME = "rho_nsh32_random0_eps0.01_nang161.json"
    CENSUS = dict(NB_LIST=32, EPS_LIST="0.05,0.01", NSTART=12, MAXIT=40, SEED=7, JOINT_TIME_LIMIT=1800)
    CENSUS_EPS_KEY = "0.01"
SMOKE_BREED = dict(NSH=4, NW=128, NANG=5, EPS=0.3, NSTART=2, MAXIT=3, MAXO=2, ROUNDS=1, START="greedy")
SMOKE_SWEEP = dict(NANG=13, NSH_LIST=4, NW=128, TAUS="0.5,1", EPS_LIST=0.5, MODES="random", SEEDS=0, SAVE_RHO=1, REGISTER="paper")
SMOKE_CENSUS = dict(NB_LIST=4, EPS_LIST="0.5", NSTART=1, MAXIT=2, SEED=7, JOINT_TIME_LIMIT=60)

def read_breed(d):
    r = json.load(open(os.path.join(d, "breed_product_results.json")))
    h = r["history"][-1]
    pj, pp = h["joint"], h["product"]
    feas = bool(pp and pp.get("feasible"))
    ratio = (pp["md"] / pj) if (feas and pj) else None
    promising = feas and ratio is not None and ratio <= PROMISING_RATIO
    txt = (f"rounds completed {h['round']}; joint " + (f"{pj/FLOOR:.3f}x floor" if pj else "infeasible") + "; product "
           + (f"{pp['md']/FLOOR:.3f}x floor" + (f" (ratio product/joint {ratio:.2f})" if ratio else "") if feas else
              ("infeasible (slack %.3f)" % pp["viol"] if pp else "LP failure"))
           + f"; PROMISING={'yes' if promising else 'no'} (rule: feasible and ratio <= {PROMISING_RATIO})")
    return promising, txt

def read_census(d):
    stem = os.path.splitext(RHO_NAME)[0]
    r = json.load(open(os.path.join(d, f"product_census_v3_{stem}.json")))
    rows = []
    for rec in r["results"]:
        j = rec["joint"]; j = f"{j:.3f}" if isinstance(j, float) else str(j)
        s = rec.get("semi"); p = rec.get("product")
        rows.append(f"  eps={rec['eps']}: joint {j} | semi {'-' if s is None else f'{s:.3f}'} | product {'-' if p is None else f'{p:.3f}'}")
    last = [rec for rec in r["results"] if str(rec["eps"]) == CENSUS_EPS_KEY]
    p = last[0].get("product") if last else None
    return p, "\n".join(rows)

# --------------------------------------------------------------------------------------- stage 0
missing = [f for f in REQUIRED if not os.path.exists(os.path.join(HERE, f))]
if missing:
    print("missing files:", missing); sys.exit(1)
status("0-SMOKE", f"pipeline start; CORES={CORES}; TINY={TINY}; promising ratio {PROMISING_RATIO}; rounds A={ROUNDS_A} B={ROUNDS_B}")
sm = fresh_dir("smoke", [(f, f) for f in REQUIRED])
run("0-SMOKE", "breed_product_v1.py", SMOKE_BREED, sm, "smoke_breed")
if not os.path.exists(os.path.join(sm, "bred_product_champion.json")):
    status("0-SMOKE", "FAILED - breed wrote no champion"); finish("FAILED: smoke breed wrote no champion"); sys.exit(1)
run("0-SMOKE", "dither_sweep_rho.py", SMOKE_SWEEP, sm, "smoke_sweep")
smoke_rho = "rho_nsh4_random0_eps0.5_nang13.json"
if not os.path.exists(os.path.join(sm, smoke_rho)):
    status("0-SMOKE", "FAILED - sweep wrote no rho file"); finish("FAILED: smoke sweep wrote no rho file"); sys.exit(1)
env = dict(SMOKE_CENSUS);
with open(os.path.join(sm, "run_census.py"), "w") as f:
    f.write(f"import runpy, sys; sys.argv=['product_census_v3.py', {smoke_rho!r}]; runpy.run_path('product_census_v3.py', run_name='__main__')\n")
run("0-SMOKE", "run_census.py", env, sm, "smoke_census")
if not os.path.exists(os.path.join(sm, "product_census_v3_" + os.path.splitext(smoke_rho)[0] + ".json")):
    status("0-SMOKE", "FAILED - census wrote no result"); finish("FAILED: smoke census wrote no result"); sys.exit(1)
status("0-SMOKE", "all three scripts ran end to end at toy scale - pipeline is healthy")

# --------------------------------------------------------------------------------------- helpers for long stages
def breed_stage(name, start, rounds):
    d = fresh_dir(name, [("breed_product_v1.py", "breed_product_v1.py"), ("greedy_results.json", "greedy_results.json")])
    env = dict(BREED); env.update(ROUNDS=rounds, START=start)
    run(name, "breed_product_v1.py", env, d, name)
    promising, txt = read_breed(d)
    status(name, txt)
    return d, promising, txt

def validate_stage(name, breed_dir):
    d = fresh_dir(name, [("flight_time.py", "flight_time.py"), ("dither_sweep_rho.py", "dither_sweep_rho.py"),
                         ("product_census_v3.py", "product_census_v3.py"),
                         (os.path.join(breed_dir, "bred_product_champion.json"), "greedy_results.json")])
    run(name, "dither_sweep_rho.py", SWEEP, d, name + "_sweep")
    if not os.path.exists(os.path.join(d, RHO_NAME)):
        status(name, "FAILED - no rho file from the full-scale sweep"); return None, "no rho file"
    with open(os.path.join(d, "run_census.py"), "w") as f:
        f.write(f"import runpy, sys; sys.argv=['product_census_v3.py', {RHO_NAME!r}]; runpy.run_path('product_census_v3.py', run_name='__main__')\n")
    run(name, "run_census.py", CENSUS, d, name + "_census")
    p, table = read_census(d)
    status(name, "full-scale product-census table:\n" + table)
    return p, table

# --------------------------------------------------------------------------------------- stages A, B, validation
if stop_requested(): finish("STOPPED by STOP.txt before stage A"); sys.exit(0)
dA, promA, txtA = breed_stage("A_breed_greedy", "greedy", ROUNDS_A)
if promA:
    if stop_requested(): finish("STOPPED by STOP.txt before validation of A", txtA); sys.exit(0)
    p, table = validate_stage("VA_validate_greedy", dA)
    verdict = ("VERDICT: a product-friendly terrain bred from the greedy start " +
               (f"attains MD/floor = {p:.3f} with a station-local census at full scale (eps={CENSUS_EPS_KEY})." if p is not None
                else f"does NOT reproduce the singlet with a station-local census at full scale (eps={CENSUS_EPS_KEY}) - breeding-scale promise did not survive."))
    finish(verdict, "Stage A: " + txtA + "\n\nFull-scale table:\n" + (table or ""))
    sys.exit(0)
if stop_requested(): finish("STOPPED by STOP.txt before stage B", "Stage A: " + txtA); sys.exit(0)
dB, promB, txtB = breed_stage("B_breed_minimal", "minimal", ROUNDS_B)
if promB:
    if stop_requested(): finish("STOPPED by STOP.txt before validation of B", "A: " + txtA + "\nB: " + txtB); sys.exit(0)
    p, table = validate_stage("VB_validate_minimal", dB)
    verdict = ("VERDICT: a product-friendly terrain bred from the minimal start " +
               (f"attains MD/floor = {p:.3f} with a station-local census at full scale (eps={CENSUS_EPS_KEY})." if p is not None
                else f"does NOT reproduce the singlet with a station-local census at full scale (eps={CENSUS_EPS_KEY}) - breeding-scale promise did not survive."))
    finish(verdict, "Stage A: " + txtA + "\nStage B: " + txtB + "\n\nFull-scale table:\n" + (table or ""))
    sys.exit(0)
finish("VERDICT: greedy terrain search finds no station-local terrain (neither the greedy nor the minimal start "
       f"produced a feasible product census within {PROMISING_RATIO}x of the joint price at breeding scale).",
       "Stage A: " + txtA + "\nStage B: " + txtB)
