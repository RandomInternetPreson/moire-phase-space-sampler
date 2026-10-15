"""c6_above_certify_v3.py (v3: asserts, on a full run with START=0 STRIDE=1, that all 112 cores are certified on both sides at exactly 1/3) (v2: subprocess no longer uses python -I, which hid user-site packages; stderr shown on failure) -- for every C6 five-core whose faithful LP value exceeds
1/(k-1) (list from the sweep logs, one per line: 'k +core,+core,...'), run the
tightened construction (exact upper bound) and the rational dual certificate
(exact lower bound).  Usage: python c6_above_certify_v3.py c6_above_cores.txt START STRIDE"""
import sys, subprocess, importlib.util, os, functools, time
from fractions import Fraction
print = functools.partial(print, flush=True)
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('cc', os.path.join(here, 'cluster_cores_v4.py')); cc = importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)
src = open(os.path.join(here, 'core_tight_bound_v3.py')).read(); src = src[:src.index("name = sys.argv[1]")]
ns = {'__file__': os.path.join(here, 'core_tight_bound_v3.py')}; exec(src, ns); tight = ns['tight_certify']
n, edges = cc.GRAPHS['C6 ring']; grp = cc.stabilizer_group(n, edges)
lines = [l.split() for l in open(sys.argv[1]) if l.strip()]
start, stride = int(sys.argv[2]), int(sys.argv[3])
n_both = 0
for idx in range(start, len(lines), stride):
    k, corestr = int(lines[idx][0]), lines[idx][1]
    core = []
    for tok in corestr.split(','):
        sg = -1 if tok[0] == '-' else 1; x, z = cc.pauli_from_str(tok.lstrip('+-')); assert grp[(x, z)] == sg; core.append((sg, x, z))
    t0 = time.time()
    m_t, lower_t, feas, faithful, F = tight(n, grp, core)
    proc = subprocess.run([sys.executable, os.path.join(here, 'core_dual_certificate_v2.py'), 'C6 ring', str(k), '--core', corestr],
                          capture_output=True, text=True)
    cert = [l for l in proc.stdout.splitlines() if l.startswith('CERTIFIED') or l.startswith('no exact')]
    msg = cert[0] if cert else 'ERROR ' + (proc.stderr.strip().splitlines() or ['(no output)'])[-1]
    print(f"idx {idx} {corestr}: tightened faithful={faithful} F={F}; dual: {msg}   [{time.time()-t0:.0f}s]")
    if faithful and F == Fraction(1, 3) and msg.startswith("CERTIFIED: F_min >= 1/3 "): n_both += 1
print(f"cores certified on both sides at exactly 1/3: {n_both}")
if start == 0 and stride == 1:
    assert len(lines) == 112 and n_both == 112, (len(lines), n_both)
    print("ALL 112 ASSERTED: faithful minimum exactly 1/3 on every core in the list")
