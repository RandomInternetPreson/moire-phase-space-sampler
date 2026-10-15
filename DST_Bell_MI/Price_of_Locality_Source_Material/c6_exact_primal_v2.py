"""c6_exact_primal_v2.py (v2: denominator ladder to 10^8, and the interior-point solution is tried when the dual-simplex vertex does not rationalize; two five-cores needed this) -- exact upper-bound certificates for the C6 cores whose
faithful minimum equals the Theorem-9 lower bound 1/(k-1) by LP value only.

For each such core: solve the collapsed-space faithful LP with the dual simplex
(a vertex solution), rationalize the solution with a ladder of denominators up
to 10^6, verify every constraint, nonnegativity and the value F = 1/(k-1) in
exact rational arithmetic.  A verified solution is an exact faithful model at
the lower bound, so F_min = 1/(k-1) exactly.
Usage: python c6_exact_primal_v2.py <sweep-log-dir> START STRIDE
  reads the 20 c6_lp_*.log files, selects cores with LP value == 1/(k-1) and
  exact-upper == None, and certifies every STRIDE-th starting at START.
Asserts, at the end of a full run (START=0, STRIDE=1), 28 four-cores and 132
five-cores certified.
"""
import sys, os, re, glob, time, itertools, functools, importlib.util
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
print = functools.partial(print, flush=True)
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('cc', os.path.join(here, 'cluster_cores_v4.py')); cc = importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)

def build(n, grp, core):
    k = len(core)
    measured = [set() for _ in range(n)]
    for (sg, x, z) in core:
        for i in range(n):
            xi, zi = (x >> i) & 1, (z >> i) & 1
            pp = {(1,0):0,(1,1):1,(0,1):2}.get((xi, zi))
            if pp is not None: measured[i].add(pp)
    slots = [(i, p) for i in range(n) for p in sorted(measured[i])]; pos = {sp: j for j, sp in enumerate(slots)}
    S = 1 << len(slots); strat = np.arange(S, dtype=np.int64)
    def parity(x, z):
        p = np.zeros(S, dtype=np.int64)
        for i in range(n):
            xi, zi = (x >> i) & 1, (z >> i) & 1
            pp = {(1,0):0,(1,1):1,(0,1):2}.get((xi, zi))
            if pp is not None: p ^= (strat >> pos[(i, pp)]) & 1
        return p
    rows_per = []
    for (sg, x, z) in core:
        supp = [i for i in range(n) if ((x | z) >> i) & 1]
        rows = [(np.ones(S, dtype=np.int64), 1)]
        for r in range(1, len(supp) + 1):
            for T in itertools.combinations(supp, r):
                mask = sum(1 << i for i in T); xT, zT = x & mask, z & mask
                qv = sg if r == len(supp) else grp.get((xT, zT), 0)
                rows.append((1 - 2 * parity(xT, zT), qv))
        rows_per.append(rows)
    return k, S, rows_per

def solve_and_certify(n, grp, core, target):
    k, S, rows_per = build(n, grp, core)
    pairs = list(itertools.combinations(range(k), 2))
    nrho = k * S; nu = len(pairs) * S; Tv = nrho + nu; nvar = Tv + 1
    er, ec, ev, beq = [], [], [], []; row = 0
    for si, rows in enumerate(rows_per):
        for vec, b in rows:
            nz = np.nonzero(vec)[0]; er += [row] * len(nz); ec += (nz + si * S).tolist(); ev += vec[nz].tolist(); beq.append(float(b)); row += 1
    A_eq = coo_matrix((ev, (er, ec)), shape=(row, nvar)).tocsr()
    ur, uc, uv, bub = [], [], [], []; row = 0
    for pi, (i, j) in enumerate(pairs):
        base = nrho + pi * S
        for c in range(S):
            for sgn in (1.0, -1.0):
                ur += [row] * 3; uc += [i * S + c, j * S + c, base + c]; uv += [sgn, -sgn, -1.0]; bub.append(0.0); row += 1
        ur += [row] * (S + 1); uc += list(range(base, base + S)) + [Tv]; uv += [1.0] * S + [-1.0]; bub.append(0.0); row += 1
    A_ub = coo_matrix((uv, (ur, uc)), shape=(row, nvar)).tocsr()
    cvec = np.zeros(nvar); cvec[Tv] = 1.0
    F_lp = None
    for method in ("highs-ds", "highs-ipm"):
        res = linprog(cvec, A_ub=A_ub, b_ub=np.array(bub), A_eq=A_eq, b_eq=np.array(beq), bounds=(0, None), method=method)
        assert res.success, res.message
        F_lp = res.x[Tv] / 2
        for den in (6, 12, 24, 48, 120, 240, 480, 1200, 5000, 20000, 100000, 1000000, 10**7, 10**8):
            rho = [[Fraction(float(v)).limit_denominator(den) for v in res.x[si * S:(si + 1) * S]] for si in range(k)]
            if any(r < 0 for rs in rho for r in rs): continue
            ok = True
            for si, rows in enumerate(rows_per):
                for vec, b in rows:
                    if sum(Fraction(int(vv)) * r for vv, r in zip(vec, rho[si]) if r != 0) != b: ok = False; break
                if not ok: break
            if not ok: continue
            F = max(sum(abs(a - b) for a, b in zip(rho[i], rho[j])) / 2 for i, j in pairs)
            if F == target: return F_lp, f"{den} ({method})", True
            break   # an exact feasible model with a larger F: try the other method
    return F_lp, None, False

if __name__ == "__main__":
    logdir, start, stride = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    rows = []
    for f in glob.glob(os.path.join(logdir, 'c6_lp_*.log')):
        for line in open(f):
            m = re.match(r"idx (\d+) size (\d+) core (.*?): LP = (\S+) \(([\d.]+)\), exact-upper = (\S+), strategies = (\d+)", line)
            if m: rows.append((int(m[1]), int(m[2]), m[3], Fraction(m[4]), m[6]))
    todo = sorted(r for r in rows if r[3] == Fraction(1, r[1] - 1) and r[4] == 'None')
    print(f"{len(todo)} LP-only cores at 1/(k-1): sizes {dict((k, sum(1 for r in todo if r[1]==k)) for k in (4,5))}; this process: {start}::{stride}")
    n, edges = cc.GRAPHS['C6 ring']; grp = cc.stabilizer_group(n, edges)
    done = {4: 0, 5: 0}
    for j in range(start, len(todo), stride):
        idx, k, corestr, lp, _ = todo[j]
        core = []
        for tok in corestr.split():
            sg = -1 if tok[0] == '-' else 1; x, z = cc.pauli_from_str(tok.lstrip('+-')); assert grp[(x, z)] == sg; core.append((sg, x, z))
        t0 = time.time()
        F_lp, den, ok = solve_and_certify(n, grp, core, Fraction(1, k - 1))
        if ok is True: done[k] += 1
        print(f"idx {idx} size {k} core {corestr}: LP {F_lp:.9f}; exact rational model at 1/(k-1): {ok} (denominator {den})   [{time.time()-t0:.0f}s]")
    print("certified exactly:", done)
    if start == 0 and stride == 1:
        assert done == {4: 28, 5: 132}, done
        print("ALL 160 LP-ONLY C6 CORES CERTIFIED EXACTLY AT 1/(k-1)")
