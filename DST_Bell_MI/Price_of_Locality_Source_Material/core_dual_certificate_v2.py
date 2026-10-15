"""
core_dual_certificate_v2.py -- exact lower bound on the FAITHFUL minimum of F
for one core, by an LP dual certificate verified in rational arithmetic.

The faithful LP (collapsed strategy space, every marginal constrained) is
solved in double precision; its dual solution is rationalized with bounded
denominators and checked EXACTLY for dual feasibility (A_eq^T y + A_ub^T w <= c,
w <= 0, in Fractions).  Weak duality then gives, for every feasible primal x,
c.x >= y.b_eq + w.b_ub, i.e. an exact lower bound on F.  Combined with an exact
faithful model of the same value (e.g. the tightened construction of
core_tight_bound_v1.py or a rationalized primal solution), F_min is certified.

Usage: python core_dual_certificate_v2.py "<graph>" <size> <index> [--expect 1/4]
       python core_dual_certificate_v2.py "<graph>" <size> --core "+IIZYYZ,..." [--expect 1/3]
(v2: --expect P/Q asserts that the certified bound equals P/Q; an assertion failure is a finding)
Graphs: "C4 chain", "C5 ring", "C6 ring" (as in cluster_cores_v3.py).
"""
import importlib.util, itertools, sys, time
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
import os
spec = importlib.util.spec_from_file_location('cc', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cluster_cores_v4.py'))
cc = importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)
expect = None
if '--expect' in sys.argv:
    i = sys.argv.index('--expect'); expect = Fraction(sys.argv[i + 1]); del sys.argv[i:i + 2]
name = sys.argv[1]; ksz = int(sys.argv[2])
n, edges = cc.GRAPHS[name]; grp = cc.stabilizer_group(n, edges)
if len(sys.argv) > 4 and sys.argv[3] == "--core":
    core = []
    for tok in sys.argv[4].split(","):
        tok = tok.strip(); sg = -1 if tok[0] == "-" else 1; s = tok.lstrip("+-")
        x, z = cc.pauli_from_str(s); assert grp[(x, z)] == sg, ("not a stabilizer with that sign", tok); core.append((sg, x, z))
    assert len(core) == ksz
else:
    which = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    found = cc.minimal_contradictions(n, grp, ksz)
    core = found[ksz][which]
print("core:", ', '.join(('+' if sg > 0 else '-') + cc.pauli_to_str(x, z, n) for (sg, x, z) in core))
k = len(core)
measured = [set() for _ in range(n)]
for (sg, x, z) in core:
    for i in range(n):
        xi, zi = (x >> i) & 1, (z >> i) & 1
        if (xi, zi) == (1, 0): measured[i].add(0)
        elif (xi, zi) == (1, 1): measured[i].add(1)
        elif (xi, zi) == (0, 1): measured[i].add(2)
slots = [(i, p) for i in range(n) for p in sorted(measured[i])]
pos = {sp: j for j, sp in enumerate(slots)}
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
pairs = list(itertools.combinations(range(k), 2))
nrho = k * S; nu = len(pairs) * S; Tv = nrho + nu; nvar = Tv + 1
# Build: min t  s.t.  A_eq x = b_eq ;  rho_i(c) - rho_j(c) - u_p(c) <= 0 ; -(rho_i - rho_j) - u_p <= 0 ; sum_c u_p(c) - t <= 0 ; x >= 0
er, ec, ev, beq = [], [], [], []; row = 0
for si, rows in enumerate(rows_per):
    for vec, b in rows:
        nz = np.nonzero(vec)[0]
        er += [row] * len(nz); ec += (nz + si * S).tolist(); ev += vec[nz].tolist(); beq.append(float(b)); row += 1
neq = row
A_eq = coo_matrix((ev, (er, ec)), shape=(neq, nvar)).tocsr()
ur, uc, uv, bub = [], [], [], []; row = 0
for pi, (i, j) in enumerate(pairs):
    base = nrho + pi * S
    for c in range(S):
        for sgn in (1.0, -1.0):
            ur += [row] * 3; uc += [i * S + c, j * S + c, base + c]; uv += [sgn, -sgn, -1.0]; bub.append(0.0); row += 1
    ur += [row] * (S + 1); uc += list(range(base, base + S)) + [Tv]; uv += [1.0] * S + [-1.0]; bub.append(0.0); row += 1
nub = row
A_ub = coo_matrix((uv, (ur, uc)), shape=(nub, nvar)).tocsr()
cvec = np.zeros(nvar); cvec[Tv] = 1.0
t0 = time.time()
res = linprog(cvec, A_ub=A_ub, b_ub=np.array(bub), A_eq=A_eq, b_eq=np.array(beq), bounds=(0, None), method="highs-ds")
print("LP:", res.x[Tv] / 2, f"[{time.time()-t0:.0f}s]")
y = res.eqlin.marginals   # free
w = res.ineqlin.marginals # <= 0
print("dual obj (float):", (np.array(beq) @ y + np.array(bub) @ w) / 2)
# Rationalize duals and verify exactly: need  A_eq^T y + A_ub^T w <= c  componentwise, w <= 0.  Then for any feasible x>=0:
#  c.x >= (A_eq^T y + A_ub^T w).x = y.b_eq + w.(A_ub x) >= y.b_eq + w.b_ub  (since w<=0 and A_ub x <= b_ub).
for den in (6, 12, 24, 48, 60, 120, 240, 480, 1200, 5000, 20000, 100000, 1000000, 10**7):
    yq = [Fraction(float(v)).limit_denominator(den) for v in y]
    wq = [min(Fraction(float(v)).limit_denominator(den), Fraction(0)) for v in w]
    # compute A^T y + A_ub^T w exactly via sparse structure
    red = [Fraction(0)] * nvar
    Aeq = A_eq.tocoo(); Aub = A_ub.tocoo()
    for r, c, v in zip(Aeq.row, Aeq.col, Aeq.data):
        if yq[r] != 0: red[c] += Fraction(int(v)) * yq[r]
    for r, c, v in zip(Aub.row, Aub.col, Aub.data):
        if wq[r] != 0: red[c] += Fraction(int(v)) * wq[r]
    ok = all(red[c] <= Fraction(int(cvec[c])) for c in range(nvar))
    bound = sum(Fraction(int(b)) * yy for b, yy in zip(beq, yq)) + sum(Fraction(int(b)) * ww for b, ww in zip(bub, wq))
    print(f"den {den}: dual feasible exactly: {ok}; certified lower bound on F = {bound/2}")
    if ok:
        print("CERTIFIED: F_min >=", bound / 2, "(exact); 1/(k-1) =", Fraction(1, k - 1), "; LP value =", Fraction(res.x[Tv] / 2).limit_denominator(400))
        if expect is not None:
            assert bound / 2 == expect, ("certified bound differs from --expect", bound / 2, expect)
            print("assertion passed: certified bound ==", expect)
        break
else:
    print("no exact dual certificate found at these denominators")
    assert expect is None, "no certificate found but --expect was given"
