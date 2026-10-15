"""core_tight_bound_v3.py (v3: optional 4th argument P/Q asserts that the tightened construction is faithful with F = P/Q on every core run) -- for cores where the ground-state construction is not
faithful: (a) Theorem 3 with each setting's sector tightened to ALL its
deterministic-target constraints (full product and the +-1 marginals), giving
an exact lower bound m'-based; (b) the 'tightened' ground-state construction:
uniform on strategies minimizing the number of tightened sectors missed,
within the tightened sector of s -- checked for faithfulness and F exactly.
Usage: python core_tight_bound_v3.py ["C6 ring"|"C4 chain"|"C5 ring"] [sizes e.g. 4,5 or 6,7] [sample size, default 40] [expected F, e.g. 1/4]"""
import itertools, functools, importlib.util, sys
from fractions import Fraction
import numpy as np
print = functools.partial(print, flush=True)
import os
spec = importlib.util.spec_from_file_location('cc', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cluster_cores_v4.py')); cc = importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)

def tight_certify(n, group, core):
    k = len(core); S = 1 << (3 * n); strat = np.arange(S, dtype=np.int64)
    def parity(x, z):
        p = np.zeros(S, dtype=np.int64)
        for i in range(n):
            xi, zi = (x >> i) & 1, (z >> i) & 1
            pp = {(1,0):0,(1,1):1,(0,1):2}.get((xi, zi))
            if pp is not None: p ^= (strat >> (3 * i + pp)) & 1
        return p
    A = []; zero_rows = []
    for (sg, x, z) in core:
        supp = [i for i in range(n) if ((x | z) >> i) & 1]
        ok = parity(x, z) == (0 if sg == 1 else 1)
        zr = []
        for r in range(1, len(supp)):
            for T in itertools.combinations(supp, r):
                mask = sum(1 << i for i in T); xT, zT = x & mask, z & mask
                qv = group.get((xT, zT), 0)
                if qv != 0: ok = ok & (parity(xT, zT) == (0 if qv == 1 else 1))
                else: zr.append(parity(xT, zT))
        A.append(ok); zero_rows.append(zr)
    A = np.array(A); miss = (~A).sum(axis=0)
    m_t = k - int(miss.min())
    lower_t = Fraction(k - m_t, k - 1)
    ground = np.nonzero(miss == miss.min())[0]
    G = A[:, ground]; g = G.sum(axis=1)
    feas = bool((g > 0).all()); faithful = feas; F = None
    if feas:
        for si in range(k):
            idx = ground[G[si]]
            for par in zero_rows[si]:
                if 2 * int(par[idx].sum()) != len(idx): faithful = False
        F = Fraction(0)
        for i, j in itertools.combinations(range(k), 2):
            inter = int((G[i] & G[j]).sum()); gi, gj = int(g[i]), int(g[j])
            d = Fraction(inter) * abs(Fraction(1, gi) - Fraction(1, gj)) + Fraction(gi - inter, gi) + Fraction(gj - inter, gj)
            F = max(F, d / 2)
    return m_t, lower_t, feas, faithful, F

name = sys.argv[1] if len(sys.argv) > 1 else 'C6 ring'
sizes = [int(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [4, 5]
sample = int(sys.argv[3]) if len(sys.argv) > 3 else 40
expect = Fraction(sys.argv[4]) if len(sys.argv) > 4 else None
n_checked = 0
n, edges = cc.GRAPHS[name]; grp = cc.stabilizer_group(n, edges)
found = cc.minimal_contradictions(n, grp, max(sizes))
rng = np.random.default_rng(2026)
for k in sizes:
    cores = found[k]; todo = cores if sample >= len(cores) else [cores[i] for i in rng.choice(len(cores), sample, replace=False)]
    for c in todo:
        r = cc.certify_core(n, grp, c)
        if r['exact']: continue
        m_t, lower_t, feas, faithful, F = tight_certify(n, grp, c)
        F_lp, F_ex, Sz = cc.faithful_lp_core(n, grp, c)
        print(f"size {k} core {', '.join(('+' if sg>0 else '-')+cc.pauli_to_str(x,z,n) for (sg,x,z) in c)}: "
              f"plain construction faithful={r['faithful']} F={r['F']}; tightened: m'={m_t} lower={lower_t}, feasible={feas}, faithful={faithful}, F={F}; "
              f"LP={Fraction(F_lp).limit_denominator(400)} exact-upper={F_ex}")
        n_checked += 1
        if expect is not None:
            assert faithful and F == expect, ("tightened construction not faithful at the expected value", faithful, F, expect)
if expect is not None:
    print(f"assertion passed on all {n_checked} cores where the plain construction is not faithful: tightened faithful with F = {expect}")
