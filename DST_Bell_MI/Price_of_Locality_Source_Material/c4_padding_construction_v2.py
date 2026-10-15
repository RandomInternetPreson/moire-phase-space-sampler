"""
c4_padding_construction_v2.py (v2: asserts the quoted split 43 faithful at 1/3 / 8 not faithful / 4 faithful at 1/2 or 1, and s = 5/6 on all 55)
c4_padding_construction_v1.py -- the C4 AVN core padded with two further
stabilizers (all 55 pairs): is the ground-state construction feasible (all
marginal correlators at their stabilizer values) and what is its F?
Since F_min can only rise when constraints are added, F_min(padded) >= 1/3
(the core's certified value); a feasible construction with F = 1/3 is
therefore a proof that F_min(padded) = 1/3 and that the construction is
optimal -- no LP required.  Exact: counts and Fractions only.
"""
import itertools, functools
from fractions import Fraction
import numpy as np
print = functools.partial(print, flush=True)

PAULI = {'I': 0, 'X': 1, 'Y': 2, 'Z': 3}
def pmul(a, b):
    if a == 0: return 1, b
    if b == 0: return 1, a
    if a == b: return 1, 0
    cyc = {(1, 2): (1j, 3), (2, 3): (1j, 1), (3, 1): (1j, 2)}
    if (a, b) in cyc: return cyc[(a, b)]
    ph, p = cyc[(b, a)]; return -ph, p
def smul(s1, s2):
    sign = s1[0] * s2[0]; out = []
    for a, b in zip(s1[1], s2[1]):
        ph, p = pmul(a, b); sign *= ph; out.append(p)
    return sign, tuple(out)
def group(gens):
    g = {}
    for bits in itertools.product((0, 1), repeat=len(gens)):
        el = (1, tuple(0 for _ in gens[0][1]))
        for b, gq in zip(bits, gens):
            if b: el = smul(el, gq)
        g[el[1]] = int(round(el[0].real)) if isinstance(el[0], complex) else el[0]
    return g
pstr = lambda t: ''.join('IXYZ'[p] for p in t)

n = 4
grp = group([(1, (1, 3, 0, 0)), (1, (3, 1, 3, 0)), (1, (0, 3, 1, 3)), (1, (0, 0, 3, 1))])
core = ['IZXZ', 'IZYY', 'ZYYZ', 'ZYXY']
others = [pstr(t) for t in grp if t != (0,)*4 and pstr(t) not in core]
strat = np.array(list(itertools.product(*[list(itertools.product((1, -1), repeat=3))] * n)))  # [4096,4,3]

def analyse(names):
    S = len(names); full = []; margs = []
    for name in names:
        t = tuple(PAULI[ch] for ch in name); sg = grp[t]
        parties = [i for i in range(n) if t[i]]
        prod = np.ones(len(strat))
        for i in parties: prod = prod * strat[:, i, t[i]-1]
        full.append(prod == sg)
        rows = []
        for r in range(1, len(parties)):
            for T in itertools.combinations(parties, r):
                sub = tuple(t[i] if i in T else 0 for i in range(n))
                pr = np.ones(len(strat))
                for i in T: pr = pr * strat[:, i, t[i]-1]
                rows.append((pr, grp.get(sub, 0)))
        margs.append(rows)
    sat = np.array(full); H = (~sat).sum(axis=0)
    s = Fraction(S - int(H.min()), S)
    ground = np.nonzero(H == H.min())[0]
    G = sat[:, ground]; g = G.sum(axis=1)
    if not (g > 0).all(): return s, False, False, None
    marg_ok = True
    for si in range(S):
        idx = ground[G[si]]
        for pr, qv in margs[si]:
            # exact: sum over support of integer +-1 values / g must equal qv
            if int(pr[idx].sum()) != qv * int(g[si]): marg_ok = False
    F = Fraction(0)
    for i, j in itertools.combinations(range(S), 2):
        inter = int((G[i] & G[j]).sum()); gi, gj = int(g[i]), int(g[j])
        d = Fraction(inter) * abs(Fraction(1, gi) - Fraction(1, gj)) + Fraction(gi - inter, gi) + Fraction(gj - inter, gj)
        F = max(F, d / 2)
    return s, True, marg_ok, F

s0, fs0, mo0, F0 = analyse(core)
print(f"core: s={s0}, support-feasible={fs0}, marginals exact={mo0}, construction F={F0}")
res = []
for a, b in itertools.combinations(others, 2):
    s, fs, mo, F = analyse(core + [a, b]); res.append((a, b, s, fs, mo, F))
from collections import Counter
print("padded pairs:", len(res))
print("s values:", dict(Counter(r[2] for r in res)))
print("support-feasible:", sum(r[3] for r in res), " marginals exact:", sum(r[4] for r in res))
print("construction F values:", dict(Counter(str(r[5]) for r in res)))
bad = [r for r in res if not (r[3] and r[4] and r[5] == Fraction(1, 3))]
print("pairs where construction is NOT (feasible with F = 1/3):", [(r[0], r[1], str(r[5]), r[3], r[4]) for r in bad] or "none")

assert len(res) == 55 and all(r[2] == Fraction(5, 6) for r in res)
n_good = sum(1 for r in res if r[3] and r[4] and r[5] == Fraction(1, 3))
n_marg = sum(1 for r in res if r[3] and not r[4])
n_bigF = sum(1 for r in res if r[3] and r[4] and r[5] in (Fraction(1, 2), Fraction(1)))
assert (n_good, n_marg, n_bigF) == (43, 8, 4), (n_good, n_marg, n_bigF)
assert s0 == Fraction(3, 4) and fs0 and mo0 and F0 == Fraction(1, 3)
print("ASSERTED: core exact at 1/3; padded: 43 faithful at 1/3, 8 not faithful, 4 faithful at 1/2 or 1; s = 5/6 on all 55")
