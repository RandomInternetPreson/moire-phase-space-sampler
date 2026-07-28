"""
Minimum measurement dependence (Hall's fraction M/2) that ANY local
deterministic model must spend to reproduce n-partite GHZ/Mermin
statistics exactly. Exact linear program: n=2 (CHSH at Tsirelson),
n=3 (should print 0.3333, confirming the hand proof), n=4 (the test).
"""
import itertools
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix

def solve_min_dependence(n, settings, targets):
    strategies = list(itertools.product([1, -1], repeat=2 * n))
    nL, nS = len(strategies), len(settings)
    pairs = list(itertools.combinations(range(nS), 2))
    nP = len(pairs)

    def out(lam, i, s):          # party i's outcome under setting s
        return lam[2 * i + s[i]]
    def r(s, l):  return s * nL + l                 # rho vars
    def dv(p, l): return nS * nL + p * nL + l       # |diff| vars
    T = nS * nL + nP * nL                           # index of t
    nvar = T + 1

    subsets = [sub for k in range(1, n)
               for sub in itertools.combinations(range(n), k)]

    A_eq = lil_matrix((nS * (2 + len(subsets)), nvar))
    b_eq = np.zeros(A_eq.shape[0]); row = 0
    for si, s in enumerate(settings):
        for l in range(nL):                          # normalization
            A_eq[row, r(si, l)] = 1.0
        b_eq[row] = 1.0; row += 1
        for l, lam in enumerate(strategies):         # full correlator
            A_eq[row, r(si, l)] = np.prod([out(lam, i, s) for i in range(n)])
        b_eq[row] = targets[si]; row += 1
        for sub in subsets:                          # marginals = 0
            for l, lam in enumerate(strategies):
                A_eq[row, r(si, l)] = np.prod([out(lam, i, s) for i in sub])
            row += 1

    A_ub = lil_matrix((2 * nP * nL + nP, nvar))
    b_ub = np.zeros(A_ub.shape[0]); row = 0
    for pi, (s1, s2) in enumerate(pairs):
        for l in range(nL):
            A_ub[row, r(s1, l)] = 1;  A_ub[row, r(s2, l)] = -1
            A_ub[row, dv(pi, l)] = -1; row += 1
            A_ub[row, r(s2, l)] = 1;  A_ub[row, r(s1, l)] = -1
            A_ub[row, dv(pi, l)] = -1; row += 1
        for l in range(nL):
            A_ub[row, dv(pi, l)] = 1
        A_ub[row, T] = -1; row += 1                  # sum|diff| <= t

    c = np.zeros(nvar); c[T] = 1.0
    res = linprog(c, A_ub=A_ub.tocsr(), b_ub=b_ub,
                  A_eq=A_eq.tocsr(), b_eq=b_eq,
                  bounds=[(0, None)] * nvar, method="highs")
    assert res.success, res.message
    return res.x[T] / 2.0

S2 = np.sqrt(2)
cases = [
    (2, [(0,0),(0,1),(1,0),(1,1)], [1/S2, 1/S2, 1/S2, -1/S2]),
    (3, [(0,0,0),(0,1,1),(1,0,1),(1,1,0)], [1,-1,-1,-1]),
    (4, [s for s in itertools.product([0,1], repeat=4) if sum(s)%2 == 0],
        [1 if sum(s)%4 == 0 else -1
         for s in itertools.product([0,1], repeat=4) if sum(s)%2 == 0]),
]
for n, ss, tt in cases:
    f = solve_min_dependence(n, ss, tt)
    print(f"n={n}: required = {f:.4f}   additive supply n*(9/64) = {n*9/64:.4f}")