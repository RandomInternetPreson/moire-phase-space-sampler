"""
quantumtestv5.py -- symmetry-reduced exact LP for Mermin measurement-
dependence floors.

Two exact reductions stacked:
 1. Class-fiber reduction (v4): strategies -> (P, gamma) in GF(2)^(n+1).
 2. NEW: S_n symmetrization. The scenario (settings, targets, constraints,
    objective) is invariant under party permutations; the feasible set is
    convex and the objective T = max_{pairs} L1 is convex and invariant, so
    group-averaging any optimal solution yields a symmetric optimal solution.
    Restricting to S_n-covariant densities is therefore EXACT (WLOG), not an
    approximation.

Covariant variables: x[k, P, a, b] where k = |s| (even), a = |gamma ^ s|,
b = |gamma ^ s-bar|. Pair orbits: (k, k', m = |s ^ s'|). Cell weights are
products of binomial coefficients. Variable count is O(n^3): n = 9 needs
~200 density variables instead of 33,000,000.

Usage: python quantumtestv5.py 3 4 5 6 7 8 9 10 11
"""
import sys, time, itertools
from math import comb
import numpy as np
from scipy.sparse import coo_matrix
from scipy.optimize import linprog
from fractions import Fraction

def solve_sym(n, verbose=True):
    ks = [k for k in range(0, n + 1, 2)]                     # even |s|
    tgt = {k: (1.0 if k % 4 == 0 else -1.0) for k in ks}

    # density variables x[k, P, a, b]
    xidx, xs = {}, []
    for k in ks:
        for P in (0, 1):
            for a in range(k + 1):
                for b in range(n - k + 1):
                    xidx[(k, P, a, b)] = len(xs); xs.append((k, P, a, b))
    nx = len(xs)

    # pair orbits (k, k', m): realizable, s != s'
    orbits = []
    for k in ks:
        for k2 in ks:
            if k2 < k: continue
            lo, hi = max(0, k + k2 - n), min(k, k2)
            for m in range(lo, hi + 1):
                if k == k2 and m == k: continue              # s == s'
                orbits.append((k, k2, m))

    # abs-value auxiliaries u[orbit, P, al, be, ga, de]
    uidx, us = {}, []
    for oi, (k, k2, m) in enumerate(orbits):
        for P in (0, 1):
            for al in range(m + 1):
                for be in range(k - m + 1):
                    for ga in range(k2 - m + 1):
                        for de in range(n - k - k2 + m + 1):
                            uidx[(oi, P, al, be, ga, de)] = nx + len(us)
                            us.append((oi, P, al, be, ga, de))
    Tvar = nx + len(us); nvar = Tvar + 1

    er, ec, ev, beq = [], [], [], []
    row = 0
    for k in ks:                                             # per s-orbit
        # normalization and correlator
        for (kind, target) in (("norm", 1.0), ("corr", tgt[k])):
            for P in (0, 1):
                for a in range(k + 1):
                    for b in range(n - k + 1):
                        w = comb(k, a) * comb(n - k, b)
                        sgn = 1.0 if kind == "norm" else (-1.0) ** ((P + a) % 2)
                        er.append(row); ec.append(xidx[(k, P, a, b)]); ev.append(w * sgn)
            beq.append(target); row += 1
    A_eq = coo_matrix((ev, (er, ec)), shape=(row, nvar)).tocsr()

    ur, uc, uv, bub = [], [], [], []
    row = 0
    for oi, (k, k2, m) in enumerate(orbits):
        for P in (0, 1):
            for al in range(m + 1):
                for be in range(k - m + 1):
                    for ga in range(k2 - m + 1):
                        for de in range(n - k - k2 + m + 1):
                            xi = xidx[(k, P, al + be, ga + de)]
                            xj = xidx[(k2, P, al + ga, be + de)]
                            ui = uidx[(oi, P, al, be, ga, de)]
                            for sgn in (1.0, -1.0):
                                ur += [row] * 3; uc += [xi, xj, ui]
                                uv += [sgn, -sgn, -1.0]; bub.append(0.0); row += 1
        # sum of weighted u <= T
        cells = [key for key in uidx if key[0] == oi]
        ur += [row] * (len(cells) + 1)
        for (_, P, al, be, ga, de) in cells:
            w = (comb(m, al) * comb(k - m, be) * comb(k2 - m, ga)
                 * comb(n - k - k2 + m, de))
            uc.append(uidx[(oi, P, al, be, ga, de)]); uv.append(float(w))
        uc.append(Tvar); uv.append(-1.0); bub.append(0.0); row += 1
    A_ub = coo_matrix((uv, (ur, uc)), shape=(row, nvar)).tocsr()

    c = np.zeros(nvar); c[Tvar] = 1.0
    t0 = time.time()
    res = linprog(c, A_ub=A_ub, b_ub=np.array(bub), A_eq=A_eq, b_eq=np.array(beq),
                  bounds=(0, None), method="highs")
    assert res.success, res.message
    F = res.x[Tvar] / 2.0
    frac = Fraction(F).limit_denominator(200)
    R = 2 ** ((n - 1) // 2)
    pred = Fraction(R, 2 * (R + 1))
    tag = "MATCHES staircase R/(2(R+1))" if abs(float(pred) - F) < 1e-8 \
          else f"!= staircase prediction {pred}"
    if verbose:
        print(f"n={n}: sym vars={nvar:,}, orbits={len(orbits)}, "
              f"solve={time.time()-t0:.2f}s -> F = {F:.9f} = {frac}  [{tag}]")
    return F

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [3, 4, 5, 6, 7]
    for n in ns:
        solve_sym(n)
