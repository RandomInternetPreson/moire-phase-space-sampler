"""
Symmetry-reduced LP for the minimum measurement dependence (Hall fraction M/2)
of n-partite GHZ/Mermin statistics.

Reduction: a deterministic strategy lambda in {+-1}^{2n} only enters every
full-correlator constraint through (P, gamma) in GF(2)^{n+1}, where
P = parity of the X-outcomes and gamma_i = (X_i != Y_i). Averaging each
distribution uniformly over the 2^{n-1}-element fiber of each (P, gamma)
class (a) preserves all full correlators, (b) makes every proper-subset
correlator exactly zero automatically, (c) does not increase any pairwise
L1 distance. Hence the LP over 4^n strategies is exactly equivalent to an
LP over 2^{n+1} classes with NO marginal constraints.

Usage: python mermin_lp_reduced.py 7
"""
import sys, time, itertools
import numpy as np
from scipy.sparse import coo_matrix
from scipy.optimize import linprog

def build_and_solve(n, verbose=True):
    if n == 2:
        settings = [(0,0),(0,1),(1,0),(1,1)]
        targets  = [1/np.sqrt(2)]*3 + [-1/np.sqrt(2)]
    else:
        settings = [s for s in itertools.product([0,1], repeat=n)
                    if sum(s) % 2 == 0]
        targets  = [1.0 if sum(s) % 4 == 0 else -1.0 for s in settings]
    nS = len(settings)
    effs = np.array(list(itertools.product([0,1], repeat=n+1)), dtype=np.int8)
    nL = len(effs)                      # 2^(n+1) classes: (P, g_1..g_n)
    P = effs[:, 0]
    G = effs[:, 1:]
    sett = np.array(settings, dtype=np.int8)
    # sign[s, c] = (-1)^(P + sum_{i in S} g_i)
    bits = (P[None, :] + sett @ G.T) % 2
    sign = np.where(bits == 1, -1.0, 1.0)

    pairs = list(itertools.combinations(range(nS), 2)); nP = len(pairs)
    def rho(si, l):  return si*nL + l
    def dvar(pi, l): return nS*nL + pi*nL + l
    Tidx = nS*nL + nP*nL
    nvar = Tidx + 1

    er, ec, ev, b_eq = [], [], [], []
    row = 0
    for si in range(nS):
        er += [row]*nL; ec += list(range(rho(si,0), rho(si,0)+nL)); ev += [1.0]*nL
        b_eq.append(1.0); row += 1
        er += [row]*nL; ec += list(range(rho(si,0), rho(si,0)+nL)); ev += list(sign[si])
        b_eq.append(targets[si]); row += 1
    A_eq = coo_matrix((ev, (er, ec)), shape=(row, nvar)).tocsr()

    ur, uc, uv, b_ub = [], [], [], []
    row = 0
    for pi, (s1, s2) in enumerate(pairs):
        for l in range(nL):
            ur += [row]*3; uc += [rho(s1,l), rho(s2,l), dvar(pi,l)]
            uv += [1.0, -1.0, -1.0]; b_ub.append(0.0); row += 1
            ur += [row]*3; uc += [rho(s2,l), rho(s1,l), dvar(pi,l)]
            uv += [1.0, -1.0, -1.0]; b_ub.append(0.0); row += 1
        ur += [row]*(nL+1)
        uc += list(range(dvar(pi,0), dvar(pi,0)+nL)) + [Tidx]
        uv += [1.0]*nL + [-1.0]; b_ub.append(0.0); row += 1
    A_ub = coo_matrix((uv, (ur, uc)), shape=(row, nvar)).tocsr()

    c = np.zeros(nvar); c[Tidx] = 1.0
    t0 = time.time()
    res = linprog(c, A_ub=A_ub, b_ub=np.array(b_ub),
                  A_eq=A_eq, b_eq=np.array(b_eq),
                  bounds=(0, None), method="highs")
    assert res.success, res.message
    val = res.x[Tidx] / 2.0
    if verbose:
        print(f"n={n}: settings={nS}, classes={nL}, vars={nvar:,}, "
              f"solve={time.time()-t0:.1f}s -> required = {val:.6f}")
    return val

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [3, 4, 5, 6, 7]
    for n in ns:
        build_and_solve(n)