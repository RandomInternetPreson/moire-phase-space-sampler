"""
Minimum measurement dependence (Hall fraction M/2) for n-partite
GHZ/Mermin statistics. Usage: python mermin_lp.py 5
"""
import sys, time, itertools
import numpy as np
from scipy.sparse import coo_matrix
from scipy.optimize import linprog

def solve(n):
    if n == 2:
        settings = [(0,0),(0,1),(1,0),(1,1)]
        targets  = [1/np.sqrt(2)]*3 + [-1/np.sqrt(2)]
    else:
        settings = [s for s in itertools.product([0,1], repeat=n)
                    if sum(s) % 2 == 0]
        targets  = [1.0 if sum(s) % 4 == 0 else -1.0 for s in settings]
    strategies = list(itertools.product([1,-1], repeat=2*n))
    strat = np.array(strategies)
    nL, nS = len(strategies), len(settings)
    pairs = list(itertools.combinations(range(nS), 2)); nP = len(pairs)
    outc = [[strat[:, 2*i + s[i]] for i in range(n)] for s in settings]

    def rho(si, l):  return si*nL + l
    def dvar(pi, l): return nS*nL + pi*nL + l
    Tidx = nS*nL + nP*nL
    nvar = Tidx + 1
    subsets = [sub for k in range(1, n)
               for sub in itertools.combinations(range(n), k)]

    er, ec, ev, b_eq, row = [], [], [], [], 0
    cols_all = list(range(nL))
    for si in range(nS):
        er += [row]*nL; ec += [rho(si,l) for l in cols_all]; ev += [1.0]*nL
        b_eq.append(1.0); row += 1
        prod = np.ones(nL)
        for i in range(n): prod = prod * outc[si][i]
        er += [row]*nL; ec += [rho(si,l) for l in cols_all]; ev += list(prod)
        b_eq.append(targets[si]); row += 1
        for sub in subsets:
            p = np.ones(nL)
            for i in sub: p = p * outc[si][i]
            er += [row]*nL; ec += [rho(si,l) for l in cols_all]; ev += list(p)
            b_eq.append(0.0); row += 1
    A_eq = coo_matrix((ev,(er,ec)), shape=(row, nvar)).tocsr()

    ur, uc, uv, b_ub, row = [], [], [], [], 0
    for pi,(s1,s2) in enumerate(pairs):
        for l in range(nL):
            ur += [row]*3; uc += [rho(s1,l), rho(s2,l), dvar(pi,l)]
            uv += [1.0,-1.0,-1.0]; b_ub.append(0.0); row += 1
            ur += [row]*3; uc += [rho(s2,l), rho(s1,l), dvar(pi,l)]
            uv += [1.0,-1.0,-1.0]; b_ub.append(0.0); row += 1
        ur += [row]*(nL+1)
        uc += [dvar(pi,l) for l in cols_all] + [Tidx]
        uv += [1.0]*nL + [-1.0]; b_ub.append(0.0); row += 1
    A_ub = coo_matrix((uv,(ur,uc)), shape=(row, nvar)).tocsr()

    c = np.zeros(nvar); c[Tidx] = 1.0
    t0 = time.time()
    res = linprog(c, A_ub=A_ub, b_ub=np.array(b_ub),
                  A_eq=A_eq, b_eq=np.array(b_eq),
                  bounds=(0, None), method="highs")
    assert res.success, res.message
    print(f"n={n}: settings={nS}, vars={nvar:,}, "
          f"solve={time.time()-t0:.1f}s -> required = {res.x[Tidx]/2:.6f}   "
          f"(1/3 = {1/3:.6f};  n*(9/64) = {n*9/64:.4f})")

if __name__ == "__main__":
    solve(int(sys.argv[1]) if len(sys.argv) > 1 else 5)