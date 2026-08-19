#!/usr/bin/env python3
"""
price_surface_lp.py - The full price surface M(S, eta).

Minimal measurement dependence for ANY faithful local model to achieve CHSH
value S at symmetric detector efficiency eta (singlet-like phenomenology:
singles rates eta, coincidence rate eta^2, unbiased marginals; the four
correlators free except for their CHSH combination).

Registered conjecture:  M(S, eta) = max{0, eta((S+2) eta - 4)/6}
Edges: S=2sqrt2 -> certified frontier; eta=1 -> Hall (S-2)/6;
       M=0 boundary -> S = 4/eta - 2 (postselection ceiling).
"""
import numpy as np
from scipy.optimize import linprog

vals = [1, -1, 0]
astrat = [(v1, v2) for v1 in vals for v2 in vals]
NS = 81
strat = [(sa, sb) for sa in astrat for sb in astrat]
pair_idx = [(0,0),(0,1),(1,0),(1,1)]
chsh_sign = [+1, -1, +1, +1]
cmb = [(i,j) for i in range(4) for j in range(i+1,4)]

def md_floor(S, eta):
    eta2 = eta*eta
    nq = 4*NS; nd = len(cmb)*NS; nvar = nq + nd + 1
    c = np.zeros(nvar); c[-1] = 1.0
    A_eq, b_eq = [], []
    chsh_row = np.zeros(nvar)
    W_list, dA_list, dB_list = [], [], []
    for k,(ia,ib) in enumerate(pair_idx):
        A = np.array([sa[ia] for sa,sb in strat], float)
        B = np.array([sb[ib] for sa,sb in strat], float)
        aA, aB = np.abs(A), np.abs(B)
        W = aA*aB
        W_list.append(W); dA_list.append(aA); dB_list.append(aB)
        row = np.zeros(nvar); row[k*NS:(k+1)*NS] = 1.0
        A_eq.append(row); b_eq.append(1.0)
        row = np.zeros(nvar); row[k*NS:(k+1)*NS] = W
        A_eq.append(row); b_eq.append(eta2)
        row = np.zeros(nvar); row[k*NS:(k+1)*NS] = aA
        A_eq.append(row); b_eq.append(eta)
        row = np.zeros(nvar); row[k*NS:(k+1)*NS] = aB
        A_eq.append(row); b_eq.append(eta)
        for vec in (A*aB, aA*B, A, B):
            row = np.zeros(nvar); row[k*NS:(k+1)*NS] = vec
            A_eq.append(row); b_eq.append(0.0)
        chsh_row[k*NS:(k+1)*NS] = chsh_sign[k]*(A*B)
    A_eq.append(chsh_row); b_eq.append(S*eta2)     # CHSH on coincidences = S
    A_ub, b_ub = [], []
    for p,(i,j) in enumerate(cmb):
        for s in range(NS):
            r = np.zeros(nvar); r[i*NS+s]=1; r[j*NS+s]=-1; r[nq+p*NS+s]=-1
            A_ub.append(r); b_ub.append(0.0)
            r = np.zeros(nvar); r[i*NS+s]=-1; r[j*NS+s]=1; r[nq+p*NS+s]=-1
            A_ub.append(r); b_ub.append(0.0)
    for p in range(len(cmb)):
        r = np.zeros(nvar); r[nq+p*NS:nq+(p+1)*NS]=1; r[-1]=-1
        A_ub.append(r); b_ub.append(0.0)
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(0,None)]*nvar, method="highs")
    return res.x[-1]/2 if res.status == 0 else np.nan

def conj(S, eta): return max(0.0, eta*((S+2)*eta - 4)/6)

print(f"{'S':>7} {'eta':>6} {'M_LP':>13} {'conjecture':>13} {'diff':>10}")
worst = 0.0
Ss = [2.1, 2.3, 2.5, 2.7, 2*np.sqrt(2), 2.9, 3.1, 3.4]
etas = [0.85, 0.90, 0.95, 1.00]
rows = []
for S in Ss:
    for eta in etas:
        m = md_floor(S, eta)
        cj = conj(S, eta)
        d = m - cj
        if not np.isnan(m):
            worst = max(worst, abs(d))
        rows.append((S, eta, m, cj))
        print(f"{S:7.4f} {eta:6.2f} {m:13.9f} {cj:13.9f} {d:10.2e}")
print(f"\nmax |LP - conjecture| over the surface: {worst:.2e}")
np.savez("/home/claude/price_surface.npz", rows=np.array(rows))
