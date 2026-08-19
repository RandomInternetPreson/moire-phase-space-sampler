#!/usr/bin/env python3
"""
certify_ladder.py - Gate G1: exact certification of the price surface at
rational (S, eta) points, in plain Fraction arithmetic.
Pipeline: float LP proposes -> Fraction snapping -> EXACT verification of
primal (all equalities, nonnegativity, objective) and dual (feasibility,
matching objective). Every accepted value is a theorem.
"""
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as Fr

vals = [1,-1,0]
astrat = [(v1,v2) for v1 in vals for v2 in vals]
NS = 81
strat = [(sa,sb) for sa in astrat for sb in astrat]
pair_idx = [(0,0),(0,1),(1,0),(1,1)]
chsh_sign = [1,-1,1,1]
cmb = [(i,j) for i in range(4) for j in range(i+1,4)]

def build_rows(S, eta):
    """Exact (Fraction) equality rows and float copies."""
    rows_F, rhs_F = [], []
    chsh_F = [Fr(0)]*(4*NS)
    for k,(ia,ib) in enumerate(pair_idx):
        A = [sa[ia] for sa,sb in strat]; B = [sb[ib] for sa,sb in strat]
        aA = [abs(v) for v in A]; aB = [abs(v) for v in B]
        base = k*NS
        def row(vec, rhs):
            r = [Fr(0)]*(4*NS)
            for s in range(NS): r[base+s] = Fr(vec[s])
            rows_F.append(r); rhs_F.append(rhs)
        row([1]*NS, Fr(1))
        row([aA[s]*aB[s] for s in range(NS)], eta*eta)
        row(aA, eta); row(aB, eta)
        row([A[s]*aB[s] for s in range(NS)], Fr(0))
        row([aA[s]*B[s] for s in range(NS)], Fr(0))
        row(A, Fr(0)); row(B, Fr(0))
        for s in range(NS):
            chsh_F[base+s] = Fr(chsh_sign[k]*A[s]*B[s])
    rows_F.append(chsh_F); rhs_F.append(S*eta*eta)
    return rows_F, rhs_F

def float_lp(S, eta):
    nq = 4*NS; nd = len(cmb)*NS; nvar = nq+nd+1
    rows_F, rhs_F = build_rows(S, eta)
    A_eq = np.zeros((len(rows_F), nvar)); b_eq = np.zeros(len(rows_F))
    for r,(row,rhs) in enumerate(zip(rows_F, rhs_F)):
        A_eq[r,:4*NS] = [float(v) for v in row]; b_eq[r] = float(rhs)
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
    res = linprog(np.eye(nvar)[-1], A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=A_eq, b_eq=b_eq, bounds=[(0,None)]*nvar, method="highs")
    return res, rows_F, rhs_F

def snapF(v, md=10**9):
    return Fr(v).limit_denominator(md)

def certify(Sf, ef, S, eta):
    target = max(Fr(0), eta*((S+2)*eta-4)/6)
    res, _, _ = float_lp(Sf, ef)
    if res.status != 0: return "LP FAIL", target
    rows_F, rhs_F = build_rows(S, eta)   # EXACT rows from the Fraction inputs
    nq = 4*NS
    q = [snapF(v) if v > 1e-9 else Fr(0) for v in res.x[:nq]]
    # primal verification
    for row, rhs in zip(rows_F, rhs_F):
        if sum(c*qv for c,qv in zip(row,q) if qv or c) != rhs:
            return "PRIMAL EQ FAIL", target
    if any(qv < 0 for qv in q): return "PRIMAL NEG FAIL", target
    tvs = [sum(abs(q[i*NS+s]-q[j*NS+s]) for s in range(NS))/2 for i,j in cmb]
    M_ex = max(tvs)
    if M_ex != target: return f"PRIMAL VALUE {M_ex} != {target}", target
    # dual verification
    y = [snapF(v) if abs(v) > 1e-9 else Fr(0) for v in res.eqlin.marginals]
    z = [snapF(v) if abs(v) > 1e-9 else Fr(0) for v in res.ineqlin.marginals]
    if any(zz > 0 for zz in z): return "DUAL SIGN FAIL", target
    obj = sum(rhs*yy for rhs,yy in zip(rhs_F,y))
    if obj != 2*target: return f"DUAL OBJ {obj} != {2*target}", target
    nvar = nq + len(cmb)*NS + 1
    red = [Fr(0)]*nvar
    for row, yy in zip(rows_F, y):
        if yy == 0: continue
        for s in range(4*NS):
            if row[s]: red[s] += row[s]*yy
    idx = 0
    for p,(i,j) in enumerate(cmb):
        for s in range(NS):
            for (si,sj) in ((1,-1),(-1,1)):
                zz = z[idx]
                if zz:
                    red[i*NS+s] += si*zz; red[j*NS+s] += sj*zz
                    red[nq+p*NS+s] += -zz
                idx += 1
    for p in range(len(cmb)):
        zz = z[idx]
        if zz:
            for s in range(NS): red[nq+p*NS+s] += zz
            red[nvar-1] += -zz
        idx += 1
    c = [Fr(0)]*nvar; c[-1] = Fr(1)
    if any(red[i] > c[i] for i in range(nvar)): return "DUAL FEAS FAIL", target
    return "CERTIFIED", target

points = [(Fr(5,2), Fr(1)), (Fr(5,2), Fr(19,20)), (Fr(27,10), Fr(9,10)),
          (Fr(29,10), Fr(19,20)), (Fr(17,5), Fr(17,20)), (Fr(21,10), Fr(1))]
print("Gate G1: rational certification ladder for M(S,eta) = eta((S+2)eta-4)/6")
passed = 0
for S, eta in points:
    verdict, target = certify(float(S), float(eta), S, eta)
    tag = f"M = {target}"
    print(f"  (S={S}, eta={eta}): {verdict:12s} {tag}")
    if verdict == "CERTIFIED": passed += 1
print(f"\n{passed}/{len(points)} points certified exactly in Q.")
