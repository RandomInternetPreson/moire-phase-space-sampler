#!/usr/bin/env python3
"""
certify_frontier.py - Exact certification of M(9/10) = (27*sqrt(2) - 33)/100.

Method (the exact_certificates.py pattern):
  1. Float LP (HiGHS) at eta = 9/10 finds the optimal solution structure.
  2. Every nonzero primal weight and dual multiplier is snapped to Q(sqrt2)
     via integer-relation detection (PSLQ over the basis {1, sqrt2}).
  3. The snapped candidates are verified in EXACT symbolic arithmetic (sympy,
     sqrt(2) kept symbolic):
       primal: q >= 0, all 36 faithfulness equalities hold exactly, and the
               max pairwise variational distance equals the target exactly
               -> an explicit local model achieving M (upper bound).
       dual:   dual feasibility A_eq^T y + A_ub^T z <= c with z <= 0, and
               objective b^T y equals the target exactly
               -> proof that NO local model does better (lower bound).
  4. Upper bound + lower bound = certified exact value.
Float arithmetic only ever PROPOSES; every accepted statement is verified
symbolically. A failed snap or failed verification is reported, not hidden.
"""
import numpy as np
from scipy.optimize import linprog
import mpmath
import sympy as sp

SQ2f = np.sqrt(2.0)
sq2 = sp.sqrt(2)
ETA = sp.Rational(9, 10)
TARGET = (27*sq2 - 33) / 100

vals = [1, -1, 0]
astrat = [(v1, v2) for v1 in vals for v2 in vals]
NS = 81
strat = [(sa, sb) for sa in astrat for sb in astrat]
pair_idx = [(0,0),(0,1),(1,0),(1,1)]
E_sign = [-1, +1, -1, -1]          # E_k = E_sign[k] * sqrt2/2
cmb = [(i,j) for i in range(4) for j in range(i+1,4)]

# ---------- exact constraint data ----------
def exact_rows():
    A_eq, b_eq = [], []
    for k,(ia,ib) in enumerate(pair_idx):
        A = [sa[ia] for sa,sb in strat]
        B = [sb[ib] for sa,sb in strat]
        aA = [abs(x) for x in A]; aB = [abs(x) for x in B]
        Ek = E_sign[k]*sq2/2
        rows = [
            ([sp.Integer(1)]*NS, sp.Integer(1)),
            ([A[s]*B[s] - Ek*aA[s]*aB[s] for s in range(NS)], sp.Integer(0)),
            ([sp.Integer(aA[s]*aB[s]) for s in range(NS)], ETA*ETA),
            ([sp.Integer(aA[s]) for s in range(NS)], ETA),
            ([sp.Integer(aB[s]) for s in range(NS)], ETA),
            ([sp.Integer(A[s]*aB[s]) for s in range(NS)], sp.Integer(0)),
            ([sp.Integer(aA[s]*B[s]) for s in range(NS)], sp.Integer(0)),
            ([sp.Integer(A[s]) for s in range(NS)], sp.Integer(0)),
            ([sp.Integer(B[s]) for s in range(NS)], sp.Integer(0)),
        ]
        for vec, rhs in rows:
            A_eq.append((k, vec)); b_eq.append(rhs)
    return A_eq, b_eq

# ---------- float LP (structure finder) ----------
def float_lp(eta):
    eta2 = eta*eta
    nq = 4*NS; nd = len(cmb)*NS; nvar = nq + nd + 1
    c = np.zeros(nvar); c[-1] = 1.0
    A_eq, b_eq = [], []
    for k,(ia,ib) in enumerate(pair_idx):
        A = np.array([sa[ia] for sa,sb in strat], float)
        B = np.array([sb[ib] for sa,sb in strat], float)
        aA, aB = np.abs(A), np.abs(B)
        Ek = E_sign[k]*SQ2f/2
        for vec, rhs in [(np.ones(NS),1.0),(A*B - Ek*aA*aB,0.0),(aA*aB,eta2),
                         (aA,eta),(aB,eta),(A*aB,0.0),(aA*B,0.0),(A,0.0),(B,0.0)]:
            r = np.zeros(nvar); r[k*NS:(k+1)*NS] = vec
            A_eq.append(r); b_eq.append(rhs)
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
    return res, np.array(A_eq), np.array(A_ub), np.array(b_eq)

# ---------- Q(sqrt2) snapping via PSLQ ----------
mpmath.mp.dps = 30
def snap(v, tol=1e-7, maxden=100000):
    if abs(v) < 5e-9:
        return sp.Integer(0)
    rel = mpmath.pslq([mpmath.mpf(1), mpmath.sqrt(2), mpmath.mpf(v)],
                      tol=mpmath.mpf(tol), maxcoeff=maxden)
    if rel is None or rel[2] == 0:
        return None
    cand = sp.Rational(-rel[0], rel[2]) + sp.Rational(-rel[1], rel[2])*sq2
    if abs(float(cand) - v) > 1e-6:
        return None
    return cand

print("=== float LP at eta = 9/10 ===")
res, Aeq_f, Aub_f, beq_f = float_lp(0.9)
t_float = res.x[-1]/2
print(f"status: {res.message}")
print(f"float MD = {t_float:.12f}   target (27*sqrt2-33)/100 = {float(TARGET):.12f}")

q_float = res.x[:4*NS]
support = [(idx, q_float[idx]) for idx in range(4*NS) if q_float[idx] > 1e-8]
print(f"primal support size: {len(support)} nonzero weights across 4 distributions")

print("\n=== snapping primal weights to Q(sqrt2) ===")
q_exact = [sp.Integer(0)]*(4*NS)
fails = 0
for idx, w in support:
    s = snap(w)
    if s is None:
        fails += 1
        print(f"  SNAP FAIL at var {idx}: {w:.12f}")
    else:
        q_exact[idx] = s
print(f"snapped {len(support)-fails}/{len(support)} weights")

if fails == 0:
    print("\n=== exact primal verification (sympy) ===")
    A_eq_sym, b_eq_sym = exact_rows()
    ok = True
    for r,((k, vec), rhs) in enumerate(zip(A_eq_sym, b_eq_sym)):
        lhs = sum(vec[s]*q_exact[k*NS+s] for s in range(NS))
        if sp.simplify(lhs - rhs) != 0:
            ok = False
            print(f"  EQUALITY {r} FAILS: lhs - rhs = {sp.simplify(lhs-rhs)}")
    neg = [i for i in range(4*NS) if q_exact[i].is_negative]
    if neg:
        ok = False; print(f"  NEGATIVE WEIGHTS at {neg}")
    # exact pairwise variational distances
    tvs = []
    for (i,j) in cmb:
        tv = sum(sp.Abs(q_exact[i*NS+s] - q_exact[j*NS+s]) for s in range(NS))/2
        tvs.append(sp.simplify(tv))
    M_exact = sp.simplify(sp.Max(*tvs))
    print(f"  all 36 equalities hold exactly: {ok}")
    print(f"  exact pairwise TVs: {[sp.nsimplify(t) for t in tvs]}")
    print(f"  exact M (max TV)   = {sp.nsimplify(M_exact)}")
    print(f"  target             = {sp.nsimplify(TARGET)}")
    print(f"  M == target: {sp.simplify(M_exact - TARGET) == 0}")
    print(f"  -> PRIMAL CERTIFICATE {'VALID' if ok and sp.simplify(M_exact-TARGET)==0 else 'INVALID'}: "
          f"an explicit exact local model achieves M(9/10) = (27*sqrt2-33)/100")

print("\n=== dual certificate (lower bound) ===")
y_f = res.eqlin.marginals          # 36 duals for equalities (free sign)
z_f = res.ineqlin.marginals        # duals for <= constraints (<= 0)
y_nz = int(np.sum(np.abs(y_f) > 5e-9)); z_nz = int(np.sum(np.abs(z_f) > 5e-9))
print(f"nonzero duals: {y_nz}/36 equality, {z_nz}/{len(z_f)} inequality")
y_exact, z_exact, dfails = [], [], 0
for v in y_f:
    s = snap(v)
    if s is None: dfails += 1; s = sp.Integer(0)
    y_exact.append(s)
for v in z_f:
    s = snap(v)
    if s is None: dfails += 1; s = sp.Integer(0)
    z_exact.append(s)
print(f"dual snap failures: {dfails}")

if dfails == 0:
    A_eq_sym, b_eq_sym = exact_rows()
    nq = 4*NS; nvar = nq + len(cmb)*NS + 1
    # exact objective: b^T y  (b_ub = 0 contributes nothing)
    obj = sp.simplify(sum(b_eq_sym[r]*y_exact[r] for r in range(36)))
    # dual objective in linprog's convention equals the primal optimum of t;
    # our M is t/2, so compare obj to 2*TARGET
    print(f"  exact dual objective b^T y = {sp.nsimplify(obj)}  (should equal 2*target = {sp.nsimplify(2*TARGET)})")
    # dual feasibility: A_eq^T y + A_ub^T z <= c, z <= 0
    z_ok = all(not z.is_positive for z in z_exact)
    # build A_ub^T z exactly using its sparse structure
    red = [sp.Integer(0)]*nvar
    for r,((k, vec), _) in enumerate(zip(A_eq_sym, b_eq_sym)):
        if y_exact[r] == 0: continue
        for s in range(NS):
            if vec[s] != 0:
                red[k*NS+s] += vec[s]*y_exact[r]
    row = 0
    for p,(i,j) in enumerate(cmb):
        for s in range(NS):
            for (sgn_i, sgn_j) in ((1,-1),(-1,1)):
                zv = z_exact[row]
                if zv != 0:
                    red[i*NS+s] += sgn_i*zv
                    red[j*NS+s] += sgn_j*zv
                    red[nq+p*NS+s] += -zv
                row += 1
    for p in range(len(cmb)):
        zv = z_exact[row]
        if zv != 0:
            for s in range(NS):
                red[nq+p*NS+s] += zv
            red[nvar-1] += -zv
        row += 1
    c_exact = [sp.Integer(0)]*nvar; c_exact[-1] = sp.Integer(1)
    viol = [i for i in range(nvar) if sp.simplify(red[i] - c_exact[i]).is_positive]
    print(f"  z <= 0 exactly: {z_ok}")
    print(f"  dual feasibility A^T(y,z) <= c: {'holds exactly' if not viol else f'VIOLATED at {len(viol)} components'}")
    valid = z_ok and not viol and sp.simplify(obj - 2*TARGET) == 0
    print(f"  -> DUAL CERTIFICATE {'VALID' if valid else 'INVALID'}: no local model achieves less")
    if valid:
        print("\n*** CERTIFIED: M(9/10) = (27*sqrt(2) - 33)/100 exactly. ***")
else:
    print("  dual snapping incomplete -> lower bound not certified this run;")
    print("  recourse: re-solve with exact rational simplex on the identified basis (workstation task)")
