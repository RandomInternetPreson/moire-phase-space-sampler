#!/usr/bin/env python3
"""
exact_frontier_lp.py - The model-independent MD/detection frontier for CHSH.

Strategy space: Alice assigns each of her 2 settings a value in {+1,-1,0}
(0 = no click) -> 9 strategies; Bob likewise -> 81 joint deterministic
strategies. ANY local hidden-variable model (with or without MD) reduces to
four distributions q_k over these 81 strategies, one per setting pair.
MD = max pairwise variational distance between the q_k.

Faithfulness constraints (what an experimenter sees must match QM with
symmetric detector efficiency eta):
  correlators on coincidences = -cos 2(a-b)     (all four)
  coincidence rate            = eta^2
  singles rates (each party)  = eta
  marginals (coincidence-conditioned and detected-singles) = 0
Objective: minimize MD. Sweep eta.

Registered conjecture: M(eta) = (1+sqrt2)/3 - 2/(3 eta)
  -> M(1) = (sqrt2-1)/3 (Hall), zero crossing at eta = 2(sqrt2-1) (Garg-Mermin)
"""
import numpy as np
from scipy.optimize import linprog

SQ2 = np.sqrt(2.0)
vals = [1, -1, 0]
astrat = [(v1, v2) for v1 in vals for v2 in vals]   # 9
NS = 81
strat = [(sa, sb) for sa in astrat for sb in astrat]

# setting pairs: (Alice idx, Bob idx), E targets with signs for CHSH (+,-,+,+)
pair_idx = [(0,0),(0,1),(1,0),(1,1)]
E_t = [-SQ2/2, +SQ2/2, -SQ2/2, -SQ2/2]

def build(eta):
    eta2 = eta*eta
    nq = 4*NS
    cmb = [(i,j) for i in range(4) for j in range(i+1,4)]
    nd = len(cmb)*NS
    nvar = nq + nd + 1
    c = np.zeros(nvar); c[-1] = 1.0
    A_eq, b_eq = [], []
    for k,(ia,ib) in enumerate(pair_idx):
        A = np.array([sa[ia] for sa,sb in strat], float)
        B = np.array([sb[ib] for sa,sb in strat], float)
        aA, aB = np.abs(A), np.abs(B)
        rows = [
            (np.ones(NS), 1.0),
            (A*B - E_t[k]*aA*aB, 0.0),
            (aA*aB, eta2),
            (aA, eta), (aB, eta),
            (A*aB, 0.0), (aA*B, 0.0),
            (A, 0.0), (B, 0.0),
        ]
        for vec, rhs in rows:
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
    return c, np.array(A_ub), np.array(b_ub), np.array(A_eq), np.array(b_eq), nvar

def md_floor(eta):
    c, Aub, bub, Aeq, beq, nvar = build(eta)
    res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq,
                  bounds=[(0,None)]*nvar, method="highs")
    return res.x[-1]/2 if res.status == 0 else np.nan

def conj(eta): return max(0.0, (1+SQ2)/3 - 2/(3*eta))

print(f"{'eta':>8} {'MD_LP':>12} {'conjecture':>12} {'diff':>11}")
etas = np.concatenate([np.linspace(0.80, 1.00, 21), [0.828427, 0.9, 0.95, 8/9, 19/20]])
etas = np.unique(np.round(etas, 6))
rows = []
for eta in etas:
    md = md_floor(eta)
    cj = conj(eta)
    rows.append((eta, md, cj))
    print(f"{eta:8.4f} {md:12.8f} {cj:12.8f} {md-cj:11.2e}")

# zero crossing by bisection
lo, hi = 0.80, 0.90
for _ in range(40):
    mid = 0.5*(lo+hi)
    if md_floor(mid) > 1e-9: hi = mid
    else: lo = mid
eta_star = 0.5*(lo+hi)
print(f"\nzero-MD crossing: eta* = {eta_star:.8f}")
print(f"2(sqrt2 - 1)      = {2*(SQ2-1):.8f}   diff = {eta_star-2*(SQ2-1):.2e}")
print(f"Hall value at eta=1: LP = {md_floor(1.0):.8f}, (sqrt2-1)/3 = {(SQ2-1)/3:.8f}")

rows = np.array(rows)
np.savez("exact_frontier.npz", rows=rows, eta_star=eta_star)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7.5,5))
ee = np.linspace(0.80, 1.0, 200)
ax.plot(ee, [conj(x) for x in ee], color="#D85A30", ls="--", lw=1.5,
        label="conjecture (1+√2)/3 − 2/(3η)")
ax.plot(rows[:,0], rows[:,1], "o", color="#1D9E75", ms=5, label="exact LP (all 81 strategies)")
ax.axvline(2*(SQ2-1), color="#888", ls=":", lw=1, label="Garg-Mermin 2(√2−1)")
ax.axhline((SQ2-1)/3, color="#888", ls=":", lw=1)
ax.set_xlabel("detector efficiency η (faithful rates: singles η, coinc η²)")
ax.set_ylabel("minimal measurement dependence (variational)")
ax.set_title("Model-independent frontier: MD floor vs detection efficiency")
ax.legend(fontsize=9); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig("/home/claude/exact_frontier.png", dpi=150)
print("saved exact frontier plot")
