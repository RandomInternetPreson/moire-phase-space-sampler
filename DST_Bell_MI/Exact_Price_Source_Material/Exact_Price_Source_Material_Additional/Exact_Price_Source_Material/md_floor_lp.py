#!/usr/bin/env python3
"""
md_floor_lp.py - Minimal measurement dependence to reproduce the four CHSH
singlet correlators with the deterministic sign model at PERFECT detection.

Variables: four conditional hidden-variable distributions q_ab(lambda), one
per setting pair. Constraints: correct correlators E_ab = -cos 2(a-b), zero
marginals (faithfulness), normalization. Objective: minimize the maximum
variational distance between any two of the conditionals (Hall-style MD
measure; both the integral |.| and half-integral conventions reported).
"""
import numpy as np
from scipy.optimize import linprog

NB = 180
lam = (np.arange(NB) + 0.5) * np.pi / NB
_deg = np.pi/180
settings = [(0.0, 22.5), (0.0, 67.5), (45.0, 22.5), (45.0, 67.5)]
targets = [-np.cos(2*(a-b)*_deg) for a, b in settings]

def A(a): return np.sign(np.cos(2*(lam - a*_deg)) + 1e-300)
def B(b): return -np.sign(np.cos(2*(lam - b*_deg)) + 1e-300)

nq = 4*NB
pairs = [(i, j) for i in range(4) for j in range(i+1, 4)]
nd = len(pairs)*NB
nvar = nq + nd + 1   # q's, abs-aux d's, t

c = np.zeros(nvar); c[-1] = 1.0

A_eq, b_eq = [], []
for k, (a, b) in enumerate(settings):
    row = np.zeros(nvar); row[k*NB:(k+1)*NB] = 1.0
    A_eq.append(row); b_eq.append(1.0)                      # normalization
    row = np.zeros(nvar); row[k*NB:(k+1)*NB] = A(a)*B(b)
    A_eq.append(row); b_eq.append(targets[k])               # correlator
    row = np.zeros(nvar); row[k*NB:(k+1)*NB] = A(a)
    A_eq.append(row); b_eq.append(0.0)                      # Alice marginal
    row = np.zeros(nvar); row[k*NB:(k+1)*NB] = B(b)
    A_eq.append(row); b_eq.append(0.0)                      # Bob marginal

A_ub, b_ub = [], []
for p, (i, j) in enumerate(pairs):
    for l in range(NB):
        r = np.zeros(nvar)
        r[i*NB+l] = 1.0; r[j*NB+l] = -1.0; r[nq+p*NB+l] = -1.0
        A_ub.append(r); b_ub.append(0.0)
        r = np.zeros(nvar)
        r[i*NB+l] = -1.0; r[j*NB+l] = 1.0; r[nq+p*NB+l] = -1.0
        A_ub.append(r); b_ub.append(0.0)
for p in range(len(pairs)):
    r = np.zeros(nvar)
    r[nq+p*NB:nq+(p+1)*NB] = 1.0; r[-1] = -1.0
    A_ub.append(r); b_ub.append(0.0)

res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
              A_eq=np.array(A_eq), b_eq=np.array(b_eq),
              bounds=[(0, None)]*nq + [(0, None)]*nd + [(0, None)],
              method="highs")
print("status:", res.message)
t = res.x[-1]
print(f"minimal MD (integral |q_i - q_j| convention):      {t:.6f}")
print(f"minimal MD (variational/half convention):          {t/2:.6f}")
print(f"Hall 2010 reference value (sqrt(2)-1)/3:           {(np.sqrt(2)-1)/3:.6f}")
print(f"as a fraction of maximal MD (=2 in integral conv): {t/2*100:.2f}%")
