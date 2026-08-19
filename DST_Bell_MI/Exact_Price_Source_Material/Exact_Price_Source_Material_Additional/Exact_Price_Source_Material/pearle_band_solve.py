#!/usr/bin/env python3
"""
pearle_band_solve.py — Promote the edge-band width to a hidden variable.

Model: each pair carries (lambda, u), both set at the source, independent of
settings (no measurement dependence). lambda ~ uniform [0, pi). u ~ rho(u),
the unknown "ambiguity width" distribution. Station rule stays deterministic
and local: outcome = sign(cos 2(lambda - a)); click iff |cos 2(lambda - a)| >= u.

For fixed u we can compute exactly (lambda-grid integration):
  Nu(D)  = pair coincidence probability at relative angle D
  Eu(D)  = coincidence correlator at D
The mixture gives  E(D) = sum_k p_k Nu_k(D) Eu_k(D) / sum_k p_k Nu_k(D).

Demanding E(D) = -cos(2D) for all D is LINEAR in p:
  sum_k p_k * c_k(D) = 0,   c_k(D) = Nu_k(D) * (Eu_k(D) + cos 2D)

Solve by NNLS over the simplex, then verify with independent Monte Carlo.
Consistency checks the solution MUST pass if it reproduces full QM statistics:
  - S at Bell angles = 2*sqrt(2) with no Tsirelson overshoot possible
  - singles efficiency <= Garg-Mermin bound 0.8284 (theorem, non-negotiable)
"""

import numpy as np
from scipy.optimize import nnls

_deg = np.pi / 180.0
S_QM = 2.0 * np.sqrt(2.0)
ETA_GM = 2.0 * (np.sqrt(2.0) - 1.0)

NLAM = 5760
LAM = np.linspace(0.0, np.pi, NLAM, endpoint=False)
THETAS = np.linspace(0.0, 180.0, 181)          # constraint grid for Delta
K = 80
US = np.linspace(0.0, 0.9875, K)               # edge-band width grid

def pair_stats(u, D):
    """Exact E and N at band u, relative angle D (deg), a=0, b=D."""
    ga = np.cos(2.0 * LAM)
    gb = np.cos(2.0 * (LAM + np.pi / 2.0 - D * _deg))
    det = (np.abs(ga) >= u) & (np.abs(gb) >= u)
    n = det.mean()
    if n == 0:
        return 0.0, 0.0
    oa = np.where(ga >= 0, 1.0, -1.0)
    ob = np.where(gb >= 0, 1.0, -1.0)
    return float(np.mean(oa[det] * ob[det]) * NLAM / det.sum() * n) / n, n

print("building constraint matrix ...")
E_tab = np.zeros((K, len(THETAS)))
N_tab = np.zeros((K, len(THETAS)))
for k, u in enumerate(US):
    for j, D in enumerate(THETAS):
        e, n = pair_stats(u, D)
        E_tab[k, j], N_tab[k, j] = e, n

C = (N_tab * (E_tab + np.cos(2.0 * THETAS * _deg))).T   # (n_theta, K)

# NNLS: minimize ||C p||  subject to p >= 0, sum p = 1 (soft constraint row)
gamma = 5.0
A = np.vstack([C, gamma * np.ones((1, K))])
b = np.concatenate([np.zeros(len(THETAS)), [gamma]])
p, res = nnls(A, b)
p = p / p.sum()

# Evaluate the mixture exactly
num = (p[:, None] * N_tab * E_tab).sum(axis=0)
den = (p[:, None] * N_tab).sum(axis=0)
E_mix = num / den
target = -np.cos(2.0 * THETAS * _deg)
dev = np.abs(E_mix - target)
print(f"NNLS residual: {res:.3e}")
print(f"max |E(D) + cos 2D| = {dev.max():.6f}")
print(f"mean deviation      = {dev.mean():.6f}")

# CHSH at Bell angles from the exact mixture
def E_at(a, b):
    ga = np.cos(2.0 * (LAM - a * _deg))
    gb = np.cos(2.0 * (LAM + np.pi / 2.0 - b * _deg))
    oa = np.where(ga >= 0, 1.0, -1.0)
    ob = np.where(gb >= 0, 1.0, -1.0)
    num = den = 0.0
    for k, u in enumerate(US):
        if p[k] < 1e-12:
            continue
        det = (np.abs(ga) >= u) & (np.abs(gb) >= u)
        num += p[k] * np.mean(oa * ob * det)
        den += p[k] * det.mean()
    return num / den

S = abs(E_at(0, 22.5) - E_at(0, 67.5) + E_at(45, 22.5) + E_at(45, 67.5))
eta = float(sum(p[k] * (1.0 - (2.0 / np.pi) * np.arcsin(min(US[k], 1.0)))
               for k in range(K)))
pair_eff = float((p[:, None] * N_tab).sum(axis=0).mean())
print(f"\nS at Bell angles     = {S:.6f}   (QM: {S_QM:.6f})")
print(f"singles efficiency   = {eta:.4f}   (Garg-Mermin bound: {ETA_GM:.4f}, "
      f"{'RESPECTED' if eta <= ETA_GM + 1e-6 else 'VIOLATED -> bug'})")
print(f"mean pair efficiency = {pair_eff:.4f}")

print("\nrho(u) support (p_k > 0.001):")
for k in range(K):
    if p[k] > 0.001:
        print(f"  u = {US[k]:.4f}   weight = {p[k]:.4f}")

# Independent Monte Carlo cross-check
rng = np.random.default_rng(7)
n = 4_000_000
lam = rng.uniform(0, np.pi, n)
u_s = rng.choice(US, size=n, p=p)
def mc_E(a, b):
    ga = np.cos(2 * (lam - a * _deg)); gb = np.cos(2 * (lam + np.pi/2 - b * _deg))
    det = (np.abs(ga) >= u_s) & (np.abs(gb) >= u_s)
    return float(np.mean(np.sign(ga)[det] * np.sign(gb)[det]))
S_mc = abs(mc_E(0,22.5) - mc_E(0,67.5) + mc_E(45,22.5) + mc_E(45,67.5))
print(f"\nMonte Carlo cross-check (4M pairs): S = {S_mc:.4f}, "
      f"E(22.5) = {mc_E(0,22.5):.4f} (QM: {-np.cos(45*_deg):.4f}), "
      f"E(30) = {mc_E(0,30):.4f} (QM: {-np.cos(60*_deg):.4f})")

# Plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.2))
ax1.bar(US, p, width=0.011, color="#7F77DD")
ax1.set_xlabel("edge-band width u"); ax1.set_ylabel("weight rho(u)")
ax1.set_title("Solved hidden-width distribution")
ax2.plot(THETAS, E_mix, color="#1D9E75", lw=2, label="mixture model")
ax2.plot(THETAS, target, color="#D85A30", ls="--", lw=1.5, label="QM: -cos 2D")
ax2.set_xlabel("relative angle (deg)"); ax2.set_ylabel("E coincidence")
ax2.set_title(f"Coincidence curve (max dev {dev.max():.4f})"); ax2.legend(fontsize=8)
ax3.plot(THETAS, E_mix - target, color="#1D9E75", lw=1.5)
ax3.axhline(0, color="#888", ls=":", lw=1)
ax3.set_xlabel("relative angle (deg)"); ax3.set_ylabel("residual")
ax3.set_title("Residual vs QM")
fig.tight_layout()
fig.savefig("/home/claude/pearle_band_solution.png", dpi=150)
np.savez("/home/claude/pearle_band_solution.npz", us=US, p=p)
print("\nsaved plot and solution arrays")
