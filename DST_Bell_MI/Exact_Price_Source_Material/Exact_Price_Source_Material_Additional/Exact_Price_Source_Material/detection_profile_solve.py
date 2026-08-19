#!/usr/bin/env python3
"""
detection_profile_solve.py - Solve for the smooth detection profile D(m).

Model: pair carries shared phase lambda ~ U[0,pi) plus INDEPENDENT per-photon
uniforms r_A, r_B (all set at source, setting-independent -> local, no MD).
Station: m = |cos 2(lambda - a)|, outcome = sign(cos 2(lambda - a)),
click iff r <= D(m). Solve for D: [0,1] -> [0,1] such that the coincidence
correlation equals -cos(2 Delta) exactly.
"""
import numpy as np
from scipy.optimize import least_squares

_deg = np.pi/180; S_QM = 2*np.sqrt(2); ETA_GM = 2*(np.sqrt(2)-1)
NLAM = 4320
LAM = np.linspace(0, np.pi, NLAM, endpoint=False)
TH = np.linspace(0.0, 180.0, 91)
M = 33
MG = np.linspace(0.0, 1.0, M)

cA = np.cos(2*LAM); mA = np.abs(cA); sA = np.sign(cA + 1e-300)

def D_of(dvals, m):
    return np.interp(m, MG, dvals)

def E_D(dvals, a, b):
    ca = np.cos(2*(LAM - a*_deg)); cb = np.cos(2*(LAM + np.pi/2 - b*_deg))
    wA = D_of(dvals, np.abs(ca)); wB = D_of(dvals, np.abs(cb))
    ww = wA*wB
    den = ww.sum()
    if den < 1e-12: return 0.0
    return float((np.sign(ca)*np.sign(cb)*ww).sum()/den)

def residuals(dvals):
    return np.array([E_D(dvals, 0.0, t) + np.cos(2*t*_deg) for t in TH])

# initialize between the two brackets: D = m^0.5
x0 = MG**0.5
sol = least_squares(residuals, x0, bounds=(0.0, 1.0), xtol=1e-14, ftol=1e-14,
                    max_nfev=6000)
D = sol.x
r = residuals(D)
print(f"converged: {sol.status}, max |E + cos2D| = {np.max(np.abs(r)):.2e}, "
      f"rms = {np.sqrt(np.mean(r**2)):.2e}")

S = abs(E_D(D,0,22.5) - E_D(D,0,67.5) + E_D(D,45,22.5) + E_D(D,45,67.5))
eta = float(np.mean(D_of(D, mA)))
pair = float(np.mean(D_of(D,mA)*D_of(D,np.abs(np.cos(2*(LAM+np.pi/2-22.5*_deg))))))
print(f"S at Bell angles = {S:.6f} (QM {S_QM:.6f})")
print(f"singles efficiency eta = {eta:.4f}  (Garg-Mermin bound {ETA_GM:.4f}: "
      f"{'respected' if eta <= ETA_GM + 1e-4 else 'VIOLATED - investigate'})")
print(f"pair efficiency at 22.5 deg = {pair:.4f}")
print("\nD(m) profile:")
for i in range(0, M, 4):
    print(f"  m = {MG[i]:.3f}   D = {D[i]:.4f}")
print(f"  m = 1.000   D = {D[-1]:.4f}")

# fit check: is D close to a power law m^p?
mask = (MG > 0.02) & (D > 1e-6)
p_fit = np.polyfit(np.log(MG[mask]), np.log(D[mask]/max(D[-1],1e-12)), 1)[0]
print(f"\npower-law exponent fit (D ~ m^p): p = {p_fit:.4f}")

# Monte Carlo verification, fully independent
rng = np.random.default_rng(42); n = 6_000_000
lam = rng.uniform(0, np.pi, n)
rA = rng.uniform(0,1,n); rB = rng.uniform(0,1,n)
def mc(a,b):
    ca = np.cos(2*(lam-a*_deg)); cb = np.cos(2*(lam+np.pi/2-b*_deg))
    det = (rA <= D_of(D,np.abs(ca))) & (rB <= D_of(D,np.abs(cb)))
    return float(np.mean(np.sign(ca)[det]*np.sign(cb)[det]))
S_mc = abs(mc(0,22.5)-mc(0,67.5)+mc(45,22.5)+mc(45,67.5))
print(f"\nMC (6M pairs): S = {S_mc:.4f}, E(22.5)={mc(0,22.5):.4f} (QM -0.7071), "
      f"E(10)={mc(0,10):.4f} (QM {-np.cos(20*_deg):.4f}), "
      f"E(70)={mc(0,70):.4f} (QM {-np.cos(140*_deg):.4f})")

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig,(a1,a2,a3)=plt.subplots(1,3,figsize=(15,4.2))
a1.plot(MG, D, color="#7F77DD", lw=2, label="solved D(m)")
a1.plot(MG, MG**max(p_fit,0.01), color="#888", ls=":", lw=1.5, label=f"m^{p_fit:.2f}")
a1.set_xlabel("alignment magnitude m = |cos 2δ|"); a1.set_ylabel("detection probability D(m)")
a1.set_title("Solved sliding profile"); a1.legend(fontsize=8)
Ecurve = np.array([E_D(D,0,t) for t in TH])
a2.plot(TH, Ecurve, color="#1D9E75", lw=2, label="model")
a2.plot(TH, -np.cos(2*TH*_deg), color="#D85A30", ls="--", lw=1.5, label="QM")
a2.set_xlabel("relative angle (deg)"); a2.set_ylabel("E coincidence")
a2.set_title("Coincidence curve"); a2.legend(fontsize=8)
a3.plot(TH, r, color="#1D9E75", lw=1.5); a3.axhline(0,color="#888",ls=":",lw=1)
a3.set_xlabel("relative angle (deg)"); a3.set_ylabel("residual")
a3.set_title(f"Residual (max {np.max(np.abs(r)):.1e})")
fig.tight_layout(); fig.savefig("/home/claude/detection_profile_solution.png", dpi=150)
np.savez("/home/claude/detection_profile.npz", m=MG, D=D)
print("saved plot and profile")
