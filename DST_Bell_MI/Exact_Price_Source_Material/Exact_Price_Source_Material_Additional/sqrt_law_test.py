#!/usr/bin/env python3
"""
sqrt_law_test.py - Consistency test of the derived edge law D(m) ~ sqrt(m).
Re-solve the exact-reproduction problem with D(m) = sqrt(m) * h(m), h smooth
on a grid, h >= 0. If the sqrt law is right, this constrained family reaches
the same ~1e-5 exactness as the free solve, with finite h(0).
"""
import numpy as np
from scipy.optimize import least_squares

_deg = np.pi/180
NLAM = 4320
LAM = np.linspace(0, np.pi, NLAM, endpoint=False)
TH = np.linspace(0.0, 180.0, 91)
M = 41
HG = np.linspace(0.0, 1.0, M)

def D_of(h, m):
    return np.sqrt(m) * np.interp(m, HG, h)

def E_D(h, delta):
    ca = np.cos(2*LAM); cb = np.cos(2*(LAM + np.pi/2 - delta*_deg))
    wa = D_of(h, np.abs(ca)); wb = D_of(h, np.abs(cb))
    ww = wa*wb
    den = ww.sum()
    return float((np.sign(ca)*np.sign(cb)*ww).sum()/den) if den > 1e-12 else 0.0

def residuals(h):
    return np.array([E_D(h, t) + np.cos(2*t*_deg) for t in TH])

sol = least_squares(residuals, np.ones(M), bounds=(0.0, 5.0),
                    xtol=1e-14, ftol=1e-14, max_nfev=4000)
h = sol.x
r = residuals(h)
print(f"constrained solve (D = sqrt(m) h(m)):")
print(f"  max |E + cos 2D| = {np.max(np.abs(r)):.2e}   rms = {np.sqrt(np.mean(r**2)):.2e}")
print(f"  (free-D benchmark was 1.3e-05; matching it confirms the sqrt law)")
print(f"  h(0) = {h[0]:.5f}  h(0.25) = {np.interp(0.25,HG,h):.5f}  "
      f"h(0.5) = {np.interp(0.5,HG,h):.5f}  h(1) = {h[-1]:.5f}")
print(f"  h monotone decreasing: {bool(np.all(np.diff(h) <= 1e-3))}")
print(f"  h(0)/h(1) = {h[0]/h[-1]:.5f}")
# simple-form probes for h (normalized to h(1)=1)
hn = h / h[-1]
for name, f in [("1/(1+a m): a from h(0)", None)]:
    pass
a_est = 1/hn[0]  # if hn = c/(1+am): hn(0)=c, hn(1)=c/(1+a)=1
print(f"  probe hn(m) ~ hn(0)/(1+a m): implied a = hn(0)-1 = {hn[0]-1:.5f}")
pred = hn[0]/(1+(hn[0]-1)*HG)
print(f"  max |hn - probe| = {np.max(np.abs(hn - pred)):.4f}")
np.savez("/home/claude/sqrt_law_h.npz", m=HG, h=h)
