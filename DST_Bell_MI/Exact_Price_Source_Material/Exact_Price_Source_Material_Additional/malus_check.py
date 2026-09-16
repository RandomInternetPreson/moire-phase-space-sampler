#!/usr/bin/env python3
"""
malus_check.py  --  Does the pure-detection profile of arXiv:2608.18886 Sec. 5
reproduce Malus's law for one photon through two sequential polarizers?

Model class (static hidden phase): lambda ~ U[0, pi); analyzer at angle a routes
by sgn cos 2(lambda - a); click probability D(m), m = |cos 2(lambda - a)|.
Two completions for a compound station (polarizer 1 at 0, polarizer 2 at theta,
detector behind polarizer 2's + port):
  variant 1: polarizers lossless, D applied once at the detector (m relative to pol 2)
  variant 2: D applied at each polarizer (undecided photon lost at the crystal)
Reference intensity is the parallel configuration (theta = 0).

Analytic fact (variant 1): D(m) = m gives cos^2 theta exactly, since
  int_{theta-pi/4}^{pi/4} cos 2(lambda-theta) dlambda = cos^2 theta.
Edge exponents: singlet (Lemma 1) forces D ~ sqrt(m); Malus forces D ~ m. Incompatible.

Run from Exact_Price_Source_Material_Additional/ (needs data/h_spectral.npz).
Prepared 2026-09-16 in session with Claude; single instance -- cold-audit before use.
"""
import numpy as np
from numpy.polynomial import chebyshev as C

c = np.load("data/h_spectral.npz")["cheb"]
h = lambda m: C.chebval(2*m - 1, c)
PROFILES = {
    "solved sqrt(m)h(m)": lambda m: np.sqrt(m) * h(m),
    "pure sqrt(m)":       lambda m: np.sqrt(m),
    "linear m":           lambda m: m,
    "constant":           lambda m: np.ones_like(m),
}

lam = np.linspace(0, np.pi, 400001)[:-1]

def R1(th, D):
    c1 = np.cos(2*lam); c2 = np.cos(2*(lam - th)); sel = (c1 > 0) & (c2 > 0)
    return np.mean(np.where(sel, D(np.abs(c2)), 0)) / np.mean(np.where(c1 > 0, D(np.abs(c1)), 0))

def R2(th, D):
    c1 = np.cos(2*lam); c2 = np.cos(2*(lam - th)); sel = (c1 > 0) & (c2 > 0)
    return (np.mean(np.where(sel, D(np.abs(c1))*D(np.abs(c2)), 0))
            / np.mean(np.where(c1 > 0, D(np.abs(c1))**2, 0)))

def singlet(D):
    """E(Delta) = -C_f(2Delta)/C_g(2Delta) and normalized coincidence rate C_g/C_g(0)."""
    x = np.linspace(0, 2*np.pi, 400001)[:-1]
    g = D(np.abs(np.cos(x))); f = np.sign(np.cos(x))*g
    G = np.fft.rfft(g); F = np.fft.rfft(f); n = len(x)
    Cg = np.fft.irfft(G*np.conj(G), n=n); Cf = np.fft.irfft(F*np.conj(F), n=n)
    idx = lambda th: int(round(2*th/(2*np.pi)*n)) % n
    E = lambda th: -Cf[idx(th)]/Cg[idx(th)]
    rate = np.array([Cg[idx(t)]/Cg[0] for t in np.linspace(0, np.pi/2, 181)])
    S = abs(E(np.pi/8) - E(3*np.pi/8) + 2*E(np.pi/8))
    return E(np.pi/8), S, rate.max() - rate.min()

ths = np.deg2rad(np.linspace(0, 90, 361))
print(f"{'profile':22s} {'var1 max|R-cos2|':>17s} {'var2 max|R-cos2|':>17s} {'E(22.5)':>8s} {'S':>6s} {'rate p2p':>9s}")
for name, D in PROFILES.items():
    d1 = max(abs(R1(t, D) - np.cos(t)**2) for t in ths)
    d2 = max(abs(R2(t, D) - np.cos(t)**2) for t in ths)
    E, S, mod = singlet(D)
    print(f"{name:22s} {d1:17.4f} {d2:17.4f} {E:8.4f} {S:6.3f} {mod:9.3f}")
print("\nQM: Malus deviation 0, E(22.5)=-0.7071, S=2.828, rate modulation 0.")
