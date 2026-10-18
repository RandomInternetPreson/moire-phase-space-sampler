#!/usr/bin/env python3
"""
solve_h_spectral.py - Solve the detection-sector master equation exactly
(to machine precision) for h(m), where D(m) = sqrt(m) h(m).

Arc form of the master equation, phi(t) = D(|sin t|):
    4 (phi * phi)(theta) = (1 - cos theta) A(theta),   theta in (0, pi)
    (phi*phi)(theta) = int_0^theta phi(t) phi(theta - t) dt
    A(theta) = 2 int_0^pi phi(t) phi(t - theta) dt   (phi even, period pi)

All integrands have exponent-1/2 endpoint singularities exactly where a phi
vanishes -> Gauss-Jacobi(1/2,1/2) quadrature is spectrally exact for the
smooth remainders. h is a Chebyshev series on m in [0,1], scale fixed h(1)=1.
"""
import numpy as np
from numpy.polynomial import chebyshev as C
from scipy.special import roots_jacobi
from scipy.optimize import least_squares

NJ = 64
xj, wj = roots_jacobi(NJ, 0.5, 0.5)     # weight (1-x)^{1/2} (1+x)^{1/2} on [-1,1]

def jacobi_int(a, b, smooth):
    """int_a^b (t-a)^{1/2} (b-t)^{1/2} smooth(t) dt via Gauss-Jacobi."""
    L = 0.5*(b - a)
    t = a + L*(xj + 1.0)
    return (L**2) * np.dot(wj, smooth(t))

def h_of(c, m):
    return C.chebval(2.0*m - 1.0, c)

def sqrt_sinc(u):
    """sqrt(sin u / u), smooth on [0, pi), safe at 0."""
    u = np.asarray(u, float)
    out = np.ones_like(u)
    nz = u > 1e-12
    out[nz] = np.sqrt(np.sin(u[nz])/u[nz])
    return out

def sqrt_sin_over_pi_minus(u):
    """sqrt(sin u / (pi - u)), smooth near pi, general u in (0, pi)."""
    return np.sqrt(np.sin(u)/(np.pi - u))

def residuals(c, thetas):
    res = []
    for th in thetas:
        # convolution term on [0, th]: sing at both ends
        def s_conv(t):
            return (sqrt_sinc(t)*h_of(c, np.sin(t)) *
                    sqrt_sinc(th - t)*h_of(c, np.sin(th - t)))
        conv = jacobi_int(0.0, th, s_conv)
        # A(theta): segment [0, th]: phi(t)~ends? phi(t) sing at t=0; phi(t-th)=phi(th-t) sing at t=th
        def s_A1(t):
            return (sqrt_sinc(t)*h_of(c, np.sin(t)) *
                    sqrt_sinc(th - t)*h_of(c, np.sin(th - t)))
        A1 = jacobi_int(0.0, th, s_A1)
        # segment [th, pi]: phi(t-th) sing at t=th; phi(t) sing at t=pi
        def s_A2(t):
            return (sqrt_sinc(t - th)*h_of(c, np.sin(t - th)) *
                    sqrt_sin_over_pi_minus(t)*h_of(c, np.sin(t)))
        A2 = jacobi_int(th, np.pi, s_A2)
        A = 2.0*(A1 + A2)
        res.append(4.0*conv - (1.0 - np.cos(th))*A)
    res.append(50.0*(h_of(c, 1.0) - 1.0))     # scale: h(1) = 1
    return np.array(res)

thetas = np.linspace(0.02, np.pi - 0.02, 81)

# validate machinery against the previously solved D (expect ~1e-5-scale residuals)
dat = np.load("sqrt_law_h.npz")
h_prev = lambda m: np.interp(m, dat["m"], dat["h"])
c0 = C.chebfit(2*dat["m"] - 1, dat["h"]/dat["h"][-1], 16)
r0 = residuals(c0, thetas)[:-1]
print(f"machinery validation on previous h: max|R| = {np.max(np.abs(r0)):.2e} "
      f"(previous solve precision was ~1e-5-ish: expected small but nonzero)")

NC = 25
cinit = np.zeros(NC); cinit[:len(c0)] = c0
sol = least_squares(lambda c: residuals(c, thetas), cinit,
                    xtol=3e-16, ftol=3e-16, gtol=3e-16, max_nfev=6000)
c = sol.x
r = residuals(c, thetas)[:-1]
print(f"\nspectral solve: max|R| = {np.max(np.abs(r)):.2e}   rms = {np.sqrt(np.mean(r**2)):.2e}")

h0 = h_of(c, 0.0)
print(f"\nh(0)/h(1) = {h0:.14f}")
print(f"3*pi/8    = {3*np.pi/8:.14f}")
print(f"difference: {h0 - 3*np.pi/8:.3e}")
verdict = "CONFIRMED" if abs(h0 - 3*np.pi/8) < 1e-8 else "REJECTED"
print(f"=> 3*pi/8 hypothesis: {verdict}")

print(f"\nh at reference points (h(1)=1):")
for mv in [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]:
    print(f"  h({mv:5.3f}) = {h_of(c, mv):.12f}")

# closed-form probes with PSLQ over small bases
import mpmath
mpmath.mp.dps = 30
def probe(val, basis, names, maxc=2000):
    rel = mpmath.pslq([mpmath.mpf(1)] + [mpmath.mpf(b) for b in basis] + [mpmath.mpf(val)],
                      tol=mpmath.mpf(1e-10), maxcoeff=maxc)
    if rel is None or rel[-1] == 0:
        return None
    terms = " + ".join(f"({-r}/{rel[-1]})*{n}" for r, n in zip(rel[:-1], ["1"]+names) if r != 0)
    return terms
for label, mv in [("h(0)", 0.0), ("h(1/2)", 0.5), ("h(1/4)", 0.25)]:
    v = h_of(c, mv)
    hit = probe(v, [np.pi, np.pi**2, np.sqrt(2), np.pi*np.sqrt(2)],
                ["pi", "pi^2", "sqrt2", "pi*sqrt2"])
    print(f"  PSLQ {label} = {v:.12f}: {hit if hit else 'no small relation found'}")

np.savez("h_spectral.npz", cheb=c, h0=h0)
print("\nsaved h_spectral.npz")

# ---------------------------------------------------------------------------
# Exact-profile harmonics of the fourfold coincidence sum (Secs. 6, 9, 10).
# R4(a,b) is the circular autocorrelation C_g(2*Delta) of g(x) = D(|cos x|),
# period pi, so R4/mean = 1 + sum_k a_k cos(4k*Delta) with a_k = g_k^2/(2 g_0^2),
# g_k the cos(2kx) Fourier coefficients of g. The sqrt endpoint singularity at
# x = pi/2 is removed by the substitution x = (pi/2) cos^2 v.
# ---------------------------------------------------------------------------
from numpy.polynomial.legendre import leggauss
_v, _w = leggauss(600)
_v = (_v + 1.0)*np.pi/4; _w = _w*np.pi/4
_x = (np.pi/2)*np.cos(_v)**2
_dx = (np.pi/2)*2*np.sin(_v)*np.cos(_v)
_g = np.sqrt(np.abs(np.cos(_x)))*h_of(c, np.abs(np.cos(_x)))
_g0 = np.sum(_w*_dx*_g)/(np.pi/2)
_gk = np.array([2*np.sum(_w*_dx*_g*np.cos(2*k*_x))/(np.pi/2) for k in range(1, 13)])
a_k = (_gk**2/2)/_g0**2
k = np.arange(1, 13)
R4 = lambda d: 1 + np.sum(a_k*np.cos(4*k*d))
print("\nexact-profile harmonics a_k of R4/mean = 1 + sum a_k cos(4k Delta):")
for i, ak in enumerate(a_k, start=1):
    print(f"  a_{i:<2d} = {ak:.6e}")
print(f"  a_1/a_2            = {a_k[0]/a_k[1]:.3f}")
print(f"  peak-to-peak       = {R4(0.0)-R4(np.pi/4):.4f}  (of the mean)")
print(f"  R4(0)/R4(45 deg)   = {R4(0.0)/R4(np.pi/4):.4f}")
print(f"  R4(22.5 deg)/mean  = {R4(np.pi/8):.4f}   (= R4(67.5 deg)/mean: {R4(3*np.pi/8):.4f})")
print(f"  cross-ratio X      = {R4(0.0)*R4(0.0)/(R4(np.pi/4)**2):.4f}")
print(f"  max-scale efficiency (mean click probability, h(1)=1) = {_g0:.4f}")
np.save("fringe_harmonics.npy", a_k)
print("saved fringe_harmonics.npy")
