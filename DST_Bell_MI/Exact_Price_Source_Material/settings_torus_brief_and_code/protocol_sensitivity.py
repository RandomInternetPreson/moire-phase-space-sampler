#!/usr/bin/env python3
"""
Sensitivity analysis for the settings-torus fringe protocol of
arXiv:2608.18886, Sec. 10.

Simulates a 24x24 scan of the fourfold coincidence sum R4(a,b) with
(i) factorizable apparatus ripple A(a)B(b), (ii) Poisson counting noise,
and optionally (iii) an injected cos4(a-b) fringe of peak-to-peak
relative amplitude df.  The estimator is the anti-diagonal (2,-2)
Fourier pixel of log R4, which factorizable responses cannot populate.

Estimator note.  The physical fringe cos4(a-b) has zero phase, so the
estimator is the *signed* projection Re F[2,-2] (not |F|).  A magnitude
estimator is biased positive under noise (|F| >= 0 always) and would
report a spurious non-zero null; the signed estimator is unbiased and its
null is centred on zero.  Significance is reported as
(signal - null) / sigma_null, i.e. against the measured null distribution.

Outputs: null-recovery test, estimator linearity, and the integration
time required for 5-sigma detection at df = 0.1% as a function of pair rate.

Requires numpy only.
"""
import numpy as np

M = 24
rng = np.random.default_rng(20260824)
grid = np.pi*np.arange(M)/M
A, B = np.meshgrid(grid, grid, indexing='ij')

def estimate(counts):
    """Signed peak-to-peak fringe amplitude from the (2,-2) Fourier pixel of log R4,
    plus the quadrature (imaginary) component of the same pixel.

    The quadrature component is the built-in systematic check: the physical fringe
    has zero phase (after relative analyzer calibration), so the quadrature must be
    null-consistent in every dataset, signal or none.  A relative analyzer zero
    offset delta rotates signal into quadrature and costs a factor cos(4*delta)
    in the in-phase amplitude (0.3% at delta = 1 degree)."""
    L = np.log(np.maximum(counts, 1))
    F = np.fft.fft2(L)
    return 4*np.real(F[2, M-2])/M**2, 4*np.imag(F[2, M-2])/M**2

def run(df, N, ripple=0.03, trials=100):
    ests, quads = [], []
    for _ in range(trials):
        epsA = ripple*np.sin(2*grid + rng.uniform(0, 2*np.pi))
        epsB = ripple*np.sin(2*grid + rng.uniform(0, 2*np.pi))
        R = np.exp(epsA[:, None] + epsB[None, :])*(1 + (df/2)*np.cos(4*(A - B)))
        counts = rng.poisson(N*R/R.mean())
        e, q = estimate(counts)
        ests.append(e); quads.append(q)
    return float(np.mean(ests)), float(np.std(ests)), float(np.mean(quads)), float(np.std(quads))

if __name__ == "__main__":
    print("Settings-torus protocol sensitivity (24 x 24, 3% factorizable ripple,\nsigned estimator with quadrature check)\n")

    m, s, qm, qs = run(0.0, 1_000_000, trials=200)
    print(f"Null recovery @ 1e6 counts/point:   in-phase ({m*100:+.4f} +/- {s*100:.4f}) %")
    print(f"                                  quadrature ({qm*100:+.4f} +/- {qs*100:.4f}) %")

    print("\nLinearity (injected vs recovered, 1e6 counts/point):")
    for df in [0.002, 0.005, 0.01, 0.05, 0.124]:
        m, s, qm, qs = run(df, 1_000_000, trials=40)
        print(f"  df = {df*100:6.2f}%  ->  {m*100:6.3f} +/- {s*100:.3f} %   (quadrature {qm*100:+.3f} +/- {qs*100:.3f} %)")

    print("\n5-sigma requirement at df = 0.1% (peak-to-peak), null-corrected:")
    for N in [2e5, 4e5, 6e5, 1e6, 2e6]:
        m0, s0, _, _ = run(0.0, int(N), trials=200)
        m1, s1, _, _ = run(0.001, int(N), trials=100)
        sig = (m1 - m0)/s0
        tot = 576*N
        print(f"  N = {N:.0e}/pt ({tot:.1e} pairs total): null {m0*100:+.4f}+/-{s0*100:.4f} %, "
              f"signal {m1*100:+.4f} %  ->  {sig:4.1f} sigma "
              f"| {tot/5e4/3600:4.1f} h @ 5e4 pairs/s "
              f"| {tot/1e6/60:4.0f} min @ 1e6 pairs/s")

    print("\nPaper's quoted requirement (Sec. 10): ~7 h @ 5e4 pairs/s, ~20 min @ 1e6.")
    print("This Poisson-only simulation reaches 5 sigma several times faster; the paper's")
    print("figure is deliberately conservative and is the one the experimental brief quotes.")
