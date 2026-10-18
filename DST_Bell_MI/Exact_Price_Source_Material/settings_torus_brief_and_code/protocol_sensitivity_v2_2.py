#!/usr/bin/env python3
"""
Sensitivity analysis for the settings-torus fringe protocol of
arXiv:2608.18886, Sec. 10 -- version 2.2 (October 2026).

Changes from v1, following the review of 13 Oct 2026:

  1. The apparatus model now has OUTCOME-DEPENDENT channel efficiencies
     eta_+(a) != eta_-(a) on each side.  For a maximally entangled source
         R4(a,b) = (1/4) [ S_A(a) S_B(b) + D_A(a) D_B(b) E(a,b) ],
     S = eta_+ + eta_-,  D = eta_+ - eta_-,  E = -V cos 2(a-b).
     The second term is local but NOT factorizable, so it leaks into the
     (2,-2) pixel of log R4: at second order for a constant relative
     imbalance delta = D/S (amplitude -V^2 delta_A^2 delta_B^2 / 4, opposite
     sign to the physical fringe) and at first order for an angle-dependent
     imbalance d1 cos(2a - p) (amplitude -V d1 e1 / 4 when the phases
     agree).  v1 injected only the factorizable ripple exp(eps_A + eps_B)
     and therefore never tested this.

  2. The CALIBRATED estimator.  For a maximally entangled source the local
     states are maximally mixed, so the two channel singles on each side
     give delta_A(a) and delta_B(b) directly.  The four coincidence
     channels, re-weighted by the singles-derived efficiencies, give the
     correlator E(a,b).  R4 is then divided by (1 + delta_A delta_B E),
     which removes the non-factorizable term, and the residual is purely
     factorizable apparatus response times the physical fringe.  The
     (2,-2) pixel of its logarithm is the estimator.  The naive v1
     estimator is kept for comparison.

  3. A coincidence-marginal CROSS-CHECK of the singles calibration: at relative
     angle 45 deg (E = 0) Alice's coincidence marginal equals delta_A exactly, so the
     coincidence table gives a second estimate; agreement with the singles estimate
     checks that singles and coincidence photons see the same channel efficiencies.

  4. (v2.2) The EQUAL-PHASE worst case: the angle-dependent imbalance on the two
     sides with the same phase (p = q), where its (2,-2) contribution is purely
     real and the quadrature check is blind.  Section E fixes p = q = 0 instead
     of drawing the phases at random; sections A-D are unchanged and reproduce
     v2.1's output exactly.

  5. The 2x2 CROSS-RATIO estimator at relative angles 0 and 45 degrees:
         X = R4(0,0) R4(45,45) / [ R4(0,45) R4(45,0) ],
     which cancels every factorizable response and equals 1.28 for the
     Sec. 5 profile.  Any passive-basis BBM92 setup records this table.
     Reported with the same calibration.

Fringe shape: the exact Sec. 5 profile, as the first ten harmonics
a_k cos 4k(a-b) of R4/mean computed from the solved D(m) = sqrt(m) h(m);
f = a_1 = 0.0577 is the full profile (12.39% peak-to-peak), smaller f
emulates its mixtures with flat detection.

Requires numpy only.
"""
import numpy as np

M = 24
rng = np.random.default_rng(20261013)
grid = np.pi*np.arange(M)/M
A, B = np.meshgrid(grid, grid, indexing='ij')
V = 1.0          # source visibility

# ---------------------------------------------------------------- apparatus
def apparatus(ripple=0.03, d0=0.0, d1=0.0, phase=None):
    """Per-side channel efficiencies eta_+(a), eta_-(a) (scale arbitrary).
    ripple: factorizable |S| ripple amplitude; d0: constant relative channel
    imbalance; d1: angle-dependent imbalance amplitude.  Phases random, unless
    phase is given: then both sides' imbalance phases equal it (p = q)."""
    def side():
        S = np.exp(ripple*np.sin(2*grid + rng.uniform(0, 2*np.pi)))
        ph = rng.uniform(0, 2*np.pi)
        if phase is not None:
            ph = phase
        delta = d0 + d1*np.cos(2*grid - ph)
        return S*(1 + delta)/2, S*(1 - delta)/2          # eta_+, eta_-
    return side(), side()

# Exact Sec. 5 fringe: relative amplitudes a_k of cos(4k Delta) in R4/mean,
# from the solved profile D(m) = sqrt(m) h(m) (solve_h_spectral.py, residual 1.5e-12).
# a_1 = 0.05773 gives the 12.39% peak-to-peak of the paper; a_1/a_2 = 5.87.
A_K = np.array([5.773052e-02, 9.834086e-03, 3.100640e-03, 1.331055e-03, 6.853649e-04,
                3.973659e-04, 2.503605e-04, 1.677163e-04, 1.177647e-04, 8.582444e-05])
F_SEC5 = A_K[0]

def fringe(f):
    """Fringe with cos-4 amplitude f and the exact Sec. 5 harmonic shape
    (f = F_SEC5 is the full profile; f < F_SEC5 emulates its mixtures with flat detection)."""
    s = f/F_SEC5
    return 1 + s*sum(ak*np.cos(4*k*(A - B)) for k, ak in enumerate(A_K, start=1))

def simulate(f, N, ripple, d0, d1, phase=None):
    """Returns Poisson counts: coincidences in the four channels (4,M,M)
    and channel singles on each side (2,M).  N = mean coincidences/point."""
    (apA, amA), (bpB, bmB) = apparatus(ripple, d0, d1, phase)
    etaA = {+1: apA, -1: amA}; etaB = {+1: bpB, -1: bmB}
    E = -V*np.cos(2*(A - B))
    P = {(i, j): (1 + i*j*E)/4 for i in (+1, -1) for j in (+1, -1)}
    F = fringe(f)
    # pairs emitted per point, normalised so mean coincidences/point = N
    R4 = sum(etaA[i][:, None]*etaB[j][None, :]*P[(i, j)]*F for (i, j) in P)
    Np = N/R4.mean()
    C = {ij: rng.poisson(Np*etaA[ij[0]][:, None]*etaB[ij[1]][None, :]*P[ij]*F) for ij in P}
    # singles: each channel sees half the pairs times its efficiency (maximally mixed local state)
    sA = {i: rng.poisson(Np*M*etaA[i]/2) for i in (+1, -1)}   # summed over partner settings
    sB = {j: rng.poisson(Np*M*etaB[j]/2) for j in (+1, -1)}
    return C, sA, sB

# ---------------------------------------------------------------- estimators
def pixel22(L):
    F = np.fft.fft2(L)
    return 2*np.real(F[2, M-2])/M**2, 2*np.imag(F[2, M-2])/M**2   # cos4 amplitude, quadrature

def naive(C):
    R4 = sum(C.values())
    return pixel22(np.log(np.maximum(R4, 1)))

def calibrated(C, sA, sB):
    dA = (sA[+1] - sA[-1])/(sA[+1] + sA[-1])            # delta_A(a) from singles
    dB = (sB[+1] - sB[-1])/(sB[+1] + sB[-1])
    wA = {i: 1/(1 + i*dA) for i in (+1, -1)}             # 1/eta up to a common factor
    wB = {j: 1/(1 + j*dB) for j in (+1, -1)}
    num = sum(i*j*C[(i, j)]*wA[i][:, None]*wB[j][None, :] for (i, j) in C)
    den = sum(C[(i, j)]*wA[i][:, None]*wB[j][None, :] for (i, j) in C)
    E = num/np.maximum(den, 1)                            # efficiency-corrected correlator
    R4 = sum(C.values())
    Rc = R4/(1 + dA[:, None]*dB[None, :]*E)               # remove the S S + D D E term
    return pixel22(np.log(np.maximum(Rc, 1)))

def marginal_check(C, sA, sB):
    """Second estimate of delta_A from the coincidence table at relative angle 45 deg
    (E = 0 there, so Alice's coincidence marginal equals delta_A exactly).  Returns the
    rms difference between the singles-based and coincidence-based delta_A over a,
    a check that singles and coincidence photons see the same channel efficiencies."""
    dA_s = (sA[+1] - sA[-1])/(sA[+1] + sA[-1])
    k = M//4
    cols = [(a, (a + k) % M) for a in range(M)]            # pairs (a, b = a + 45 deg)
    num = np.array([sum(i*C[(i, j)][a, b] for (i, j) in C) for a, b in cols])
    den = np.array([sum(C[(i, j)][a, b] for (i, j) in C) for a, b in cols])
    dA_c = num/np.maximum(den, 1)
    return float(np.sqrt(np.mean((dA_s - dA_c)**2)))

def cross_ratio(C, sA, sB):
    """2x2 cross-ratio from the torus points a,b in {0, 45 deg}, calibrated.
    Returns X - 1, which equals 4 f (1 - f/6 ...) ~ 4f for small f."""
    dA = (sA[+1] - sA[-1])/(sA[+1] + sA[-1]); dB = (sB[+1] - sB[-1])/(sB[+1] + sB[-1])
    wA = {i: 1/(1 + i*dA) for i in (+1, -1)}; wB = {j: 1/(1 + j*dB) for j in (+1, -1)}
    num = sum(i*j*C[(i, j)]*wA[i][:, None]*wB[j][None, :] for (i, j) in C)
    den = sum(C[(i, j)]*wA[i][:, None]*wB[j][None, :] for (i, j) in C)
    E = num/np.maximum(den, 1)
    Rc = sum(C.values())/(1 + dA[:, None]*dB[None, :]*E)
    k = M//4                                              # index of 45 degrees
    return Rc[0, 0]*Rc[k, k]/(Rc[0, k]*Rc[k, 0]) - 1

def run(f, N, ripple=0.03, d0=0.0, d1=0.0, trials=100, phase=None):
    out = {'naive': [], 'naive_q': [], 'cal': [], 'cal_q': [], 'x': []}
    for _ in range(trials):
        C, sA, sB = simulate(f, N, ripple, d0, d1, phase)
        e, q = naive(C); out['naive'].append(e); out['naive_q'].append(q)
        e, q = calibrated(C, sA, sB); out['cal'].append(e); out['cal_q'].append(q)
        out['x'].append(cross_ratio(C, sA, sB))
        out.setdefault('mchk', []).append(marginal_check(C, sA, sB))
    return {k: (float(np.mean(v)), float(np.std(v))) for k, v in out.items()}

def pct(ms):
    return f"{ms[0]*100:+.4f} +/- {ms[1]*100:.4f} %"

if __name__ == "__main__":
    print("Settings-torus protocol sensitivity v2.2 (24 x 24; cos4 amplitude f, peak-to-peak = 2f)\n")

    print("A. Null recovery at 1e6 coincidences/point, 3% factorizable ripple:")
    for d0, d1, label in [(0.0, 0.0, "no imbalance          "),
                          (0.20, 0.0, "20% const imbalance   "),
                          (0.05, 0.03, "5% + 3% ripple imbal. "),
                          (0.20, 0.03, "20% + 3% ripple imbal.")]:
        r = run(0.0, 1_000_000, d0=d0, d1=d1, trials=100)
        print(f"  {label}  naive {pct(r['naive'])}   calibrated {pct(r['cal'])}   marginal-check rms {r['mchk'][0]:.4f}"
              f"   (predicted naive bias: const {-d0**4/4*100:+.4f}%, ripple ~{-d1**2/4*100:+.4f}% at equal phase)")

    print("\nB. Linearity, calibrated estimator, 20% + 3% imbalance, 1e6/point:")
    for f in [0.0005, 0.001, 0.005, F_SEC5]:
        r = run(f, 1_000_000, d0=0.20, d1=0.03, trials=40)
        print(f"  f = {f*100:6.3f}%  ->  {pct(r['cal'])}   quadrature {pct(r['cal_q'])}   cross-ratio-1 = {r['x'][0]:+.5f} (4f = {4*f:+.5f})")

    print("\nC. 5-sigma requirement at f = 0.05% (0.1% peak-to-peak), calibrated, 20% + 3% imbalance:")
    for N in [2e5, 4e5, 6e5, 1e6]:
        r0 = run(0.0, int(N), d0=0.20, d1=0.03, trials=200)
        r1 = run(0.0005, int(N), d0=0.20, d1=0.03, trials=100)
        sig = (r1['cal'][0] - r0['cal'][0])/r0['cal'][1]
        sigx = (r1['x'][0] - r0['x'][0])/r0['x'][1]
        tot = 576*N
        print(f"  N = {N:.0e}/pt ({tot:.1e} coinc. total): torus {sig:4.1f} sigma | "
              f"2x2 from the same 4 points {sigx:4.1f} sigma | {tot/5e4/3600:4.1f} h @ 5e4/s")

    print("\nD. Dedicated 2x2 (four settings only), 20% + 3% imbalance:")
    for N in [1e7, 2.5e7]:
        # emulate: only the four corner points matter; run with that N per point
        r0 = run(0.0, int(N), d0=0.20, d1=0.03, trials=200)
        r1 = run(0.0005, int(N), d0=0.20, d1=0.03, trials=100)
        sigx = (r1['x'][0] - r0['x'][0])/r0['x'][1]
        print(f"  N = {N:.1e}/setting (4N = {4*N:.1e} total): {sigx:4.1f} sigma | {4*N/5e4/3600:4.1f} h @ 5e4/s")

    print("\nE. Equal-phase worst case (p = q = 0), 1e6 coincidences/point, 3% factorizable ripple:")
    for d0, d1, label in [(0.05, 0.03, "5% + 3% ripple, p = q  "),
                          (0.20, 0.03, "20% + 3% ripple, p = q "),
                          (0.20, 0.01, "20% + 1% ripple, p = q ")]:
        r = run(0.0, 1_000_000, d0=d0, d1=d1, trials=100, phase=0.0)
        print(f"  {label} naive {pct(r['naive'])} (quadrature {pct(r['naive_q'])})   "
              f"calibrated {pct(r['cal'])} (quadrature {pct(r['cal_q'])})   "
              f"(predicted naive: {(-d0**4/4 - d1**2/4)*100:+.4f}%)")

    print("\nNotes: (i) the naive estimator's bias is the non-factorizable S S + D D E term; "
          "(ii) the calibrated\nestimator removes it using the channel singles; (iii) the 2x2 cross-ratio "
          "needs ~1e8 coincidences\nfor 5 sigma at 0.1% against 2.3e8 for the full torus.")
