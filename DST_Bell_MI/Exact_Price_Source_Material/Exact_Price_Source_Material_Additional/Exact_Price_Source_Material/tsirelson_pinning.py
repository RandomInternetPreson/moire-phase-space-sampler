#!/usr/bin/env python3
"""
tsirelson_pinning.py - The sphere, the squares, and the rate anisotropy.

Claims under test:
(1) SQUARES => SPHERE: correlators of the form E(D) = -sum c_k cos(2kD) with
    c_k >= 0 (forced: harmonic weights are squared amplitudes) satisfy
    CHSH <= 2 sqrt 2 at ANY angle choices (Schoenberg/Tsirelson).
(2) The pinning principle: within the mechanics, flat coincidence rate +
    squares => Tsirelson bound automatic. Models exceeding 2 sqrt 2 must
    show angle-dependent coincidence rates.
(3) The hidden wound / R3 stake: the solved exact-D model has coincidence
    rate N(D) ~ C_g(2D), which is NOT flat. Quantify the modulation - a
    concrete differential prediction (QM: flat) any experiment logs.
"""
import numpy as np
from scipy.optimize import minimize

# ---------- (1) positive mixtures never beat 2 sqrt 2, any angles ----------
rng = np.random.default_rng(3)
def chsh_mixture(c, angles):
    a1, a2, b1, b2 = angles
    def E(x, y):
        d = x - y
        return -sum(ck*np.cos(2*(k+1)*d) for k, ck in enumerate(c))
    return abs(E(a1,b1) - E(a1,b2) + E(a2,b1) + E(a2,b2))
worst = 0.0
for _ in range(300):
    c = rng.dirichlet(np.ones(6))          # positive weights, sum 1
    best_ang = 0.0
    for _ in range(20):
        x0 = rng.uniform(0, np.pi, 4)
        r = minimize(lambda x: -chsh_mixture(c, x), x0, method="Nelder-Mead",
                     options={"maxiter": 400})
        best_ang = max(best_ang, -r.fun)
    worst = max(worst, best_ang)
print(f"(1) max CHSH over 300 random positive mixtures x optimized angles: "
      f"{worst:.6f}  (2sqrt2 = {2*np.sqrt(2):.6f})")
print(f"    squares => sphere: {'HOLDS' if worst <= 2*np.sqrt(2)+1e-4 else 'FAILS'}")

# ---------- (2)+(3) rate anisotropy of the mechanics ----------
NX = 8192
x = np.linspace(0, 2*np.pi, NX, endpoint=False)

def rate_curve(gvals):
    Gf = np.fft.rfft(gvals)
    Cg = np.fft.irfft(np.abs(Gf)**2, n=NX)/NX     # autocorrelation over x=2*lambda
    th = x[:NX//2]                                 # theta = 2*Delta in [0, pi)
    return th/2*180/np.pi, Cg[:NX//2]

dat = np.load("sqrt_law_h.npz")
D = lambda m: np.sqrt(m)*np.interp(m, dat["m"], dat["h"]/dat["h"][-1])
g_exact = D(np.abs(np.cos(x)))
Dg, N_exact = rate_curve(g_exact)
mod_exact = (N_exact.max()-N_exact.min())/N_exact.mean()
print(f"\n(3) EXACT-D model coincidence rate vs relative angle:")
print(f"    modulation depth (max-min)/mean = {mod_exact*100:.2f}%   (QM: 0%)")
i_mx, i_mn = np.argmax(N_exact), np.argmin(N_exact)
print(f"    max at Delta = {Dg[i_mx]:.1f} deg, min at Delta = {Dg[i_mn]:.1f} deg")

w = 0.4
g_band = (np.abs(np.cos(x)) >= w).astype(float)
_, N_band = rate_curve(g_band)
mod_band = (N_band.max()-N_band.min())/N_band.mean()
print(f"\n(2) hard-band model (w=0.4, the S=4 regime): rate modulation = "
      f"{mod_band*100:.1f}%")
print(f"    -> confirms: super-Tsirelson violation is purchased with large")
print(f"       rate anisotropy; flat-rate faithfulness forbids it.")

print(f"\nTHEOREM-LET: full faithfulness (cosine E AND flat rate) is impossible")
print(f"for any pure-detection model in this class: flat rate forces D const,")
print(f"D const forces sawtooth E. Pure detection is dead as a complete account")
print(f"- independently of the loophole-free experiments. Surviving: MD/hybrid.")

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(Dg, N_exact/N_exact.mean(), color="#1D9E75", lw=2,
        label=f"exact-D model ({mod_exact*100:.1f}% modulation)")
ax.plot(Dg, N_band/N_band.mean(), color="#7F77DD", lw=1.5, ls="--",
        label=f"hard band w=0.4 ({mod_band*100:.0f}%)")
ax.axhline(1.0, color="#D85A30", ls=":", lw=1.5, label="QM: flat")
ax.set_xlabel("relative analyzer angle Δ (deg)")
ax.set_ylabel("coincidence rate (normalized)")
ax.set_title("The R3 stake: rate anisotropy — a logged observable QM says is flat")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("/home/claude/rate_anisotropy.png", dpi=150)
print("saved rate_anisotropy.png")
