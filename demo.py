"""
demo.py
-------
Generate all figures for the moiré Wigner reconstruction project.

Usage:
    python demo.py

Outputs (saved to ./figures/):
    reconstruction_fock1.png
    reconstruction_fock2.png
    reconstruction_cat.png
    reconstruction_coherent.png
    reconstruction_squeezed.png
    convergence.png
    chi_structure.png

Each reconstruction figure shows:
    True W(q,p) | |χ(s,τ)| | Moiré T | Sampled χ bins | Reconstructed Ŵ | Error
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')   # no display required

from moire_wigner import (
    STATES,
    plot_reconstruction,
    plot_convergence,
    plot_chi_structure,
    make_grid,
    wigner_fock,
    characteristic_function,
    sampling_mask,
    reconstruct_wigner,
    reconstruction_stats,
    beat_to_bin,
)

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

os.makedirs('figures', exist_ok=True)

N    = 64
QMAX = 3.5
F1   = 1.5      # fine screen spatial frequency
F2   = 1.0      # coarse screen spatial frequency
NP   = 15       # number of screen pair configurations (chi coverage)

q, p, Q, dq = make_grid(N, QMAX)

print('=' * 60)
print('MOIRÉ WIGNER RECONSTRUCTION  ·  DEMO')
print('=' * 60)
print(f'Grid:    {N}×{N},  q ∈ [{-QMAX}, {QMAX}],  dq = {dq:.4f}')
print(f'Screens: f₁={F1}, f₂={F2}, Δf={F1-F2:.2f}')
print(f'         Beat bin radius = {beat_to_bin(F1-F2, N, dq):.1f} DFT bins')
print(f'N pairs: {NP}  →  disk r ≤ {1+0.5*NP:.1f} bins')
print()

# ---------------------------------------------------------------------------
# 1. Reconstruction figures — one per state
# ---------------------------------------------------------------------------

for key in STATES:
    print(f'Generating reconstruction: {key} ...')
    fig, stats = plot_reconstruction(
        state_key    = key,
        N=N, qmax=QMAX,
        f1=F1, f2=F2,
        n_pairs      = NP,
        use_ring_mask= False,
        savepath     = f'figures/reconstruction_{key}.png',
        show         = False,
    )
    import matplotlib.pyplot as plt
    plt.close(fig)

    print(f'  RMSE = {stats["rmse"]:.5f}')
    if not np.isnan(stats['neg_fidelity']):
        print(f'  neg-vol recovery = {100*stats["neg_fidelity"]:.1f}%')
    else:
        print(f'  no negative volume (classical state)')
    print()

# ---------------------------------------------------------------------------
# 2. Convergence curves
# ---------------------------------------------------------------------------

print('Generating convergence curves ...')
import matplotlib.pyplot as plt
fig = plot_convergence(
    n_pairs_range = range(1, 51, 2),
    savepath      = 'figures/convergence.png',
    show          = False,
)
plt.close(fig)
print('  Done.')
print()

# ---------------------------------------------------------------------------
# 3. Characteristic function structure
# ---------------------------------------------------------------------------

print('Generating chi structure panel ...')
fig = plot_chi_structure(
    savepath = 'figures/chi_structure.png',
    show     = False,
)
plt.close(fig)
print('  Done.')
print()

# ---------------------------------------------------------------------------
# 4. Key numbers printed to terminal
# ---------------------------------------------------------------------------

print('=' * 60)
print('KEY NUMBERS  (validating against experimental literature)')
print('=' * 60)

W_fock1 = wigner_fock(q, p, n=1)
origin_val = wigner_fock(
    np.array([[0.0]]), np.array([[0.0]]), n=1
)[0, 0]

print(f'\nFock |1⟩:')
print(f'  W(0,0) = {origin_val:.6f}')
print(f'  Theory: -1/π = {-1/np.pi:.6f}')
print(f'  Experiment (Lvovsky 2001): -0.318 ± 0.030')
print(f'  Match: {"YES ✓" if abs(origin_val - (-1/np.pi)) < 0.01 else "NO ✗"}')

print(f'\nFock |1⟩ reconstruction with {NP} screen pairs:')
chi   = characteristic_function(W_fock1)
mask  = sampling_mask(N, NP)
W_rec = reconstruct_wigner(chi, mask)
stats = reconstruction_stats(W_fock1, W_rec, dq)
print(f'  True neg-vol  = {stats["neg_vol_true"]:.4f}')
print(f'  Rec  neg-vol  = {stats["neg_vol_rec"]:.4f}')
print(f'  Recovery      = {100*stats["neg_fidelity"]:.1f}%')
print(f'  RMSE          = {stats["rmse"]:.5f}')

print(f'\nCoverage at N_pairs={NP}:')
print(f'  χ bins sampled = {mask.sum()} / {N*N}')
print(f'  Coverage       = {100*mask.sum()/N**2:.1f}%')
print(f'  r_max          = {1+0.5*NP:.1f} DFT bins')

print('\n' + '=' * 60)
print('All figures saved to ./figures/')
print('=' * 60)
