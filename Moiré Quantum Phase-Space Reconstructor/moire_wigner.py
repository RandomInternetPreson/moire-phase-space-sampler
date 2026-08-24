"""
moire_wigner.py
---------------
Moiré-based quantum phase-space reconstruction.

Core physics:
    The moiré beat between two periodic screens samples the Wigner
    characteristic function χ(s,τ) at the beat frequency. N screen
    pairs at different (Δf, θ) sample N arcs in χ-space. Sparse iFFT
    recovers the Wigner function W(q,p).

    The key identity (not an analogy):
        I_moiré ∝ Re[ψ̃*(p) · ψ̃(p + ℏΔf)] = χ(0, 2πΔf)

    See DERIVATION.md for the full proof.

Wigner functions are validated against experimental measurements:
    Fock |1⟩:  W(0,0) = -1/π ≈ -0.318  [Lvovsky & Babichev, PRL 87 (2001)]
    Fock |2⟩:  negative annular ring     [Ourjoumtsev et al., Science 312 (2006)]
    Cat state: interference fringes W<0  [Deleglise et al., Nature 455 (2008)]
    Coherent:  Gaussian, W≥0            [Smithey et al., PRL 70 (1993)]
    Squeezed:  elliptical Gaussian, W≥0  [Breitenbach et al., Nature 387 (1997)]
"""

import numpy as np
from scipy.special import genlaguerre
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec


# ---------------------------------------------------------------------------
# Grid setup
# ---------------------------------------------------------------------------

def make_grid(N=64, qmax=3.5):
    """
    Create phase-space coordinate grids.
    
    Parameters
    ----------
    N    : grid size (N×N)
    qmax : extent in both q and p directions
    
    Returns
    -------
    q, p : 2D meshgrids (shape N×N)
    Q    : 1D coordinate array (length N)
    dq   : grid spacing
    """
    dq = 2 * qmax / N
    Q  = np.linspace(-qmax + dq/2, qmax - dq/2, N)
    q, p = np.meshgrid(Q, Q, indexing='ij')
    return q, p, Q, dq


# ---------------------------------------------------------------------------
# Wigner functions  (ħ = ½ convention: W_vacuum peak = 1/π)
# ---------------------------------------------------------------------------

def wigner_coherent(q, p, alpha=1.5):
    """
    Wigner function for coherent state |α⟩ with real displacement α.

        W(q,p) = (1/π) exp(-(q - √2 α)² - p²)

    Always non-negative. Classical benchmark.
    Ref: Smithey et al., PRL 70, 1244 (1993).
    """
    c = np.sqrt(2) * alpha
    return (1/np.pi) * np.exp(-(q - c)**2 - p**2)


def wigner_fock(q, p, n=1):
    """
    Wigner function for Fock (number) state |n⟩.

        W_n(q,p) = (1/π) (-1)^n  L_n(2r²) exp(-r²)

    where r² = q² + p² and L_n is the Laguerre polynomial.

    n=1: W(0,0) = -1/π ≈ -0.318  (experimentally confirmed negative peak)
    n=2: negative annular ring with positive outer lobe

    Refs:
        n=1: Lvovsky & Babichev, PRL 87, 050402 (2001)
        n=2: Ourjoumtsev, Tualle-Brouri & Grangier, Science 312 (2006)
    """
    r2 = q**2 + p**2
    Ln = genlaguerre(n, 0)(2 * r2)
    return (1/np.pi) * ((-1)**n) * Ln * np.exp(-r2)


def wigner_cat(q, p, alpha=2.0):
    """
    Wigner function for even cat state (|α⟩ + |-α⟩) / N  with real α.

        W_cat = [W_α + W_{-α} + (2/π) exp(-r²) cos(2√2 α p)] / (2(1 + e^{-2α²}))

    The cross-interference term (2/π)exp(-r²)cos(2√2αp) produces fringes
    between the two coherent blobs with W < 0 in the troughs.
    Fringe period in p: 2π / (2√2 α).

    Ref: Deleglise, Dotsenko et al., Nature 455, 510 (2008).
    """
    c   = np.sqrt(2) * alpha
    r2  = q**2 + p**2
    nm  = 2 * (1 + np.exp(-2 * alpha**2))
    W1  = (1/np.pi) * np.exp(-(q - c)**2 - p**2)
    W2  = (1/np.pi) * np.exp(-(q + c)**2 - p**2)
    Wc  = (2/np.pi) * np.exp(-r2) * np.cos(2 * c * p)
    return (W1 + W2 + Wc) / nm


def wigner_squeezed(q, p, r=1.0):
    """
    Wigner function for squeezed vacuum with squeezing parameter r.

        W(q,p) = (1/π) exp(-q² e^{2r} - p² e^{-2r})

    Squeezed below shot noise in q, above in p. Always W ≥ 0.
    Ref: Breitenbach, Schiller & Mlynek, Nature 387, 471 (1997).
    """
    return (1/np.pi) * np.exp(-q**2 * np.exp(2*r) - p**2 * np.exp(-2*r))


# Catalogue used by demo
STATES = {
    'fock1': {
        'label':  'Fock |1⟩',
        'fn':     lambda q, p: wigner_fock(q, p, n=1),
        'vmax':   0.35,
        'ref':    'Lvovsky & Babichev, PRL 87, 050402 (2001)',
        'note':   'W(0,0) = -1/π ≈ -0.318  ← experimentally confirmed negative peak',
        'color':  '#ff6060',
    },
    'fock2': {
        'label':  'Fock |2⟩',
        'fn':     lambda q, p: wigner_fock(q, p, n=2),
        'vmax':   0.35,
        'ref':    'Ourjoumtsev et al., Science 312 (2006)',
        'note':   'Negative annular ring with positive outer lobe',
        'color':  '#ff9030',
    },
    'cat': {
        'label':  'Cat |α=2⟩',
        'fn':     lambda q, p: wigner_cat(q, p, alpha=2.0),
        'vmax':   0.18,
        'ref':    'Deleglise et al., Nature 455 (2008)',
        'note':   'Interference fringes with W<0 between coherent blobs',
        'color':  '#60d0ff',
    },
    'coherent': {
        'label':  'Coherent |α=1.5⟩',
        'fn':     lambda q, p: wigner_coherent(q, p, alpha=1.5),
        'vmax':   0.35,
        'ref':    'Smithey et al., PRL 70, 1244 (1993)',
        'note':   'Gaussian, W≥0 everywhere  ← classical baseline',
        'color':  '#60ff60',
    },
    'squeezed': {
        'label':  'Squeezed Vacuum (r=1)',
        'fn':     lambda q, p: wigner_squeezed(q, p, r=1.0),
        'vmax':   0.35,
        'ref':    'Breitenbach, Schiller & Mlynek, Nature 387 (1997)',
        'note':   'Elliptical Gaussian, noise below shot noise in q',
        'color':  '#d0ff60',
    },
}


# ---------------------------------------------------------------------------
# Characteristic function  χ(s,τ) = FFT2D(W)
# ---------------------------------------------------------------------------

def characteristic_function(W):
    """
    Compute the 2D Wigner characteristic function.

        χ(s,τ) = ∬ W(q,p) e^{i(sq+τp)} dq dp

    Uses numpy FFT2 (O(N² log N) — fast).
    Returns complex array, same shape as W.
    fftshift is NOT applied here; use np.fft.fftshift for display.
    """
    return np.fft.fft2(W)


def chi_magnitude_display(chi):
    """Log-scaled magnitude of χ, fftshifted for display (DC at centre)."""
    mag = np.abs(np.fft.fftshift(chi))
    mag = np.where(mag < 1e-12, 1e-12, mag)
    lmag = np.log10(mag)
    lmag -= lmag.max()          # normalise to [log_floor, 0]
    return np.clip((lmag + 5) / 5, 0, 1)   # map [-5, 0] → [0, 1]


# ---------------------------------------------------------------------------
# Moiré screen model
# ---------------------------------------------------------------------------

def moire_transmission(q, p, f1=1.5, f2=1.0):
    """
    Combined transmission of two 2D square-weave screens.

        T(q,p) = [(1 + cos(2πf₁q))(1 + cos(2πf₁p)) ×
                  (1 + cos(2πf₂q))(1 + cos(2πf₂p))] / 16

    Range: [0, 1].  Beat frequency Δf = f₁ - f₂ sets which χ bins are sampled.
    """
    T1 = (1 + np.cos(2*np.pi*f1*q)) * (1 + np.cos(2*np.pi*f1*p)) / 4
    T2 = (1 + np.cos(2*np.pi*f2*q)) * (1 + np.cos(2*np.pi*f2*p)) / 4
    return T1 * T2


def beat_to_bin(delta_f, N, dq):
    """
    Convert a physical screen beat frequency to a DFT bin radius.

    Physical beat: s = 2π·Δf  (rad / unit)
    DFT bin:       m = s · N·dq / (2π) = Δf · N·dq

    With N=64, dq=7/64:  bin = Δf × 7
    """
    return delta_f * N * dq


# ---------------------------------------------------------------------------
# Sampling mask  (disk in χ-space)
# ---------------------------------------------------------------------------

def sampling_mask(N, n_pairs):
    """
    Build a binary sampling mask in DFT (χ) space.

    Represents the cumulative coverage from n_pairs screen configurations,
    each with a different (Δf, θ).  Models the disk of DFT bins accessible
    when you vary both the beat frequency and the screen orientation.

    Disk radius: r_max = 1 + 0.5 × n_pairs  (in DFT bin units)

    Physically: the moiré beat at Δf samples DFT bin r = Δf × N·dq.
    Varying Δf from near-zero to Δf_max traces out radii 0 → r_max.
    Varying θ over [0,π) fills each annulus completely.

    Parameters
    ----------
    N       : grid size
    n_pairs : number of screen pair configurations

    Returns
    -------
    mask    : bool array, shape (N, N), True = sampled
    """
    r_max = 1.0 + 0.5 * n_pairs
    mask  = np.zeros((N, N), dtype=bool)
    for m in range(N):
        for n_ in range(N):
            mf = m if m <= N//2 else m - N
            nf = n_ if n_ <= N//2 else n_ - N
            if np.sqrt(mf**2 + nf**2) <= r_max:
                mask[m, n_] = True
    return mask


def sampling_mask_rings(N, dq, f1=1.5, f2=1.0, n_pairs=15, n_harmonics=3):
    """
    Physically motivated ring mask: samples χ at the actual moiré beat
    frequencies and their harmonics, at n_pairs different angles.

    This is the 'honest' mask — it only includes (s,τ) points that a
    real screen pair at (f1, f2) would actually sample.

    Parameters
    ----------
    N           : grid size
    dq          : grid spacing
    f1, f2      : screen spatial frequencies (cycles per unit)
    n_pairs     : number of angular orientations
    n_harmonics : how many harmonics to include (1=fundamental only)
    """
    mask = np.zeros((N, N), dtype=bool)
    mask[0, 0] = True   # DC always sampled (uniform transmission)

    df      = abs(f1 - f2)
    freqs   = [df * h for h in range(1, n_harmonics+1)]
    freqs  += [f1, f2]
    freqs   = [f for f in freqs if 0 < beat_to_bin(f, N, dq) < N//2 - 0.5]

    for k in range(n_pairs):
        theta = (k + 0.5) * np.pi / n_pairs
        ct, st = np.cos(theta), np.sin(theta)
        for freq in freqs:
            r = beat_to_bin(freq, N, dq)
            for sign in (+1, -1):
                m  = int(round(sign * r * ct))
                n_ = int(round(sign * r * st))
                mask[m % N, n_ % N] = True

    return mask


# ---------------------------------------------------------------------------
# Reconstruction
# ---------------------------------------------------------------------------

def reconstruct_wigner(chi, mask):
    """
    Reconstruct W from a sparse subset of χ samples.

        Ŵ = iFFT2(χ · mask)

    When mask is all-True this is exact.  Partial mask gives a
    low-pass filtered reconstruction whose fidelity grows with coverage.

    Parameters
    ----------
    chi  : complex array (N×N), output of characteristic_function()
    mask : bool array (N×N), from sampling_mask() or sampling_mask_rings()

    Returns
    -------
    W_rec : real float array (N×N)
    """
    chi_masked = np.where(mask, chi, 0.0 + 0.0j)
    W_rec      = np.fft.ifft2(chi_masked).real
    return W_rec


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def reconstruction_stats(W_true, W_rec, dq):
    """
    Compute reconstruction quality metrics.

    Returns dict with:
        rmse        : RMS error
        neg_vol_true: negative Wigner volume of true state
        neg_vol_rec : negative Wigner volume of reconstruction
        neg_fidelity: fraction of true negative volume recovered
    """
    rmse         = float(np.sqrt(np.mean((W_true - W_rec)**2)))
    neg_vol_true = float(np.sum(W_true[W_true < 0] * (-dq**2)))
    neg_vol_rec  = float(np.sum(W_rec[W_rec   < 0] * (-dq**2)))
    neg_fidelity = neg_vol_rec / neg_vol_true if neg_vol_true > 0 else float('nan')
    return {
        'rmse':          rmse,
        'neg_vol_true':  neg_vol_true,
        'neg_vol_rec':   neg_vol_rec,
        'neg_fidelity':  neg_fidelity,
    }


def convergence_curve(W_true, chi, dq, n_pairs_range=None):
    """
    Compute RMSE and negative-volume recovery as a function of n_pairs.

    Parameters
    ----------
    W_true       : true Wigner function
    chi          : its characteristic function
    dq           : grid spacing
    n_pairs_range: iterable of n_pairs values to evaluate

    Returns
    -------
    results : list of dicts (one per n_pairs value)
    """
    if n_pairs_range is None:
        n_pairs_range = range(1, 51, 2)

    N = W_true.shape[0]
    results = []
    for np_ in n_pairs_range:
        mask  = sampling_mask(N, np_)
        W_rec = reconstruct_wigner(chi, mask)
        stats = reconstruction_stats(W_true, W_rec, dq)
        stats['n_pairs']   = np_
        stats['n_bins']    = int(mask.sum())
        stats['coverage']  = mask.sum() / N**2
        results.append(stats)
    return results


# ---------------------------------------------------------------------------
# Colormaps
# ---------------------------------------------------------------------------

def wigner_cmap():
    """
    Diverging colormap for Wigner functions.
    Blue = W < 0 (quantum),  near-black = W = 0,  amber = W > 0.
    """
    colors = [
        (0.00, '#0022cc'),   # deep blue  (most negative)
        (0.35, '#001444'),   # dark blue
        (0.50, '#060a10'),   # near-black (zero)
        (0.65, '#331400'),   # dark amber
        (1.00, '#ff8800'),   # bright amber (most positive)
    ]
    positions = [c[0] for c in colors]
    hex_colors = [c[1] for c in colors]
    return mcolors.LinearSegmentedColormap.from_list(
        'wigner', list(zip(positions, hex_colors)))


def chi_cmap():
    """Hot colormap for |χ| magnitude: black → red → yellow → white."""
    return plt.cm.hot


# ---------------------------------------------------------------------------
# High-level plotting
# ---------------------------------------------------------------------------

def plot_reconstruction(state_key='fock1', N=64, qmax=3.5, f1=1.5, f2=1.0,
                        n_pairs=15, use_ring_mask=False,
                        savepath=None, show=True):
    """
    Six-panel figure showing the full reconstruction pipeline for one state.

    Panels: True W | |χ| | Moiré | Sampled χ bins | Reconstructed Ŵ | Error

    Parameters
    ----------
    state_key      : key in STATES dict
    N, qmax        : grid parameters
    f1, f2         : screen spatial frequencies
    n_pairs        : number of screen pair configurations
    use_ring_mask  : if True, use physically-motivated ring mask instead of disk
    savepath       : if set, save figure to this path
    show           : if True, call plt.show()
    """
    state = STATES[state_key]
    q, p, Q, dq = make_grid(N, qmax)

    W_true  = state['fn'](q, p)
    chi     = characteristic_function(W_true)
    T       = moire_transmission(q, p, f1, f2)

    if use_ring_mask:
        mask = sampling_mask_rings(N, dq, f1, f2, n_pairs)
    else:
        mask = sampling_mask(N, n_pairs)

    W_rec   = reconstruct_wigner(chi, mask)
    stats   = reconstruction_stats(W_true, W_rec, dq)

    n_bins  = int(mask.sum())
    coverage = 100 * n_bins / N**2

    # --- colormaps ---
    wcmap = wigner_cmap()
    ccmap = chi_cmap()
    vmax  = state['vmax']

    fig = plt.figure(figsize=(14, 9), facecolor='#05080f')
    fig.patch.set_facecolor('#05080f')
    gs  = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.25)

    panel_data = [
        (W_true,                          wcmap,  -vmax, vmax, 'True  W(q,p)',       'Analytical · ground truth'),
        (chi_magnitude_display(chi),       ccmap,  0,    1,    '|χ(s,τ)|',           'Characteristic function · log scale · DC centred'),
        (T,                                'gray', 0,    1,    f'Moiré  T(q,p)²',    f'f₁={f1}  f₂={f2}  Δf={abs(f1-f2):.2f}'),
        (np.fft.fftshift(mask).astype(float), None, 0, 1, 'Sampled χ bins',     f'{n_bins}/{N*N} ({coverage:.1f}%)'),
        (W_rec,                            wcmap,  -vmax, vmax, 'Reconstructed  Ŵ',  f'RMSE = {stats["rmse"]:.5f}'),
        (np.abs(W_true - W_rec),           ccmap,  0, None,    '|W − Ŵ|  Residual', 'Error magnitude'),
    ]

    axes = []
    for idx, (data, cmap, vmin, vmx, title, sub) in enumerate(panel_data):
        ax = fig.add_subplot(gs[idx // 3, idx % 3])
        ax.set_facecolor('#05080f')

        if cmap is None:   # mask panel — cyan / dark
            rgb = np.zeros((*data.shape, 3))
            rgb[data > 0.5] = [0, 0.86, 1.0]
            rgb[data < 0.5] = [0.04, 0.04, 0.09]
            ax.imshow(rgb.transpose(1, 0, 2), origin='lower',
                      extent=[-qmax, qmax, -qmax, qmax])
        else:
            im = ax.imshow(data.T, origin='lower', cmap=cmap,
                           vmin=vmin, vmax=vmx,
                           extent=[-qmax, qmax, -qmax, qmax],
                           interpolation='nearest')
            plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02,
                         label='', format='%.2f').ax.tick_params(
                         labelsize=6, colors='#556677')

        ax.set_title(title, color=state['color'], fontsize=8,
                     fontfamily='monospace', pad=4)
        ax.set_xlabel('q  (position)', color='#334455', fontsize=7)
        ax.set_ylabel('p  (momentum)', color='#334455', fontsize=7)
        ax.tick_params(colors='#334455', labelsize=6)
        for spine in ax.spines.values():
            spine.set_edgecolor('#1a2535')
        ax.text(0.5, -0.16, sub, transform=ax.transAxes,
                color='#2a4050', fontsize=6.5, ha='center',
                fontfamily='monospace')
        axes.append(ax)

    # Title block
    title_str  = f'MOIRÉ WIGNER RECONSTRUCTION  ·  {state["label"]}'
    stats_str  = (f'neg-vol (true) = {stats["neg_vol_true"]:.4f}   '
                  f'neg-vol (rec) = {stats["neg_vol_rec"]:.4f}   '
                  f'recovery = {100*stats["neg_fidelity"]:.1f}%')
    ref_str    = f'Ref: {state["ref"]}   |   {state["note"]}'

    fig.text(0.5, 0.97, title_str, ha='center', va='top',
             color='#00d4ff', fontsize=12, fontfamily='monospace',
             fontweight='bold')
    fig.text(0.5, 0.93, stats_str, ha='center', va='top',
             color='#446688', fontsize=8, fontfamily='monospace')
    fig.text(0.5, 0.005, ref_str, ha='center', va='bottom',
             color='#2a4050', fontsize=7.5, fontfamily='monospace',
             style='italic')

    if savepath:
        fig.savefig(savepath, dpi=150, bbox_inches='tight',
                    facecolor='#05080f')
        print(f'Saved: {savepath}')
    if show:
        plt.show()
    return fig, stats


def plot_convergence(state_keys=None, N=64, qmax=3.5,
                     n_pairs_range=None, savepath=None, show=True):
    """
    Plot RMSE and negative-volume recovery vs n_pairs for all states.

    Shows the compressed sensing argument: quantum states require far fewer
    measurements than a classical uniform sampling strategy.
    """
    if state_keys is None:
        state_keys = list(STATES.keys())
    if n_pairs_range is None:
        n_pairs_range = list(range(1, 51, 2))

    q, p, Q, dq = make_grid(N, qmax)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), facecolor='#05080f')
    for ax in (ax1, ax2):
        ax.set_facecolor('#05080f')
        ax.tick_params(colors='#556677', labelsize=8)
        for spine in ax.spines.values():
            spine.set_edgecolor('#1a2535')
        ax.grid(True, color='#0d1a26', linewidth=0.5)

    for key in state_keys:
        state   = STATES[key]
        W_true  = state['fn'](q, p)
        chi     = characteristic_function(W_true)
        results = convergence_curve(W_true, chi, dq, n_pairs_range)

        coverages  = [r['coverage'] * 100 for r in results]
        rmses      = [r['rmse']           for r in results]
        neg_fids   = [r['neg_fidelity'] * 100
                      if not np.isnan(r['neg_fidelity']) else 0
                      for r in results]

        ax1.plot(coverages, rmses,      color=state['color'], lw=1.5,
                 label=state['label'])
        ax2.plot(coverages, neg_fids,   color=state['color'], lw=1.5,
                 label=state['label'])

    ax1.set_xlabel('χ-space coverage (%)',    color='#556677', fontsize=9)
    ax1.set_ylabel('RMSE',                    color='#556677', fontsize=9)
    ax1.set_title('Reconstruction Error vs Coverage',
                  color='#00d4ff', fontsize=10, fontfamily='monospace')
    ax1.legend(fontsize=7, facecolor='#0a0f18', edgecolor='#1a2535',
               labelcolor='#aabbd0')

    ax2.set_xlabel('χ-space coverage (%)',           color='#556677', fontsize=9)
    ax2.set_ylabel('Negative volume recovered (%)',  color='#556677', fontsize=9)
    ax2.set_title('Quantum Signature Recovery vs Coverage',
                  color='#00d4ff', fontsize=10, fontfamily='monospace')
    ax2.axhline(100, color='#334455', ls='--', lw=0.8, label='perfect')
    ax2.legend(fontsize=7, facecolor='#0a0f18', edgecolor='#1a2535',
               labelcolor='#aabbd0')

    fig.text(0.5, 0.97,
             'MOIRÉ WIGNER RECONSTRUCTION  ·  CONVERGENCE',
             ha='center', color='#00d4ff', fontsize=11,
             fontfamily='monospace', fontweight='bold')
    fig.text(0.5, 0.01,
             'Classical states (coherent, squeezed) have no negative volume to recover.',
             ha='center', color='#2a4050', fontsize=8,
             fontfamily='monospace', style='italic')

    fig.tight_layout(rect=[0, 0.04, 1, 0.94])

    if savepath:
        fig.savefig(savepath, dpi=150, bbox_inches='tight',
                    facecolor='#05080f')
        print(f'Saved: {savepath}')
    if show:
        plt.show()
    return fig


def plot_chi_structure(state_keys=None, N=64, qmax=3.5,
                       savepath=None, show=True):
    """
    Show |χ(s,τ)| for each state to illustrate why sparse sampling works.

    The characteristic function decays rapidly for sparse quantum states,
    concentrating information near the origin. This is the compressed
    sensing prior — and it's physically guaranteed, not assumed.
    """
    if state_keys is None:
        state_keys = list(STATES.keys())

    q, p, Q, dq = make_grid(N, qmax)
    n            = len(state_keys)
    fig, axes    = plt.subplots(1, n, figsize=(3*n, 3.5), facecolor='#05080f')
    if n == 1:
        axes = [axes]

    for ax, key in zip(axes, state_keys):
        state  = STATES[key]
        W      = state['fn'](q, p)
        chi    = characteristic_function(W)
        mag    = chi_magnitude_display(chi)

        ax.imshow(mag.T, origin='lower', cmap='hot',
                  vmin=0, vmax=1, interpolation='bilinear',
                  extent=[-N//2, N//2, -N//2, N//2])
        ax.set_title(state['label'], color=state['color'],
                     fontsize=9, fontfamily='monospace')
        ax.set_xlabel('s-bin', color='#334455', fontsize=7)
        ax.set_ylabel('τ-bin', color='#334455', fontsize=7)
        ax.tick_params(colors='#334455', labelsize=6)
        ax.set_facecolor('#05080f')
        for spine in ax.spines.values():
            spine.set_edgecolor('#1a2535')

    fig.text(0.5, 0.97, '|χ(s,τ)|  CHARACTERISTIC FUNCTIONS  (log scale)',
             ha='center', color='#00d4ff', fontsize=11,
             fontfamily='monospace', fontweight='bold')
    fig.text(0.5, 0.01,
             'All information is concentrated near the DC bin. '
             'This is why sparse sampling works.',
             ha='center', color='#2a4050', fontsize=8,
             fontfamily='monospace', style='italic')

    fig.tight_layout(rect=[0, 0.05, 1, 0.94])

    if savepath:
        fig.savefig(savepath, dpi=150, bbox_inches='tight',
                    facecolor='#05080f')
        print(f'Saved: {savepath}')
    if show:
        plt.show()
    return fig
