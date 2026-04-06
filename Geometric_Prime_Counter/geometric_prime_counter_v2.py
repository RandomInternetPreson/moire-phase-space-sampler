#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════╗
║           GEOMETRIC PRIME COUNTER  (standalone)                  ║
║                                                                  ║
║  Count primes from the shape of a surface.                       ║
║  No sieve. No lookup tables. No primes as input.                 ║
║                                                                  ║
║  Pipeline:                                                       ║
║    Modular surface  →  Laplacian eigenvalues                     ║
║    →  approx ζ zeros  →  Newton refinement                       ║
║    →  Chebyshev explicit formula  →  π(x)                        ║
║                                                                  ║
║  Install:  pip install numpy scipy                               ║
║  Run:      python geometric_prime_counter_v2.py                  ║
║                                                                  ║
║  Theory:   Colin de Verdière (1983), Selberg (1956)              ║
╚══════════════════════════════════════════════════════════════════╝
"""

import numpy as np
from math import log, sqrt, pi, isqrt, exp
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
import time
import sys

# ════════════════════════════════════════════════════════════════
# PART 1 — MODULAR SURFACE LAPLACIAN
# Discretise SL(2,Z)\H² truncated at height Y.
# Eigenvalues λ = ¼ + t² → spectral params t → approx ζ zeros 2t.
# ════════════════════════════════════════════════════════════════

def build_surface(Nx=50, Nu=300, Y=50.0):
    """
    Discretise the hyperbolic Laplacian on the truncated modular surface.

    Coordinates: x ∈ (0, ½),  u = log y,  y = e^u.
    Boundary conditions:
      - Dirichlet on the arc  x² + e^{2u} = 1  (lower boundary)
      - Dirichlet at  u = log Y                 (cusp wall)
      - Neumann at  x = 0  and  x = ½           (mod-group identifications)

    The Dirichlet wall at height Y converts the continuous spectrum of the
    open cusp into a discrete spectrum.  This is inspired by Colin de
    Verdière's pseudo-Laplacian construction (1983), though the code uses
    the simpler Dirichlet truncation rather than the pseudo-Laplacian
    itself.

    The spectral parameter t from λ = ¼ + t² relates to zeta zeros via the
    scattering phase condition and Weyl law for the truncated surface.  The
    leading relationship is γ ≈ 2t, accurate enough to seed Newton
    refinement on the Riemann-Siegel Z function.

    Returns dict with keys: eigenvalues, spectral_params, approx_zeros,
                            n_points, build_time_s
    """
    t0 = time.perf_counter()

    u_max = log(Y) - 0.5
    t_p   = np.linspace(-1, 1, Nu)
    u_grid = 0.5*(u_max + 0.14)*(1 + np.tanh(2.0*t_p)/np.tanh(2.0)) - 0.14
    u_grid = np.sort(np.unique(np.clip(u_grid, -0.14, u_max)))
    Nu     = len(u_grid)

    x_grid = np.linspace(0.003, 0.497, Nx)
    dx     = x_grid[1] - x_grid[0]

    # Map (i,j) grid index → flat index, keeping only interior domain
    point_map = {}
    idx_map   = {}
    on_arc    = set()
    idx = 0
    for i in range(Nx):
        for j in range(Nu):
            x, u = x_grid[i], u_grid[j]
            if x**2 + exp(2*u) >= 1.0 - 1e-6:
                point_map[(i, j)] = idx
                idx_map[idx]      = (i, j)
                idx += 1
                if j > 0 and (x_grid[i]**2 + exp(2*u_grid[j-1])) < 1.0:
                    on_arc.add((i, j))

    N = idx
    H = lil_matrix((N, N), dtype=float)

    for k in range(N):
        i, j = idx_map[k]
        u = u_grid[j]

        # Non-uniform u-grid differences
        dp = u_grid[j+1] - u_grid[j]   if 0 < j < Nu-1 else u_grid[1]  - u_grid[0]
        dm = u_grid[j]   - u_grid[j-1] if 0 < j < Nu-1 else u_grid[1]  - u_grid[0]
        if j == 0:   dm = dp
        if j == Nu-1: dp = dm
        ds = dp + dm

        # Finite-difference coefficients for d²/du² and d/du
        c_p = 2/(dp*ds);  c_0 = -2/(dp*dm);  c_m = 2/(dm*ds)
        d_p = dm**2/(dp*dm*ds);  d_m = -dp**2/(dp*dm*ds)
        d_0 = (dp**2 - dm**2)/(dp*dm*ds)

        euu = exp(2*u)  # coefficient of d²/dx² in −Δ_hyp

        # Diagonal
        H[k, k] = -(c_0 + (-1)*d_0) + 2*euu/dx**2

        # u-neighbours
        for dj, cp, dp_ in [(+1, c_p, d_p), (-1, c_m, d_m)]:
            nb = (i, j+dj)
            if nb in point_map:
                H[k, point_map[nb]] = -(cp + (-1)*dp_)
            elif dj == -1 and (i, j) in on_arc and (i, j+1) in point_map:
                H[k, point_map[(i, j+1)]] += -(cp + (-1)*dp_)

        # x-neighbours (Neumann at x=0 and x=½: reflect)
        for di in [-1, +1]:
            nb = (i+di, j)
            if nb in point_map:
                H[k, point_map[nb]] += -euu/dx**2
            elif (di == +1 and i == Nx-1) or (di == -1 and i == 0):
                refl = (i-di, j)
                if refl in point_map:
                    H[k, point_map[refl]] += -euu/dx**2

    H_sym = (csr_matrix(H) + csr_matrix(H).T) / 2

    try:
        evals, _ = eigsh(H_sym, k=80, sigma=0.3, which='LM')
    except Exception:
        evals, _ = eigsh(H_sym, k=80, which='SM')

    evals = np.sort(np.real(evals))
    phys  = evals[evals > 0.26]
    t_vals = np.sqrt(np.maximum(phys - 0.25, 0))

    return {
        'eigenvalues':   phys,
        'spectral_params': t_vals,
        'approx_zeros':  2.0 * t_vals,   # 2t ≈ imaginary part of ζ zero
        'n_points':      N,
        'build_time_s':  time.perf_counter() - t0,
    }


# ════════════════════════════════════════════════════════════════
# PART 2 — NEWTON REFINEMENT (Riemann-Siegel Z function)
# Snap each rough zero candidate to a true zero of ζ(½+it).
# ════════════════════════════════════════════════════════════════

def _rs_theta(t):
    """Riemann-Siegel phase θ(t).  Accurate for t > 10."""
    return (t/2*log(t/(2*pi)) - t/2 - pi/8
            + 1/(48*t) + 7/(5760*t**3))

def rs_Z(t):
    """Real-valued Z(t); zeros coincide with nontrivial zeros of ζ(½+it)."""
    if t < 10:
        return 0.0
    N     = int(sqrt(t/(2*pi)))
    theta = _rs_theta(t)
    total = sum(np.cos(theta - t*log(n)) / sqrt(n) for n in range(1, N+1))
    frac  = sqrt(t/(2*pi)) - N
    corr  = np.cos(2*pi*(frac**2 - frac - 1/16)) / np.cos(2*pi*frac)
    return 2*total + (-1)**(N-1) * corr*(2*pi/t)**0.25

def refine_zero(t0, search_radius=2.0, step=0.1):
    """
    Find the true zero of Z(t) nearest to t0.
    Returns refined t, or None if no bracket found.
    """
    if t0 < 10:
        return None

    # Bracket search
    Z0 = rs_Z(t0)
    t_lo = t_hi = None
    for direction in [1, -1]:
        t, Zprev = t0, Z0
        for _ in range(int(search_radius/step) + 1):
            t += direction*step
            if t < 10: break
            Zcurr = rs_Z(t)
            if Zprev * Zcurr < 0:
                t_lo, t_hi = (t-direction*step, t) if direction == 1 else (t, t-direction*step)
                if t_lo > t_hi: t_lo, t_hi = t_hi, t_lo
                break
            Zprev = Zcurr
        if t_lo is not None:
            break

    if t_lo is None:
        return None

    # Bisect
    for _ in range(60):
        tm = (t_lo + t_hi)/2
        if t_hi - t_lo < 1e-9: break
        if rs_Z(t_lo)*rs_Z(tm) < 0: t_hi = tm
        else:                        t_lo = tm

    # Newton polish
    t = (t_lo + t_hi)/2
    for _ in range(20):
        Zt  = rs_Z(t)
        dZt = (rs_Z(t+1e-5) - rs_Z(t-1e-5)) / 2e-5
        if abs(dZt) < 1e-15: break
        dt = -Zt/dZt
        t += dt
        if abs(dt) < 1e-12: break

    return t if abs(rs_Z(t)) < 0.01 else None

def refine_all(candidates, verbose=True):
    """Refine a list of candidate zeros. Returns sorted deduplicated list."""
    refined = []
    n = len(candidates)
    for i, t0 in enumerate(candidates):
        if verbose:
            sys.stdout.write(f"\r  Refining {i+1}/{n} (t={t0:.2f}) ...    ")
            sys.stdout.flush()
        tr = refine_zero(t0)
        if tr is not None and tr > 10:
            refined.append(tr)
    if verbose:
        sys.stdout.write("\r" + " "*50 + "\r")

    refined = sorted(set(round(t, 9) for t in refined))
    deduped = [refined[0]] if refined else []
    for t in refined[1:]:
        if t - deduped[-1] > 0.5:
            deduped.append(t)
    return deduped


# ════════════════════════════════════════════════════════════════
# PART 3 — PRIME COUNTING (Chebyshev explicit formula)
# ψ(x) = x − log 2π − Σ_γ  oscillating correction
# π(x) ≈ ψ(x) / log x
# ════════════════════════════════════════════════════════════════

def count_primes(x, zeros, max_zeros=None):
    """
    Estimate π(x) from imaginary parts of ζ zeros via the explicit formula.

    Parameters
    ----------
    x        : positive real, the bound to count primes up to
    zeros    : iterable of imaginary parts γ of nontrivial ζ zeros
    max_zeros: use only the first max_zeros zeros (fewer = faster, less accurate)

    Returns
    -------
    Estimated π(x) as a float.
    """
    if x < 2:
        return 0.0
    gammas = list(zeros)[:max_zeros] if max_zeros else list(zeros)

    psi = x - log(2*pi)
    lx  = log(x)
    sqx = sqrt(x)
    for g in gammas:
        if g < 10: continue
        psi -= 2*sqx*(0.5*np.cos(g*lx) + g*np.sin(g*lx)) / (0.25 + g**2)

    return max(0.0, psi / lx)

def sieve_exact(N):
    """Exact prime count via Sieve of Eratosthenes."""
    if N < 2: return 0
    s = bytearray(b'\x01') * (N+1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(N)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return sum(s)


# ════════════════════════════════════════════════════════════════
# MAIN — build, refine, compare
# ════════════════════════════════════════════════════════════════

def main():
    print()
    print("╔" + "═"*60 + "╗")
    print("║" + "  GEOMETRIC PRIME COUNTER".center(60) + "║")
    print("║" + "  Primes from the shape of a surface".center(60) + "║")
    print("╚" + "═"*60 + "╝")
    print()

    # ── Step 1: build surface ────────────────────────────────────
    print("  Step 1 — Building modular surface SL(2,Z)\\H²")
    surf = build_surface(Nx=50, Nu=300, Y=50.0)
    print(f"  Done in {surf['build_time_s']*1000:.0f}ms  "
          f"({surf['n_points']} grid points, "
          f"{len(surf['approx_zeros'])} candidate zeros)")
    print()

    # ── Step 2: Newton refinement ────────────────────────────────
    print("  Step 2 — Newton refinement (snapping to true ζ zeros)")
    t0       = time.perf_counter()
    refined  = refine_all(surf['approx_zeros'], verbose=True)
    refine_s = time.perf_counter() - t0
    n_cand   = len(surf['approx_zeros'])
    n_conv   = len(refined)
    pct_fail = (1 - n_conv/n_cand)*100 if n_cand > 0 else 0
    print(f"  Done in {refine_s*1000:.0f}ms  "
          f"({n_conv} zeros converged, {pct_fail:.0f}% of candidates failed)")
    print()

    # ── Step 3: prime counting ───────────────────────────────────
    print("  Step 3 — Prime counting")
    print()
    test_x = [100, 1_000, 10_000, 100_000, 1_000_000]

    print(f"  {'x':>10}  {'π(x) exact':>11}  "
          f"{'Geometric':>11}  {'+ Newton':>10}  "
          f"{'Geo err':>8}  {'Newton err':>10}")
    print("  " + "─"*70)

    for x in test_x:
        exact   = sieve_exact(x)
        geo     = count_primes(x, surf['approx_zeros'])
        newton  = count_primes(x, refined) if refined else 0.0
        ge      = abs(geo    - exact)/exact*100
        ne      = abs(newton - exact)/exact*100 if refined else 999
        flag    = "  ✓" if ne < 2.0 else ""
        print(f"  {x:>10,}  {exact:>11,}  "
              f"{geo:>11.0f}  {newton:>10.0f}  "
              f"{ge:>7.1f}%  {ne:>9.1f}%{flag}")

    print()

    # ── Timing ───────────────────────────────────────────────────
    N_REP = 1000
    x_time = 1_000_000

    t0 = time.perf_counter()
    for _ in range(N_REP): sieve_exact(x_time)
    sieve_ms = (time.perf_counter()-t0)/N_REP*1000

    t0 = time.perf_counter()
    for _ in range(N_REP): count_primes(x_time, refined or surf['approx_zeros'])
    geo_ms = (time.perf_counter()-t0)/N_REP*1000

    print(f"  Timing (x = {x_time:,}):")
    print(f"    Sieve (exact, per call):  {sieve_ms:.3f} ms")
    print(f"    Geometric query:          {geo_ms:.3f} ms")
    total_setup_ms = (surf['build_time_s'] + refine_s)*1000
    if sieve_ms > geo_ms:
        breakeven = int(total_setup_ms / (sieve_ms - geo_ms))
        print(f"    One-time setup:           {total_setup_ms:.0f} ms")
        print(f"    Break-even:               ~{breakeven:,} queries")
    print()

    # ── What happened ────────────────────────────────────────────
    print("  " + "─"*60)
    print("  The pipeline in plain language:")
    print()
    print("  1. Construct the modular surface SL(2,Z)\\H² — the")
    print("     quotient of the hyperbolic plane by the modular group.")
    print("     Place a Dirichlet wall in the cusp at height Y.")
    print()
    print("  2. Solve for the eigenvalues of the hyperbolic Laplacian.")
    print("     By the scattering phase mechanism (Colin de Verdière,")
    print("     1983), these approximate the imaginary parts of")
    print("     Riemann zeta zeros.")
    print()
    print("  3. Use Newton's method on the Riemann-Siegel Z function")
    print("     to snap each rough geometric zero to a true zero.")
    print("     The Z function requires no primes — only calculus.")
    print()
    print("  4. Feed the refined zeros into Chebyshev's explicit formula")
    print("     to estimate π(x).  No primes entered at any step.")
    print()
    print("  References:")
    print("    Colin de Verdière (1983)  Ann. Inst. Fourier 33, 87")
    print("    Selberg (1956)            J. Indian Math. Soc. 20, 47")
    print()


if __name__ == "__main__":
    main()
