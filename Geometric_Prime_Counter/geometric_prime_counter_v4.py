#!/usr/bin/env python3
r"""
╔══════════════════════════════════════════════════════════════════╗
║        PRIME-FREE PRIME COUNTER  v4  (standalone)                ║
║                                                                  ║
║  Count primes with no sieve, no lookup tables, no primes as      ║
║  input — and test, with controls, whether a truncated modular    ║
║  surface contributes anything to the zeros that do the work.     ║
║                                                                  ║
║  Pipeline that works (prime-free):                               ║
║    Riemann–Siegel Z(t)  →  ζ zeros  →  Riemann's J(x)  →  π(x)   ║
║                                                                  ║
║  Pipeline under test (v1–v3 claim, now a NEGATIVE result):       ║
║    Modular surface  →  Laplacian eigenvalues  →  seeds for the   ║
║    above.  Controls in PART 4 show the surface eigenvalues do    ║
║    not locate zeros better than a ruler laid across the same     ║
║    interval; Newton + Riemann–Siegel finds the zeros from any    ║
║    seed within reach.                                            ║
║                                                                  ║
║  Install:  pip install numpy scipy        (mpmath optional)      ║
║  Run:      python geometric_prime_counter_v4.py                  ║
║            python geometric_prime_counter_v4.py --controls       ║
║                                                                  ║
║  Theory:   Riemann (1859), Siegel (1932), Selberg (1956),        ║
║            Faddeev–Pavlov (1972), Lax–Phillips (1976),           ║
║            Colin de Verdière (1983)                              ║
╚══════════════════════════════════════════════════════════════════╝

CHANGES FROM v3 (v4, 6 September 2026)
  • The claim that surface eigenvalues approximate ζ zeros is withdrawn.
    PART 4 adds the control experiments (distinct-zero count against
    evenly spaced and random seeds; raw nearest-zero distance; Newton
    search-radius sweep; wall-height sweep with a significance estimate;
    quantisation roots of the wall condition).  The surface ties or loses
    to uninformed seeds at every radius.  A weak signal remains (raw
    candidates sit ~20% closer to zeros than random points, 1.4–1.9σ
    one-sided at Y = 20, 50, 100) — real, not usable.
  • The theory in build_surface() is corrected.  The Dirichlet wall's
    continuum eigenvalues obey  Y^{2it} = −φ(½+it), a near-uniform comb
    of density (log Y + log(t/π))/π.  The ζ zeros are not eigenvalues of
    this operator; they are poles of φ(s) — scattering resonances.  The
    pseudo-Laplacian (Colin de Verdière) imposes the same condition on
    the constant term and yields the same comb.
  • "60–70% of candidates fail to bracket" (v3 README) was wrong: 47 of
    48 converge and funnel onto 17 distinct zeros.  main() now reports
    funnelling, not failure.
  • Precision statement corrected.  Newton converges to 1e-12 on the
    zero of the ONE-TERM Riemann–Siegel approximation used in rs_Z(); that
    zero differs from the true γ by up to 8e-3 for γ < 70 (median 3e-4
    over the first 108).  The effect on π(10⁶) is about one count.  v1–v3
    said "machine precision"; that was true of the iteration, not of γ.
  • Boundary-condition remark.  build_surface() imposes Dirichlet on the
    arc |z| = 1 together with Neumann at x = 0.  The S-gluing of the
    modular surface gives Neumann on the arc for even forms and Dirichlet
    only for odd forms, so the discretised operator is not exactly the
    modular Laplacian.  This does not affect the controls (which test the
    seeds as numbers) or Control 3 (which uses the exact φ of SL(2,Z)\H²);
    it is recorded so the operator is described correctly.
  • Nothing about the prime-free pipeline itself changes.  Riemann's J
    still beats ψ(x)/log x by two orders of magnitude (v3 result).
"""

import numpy as np
from math import log, sqrt, pi, isqrt, exp
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.special import exp1, expi
import time
import sys

# ════════════════════════════════════════════════════════════════
# PART 1 — TRUNCATED HYPERBOLIC LAPLACIAN  (seed generator)
# Discretise SL(2,Z)\H² truncated at height Y (see boundary-condition
# remark in build_surface).  Eigenvalues λ = ¼ + t² → seeds γ = 2t.
# ════════════════════════════════════════════════════════════════

def build_surface(Nx=50, Nu=300, Y=50.0):
    r"""
    Discretise the hyperbolic Laplacian on the truncated modular surface.

    Coordinates: x ∈ (0, ½),  u = log y,  y = e^u.
    Boundary conditions:
      - Dirichlet on the arc  x² + e^{2u} = 1  (lower boundary)
      - Dirichlet at  u = log Y                 (cusp wall)
      - Neumann at  x = 0  and  x = ½           (x ↦ −x symmetry; T-periodicity)

    BOUNDARY-CONDITION REMARK (v4).  Neumann at x = 0 selects functions
    even under x ↦ −x.  The inversion S: z ↦ −1/z maps the arc to itself
    by (x, y) ↦ (−x, y), so on the modular surface an even eigenfunction
    satisfies a NEUMANN condition on the arc and only an odd one satisfies
    Dirichlet.  The Dirichlet-arc/Neumann-x combination used here is
    therefore not the modular-surface gluing; the interior modes it
    produces are not Maass cusp forms of SL(2,Z).  The cusp-end analysis
    below (the comb) is unchanged in form, but the phase wobble on it is
    the scattering phase of THIS mixed problem, not exactly arg φ.  None
    of the controls in PART 4 depend on this: Controls 1, 2, 2b test the
    seeds as numbers, and Control 3 uses the exact φ of SL(2,Z)\H².

    WHAT THIS SPECTRUM IS (corrected, v4).
    The wall discretises the cusp's continuous spectrum.  For the constant
    Fourier term  y^s + φ(s) y^{1−s}  (φ = Λ(2s−1)/Λ(2s) the scattering
    matrix of SL(2,Z)\H²) the Dirichlet condition at y=Y reads
                    Y^{2it} = −φ(½+it),
    i.e.  2t·log Y = arg(−φ) + 2πn.  The 2t·log Y term dominates, so the
    continuum eigenvalues form a near-uniform comb in t with density
    (log Y + log(t/π))/π; the ζ zeros enter only as a phase wobble.  Those
    comb roots are no closer to ζ zeros than random points (PART 4).  The
    remaining eigenvalues are interior modes (Weyl count ~t²/24 on this
    half domain; on the true surface these would be Maass cusp forms),
    which have nothing to do with ζ.  The ζ zeros themselves are
    the POLES of φ(s), at s = ρ/2 — resonances, not L² eigenvalues of any
    self-adjoint truncation (Faddeev–Pavlov 1972; Lax–Phillips 1976).
    Colin de Verdière's pseudo-Laplacian imposes the same constant-term
    condition and produces the same comb; it is not a fix.

    The returned 'approx_zeros' = 2t are therefore SEEDS spanning an
    interval, not approximations to zeros.  They are kept so the control
    experiments in PART 4 can be reproduced.

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
# PART 2 — NEWTON REFINEMENT (Riemann–Siegel Z function)
# Locate the zero of Z(t) nearest each seed.  Precision: 1e-12 on the
# one-term Riemann–Siegel approximation below, which places the zero
# within ~1e-2 (γ ≈ 25) to ~1e-4 (γ ≈ 250) of the true γ.
# ════════════════════════════════════════════════════════════════

def _rs_theta(t):
    """Riemann-Siegel phase θ(t).  Accurate for t > 10."""
    return (t/2*log(t/(2*pi)) - t/2 - pi/8
            + 1/(48*t) + 7/(5760*t**3))

def rs_Z(t):
    """
    Real-valued Z(t) = e^{iθ(t)} ζ(½+it) by the Riemann–Siegel formula with
    the first remainder term C₀ only.  Its zeros approximate the γ:
    measured against mpmath.zetazero, the located zeros differ from the
    true γ by at most 7.5e-3 for γ < 70 (worst case near γ = 25.01) and by
    a median 2.6e-4 over the first 108 zeros.  That is the precision floor
    of everything downstream; the Newton iteration in refine_zero()
    converges to 1e-12 on THIS function, not on ζ.
    """
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
    Find the zero of rs_Z nearest to t0 (bracket → bisect → Newton).
    Returns refined t, or None if no sign change lies within search_radius.
    Precision is set by rs_Z (one-term Riemann–Siegel); see its docstring.
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
# PART 2b — GLOBAL ZERO SCAN (optional, prime-free)
# Walk Z(t) along the critical line and collect every sign change.
# Used to show how accuracy improves with zero count.
# ════════════════════════════════════════════════════════════════

def scan_zeros(t_max=250.0, step=0.05):
    """Find all zeros of Z(t) in (10, t_max] by sign-change scanning + refinement."""
    zeros = []
    t, Zp = 10.0, rs_Z(10.0)
    while t < t_max:
        t2 = t + step; Z2 = rs_Z(t2)
        if Zp * Z2 < 0:
            r = refine_zero((t + t2)/2, search_radius=2*step, step=step/2)
            if r is not None: zeros.append(r)
        t, Zp = t2, Z2
    return sorted(set(round(z, 9) for z in zeros))


# ════════════════════════════════════════════════════════════════
# PART 3 — PRIME COUNTING (Riemann's explicit formula, 1859)
#
#   J(x) = li(x) − Σ_ρ li(x^ρ) − log 2 + ∫_x^∞ dt / (t(t²−1) log t)
#   π(x) = J(x) − π(x^{1/2})/2 − π(x^{1/3})/3 − …      (recursive)
#
# where J(x) = Σ_{p^k ≤ x} 1/k counts prime powers with weight 1/k.
# Everything here is li, Ei and E₁ at real or complex arguments:
# no primes, no factorisation, no Möbius function.
# ════════════════════════════════════════════════════════════════

def riemann_J(x, zeros):
    """
    Riemann's prime-power counting function J(x) from imaginary parts of ζ zeros.

    li(x^ρ) = Ei(ρ log x); pairing ρ with its conjugate gives 2 Re Ei(ρ log x),
    and Re Ei(z) = −Re E₁(−z), which scipy evaluates for complex z.
    The tail integral equals Σ_{k≥1} E₁(2k log x).
    """
    L   = log(x)
    rho = 0.5 + 1j*np.asarray(zeros, dtype=float)
    osc = -2.0 * np.real(exp1(-rho*L)).sum()
    k   = np.arange(1, 40)
    tail = exp1(2*k*L).sum()
    return expi(L) - log(2) - osc + tail

def count_primes(x, zeros, max_zeros=None):
    """
    Estimate π(x) from imaginary parts of ζ zeros via Riemann's explicit formula.

    π(x) = J(x) − Σ_{k≥2} π(x^{1/k}) / k, applied recursively until x^{1/k} < 2.
    Equivalent to Möbius inversion but needs no arithmetic function.

    Parameters
    ----------
    x        : positive real, the bound to count primes up to
    zeros    : iterable of imaginary parts γ of nontrivial ζ zeros
    max_zeros: use only the first max_zeros zeros
    """
    if x < 2:
        return 0.0
    gammas = list(zeros)[:max_zeros] if max_zeros else list(zeros)
    val = riemann_J(x, gammas)
    k = 2
    while x**(1.0/k) >= 2:
        val -= count_primes(x**(1.0/k), gammas) / k
        k += 1
    return max(0.0, val)

def count_primes_psi(x, zeros):
    """
    The v2 conversion, kept for comparison:  π(x) ≈ ψ(x)/log x  with
    ψ(x) = x − log 2π − Σ_γ 2√x [½cos(γ log x) + γ sin(γ log x)]/(¼+γ²).

    This discards the x/log²x term of π(x) and reduces to the prime number
    theorem's leading order: the error is ≈ 1/log x whatever zeros are supplied.
    """
    if x < 2:
        return 0.0
    psi = x - log(2*pi); lx = log(x); sqx = sqrt(x)
    for g in zeros:
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
# PART 4 — CONTROL EXPERIMENTS  (new in v4)
# Does the surface locate ζ zeros better than seeds that know nothing
# about the surface?  Reference zeros come from scan_zeros() — the
# prime-free Riemann–Siegel sweep — so nothing here needs a table.
# ════════════════════════════════════════════════════════════════

def _distinct(seeds, search_radius=2.0, step=None):
    """Refine every seed; return the sorted list of DISTINCT zeros hit."""
    if step is None: step = min(0.1, search_radius/8)
    hits = [refine_zero(s, search_radius=search_radius, step=step) for s in seeds]
    hits = sorted(set(round(h, 9) for h in hits if h))
    out  = [hits[0]] if hits else []
    for h in hits[1:]:
        if h - out[-1] > 0.5: out.append(h)
    return out

def nearest_zero_distances(seeds, true_zeros):
    z = np.asarray(true_zeros)
    return np.array([np.min(np.abs(z - s)) for s in seeds])

def control_table(cand, true_zeros, n_random=20, seed=0, search_radius=2.0):
    """
    Three rows: surface candidates / evenly spaced / uniform random, all
    with the same count and span.  Reports distinct zeros recovered and the
    raw distance of each seed to the nearest true zero (no Newton).
    """
    cand = np.sort(np.asarray(cand)); lo, hi, n = cand.min(), cand.max(), len(cand)
    rng  = np.random.default_rng(seed)
    even = np.linspace(lo, hi, n)
    rows = []
    for name, seeds in [('surface eigenvalues', cand), ('evenly spaced, same span', even)]:
        d = nearest_zero_distances(seeds, true_zeros)
        rows.append((name, len(_distinct(seeds, search_radius)), np.median(d), d.mean()))
    dist, meds, means = [], [], []
    for _ in range(n_random):
        r = rng.uniform(lo, hi, n); d = nearest_zero_distances(r, true_zeros)
        dist.append(len(_distinct(r, search_radius))); meds.append(np.median(d)); means.append(d.mean())
    rows.append((f'uniform random, same span (×{n_random})', np.mean(dist), np.mean(meds), np.mean(means)))
    nz  = np.sum((np.asarray(true_zeros) >= lo) & (np.asarray(true_zeros) <= hi))
    gap = (hi - lo) / max(nz, 1)
    return {'rows': rows, 'span': (lo, hi), 'n': n, 'true_in_span': int(nz), 'gap_over_4': gap/4}

def radius_sweep(cand, radii=(0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0), n_random=6, seed=1):
    """Distinct zeros recovered vs Newton search radius, surface / even / random."""
    cand = np.sort(np.asarray(cand)); lo, hi, n = cand.min(), cand.max(), len(cand)
    even = np.linspace(lo, hi, n); rng = np.random.default_rng(seed)
    out = []
    for R in radii:
        rs = [len(_distinct(rng.uniform(lo, hi, n), R)) for _ in range(n_random)]
        out.append((R, len(_distinct(cand, R)), len(_distinct(even, R)), float(np.mean(rs))))
    return out

def proximity_significance(cand, true_zeros, n_draws=3000, seed=11):
    """
    How unusual is the surface seeds' median nearest-zero distance?
    Draw n_draws sets of len(cand) uniform points on the same span, take
    each set's median nearest-zero distance, and report the surface's
    median against that distribution (z-score and one-sided fraction of
    random sets doing at least as well).  Raw distances only; no Newton.
    """
    cand = np.sort(np.asarray(cand)); lo, hi, n = cand.min(), cand.max(), len(cand)
    rng  = np.random.default_rng(seed)
    obs  = float(np.median(nearest_zero_distances(cand, true_zeros)))
    meds = np.array([np.median(nearest_zero_distances(rng.uniform(lo, hi, n), true_zeros))
                     for _ in range(n_draws)])
    return {'surface_median': obs, 'random_median_mean': float(meds.mean()),
            'random_median_sd': float(meds.std()),
            'z': float((meds.mean() - obs) / meds.std()),
            'p_one_sided': float((meds <= obs).mean()), 'n': n, 'span': (lo, hi)}

def wall_height_sweep(true_zeros, heights=(20.0, 50.0, 100.0), Nx=50, Nu=300, n_draws=3000):
    """Repeat Control 1's raw-distance test and proximity_significance at several Y."""
    out = []
    for Y in heights:
        s  = build_surface(Nx=Nx, Nu=Nu, Y=Y)
        ps = proximity_significance(s['approx_zeros'], true_zeros, n_draws=n_draws)
        nd = len(_distinct(s['approx_zeros'], 2.0))
        ne = len(_distinct(np.linspace(*ps['span'], ps['n']), 2.0))
        out.append((Y, ps['n'], ps['span'], nd, ne, ps))
    return out

def quantization_roots(Y, t_lo=4.5, t_hi=37.0, n_grid=6500):
    """
    Roots of the Dirichlet-wall constant-term condition  Y^{2it} = −φ(½+it),
    φ(s) = Λ(2s−1)/Λ(2s).  Requires mpmath (ζ is used here for DIAGNOSIS of
    the operator only; it plays no part in the prime-counting pipeline).
    Returns the roots in t and the mean spacing.
    """
    try:
        import mpmath as mp
    except ImportError:
        return None
    mp.mp.dps = 20
    lY = log(Y)
    Lam = lambda s: mp.pi**(-s/2) * mp.gamma(s/2) * mp.zeta(s)
    phi = lambda t: Lam(2*mp.mpc(0.5, t) - 1) / Lam(2*mp.mpc(0.5, t))
    ts = np.linspace(t_lo, t_hi, n_grid); F = np.zeros_like(ts); prev = None; acc = 0.0
    for i, t in enumerate(ts):
        a = float(mp.arg(-phi(t)))
        if prev is not None:
            d = a - prev
            if d >  np.pi: acc -= 2*np.pi
            if d < -np.pi: acc += 2*np.pi
        prev = a; F[i] = 2*t*lY - (a + acc)
    roots = []
    for k in range(int(F[0]/(2*np.pi)) - 1, int(F[-1]/(2*np.pi)) + 2):
        tgt = 2*np.pi*k; idx = np.where(np.diff(np.sign(F - tgt)))[0]
        for i in idx:
            roots.append(ts[i] + (tgt - F[i])*(ts[i+1]-ts[i])/(F[i+1]-F[i]))
    roots = np.array(sorted(roots))
    return {'roots': roots, 'mean_spacing': float(np.mean(np.diff(roots))),
            'predicted_spacing_pi_over_logY': pi/lY}

def run_controls(surf, true_zeros, full=False):
    cand = surf['approx_zeros']
    print("  Control 1 — do the surface eigenvalues beat uninformed seeds?")
    ct = control_table(cand, true_zeros)
    lo, hi = ct['span']
    print(f"  {ct['n']} seeds spanning γ = {lo:.2f} … {hi:.2f};  "
          f"{ct['true_in_span']} true zeros in span;  gap/4 = {ct['gap_over_4']:.3f}")
    print(f"  {'seeds':<36} {'distinct zeros':>14} {'median dist':>12} {'mean dist':>10}")
    print("  " + "─"*76)
    for name, nd, md, mn in ct['rows']:
        print(f"  {name:<36} {nd:>14.1f} {md:>12.3f} {mn:>10.3f}")
    print("  A random point sits a median gap/4 from the nearest zero; compare the surface row.")
    print()
    if not full:
        print("  (run with --controls for the radius sweep and the quantization-root test)")
        return
    print("  Control 2 — Newton search radius sweep (distinct zeros recovered)")
    print(f"  {'radius':>8} {'surface':>9} {'even':>7} {'random':>8}")
    for R, s, e, r in radius_sweep(cand):
        print(f"  {R:>8.2f} {s:>9d} {e:>7d} {r:>8.1f}")
    print()
    print("  Control 2b — wall-height sweep: raw proximity signal and its significance")
    print(f"  {'Y':>5} {'n':>3} {'span':>15} {'distinct':>9} {'even':>5} {'surf med':>9} {'rand med':>9} {'z':>5} {'p(1-sided)':>10}")
    for Y, n, (lo, hi), nd, ne, ps in wall_height_sweep(true_zeros):
        print(f"  {Y:>5.0f} {n:>3d} {lo:>6.2f} … {hi:>6.2f} {nd:>9d} {ne:>5d} "
              f"{ps['surface_median']:>9.3f} {ps['random_median_mean']:>9.3f} {ps['z']:>5.2f} {ps['p_one_sided']:>10.3f}")
    print("  'rand med' = mean over 3000 random same-span seed sets of their median nearest-zero distance;")
    print("  z and p compare the surface row to that distribution.  A weak signal, present at every Y.")
    print()
    print("  Control 3 — what the Dirichlet wall actually quantises  (needs mpmath)")
    q = quantization_roots(50.0)
    if q is None:
        print("  mpmath not installed; skipped."); return
    tz = np.asarray(true_zeros)/2.0
    dz = nearest_zero_distances(q['roots'], tz)
    inside = tz[(tz > 5) & (tz < 36)]
    print(f"  roots of Y^(2it) = −φ(½+it) on t∈[5,36]:  {len(q['roots'])}")
    print(f"  mean spacing {q['mean_spacing']:.3f}   "
          f"(π/log Y = {q['predicted_spacing_pi_over_logY']:.3f};  with the phase term "
          f"(log Y + log(t/π))/π gives ≈ {1/((log(50)+log(20/pi))/pi):.3f})")
    print(f"  median distance root → nearest ζ-zero/2: {np.median(dz):.3f}   "
          f"vs gap/4 for random points: {np.mean(np.diff(inside))/4:.3f}")
    print("  The comb roots are no closer to ζ zeros than random points: the wall does not target zeros.")
    print()


# ════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════

def main(controls=False):
    print()
    print("╔" + "═"*60 + "╗")
    print("║" + "  PRIME-FREE PRIME COUNTER  v4".center(60) + "║")
    print("║" + "  Riemann–Siegel → ζ zeros → Riemann's J → π(x)".center(60) + "║")
    print("║" + "  with controls on the surface's contribution".center(60) + "║")
    print("╚" + "═"*60 + "╝")
    print()

    # ── Step 1: surface (kept for the controls) ──────────────────
    print("  Step 1 — Truncated hyperbolic surface  (seed generator; see PART 4)")
    surf = build_surface(Nx=50, Nu=300, Y=50.0)
    print(f"  Done in {surf['build_time_s']*1000:.0f}ms  ({surf['n_points']} grid points, "
          f"{len(surf['approx_zeros'])} eigenvalue seeds)")
    print()

    # ── Step 2: Newton on Riemann–Siegel ─────────────────────────
    print("  Step 2 — Newton refinement on Riemann–Siegel Z(t)")
    t0 = time.perf_counter()
    per_seed = [refine_zero(c) for c in surf['approx_zeros']]
    refined  = refine_all(surf['approx_zeros'], verbose=False)
    n_cand, n_conv, n_dist = len(per_seed), sum(1 for r in per_seed if r), len(refined)
    print(f"  Done in {(time.perf_counter()-t0)*1000:.0f}ms  "
          f"({n_conv}/{n_cand} seeds converge, funnelling onto {n_dist} distinct zeros)")
    print()

    # ── Step 2b: prime-free global scan ──────────────────────────
    print("  Step 2b — Global Riemann–Siegel scan to t = 250  (this is where the zeros really come from)")
    t0 = time.perf_counter(); scanned = scan_zeros(250.0)
    print(f"  Done in {(time.perf_counter()-t0)*1000:.0f}ms  ({len(scanned)} zeros)")
    print()

    # ── Step 3: prime counting ───────────────────────────────────
    print("  Step 3 — Prime counting")
    test_x = [100, 1_000, 10_000, 100_000, 1_000_000]
    print(f"  {'x':>10}  {'π(x) exact':>11}  {'ψ/log x':>9}  {'Riemann':>9}  {'err':>7}  {'Riemann':>9}  {'err':>7}")
    print(f"  {'':>10}  {'':>11}  {'(17 z)':>9}  {'(17 z)':>9}  {'':>7}  {f'({len(scanned)} z)':>9}  {'':>7}")
    print("  " + "─"*74)
    for x in test_x:
        exact = sieve_exact(x)
        old, r17, rall = count_primes_psi(x, refined), count_primes(x, refined), count_primes(x, scanned)
        print(f"  {x:>10,}  {exact:>11,}  {old:>9.0f}  {r17:>9.1f}  {abs(r17-exact)/exact*100:>6.2f}%  "
              f"{rall:>9.1f}  {abs(rall-exact)/exact*100:>6.2f}%")
    print()

    # ── Step 4: controls ─────────────────────────────────────────
    print("  Step 4 — Controls: did the surface contribute?")
    run_controls(surf, scanned, full=controls)

    print("  Summary")
    print("  • The pipeline is prime-free at every step and Riemann's J recovers π(x) to ~0.02%.")
    print("  • The zeros are found by Riemann–Siegel from any seed within reach of a zero,")
    print("    to ~1e-3 in γ (one-term Riemann–Siegel remainder): about one count at x = 10⁶.")
    print("  • The surface supplies an interval to search, plus a weak (~20%) proximity signal.")
    print("  • It does not supply the zeros.  See PART 4 and README §4.3, §5.")
    print()


if __name__ == "__main__":
    main(controls=('--controls' in sys.argv))
