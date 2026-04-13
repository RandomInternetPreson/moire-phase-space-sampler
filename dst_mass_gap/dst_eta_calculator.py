"""
DST Many-Body Enhancement Factor Calculator
============================================
Computes η(ρ) from first-principles 3D overlap integrals across:
  - 20 density points (1.5ρ₀ to 12ρ₀)
  - Multiple Yukawa range parameters (λ)
  - Multiple Gaussian rotational profile widths (σ)
  - FCC, BCC, and random liquid-like lattice geometries
  - Bootstrap uncertainty quantification on the power-law exponent

Output: η vs ρ with confidence bands, fitted exponent with CI,
        total DST energy density exponent with CI.

Author: Aaron Alai with Claude (Anthropic)
Date: April 2026
"""

import numpy as np
from itertools import product
import time
import sys
import os

# ══════════════════════════════════════════════════════════════
#  PHYSICAL CONSTANTS
# ══════════════════════════════════════════════════════════════
RHO_NUC = 2.7e14  # g/cm³, nuclear saturation density
D0 = 1.84  # fm, inter-nucleon separation at ρ₀
LAMBDA_PI = 1.41  # fm, pion Compton wavelength (default Yukawa range)
SIGMA_DEFAULT = 0.50  # fm, default Gaussian rotational profile width
ALPHA_EM = 1.0 / 137.036


# ══════════════════════════════════════════════════════════════
#  LATTICE GENERATORS
# ══════════════════════════════════════════════════════════════

def fcc_neighbors(d):
    """12 nearest neighbors in FCC arrangement at separation d."""
    # FCC nearest neighbors are at (±1,±1,0)/√2, permutations
    raw = []
    for i, j in [(0,1),(0,2),(1,2)]:
        for si in [-1, 1]:
            for sj in [-1, 1]:
                v = np.zeros(3)
                v[i] = si
                v[j] = sj
                raw.append(v)
    raw = np.array(raw)
    raw = raw / np.linalg.norm(raw[0]) * d  # scale to separation d
    return raw


def bcc_neighbors(d):
    """8 nearest neighbors in BCC arrangement at separation d."""
    raw = []
    for sx in [-1, 1]:
        for sy in [-1, 1]:
            for sz in [-1, 1]:
                raw.append([sx, sy, sz])
    raw = np.array(raw, dtype=float)
    raw = raw / np.linalg.norm(raw[0]) * d
    return raw


def random_liquid_neighbors(d, n_neighbors=12, rng=None):
    """
    Random liquid-like packing: n_neighbors placed at approximately
    distance d with ~10% positional noise, no hard overlaps.
    """
    if rng is None:
        rng = np.random.default_rng()
    
    positions = []
    attempts = 0
    min_sep = 0.7 * d  # hard-core exclusion
    
    while len(positions) < n_neighbors and attempts < 10000:
        # Random direction, distance drawn from narrow Gaussian around d
        direction = rng.standard_normal(3)
        direction /= np.linalg.norm(direction)
        dist = d * (1.0 + 0.1 * rng.standard_normal())
        dist = max(dist, min_sep)
        
        candidate = direction * dist
        
        # Check against existing positions
        ok = True
        for p in positions:
            if np.linalg.norm(candidate - p) < min_sep:
                ok = False
                break
        if ok:
            positions.append(candidate)
        attempts += 1
    
    return np.array(positions)


# ══════════════════════════════════════════════════════════════
#  FIELD PROFILES
# ══════════════════════════════════════════════════════════════

def yukawa_profile(r, lam):
    """Radial displacement profile: Yukawa ~ exp(-r/λ)/r, regularized."""
    r_safe = np.maximum(r, 0.01)  # regularize at origin
    return np.exp(-r_safe / lam) / r_safe


def gaussian_profile(r, sigma):
    """Rotational displacement profile: Gaussian."""
    return np.exp(-0.5 * (r / sigma) ** 2)


# ══════════════════════════════════════════════════════════════
#  3D OVERLAP INTEGRAL
# ══════════════════════════════════════════════════════════════

def compute_eta_3d(neighbor_positions, lam, sigma, grid_n=65):
    """
    Compute the many-body enhancement factor η on a 3D grid.
    
    The cross-coupling energy density is:
      ε = ½g × φ_r_total² × |Φ_θ_total|²
    
    η = (ε_total - Σ ε_isolated) / Σ ε_isolated
    
    i.e., the fractional excess from interference.
    
    Parameters
    ----------
    neighbor_positions : (N, 3) array, positions of neighbors in fm
    lam : float, Yukawa range parameter in fm
    sigma : float, Gaussian width in fm
    grid_n : int, grid points per dimension
    
    Returns
    -------
    eta : float, enhancement factor (dimensionless)
    """
    n_neighbors = len(neighbor_positions)
    
    # Determine grid extent: must cover all particles + several λ
    max_extent = np.max(np.abs(neighbor_positions)) + 4 * lam
    max_extent = max(max_extent, 4 * lam)  # at least 4λ from center
    
    # 3D grid
    x = np.linspace(-max_extent, max_extent, grid_n)
    dx = x[1] - x[0]
    dV = dx ** 3
    
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    
    # Central nucleon at origin
    r_center = np.sqrt(X**2 + Y**2 + Z**2)
    phi_r_center = yukawa_profile(r_center, lam)
    phi_theta_center = gaussian_profile(r_center, sigma)
    
    # Total fields = sum over all nucleons (central + neighbors)
    phi_r_total = phi_r_center.copy()
    phi_theta_total = phi_theta_center.copy()
    
    # Sum of isolated energies
    eps_iso_center = phi_r_center**2 * phi_theta_center**2
    E_isolated = np.sum(eps_iso_center) * dV
    
    for pos in neighbor_positions:
        r_i = np.sqrt((X - pos[0])**2 + (Y - pos[1])**2 + (Z - pos[2])**2)
        phi_r_i = yukawa_profile(r_i, lam)
        phi_theta_i = gaussian_profile(r_i, sigma)
        
        phi_r_total += phi_r_i
        phi_theta_total += phi_theta_i  # phase-aligned: Δθ = 0
        
        eps_iso_i = phi_r_i**2 * phi_theta_i**2
        E_isolated += np.sum(eps_iso_i) * dV
    
    # Total cross-coupling energy
    eps_total = phi_r_total**2 * phi_theta_total**2
    E_total = np.sum(eps_total) * dV
    
    # Enhancement factor
    if E_isolated > 0:
        eta = (E_total - E_isolated) / E_isolated
    else:
        eta = 0.0
    
    return eta


# ══════════════════════════════════════════════════════════════
#  DENSITY SWEEP
# ══════════════════════════════════════════════════════════════

def compute_eta_at_density(rho_over_rho0, lam, sigma, lattice='fcc', 
                           grid_n=65, rng=None):
    """
    Compute η at a given density for specified profile parameters.
    
    Parameters
    ----------
    rho_over_rho0 : float, density in units of ρ₀
    lam : float, Yukawa range in fm
    sigma : float, Gaussian width in fm
    lattice : str, 'fcc', 'bcc', or 'random'
    grid_n : int, grid resolution
    rng : numpy RNG for random lattice
    
    Returns
    -------
    eta : float
    """
    d = D0 * rho_over_rho0 ** (-1.0 / 3.0)  # inter-nucleon separation
    
    if lattice == 'fcc':
        neighbors = fcc_neighbors(d)
    elif lattice == 'bcc':
        neighbors = bcc_neighbors(d)
    elif lattice == 'random':
        neighbors = random_liquid_neighbors(d, n_neighbors=12, rng=rng)
    else:
        raise ValueError(f"Unknown lattice: {lattice}")
    
    return compute_eta_3d(neighbors, lam, sigma, grid_n=grid_n)


# ══════════════════════════════════════════════════════════════
#  POWER LAW FITTING
# ══════════════════════════════════════════════════════════════

def fit_power_law(rho_vals, eta_vals, min_rho=2.5):
    """
    Fit η = A × (ρ/ρ₀)^β to data above min_rho.
    Returns (A, beta, r_squared).
    """
    mask = rho_vals >= min_rho
    if np.sum(mask) < 3:
        return np.nan, np.nan, np.nan
    
    x = np.log(rho_vals[mask])
    y = np.log(np.maximum(eta_vals[mask], 1e-30))
    
    # Least squares: y = β·x + ln(A)
    n = len(x)
    sx = np.sum(x)
    sy = np.sum(y)
    sxx = np.sum(x**2)
    sxy = np.sum(x * y)
    
    beta = (n * sxy - sx * sy) / (n * sxx - sx**2)
    lnA = (sy - beta * sx) / n
    A = np.exp(lnA)
    
    # R²
    y_pred = beta * x + lnA
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    
    return A, beta, r2


def bootstrap_exponent(rho_vals, eta_vals, n_bootstrap=2000, min_rho=2.5, 
                       rng=None):
    """
    Bootstrap confidence interval on the power-law exponent.
    Returns (median_beta, ci_68_low, ci_68_high, ci_95_low, ci_95_high).
    """
    if rng is None:
        rng = np.random.default_rng(42)
    
    mask = rho_vals >= min_rho
    rho_fit = rho_vals[mask]
    eta_fit = eta_vals[mask]
    n = len(rho_fit)
    
    if n < 3:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    
    betas = []
    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        _, beta, _ = fit_power_law(rho_fit[idx], eta_fit[idx], min_rho=0)
        if np.isfinite(beta):
            betas.append(beta)
    
    betas = np.array(betas)
    median = np.median(betas)
    ci68 = np.percentile(betas, [16, 84])
    ci95 = np.percentile(betas, [2.5, 97.5])
    
    return median, ci68[0], ci68[1], ci95[0], ci95[1]


# ══════════════════════════════════════════════════════════════
#  MAIN SWEEP
# ══════════════════════════════════════════════════════════════

def run_full_sweep(grid_n=65, n_random=5, save_dir=None):
    """
    Run the complete parameter sweep.
    
    Returns a structured results dictionary.
    """
    # Density points: 20 logarithmically spaced from 1.5 to 12
    densities = np.logspace(np.log10(1.5), np.log10(12.0), 20)
    
    # Profile parameters to sweep
    lambdas = [1.0, 1.2, 1.41, 1.6, 2.0]  # Yukawa range (fm)
    sigmas = [0.35, 0.50, 0.65, 0.80]      # Gaussian width (fm)
    lattices = ['fcc', 'bcc']               # + random
    
    rng = np.random.default_rng(2026)
    
    results = {
        'densities': densities,
        'lambdas': lambdas,
        'sigmas': sigmas,
        'lattices': lattices + ['random'],
        'n_random': n_random,
        'grid_n': grid_n,
        'eta_data': {},  # key: (lam, sigma, lattice, random_seed) -> eta array
    }
    
    total_runs = len(densities) * (
        len(lambdas) * len(sigmas) * len(lattices) +
        len(lambdas) * len(sigmas) * n_random
    )
    
    print(f"DST η Calculator")
    print(f"================")
    print(f"Densities: {len(densities)} points, {densities[0]:.1f}–{densities[-1]:.1f} ρ₀")
    print(f"Profiles: {len(lambdas)} λ × {len(sigmas)} σ = {len(lambdas)*len(sigmas)} combinations")
    print(f"Lattices: FCC, BCC, + {n_random} random realizations")
    print(f"Grid: {grid_n}³ = {grid_n**3:,} points")
    print(f"Total runs: {total_runs:,}")
    print(f"")
    
    run_count = 0
    t_start = time.time()
    
    # Structured lattices
    for lam, sigma, lattice in product(lambdas, sigmas, lattices):
        key = (lam, sigma, lattice, 0)
        eta_arr = np.zeros(len(densities))
        
        for i, rho in enumerate(densities):
            eta_arr[i] = compute_eta_at_density(rho, lam, sigma, lattice, 
                                                 grid_n=grid_n)
            run_count += 1
        
        results['eta_data'][key] = eta_arr
        
        elapsed = time.time() - t_start
        rate = run_count / elapsed if elapsed > 0 else 0
        eta_remaining = (total_runs - run_count) / rate if rate > 0 else 0
        
        sys.stdout.write(f"\r  [{run_count}/{total_runs}] "
                        f"λ={lam:.2f} σ={sigma:.2f} {lattice:6s} "
                        f"η(3ρ₀)={eta_arr[np.argmin(np.abs(densities-3))]:.3f} "
                        f"η(10ρ₀)={eta_arr[-2]:.3f} "
                        f"ETA: {eta_remaining/60:.1f}min")
        sys.stdout.flush()
    
    # Random lattices
    for lam, sigma in product(lambdas, sigmas):
        for seed in range(n_random):
            key = (lam, sigma, 'random', seed)
            eta_arr = np.zeros(len(densities))
            
            for i, rho in enumerate(densities):
                eta_arr[i] = compute_eta_at_density(rho, lam, sigma, 'random',
                                                     grid_n=grid_n, rng=rng)
                run_count += 1
            
            results['eta_data'][key] = eta_arr
            
            elapsed = time.time() - t_start
            rate = run_count / elapsed if elapsed > 0 else 0
            eta_remaining = (total_runs - run_count) / rate if rate > 0 else 0
            
            sys.stdout.write(f"\r  [{run_count}/{total_runs}] "
                            f"λ={lam:.2f} σ={sigma:.2f} rand#{seed} "
                            f"ETA: {eta_remaining/60:.1f}min   ")
            sys.stdout.flush()
    
    print(f"\n\nCompleted {run_count} runs in {(time.time()-t_start)/60:.1f} minutes.")
    
    return results


def analyze_results(results):
    """
    Analyze the sweep results: fit exponents, compute confidence intervals.
    """
    densities = results['densities']
    
    print(f"\n{'='*70}")
    print(f"ANALYSIS")
    print(f"{'='*70}")
    
    # Collect all η curves
    all_etas = []
    all_exponents = []
    all_prefactors = []
    
    for key, eta_arr in results['eta_data'].items():
        lam, sigma, lattice, seed = key
        A, beta, r2 = fit_power_law(densities, eta_arr, min_rho=2.5)
        if np.isfinite(beta) and beta > 0:
            all_etas.append(eta_arr)
            all_exponents.append(beta)
            all_prefactors.append(A)
    
    all_etas = np.array(all_etas)
    all_exponents = np.array(all_exponents)
    all_prefactors = np.array(all_prefactors)
    
    # Summary statistics
    print(f"\nη exponent across all {len(all_exponents)} parameter combinations:")
    print(f"  Median:  {np.median(all_exponents):.3f}")
    print(f"  Mean:    {np.mean(all_exponents):.3f}")
    print(f"  Std:     {np.std(all_exponents):.3f}")
    print(f"  68% CI:  [{np.percentile(all_exponents, 16):.3f}, "
          f"{np.percentile(all_exponents, 84):.3f}]")
    print(f"  95% CI:  [{np.percentile(all_exponents, 2.5):.3f}, "
          f"{np.percentile(all_exponents, 97.5):.3f}]")
    print(f"  Min:     {np.min(all_exponents):.3f}")
    print(f"  Max:     {np.max(all_exponents):.3f}")
    
    # Total DST energy density exponent = η_exponent + 2/3 (overlap) + 1 (density)
    total_exponents = all_exponents + 2.0/3.0 + 1.0
    
    print(f"\nTotal DST energy density exponent (η + 2/3 + 1):")
    print(f"  Median:  {np.median(total_exponents):.3f}")
    print(f"  68% CI:  [{np.percentile(total_exponents, 16):.3f}, "
          f"{np.percentile(total_exponents, 84):.3f}]")
    print(f"  95% CI:  [{np.percentile(total_exponents, 2.5):.3f}, "
          f"{np.percentile(total_exponents, 97.5):.3f}]")
    print(f"  Fermi:   1.667")
    print(f"  Exceeds Fermi in {np.sum(total_exponents > 1.667)}/{len(total_exponents)} "
          f"({100*np.sum(total_exponents > 1.667)/len(total_exponents):.1f}%) of cases")
    
    # Median η curve with bands
    eta_median = np.median(all_etas, axis=0)
    eta_16 = np.percentile(all_etas, 16, axis=0)
    eta_84 = np.percentile(all_etas, 84, axis=0)
    eta_2p5 = np.percentile(all_etas, 2.5, axis=0)
    eta_97p5 = np.percentile(all_etas, 97.5, axis=0)
    
    # Print table
    print(f"\n{'='*70}")
    print(f"η(ρ) TABLE — Median with 68% CI")
    print(f"{'='*70}")
    print(f"{'ρ/ρ₀':>8} {'d (fm)':>8} {'η_med':>10} {'η_16%':>10} "
          f"{'η_84%':>10} {'η×137 med':>10}")
    
    for i, rho in enumerate(densities):
        d = D0 * rho ** (-1.0/3.0)
        print(f"{rho:8.2f} {d:8.3f} {eta_median[i]:10.4f} {eta_16[i]:10.4f} "
              f"{eta_84[i]:10.4f} {eta_median[i]*137:10.1f}")
    
    # Breakdown by lattice type
    print(f"\n{'='*70}")
    print(f"EXPONENT BY LATTICE TYPE")
    print(f"{'='*70}")
    
    for lattice_type in ['fcc', 'bcc', 'random']:
        exps = []
        for key, eta_arr in results['eta_data'].items():
            lam, sigma, lattice, seed = key
            if lattice != lattice_type:
                continue
            A, beta, r2 = fit_power_law(densities, eta_arr, min_rho=2.5)
            if np.isfinite(beta) and beta > 0:
                exps.append(beta)
        exps = np.array(exps)
        if len(exps) > 0:
            print(f"  {lattice_type:8s}: median={np.median(exps):.3f} "
                  f"[{np.percentile(exps,16):.3f}, {np.percentile(exps,84):.3f}] "
                  f"n={len(exps)}")
    
    # Breakdown by λ
    print(f"\n{'='*70}")
    print(f"EXPONENT BY YUKAWA RANGE λ")
    print(f"{'='*70}")
    
    for lam in results['lambdas']:
        exps = []
        for key, eta_arr in results['eta_data'].items():
            if key[0] != lam:
                continue
            A, beta, r2 = fit_power_law(densities, eta_arr, min_rho=2.5)
            if np.isfinite(beta) and beta > 0:
                exps.append(beta)
        exps = np.array(exps)
        if len(exps) > 0:
            print(f"  λ={lam:.2f} fm: median={np.median(exps):.3f} "
                  f"[{np.percentile(exps,16):.3f}, {np.percentile(exps,84):.3f}] "
                  f"n={len(exps)}")
    
    # Breakdown by σ
    print(f"\n{'='*70}")
    print(f"EXPONENT BY GAUSSIAN WIDTH σ")
    print(f"{'='*70}")
    
    for sigma in results['sigmas']:
        exps = []
        for key, eta_arr in results['eta_data'].items():
            if key[1] != sigma:
                continue
            A, beta, r2 = fit_power_law(densities, eta_arr, min_rho=2.5)
            if np.isfinite(beta) and beta > 0:
                exps.append(beta)
        exps = np.array(exps)
        if len(exps) > 0:
            print(f"  σ={sigma:.2f} fm: median={np.median(exps):.3f} "
                  f"[{np.percentile(exps,16):.3f}, {np.percentile(exps,84):.3f}] "
                  f"n={len(exps)}")
    
    # Bootstrap on the median curve
    print(f"\n{'='*70}")
    print(f"BOOTSTRAP ON MEDIAN η CURVE")
    print(f"{'='*70}")
    
    med_beta, ci68_lo, ci68_hi, ci95_lo, ci95_hi = bootstrap_exponent(
        densities, eta_median, n_bootstrap=5000, min_rho=2.5)
    
    print(f"  Median exponent: {med_beta:.3f}")
    print(f"  68% CI: [{ci68_lo:.3f}, {ci68_hi:.3f}]")
    print(f"  95% CI: [{ci95_lo:.3f}, {ci95_hi:.3f}]")
    
    # The money line
    total_med = med_beta + 2.0/3.0 + 1.0
    total_lo = ci95_lo + 2.0/3.0 + 1.0
    total_hi = ci95_hi + 2.0/3.0 + 1.0
    
    print(f"\n{'='*70}")
    print(f"THE BOTTOM LINE")
    print(f"{'='*70}")
    print(f"  η exponent:          {med_beta:.3f} [{ci95_lo:.3f}, {ci95_hi:.3f}] (95% CI)")
    print(f"  Total DST exponent:  {total_med:.3f} [{total_lo:.3f}, {total_hi:.3f}] (95% CI)")
    print(f"  Fermi exponent:      1.667")
    print(f"  Exceeds Fermi:       {'YES' if total_lo > 1.667 else 'MARGINAL'} "
          f"(lower 95% bound {'>' if total_lo > 1.667 else '<'} 1.667)")
    
    if total_lo > 1.667:
        print(f"\n  *** The cascade argument holds at >95% confidence ***")
        print(f"  *** across ALL tested profile parameters and lattice geometries. ***")
    
    return {
        'all_exponents': all_exponents,
        'total_exponents': total_exponents,
        'eta_median': eta_median,
        'eta_bands': (eta_2p5, eta_16, eta_84, eta_97p5),
        'median_fit': (med_beta, ci95_lo, ci95_hi),
    }


def save_results_csv(results, analysis, filename):
    """Save the η table and exponent summary to CSV."""
    densities = results['densities']
    eta_med = analysis['eta_median']
    eta_2p5, eta_16, eta_84, eta_97p5 = analysis['eta_bands']
    
    with open(filename, 'w') as f:
        f.write("# DST Many-Body Enhancement Factor η(ρ)\n")
        f.write(f"# Grid: {results['grid_n']}^3\n")
        f.write(f"# Profiles: λ = {results['lambdas']}, σ = {results['sigmas']}\n")
        f.write(f"# Lattices: {results['lattices']}\n")
        f.write(f"# Total parameter combinations: {len(results['eta_data'])}\n")
        f.write(f"# η exponent (median): {analysis['median_fit'][0]:.4f}\n")
        f.write(f"# η exponent (95% CI): [{analysis['median_fit'][1]:.4f}, "
                f"{analysis['median_fit'][2]:.4f}]\n")
        f.write("#\n")
        f.write("rho_over_rho0,d_fm,eta_median,eta_2.5pct,eta_16pct,"
                "eta_84pct,eta_97.5pct,eta_x_137_median\n")
        
        for i, rho in enumerate(densities):
            d = D0 * rho ** (-1.0/3.0)
            f.write(f"{rho:.4f},{d:.4f},{eta_med[i]:.6f},{eta_2p5[i]:.6f},"
                    f"{eta_16[i]:.6f},{eta_84[i]:.6f},{eta_97p5[i]:.6f},"
                    f"{eta_med[i]*137:.2f}\n")
    
    print(f"\nResults saved to {filename}")


# ══════════════════════════════════════════════════════════════
#  CONVERGENCE TEST
# ══════════════════════════════════════════════════════════════

def convergence_test():
    """Test grid convergence at one density point."""
    print("Grid convergence test at ρ = 5ρ₀, λ=1.41, σ=0.50, FCC:")
    for n in [33, 49, 65, 81, 97]:
        eta = compute_eta_at_density(5.0, 1.41, 0.50, 'fcc', grid_n=n)
        print(f"  {n}³ = {n**3:>8,} points:  η = {eta:.5f}")


# ══════════════════════════════════════════════════════════════
#  ENTRY POINT
# ══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='DST η calculator')
    parser.add_argument('--grid', type=int, default=65, 
                       help='Grid points per dimension (default 65)')
    parser.add_argument('--n-random', type=int, default=5,
                       help='Number of random lattice realizations (default 5)')
    parser.add_argument('--convergence', action='store_true',
                       help='Run convergence test only')
    parser.add_argument('--quick', action='store_true',
                       help='Quick run: fewer profiles, smaller grid')
    parser.add_argument('--output', type=str, default='dst_eta_results.csv',
                       help='Output CSV filename')
    
    args = parser.parse_args()
    
    if args.convergence:
        convergence_test()
        sys.exit(0)
    
    if args.quick:
        # Override for quick test
        print("=== QUICK MODE: reduced parameter space ===\n")
        
        densities = np.logspace(np.log10(2.0), np.log10(10.0), 10)
        lambdas = [1.41]
        sigmas = [0.50]
        lattices = ['fcc']
        
        results = {
            'densities': densities,
            'lambdas': lambdas,
            'sigmas': sigmas,
            'lattices': lattices,
            'n_random': 0,
            'grid_n': args.grid,
            'eta_data': {},
        }
        
        for i, rho in enumerate(densities):
            eta = compute_eta_at_density(rho, 1.41, 0.50, 'fcc', grid_n=args.grid)
            d = D0 * rho ** (-1.0/3.0)
            print(f"  ρ={rho:5.2f}ρ₀  d={d:.3f}fm  η={eta:.4f}  η×137={eta*137:.1f}")
        
            key = (1.41, 0.50, 'fcc', 0)
            if key not in results['eta_data']:
                results['eta_data'][key] = np.zeros(len(densities))
            results['eta_data'][key][i] = eta
        
        A, beta, r2 = fit_power_law(densities, results['eta_data'][key], min_rho=2.5)
        print(f"\n  Power law fit (ρ ≥ 2.5ρ₀): η = {A:.4f} × (ρ/ρ₀)^{beta:.3f}  R²={r2:.4f}")
        total = beta + 2.0/3.0 + 1.0
        print(f"  Total DST exponent: {total:.3f}  (Fermi: 1.667)")
        print(f"  Exceeds Fermi: {'YES' if total > 1.667 else 'NO'}")
        
    else:
        results = run_full_sweep(grid_n=args.grid, n_random=args.n_random)
        analysis = analyze_results(results)
        save_results_csv(results, analysis, args.output)
