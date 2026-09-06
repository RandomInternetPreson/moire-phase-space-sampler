# Prime-Free Prime Counting, and a Negative Result on Geometric Provenance
https://doi.org/10.5281/zenodo.22073845

**A prime-free pipeline: Riemann–Siegel zeros → Riemann's J(x) → π(x) — and controlled tests showing that a truncated modular surface does not supply the zeros**

Aaron Alai · Baltimore, MD · ORCID [0000-0003-4704-7271](https://orcid.org/0000-0003-4704-7271)

Version 4 — 6 September 2026. v1 April 2026; v3 August 2026. This version withdraws the geometric-provenance claim of v1–v3.

Code: [github.com/RandomInternetPreson/moire-phase-space-sampler/tree/main/Geometric_Prime_Counter](https://github.com/RandomInternetPreson/moire-phase-space-sampler/tree/main/Geometric_Prime_Counter)

---

## Abstract

I present a self-contained implementation of prime counting that uses no primes, no sieve, and no number-theoretic lookup tables as input. Nontrivial zeros of ζ(*s*) are located by Newton's method on the Riemann–Siegel *Z* function, which sums over integers and never touches a prime; Riemann's explicit formula for *J*(*x*) then yields π(*x*) by a prime-free recursion. With 17 zeros, π(*x*) is recovered to within 0.25% for 100 ≤ *x* ≤ 10⁶ and within 0.1% for *x* ≥ 10⁴; the shortcut π(*x*) ≈ ψ(*x*)/log *x* is shown to discard the information the zeros carry.

In versions 1–3 of this note I claimed a further step: that eigenvalues of the hyperbolic Laplacian on the modular surface SL(2,ℤ)\ℍ², truncated by a Dirichlet wall in the cusp, approximate the zeta zeros and therefore give the pipeline a geometric origin. In this version I withdraw that claim. Control experiments (§4.3) show that the surface eigenvalues locate zeros no better than evenly spaced or random seeds laid across the same interval, at every Newton search radius tested; the raw eigenvalues sit only ~20% closer to zeros than random points (1.4–1.9σ one-sided across three wall heights), a real but unusable signal. I also correct the theoretical premise (§2.2): the wall's continuum eigenvalues obey *Y*^(2*it*) = −φ(½+*it*) and form a near-uniform comb whose roots are no closer to zeros than random; the zeros of ζ are the *poles* of the scattering matrix φ — resonances of the surface — not eigenvalues of any self-adjoint truncation. The pipeline remains prime-free. It is not geometric.

---

## 1. Introduction

Prime counting is one of the oldest problems in mathematics. Given *x*, how many primes are ≤ *x*? The Sieve of Eratosthenes answers exactly in O(*n* log log *n*); Meissel–Lehmer and Lagarias–Miller–Odlyzko extend exact counting to enormous *x* in sub-linear time. All of these are arithmetic: divisibility is the primary input. The analytic approximations — li(*x*), Riemann's R(*x*) — are typically initialised from tabulated zeta zeros, themselves computed by arithmetic algorithms.

I ask two questions here, and they receive opposite answers.

**Can π(*x*) be computed with no prime as input?** Yes. Riemann–Siegel locates the zeros from integers alone; Riemann's 1859 formula converts zeros to *J*(*x*) and thence to π(*x*). §4.1 demonstrates this with runnable code, to 0.02% at *x* = 10⁶.

**Can those zeros be obtained from geometry — from the shape of a surface?** In earlier versions I said yes, via the Dirichlet-truncated modular Laplacian. Controlled tests (§4.3) say no: the surface's contribution is an interval to search, plus a weak proximity signal that does not survive any reasonable test of usefulness. I report the negative result in full, with the diagnosis, because a clean failure is more useful than an unexamined success.

**A note on "prime-free" vs "arithmetic-free."** The Riemann–Siegel sum runs over integers *n* = 1, …, *N*, so integers and their logarithms appear. What does not appear is any use of primes, factorisation, or divisibility — the structures that make prime counting hard. The pipeline is therefore *prime-free*, not arithmetic-free.

---

## 2. Background

### 2.1 The Modular Surface

The modular surface is SL(2,ℤ)\ℍ², the quotient of the hyperbolic upper half-plane by the modular group. Its fundamental domain {|*z*| ≥ 1, |Re *z*| ≤ ½} is a horn: compact at the arc, opening into a cusp as *y* → ∞. The hyperbolic Laplacian is −Δ = −*y*²(∂²ₓ + ∂²ᵧ). Its L² spectrum on the full surface consists of the Maaß cusp forms (discrete, eigenvalues ¼ + *t*², obeying Weyl's law *N*(*T*) ≈ *T*²/12) and a continuous spectrum on [¼, ∞) carried by the Eisenstein series *E*(*z*, *s*), whose constant Fourier term is *y*^*s* + φ(*s*)*y*^(1−*s*) with scattering matrix

> φ(*s*) = Λ(2*s* − 1) / Λ(2*s*),  Λ(*s*) = π^(−*s*/2) Γ(*s*/2) ζ(*s*).

### 2.2 Where the zeta zeros actually live on the surface (corrected)

The zeros of ζ enter through φ. Writing ρ = ½ + *i*γ for a nontrivial zero, φ has a pole where Λ(2*s*) = 0, i.e. at *s* = ρ/2 = ¼ + *i*γ/2, and a zero at *s* = (1+ρ)/2 = ¾ + *i*γ/2. On the critical line Re *s* = ½, |φ| = 1 and φ has neither zeros nor poles. The zeta zeros are therefore *resonances* of the modular surface — poles of the scattering matrix off the spectral line — and equivalently eigenvalues of the non-self-adjoint Lax–Phillips semigroup generator (Faddeev–Pavlov 1972; Lax–Phillips 1976). They are not L² eigenvalues of the Laplacian, and no self-adjoint truncation makes them so.

**What a Dirichlet wall does.** Placing a wall at *y* = *Y* discretises the continuous spectrum. The constant term must vanish there: *Y*^*s* + φ(*s*)*Y*^(1−*s*) = 0, which for *s* = ½ + *it* reads

> *Y*^(2*it*) = −φ(½ + *it*),  i.e.  2*t* log *Y* = arg(−φ(½+*it*)) + 2π*n*.

The term 2*t* log *Y* dominates. The resulting eigenvalues form a comb in *t* with density (log *Y* + log(*t*/π))/π: at *Y* = 50 the predicted spacing is 0.545 and the measured spacing of the roots is 0.551 (§4.3, Control 3). The scattering phase contributes only a wobble on this comb. The roots of the condition sit a median 0.410 from the nearest γ/2, where a point placed at random would sit gap/4 = 0.433: *the wall does not target zeros*. Increasing *Y* makes the comb denser, so some root is always near any zero, but no root is attracted to one. In versions 1–3 I stated that these eigenvalues "cluster near values related to zeta zeros" and that "the approximation improves as *Y* → ∞"; I withdraw both statements.

**Colin de Verdière's pseudo-Laplacian.** The pseudo-Laplacian Δ_*a* (Colin de Verdière 1983) forces the constant Fourier term to vanish for *y* > *a* while leaving the other Fourier modes free. Its continuum-derived eigenvalues satisfy the same condition *a*^(2*s*−1) = −φ(*s*) and form the same comb; it differs from the Dirichlet wall only in how it treats the non-constant modes (which decay in the cusp regardless). In versions 1–3 I described the pseudo-Laplacian as an operator whose eigenvalues "coincide exactly with the nontrivial zeros of ζ." That description is not correct and I withdraw it; the reader should consult the primary source. The rigorous connection between the modular surface and the zeros runs through the scattering matrix's poles, not through the L² spectrum of any truncation.

**Consequence for the pipeline.** The surface eigenvalues, converted by γ = 2*t*, are seeds spanning an interval on the critical line. Whether they carry any information about zero *positions* beyond that interval is an empirical question, answered in §4.3.

### 2.3 The Riemann–Siegel Z Function

*Z*(*t*) = *e*^(*i*θ(*t*)) ζ(½ + *it*) with θ(*t*) = Im log Γ(¼ + *it*/2) − (*t*/2) log π is real for real *t*, and its zeros are the γ. The Riemann–Siegel formula gives

> *Z*(*t*) = 2 Σ\_{*n*=1}^{*N*} cos(θ(*t*) − *t* log *n*) / √*n* + remainder,  *N* = ⌊√(*t*/2π)⌋.

Elementary functions at integer arguments; no primes anywhere. The code keeps only the first remainder term *C*₀, so the *Z* it evaluates is an approximation whose zeros differ from the true γ by at most 7.5 × 10⁻³ for γ < 70 (worst case near γ = 25.01) and by a median 2.6 × 10⁻⁴ over the first 108 zeros, measured against `mpmath.zetazero`. Newton's method converges to 10⁻¹² on that approximation from any starting point within reach of a sign change; earlier versions of this note said "to machine precision," which was true of the iteration but not of γ. The effect on π(10⁶) with 17 zeros is about one count (78 514.6 against 78 515.9 with exact zeros).

### 2.4 Riemann's Explicit Formula

Riemann counts prime powers with weight 1/*k*:

> *J*(*x*) = Σ\_{*p*^*k* ≤ *x*} 1/*k* = π(*x*) + ½ π(*x*^(1/2)) + ⅓ π(*x*^(1/3)) + …

and expresses *J* through the zeros:

> *J*(*x*) = li(*x*) − Σ\_ρ li(*x*^ρ) − log 2 + ∫\_*x*^∞ d*t* / (*t*(*t*² − 1) log *t*).

Inverting, π(*x*) = *J*(*x*) − ½ π(√*x*) − ⅓ π(*x*^(1/3)) − …, which terminates once *x*^(1/*k*) < 2. Every ingredient is analysis. The truncation error from using zeros with γ ≤ *T* is O(*x* log² *x* / *T*) in *J*.

**Why not ψ(*x*)/log *x*.** The shortcut through Chebyshev's ψ drops the *x*/log² *x* term of π(*x*) and reduces the pipeline to the prime number theorem's leading order, with error ≈ 1/log *x* — 7.8% at 10⁶ — *regardless of which zeros are supplied*. Any pipeline that ends with this step is measuring the conversion, not the zeros. §4.1 reports both conversions side by side.

---

## 3. The Pipeline

```
  1.  (optional) Discretise SL(2,Z)\H² with a Dirichlet wall at height Y;
      take eigenvalue seeds γ = 2t                              — PART 1
  2.  Newton's method on Riemann–Siegel Z(t) from each seed      — PART 2
  2b. Or simply scan Z(t) for sign changes: no seeds needed      — PART 2b
  3.  Riemann's explicit formula: {γ} → J(x) → π(x)              — PART 3
  4.  Controls: do the step-1 seeds beat a ruler?                — PART 4

  Prime-free at every step.  Step 1 is retained for the controls.
```

### 3.1 Discretisation

Coordinates (*x*, *u*) with *y* = *e*^*u*; domain *x* ∈ (0, ½), *u* ∈ (*u*\_min, log *Y*), a non-uniform grid clustered at both boundaries. −Δ = −(∂²ᵤ − ∂ᵤ) − *e*^(2*u*) ∂²ₓ with Dirichlet on the arc and the wall, Neumann at *x* = 0 and ½. The sparse matrix is symmetrised and its 80 lowest eigenvalues found by ARPACK. At *Y* = 50, Nx = 50, Nu = 300 this yields 48 seeds spanning γ = 12.65 … 69.75. That count is consistent with Weyl's area term for the half domain (*T*²/24 ≈ 51 at *T* = 35): the seeds are predominantly interior modes and comb roots, as §2.2 predicts.

**A remark on the boundary conditions (new in v4).** Neumann at *x* = 0 selects eigenfunctions even under *x* ↦ −*x*. The inversion *S*: *z* ↦ −1/*z* maps the arc |*z*| = 1 to itself by (*x*, *y*) ↦ (−*x*, *y*), so on the modular surface an even eigenfunction satisfies a *Neumann* condition on the arc; only odd eigenfunctions satisfy Dirichlet there. The Dirichlet-arc/Neumann-*x* combination in the code — inherited unchanged from v1 — is therefore not the SL(2,ℤ) gluing, and its interior modes are not Maaß cusp forms of the modular group. I record this so that the operator is described correctly, not because it bears on the result: the comb analysis of §2.2 depends only on the cusp end and is unchanged in form (the phase wobble on the comb is the scattering phase of this mixed problem rather than exactly arg φ); Controls 1, 2 and 2b test the seeds as numbers, whatever operator produced them; and Control 3 is computed from the exact φ of SL(2,ℤ)\ℍ². If anything the remark strengthens the negative result: the operator that was supposed to know about ζ was not even the right operator, and the controls would have caught the claim either way.

### 3.2 Newton Refinement (corrected)

For each seed: bracket a sign change of *Z* within radius 2.0; bisect; Newton-polish to 10⁻¹² on the one-term Riemann–Siegel *Z* (§2.3). **47 of the 48 seeds converge.** They funnel onto **17 distinct zeros** — between one and five seeds per zero — which are the 17 nontrivial zeros below γ = 70 (γ₁₇ = 69.55; γ₁₈ = 72.07 lies outside the seed span), each located to within 7.5 × 10⁻³. In versions 1–3 I reported that "60–70% of candidates fail to bracket"; that was a misreading of the deduplicated count. Nothing fails; the seeds are redundant. Refinement takes ~30 ms.

### 3.3 Prime Counting

*J*(*x*) is evaluated with `scipy.special.expi` for li(*x*) and `exp1` at complex argument for the zero sum; the tail integral is Σ\_{*k*≥1} E₁(2*k* log *x*). The recursion for π(*x*) descends through *x*^(1/*k*) until the argument drops below 2. A query costs ~2 ms.

### 3.4 Global Scan

`scan_zeros` walks *Z*(*t*) from *t* = 10 to 250 in steps of 0.05, brackets every sign change, and refines it: 108 zeros in ~40 ms, no seeds required, equally prime-free. **This is the step that actually supplies the zeros.** It also provides the reference set for the controls in §4.3, so no external zero table is used anywhere.

---

## 4. Results

### 4.1 Accuracy (unchanged from v3)

| *x* | π(*x*) exact | ψ/log *x*, 17 zeros | Riemann, 17 zeros | Riemann, 108 zeros | li(*x*) alone |
|---:|---:|---:|---:|---:|---:|
| 100 | 25 | 21 (17.1%) | 25.3 (1.04%) | 25.0 (0.19%) | 30.1 (20.5%) |
| 1,000 | 168 | 144 (14.1%) | 168.0 (0.01%) | 167.6 (0.24%) | 177.6 (5.7%) |
| 10,000 | 1,229 | 1,088 (11.5%) | 1,229.9 (0.07%) | 1,229.7 (0.06%) | 1,246.1 (1.4%) |
| 100,000 | 9,592 | 8,689 (9.4%) | 9,591.0 (0.01%) | 9,589.7 (0.02%) | 9,629.8 (0.39%) |
| 1,000,000 | 78,498 | 72,369 (7.8%) | 78,514.6 (0.02%) | 78,502.8 (0.01%) | 78,627.5 (0.17%) |

*Table 1. The ψ/log x column is the prime number theorem to within a few counts, independent of the zeros. Riemann's formula with the same 17 zeros beats li(x) at every x. The 0.01% entries are partly fortunate; the fair summary is the bound.*

### 4.2 Speed

| Method | *x* = 1M, per call | Notes |
|--------|---:|------|
| Sieve (exact) | ~6.5 ms | Rebuilds each call |
| li(*x*) | ~0.01 ms | No zeros; 0.17% error |
| Riemann query, 17 zeros | ~2 ms | After ~1 s setup; 0.02% error |

*Table 2. Break-even against the sieve at x = 1M occurs at roughly 150–300 queries. This is not a production algorithm (§5).*

### 4.3 Control experiments (new in v4)

The question: do the 48 surface eigenvalues locate zeros better than 48 numbers that know nothing about the surface? Reference zeros are those found by the prime-free scan of §3.4. All controls are in `PART 4` of the code; Control 1 runs by default and the rest with `python geometric_prime_counter_v4.py --controls`.

**Control 1 — distinct zeros recovered, and raw distance to the nearest zero.** Same count, same span (γ = 12.65 … 69.75, containing 17 zeros, so gap/4 = 0.840).

| seeds | distinct zeros (radius 2.0) | median distance to nearest zero | mean |
|---|---:|---:|---:|
| surface eigenvalues | 17 | 0.692 | 0.819 |
| evenly spaced, same span | 17 | 0.924 | 1.006 |
| uniform random, same span (20 draws) | 15.7 | 0.895 | 1.018 |

A ruler recovers every zero the surface recovers. The raw eigenvalues do sit closer to zeros than random points — 0.69 against 0.84–0.90 — by about 20%. Control 2b below asks how unusual that is.

**Control 2 — Newton search radius sweep.** With mean zero spacing 3.36, a radius of 2.0 sweeps 4.0 units and finds a zero from almost any seed. Shrinking the radius asks whether the surface knows where the zeros are.

| radius | surface | evenly spaced | random (6 draws) |
|---:|---:|---:|---:|
| 0.10 | 3 | 1 | 3.5 |
| 0.20 | 8 | 5 | 4.5 |
| 0.35 | 8 | 12 | 8.2 |
| 0.50 | 13 | 16 | 10.8 |
| 0.75 | 15 | 17 | 12.8 |
| 1.00 | 15 | 17 | 14.7 |
| 1.50 | 17 | 17 | 15.5 |
| 2.00 | 17 | 17 | 15.5 |

The surface edges ahead of both controls only at radius ≤ 0.2, where it recovers 8 of 17 — the weak signal of Control 1 made visible. At every radius from 0.35 up it ties or loses to a ruler. There is no regime in which the surface is the reason the zeros are found.

**Control 2b — wall-height sweep and significance of the proximity signal.** For each of *Y* = 20, 50, 100 I rebuild the surface, take the median nearest-zero distance of its seeds, and compare it with the distribution of that same statistic over 3000 sets of uniformly random seeds on the same span (`proximity_significance` in the code).

| *Y* | seeds | seed span (γ) | distinct zeros, surface / even | surface median | random-set median | *z* | one-sided *p* |
|---:|---:|------------|-------:|------:|------:|---:|-----:|
| 20 | 50 | 5.34–72.08 | 17 / 18 | 0.686 | 0.975 | 1.88 | 0.024 |
| 50 | 48 | 12.65–69.75 | 17 / 17 | 0.692 | 0.874 | 1.37 | 0.082 |
| 100 | 50 | 9.15–72.16 | 18 / 18 | 0.686 | 0.913 | 1.61 | 0.045 |

The seeds sit closer to zeros than random points at every wall height, by 1.4–1.9σ one-sided, and the three heights are not independent samples (the same low-lying zeros dominate each span). A weak signal is present. It is not enough to bracket a zero without a search radius comparable to the zero spacing, and at *Y* = 20 the surface recovers one zero fewer than the ruler.

**Control 3 — what the wall quantises.** The roots of *Y*^(2*it*) = −φ(½+*it*) on *t* ∈ [5, 36], computed with mpmath (ζ used here for diagnosis of the operator only): 58 roots, mean spacing 0.551 against the density prediction 0.545. Median distance from a root to the nearest γ/2: 0.410, against gap/4 = 0.433 for random points. The comb does not target zeros. This is the operator-level version of Controls 1, 2 and 2b — it uses the exact scattering matrix of SL(2,ℤ)\ℍ², so it is independent of the boundary-condition remark in §3.1 — and it confirms that the negative result is structural, not a grid artefact.

---

## 5. Discussion

### 5.1 What survives

The prime-free pipeline stands: Riemann–Siegel supplies zeros from integers, Riemann's formula converts them to π(*x*), and the demonstration that this conversion beats ψ(*x*)/log *x* by two orders of magnitude is unchanged. The code is self-contained and reproduces Table 1 in about one second.

### 5.2 What does not survive

The claim of geometric provenance. Measured four ways, the truncated Laplacian's eigenvalues locate zeros no better than a ruler. The Newton step, which I described in earlier versions as "snapping" geometric candidates to true zeros, is in fact doing the finding; the search radius of 2.0 against a zero spacing of 3.36 guarantees success from essentially any seed. The zeros in Table 1 came from Riemann–Siegel. The geometry contributed the interval [12.65, 69.75].

The theoretical premise does not survive either (§2.2). The wall discretises the continuum into a comb governed by 2*t* log *Y*; the zeros are the poles of φ and are not eigenvalues of the truncated operator, nor of the pseudo-Laplacian. The 20% proximity signal is the scattering-phase wobble riding on the comb, and it is the right size for that.

### 5.3 What the experience is worth

Two things. First, the structural fact: the fraction of the surface's L² spectrum that is zeta-related decays like log *T* / *T*, because cusp forms obey Weyl's *T*² law while the continuum contributes *T* log *T*. Scaling this construction up cannot help. Second, the method: the controls of §4.3 are the tests any claim of geometric provenance must pass, and they cost a few seconds to run. They should have been in v1.

### 5.4 Not a production algorithm

Lagarias–Miller–Odlyzko and Deléglise–Rivat count exactly in sub-linear time; the analytic method of Lagarias–Odlyzko achieves *x*^(1/2+ε). Computing zeros to height *T* costs *T*^(1+ε) by Odlyzko–Schönhage, and a sparse eigensolve on a grid resolving oscillations up to *t* is not competitive with any of these. Suggestions that a "phase-locked" geometric solver could beat combinatorial methods founder on the same point as §2.2: the phase to be locked to is arg φ(½+*it*), and evaluating φ is evaluating ζ.

### 5.5 Open directions

If a genuinely geometric route to the zeros is wanted, it runs through the resonances — the poles of φ — rather than through any L² spectrum. That is a scattering problem, non-self-adjoint, and its numerics are a research problem in their own right; the Lax–Phillips generator is the object. Whether it can be discretised without computing ζ is exactly the question, and I do not answer it here.

---

## 6. Code

`geometric_prime_counter_v4.py`, single file, numpy + scipy (mpmath optional for Control 3).

```bash
pip install numpy scipy
python geometric_prime_counter_v4.py             # pipeline + Control 1
python geometric_prime_counter_v4.py --controls  # + radius sweep + quantisation roots
```

Functions: `build_surface`, `refine_zero`, `refine_all`, `scan_zeros`, `riemann_J`, `count_primes`, `count_primes_psi`, `sieve_exact`; new in v4: `control_table`, `radius_sweep`, `proximity_significance`, `wall_height_sweep`, `quantization_roots`, `run_controls`. The default run takes about 2 s; `--controls` takes about 20 s (three further surface builds and the mpmath evaluation of φ).

---

## 7. Changes from v3 (v4, 6 September 2026)

| | v3 said | v4 says | basis |
|---|---|---|---|
| Surface eigenvalues | approximate ζ zeros via the scattering phase | form a comb of density (log *Y* + log(*t*/π))/π; no closer to zeros than random | §2.2, Control 3 |
| Pseudo-Laplacian | eigenvalues coincide exactly with ζ zeros | same constant-term condition, same comb; zeros are poles of φ, not eigenvalues | §2.2 |
| Newton step | 60–70% of candidates fail to bracket | 47/48 converge, funnelling onto 17 distinct zeros | §3.2 |
| Provenance | zeros came from geometry | zeros came from Riemann–Siegel; geometry supplied an interval and a ~20% proximity signal (1.4–1.9σ) | §4.3 |
| Zero precision | "machine precision," 10⁻¹² | 10⁻¹² on the one-term Riemann–Siegel *Z*; ≤ 7.5 × 10⁻³ in γ below 70; ~1 count at *x* = 10⁶ | §2.3, §3.2 |
| Boundary conditions | Neumann at *x* = 0, ½ "the mod-group identifications" | Dirichlet arc + Neumann *x* is not the *S*-gluing; interior modes are not Maaß forms; controls unaffected | §3.1 |
| Scaling | more zeros need larger *Y* and finer grid | zeta fraction of the spectrum decays like log *T*/*T*; scaling cannot help | §5.3 |
| Pipeline accuracy, prime-free status, ψ vs *J* result | — | unchanged | §4.1 |

The v3 code and README are retained in the repository as `geometric_prime_counter_v3.py` and `README_v3.md`. The error was caught by a control experiment that should have been run before v1: seed the Newton step with a ruler and see whether anything is lost. Nothing was.

---

## References

1. B. Riemann, *Über die Anzahl der Primzahlen unter einer gegebenen Grösse*, Monatsberichte der Berliner Akademie (1859).
2. C. L. Siegel, *Über Riemanns Nachlass zur analytischen Zahlentheorie*, Quellen und Studien zur Geschichte der Mathematik **2** (1932) 45–80.
3. A. Selberg, *Harmonic analysis and discontinuous groups in weakly symmetric Riemannian spaces with applications to Dirichlet series*, J. Indian Math. Soc. **20** (1956) 47–87.
4. L. D. Faddeev and B. S. Pavlov, *Scattering theory and automorphic functions*, Proc. Steklov Inst. Math. **27** (1972) 161–193.
5. P. D. Lax and R. S. Phillips, *Scattering Theory for Automorphic Functions*, Annals of Mathematics Studies **87**, Princeton (1976).
6. Y. Colin de Verdière, *Pseudo-laplaciens II*, Annales de l'Institut Fourier **33** (1983) 87–113.
7. J. C. Lagarias and A. M. Odlyzko, *Computing π(x): an analytic method*, J. Algorithms **8** (1987) 173–191.
8. A. M. Odlyzko and A. Schönhage, *Fast algorithms for multiple evaluations of the Riemann zeta function*, Trans. AMS **309** (1988) 797–809.
9. H. M. Edwards, *Riemann's Zeta Function*, Academic Press (1974); Dover (2001).
10. M. V. Berry and J. P. Keating, *The Riemann zeros and eigenvalue asymptotics*, SIAM Review **41** (1999) 236–266.

---

## AI-assisted preparation

Large-language-model assistance (Claude, Anthropic) was used during this research for checking derivations, drafting text, and writing verification code, under the sole direction of the author, who takes full responsibility for the entire content. The control experiments of §4.3 were proposed and first run by the assistant in September 2026 while auditing v3; a second, independent assistant instance then audited the v4 draft and code, reproduced every table, and identified the precision overstatement (§2.3) and the boundary-condition issue (§3.1) recorded in this version. I reran everything on my own machine before release. All numerical claims are reproducible from the released script. The reader is asked to verify the statement of Colin de Verdière's theorem in §2.2 against the primary source; the correction was reconstructed from the operator's definition rather than from the paper.
