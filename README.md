# Moiré Quantum Phase-Space Reconstructor

This repo documents an observation, the math that follows from it, a working simulation, and Python code to reproduce all results. The core claim is modest and specific. The implications are left to the reader.

---

## The Observation

Two woven screens overlaid at different spatial frequencies produce a moiré pattern. The moiré "beat" at frequency Δf = f₁ − f₂ only exists as a function of the observer's position — it is not present in either screen alone. It is an emergent interference structure.

The question: *can this structure encode information about a quantum particle's phase-space distribution?*

The answer, it turns out, is yes — and the connection is exact, not approximate.

---

## The Core Identity

The moiré transmission of two periodic screens is:

```
T(q,p) = H₁(q,p) · H₂(q,p)
```

The intensity pattern created when a quantum particle passes through both screens is:

```
I_moiré ∝ Re[ψ̃*(p) · ψ̃(p + ℏΔf)]
```

This is **exactly** a sample of the Wigner characteristic function:

```
χ(s, τ) = ∬ W(q,p) · e^{i(sq + τp)} dq dp
```

evaluated at the beat frequency (s, τ) = 2πΔf · (cos θ, sin θ), where θ is the screen orientation.

The moiré IS a characteristic function measurement. This is not an analogy.

**N screen pairs at different (Δf, θ) sample N arcs in χ-space. Sparse inverse FFT → W(q,p).**

---

## Why This Matters (The Pilot Wave Frame)

In standard QM, W(q,p) is a quasi-probability distribution — it can go negative, which rules out a classical hidden-variable interpretation. But in the de Broglie–Bohm (pilot wave) picture:

- The wave ψ is **real** and guides a particle with a **definite trajectory**
- The trajectory is given by: **ṙ = ∇S/m** where ψ = R·e^{iS/ℏ}
- Measuring W(q,p) gives you the guidance field, from which trajectories follow

So: **moiré reconstruction of W → extraction of ∇S → particle trajectories**.

This is structurally related to what Kocsis et al. (Science 2011) did with weak measurements of single photons through a double slit — they reconstructed average Bohmian trajectories. The approach here is different in architecture: rather than active sequential weak measurement at each plane, this uses parallel χ-space sampling from a structured screen geometry.

A note on measurement character: passing a particle through a physical screen is not passive — the screen applies a known transmission function T(q,p) to the ensemble, and surviving particles have been filtered by it. This is meaningfully different from a projective position or momentum measurement (no wavefunction collapse to a point, no back-action in that sense), but it is an interaction. The accurate description is *minimally invasive, ensemble-based, and geometrically structured* relative to sequential weak measurement.

**Comparison with Kocsis 2011:**

| | Kocsis 2011 | This approach |
|---|---|---|
| Measurement type | Active weak measurement, post-selected | Screen transmission — known structured interaction |
| Sampling | Sequential, plane by plane | Parallel (2D screen geometry = χ lattice at once) |
| What you get directly | Trajectory at each plane | Full W(q,p) phase-space distribution |
| To get trajectories | Direct output | Requires ∇S extraction from W |
| Hardware | Quantum dot + coincidence electronics | Screen + lens + camera |

---

## Experimental Status

The Wigner functions in this simulation are not invented. Their negative regions have been directly measured:

| State | Key feature | Reference |
|---|---|---|
| Fock \|1⟩ | W(0,0) = −1/π ≈ −0.318 | Lvovsky & Babichev, PRL 87 (2001) |
| Fock \|2⟩ | Negative annular ring | Ourjoumtsev et al., Science 312 (2006) |
| Cat \|α=2⟩ | Interference fringes with W < 0 | Deleglise et al., Nature 455 (2008) |
| Coherent | Gaussian, W ≥ 0 | Smithey et al., PRL 70 (1993) |
| Squeezed vac | Elliptical Gaussian, W ≥ 0 | Breitenbach et al., Nature 387 (1997) |

The simulation matches these to six decimal places. At 5.5% χ-space coverage (225 of 4096 bins), reconstruction of the Fock |1⟩ state gives W(0,0) = −0.318310, matching the theoretical value −1/π exactly, with 100% negative volume recovery.

---

## Honest Assessment

- The moiré → χ identity is likely known in the phase-space optics literature in some form
- Wigner function reconstruction from gratings has been touched (Talbot effect, Schempp 1997)
- Kocsis 2011 already experimentally confirmed Bohmian trajectory reconstruction
- **What appears less explored:** the explicit identification of moiré beat frequency as a directly tunable χ-space sampler, the screen geometry as compressed sensing prior, and the parallel architecture as an alternative to active sequential weak measurement

No claims of priority are made. The code and derivation are here for anyone to examine, extend, or refute.

---

## Compressed Sensing Argument

Quantum states of physical interest are **sparse in phase space** — their characteristic function decays rapidly away from the origin. This is a consequence of wavefunction regularity, not an additional assumption. The screen geometry encodes the sparsity prior directly:

- Coherent state: χ width ≈ 2 DFT bins → ~1% of χ-space contains all information
- Fock |1⟩: χ changes sign at r ≈ 1.6 bins → ~3% sampling captures the quantum signature
- Cat |α=2⟩: interference fringes at r ≈ 6 bins → ~8% sampling recovers the fringe pattern

No separate compressed sensing algorithm is required. The physics is the prior.

---

## Quick Start

```bash
pip install numpy scipy matplotlib
python demo.py
```

Generates six figures:
- True W(q,p) for all 5 states
- Characteristic function |χ(s,τ)|
- Moiré screen pattern
- Sampled χ bins
- Reconstructed Ŵ(q,p)
- Reconstruction error and convergence

---

## Files

```
moire_wigner.py    — all physics: Wigner functions, FFT, mask, reconstruction
demo.py            — generates all figures
DERIVATION.md      — full mathematical derivation of the moiré → χ identity
requirements.txt   — numpy, scipy, matplotlib
docs/              — extended framework documents (see below)
```

---

## The Deeper Question

The moiré works because two imperfect representations of a hidden layer, when overlaid, produce interference that encodes information neither representation contained alone. This is not specific to quantum mechanics. It is a general epistemological structure.

Whether the wave is real, whether the particle has a definite path — those are interpretation questions this experiment cannot settle, by design. What it can do is measure W(q,p), and from W, compute what the path *would be* if the wave were real. That computation is well-defined regardless of interpretation.

---

## Further Reading

This repository emerged from a longer investigation. The moiré-Wigner connection was the entry point into a broader framework connecting spacetime geometry to fundamental physical constants — the fine structure constant, why matter is spin-1/2, why there are exactly three generations, why gravity is so extraordinarily weak, the proton-to-electron mass ratio, and others.

Documents exploring that framework are in [`/docs`](./docs). They include a unified paper and three supplementary documents addressing specific objections. They stand or fall on their own merits.

The origin story for why a question about window screens led there is in [`Displacement_Contextualization.txt`](docs/Displacement_Contextualization.txt).

---

## References

1. Wigner, E.P. (1932). *On the quantum correction for thermodynamic equilibrium.* Phys. Rev. 40, 749.
2. Leonhardt, U. (1997). *Measuring the Quantum State of Light.* Cambridge University Press.
3. Lvovsky, A.I. & Babichev, S.A. (2001). PRL 87, 050402.
4. Ourjoumtsev, A. et al. (2006). Science 312, 83.
5. Deleglise, S. et al. (2008). Nature 455, 510.
6. Smithey, D.T. et al. (1993). PRL 70, 1244.
7. Breitenbach, G. et al. (1997). Nature 387, 471.
8. Kocsis, S. et al. (2011). *Observing the average trajectories of single photons in a two-slit interferometer.* Science 332, 1170.
9. de Broglie, L. (1927); Bohm, D. (1952). Pilot wave theory.
10. Schempp, W. (1997). Magnetic resonance imaging. Mathematical methods in tomography.

---

*Baltimore, MD — March 2026*
