"""
The vise toy (item #18 baseline): naive realized-context Gibbs family on the
n = 2 reduced classes, CHSH at the Tsirelson settings.

Family: rho_s(class) ~ exp(beta * tau_s * sign_s(class)) over the 8 classes
(P, g1, g2) in GF(2)^3, sign_s = (-1)^(P + s.g), tau = (+,+,+,-) the sign
pattern of the CHSH targets (magnitude 1/sqrt2). Marginals vanish by the
fiber-averaged class construction.

Closed forms to verify (the baseline any derived ring functional must beat):
  E(s; beta) = tau_s * tanh(beta)
  S(beta)    = 4 * tanh(beta)      -- local (0) -> Tsirelson -> PR-box (4)
  F(beta)    = (1/2) * tanh(beta)  -- Hall fraction paid

Landmarks:
  beta* = arctanh(1/sqrt2) ~ 0.881374 : S(beta*) = 2*sqrt2 exactly
  F(beta*) = 1/(2*sqrt2) ~ 0.353553   : 2.561x Hall's floor (sqrt2-1)/3
  S > 2*sqrt2 for beta > beta*        : excluded at ~30 sigma
                                        (Poh et al., PRL 115, 180408 (2015):
                                         S = 2.82759 +/- 0.00051)

The mechanism owes: cos-saturation of the ring functional, OR beta-pinning
<= beta*. Generic free beta passes THROUGH the quantum bound. Pre-registered
standard: a generic result gets filed at face value.
"""
import numpy as np, itertools

effs = np.array(list(itertools.product([0, 1], repeat=3)))
P, G = effs[:, 0], effs[:, 1:]
settings = np.array([(0, 0), (0, 1), (1, 0), (1, 1)])
tau = np.array([1, 1, 1, -1], dtype=float)
sign = np.where((P[None, :] + settings @ G.T) % 2 == 1, -1.0, 1.0)  # [s, class]

def family(beta):
    rhos = []
    for si in range(4):
        w = np.exp(beta * tau[si] * sign[si])
        rhos.append(w / w.sum())
    return np.array(rhos)

def observables(beta):
    r = family(beta)
    E = np.array([(r[si] * sign[si]).sum() for si in range(4)])
    S = E[0] + E[1] + E[2] - E[3]
    F = max(0.5 * np.abs(r[i] - r[j]).sum()
            for i in range(4) for j in range(i + 1, 4))
    return E, S, F

print(f"{'beta':>10} {'E/tau':>10} {'S':>10} {'F':>10}")
for beta in [0.0, 0.25, 0.5, 0.881374, 1.5, 3.0, 8.0]:
    E, S, F = observables(beta)
    print(f"{beta:>10.6f} {E[0]:>10.6f} {S:>10.6f} {F:>10.6f}")
    # closed forms, exact
    assert np.allclose(E, tau * np.tanh(beta), atol=1e-12)
    assert abs(S - 4 * np.tanh(beta)) < 1e-12
    assert abs(F - 0.5 * np.tanh(beta)) < 1e-12

bstar = np.arctanh(1 / np.sqrt(2))
E, S, F = observables(bstar)
hall = (np.sqrt(2) - 1) / 3
print(f"\nbeta* = arctanh(1/sqrt2) = {bstar:.6f}")
print(f"S(beta*) = {S:.12f}  vs 2*sqrt2 = {2*np.sqrt(2):.12f}  "
      f"(diff {abs(S-2*np.sqrt(2)):.1e})")
print(f"F(beta*) = {F:.6f}  vs Hall floor {hall:.6f}  ratio {F/hall:.3f}")
assert abs(S - 2 * np.sqrt(2)) < 1e-12
assert abs(F - 1 / (2 * np.sqrt(2))) < 1e-12

# overshoot check: the family DOES pass through Tsirelson
E, S, F = observables(bstar * 1.5)
print(f"\nS(1.5*beta*) = {S:.6f} > 2*sqrt2: the vise is real "
      f"(nothing in the naive family pins beta)")
assert S > 2 * np.sqrt(2)

print("\nALL ASSERTIONS PASSED")
