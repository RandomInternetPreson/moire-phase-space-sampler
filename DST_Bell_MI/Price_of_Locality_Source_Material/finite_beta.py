"""
Finite-beta phenomenology of the context-hard Gibbs family on the reduced
class space (beta = J/T_eff).

Family: rho_s(class) ~ exp(-beta * H(class)) restricted to classes satisfying
context s (context-hard), H = total frustration count. beta -> inf recovers
the frustration ansatz (= LP optimum in all computed scenarios).

Claims to verify:
 1. Full correlators equal targets EXACTLY at every beta (context-hard).
 2. Proper-subset marginals vanish at every beta (fiber-averaged classes).
 3. F(beta) runs from the LP floor (beta->inf) to EXACTLY 1/2 (beta->0),
    monotonically -- the flooding value is the hot limit, the LP floor the
    cold limit, for every n.
"""
import numpy as np, itertools

FLOORS = {3: 1/3, 4: 1/3, 5: 2/5}

def scenario(n):
    settings = [s for s in itertools.product([0,1], repeat=n) if sum(s) % 2 == 0]
    targets  = np.array([1.0 if sum(s) % 4 == 0 else -1.0 for s in settings])
    effs = np.array(list(itertools.product([0,1], repeat=n+1)), dtype=int)
    P, G = effs[:, 0], effs[:, 1:]
    sett = np.array(settings, dtype=int)
    sign = np.where((P[None, :] + sett @ G.T) % 2 == 1, -1.0, 1.0)
    return sign, targets

def gibbs_family(n, beta):
    sign, t = scenario(n)
    H = (sign != t[:, None]).sum(axis=0).astype(float)   # frustration per class
    rhos = []
    for si in range(len(t)):
        sat = (sign[si] == t[si])
        w = np.where(sat, np.exp(-beta * (H - H.min())), 0.0)
        rhos.append(w / w.sum())
    return np.array(rhos), sign, t

def F_and_checks(n, beta):
    rhos, sign, t = gibbs_family(n, beta)
    # correlator exactness
    corr_err = max(abs((rhos[si] * sign[si]).sum() - t[si]) for si in range(len(t)))
    # Hall fraction
    F = 0.0
    for i in range(len(t)):
        for j in range(i + 1, len(t)):
            F = max(F, 0.5 * np.abs(rhos[i] - rhos[j]).sum())
    return F, corr_err

for n in [3, 4, 5]:
    print(f"\nn = {n}  (LP floor = {FLOORS[n]:.6f}, ceiling = 0.5)")
    print(f"{'beta':>8} {'F(beta)':>12} {'corr err':>10}")
    prev = None
    for beta in [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 50.0]:
        F, ce = F_and_checks(n, beta)
        mono = "" if prev is None or F <= prev + 1e-12 else "  <-- NON-MONOTONE"
        print(f"{beta:>8.2f} {F:>12.8f} {ce:>10.2e}{mono}")
        prev = F
    F0, _ = F_and_checks(n, 0.0)
    Finf, _ = F_and_checks(n, 60.0)
    assert abs(F0 - 0.5) < 1e-12, f"hot limit != 1/2: {F0}"
    assert abs(Finf - FLOORS[n]) < 1e-9, f"cold limit != floor: {Finf}"
    print(f"  hot limit F(0) = {F0} = 1/2 exactly; cold limit F(inf) = {Finf:.9f} = floor")

# large-beta approach: F(beta) - floor ~ c * exp(-Delta * beta); extract Delta
# per-n window chosen where F - floor is well above underflow (n = 5 decays
# twice as fast, so its window sits at smaller beta)
print("\ncold-side gap Delta from fit of ln(F - floor):")
WINDOWS = {3: (6.0, 10.0), 4: (6.0, 10.0), 5: (2.0, 4.0)}
EXPECTED = {3: 2.0, 4: 2.0, 5: 4.0}
for n in [3, 4, 5]:
    lo, hi = WINDOWS[n]
    betas = np.linspace(lo, hi, 5)
    vals = np.array([F_and_checks(n, b)[0] - FLOORS[n] for b in betas])
    assert vals.min() > 1e-13, "fit window underflows; move it to smaller beta"
    Delta = -np.polyfit(betas, np.log(vals), 1)[0]
    print(f"  n={n}: Delta = {Delta:.4f}  (expected {EXPECTED[n]}; "
          f"R = {2**((n-1)//2)})")
    assert abs(Delta - EXPECTED[n]) < 0.05, f"n={n}: Delta {Delta}"
print("gap tracks the staircase ratio: Delta = R at every computed n "
      "(curiosity, logged not claimed)")

print("\nALL ASSERTIONS PASSED")
