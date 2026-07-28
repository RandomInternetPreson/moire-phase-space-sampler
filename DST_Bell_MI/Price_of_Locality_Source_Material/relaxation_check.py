"""
relaxation_check.py -- machine verification of the relaxation-dynamics claims
of "The glassy landscape" subsection: single-site relaxation dynamics on the
class-space landscape are non-ergodic on constraint sectors, so any dynamical
account requires collective (>= 2-site) moves; with collective moves, a local
stochastic automaton relaxes each context to the uniform contextual ground
shell and regenerates the exact floor F = 1/3 at n = 3 to sampling accuracy
[TOY].

The check has two parts.

PART A (exact, no randomness). On the fiber-reduced class space GF(2)^(n+1)
with frustration energy H = number of violated Mermin constraints:

  * The constraint sector of a class is the set of constraints it violates.
  * Zero-temperature single-site dynamics accepts a one-bit flip iff H does
    not increase. Non-ergodicity claim: the ground set decomposes into
    multiple dynamically disconnected components under this rule.
      n = 4: every one-bit flip from a ground class raises H by >= 2, so all
             8 ground classes are isolated fixed points (8 components) and
             the violated sector is strictly conserved: single-site dynamics
             can never move between constraint sectors.
      n = 3: 3 components; the two ground classes violating the all-X
             constraint are isolated singletons, unreachable from the other
             sector.
  * Collective moves repair this: 2-site flips create equal-energy,
    sector-changing moves inside the ground set (components drop), and
    3-site flips connect the ground set completely at both sizes.

PART B (seeded simulation). Contextual Metropolis relaxation at n = 3: for
each setting s the energy is E_s = H + kappa * [constraint s violated]
(kappa = 3), T = 0.30, proposals drawn from 1- and 2-site flips. The chain
relaxes to the uniform measure on the contextual ground shell
(ground classes satisfying s) -- the frustration-ansatz density of the
paper. Certified: per-context spread from uniformity <= 0.004 (typical run
<= 0.003), and the pairwise variational distance of the empirical contextual
densities gives F within 0.01 of the exact value 1/3. A zero-temperature
single-site control run visits only its own component's classes (6 of 8
ground classes from the generic start), exhibiting the frozen sectors of
Part A directly.

Runtime: ~1-2 min. Requires numpy.
"""
import itertools
import numpy as np
from collections import deque

# ---------------------------------------------------------------------------
# Common scenario machinery
# ---------------------------------------------------------------------------


def scenario(n):
    settings = [s for s in itertools.product([0, 1], repeat=n)
                if sum(s) % 2 == 0]
    targets = [1 if sum(s) % 4 == 0 else -1 for s in settings]
    classes = list(itertools.product([0, 1], repeat=n + 1))

    def satisfies(c, si):
        P, g = c[0], c[1:]
        s, t = settings[si], targets[si]
        return ((-1) ** ((P + sum(a * b for a, b in zip(g, s))) % 2)) == t

    def sector(c):
        return tuple(si for si in range(len(settings)) if not satisfies(c, si))

    H = {c: len(sector(c)) for c in classes}
    return settings, classes, H, satisfies, sector


def components(nodes, adjacency):
    seen, comps = set(), []
    for v in nodes:
        if v in seen:
            continue
        comp, dq = {v}, deque([v])
        while dq:
            u = dq.popleft()
            for w in adjacency(u):
                if w in comp:
                    continue
                comp.add(w)
                dq.append(w)
        seen |= comp
        comps.append(comp)
    return comps


def k_flips(c, k, nbits):
    for idx in itertools.combinations(range(nbits), k):
        l = list(c)
        for i in idx:
            l[i] ^= 1
        yield tuple(l)


# ---------------------------------------------------------------------------
# PART A: exact non-ergodicity / collectivity structure
# ---------------------------------------------------------------------------
print("PART A: exact zero-temperature structure")
expected_components = {3: {1: 3, 2: 2, 3: 1}, 4: {1: 8, 2: 3, 3: 1}}

for n in (3, 4):
    settings, classes, H, satisfies, sector = scenario(n)
    Hmin = min(H.values())
    ground = [c for c in classes if H[c] == Hmin]
    gset = set(ground)
    nbits = n + 1

    print(f"\nn = {n}: N = {len(settings)}, Hmin = {Hmin}, "
          f"|ground| = {len(ground)}")

    # barrier of single-site moves out of the ground set
    dH1 = sorted({H[x] - Hmin for g in ground for x in k_flips(g, 1, nbits)})
    print(f"  single-site dH values from ground: {dH1}")
    if n == 4:
        assert min(dH1) == 2, "n=4 single-site freeze-out failed"
        # sector conservation: no single flip within the ground set at all
        assert all(x not in gset
                   for g in ground for x in k_flips(g, 1, nbits))

    for k in (1, 2, 3):
        def adj(u, k=k):
            return [x for x in k_flips(u, k, nbits) if x in gset]
        comps = components(ground, adj)
        mixes = any(sector(a) != sector(b)
                    for a in ground for b in k_flips(a, k, nbits)
                    if b in gset and b != a)
        print(f"  k = {k}-site moves at T = 0: {len(comps)} ground "
              f"component(s); sector-changing in-ground move exists: {mixes}")
        assert len(comps) == expected_components[n][k], \
            f"component structure changed at n={n}, k={k}"
        if k == 1:
            assert len(comps) > 1, "single-site dynamics is ergodic -- " \
                                   "claim falsified"
        if k >= 2 and n == 4:
            assert mixes, "collective moves fail to change sectors"

print("\n  => single-site zero-temperature dynamics is non-ergodic on "
      "constraint sectors at n = 3 and n = 4 (exact); collective moves "
      "restore connectivity.")

# ---------------------------------------------------------------------------
# PART B: contextual Metropolis relaxation at n = 3
# ---------------------------------------------------------------------------
print("\nPART B: contextual Metropolis relaxation (n = 3, seeded)")
n = 3
settings, classes, H, satisfies, sector = scenario(n)
N = len(settings)
Hmin = min(H.values())
ground = [c for c in classes if H[c] == Hmin]
nbits = n + 1

rng = np.random.default_rng(20260728)
KAPPA, T = 3.0, 0.30
STEPS, BURN = 1_200_000, 200_000


def run_context(si, move_sizes):
    def E(c):
        return H[c] + KAPPA * (0 if satisfies(c, si) else 1)

    c = next(g for g in ground if satisfies(g, si))
    occ = {}
    for step in range(STEPS):
        k = move_sizes[rng.integers(len(move_sizes))]
        idx = rng.choice(nbits, size=k, replace=False)
        cn = list(c)
        for i in idx:
            cn[i] ^= 1
        cn = tuple(cn)
        dE = E(cn) - E(c)
        if dE <= 0 or rng.random() < np.exp(-dE / T):
            c = cn
        if step >= BURN:
            occ[c] = occ.get(c, 0) + 1
    tot = sum(occ.values())
    emp = {k_: v / tot for k_, v in occ.items()}
    support = [g for g in ground if satisfies(g, si)]
    mass = sum(emp.get(g, 0.0) for g in support)
    cond = {g: emp.get(g, 0.0) / mass for g in support}
    spread = max(abs(v - 1.0 / len(support)) for v in cond.values())
    return cond, spread


conds, spreads = [], []
for si in range(N):
    cond, spread = run_context(si, move_sizes=(1, 2))
    conds.append(cond)
    spreads.append(spread)
    print(f"  context {si}: support {len(cond)}, spread from uniform "
          f"= {spread:.4f}")

F_emp = 0.0
for a, b in itertools.combinations(range(N), 2):
    keys = set(conds[a]) | set(conds[b])
    F_emp = max(F_emp, 0.5 * sum(abs(conds[a].get(k_, 0) - conds[b].get(k_, 0))
                                 for k_ in keys))
print(f"  empirical F = {F_emp:.4f}   (exact value 1/3 = {1/3:.4f})")

assert max(spreads) <= 0.004, "contextual shells not uniform"
assert abs(F_emp - 1.0 / 3.0) <= 0.01, "empirical floor off target"

# zero-temperature single-site control: frozen sectors, directly observed
c = ground[0]
visited = {c}
for _ in range(200_000):
    i = rng.integers(nbits)
    cn = list(c)
    cn[i] ^= 1
    cn = tuple(cn)
    if H[cn] <= H[c]:
        c = cn
        visited.add(c)
vg = sum(1 for v in visited if H[v] == Hmin)
print(f"  T = 0 single-site control: visited {vg} of {len(ground)} ground "
      f"classes (trapped in its component)")
assert vg < len(ground), "single-site control escaped its component"

print("\nALL ASSERTIONS PASSED: single-site dynamics non-ergodic on "
      "constraint sectors; collective (>= 2-site) moves relax each context "
      "to the uniform contextual ground shell and regenerate F = 1/3 to "
      "sampling accuracy.")
