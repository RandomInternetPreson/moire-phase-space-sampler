"""
ultrametric_check.py -- machine verification of the landscape-geometry claims
of "The glassy landscape" subsection: no bit-hierarchy ultrametric grades the
violation shells at n = 4 (0 of 120 orderings, exhaustive), and Hamming
proximity also fails (maximal-violation classes at Hamming distance one from
the ground set).

Definitions (stated precisely so the check is self-contained):

  Class space: (P, gamma) in GF(2)^(n+1), the fiber-reduced strategy classes.
  Frustration H(P, gamma): number of violated Mermin constraints (of N = 8
  at n = 4). The violation shells are the distinct values of H.

  Bit-hierarchy ultrametric: choose an ordering pi of the n+1 bits. This
  induces the standard dyadic tree on GF(2)^(n+1): two classes lie in the
  same depth-k subtree iff they agree on the first k bits under pi, and
  d_pi(c, c') = 2^(-k) with k the length of the longest common prefix.
  Distance to the ground set: D_pi(c) = min over ground classes g of
  d_pi(c, g), with D_pi = 0 on the ground set itself.

  "Grades the violation shells": H is a strictly increasing function of
  D_pi -- i.e. each distance level of D_pi is H-homogeneous and H increases
  with distance. This is the property required for frustration to be
  readable as ultrametric distance from an invariant set.

  There are (n+1)! = 120 bit orderings at n = 4; the check is exhaustive.

  A strictly weaker property -- monotonicity only (min H at each level not
  below the running max of previous levels, with mixed levels allowed) --
  is also reported for transparency. Monotone-only orderings do NOT grade
  the shells, since distance then fails to determine frustration.

Assertions certify: shells = {2, 4, 6}; exact-grading orderings = 0 of 120;
some maximal-violation class at Hamming distance 1 from the ground set
(the parity flip P -> P xor 1 maps H -> 8 - H, so a ground class's parity
flip is maximally violating at distance one).

Runtime: < 5 s. No dependencies beyond the standard library.
"""
import itertools
from collections import defaultdict

n = 4
settings = [s for s in itertools.product([0, 1], repeat=n) if sum(s) % 2 == 0]
targets = [1 if sum(s) % 4 == 0 else -1 for s in settings]
N = len(settings)
assert N == 8

classes = list(itertools.product([0, 1], repeat=n + 1))   # (P, g1..gn)


def H(c):
    P, g = c[0], c[1:]
    h = 0
    for s, t in zip(settings, targets):
        sign = (-1) ** ((P + sum(a * b for a, b in zip(g, s))) % 2)
        if sign != t:
            h += 1
    return h


Hs = {c: H(c) for c in classes}
shells = sorted(set(Hs.values()))
Hmin, Hmax = min(shells), max(shells)
ground = [c for c in classes if Hs[c] == Hmin]

print(f"n = {n}: N = {N} constraints, {len(classes)} classes")
print(f"violation shells: {shells}   Hmin = {Hmin}, |ground| = {len(ground)}")
assert shells == [2, 4, 6], "shell structure changed -- claim must be re-derived"
assert len(ground) == 8

# ---- Hamming proximity failure -------------------------------------------
maxcls = [c for c in classes if Hs[c] == Hmax]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


min_d_max_to_ground = min(hamming(mc, g) for mc in maxcls for g in ground)
print(f"minimum Hamming distance of a maximal-violation class to the ground "
      f"set: {min_d_max_to_ground}")
assert min_d_max_to_ground == 1, "Hamming-proximity claim failed"
# the witness is structural: flipping P alone maps H -> N - H
for g in ground:
    pf = (g[0] ^ 1,) + g[1:]
    assert Hs[pf] == N - Hmin == Hmax

# ---- bit-hierarchy ultrametric exhaustion --------------------------------


def levels_for_order(order):
    """Return {distance level: [H values]} for D_pi under bit ordering."""
    lev = defaultdict(list)
    gset = set(ground)
    for c in classes:
        if c in gset:
            d = 0.0
        else:
            best_k = -1
            for g in ground:
                k = 0
                for b in order:
                    if c[b] == g[b]:
                        k += 1
                    else:
                        break
                best_k = max(best_k, k)
            d = 2.0 ** (-best_k)
        lev[d].append(Hs[c])
    return lev


def grades(order):
    """Exact grading: every level H-homogeneous, strictly increasing in distance."""
    lev = levels_for_order(order)
    prev = -1
    for d in sorted(lev):
        hs = set(lev[d])
        if len(hs) > 1:
            return False
        h = hs.pop()
        if h <= prev:
            return False
        prev = h
    return True


def monotone_only(order):
    lev = levels_for_order(order)
    running = -1
    for d in sorted(lev):
        if min(lev[d]) < running:
            return False
        running = max(running, max(lev[d]))
    return True


orders = list(itertools.permutations(range(n + 1)))
assert len(orders) == 120
n_grade = sum(grades(o) for o in orders)
n_mono = sum(monotone_only(o) for o in orders)

print(f"bit orderings that GRADE the shells (exact): {n_grade} of {len(orders)}")
print(f"  (weaker, monotone-only, mixed levels allowed: {n_mono} of "
      f"{len(orders)} -- these do not constitute a grading)")
assert n_grade == 0, "an ultrametric grading exists -- claim falsified"

print("\nALL ASSERTIONS PASSED: no bit-hierarchy ultrametric grades the "
      "violation shells (0 of 120, exhaustive); maximal-violation classes "
      "lie at Hamming distance one from the ground set.")
