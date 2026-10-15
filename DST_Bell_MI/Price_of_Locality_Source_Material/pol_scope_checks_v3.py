"""
pol_scope_checks_v3.py (v3: asserts the quoted numbers -- Mermin-5 minus one: LP 2/5, B = 4/11, 180 cores of sizes 4 (60) and 6 (120), construction feasible with F = 1/2 on 12 ground classes; C4 core: LP 1/3, s = 3/4, construction faithful and optimal)
pol_scope_checks_v2.py (v2: part (C) marginals checked in integer arithmetic; dead assertion removed; prints flushed)
pol_scope_checks_v1.py -- direct computations for the three scope questions
raised by the v6 cold audit of arXiv:2608.00124 (Price of Locality).

  (D) Mermin-5 minus one constraint.  v5 reports: LP minimum stays 2/5,
      Theorem 3 bound B falls to 4/11.  Questions: does the ground-state
      construction attain 2/5 there (Conjecture 2 on a non-core scenario)?
      Does any embedded core have minimum 2/5 (the X.B "inherit from the
      tightest embedded core" clause)?  All three deletion orbits (|S| = 0,
      2, 4) are run.

  (C) C4 linear-cluster AVN core {IZXZ, IZYY, ZYYZ, ZYXY}, signs (+,+,+,-),
      alone and padded with two further stabilizers of the C4 group (all 55
      pairs).  Strategy-level faithful LP with every marginal correlator of
      the measured observables constrained to its stabilizer value (+-1 on
      the group, 0 off it).  Questions: is the ground-state construction
      feasible (does it satisfy the marginals) and optimal on the core and
      on each padded scenario?  What is s for each?

  (B) is implicit in (C): "faithful" for a cluster scenario means every
      marginal correlator at its quantum value, not "vanishing".

LP values are double precision (HiGHS) and reported with their rational
reconstruction; construction values are exact Fractions of integer counts.
Theorem 3 bounds are exact.

Usage: python pol_scope_checks_v1.py
"""
import itertools, sys, time, functools
print = functools.partial(print, flush=True)
from fractions import Fraction
from math import comb
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


# --------------------------------------------------------------------------
# generic faithful LP: minimise max pairwise L1 distance between densities
# --------------------------------------------------------------------------
def solve_minmax_l1(nvar_per_setting, settings_eq):
    """settings_eq: list over settings of (A_s, b_s) with A_s dense/sparse
    (rows x nvar_per_setting).  Variables: rho_s (one block per setting),
    u_{ss'} (one block per pair), T.  Returns (T/2, x)."""
    S = len(settings_eq); V = nvar_per_setting
    pairs = list(itertools.combinations(range(S), 2))
    nrho = S * V; nu = len(pairs) * V; Tvar = nrho + nu; nvar = Tvar + 1
    er, ec, ev, beq = [], [], [], []
    row = 0
    for si, (A, b) in enumerate(settings_eq):
        A = coo_matrix(A)
        er += (A.row + row).tolist(); ec += (A.col + si * V).tolist(); ev += A.data.tolist()
        beq += list(b); row += A.shape[0]
    A_eq = coo_matrix((ev, (er, ec)), shape=(row, nvar)).tocsr()
    ur, uc, uv, bub = [], [], [], []
    row = 0
    for pi, (i, j) in enumerate(pairs):
        base_u = nrho + pi * V
        for k in range(V):
            for sgn in (1.0, -1.0):
                ur += [row] * 3; uc += [i * V + k, j * V + k, base_u + k]
                uv += [sgn, -sgn, -1.0]; bub.append(0.0); row += 1
        ur += [row] * (V + 1); uc += list(range(base_u, base_u + V)) + [Tvar]
        uv += [1.0] * V + [-1.0]; bub.append(0.0); row += 1
    A_ub = coo_matrix((uv, (ur, uc)), shape=(row, nvar)).tocsr()
    c = np.zeros(nvar); c[Tvar] = 1.0
    res = linprog(c, A_ub=A_ub, b_ub=np.array(bub), A_eq=A_eq, b_eq=np.array(beq),
                  bounds=(0, None), method="highs")
    assert res.success, res.message
    return res.x[Tvar] / 2.0, res.x


def frac(x, den=400):
    return Fraction(x).limit_denominator(den)


# --------------------------------------------------------------------------
# (D) Mermin-5 minus one constraint, class level (Theorem 1 applies to any
#     subset of Mermin settings: the proof uses only the Mermin functionals)
# --------------------------------------------------------------------------
def mermin_minus_one(n=5):
    print(f"\n=== (D) Mermin-{n} minus one constraint ===")
    masks = [m for m in range(1 << n) if bin(m).count("1") % 2 == 0]
    tgt = {m: (0 if bin(m).count("1") % 4 == 0 else 1) for m in masks}
    classes = [(P, g) for P in (0, 1) for g in range(1 << n)]
    def sat(m, c):
        P, g = c
        return (P ^ (bin(g & m).count("1") % 2)) == tgt[m]
    for k_del in sorted(set(bin(m).count("1") for m in masks)):
        s_del = next(m for m in masks if bin(m).count("1") == k_del)
        sett = [m for m in masks if m != s_del]
        N = len(sett)
        satmat = np.array([[sat(m, c) for c in classes] for m in sett])   # [setting, class]
        # LP
        eqs = []
        for si, m in enumerate(sett):
            A = np.zeros((2, len(classes))); A[0, :] = 1.0
            A[1, :] = np.where(satmat[si], 1.0, -1.0); eqs.append((A, [1.0, 1.0]))
        F_lp, _ = solve_minmax_l1(len(classes), eqs)
        # Theorem 3 bound B over ALL subsets (exact), and over cores only
        satbits = [sum(1 << si for si in range(N) if satmat[si, ci]) for ci in range(len(classes))]
        B = Fraction(0); B_sub = None
        vec = {m: (1 << n) | m for m in sett}
        xor_v = np.array([vec[m] for m in sett]); tg = np.array([tgt[m] for m in sett])
        sb_arr = np.array(satbits, dtype=np.int64)
        POP = np.array([bin(i).count("1") for i in range(1 << 16)])
        allmasks = np.arange(1 << N, dtype=np.int64)
        r_of = POP[allmasks]
        # m(S') for every subset mask, vectorised over classes
        mS_all = np.zeros(1 << N, dtype=np.int64)
        for sb in satbits:
            mS_all = np.maximum(mS_all, POP[allmasks & sb])
        valid = r_of >= 2
        num = (r_of - mS_all); den = np.maximum(r_of - 1, 1)
        # best bound: compare fractions exactly via cross-multiplication
        best = None
        for mask in np.nonzero(valid)[0]:
            v = Fraction(int(num[mask]), int(den[mask]))
            if best is None or v > best[0]: best = (v, int(mask))
        B, B_sub = best
        # contradictions: XOR of functionals 0, XOR of targets 1, unsatisfiable (mS = r-1)
        xr = np.zeros(1 << N, dtype=np.int64); tr = np.zeros(1 << N, dtype=np.int64)
        for si in range(N):
            bit = allmasks >> si & 1
            xr ^= bit * int(xor_v[si]); tr ^= bit * int(tg[si])
        contr = valid & (xr == 0) & (tr == 1)
        # has[mask] = some submask (incl. itself) is a contradiction
        has = contr.copy()
        for si in range(N):
            bit = 1 << si
            idx = np.nonzero(allmasks & bit)[0]
            has[idx] |= has[idx ^ bit]
        minimal = contr.copy()
        for si in range(N):
            bit = 1 << si
            idx = np.nonzero(allmasks & bit)[0]
            minimal[idx] &= ~has[idx ^ bit]
        cores = np.nonzero(minimal)[0]
        core_sizes = sorted(set(r_of[cores].tolist()))
        best_core = (Fraction(0), None)
        for mask in cores:
            v = Fraction(int(num[mask]), int(den[mask]))
            if v > best_core[0]: best_core = (v, int(mask))
        B_sub = [si for si in range(N) if B_sub >> si & 1]
        best_core = (best_core[0], [si for si in range(N) if best_core[1] >> si & 1] if best_core[1] is not None else None)
        # ground-state construction, exact
        H = (~satmat).sum(axis=0)
        ground = np.nonzero(H == H.min())[0]
        G = satmat[:, ground]
        g = G.sum(axis=1)
        feasible = bool((g > 0).all())
        Fc = None
        if feasible:
            Fc = Fraction(0)
            for i, j in itertools.combinations(range(N), 2):
                # L1/2 between uniform-on-G_i and uniform-on-G_j (may have different g)
                inter = int((G[i] & G[j]).sum())
                # 1/2 sum |1/g_i 1_{G_i} - 1/g_j 1_{G_j}|
                gi, gj = int(g[i]), int(g[j])
                d = Fraction(inter, 1) * abs(Fraction(1, gi) - Fraction(1, gj)) \
                    + Fraction(gi - inter, gi) + Fraction(gj - inter, gj)
                Fc = max(Fc, d / 2)
        assert abs(F_lp - 0.4) < 1e-7 and B == Fraction(4, 11) and len(cores) == 180 and core_sizes == [4, 6], (F_lp, B, len(cores), core_sizes)
        assert sum(1 for c in cores if r_of[c] == 4) == 60 and feasible and Fc == Fraction(1, 2) and len(ground) == 12 and int(H.min()) == 5
        print(f"deleted |S|={k_del}: N={N}  LP F_min = {F_lp:.9f} = {frac(F_lp)}   "
              f"Theorem-3 B over all subsets = {B} (|S'|={len(B_sub)})   "
              f"cores: {len(cores)} (sizes {core_sizes}), best core bound = {best_core[0]} (|core|={len(best_core[1]) if best_core[1] else None})   "
              f"construction: feasible={feasible}, F = {Fc}, |ground|={len(ground)}, Hmin={int(H.min())}, "
              f"sector sizes={sorted(set(g.tolist()))}")


# --------------------------------------------------------------------------
# (C) C4 linear cluster, strategy level
# --------------------------------------------------------------------------
PAULI = {'I': 0, 'X': 1, 'Y': 2, 'Z': 3}
def pmul(a, b):
    """single-qubit Pauli product: returns (phase, pauli) with phase in {1,-1,1j,-1j}"""
    if a == 0: return 1, b
    if b == 0: return 1, a
    if a == b: return 1, 0
    # XY = iZ, YZ = iX, ZX = iY
    cyc = {(1, 2): (1j, 3), (2, 3): (1j, 1), (3, 1): (1j, 2)}
    if (a, b) in cyc: return cyc[(a, b)]
    ph, p = cyc[(b, a)]; return -ph, p

def smul(s1, s2):
    """product of signed Pauli strings (sign, tuple)"""
    sign = s1[0] * s2[0]; out = []
    for a, b in zip(s1[1], s2[1]):
        ph, p = pmul(a, b); sign *= ph; out.append(p)
    return sign, tuple(out)

def stabilizer_group(gens):
    group = {}
    for bits in itertools.product((0, 1), repeat=len(gens)):
        el = (1, tuple(0 for _ in gens[0][1]))
        for b, gq in zip(bits, gens):
            if b: el = smul(el, gq)
        assert el[0] in (1, -1), el
        group[el[1]] = int(el[0].real) if isinstance(el[0], complex) else el[0]
    return group

def pstr(t): return ''.join('IXYZ'[p] for p in t)

def c4_scenario(extra_pairs=None):
    print("\n=== (C) C4 linear cluster: AVN core and padded scenarios ===")
    n = 4
    gens = [(1, (1, 3, 0, 0)), (1, (3, 1, 3, 0)), (1, (0, 3, 1, 3)), (1, (0, 0, 3, 1))]
    grp = stabilizer_group(gens)
    core = [('IZXZ', +1), ('IZYY', +1), ('ZYYZ', +1), ('ZYXY', -1)]
    for name, sg in core:
        t = tuple(PAULI[ch] for ch in name)
        assert t in grp and grp[t] == sg, (name, grp.get(t))
    others = [pstr(t) for t in grp if t != (0, 0, 0, 0) and pstr(t) not in [c for c, _ in core]]
    print("C4 stabilizer group (signed):", ', '.join(('+' if grp[t] > 0 else '-') + pstr(t) for t in sorted(grp)))
    print("core verified in group with signs", core)

    # strategies: per party outcomes (x, y, z) in {+-1}^3 -> 8 per party, 4096 total
    strat = np.array(list(itertools.product(*[list(itertools.product((1, -1), repeat=3))] * n)))  # [4096, 4, 3]

    def run(settings):
        """settings: list of (name, sign).  returns LP F, s, construction feasibility & F."""
        S = len(settings)
        eqs = []; full_sat = []
        for name, sg in settings:
            t = tuple(PAULI[ch] for ch in name)
            parties = [i for i in range(n) if t[i] != 0]
            rows = [np.ones(len(strat))]; rhs = [1.0]
            for r in range(1, len(parties) + 1):
                for T in itertools.combinations(parties, r):
                    sub = tuple(t[i] if i in T else 0 for i in range(n))
                    qv = grp.get(sub, 0)        # stabilizer value: +-1 on the group, 0 off it
                    prod = np.ones(len(strat))
                    for i in T: prod = prod * strat[:, i, t[i] - 1]
                    rows.append(prod); rhs.append(float(qv))
                    if r == len(parties):
                        assert qv == sg; full_sat.append(prod == sg)
            eqs.append((np.array(rows), rhs))
        F_lp, _ = solve_minmax_l1(len(strat), eqs)
        satmat = np.array(full_sat)                      # [setting, strategy]
        H = (~satmat).sum(axis=0)
        s_val = Fraction(S - int(H.min()), S)
        ground = np.nonzero(H == H.min())[0]
        G = satmat[:, ground]; g = G.sum(axis=1)
        feasible_support = bool((g > 0).all())
        # check marginals of the construction exactly
        marg_ok = True; Fc = None
        if feasible_support:
            for si, (A, b) in enumerate(eqs):
                idx = ground[G[si]]
                # exact: each row of A is +-1 (or 1 for normalization); sum over support must equal b * g
                for row, bb in zip(A, b):
                    if int(np.rint(row[idx]).sum()) != int(bb) * int(g[si]): marg_ok = False
            Fc = Fraction(0)
            for i, j in itertools.combinations(range(S), 2):
                inter = int((G[i] & G[j]).sum()); gi, gj = int(g[i]), int(g[j])
                d = Fraction(inter) * abs(Fraction(1, gi) - Fraction(1, gj)) + Fraction(gi - inter, gi) + Fraction(gj - inter, gj)
                Fc = max(Fc, d / 2)
        return F_lp, s_val, feasible_support, marg_ok, Fc, len(ground)

    F_lp, s_val, fs, mo, Fc, ng = run(core)
    assert abs(F_lp - 1/3) < 1e-7 and s_val == Fraction(3, 4) and fs and mo and Fc == Fraction(1, 3), (F_lp, s_val, fs, mo, Fc)
    print(f"core alone: LP F_min = {F_lp:.9f} = {frac(F_lp)}, s = {s_val}, F*s = {frac(F_lp)*s_val};  "
          f"construction: support-feasible={fs}, marginals satisfied={mo}, F = {Fc}, |ground|={ng}")

    if '--padded' not in sys.argv:
        print("(the 55 padded scenarios are certified by c4_padding_construction_v2.py and c4_bad_pairs_lp_v2.py; pass --padded to also LP all 55 here, ~5 min each)")
        print("ASSERTED: Mermin-5 minus one (all three orbits) and the C4 core match the quoted values")
        return
    results = []
    for a, b in itertools.combinations(others, 2):
        sett = core + [(a, grp[tuple(PAULI[ch] for ch in a)]), (b, grp[tuple(PAULI[ch] for ch in b)])]
        F_lp, s_val, fs, mo, Fc, ng = run(sett)
        results.append((a, b, frac(F_lp), s_val, fs, mo, Fc, ng))
    print(f"\ncore + 2 further stabilizers, all {len(results)} pairs:")
    print(f"{'pad1':>6} {'pad2':>6} {'LP F':>6} {'s':>6} {'F*s':>6}  support  marg  constr F  |ground|")
    for a, b, Fl, sv, fs, mo, Fc, ng in results:
        print(f"{a:>6} {b:>6} {str(Fl):>6} {str(sv):>6} {str(Fl*sv):>6}  {str(fs):>7}  {str(mo):>4}  {str(Fc):>8}  {ng}")
    s56 = [r for r in results if r[3] == Fraction(5, 6)]
    print(f"\npairs with s = 5/6: {len(s56)}; of these LP F = 1/3: {sum(1 for r in s56 if r[2] == Fraction(1,3))}, "
          f"construction feasible+optimal: {sum(1 for r in s56 if r[4] and r[5] and r[6] == Fraction(1,3))}")


if __name__ == "__main__":
    t0 = time.time()
    mermin_minus_one(5)
    c4_scenario()
    print(f"\n[{time.time()-t0:.1f}s]")
