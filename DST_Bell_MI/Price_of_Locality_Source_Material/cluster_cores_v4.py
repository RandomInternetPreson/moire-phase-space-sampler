"""
cluster_cores_v4.py -- exact certificates for the cluster-state rows of Table I
of arXiv:2608.00124 (v8), replacing LP-only evidence.

For each of the graph states C4 (linear chain), C5 (ring), C6 (ring):
  1. Build the signed stabilizer group (2^n elements).
  2. Enumerate the minimal contradictions ("cores") of sizes up to --max-size
     (default 5) among the nontrivial stabilizers: a set of signed stabilizers
     is a contradiction iff every (party, Pauli) occurrence count is even (the
     GF(2) functionals on the outcome bits sum to zero) and the signs multiply
     to -1; a core is a contradiction with no proper sub-contradiction.  This
     is a census of sizes <= max-size only; core_census_v1.py gives the complete
     census (all sizes) for C4 and C5 by kernel enumeration.
  3. For every core (or, for C6's 1,920 five-cores, every core unless --sample
     is given): take it as a scenario on its own, at strategy level (8 pre-set
     outcomes per party: one for each of X, Y, Z; in the setting of stabilizer
     g the parties in its support measure the indicated Pauli and the others
     are not measured).  Compute exactly:
       - s = (k - 1)/k and the Theorem 3 lower bound 1/(k - 1)   [exact];
       - the ground-state construction (rho_s uniform on strategies that
         violate exactly one constraint, within the sector of s): its F, and
         whether every marginal correlator over a subset T of the support
         equals +1 / -1 / 0 according to whether +g|_T / -g|_T / neither lies
         in the stabilizer group  [exact integer counts];
     If the construction is faithful with F = 1/(k - 1), the faithful minimum
     is certified exactly: lower bound Theorem 3, upper bound the construction.
  4. Report, per graph state and core size: number of cores, how many are
     certified exactly, and the F values seen.

Also prints the counts that the paper quotes: C6 has no 3-element contradiction
and 1,920 minimal 5-element contradictions.

  5. For cores on which the construction is NOT faithful, solve the faithful
     LP on the collapsed strategy space (only the Paulis each party actually
     measures in the core), then try to certify the LP optimum exactly:
     rationalize the solution (limit_denominator) and verify every constraint
     and the value in Fractions.  'LP+exact' = certified; 'LP' = double
     precision only.
Usage: python cluster_cores_v4.py [--sample K] [--max-size M] [--no-lp] [--lp-only START STRIDE]
(v4 = v3 renamed for the v2/v3 dependents; asserts the counts quoted in the paper when run WITHOUT --sample;
 with --sample the assertions are skipped and the run says so)   (K random five-cores of C6
       instead of all 1,920; default: all)
"""
import itertools, sys, time, argparse, functools
from fractions import Fraction
import numpy as np
print = functools.partial(print, flush=True)

# ---------------------------------------------------------------- Paulis --
# encode a Pauli string as (x bits, z bits) with X=(1,0), Z=(0,1), Y=(1,1)
def pauli_from_str(s):
    x = z = 0
    for i, ch in enumerate(s):
        if ch in 'XY': x |= 1 << i
        if ch in 'ZY': z |= 1 << i
    return x, z

def pauli_to_str(x, z, n):
    return ''.join('IXZY'[((x >> i) & 1) | (((z >> i) & 1) << 1)] for i in range(n))

def pauli_mul(a, b, n):
    """(sign, pauli) of a*b for signed Paulis a=(s,x,z), b=(s',x',z').  Phase
    tracked exactly; result sign in {1,-1} iff the product is Hermitian (always
    true for commuting stabilizers)."""
    sa, xa, za = a; sb, xb, zb = b
    # per-qubit phase: X*Z = -iY, Z*X = iY, X*Y = iZ, Y*X = -iZ, Y*Z = iX, Z*Y = -iX
    phase = 0  # exponent of i
    for i in range(n):
        pa = ((xa >> i) & 1, (za >> i) & 1); pb = ((xb >> i) & 1, (zb >> i) & 1)
        if pa == (0, 0) or pb == (0, 0) or pa == pb: continue
        table = {((1,0),(0,1)): -1, ((0,1),(1,0)): 1, ((1,0),(1,1)): 1, ((1,1),(1,0)): -1,
                 ((1,1),(0,1)): 1, ((0,1),(1,1)): -1}
        phase += table[(pa, pb)]
    phase %= 4
    assert phase in (0, 2), "non-Hermitian product"
    return (sa * sb * (1 if phase == 0 else -1), xa ^ xb, za ^ zb)

def stabilizer_group(n, edges):
    gens = []
    for i in range(n):
        s = ['I'] * n; s[i] = 'X'
        for (a, b) in edges:
            if a == i: s[b] = 'Z'
            if b == i: s[a] = 'Z'
        x, z = pauli_from_str(''.join(s)); gens.append((1, x, z))
    group = {}
    for bits in itertools.product((0, 1), repeat=n):
        el = (1, 0, 0)
        for bt, g in zip(bits, gens):
            if bt: el = pauli_mul(el, g, n)
        group[(el[1], el[2])] = el[0]
    assert len(group) == 2 ** n
    return group

# ------------------------------------------------------- contradictions --
def occurrence_vector(x, z, n):
    """3n-bit vector with one bit per (party, Pauli) occurring in the string."""
    v = 0
    for i in range(n):
        xi, zi = (x >> i) & 1, (z >> i) & 1
        if (xi, zi) == (1, 0): v |= 1 << (3 * i)       # X
        elif (xi, zi) == (1, 1): v |= 1 << (3 * i + 1) # Y
        elif (xi, zi) == (0, 1): v |= 1 << (3 * i + 2) # Z
    return v

def minimal_contradictions(n, group, max_size):
    """All minimal classical contradictions (all-versus-nothing sets) of size
    <= max_size among nontrivial stabilizers: every (party, Pauli) occurrence
    count is even (the GF(2) functionals on the pre-set outcome bits sum to
    zero) while the signs multiply to -1.  Minimal: no proper subset is one."""
    elems = [(sg, x, z) for (x, z), sg in group.items() if (x, z) != (0, 0)]
    vec = [occurrence_vector(x, z, n) for (sg, x, z) in elems]
    sgn = [sg for (sg, x, z) in elems]
    idx = list(range(len(elems)))
    contr = []                      # list of frozensets, minimal ones only
    found = {k: [] for k in range(2, max_size + 1)}
    for k in range(2, max_size + 1):
        for comb in itertools.combinations(idx, k):
            v = 0; s = 1
            for j in comb: v ^= vec[j]; s *= sgn[j]
            if v == 0 and s == -1:
                fs = frozenset(comb)
                if not any(c < fs for c in contr):
                    contr.append(fs); found[k].append([elems[j] for j in comb])
    return found

# ------------------------------------------------------------ scenario --
def certify_core(n, group, core):
    """core: list of signed stabilizers (sign, x, z).  Returns dict."""
    k = len(core)
    S = 1 << (3 * n)                              # strategies: 3 bits per party (X,Y,Z outcomes)
    strat = np.arange(S, dtype=np.int64)
    def bit(party, pauli):                        # pauli index 0=X,1=Y,2=Z
        return (strat >> (3 * party + pauli)) & 1
    def product_parity(x, z):
        """parity bit of the product of outcomes over the support of Pauli (x,z)"""
        p = np.zeros(S, dtype=np.int64)
        for i in range(n):
            xi, zi = (x >> i) & 1, (z >> i) & 1
            if (xi, zi) == (1, 0): p ^= bit(i, 0)
            elif (xi, zi) == (1, 1): p ^= bit(i, 1)
            elif (xi, zi) == (0, 1): p ^= bit(i, 2)
        return p
    sat = []
    for (sg, x, z) in core:
        tb = 0 if sg == 1 else 1
        sat.append(product_parity(x, z) == tb)
    sat = np.array(sat)                           # [k, S]
    H = (~sat).sum(axis=0)
    s = Fraction(k - int(H.min()), k)
    lower = Fraction(k - (k - int(H.min())), k - 1)
    ground = np.nonzero(H == H.min())[0]
    G = sat[:, ground]; g = G.sum(axis=1)
    feasible = bool((g > 0).all())
    faithful = feasible
    # marginals, exact: for each setting, each proper subset T of its support
    if feasible:
        for si, (sg, x, z) in enumerate(core):
            supp = [i for i in range(n) if ((x | z) >> i) & 1]
            idx = ground[G[si]]
            for r in range(1, len(supp)):
                for T in itertools.combinations(supp, r):
                    mask = sum(1 << i for i in T)
                    xT, zT = x & mask, z & mask
                    qv = group.get((xT, zT), 0)   # +1/-1 if in group, 0 if not
                    par = product_parity(xT, zT)[idx]
                    val = int(len(idx)) - 2 * int(par.sum())   # sum of +-1 outcomes
                    if val != qv * int(len(idx)): faithful = False
    F = None
    if feasible:
        F = Fraction(0)
        for i, j in itertools.combinations(range(k), 2):
            inter = int((G[i] & G[j]).sum()); gi, gj = int(g[i]), int(g[j])
            d = Fraction(inter) * abs(Fraction(1, gi) - Fraction(1, gj)) + Fraction(gi - inter, gi) + Fraction(gj - inter, gj)
            F = max(F, d / 2)
    return dict(k=k, s=s, lower=lower, feasible=feasible, faithful=faithful, F=F,
                exact=(feasible and faithful and F == lower), ground=len(ground))



# ------------------------------------------------------ faithful LP path --
def faithful_lp_core(n, group, core, max_den=5000):
    """Faithful LP on the collapsed strategy space; returns (F_lp, F_exact or None)."""
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    k = len(core)
    measured = [set() for _ in range(n)]
    for (sg, x, z) in core:
        for i in range(n):
            xi, zi = (x >> i) & 1, (z >> i) & 1
            if (xi, zi) == (1, 0): measured[i].add(0)
            elif (xi, zi) == (1, 1): measured[i].add(1)
            elif (xi, zi) == (0, 1): measured[i].add(2)
    slots = [(i, p) for i in range(n) for p in sorted(measured[i])]
    pos = {sp: j for j, sp in enumerate(slots)}
    S = 1 << len(slots); strat = np.arange(S, dtype=np.int64)
    def parity(x, z):
        p = np.zeros(S, dtype=np.int64)
        for i in range(n):
            xi, zi = (x >> i) & 1, (z >> i) & 1
            pp = {(1,0):0,(1,1):1,(0,1):2}.get((xi, zi))
            if pp is not None: p ^= (strat >> pos[(i, pp)]) & 1
        return p
    # constraint rows per setting: normalization, full product, every proper marginal
    rows_per = []
    for (sg, x, z) in core:
        supp = [i for i in range(n) if ((x | z) >> i) & 1]
        rows = [(np.ones(S), 1)]
        for r in range(1, len(supp) + 1):
            for T in itertools.combinations(supp, r):
                mask = sum(1 << i for i in T); xT, zT = x & mask, z & mask
                qv = sg if r == len(supp) else group.get((xT, zT), 0)
                rows.append((1 - 2 * parity(xT, zT), qv))
        rows_per.append(rows)
    pairs = list(itertools.combinations(range(k), 2))
    nrho = k * S; nu = len(pairs) * S; Tv = nrho + nu; nvar = Tv + 1
    er, ec, ev, beq = [], [], [], []; row = 0
    for si, rows in enumerate(rows_per):
        for vec, b in rows:
            nz = np.nonzero(vec)[0]
            er += [row] * len(nz); ec += (nz + si * S).tolist(); ev += vec[nz].tolist(); beq.append(float(b)); row += 1
    A_eq = coo_matrix((ev, (er, ec)), shape=(row, nvar)).tocsr()
    ur, uc, uv, bub = [], [], [], []; row = 0
    for pi, (i, j) in enumerate(pairs):
        base = nrho + pi * S
        for c in range(S):
            for sgn in (1.0, -1.0):
                ur += [row] * 3; uc += [i * S + c, j * S + c, base + c]; uv += [sgn, -sgn, -1.0]; bub.append(0.0); row += 1
        ur += [row] * (S + 1); uc += list(range(base, base + S)) + [Tv]; uv += [1.0] * S + [-1.0]; bub.append(0.0); row += 1
    A_ub = coo_matrix((uv, (ur, uc)), shape=(row, nvar)).tocsr()
    cvec = np.zeros(nvar); cvec[Tv] = 1.0
    res = linprog(cvec, A_ub=A_ub, b_ub=np.array(bub), A_eq=A_eq, b_eq=np.array(beq), bounds=(0, None), method="highs")
    assert res.success, res.message
    F_lp = res.x[Tv] / 2
    # exact verification of a rationalized solution
    rho = [[Fraction(float(v)).limit_denominator(max_den) for v in res.x[si * S:(si + 1) * S]] for si in range(k)]
    ok = True
    for si, rows in enumerate(rows_per):
        for vec, b in rows:
            tot = sum(Fraction(int(vv)) * r for vv, r in zip(vec, rho[si]) if r != 0)
            if tot != b: ok = False; break
        if not ok: break
    F_exact = None
    if ok and all(r >= 0 for rs in rho for r in rs):
        F_exact = max(sum(abs(a - b) for a, b in zip(rho[i], rho[j])) / 2 for i, j in pairs)
    return F_lp, F_exact, S

GRAPHS = {
    'C4 chain': (4, [(0, 1), (1, 2), (2, 3)]),
    'C5 ring':  (5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
    'C6 ring':  (6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]),
}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--sample", type=int, default=None)
    ap.add_argument("--max-size", type=int, default=5); ap.add_argument("--no-lp", action="store_true")
    ap.add_argument("--lp-only", nargs=2, type=int, metavar=("START","STRIDE"), default=None,
                    help="LP-only mode for parallel runs: for C6, LP every STRIDE-th non-faithful core starting at START, one line per core; e.g. 20 processes with START=0..19 STRIDE=20")
    a = ap.parse_args()
    rng = np.random.default_rng(2026)
    if a.lp_only is not None:
        start, stride = a.lp_only
        n, edges = GRAPHS['C6 ring']; grp = stabilizer_group(n, edges)
        found = minimal_contradictions(n, grp, a.max_size)
        allcores = [(k, c) for k in sorted(found) for c in found[k]]
        bad = [(k, c) for (k, c) in allcores if not certify_core(n, grp, c)['exact']]
        print(f"C6: {len(bad)} cores where the construction is not faithful; this process handles indices {start}::{stride}")
        for idx in range(start, len(bad), stride):
            k, c = bad[idx]; t0 = time.time()
            F_lp, F_ex, Sz = faithful_lp_core(n, grp, c)
            print(f"idx {idx} size {k} core {' '.join(('+' if sg>0 else '-')+pauli_to_str(x,z,n) for (sg,x,z) in c)}: "
                  f"LP = {Fraction(F_lp).limit_denominator(400)} ({F_lp:.9f}), exact-upper = {F_ex}, strategies = {Sz}  [{time.time()-t0:.0f}s]")
        sys.exit(0)
    if a.sample:
        print("NOTE: --sample given; the count assertions are skipped (run without --sample to assert the quoted counts)")
    for name, (n, edges) in GRAPHS.items():
        t0 = time.time()
        grp = stabilizer_group(n, edges)
        found = minimal_contradictions(n, grp, a.max_size)
        if name == 'C6 ring' and not a.sample:
            assert {k: len(v) for k, v in found.items() if k <= 5} == {2: 0, 3: 0, 4: 360, 5: 1920}
        print(f"\n=== {name}: {len(grp)} stabilizers; minimal contradictions by size: "
              f"{ {k: len(v) for k, v in found.items()} }   [{time.time()-t0:.1f}s]")
        if name == 'C4 chain':
            core_names = {pauli_to_str(x, z, n) for (sg, x, z) in found[4][0]} if found[4] else set()
        for k, cores in found.items():
            if not cores: continue
            todo = cores
            if a.sample and len(cores) > a.sample:
                todo = [cores[i] for i in rng.choice(len(cores), a.sample, replace=False)]
            t1 = time.time(); res = [certify_core(n, grp, c) for c in todo]
            EXPECT = {('C4 chain', 4): (8, 8), ('C4 chain', 5): (8, 8), ('C5 ring', 4): (60, 60), ('C5 ring', 5): (180, 180),
                      ('C6 ring', 4): (360, 280), ('C6 ring', 5): (1920, 1616)}
            if (name, k) in EXPECT and not a.sample:
                assert (len(cores), sum(r['exact'] for r in res)) == EXPECT[(name, k)], (name, k, len(cores), sum(r['exact'] for r in res))
                assert all(r['feasible'] for r in res) and all(r['F'] == Fraction(1, k - 1) for r in res)
            from collections import Counter
            print(f"  size {k}: {len(todo)} of {len(cores)} cores checked; "
                  f"s = {set(str(r['s']) for r in res)}, Theorem-3 bound = {set(str(r['lower']) for r in res)}; "
                  f"construction feasible: {sum(r['feasible'] for r in res)}, faithful (all marginals exact): {sum(r['faithful'] for r in res)}, "
                  f"F values: {dict(Counter(str(r['F']) for r in res))}; "
                  f"EXACTLY CERTIFIED F_min = 1/(k-1): {sum(r['exact'] for r in res)} of {len(todo)}   [{time.time()-t1:.1f}s]")
            bad = [c for c, r in zip(todo, res) if not r['exact']]
            if bad and not a.no_lp:
                t2 = time.time(); lp_out = []
                for c in bad:
                    F_lp, F_ex, Sz = faithful_lp_core(n, grp, c); lp_out.append((F_lp, F_ex, Sz))
                from collections import Counter
                print(f"    faithful LP on the {len(bad)} cores where the construction is not faithful "
                      f"(collapsed strategy spaces {sorted(set(o[2] for o in lp_out))}): LP values "
                      f"{dict(Counter(str(Fraction(o[0]).limit_denominator(400)) for o in lp_out))}; "
                      f"exactly certified by rationalized LP solution: {sum(1 for o in lp_out if o[1] is not None)} "
                      f"with F = {set(str(o[1]) for o in lp_out if o[1] is not None)}   [{time.time()-t2:.1f}s]")
            ex = next((c for c in todo), None)
            print("    example core:", ', '.join(('+' if sg > 0 else '-') + pauli_to_str(x, z, n) for (sg, x, z) in ex))
