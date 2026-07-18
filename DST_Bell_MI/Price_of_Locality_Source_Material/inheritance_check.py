"""
Machine verification for Theorem 6 (coherence-support inheritance).

Part A: For stabilizer states (GHZ n=3,4,5; cluster chain C4; cluster ring C5),
        the set of Pauli operators with nonzero expectation is EXACTLY the
        stabilizer group. In particular, for GHZ every X/Y-containing element
        is full weight (all pairwise XX/XY/YY channels vanish) -- the Class 3
        channels are unpopulated.

Part B: Emission model. Post-emission superselected record
        |Psi> = (|0...0>|N> + |1...1>|N-n>)/sqrt(2).
        Sweep ALL channel assignments (eps_1..eps_n, q):
        eps_i in {0,+1,-1} (system winding transfer per site),
        q in {-(n+2)..(n+2)} (condensate winding transfer).
        Claim: expectation nonzero ONLY for
          (eps = 0, q = 0)          [diagonal / Class 1]
          (eps = +1^n, q = +n)      [full collective flip, imprint channel]
          (eps = -1^n, q = -n)      [conjugate]
        and the imprint coherence amplitude (normalized) = 1/2 exactly,
        independent of N.
"""
import numpy as np, itertools

I2 = np.eye(2); X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]]); Z = np.diag([1.,-1.]).astype(complex)
P = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

def kron_list(ops):
    M = np.array([[1.+0j]])
    for o in ops: M = np.kron(M, o)
    return M

def pauli_support(psi, n, tol=1e-9):
    sup = {}
    for s in itertools.product('IXYZ', repeat=n):
        M = kron_list([P[c] for c in s])
        v = np.vdot(psi, M @ psi)
        if abs(v) > tol:
            sup[''.join(s)] = round(v.real, 6) if abs(v.imag) < tol else v
    return sup

def ghz(n):
    v = np.zeros(2**n, dtype=complex); v[0] = 1; v[-1] = 1
    return v / np.sqrt(2)

def cluster(n, edges):
    v = np.ones(2**n, dtype=complex) / 2**(n/2)
    for (i,j) in edges:
        for b in range(2**n):
            if (b >> (n-1-i)) & 1 and (b >> (n-1-j)) & 1:
                v[b] *= -1
    return v

def stab_group_from_generators(gens, n):
    """gens: list of (pauli_string, sign). Multiply out the full group."""
    def mul(a, b):
        # multiply two pauli strings, return (string, phase)
        table = {('I','I'):('I',1),('I','X'):('X',1),('I','Y'):('Y',1),('I','Z'):('Z',1),
                 ('X','I'):('X',1),('X','X'):('I',1),('X','Y'):('Z',1j),('X','Z'):('Y',-1j),
                 ('Y','I'):('Y',1),('Y','X'):('Z',-1j),('Y','Y'):('I',1),('Y','Z'):('X',1j),
                 ('Z','I'):('Z',1),('Z','X'):('Y',1j),('Z','Y'):('X',-1j),('Z','Z'):('I',1)}
        s, ph = [], 1
        for ca, cb in zip(a, b):
            c, p = table[(ca, cb)]; s.append(c); ph *= p
        return ''.join(s), ph
    group = {('I'*n): 1}
    for bits in itertools.product([0,1], repeat=len(gens)):
        s, ph = 'I'*n, 1
        for bit, (g, sg) in zip(bits, gens):
            if bit:
                s2, p2 = mul(s, g); s, ph = s2, ph * p2 * sg
        group[s] = ph
    return group

print("="*72)
print("PART A: Pauli support = stabilizer group")
print("="*72)

# GHZ n = 3, 4, 5
for n in [3, 4, 5]:
    sup = pauli_support(ghz(n), n)
    xy_weights = sorted({sum(c in 'XY' for c in s) for s in sup})
    pairwise = [s for s in sup if 0 < sum(c in 'XY' for c in s) < n]
    print(f"\nGHZ n={n}: |support| = {len(sup)} (expect 2^{n} = {2**n})")
    print(f"  X/Y-weights present: {xy_weights}  (claim: only 0 and {n})")
    print(f"  partial-weight X/Y elements (Class 2/3 channels): {len(pairwise)} (claim: 0)")
    assert len(sup) == 2**n
    assert xy_weights == [0, n]
    assert len(pairwise) == 0

# C4 chain cluster: generators Z X Z I pattern
n = 4
edges = [(0,1),(1,2),(2,3)]
psi = cluster(n, edges)
sup = pauli_support(psi, n)
gens = [('XZII', 1), ('ZXZI', 1), ('IZXZ', 1), ('IIZX', 1)]
grp = stab_group_from_generators(gens, n)
match = set(sup.keys()) == set(grp.keys()) and all(abs(sup[k] - grp[k].real) < 1e-9 for k in grp)
wts = sorted({sum(c in 'XY' for c in s) for s in sup})
print(f"\nC4 chain cluster: |support| = {len(sup)}, equals stabilizer group: {match}")
print(f"  X/Y-weights present: {wts}  (mixed weights EXPECTED here -- support matches ITS OWN group)")
assert match

# C5 ring cluster
n = 5
edges = [(0,1),(1,2),(2,3),(3,4),(4,0)]
psi = cluster(n, edges)
sup = pauli_support(psi, n)
gens = [('XZIIZ',1), ('ZXZII',1), ('IZXZI',1), ('IIZXZ',1), ('ZIIZX',1)]
grp = stab_group_from_generators(gens, n)
match = set(sup.keys()) == set(grp.keys()) and all(abs(sup[k] - grp[k].real) < 1e-9 for k in grp)
print(f"\nC5 ring cluster: |support| = {len(sup)}, equals stabilizer group: {match}")
assert match

print("\n" + "="*72)
print("PART B: emission-record channel sweep")
print("="*72)

def emission_sweep(n, N):
    dC = N + n + 2               # condensate Fock space 0..N+n+1
    dS = 2**n
    psi = np.zeros(dS * dC, dtype=complex)
    psi[0*dC + N] = 1                 # |0...0>|N>
    psi[(dS-1)*dC + (N-n)] = 1        # |1...1>|N-n>
    psi /= np.sqrt(2)

    sp = np.array([[0,0],[1,0]], dtype=complex)   # sigma^+ : |0> -> |1>
    sm = sp.T.conj()
    proj = {0: I2, +1: sp, -1: sm}

    # condensate winding-shift (Susskind-Glogower style, unit matrix elements
    # so the reported amplitude is the NORMALIZED coherence)
    def shift(q):
        W = np.zeros((dC, dC), dtype=complex)
        for m in range(dC):
            if 0 <= m - q < dC:
                W[m - q, m] = 1     # removes q units of winding: |m> -> |m-q>
        return W

    hits = []
    for eps in itertools.product([0, 1, -1], repeat=n):
        Osys = kron_list([proj[e] for e in eps])
        for q in range(-(n+2), n+3):
            O = np.kron(Osys, shift(q))
            v = np.vdot(psi, O @ psi)
            if abs(v) > 1e-12:
                hits.append((eps, q, v.real if abs(v.imag) < 1e-12 else v))
    return hits

for n, N in [(3, 12), (4, 12)]:
    hits = emission_sweep(n, N)
    print(f"\nn={n}, N={N}: swept {3**n * (2*n+5)} channels; nonzero:")
    for eps, q, v in hits:
        print(f"   eps={eps}, condensate winding q={q:+d}: <O> = {v}")
    allowed = {((0,)*n, 0), ((1,)*n, n), ((-1,)*n, -n)}
    got = {(e, q) for e, q, _ in hits}
    assert got == allowed, f"UNEXPECTED CHANNELS: {got - allowed}"
    amp = [abs(v) for e, q, v in hits if e == (1,)*n][0]
    print(f"   imprint winding-{n} coherence amplitude = {amp}  (claim: 0.5 exactly)")
    assert abs(amp - 0.5) < 1e-12

# N-independence of the amplitude
amps = []
for N in [6, 10, 20, 40]:
    hits = emission_sweep(3, N)
    amps.append([abs(v) for e, q, v in hits if e == (1,1,1)][0])
print(f"\nN-independence (n=3, N=6,10,20,40): amplitudes = {amps}")
assert all(abs(a - 0.5) < 1e-12 for a in amps)

print("\nALL ASSERTIONS PASSED")
