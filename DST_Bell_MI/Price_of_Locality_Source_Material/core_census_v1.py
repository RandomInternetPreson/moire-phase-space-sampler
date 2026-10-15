"""
core_census_v1.py -- COMPLETE census of the minimal contradictions (cores) of a
stabilizer state, at all sizes, by GF(2) linear algebra.

A set of signed stabilizers is a contradiction iff the (party, Pauli)
occurrence vectors sum to zero over GF(2) and the signs multiply to -1; a core
is a minimal one.  Any proper dependent subset of a contradiction, or its
complement, is a smaller contradiction, so cores are exactly the circuits
(minimal dependent sets) of the occurrence vectors with sign product -1.  The
circuits are found by enumerating the kernel of the occurrence map, which is
feasible when its dimension d = N - rank is small (C4: d = 8, C5: d = 22;
C6: d = 52, not feasible this way -- use cluster_cores_v4.py --max-size for C6).

Usage: python core_census_v1.py "C4 chain" | "C5 ring"
Prints cores by size and asserts the counts quoted in the paper.
"""
import itertools, sys, importlib.util, time
import numpy as np
import os
spec = importlib.util.spec_from_file_location('cc', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cluster_cores_v4.py'))
cc = importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)
def setup(name):
    n, edges = cc.GRAPHS[name]; grp = cc.stabilizer_group(n, edges)
    elems = [(sg, x, z) for (x, z), sg in grp.items() if (x, z) != (0, 0)]
    vec = [cc.occurrence_vector(x, z, n) for (sg, x, z) in elems]
    sgn = [sg for (sg, x, z) in elems]
    return n, elems, vec, sgn
def gf2_rank(vecs):
    rows = list(vecs); r = 0
    piv = {}
    for v in rows:
        for p, pv in piv.items():
            if (v >> p) & 1: v ^= pv
        if v:
            p = v.bit_length() - 1
            # reduce other pivots
            for q in list(piv):
                if (piv[q] >> p) & 1: piv[q] ^= v
            piv[p] = v; r += 1
    return r
def kernel_basis(vec):
    # kernel of the map GF(2)^N -> GF(2)^{3n}, e_j -> vec[j]
    N = len(vec); 
    # Gaussian elimination tracking combinations
    rows = [(vec[j], 1 << j) for j in range(N)]
    piv = {}
    basis = []
    for v, comb in rows:
        for p, (pv, pc) in piv.items():
            if (v >> p) & 1: v ^= pv; comb ^= pc
        if v:
            p = v.bit_length() - 1
            for q in list(piv):
                if (piv[q][0] >> p) & 1: piv[q] = (piv[q][0] ^ v, piv[q][1] ^ comb)
            piv[p] = (v, comb)
        else:
            basis.append(comb)
    return basis
name = sys.argv[1]
n, elems, vec, sgn = setup(name)
N = len(vec); print(name, "N =", N, "rank =", gf2_rank(vec))
kb = kernel_basis(vec); d = len(kb); print("kernel dim", d)
if d > 24: sys.exit("kernel too large to enumerate; use cluster_cores_v4.py --max-size")
# enumerate all kernel vectors (dependent sets), gray code
deps = []
cur = 0
for i in range(1, 1 << d):
    # gray code: flip bit of lowest set bit of i
    b = (i & -i).bit_length() - 1
    cur ^= kb[b]
    deps.append(cur)
deps.sort(key=lambda m: bin(m).count('1'))
# minimal dependent sets (circuits): no proper subset dependent
circuits = []
for m in deps:
    if any((c & m) == c for c in circuits): continue
    circuits.append(m)
from collections import Counter
sizes = Counter()
for m in circuits:
    s = 1
    for j in range(N):
        if (m >> j) & 1: s *= sgn[j]
    if s == -1: sizes[bin(m).count('1')] += 1
print("minimal contradictions by size:", dict(sorted(sizes.items())))
print("circuits by size:", dict(sorted(Counter(bin(m).count('1') for m in circuits).items())))

EXPECTED = {"C4 chain": {4: 8, 5: 8, 6: 8, 7: 16}, "C5 ring": {4: 60, 5: 180, 6: 572, 7: 1980, 8: 5220, 9: 9640, 10: 8912}}
if name in EXPECTED:
    assert dict(sorted(sizes.items())) == EXPECTED[name], (dict(sizes), EXPECTED[name])
    print("census matches the counts quoted in the paper; total cores =", sum(sizes.values()))
