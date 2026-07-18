"""
Delta = R: from fitted curiosity to spectral theorem (July 15, late session).

History: finite_beta.py's cold-side fits gave Delta = 2, 2, 4 at n = 3, 4, 5
(after the n=5 underflow was chased down). Discriminator question: is Delta
the staircase ratio R, or the spectral gap of the frustration count H on
satisfying sectors? Answer: both, because the gap IS R -- provably, all n.

THE LEMMA (Mermin-1990 identity, repurposed):
  W(class) := sum_s t_s sign_s(class) = Re[ prod_j (x_j + i y_j) ]
(one-line expansion: the real part of the product keeps exactly the
even-weight settings with sign i^|S| -> the Mermin target pattern).
Each factor = sqrt2 * e^(i k pi/4), k odd  ==>  product has modulus 2^(n/2),
phase quantized in pi/4, parity K = n mod 2.  Hence:
  odd n:  W in {+-R},        R = 2^((n-1)/2)   (two satisfying levels)
  even n: W in {0, +-2R},    R = 2^(n/2 - 1)   (three levels)
H = (n_c - W)/2  ==>  spectrum on an R-lattice, satisfying-sector gap = R,
ground = (n_c - W_max)/2, recovering s = 1/2 + W_max/2^n.

STATUS: spectral gap = R  [PROVEN via identity; identity machine-verified
exhaustively over all 4^n assignments, n = 3..9].
F(beta) - floor ~ e^(-R beta)  [machine-verified n = 3-5 in finite_beta.py;
the nonzero-coefficient step is generic but not separately proven].
Does NOT yet explain the n=3/n=4 identical-F(beta)-curves curiosity
(same gap, different level counts -- multiplicities must conspire; open).
"""
import numpy as np, itertools

for n in range(3, 10):
    settings = np.array([s for s in itertools.product([0,1], repeat=n) if sum(s)%2==0])
    t = np.array([1.0 if s.sum()%4==0 else -1.0 for s in settings])
    nc = len(t); R = 2**((n-1)//2)

    # exhaustive identity check over all 4^n deterministic assignments
    m = 4**n
    bits = ((np.arange(m)[:,None] >> np.arange(2*n)) & 1)
    x = 1 - 2*bits[:, :n]; y = 1 - 2*bits[:, n:]
    W = np.zeros(m)
    for s, tt in zip(settings, t):
        W += tt * np.where(s == 0, x, y).prod(axis=1)
    prod = (x + 1j*y).prod(axis=1)
    assert np.allclose(W, prod.real), f"n={n}: identity W = Re[prod] FAILS"

    Wvals = np.unique(np.round(W).astype(int))
    expected = np.array([-R, R]) if n % 2 else np.array([-2*R, 0, 2*R])
    assert (Wvals == expected).all(), f"n={n}: W values {Wvals} != {expected}"

    # spectrum + satisfying-sector gap on the reduced classes
    effs = np.array(list(itertools.product([0,1], repeat=n+1)), dtype=int)
    sign = np.where((effs[:,0][None,:] + settings @ effs[:,1:].T) % 2 == 1, -1., 1.)
    H = (sign != t[:, None]).sum(axis=0)
    gaps, grounds = [], []
    for si in range(nc):
        lv = np.unique(H[sign[si] == t[si]])
        gaps.append(lv[1] - lv[0]); grounds.append(lv[0])
    assert min(gaps) == max(gaps) == R, f"n={n}: gap != R"
    g = grounds[0]
    assert all(gr == g for gr in grounds)
    assert g == (nc - (R if n % 2 else 2*R)) // 2
    s_max = 1 - g / nc
    print(f"n={n}: identity OK on {m} assignments; W in {Wvals.tolist()}; "
          f"gap = {R} = R; ground = {g} (s = {s_max})")

print("\nALL ASSERTIONS PASSED -- Delta = R is spectral, exact, all computed n")
