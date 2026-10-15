"""
exact_certificates_fwht_v2.py -- float-free certification of the Mermin
staircase F_min(n) = R/(2(R+1)), for n well beyond 13.

Same certificate as exact_certificates.py (lower bound = Theorem 3 overlap
bound with exact satisfiability counts; upper bound = ground-state
construction with exact overlap counts), but the violation count H is
obtained by a Walsh-Hadamard transform instead of an N x 2^n satisfaction
matrix, so memory is O(2^(n+1)) integers instead of O(2^(2n-1)) bytes.

  H(P, gamma) = 2^(n-2) - (1/2) (-1)^P T(gamma),      T = WHT[t-hat]
  where t-hat(S) = E(s) = +1 if |S| = 0 mod 4, -1 if |S| = 2 mod 4, and 0
  for odd |S|  (Eq. (6) of the paper, Appendix E).

Everything is integer arithmetic: the WHT of a {-1,0,+1} vector of length
2^n is exact in int64 (|T| <= 2^(n-1)). Fractions for the final ratios.

Checks performed at each n (all are hard assertions):
  1. m = N - min H equals (R+1)N/(2R)                 [Mermin's local bound]
  2. odd n : (N - m)/(N - 1) == R/(2(R+1))            [Theorem 3, full set]
     even n: same on the embedded Mermin_{n-1} (X on the last party)
  3. every sector A_s contains the same number g > 0 of ground states
     (feasibility of the ground-state construction)
  4. max over S_n-orbits of setting pairs of 1 - |G_s & G_s'|/g == R/(2(R+1))
     (value of the construction; pairs enumerated per orbit (k, k', |s & s'|),
      which is exact by Theorem 2; brute force over all pairs at n <= 8 as a
      check of the orbit reduction)
  5. the closed-form counts of Lemma 11 of the paper, asserted (not just
     printed): |ground| = 2^n (odd n) / 2^(n-1) (even n); g = R(R+1); the
     minimum pair overlap is R(R+2)/2; every non-complementary pair of
     settings has overlap exactly R(R+2)/2 (equidistance); at even n each
     complementary pair (S, [n] minus S) has overlap g (identical densities).
     (v2: these assertions added; v1 printed the counts but asserted ratios only.)

Usage:  python exact_certificates_fwht_v2.py 3 4 5 ... 20
        python exact_certificates_fwht_v2.py --upto 22
Memory: ~ 2^(n+1) * 8 bytes for H plus a few 2^n arrays; n = 24 needs ~0.5 GB.
Time:   dominated by O(n^3) orbit pairs x O(2^n) per pair.
"""
import sys, time, argparse
from fractions import Fraction
import numpy as np


def popcount64(x):
    x = x.astype(np.int64)
    c = np.zeros_like(x)
    while True:
        nz = x != 0
        if not nz.any():
            break
        c += (x & 1)
        x = x >> 1
    return c


def fwht_int(a):
    """In-place integer Walsh-Hadamard transform (unnormalized). Exact."""
    a = a.copy()
    n = a.shape[0]
    h = 1
    while h < n:
        a = a.reshape(-1, 2 * h)
        x = a[:, :h].copy()
        y = a[:, h:].copy()
        a[:, :h] = x + y
        a[:, h:] = x - y
        a = a.reshape(-1)
        h *= 2
    return a


def certify(n, brute_pairs=False, quiet=False):
    t0 = time.time()
    R = 2 ** ((n - 1) // 2)
    staircase = Fraction(R, 2 * (R + 1))
    N = 1 << (n - 1)

    masks = np.arange(1 << n, dtype=np.int64)
    pc = popcount64(masks)
    # t-hat on all masks: 0 for odd |S|; +1 / -1 for |S| = 0 / 2 mod 4
    that = np.where(pc % 2 == 1, 0, np.where(pc % 4 == 0, 1, -1)).astype(np.int64)
    T = fwht_int(that)                      # T(gamma) = sum_S t-hat(S) (-1)^{|gamma & S|}
    # H(P=0,gamma) = N/2 - T/2 ; H(P=1,gamma) = N/2 + T/2
    assert ((N - T) % 2 == 0).all()
    H0 = (N - T) // 2
    H1 = (N + T) // 2
    Hmin = int(min(H0.min(), H1.min()))
    m = N - Hmin
    assert Fraction(m, N) == Fraction(R + 1, 2 * R), (n, Fraction(m, N))

    # ---- lower bound (Theorem 3) ------------------------------------------
    if n % 2 == 1:
        lower = Fraction(N - m, N - 1)
        lower_desc = f"full set, m={m}"
    else:
        # embedded Mermin_{n-1}: settings with X on the last party (bit n-1 clear)
        that_sub = that.copy()
        that_sub[(masks >> (n - 1)) & 1 == 1] = 0
        Tsub = fwht_int(that_sub)
        Nsub = N // 2
        Hsub_min = int(min(((Nsub - Tsub) // 2).min(), ((Nsub + Tsub) // 2).min()))
        m_sub = Nsub - Hsub_min
        lower = Fraction(Nsub - m_sub, Nsub - 1)
        lower_desc = f"embedded Mermin_{n-1}, m'={m_sub} of N'={Nsub}"
    assert lower == staircase, (n, lower, staircase)

    # ---- ground set ---------------------------------------------------------
    g0 = np.nonzero(H0 == Hmin)[0]          # gammas with P = 0 in ground set
    g1 = np.nonzero(H1 == Hmin)[0]          # gammas with P = 1 in ground set
    n_ground = len(g0) + len(g1)
    assert n_ground == (1 << n if n % 2 == 1 else 1 << (n - 1)), (n, n_ground)   # Lemma 10

    sett = masks[pc % 2 == 0]
    t_bit = ((pc[pc % 2 == 0] % 4) != 0).astype(np.int64)   # GF(2) target
    assert len(sett) == N

    def ground_sat(si):
        """Boolean vector over the ground set: does ground class j satisfy setting si?"""
        s = int(sett[si]); tb = int(t_bit[si])
        p0 = popcount64(g0 & s) % 2         # P = 0: parity = popcount(gamma & S)
        p1 = (popcount64(g1 & s) + 1) % 2   # P = 1
        return np.concatenate([p0 == tb, p1 == tb])

    # sector ground counts: by covariance only |S| matters, but check every k and
    # a random sample of settings as well
    ks = sorted(set(pc[pc % 2 == 0].tolist()))
    rng = np.random.default_rng(12345)
    sample = sorted(set([int(np.nonzero(pc[pc % 2 == 0] == k)[0][0]) for k in ks]
                        + rng.integers(0, N, size=min(N, 64)).tolist()))
    g_counts = {si: int(ground_sat(si).sum()) for si in sample}
    g = g_counts[sample[0]]
    assert g > 0 and all(v == g for v in g_counts.values()), "sector ground counts differ"
    assert g == R * (R + 1), (n, g)                                              # Lemma 11

    # ---- upper bound: construction overlaps -------------------------------
    sett_list = sett.tolist()
    index = {s: i for i, s in enumerate(sett_list)}
    inter_min = None
    noncomp = set(); comp = set()          # overlaps seen for non-complementary / complementary pairs
    if brute_pairs:
        cache = [ground_sat(i) for i in range(N)]
        full = (1 << n) - 1
        for i in range(N):
            for j in range(i + 1, N):
                v = int((cache[i] & cache[j]).sum())
                (comp if sett_list[i] ^ sett_list[j] == full else noncomp).add(v)
                inter_min = v if inter_min is None else min(inter_min, v)
    else:
        for k in ks:
            sat1 = ground_sat(index[(1 << k) - 1])
            for k2 in ks:
                if k2 < k:
                    continue
                lo, hi = max(0, k + k2 - n), min(k, k2)
                for ov in range(lo, hi + 1):
                    if k == k2 and ov == k:
                        continue
                    s2 = ((1 << ov) - 1) | (((1 << (k2 - ov)) - 1) << k)
                    v = int((sat1 & ground_sat(index[s2])).sum())
                    (comp if (ov == 0 and k + k2 == n) else noncomp).add(v)
                    inter_min = v if inter_min is None else min(inter_min, v)
    upper = 1 - Fraction(inter_min, g)
    assert upper == staircase, (n, upper, staircase)
    assert inter_min == R * (R + 2) // 2, (n, inter_min)                        # Lemma 11
    assert noncomp == {R * (R + 2) // 2}, (n, noncomp)                           # equidistance
    assert comp == (set() if n % 2 == 1 else {g}), (n, comp)                      # complementary pairs

    if not quiet:
        print(f"n={n:>2}: N={N:>7}, Hmin={Hmin:>6}, m={m:>7} (s={Fraction(m, N)}), "
              f"|ground|={n_ground:>8}, g={g:>7}, min-overlap={inter_min:>7}  "
              f"lower={lower} [{lower_desc}]  upper={upper}  == R/(2(R+1)) = {staircase}"
              f"   [{time.time() - t0:.1f}s]")
    return dict(n=n, N=N, Hmin=Hmin, m=m, ground=n_ground, g=g, inter_min=inter_min,
                lower=lower, upper=upper, staircase=staircase)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("ns", nargs="*", type=int)
    ap.add_argument("--upto", type=int, default=None)
    ap.add_argument("--brute-upto", type=int, default=8,
                    help="brute-force all setting pairs for n <= this (checks the orbit reduction)")
    a = ap.parse_args()
    ns = a.ns or list(range(3, (a.upto or 13) + 1))
    for n in ns:
        certify(n, brute_pairs=(n <= a.brute_upto))
    print("\nALL VALUES CERTIFIED EXACTLY, Lemma-11 counts asserted at every n (integer WHT + Fractions; no floating point)")
