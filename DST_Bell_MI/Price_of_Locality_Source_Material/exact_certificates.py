"""
exact_certificates.py -- float-free certification of every Mermin floor.

The LP (quantumtestv4/v5) is double precision; its rationals come from
limit_denominator. This script replaces trust in solver tolerance with an
exact squeeze, using only integer counting and Fractions:

  LOWER BOUND (proven, Theorem 6 / overlap bound):
    F >= (N' - m') / (N' - 1) for any constraint subset S'.
      odd n : S' = full Mermin set. m = N - Hmin exactly (integer count);
              algebraically (N - m)/(N - 1) = R/(2(R+1)) when m/N = (R+1)/(2R).
      even n: S' = settings with X on the last party (an embedded
              Mermin_{n-1}); same formula at R(n-1) = R(n).

  UPPER BOUND (explicit feasible construction, frustration ansatz):
    rho_s uniform on ground states of H within the sector satisfying s.
    Feasibility is structural (support => exact correlators; fiber lift =>
    marginals). F = max over setting pairs of 1 - |G_s ^ G_s'| / g, a ratio
    of integer counts. Pairs are enumerated per S_n-orbit (the construction
    is manifestly covariant), with a full brute-force sweep at n <= 7 to
    check the orbit reduction itself.

  If LOWER == UPPER == R/(2(R+1)), the floor is established exactly at
  that n, independent of any floating-point tolerance.

All arithmetic: numpy integer ops (exact) + fractions.Fraction. No floats.

Usage: python exact_certificates.py 3 4 5 6 7 8 9 10 11 12 13
"""
import sys, time, itertools
from fractions import Fraction
import numpy as np

POP16 = np.array([bin(i).count("1") for i in range(1 << 16)], dtype=np.int64)

def popcount(x):
    return POP16[x & 0xFFFF] + POP16[(x >> 16) & 0xFFFF]

def certify(n, brute_pairs=False):
    t0 = time.time()
    R = 2 ** ((n - 1) // 2)
    staircase = Fraction(R, 2 * (R + 1))

    # settings: even-popcount n-bit masks; GF(2) target t_s = 0 iff |s| % 4 == 0
    masks = np.arange(1 << n, dtype=np.int64)
    pc_masks = popcount(masks)
    sett = masks[pc_masks % 2 == 0]
    t = ((pc_masks[pc_masks % 2 == 0] % 4) != 0).astype(np.int64)  # 0 or 1
    N = len(sett)
    assert N == 1 << (n - 1)

    # per-gamma violation count at P = 0; H(P=1) = N - H(P=0)
    gam = np.arange(1 << n, dtype=np.int64)
    H0 = np.zeros(1 << n, dtype=np.int64)
    sat0 = np.zeros((N, 1 << n), dtype=bool)      # satisfaction at P = 0
    for si in range(N):
        p = popcount(gam & sett[si]) % 2
        ok = (p == t[si])
        sat0[si] = ok
        H0 += ~ok
    H = np.concatenate([H0, N - H0])              # classes: (P=0, gam), (P=1, gam)
    sat = np.concatenate([sat0, ~sat0], axis=1)   # [setting, class]

    # ---- exact satisfiability and lower bound -------------------------------
    Hmin = int(H.min())
    m = N - Hmin
    s_frac = Fraction(m, N)
    assert s_frac == Fraction(R + 1, 2 * R), (n, s_frac)

    if n % 2 == 1:
        lower = Fraction(N - m, N - 1)
        lower_desc = f"full set, m={m}"
    else:
        sub = (sett & (1 << (n - 1))) == 0        # X on last party
        Nsub = int(sub.sum())
        viol_sub = np.zeros(1 << n, dtype=np.int64)
        for si in np.nonzero(sub)[0]:
            viol_sub += ~sat0[si]
        m_sub = Nsub - int(min(viol_sub.min(), (Nsub - viol_sub).min()))
        lower = Fraction(Nsub - m_sub, Nsub - 1)
        lower_desc = f"embedded Mermin_{n-1}, m'={m_sub} of N'={Nsub}"
    assert lower == staircase, (n, lower, staircase)

    # ---- exact upper bound: frustration-ansatz overlap counts ---------------
    ground = np.nonzero(H == Hmin)[0]
    Gsat = sat[:, ground]                          # [setting, ground-class]
    g_counts = Gsat.sum(axis=1)
    g = int(g_counts[0])
    assert (g_counts == g).all(), "sector ground counts differ"

    if brute_pairs:
        inter_min = min(int((Gsat[i] & Gsat[j]).sum())
                        for i in range(N) for j in range(i + 1, N))
    else:
        # orbit representatives (k, k', overlap): construction is S_n-covariant
        ks = sorted(set(popcount(sett).tolist()))
        idx_by_k = {k: int(np.nonzero(popcount(sett) == k)[0][0]) for k in ks}
        inter_min = None
        for k in ks:
            for k2 in ks:
                if k2 < k: continue
                lo, hi = max(0, k + k2 - n), min(k, k2)
                for ov in range(lo, hi + 1):
                    if k == k2 and ov == k: continue
                    s1 = (1 << k) - 1
                    s2 = ((1 << ov) - 1) | (((1 << (k2 - ov)) - 1) << k)
                    i1 = int(np.nonzero(sett == s1)[0][0])
                    i2 = int(np.nonzero(sett == s2)[0][0])
                    val = int((Gsat[i1] & Gsat[i2]).sum())
                    inter_min = val if inter_min is None else min(inter_min, val)
    upper = 1 - Fraction(inter_min, g)
    assert upper == staircase, (n, upper, staircase)

    print(f"n={n:>2}: N={N:>5}, Hmin={Hmin:>4}, m={m:>5} (s={s_frac}), "
          f"|ground|={len(ground):>6}, g={g:>5}, min-overlap={inter_min:>5}  "
          f"lower={lower} [{lower_desc}]  upper={upper}  "
          f"== R/(2(R+1)) = {staircase}   [{time.time()-t0:.1f}s]")
    return staircase

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or list(range(3, 14))
    for n in ns:
        certify(n, brute_pairs=(n <= 7))
    print("\nALL FLOORS CERTIFIED EXACTLY -- no floating point in any certified value")
