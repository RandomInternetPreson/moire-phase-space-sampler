#!/usr/bin/env python3
"""
cross_ratio_v1.py — the two-setting cross-ratio X of the fourfold coincidence sum
from passive-basis BBM92 / two-channel Bell data.  Companion to arXiv:2608.18886
(Sec. 10, Design A) and the settings-torus brief v2.4.  numpy only.

    X = R4(HV,HV) * R4(DA,DA) / [ R4(HV,DA) * R4(DA,HV) ]

R4(a,b) is the sum of the four coincidence channels for one basis combination.
Any factorizable apparatus response A(a)B(b) — splitter ratios, path losses,
detector efficiencies — cancels in X.  For quantum mechanics with a factorizable
response X = 1.  The Sec. 5 detection profile gives X = 1.278; a mixture with
detection share s gives X(s) (function below).

The one non-factorizable term is outcome-dependent channel efficiency:
R4 = 1/4 [ S_A S_B + D_A D_B E ],  S = eta+ + eta-,  D = eta+ - eta-  (paper Eq. 6).
It is removed with delta = D/S from the channel singles (local states are maximally
mixed for a maximally entangled source) and E from the efficiency-reweighted table.

TWO INPUT MODES
---------------
Mode "table"  (four detectors per side; Kurtsiefer / Jennewein / Vienna / Aliro):
  JSON with, for each basis combination "HV_HV", "HV_DA", "DA_HV", "DA_DA",
  a 2x2 list of coincidence counts [[N++, N+-], [N-+, N--]] (A outcome rows,
  B outcome columns), and "singles": {"A_HV": [n+, n-], "A_DA": [...],
  "B_HV": [...], "B_DA": [...]}.  Optional: "duration_s", "window_s" for an
  accidental-coincidence subtraction.

Mode "peaks"  (time-multiplexed basis choice, two detectors per user; Bristol):
  JSON with "users": [names], "pairs": [{"a": name, "b": name,
  "central": 2x2, "side_ab": 2x2, "side_ba": 2x2}], where "central" holds the
  same-basis coincidences (HV_HV + DA_DA), "side_ab" the (A short, B long) peak,
  "side_ba" the (A long, B short) peak, each split by detector pair, and
  "singles": {name: [n+, n-]}.  With three or more users the per-user
  long/short efficiency ratios x_u are solved from the data themselves
  (side-peak ratios give x_b/x_a; central peaks overdetermine the x_u).
  With two users X is bounded, not measured, unless x_a, x_b are supplied.

USAGE
-----
  python cross_ratio_v1.py data.json            # analyse
  python cross_ratio_v1.py --selftest           # synthetic data, both modes
  python cross_ratio_v1.py --example table      # print an example input file
  python cross_ratio_v1.py --example peaks

Aaron Alai, November 2026.  ORCID 0000-0003-4704-7271.
"""
import json, sys
import numpy as np

# Exact-profile harmonics a_k of R4/mean = 1 + sum_k a_k cos(4k Delta), from
# solve_h_spectral.py (Sec. 5 profile).  Twelve harmonics; X(1) = 1.2779.
A_K = np.array([5.773047e-02, 9.834067e-03, 3.100629e-03, 1.331048e-03,
                6.853601e-04, 3.973622e-04, 2.503575e-04, 1.677139e-04,
                1.177627e-04, 8.582274e-05, 6.445952e-05, 4.963526e-05])

def fringe(delta_rad, share=1.0):
    k = np.arange(1, len(A_K) + 1)
    return 1.0 + share * np.sum(A_K * np.cos(4 * k * delta_rad))

def X_of_share(s):
    """Cross-ratio predicted for a mixture with detection share s."""
    return (fringe(0.0, s) / fringe(np.pi / 4, s)) ** 2

def share_of_X(X):
    """Invert X_of_share by bisection (X in [1, 1.278])."""
    if X <= 1.0:
        return 0.0
    lo, hi = 0.0, 1.0
    if X >= X_of_share(1.0):
        return 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if X_of_share(mid) < X:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

def leading_amplitude(X):
    """Equivalent single-harmonic amplitude a with (1+a)/(1-a) = sqrt(X)."""
    r = np.sqrt(X)
    return (r - 1) / (r + 1)

# --------------------------------------------------------------------------- core
def cross_ratio(N00, N45a, N45b, N4545):
    """X and the standard error of ln X from four Poisson counts."""
    counts = np.array([N00, N4545, N45a, N45b], float)
    if np.any(counts <= 0):
        raise ValueError("every fourfold sum must be positive")
    X = (N00 * N4545) / (N45a * N45b)
    sig_lnX = np.sqrt(np.sum(1.0 / counts))
    return X, sig_lnX

def report(label, X, sig_lnX):
    a = leading_amplitude(X)
    s = share_of_X(X)
    # error on share: propagate through dX/ds numerically
    ds = 1e-4
    dXds = (X_of_share(min(s + ds, 1)) - X_of_share(max(s - ds, 0))) / (2 * ds) if 0 < s < 1 else \
           (X_of_share(ds) - X_of_share(0)) / ds
    sig_s = X * sig_lnX / dXds
    print(f"  {label}")
    print(f"    X            = {X:.6f}   (ln X = {np.log(X):+.6f} +/- {sig_lnX:.6f})")
    print(f"    deviation    = {(X - 1) / (X * sig_lnX):+.2f} sigma from X = 1")
    print(f"    equiv. cos4  = {a:+.5f}  (peak-to-peak {2 * a:+.4f} of mean)")
    print(f"    Sec.5 share  = {s:.4f} +/- {sig_s:.4f}   (X = 1.278 at share 1)")
    print(f"    95% upper bound on share: {max(0.0, s) + 1.645 * sig_s:.4f}")

# --------------------------------------------------------------- imbalance calib.
def delta_from_singles(n_plus, n_minus):
    return (n_plus - n_minus) / (n_plus + n_minus)

def reweighted_E(table2x2, etaA, etaB):
    """Correlator from a 2x2 coincidence table with per-channel efficiencies."""
    t = np.asarray(table2x2, float)
    w = t / np.outer(etaA, etaB)
    sgn = np.array([[1, -1], [-1, 1]], float)
    return np.sum(sgn * w) / np.sum(w)

def coincidence_marginal_delta(table2x2):
    """delta_A, delta_B from coincidence marginals (exact at relative angle 45)."""
    t = np.asarray(table2x2, float)
    dA = (t[0].sum() - t[1].sum()) / t.sum()
    dB = (t[:, 0].sum() - t[:, 1].sum()) / t.sum()
    return dA, dB

def accidentals(singles_a, singles_b, duration_s, window_s):
    """Expected accidental coincidences per detector pair over the run."""
    ra = np.asarray(singles_a, float) / duration_s
    rb = np.asarray(singles_b, float) / duration_s
    return np.outer(ra, rb) * window_s * duration_s

# ----------------------------------------------------------------- mode: table
def analyse_table(d):
    combos = ["HV_HV", "HV_DA", "DA_HV", "DA_DA"]
    tabs = {c: np.asarray(d[c], float) for c in combos}
    sing = {k: np.asarray(v, float) for k, v in d["singles"].items()}
    if "duration_s" in d and "window_s" in d:
        for c in combos:
            ab, bb = c.split("_")
            acc = accidentals(sing["A_" + ab], sing["B_" + bb], d["duration_s"], d["window_s"])
            tabs[c] = tabs[c] - acc
        print("  accidentals subtracted (duration and window given)")
    R4 = {c: tabs[c].sum() for c in combos}
    print("\nRAW (no calibration)")
    X, s = cross_ratio(R4["HV_HV"], R4["HV_DA"], R4["DA_HV"], R4["DA_DA"])
    report("cross-ratio from the four fourfold sums", X, s)

    # imbalance calibration from singles
    delta = {k: delta_from_singles(*v) for k, v in sing.items()}
    print("\nIMBALANCE delta = (eta+ - eta-)/(eta+ + eta-) from channel singles")
    for k in ["A_HV", "A_DA", "B_HV", "B_DA"]:
        print(f"    {k}: {delta[k]:+.4f}")
    print("  cross-check from coincidence marginals (exact where E = 0, i.e. cross-basis):")
    for c in ["HV_DA", "DA_HV"]:
        dA, dB = coincidence_marginal_delta(tabs[c])
        ab, bb = c.split("_")
        print(f"    {c}: delta_A {dA:+.4f} (singles {delta['A_' + ab]:+.4f})   "
              f"delta_B {dB:+.4f} (singles {delta['B_' + bb]:+.4f})")
    R4c = {}
    for c in combos:
        ab, bb = c.split("_")
        etaA = sing["A_" + ab] / sing["A_" + ab].sum()
        etaB = sing["B_" + bb] / sing["B_" + bb].sum()
        E = reweighted_E(tabs[c], etaA, etaB)
        corr = 1.0 + delta["A_" + ab] * delta["B_" + bb] * E
        R4c[c] = R4[c] / corr
        print(f"    {c}: E = {E:+.4f}, correction factor {corr:.6f}")
    print("\nCALIBRATED (Eq. 6 term removed)")
    X, s = cross_ratio(R4c["HV_HV"], R4c["HV_DA"], R4c["DA_HV"], R4c["DA_DA"])
    report("cross-ratio from the corrected fourfold sums", X, s)

# ----------------------------------------------------------------- mode: peaks
def analyse_peaks(d):
    users = d["users"]
    pairs = d["pairs"]
    sing = {k: np.asarray(v, float) for k, v in d["singles"].items()}
    delta = {u: delta_from_singles(*sing[u]) for u in users}
    idx = {u: i for i, u in enumerate(users)}
    n = len(users)
    print(f"\n{n} users, {len(pairs)} pairs")
    print("  per-user imbalance from singles: " +
          ", ".join(f"{u} {delta[u]:+.4f}" for u in users))

    # per pair: C (same basis, both HV_HV and DA_DA), S_ab, S_ba; calibrate C for the
    # Eq. 6 term (E != 0 in the same-basis peak), sides have E ~ 0.
    rows = []
    for p in pairs:
        a, b = p["a"], p["b"]
        C = np.asarray(p["central"], float); Sab = np.asarray(p["side_ab"], float)
        Sba = np.asarray(p["side_ba"], float)
        etaA = sing[a] / sing[a].sum(); etaB = sing[b] / sing[b].sum()
        Ec = reweighted_E(C, etaA, etaB)
        Cc = C.sum() / (1.0 + delta[a] * delta[b] * Ec)
        rows.append((a, b, C.sum(), Cc, Sab.sum(), Sba.sum()))
        r = Cc / (Sab.sum() + Sba.sum())
        print(f"  pair {a}-{b}: C {C.sum():.0f} (E {Ec:+.3f}, calibrated {Cc:.0f}), "
              f"sides {Sab.sum():.0f} / {Sba.sum():.0f};  r = C/(S_ab+S_ba) = {r:.5f}, "
              f"q = S_ab/S_ba = {Sab.sum() / Sba.sum():.5f}")

    # Model: C_ab = K (1 + x_a x_b) rho ;  S_ab = K x_b ;  S_ba = K x_a
    # with rho = F(0)/F(45) the physical ratio, X = rho^2, x_u = long/short
    # efficiency ratio of user u.  Unknowns: ln rho and ln x_u.
    if "x" in d:                       # externally calibrated path ratios
        x = {u: float(d["x"][u]) for u in users}
        print("  using supplied long/short efficiency ratios x_u")
        lnrho, w = [], []
        for a, b, C, Cc, Sab, Sba in rows:
            r = Cc / (Sab + Sba)
            lnrho.append(np.log(r * (x[a] + x[b]) / (1 + x[a] * x[b])))
            w.append(1.0 / (1 / Cc + 1 / (Sab + Sba)))
        lnrho = np.array(lnrho); w = np.array(w)
        m = np.sum(w * lnrho) / np.sum(w); sig = 1 / np.sqrt(np.sum(w))
        report("X = rho^2 from calibrated pairs", np.exp(2 * m), 2 * sig)
        return
    if n < 3:
        a, b, C, Cc, Sab, Sba = rows[0]
        r = Cc / (Sab + Sba)
        sig = np.sqrt(1 / Cc + 1 / (Sab + Sba))
        print("\n  Two users: x_a, x_b not identifiable from the data alone.")
        print("  If both long paths are lossier than the short ones (x_a, x_b <= 1),")
        print("  then (1 + x_a x_b)/(x_a + x_b) >= 1 and r >= rho, so X <= r^2:")
        report("UPPER BOUND X <= r^2 (equality if either x = 1)", r * r, 2 * sig)
        print("  Supply {\"x\": {user: ratio}} from a classical path calibration for a measurement.")
        return
    # weighted nonlinear least squares in log space (Gauss-Newton)
    def residuals(theta):
        lnrho, lnx = theta[0], theta[1:]
        x = np.exp(lnx); res, wts = [], []
        for a, b, C, Cc, Sab, Sba in rows:
            xa, xb = x[idx[a]], x[idx[b]]
            res.append(np.log(Cc) - np.log(Sab + Sba) - lnrho - np.log((1 + xa * xb) / (xa + xb)))
            wts.append(1.0 / np.sqrt(1 / Cc + 1 / (Sab + Sba)))
            res.append(np.log(Sab / Sba) - (lnx[idx[b]] - lnx[idx[a]]))
            wts.append(1.0 / np.sqrt(1 / Sab + 1 / Sba))
        return np.array(res) * np.array(wts)
    # initial point: ln x ratios from the side peaks (least squares), then a scan
    # over the common scale of x against the central-peak residuals.  The model's
    # gradient in x vanishes at x = 1, so a start at x = 1 cannot be used.
    Aq = np.zeros((len(rows) + 1, n)); bq = np.zeros(len(rows) + 1)
    for i, (a, b, C, Cc, Sab, Sba) in enumerate(rows):
        Aq[i, idx[b]] = 1; Aq[i, idx[a]] = -1; bq[i] = np.log(Sab / Sba)
    Aq[-1, :] = 1.0                                     # fixes the mean of ln x to 0
    lnx0, *_ = np.linalg.lstsq(Aq, bq, rcond=None)
    best = None
    for scale in np.linspace(np.log(0.3), np.log(1.5), 241):
        for lnrho in (0.0,):
            th = np.concatenate([[lnrho], lnx0 + scale])
            # profile out ln rho analytically: it enters the r-residuals linearly
            x = np.exp(th[1:]); num = den = 0.0
            for a, b, C, Cc, Sab, Sba in rows:
                wgt = 1.0 / (1 / Cc + 1 / (Sab + Sba))
                num += wgt * (np.log(Cc) - np.log(Sab + Sba) - np.log((1 + x[idx[a]] * x[idx[b]]) / (x[idx[a]] + x[idx[b]])))
                den += wgt
            th[0] = num / den
            val = float(residuals(th) @ residuals(th))
            if best is None or val < best[0]:
                best = (val, th.copy())
    theta = best[1]
    lam = 1e-3                                           # Levenberg-Marquardt
    for _ in range(200):
        r0 = residuals(theta)
        J = np.zeros((len(r0), n + 1)); h = 1e-7
        for j in range(n + 1):
            t = theta.copy(); t[j] += h
            J[:, j] = (residuals(t) - r0) / h
        H = J.T @ J; g = J.T @ r0
        step = np.linalg.solve(H + lam * np.diag(np.diag(H) + 1e-12), -g)
        step = np.clip(step, -0.2, 0.2)
        trial = theta + step
        trial[1:] = np.clip(trial[1:], np.log(0.1), np.log(10.0))   # physical range of x_u
        step = trial - theta
        r1 = residuals(trial)
        if r1 @ r1 < r0 @ r0:
            theta += step; lam = max(lam / 3, 1e-9)
        else:
            lam *= 10
        if np.max(np.abs(step)) < 1e-10 or lam > 1e6:
            break
    r0 = residuals(theta)
    cov = np.linalg.inv(J.T @ J)
    chi2 = float(r0 @ r0); dof = len(r0) - (n + 1)
    print(f"\n  self-calibration fit: chi2 = {chi2:.2f} for {dof} dof")
    if np.any(np.abs(theta[1:]) > np.log(9.0)):
        print("  WARNING: a path ratio sits at the edge of its allowed range (0.1..10); the common")
        print("  scale of the x_u is not determined by these data. More pairs, more counts, or an")
        print("  external calibration of one user's long/short transmission is needed.")
    for u in users:
        print(f"    x_{u} = {np.exp(theta[1 + idx[u]]):.5f} +/- {np.exp(theta[1 + idx[u]]) * np.sqrt(cov[1 + idx[u], 1 + idx[u]]):.5f}")
    report("X = rho^2 from the network self-calibration", np.exp(2 * theta[0]), 2 * np.sqrt(cov[0, 0]))

# -------------------------------------------------------------------- synthetic
def synth_table(N0, share, eta, w_basis, V=0.98, rng=None):
    """eta: {'A_HV':[e+,e-], ...}; w_basis: {'A':[wHV,wDA], 'B':[...]}."""
    rng = rng or np.random.default_rng(1)
    rel = {"HV_HV": 0.0, "DA_DA": 0.0, "HV_DA": np.pi / 4, "DA_HV": -np.pi / 4}
    out = {}
    sgn = np.array([[1, -1], [-1, 1]], float)
    for c, dl in rel.items():
        ab, bb = c.split("_")
        E = -V * np.cos(2 * dl)                          # singlet correlator
        eA = np.asarray(eta["A_" + ab]); eB = np.asarray(eta["B_" + bb])
        wa = w_basis["A"][0 if ab == "HV" else 1]; wb = w_basis["B"][0 if bb == "HV" else 1]
        mean = N0 * wa * wb * 0.25 * np.outer(eA, eB) * (1 + sgn * E) * fringe(dl, share)
        out[c] = rng.poisson(mean).tolist()
    Ns = 40 * N0
    out["singles"] = {k: rng.poisson(Ns * w_basis[k[0]][0 if k[2:] == "HV" else 1] * 0.5 * np.asarray(eta[k])).tolist()
                      for k in eta}
    return out

def synth_peaks(N0, share, users, eta, xlong, V=0.98, rng=None):
    rng = rng or np.random.default_rng(2)
    sgn = np.array([[1, -1], [-1, 1]], float)
    pairs = []
    for i, a in enumerate(users):
        for b in users[i + 1:]:
            eA = np.asarray(eta[a]); eB = np.asarray(eta[b])
            base = N0 * 0.25 * np.outer(eA, eB)
            C = base * (1 + xlong[a] * xlong[b]) * (1 - sgn * V) * fringe(0.0, share)
            Sab = base * xlong[b] * fringe(np.pi / 4, share)          # A short, B long: E = 0
            Sba = base * xlong[a] * fringe(np.pi / 4, share)
            pairs.append({"a": a, "b": b, "central": rng.poisson(C).tolist(),
                          "side_ab": rng.poisson(Sab).tolist(), "side_ba": rng.poisson(Sba).tolist()})
    singles = {u: rng.poisson(40 * N0 * 0.5 * (1 + xlong[u]) * np.asarray(eta[u])).tolist() for u in users}
    return {"users": users, "pairs": pairs, "singles": singles}

def selftest():
    np.set_printoptions(precision=5)
    print("=" * 78); print("SELF-TEST, mode table: 20% constant imbalance per side, unequal basis")
    print("splitters, V = 0.98, 2.5e7 coincidences per setting")
    eta = {"A_HV": [0.6, 0.4], "A_DA": [0.6, 0.4], "B_HV": [0.4, 0.6], "B_DA": [0.4, 0.6]}
    w = {"A": [0.55, 0.45], "B": [0.48, 0.52]}
    for share in (0.0, 0.01, 0.1, 1.0):
        print("-" * 78); print(f"planted detection share {share}  ->  expected X = {X_of_share(share):.6f}")
        d = synth_table(1.0e8, share, eta, w)
        analyse_table(d)
    print("=" * 78); print("SELF-TEST, mode peaks: 4 users, long-path ratios 0.75..0.9, 10% imbalance")
    users = ["A", "B", "C", "D"]
    eta = {"A": [0.55, 0.45], "B": [0.5, 0.5], "C": [0.45, 0.55], "D": [0.52, 0.48]}
    xl = {"A": 0.75, "B": 0.9, "C": 0.8, "D": 0.85}
    for share in (0.0, 0.05, 1.0):
        print("-" * 78); print(f"planted detection share {share}  ->  expected X = {X_of_share(share):.6f}; "
                               f"planted x = {xl}")
        d = synth_peaks(2.0e7, share, users, eta, xl)
        analyse_peaks(d)
    print("-" * 78); print("two users only (bound, not measurement), share 0:")
    d = synth_peaks(2.0e7, 0.0, ["A", "B"], {"A": eta["A"], "B": eta["B"]}, {"A": 0.75, "B": 0.9})
    analyse_peaks(d)

EXAMPLES = {
"table": {"HV_HV": [[120000, 4000], [4100, 119000]], "HV_DA": [[62000, 61000], [61500, 62500]],
          "DA_HV": [[61800, 62200], [61900, 62100]], "DA_DA": [[121000, 3900], [4000, 120000]],
          "singles": {"A_HV": [2500000, 2400000], "A_DA": [2300000, 2350000],
                      "B_HV": [2600000, 2550000], "B_DA": [2450000, 2500000]},
          "duration_s": 3600, "window_s": 1.0e-9},
"peaks": {"users": ["A", "B", "C"], "pairs": [{"a": "A", "b": "B", "central": [[8978, 886258], [816811, 8341]], "side_ab": [[233719, 233641], [216006, 216163]], "side_ba": [[207596, 208147], [191958, 192162]]}, {"a": "A", "b": "C", "central": [[7624, 851794], [696521, 8034]], "side_ab": [[170791, 192932], [158080, 177407]], "side_ba": [[195530, 220292], [179609, 203287]]}, {"a": "B", "b": "C", "central": [[7750, 855538], [758888, 8623]], "side_ab": [[164387, 186796], [164703, 186068]], "side_ba": [[210909, 238596], [211050, 238572]]}], "singles": {"A": [74892009, 69121635], "B": [75997799, 75999213], "C": [63915848, 72053407]}}}

if __name__ == "__main__":
    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--help"):
        print(__doc__); sys.exit(0)
    if args[0] == "--selftest":
        selftest(); sys.exit(0)
    if args[0] == "--example":
        print(json.dumps(EXAMPLES[args[1]], indent=1)); sys.exit(0)
    with open(args[0]) as f:
        d = json.load(f)
    if "pairs" in d:
        analyse_peaks(d)
    else:
        analyse_table(d)
