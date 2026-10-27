#!/usr/bin/env python3
"""
census_cells_v1.py -- is the census the SAME offset density in every (tau, scale) cell, up to a weight?

Why.  Every birth model so far (birth_model_v1/v2) assumes rho(xi_A, xi_B, cell) = W[cell] * f(xi_A, xi_B): one shape over
the offsets, reweighted per cell.  Warm-started fits to the K = 6 census leave a 70% residual at J = 2, 3, 4, so something
the form cannot carry is large.  This script measures, on the LP censuses themselves, how far they are from that
cell-separable form -- and what a cell-separable census costs in correlator miss.

For each census file given (LP rho json, or a census_bandlimit_v4 json with a K entry via K6_K):
  (1) the 15 x NV matrix R[cell, offset]; its singular values; the share of ||R||^2 carried by rank 1 (separable),
      rank 2, rank 3 -- once as is, once with every row normalized to unit mass (shape only, weights removed)
  (2) per-cell geometry: mass, |F(2,2)| and its phase (the sum lock), preferred S (mod 90), |F(2,0)|, |F(0,2)|,
      the aligned/crossed split of D (aligned: D within 45 deg of 18; crossed: within 45 deg of 128)
  (3) shape similarity between cells: cosine similarity of the unit-mass rows, averaged within the same tau, within the
      same scale, and across everything
  (4) if the banks are next to the file: the correlator miss and price of the best rank-1 (cell-separable) approximation
      of the census -- i.e. what a birth that stamps the same offsets at every dial setting can reach at best, given
      this census's shape
Usage:  python census_cells_v1.py rho_nsh32_random0_eps0.01_nang161.json census_bandlimit_v4_r0_e01_k6_v4_rho_nsh32_random0_eps0.01_nang161.json
Env:    K6_K (default 6) for census_bandlimit jsons.   numpy only; seconds.   Writes census_cells_v1_out.json next to the first file.
"""
import os, sys, json
import numpy as np

TAUS = [0.5, 1.0, 2.0, 4.0, 8.0]; SCALES = [1.0, 1.8, 2.6]; K6_K = int(os.environ.get("K6_K", 6))
_deg = np.pi / 180; G = 256; DL = np.pi / G; FLOOR = (np.sqrt(2.0) - 1.0) / 3.0; S_QM = 2 * np.sqrt(2)

def basis(offs, K):
    xa, xb = offs[:, 0], offs[:, 1]; cols = [np.ones(len(offs))]
    for ka in range(-K, K + 1):
        for kb in range(-K, K + 1):
            if ka < 0 or (ka == 0 and kb <= 0): continue
            arg = 2 * (ka * xa + kb * xb); cols.append(np.cos(arg)); cols.append(np.sin(arg))
    return np.column_stack(cols)

def load_census(path):
    d = json.load(open(path))
    if "results" in d:                        # band-limited census: rebuild from coefficients; offsets/meta from its rho file
        dr = json.load(open(os.path.join(os.path.dirname(path) or ".", d["file"]))); offs = np.asarray(dr["offsets"], float); NV = len(offs)
        r = next((r for r in d["results"] if r.get("K") == K6_K and "coeff" in r), None)
        if r is None: sys.exit(f"{path}: no feasible K={K6_K} entry with coefficients")
        c = np.asarray(r["coeff"], float); Phi = basis(offs, K6_K); nflat = c.shape[0]
        rho = np.clip(np.concatenate([Phi @ c[ts] for ts in range(nflat)]), 0, None); name = f"K={K6_K} census ({os.path.basename(path)})"; d = dr
    else:
        offs = np.asarray(d["offsets"], float); NV = len(offs)
        rho = np.asarray(d["rho"], float); name = f"LP census ({os.path.basename(path)})"
    rho = rho / rho.sum(); nflat = rho.size // NV
    return d, name, offs, rho.reshape(nflat, NV)

def load_banks(path, d):
    f = os.path.join(os.path.dirname(path) or ".", f"banks_nsh{d['nsh']}_{d['mode']}{d['seed']}_nang{d['nang']}_nw{d['nw']}.npz")
    if not os.path.exists(f): return None
    z = np.load(f); nflat = int(z["nflat"])
    Qs = [np.vstack([z[f"q{k}_{p}"] for k in range(nflat)]) for p in range(4)]
    Svec = np.concatenate([z[f"s{k}"] for k in range(nflat)]); Ecoef = np.hstack([z[f"e{k}"] for k in range(nflat)])
    angs = np.linspace(2, 178, int(d["nang"])); return Qs, Svec, Ecoef, -np.cos(2 * angs * _deg)

def misses_price(banks, rho):
    Qs, Svec, Ecoef, tgt = banks
    E = Ecoef @ rho; qm = [Qs[p].T @ rho for p in range(4)]
    md = max(0.5 * np.sum(np.abs(qm[a] - qm[b])) * DL for a in range(4) for b in range(a + 1, 4))
    return float(np.abs(E - tgt).max()), float(abs(Svec @ rho + S_QM)), float(md / FLOOR)

def pixel(w, xa, xb, ka, kb):
    F = np.sum(w * np.exp(2j * (ka * xa + kb * xb))) / w.sum(); return abs(F), np.degrees(np.angle(F))

def analyse(path):
    d, name, offs, R = load_census(path); nflat, NV = R.shape; xa, xb = offs[:, 0], offs[:, 1]
    S = (xa + xb) % np.pi; D = (xa - xb) % np.pi
    print(f"\n===== {name}: {nflat} cells x {NV} offsets =====")
    out = dict(file=os.path.basename(path), name=name)
    # (1) separability
    for label, M in (("as is (weights included)", R), ("rows normalized (shape only)", R / np.maximum(R.sum(axis=1, keepdims=True), 1e-300))):
        s = np.linalg.svd(M, compute_uv=False); tot = np.sum(s ** 2)
        shares = np.cumsum(s ** 2) / tot
        print(f"(1) rank shares of ||R||^2, {label}: rank1 {shares[0]:.3f}  rank2 {shares[1]:.3f}  rank3 {shares[2]:.3f}  rank5 {shares[4]:.3f}   (1.000 = cell-separable)")
        out[f"rank_shares_{'raw' if 'as is' in label else 'shape'}"] = shares[:8].tolist()
    # (2) per-cell geometry
    print("(2) per cell:  mass   |F22|  phi22   S_pref(mod90)  |F20|  |F02|   aligned  crossed")
    rows = []
    aligned = (np.abs(((D - 18 * _deg + np.pi / 2) % np.pi) - np.pi / 2) < 45 * _deg)
    crossed = (np.abs(((D - 128 * _deg + np.pi / 2) % np.pi) - np.pi / 2) < 45 * _deg)
    for c in range(nflat):
        w = R[c]; mass = w.sum()
        if mass < 1e-9:
            print(f"    cell {c:2d} (tau={TAUS[c//3]:<3} s={SCALES[c%3]:<3}): {mass:.3f}   (empty)"); rows.append(dict(cell=c, mass=0.0)); continue
        a22, p22 = pixel(w, xa, xb, 2, 2); a20, _ = pixel(w, xa, xb, 2, 0); a02, _ = pixel(w, xa, xb, 0, 2)
        z = np.sum(w * np.exp(4j * S)) / mass; spref = (np.degrees(np.angle(z)) / 4) % 90
        al = w[aligned].sum() / mass; cr = w[crossed].sum() / mass
        print(f"    cell {c:2d} (tau={TAUS[c//3]:<3} s={SCALES[c%3]:<3}): {mass:.3f}   {a22:.3f}  {p22:+7.1f}   {spref:6.1f}        {a20:.3f}  {a02:.3f}    {al:.2f}     {cr:.2f}")
        rows.append(dict(cell=c, tau=TAUS[c // 3], scale=SCALES[c % 3], mass=float(mass), F22=float(a22), phi22=float(p22), S_pref=float(spref), F20=float(a20), F02=float(a02), aligned=float(al), crossed=float(cr)))
    out["cells"] = rows
    # (3) shape similarity
    U = R / np.maximum(np.linalg.norm(R, axis=1, keepdims=True), 1e-300); C = U @ U.T
    live = R.sum(axis=1) > 1e-9
    def avg(mask_pairs):
        v = [C[i, j] for i, j in mask_pairs if live[i] and live[j]]; return (np.mean(v) if v else float("nan")), len(v)
    pairs = [(i, j) for i in range(nflat) for j in range(i + 1, nflat)]
    same_tau = [(i, j) for i, j in pairs if i // 3 == j // 3]; same_sc = [(i, j) for i, j in pairs if i % 3 == j % 3]
    a_all, n_all = avg(pairs); a_tau, n_tau = avg(same_tau); a_sc, n_sc = avg(same_sc)
    print(f"(3) cosine similarity of cell shapes: all pairs {a_all:.3f} (n={n_all})   same tau {a_tau:.3f} (n={n_tau})   same scale {a_sc:.3f} (n={n_sc})   (1 = identical shape; random censuses of this sparsity ~ 0.05)")
    out["similarity"] = dict(all=float(a_all), same_tau=float(a_tau), same_scale=float(a_sc))
    # (4) the cost of separability
    banks = load_banks(path, d)
    if banks is not None:
        mE, mS, pr = misses_price(banks, R.reshape(-1)); print(f"(4) this census: band miss {mE:.5f}  CHSH miss {mS:.2e}  price {pr:.4f}")
        u, s, vt = np.linalg.svd(R, full_matrices=False)
        for r in (1, 2, 3):
            A = (u[:, :r] * s[:r]) @ vt[:r]; A = np.clip(A, 0, None); A /= A.sum()
            mE, mS, pr = misses_price(banks, A.reshape(-1)); out[f"rank{r}_approx"] = dict(band_miss=mE, chsh_miss=mS, price=pr)
            print(f"    best rank-{r} approximation{' (cell-separable: one shape, 15 weights)' if r == 1 else ''}: band miss {mE:.5f}  CHSH miss {mS:.5f}  price {pr:.4f}")
    else:
        print("(4) banks not found next to the file; skipping the miss/price of the separable approximation")
    return out

if __name__ == "__main__":
    res = [analyse(p) for p in sys.argv[1:]]
    json.dump(res, open(os.path.join(os.path.dirname(sys.argv[1]) or ".", "census_cells_v1_out.json"), "w"), indent=1)
