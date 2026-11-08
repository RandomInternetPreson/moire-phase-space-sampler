#!/usr/bin/env python3
"""
census_fuzz_v3.py -- the disc's edge: a birth distribution of the pair's phase lambda_0 that is NOT uniform.

The paper starts every pair with lambda_0 uniform on [0, pi), independent of the offsets and the dial, and the census rho lives on
(xi_A, xi_B, tau, s) alone.  Aaron (1 Nov): the disc may have a thickness, a knife-edge or a blunt rim, a multimodal rim -- in the
model, that is the birth fuzz p(lambda_0), and the richest version is a JOINT birth over (lambda_0, xi_A, xi_B, tau, s).
Test: with the paper's step reader, let the joint birth be  w(atom, lambda_0) = sum_l w_l(atom) f_l(lambda_0),  f = {1, cos 2 l lambda_0,
sin 2 l lambda_0, l = 1..L}  (L = 0 is the paper's census), w_l(atom) a trigonometric polynomial of degree K over the offsets (the usual
band limit), positivity of the joint density at a lambda_0 grid, mass, CHSH exact, band eps, MD rows on the lambda(tau) marginal (the
reader reads lambda(tau) only, so the paper's marginal MD is the right measure).  One LP, exact optimum.  Two scopes:
  cell : the fuzz depends on the dial only (w_l for l >= 1 constant over the offsets within each (tau, s) cell)
  atom : the fuzz depends on the approach phases too (full joint; w_l band-limited in xi like w_0)
Question: does a shaped rim let a SMOOTHER census (lower K) reach the floor?  (Control: L = 0 must give 1.700 / 1.429 / 1.258 at K = 6/8/10.)

Needs the PARTICLE store from pipeline_fuzz_v1.py (final lambda-bin of every initial particle, every setting, every atom, every tau).
Two stages, both cached:  contract (once, with --Lmax harmonics) -> fuzz_contracted_Lmax<Lmax>_<store>.npz   then the LP for one (K, L, scope).
Usage:  python census_fuzz_v3.py part_nsh32_grid0_nang161_nw2048.npz --L 2 --K 6 --scope atom --eps 0.01
        (K = 0 means the full census: identity basis over the offsets)
Writes census_fuzz_v3_L<L>_K<K>_<scope>.json; prints a SUMMARY block.
v3 (5 Nov) makes the atom-scope LPs tractable (v1/v2's L >= 1, K >= 8 interior-point solves ran 40-100 h and one returned "numerical
difficulties"): (a) --Krim: the rim harmonics (l >= 1) get their own, lower band limit over the offsets (the L = 1 rim's skeleton is
low-order: (1,0), (1,2), (3,2) -- degree 3), while the census w_0 keeps K; (b) --Mpl: positivity points per harmonic order, default 16
instead of 64 (the LP is then a slight relaxation; the clipped census is re-measured and gives the matching upper bound, so the
price is bracketed, not assumed); (c) HiGHS tolerances 1e-7 and an iteration-free dual-simplex retry on failure.  Log: fuzz3_*.txt.
v2 (2 Nov) adds what v1 could not show after the fact: (a) the solution is saved (census_fuzz_v3_L<L>_K<K>_<scope>_solution.npz:
u = joint-density harmonics per atom, c = basis coefficients) so the rim can be read later; (b) a CLIPPED re-verification -- the joint
density is clipped at zero on a fine lambda_0 grid, renormalized, and its MD and correlator miss re-measured: a valid census whose
price is an UPPER bound, bracketing the LP's slight lower bound (the dips between positivity points); (c) per-atom rim statistics:
the mass-weighted modulation depth of each harmonic (the pooled rim of v1 cancels by the +90 symmetry and says nothing for the atom
scope), and the circular correlation of the rim's preferred direction with the approach-phase sum, difference and each station.
"""
import os, sys, json, time, argparse
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
import census_bandlimit_v4 as cb                      # basis, run_lp, STATUS

_deg = np.pi / 180; S_QM = 2 * np.sqrt(2); FLOOR = (np.sqrt(2.0) - 1.0) / 3.0
settings = [(0.0, 22.5), (0.0, 67.5), (45.0, 22.5), (45.0, 67.5)]; signs = [+1, -1, +1, +1]

def fbasis(lam0, L):
    cols = [np.ones_like(lam0)]
    for l in range(1, L + 1): cols += [np.cos(2 * l * lam0), np.sin(2 * l * lam0)]
    return np.column_stack(cols)                                                  # (n, 2L+1)

def contract(store, L, out):
    """E[ai, atom, l], Ecs[p, atom, l], Svec[atom, l], Q[p][atom, g, l]  for the step reader; cached to `out`."""
    t0 = time.time(); z = np.load(store)
    Bpair, Bang = z["Bpair"], z["Bang"]                     # (nflat, 4, NV, NW), (nflat, nang, NV, NW) uint8 final bins
    offs = z["offsets"]; angs = z["angs"]; NW = int(z["nw"]); G = int(z["G"])
    nflat, nang, NV = Bpair.shape[0], Bang.shape[1], Bpair.shape[2]; NT = nflat * NV; DL = np.pi / G
    lam_c = (np.arange(G) + 0.5) * np.pi / G; lam0 = np.arange(NW) * np.pi / NW
    F = fbasis(lam0, L); nb0 = F.shape[1]
    SGA0 = np.sign(np.cos(2 * lam_c))
    print(f"contracting {os.path.basename(store)} with L={L}: {nflat} cells x {NV} offsets x ({nang}+4 settings) x {NW} particles -> {nb0} lambda_0 harmonics", flush=True)
    E = np.zeros((nang, NT, nb0))
    for ai, tdeg in enumerate(angs):
        kern = SGA0 * np.sign(np.cos(2 * (lam_c + np.pi / 2 - tdeg * _deg)))          # (G,)
        E[ai] = kern[Bang[:, ai].reshape(NT, NW)] @ F / NW                            # gather then (NT, NW) @ (NW, nb0)
    del Bang
    Ecs = np.zeros((4, NT, nb0)); Q = []
    idx = (np.arange(NT)[:, None] * G)
    for p, (a, b) in enumerate(settings):
        kern = np.sign(np.cos(2 * (lam_c - a * _deg))) * np.sign(np.cos(2 * (lam_c + np.pi / 2 - b * _deg)))
        Bp = Bpair[:, p].reshape(NT, NW)
        Ecs[p] = kern[Bp] @ F / NW
        Qp = np.zeros((NT, G, nb0))
        flat = (idx + Bp).ravel()
        for l in range(nb0):
            Qp[:, :, l] = np.bincount(flat, weights=np.broadcast_to(F[:, l], (NT, NW)).ravel(), minlength=NT * G).reshape(NT, G) / (NW * DL)
        Q.append(Qp)
    Svec = np.einsum("p,pal->al", np.array(signs, float), Ecs)
    np.savez(out, E=E, Ecs=Ecs, Svec=Svec, Q=np.stack(Q), offsets=offs, angs=angs, nflat=nflat, NV=NV, G=G, NW=NW, L=L)
    print(f"  contracted in {time.time()-t0:.0f}s -> {out}", flush=True)

def build_and_solve(C, L, K, scope, eps, mode, method, tlim, M_per_L=16, Krim=None, tol=1e-7):
    """one LP.  mode 'price': min MD;  mode 'miss': min sigma with the band widened by sigma (certificate by value).
       L <= the contraction's Lmax: the first 2L+1 harmonic columns are used."""
    nb0 = 2 * L + 1; assert nb0 <= 2 * int(C["L"]) + 1, "L exceeds the contraction's Lmax"
    E, Ecs, Svec, Q = C["E"][:, :, :nb0], C["Ecs"][:, :, :nb0], C["Svec"][:, :nb0], C["Q"][:, :, :, :nb0]; offs = C["offsets"]; angs = C["angs"]
    nflat, NV, G, NW = int(C["nflat"]), int(C["NV"]), int(C["G"]), int(C["NW"]); NT = nflat * NV
    nang = len(angs); DL = np.pi / G; RED = G
    Phi0 = np.eye(NV) if K == 0 else cb.basis(offs, K)[0]
    Kr = K if Krim is None else Krim
    Phir = Phi0 if Kr == K else (np.eye(NV) if Kr == 0 else cb.basis(offs, Kr)[0])
    Phil = Phir if scope == "atom" else Phi0[:, :1]                 # cell scope: the fuzz harmonics are constant over the offsets in a cell
    bases = [Phi0] + [Phil] * (2 * L); nbs = [B.shape[1] for B in bases]; ncs = [nflat * nb for nb in nbs]; nc = sum(ncs)
    off_c = np.cumsum([0] + ncs)                                     # column offsets of c_l
    n_u = NT * nb0; off_u = nc; cmb = [(i, j) for i in range(4) for j in range(i + 1, 4)]; n_s = len(cmb) * RED
    off_t = off_u + n_u; off_s = off_t + 1; nv = off_s + n_s        # [c | u | t(or sigma) | slacks]
    def red(Xl, l):      # reduce an (NT, ...) per-atom array on harmonic l onto its basis: returns (ncs[l], ...)
        B = bases[l]; return np.concatenate([np.tensordot(B.T, Xl[ts * NV:(ts + 1) * NV], axes=(1, 0)) for ts in range(nflat)], axis=0)
    # --- equalities: u_l = blockdiag(B_l) c_l   (positivity is written on u)
    rows_eq, b_eq = [], []
    for l in range(nb0):
        Bl = sp.block_diag([sp.csr_matrix(bases[l])] * nflat).tocsr()      # (NT, ncs[l])
        A = sp.hstack([sp.csr_matrix((NT, off_c[l])), Bl, sp.csr_matrix((NT, nc - off_c[l + 1])), sp.csr_matrix((NT, l * NT)), -sp.identity(NT, format="csr"),
                       sp.csr_matrix((NT, n_u - (l + 1) * NT + 1 + n_s))]).tocsr()
        rows_eq.append(A); b_eq += [0.0] * NT
    # mass (f_l, l >= 1, integrate to zero over the uniform particle grid): sum_atoms u_0 = 1
    r = sp.lil_matrix((1, nv)); r[0, off_u:off_u + NT] = 1.0; rows_eq.append(r.tocsr()); b_eq.append(1.0)
    # CHSH exact, in c
    Sr = np.concatenate([red(Svec[:, l], l) for l in range(nb0)])
    r = sp.lil_matrix((1, nv)); r[0, :nc] = Sr; rows_eq.append(r.tocsr()); b_eq.append(-S_QM)
    Aeq = sp.vstack(rows_eq).tocsr(); beq = np.array(b_eq)
    # --- inequalities
    blocks, bub = [], []
    Er = np.hstack([red(E[:, :, l].T, l).T for l in range(nb0)])    # (nang, nc)
    tgt = -np.cos(2 * angs * _deg)
    if mode == "price":
        # MD rows: |sum_l Qr_l[p] - Qr_l[p']| per bin <= slack; 0.5*sum slack <= t
        Qr = [np.vstack([red(Q[p][:, :, l], l) for l in range(nb0)]) * DL for p in range(4)]      # (nc, G)
        for ci, (i, j) in enumerate(cmb):
            Dif = (Qr[i] - Qr[j]).T                                                               # (G, nc)
            S_ = sp.lil_matrix((RED, nv));
            for g in range(RED): S_[g, off_s + ci * RED + g] = -1.0
            S_ = S_.tocsr()
            blocks += [sp.hstack([sp.csr_matrix(Dif), sp.csr_matrix((RED, nv - nc))]).tocsr() + S_, sp.hstack([sp.csr_matrix(-Dif), sp.csr_matrix((RED, nv - nc))]).tocsr() + S_]
            bub += [0.0] * (2 * RED)
        T = sp.lil_matrix((len(cmb), nv))
        for ci in range(len(cmb)): T[ci, off_s + ci * RED:off_s + (ci + 1) * RED] = 0.5; T[ci, off_t] = -1.0
        blocks.append(T.tocsr()); bub += [0.0] * len(cmb)
        Ep = sp.hstack([sp.csr_matrix(Er), sp.csr_matrix((nang, nv - nc))]).tocsr()
        blocks += [Ep, -Ep]; bub += list(tgt + eps) + list(-tgt + eps)
    else:
        Ep = sp.hstack([sp.csr_matrix(Er), sp.csr_matrix((nang, n_u)), sp.csr_matrix(-np.ones((nang, 1))), sp.csr_matrix((nang, n_s))]).tocsr()
        En = sp.hstack([sp.csr_matrix(-Er), sp.csr_matrix((nang, n_u)), sp.csr_matrix(-np.ones((nang, 1))), sp.csr_matrix((nang, n_s))]).tocsr()
        blocks += [Ep, En]; bub += list(tgt + eps) + list(-tgt + eps)
    # positivity of the joint density at M lambda_0 grid points:  -sum_l F[j, l] u_l(atom) <= 0
    M = max(8, M_per_L * L) if L > 0 else 1
    Fg = fbasis((np.arange(M) + 0.37) * np.pi / M, L)                                             # (M, nb0), off-grid points
    Pu = sp.hstack([sp.kron(sp.csr_matrix(-Fg[:, l:l + 1]), sp.identity(NT, format="csr")) for l in range(nb0)]).tocsr()   # (M*NT, n_u)
    blocks.append(sp.hstack([sp.csr_matrix((M * NT, nc)), Pu, sp.csr_matrix((M * NT, 1 + n_s))]).tocsr()); bub += [0.0] * (M * NT)
    Aub = sp.vstack(blocks).tocsr(); bub = np.array(bub)
    cobj = np.zeros(nv); cobj[off_t] = 1.0
    bounds = [(None, None)] * (nc + n_u) + [(0, None)] * (1 + n_s)
    print(f"    LP: {nv} variables ({nc} coefficients, {n_u} joint-density values), {Aub.shape[0]} inequality rows, {Aeq.shape[0]} equalities, {Aub.nnz + Aeq.nnz:,} nonzeros", flush=True)
    opts = {"primal_feasibility_tolerance": tol, "dual_feasibility_tolerance": tol, "ipm_optimality_tolerance": tol}
    if tlim > 0: opts["time_limit"] = tlim
    t = time.time(); res = linprog(cobj, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds, method=method, options=opts)
    print(f"    [{method}] status {res.status} = {cb.STATUS.get(res.status, '?')}; {res.message.strip()}   [{time.time()-t:.0f}s]", flush=True)
    if res.status not in (0, 2) and method != "highs-ds":
        print("    retrying with dual simplex", flush=True)
        t = time.time(); res = linprog(cobj, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds, method="highs-ds", options=opts)
        print(f"    [highs-ds] status {res.status} = {cb.STATUS.get(res.status, '?')}; {res.message.strip()}   [{time.time()-t:.0f}s]", flush=True)
    return res, dict(nc=nc, off_u=off_u, off_t=off_t, NT=NT, nb0=nb0, Fg=Fg, Er=Er, tgt=tgt, Q=Q, bases=bases, off_c=off_c, nflat=nflat, NV=NV, G=G, DL=DL, L=L, offs=offs, E_atom=E, M=M, Krim=Kr)

def verify(res, meta, eps):
    """everything re-measured from the solution: joint density minimum on a fine lambda_0 grid, mass, band, MD of the lambda(tau) marginal."""
    x = res.x; NT, nb0 = meta["NT"], meta["nb0"]; u = x[meta["off_u"]:meta["off_u"] + NT * nb0].reshape(nb0, NT)   # u[l, atom]
    Ff = fbasis((np.arange(512) + 0.5) * np.pi / 512, meta["L"])                                    # fine grid
    joint = Ff @ u                                                                                     # (512, NT)
    v = dict(joint_min=float(joint.min()), mass=float(u[0].sum()), joint_min_rel=float(joint.min() / max(u[0].mean(), 1e-300)))   # relative to the mean atom density
    E = meta["Er"] @ x[:meta["nc"]]; v["band_miss"] = float(np.abs(E - meta["tgt"]).max() - eps)
    qm = [np.einsum("agl,la->g", meta["Q"][p], u) for p in range(4)]                                   # lambda(tau) marginal per setting
    v["md"] = float(max(0.5 * np.sum(np.abs(qm[a] - qm[b])) * meta["DL"] for a in range(4) for b in range(a + 1, 4))); v["price"] = v["md"] / FLOOR
    v["ok"] = v["joint_min_rel"] > -1e-3 and abs(v["mass"] - 1) < 1e-7 and v["md"] >= FLOOR - 1e-7 and v["band_miss"] <= 1e-6
    print(f"    verify: joint density min {v['joint_min']:+.2e} = {100*v['joint_min_rel']:+.3f}% of the mean atom density (512-point lambda_0 grid; positivity was imposed at {meta['M']} points)  "
          f"mass {v['mass']:.9f}  band miss over eps {v['band_miss']:+.2e}  MD/floor {v['price']:.4f}   {'OK' if v['ok'] else '** VIOLATED **'}", flush=True)
    # the rim: per cell, the fuzz's harmonics relative to its mass
    rim = []
    for ts in range(meta["nflat"]):
        sl = slice(ts * meta["NV"], (ts + 1) * meta["NV"]); m0 = u[0, sl].sum()
        h = []
        for l in range(1, meta["L"] + 1):
            a_, b_ = u[2 * l - 1, sl].sum(), u[2 * l, sl].sum()
            h.append((float(np.hypot(a_, b_) / max(m0, 1e-12)), float(np.degrees(np.arctan2(b_, a_)) / (2 * l))))
        rim.append(dict(mass=float(m0), harmonics=h))
    v["rim"] = rim
    # --- per-atom rim statistics (the atom scope): depth and direction of the lambda_0 modulation, atom by atom
    if meta["L"] > 0:
        offs = meta["offs"]; xa = np.tile(offs[:, 0], meta["nflat"]); xb = np.tile(offs[:, 1], meta["nflat"]); w0 = np.clip(u[0], 0, None)
        depth = {}; corr = {}
        for l in range(1, meta["L"] + 1):
            a_, b_ = u[2 * l - 1], u[2 * l]; amp = np.hypot(a_, b_)
            depth[f"h{l}"] = float(amp.sum() / max(w0.sum(), 1e-300))                 # mass-weighted modulation depth (0 uniform, 1 full)
            psi2 = np.arctan2(b_, a_)                                                  # 2 l psi: the rim's preferred direction in lambda_0 (x l)
            for name, X in (("sum", xa + xb), ("diff", xa - xb), ("xiA", xa), ("xiB", xb)):
                for mult in (1, 2):                                                    # does 2l*psi track mult*2*X ?
                    z = np.sum(amp * np.exp(1j * (psi2 - mult * 2 * X))) / max(amp.sum(), 1e-300)
                    corr[f"h{l}~{mult}x{name}"] = (float(abs(z)), float(np.degrees(np.angle(z))))
        v["depth"] = depth; v["corr"] = corr
        print("    rim per atom (mass-weighted modulation depth; 0 = uniform birth): " + "  ".join(f"{k} {d:.3f}" for k, d in depth.items()))
        best = sorted(corr.items(), key=lambda kv: -kv[1][0])[:4]
        print("    rim direction tracks (circular correlation |R|, lag):  " + "  ".join(f"{k}: {r:.2f} @ {ph:+.0f}" for k, (r, ph) in best) + "   (|R| = 1: the rim's lean is a fixed function of that phase)")
    # --- clipped census: an UPPER bound to bracket the LP's slight lower bound
    jc = np.clip(joint, 0, None); jc *= 1.0 / (jc.mean(axis=0).sum())                 # renormalize so that the mass is 1 (mean over lambda_0 x sum over atoms)
    uc = (Ff.T @ jc) / Ff.shape[0] * np.array([1] + [2] * (2 * meta["L"]))[:, None]  # back to harmonics (orthogonality on the fine grid)
    qm = [np.einsum("agl,la->g", meta["Q"][p], uc) for p in range(4)]
    v["clip_md_over_floor"] = float(max(0.5 * np.sum(np.abs(qm[a] - qm[b])) * meta["DL"] for a in range(4) for b in range(a + 1, 4)) / FLOOR)
    Ec_u = np.einsum("atl,lt->a", meta["E_atom"], uc) if "E_atom" in meta else None
    if Ec_u is not None:
        v["clip_band_miss"] = float(np.abs(Ec_u - meta["tgt"]).max())
        print(f"    clipped census (joint density clipped at 0, renormalized): MD/floor {v['clip_md_over_floor']:.4f} (upper bound), max |E - target| {v['clip_band_miss']:.5f} (eps {eps})")
    return u, v

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("store"); ap.add_argument("--L", type=int, default=1); ap.add_argument("--K", type=int, default=6)
    ap.add_argument("--scope", default="atom"); ap.add_argument("--eps", type=float, default=0.01); ap.add_argument("--method", default="highs-ipm"); ap.add_argument("--tlim", type=float, default=0)
    ap.add_argument("--Lmax", type=int, default=3, help="harmonics in the cached contraction (one contraction serves every L <= Lmax)")
    ap.add_argument("--save", type=int, default=1, help="save the solution (u, c) to an npz")
    ap.add_argument("--Krim", type=int, default=-1, help="band limit of the rim harmonics over the offsets (-1 = same as K)")
    ap.add_argument("--Mpl", type=int, default=16, help="positivity points per lambda_0 harmonic order")
    ap.add_argument("--tol", type=float, default=1e-7, help="HiGHS feasibility / optimality tolerances")
    ap.add_argument("--log", default="auto", help="also write everything printed to this file ('auto' = fuzz3_L<L>_K<K>_<scope>.txt, '' = none)")
    A = ap.parse_args(); t0 = time.time(); Lmax = max(A.Lmax, A.L)
    rtag = "" if A.Krim < 0 else f"_r{A.Krim}"
    logf = f"fuzz3_L{A.L}_K{A.K}{rtag}_{A.scope}.txt" if A.log == "auto" else A.log
    if logf:
        class Tee:
            def __init__(self, path): self.f = open(path, "a"); self.o = sys.__stdout__
            def write(self, x): self.o.write(x); self.f.write(x)
            def flush(self): self.o.flush(); self.f.flush()
        sys.stdout = Tee(logf)
    cf = f"fuzz_contracted_Lmax{Lmax}_{os.path.splitext(os.path.basename(A.store))[0]}.npz"
    if not os.path.exists(cf): contract(A.store, Lmax, cf)
    C = np.load(cf)
    print(f"census_fuzz_v3  L={A.L} (lambda_0 harmonics: {'uniform birth, the paper' if A.L == 0 else 'cos/sin 2..' + str(2*A.L) + ' lambda_0'})  K={A.K}{' (full census)' if A.K == 0 else ''}  scope={A.scope}  eps={A.eps}")
    print(f"question: does a shaped rim (birth fuzz in lambda_0, {'dial-dependent only' if A.scope == 'cell' else 'depending on the approach phases too'}) let the census of bandwidth K reach the floor more cheaply than L=0?")
    print("  price LP:", flush=True)
    Krim = None if A.Krim < 0 else A.Krim
    print(f"rim band limit Krim = {A.K if Krim is None else Krim}; positivity at {max(8, A.Mpl*A.L) if A.L else 1} lambda_0 points; tolerances {A.tol:g}")
    res, meta = build_and_solve(C, A.L, A.K, A.scope, A.eps, "price", A.method, A.tlim, M_per_L=A.Mpl, Krim=Krim, tol=A.tol)
    rec = dict(L=A.L, K=A.K, scope=A.scope, eps=A.eps)
    if res.status == 0:
        u, v = verify(res, meta, A.eps)
        if A.save:
            sp_ = f"census_fuzz_v3_L{A.L}_K{A.K}{rtag}_{A.scope}_solution.npz"; np.savez_compressed(sp_, u=u, c=res.x[:meta["nc"]], offsets=meta["offs"], L=A.L, K=A.K, scope=A.scope); print(f"  solution saved: {sp_}")
        st = "feasible" if v["ok"] else (f"feasible (joint density dips {100*v['joint_min_rel']:+.2f}% between positivity points: the price is a slight lower bound)" if v["joint_min_rel"] > -0.05 and abs(v["mass"] - 1) < 1e-7 and v["band_miss"] <= 1e-6 else "feasible-unverified")
        rec.update(status=st, price=v["price"], verify={k: v[k] for k in ("joint_min", "joint_min_rel", "mass", "band_miss", "md", "price", "ok")}, rim=v["rim"], depth=v.get("depth"), corr=v.get("corr"), clip_md_over_floor=v.get("clip_md_over_floor"), clip_band_miss=v.get("clip_band_miss"))
        print(f"  MD/floor = {rec['price']:.4f}   {rec['status']}   [{time.time()-t0:.0f}s]")
        if A.L > 0:
            print("  the rim per cell (amplitude of each lambda_0 harmonic relative to the cell's mass, and its phase in deg):")
            for ts, r in enumerate(v["rim"]):
                print(f"    cell {ts:2d}  mass {r['mass']:.4f}  " + "  ".join(f"h{l+1}: {amp:.3f} @ {ph:+6.1f}" for l, (amp, ph) in enumerate(r["harmonics"])))
            pooled = u.sum(axis=1); print("  pooled rim: " + "  ".join(f"h{l}: {np.hypot(pooled[2*l-1], pooled[2*l])/pooled[0]:.3f} @ {np.degrees(np.arctan2(pooled[2*l], pooled[2*l-1]))/(2*l):+6.1f} deg" for l in range(1, A.L + 1)) + "   (0 = uniform birth; 1 = fully modulated)")
    else:
        print("  min-miss LP (certificate by value):", flush=True)
        rm, meta2 = build_and_solve(C, A.L, A.K, A.scope, A.eps, "miss", A.method, A.tlim, M_per_L=A.Mpl, Krim=Krim, tol=A.tol)
        if rm.status == 0:
            u, v = verify(rm, meta2, 0.0); rec.update(status="infeasible (certified by value)" if v["band_miss"] > A.eps + 1e-6 else "feasible; price LP failed numerically", sigma=float(rm.x[meta2["off_t"]]), min_miss=v["band_miss"])
            print(f"  min miss {v['band_miss']:.6f} (eps {A.eps})  {rec['status']}")
        else: rec.update(status="undetermined (both LPs failed)")
    rec["Krim"] = A.K if Krim is None else Krim; jp = f"census_fuzz_v3_L{A.L}_K{A.K}{rtag}_{A.scope}.json"; json.dump(rec, open(jp, "w"), default=float)
    print(f"\nSUMMARY  census_fuzz_v3  L={A.L}  K={A.K}  Krim={rec['Krim']}  scope={A.scope}  eps={A.eps}  (positivity points {max(8, A.Mpl*A.L) if A.L else 1})")
    if "price" in rec: print(f"  MD/floor {rec['price']:.4f}   {rec['status']}" + (f"   clipped upper bound {rec['clip_md_over_floor']:.4f} (band miss {rec['clip_band_miss']:.5f})" if rec.get("clip_md_over_floor") else "") + (f"   rim depth " + " ".join(f"{k} {d:.3f}" for k, d in rec["depth"].items()) if rec.get("depth") else ""))
    else: print(f"  -   {rec['status']}" + (f"; min miss {rec['min_miss']:.5f}" if "min_miss" in rec else ""))
    print(f"  reference L=0 (paper, lattice): K=6 1.700, K=8 1.429, K=10 1.258, full 1.000      total {time.time()-t0:.0f}s   json {jp}")

if __name__ == "__main__":
    main()
