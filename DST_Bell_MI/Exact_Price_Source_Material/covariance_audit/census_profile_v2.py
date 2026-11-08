#!/usr/bin/env python3
"""
census_profile_v2.py -- re-contract the stored lambda-histograms with a chosen DETECTION PROFILE (the bulge) and price the census.

The model's station is A = sign(cos 2(lambda - a)), B = sign(cos 2(lambda + pi/2 - b)): a step.  Here the step is replaced by
a profile D(theta) with |D| <= 1 and D(theta + pi/2) = -D(theta) (the other port), the pair's lambda-histograms (from
pipeline_bulge_v1.py) are contracted with D_A(lambda - a) D_B(lambda + pi/2 - b) at every stored angle, and the usual LPs run:
  * S_max   : the largest CHSH value any census can reach with this bulge (max over atoms of |Svec|); below 2 sqrt 2 the
              singlet is unreachable whatever the birth does;
  * miss    : the relaxed min-miss LP over the full census -- the smallest uniform miss of the 161 correlators and the four
              CHSH correlators together (no CHSH equality), always feasible;
  * full    : min_md_full with this bulge (identity basis): the floor price of the bulge itself;
  * K-curve : the band-limited price LPs of census_bandlimit_v4 (positivity at atoms, mass, CHSH exact, band eps);
              infeasible entries certified by value (min-miss LP), as in v4;
  * marginals: <A>(t) and <B>(t) under the priced census at all 161 angles (first time printed in this program; a diagnostic).
Profiles:  sign | smooth | tanh:B | notch:J:E | anotch:J:E   (see pipeline_bulge_v1.py)
           v2 adds  dome:P  -- the bead with a vertical wall and a curved cap: sign(cos 2theta) |cos 2theta|^P  (P -> 0 is the step,
           P = 1 the smooth bulge; the wall at 45 deg stays sharp for every P, only the cap bends).  Same outputs, json/log names carry the tag.  'sign' must reproduce the banks
file next to the histogram file bit-for-bit (checked and printed when that file exists).

Usage:  python census_profile_v2.py hist_nsh32_grid0_nang161_nw2048.npz tanh:4 --klist 6,8,10,12,14 --eps 0.01 --full 1
Writes census_profile_v2_<ptag>.json; prints a SUMMARY block.
"""
import os, sys, json, time, argparse
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
import census_bandlimit_v4 as cb                      # basis, reduce, price_lp, miss_lp, verify, md_of, run_lp

_deg = np.pi / 180; S_QM = 2 * np.sqrt(2); FLOOR = (np.sqrt(2.0) - 1.0) / 3.0
settings = [(0.0, 22.5), (0.0, 67.5), (45.0, 22.5), (45.0, 67.5)]; signs = [+1, -1, +1, +1]

def ptag(spec):
    f = spec.split(":")
    if f[0] in ("sign", "smooth"): return f[0]
    if f[0] == "tanh": return f"tanh{f[1].replace('.', 'p')}"
    if f[0] in ("notch", "anotch"): return f"{f[0]}{f[1]}_e{int(round(float(f[2]) * 100)):02d}"
    if f[0] == "dome": return f"dome{f[1].replace('.', 'p')}"
    raise ValueError(spec)

def profile(spec):
    """returns (D_A, D_B): functions of the station-relative angle theta (radians); both in [-1, 1], antisymmetric under +pi/2."""
    f = spec.split(":"); step = lambda th: np.sign(np.cos(2 * th))
    if f[0] == "sign": return step, step, "the paper's step  sign(cos 2 theta)"
    if f[0] == "smooth": return (lambda th: np.cos(2 * th)), (lambda th: np.cos(2 * th)), "the smooth bulge  cos 2 theta"
    if f[0] == "tanh":
        B = float(f[1]); fn = lambda th: np.tanh(B * np.cos(2 * th)) / np.tanh(B)
        return fn, fn, f"blunted step  tanh({B:g} cos 2 theta)/tanh({B:g})"
    if f[0] == "dome":
        Pw = float(f[1]); fn = lambda th: step(th) * np.abs(np.cos(2 * th)) ** Pw
        return fn, fn, f"bead with a vertical wall and a curved cap  sign(cos 2 theta) |cos 2 theta|^{Pw:g}"
    if f[0] in ("notch", "anotch"):
        J = int(f[1]); E = float(f[2])
        gA = lambda th: 1 - E * (1 + np.cos(4 * J * th)) / 2          # notches of depth E at theta = 0 mod pi/(2J)
        gB = (lambda th: 1 - E * (1 - np.cos(4 * J * th)) / 2) if f[0] == "anotch" else gA   # anti: notches between Alice's
        return (lambda th: step(th) * gA(th)), (lambda th: step(th) * gB(th)), \
               f"{'anti-' if f[0]=='anotch' else ''}notched step: sign(cos 2 theta) x (1 - {E:g}(1 {'-' if f[0]=='anotch' else '+'} cos {4*J} theta)/2) on {'Bob' if f[0]=='anotch' else 'both'}"
    raise ValueError(spec)

def relaxed_miss_lp(Er_all, Mr, P, tgt_all, nc, method, tlim):
    """min sigma  s.t.  |E_all - target_all| <= sigma (161 angles + the 4 CHSH correlators), mass 1, positivity.  Always feasible."""
    n = Er_all.shape[0]
    Ep = sp.hstack([sp.csr_matrix(Er_all), sp.csr_matrix(-np.ones((n, 1)))]).tocsr()
    En = sp.hstack([sp.csr_matrix(-Er_all), sp.csr_matrix(-np.ones((n, 1)))]).tocsr()
    Pp = sp.hstack([P, sp.csr_matrix((P.shape[0], 1))]).tocsr()
    Aub = sp.vstack([Ep, En, Pp]).tocsr(); bub = np.concatenate([tgt_all, -tgt_all, np.zeros(P.shape[0])])
    Aeq = sp.csr_matrix(np.concatenate([Mr, [0]])[None, :])
    cobj = np.zeros(nc + 1); cobj[nc] = 1
    return cb.run_lp(cobj, Aub, bub, Aeq, [1.0], [(None, None)] * nc + [(0, None)], method, tlim)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hist"); ap.add_argument("profile")
    ap.add_argument("--klist", default="6,8,10,12,14"); ap.add_argument("--eps", type=float, default=0.01); ap.add_argument("--full", type=int, default=1)
    ap.add_argument("--method", default="highs-ipm"); ap.add_argument("--method2", default="highs-ds"); ap.add_argument("--tlim", type=float, default=0)
    A = ap.parse_args(); t0 = time.time()
    DA, DB, desc = profile(A.profile); tag = ptag(A.profile)
    z = np.load(A.hist)
    Hpair, Hang = z["Hpair"], z["Hang"]                       # (nflat, 4, NV, G), (nflat, nang, NV, G) uint16 counts
    offs = z["offsets"]; angs = z["angs"]; nw = float(z["nw"]); G = int(z["G"])
    nflat, nang, NV = Hpair.shape[0], Hang.shape[1], Hpair.shape[2]; NT = nflat * NV
    DL = np.pi / G; lam_c = (np.arange(G) + 0.5) * np.pi / G
    print(f"census_profile_v2  profile {A.profile} [{tag}]: {desc}")
    print(f"histograms {os.path.basename(A.hist)}: {nflat} cells x {NV} offsets x ({nang} angles + 4 CHSH settings) x {G} bins, {int(nw)} counts each")
    print(f"K list {A.klist}  eps {A.eps}  full {A.full}  solver {A.method} (retry {A.method2})")
    # --- the profile itself, printed so a cold reader can see the bulge ---
    th = np.linspace(0, np.pi, 13, endpoint=False)
    print("profile on theta = 0..165 deg step 15:  D_A " + " ".join(f"{v:+.3f}" for v in DA(th)) + "\n" + " " * 39 + "D_B " + " ".join(f"{v:+.3f}" for v in DB(th)))
    assert np.all(np.abs(DA(lam_c)) <= 1 + 1e-12) and np.all(np.abs(DB(lam_c)) <= 1 + 1e-12), "profile exceeds [-1, 1]"

    # ---------- contraction ----------
    t1 = time.time()
    Ecoef = np.zeros((nang, NT)); mA = np.zeros((nang, NT)); mB = np.zeros((nang, NT))
    for ai, tdeg in enumerate(angs):
        kA = DA(lam_c); kB = DB(lam_c + np.pi / 2 - tdeg * _deg)
        Hf = Hang[:, ai].astype(np.float64).reshape(NT, G)
        Ecoef[ai] = Hf @ (kA * kB) / nw; mA[ai] = Hf @ kA / nw; mB[ai] = Hf @ kB / nw
    del Hang
    Qs = [Hpair[:, p].astype(np.float64).reshape(NT, G) / (nw * DL) for p in range(4)]           # densities, as the banks store them
    Ecs = np.zeros((4, NT)); Svec = np.zeros(NT)
    for p, (a, b) in enumerate(settings):
        kern = DA(lam_c - a * _deg) * DB(lam_c + np.pi / 2 - b * _deg)
        Ecs[p] = (Qs[p] @ kern) * DL; Svec += signs[p] * Ecs[p]
    print(f"contraction done in {time.time()-t1:.0f}s", flush=True)
    # --- check against the banks (only the step can match them) ---
    bankf = os.path.join(os.path.dirname(A.hist) or ".", os.path.basename(A.hist).replace("hist_", "banks_", 1))
    if os.path.exists(bankf):
        zb = np.load(bankf); nb_ = int(zb["nflat"])
        Sb = np.concatenate([zb[f"s{k}"] for k in range(nb_)]); Eb = np.hstack([zb[f"e{k}"] for k in range(nb_)])
        Qb = np.vstack([zb[f"q{k}_0"] for k in range(nb_)])
        print(f"against {os.path.basename(bankf)}:  max|Svec - banks| {np.abs(Svec-Sb).max():.2e}   max|Ecoef - banks| {np.abs(Ecoef-Eb).max():.2e}   "
              f"max|Q0 - banks| {np.abs(Qs[0]-Qb).max():.2e}   " + ("(the step must match to ~1e-15; the Q check must always match)" if A.profile == "sign" else "(a shaped bulge is expected to differ in Svec/Ecoef; Q must match)"))
    # ---------- S_max and the relaxed miss ----------
    smax = float(-Svec.min()); chsh_ok = Svec.min() <= -S_QM + 1e-9
    print(f"\nS_max with this bulge = {smax:.6f} = {smax/S_QM:.6f} x 2 sqrt 2   ->  the singlet's CHSH value is {'reachable' if chsh_ok else 'NOT reachable: no census attains S = 2 sqrt 2 with this bulge'}")
    tgt = -np.cos(2 * angs * _deg); tgt_c = np.array([-np.cos(2 * (a - b) * _deg) for a, b in settings])
    out = dict(profile=A.profile, tag=tag, desc=desc, hist=os.path.basename(A.hist), eps=A.eps, smax=smax, chsh_ok=bool(chsh_ok), results=[])
    jpath = f"census_profile_v2_{tag}.json"
    Phi_I = np.eye(NV)
    Qr, Sr, Er, Mr, P = cb.reduce(Qs, Svec, Ecoef, Phi_I, nflat, NV)
    Er_c = np.hstack([Ecs[:, ts * NV:(ts + 1) * NV] @ Phi_I for ts in range(nflat)])
    print(f"\nrelaxed min-miss LP over the full census ({nang} angles + 4 CHSH correlators, no CHSH equality):", flush=True)
    rm = relaxed_miss_lp(np.vstack([Er, Er_c]), Mr, P, np.concatenate([tgt, tgt_c]), NT, A.method, A.tlim)
    rho_miss = None
    if rm.status == 0:
        c = rm.x[:NT]; rho_miss = np.clip(np.concatenate([Phi_I @ c.reshape(nflat, NV)[ts] for ts in range(nflat)]), 0, None); rho_miss /= rho_miss.sum()
        E_all = np.concatenate([Ecoef @ rho_miss, Ecs @ rho_miss]); miss = float(np.abs(E_all - np.concatenate([tgt, tgt_c])).max())
        out["miss_relaxed"] = miss; out["miss_relaxed_md_over_floor"] = float(cb.md_of(Qs, rho_miss) / FLOOR)
        print(f"  the closest this bulge gets to the singlet, any census: max |E - target| = {miss:.6f} over {nang+4} correlators (that census's MD/floor {out['miss_relaxed_md_over_floor']:.4f})")
    else:
        print("  relaxed min-miss LP failed")

    # ---------- the price LPs ----------
    KL = [int(x) for x in A.klist.split(",") if x]
    plan = ([("full", None)] if A.full else []) + [(f"K={K}", K) for K in KL]
    rho_full = None
    for label, K in plan:
        if K is None: Phi = Phi_I; Qr_, Sr_, Er_, Mr_, P_ = Qr, Sr, Er, Mr, P
        else:
            Phi, _ = cb.basis(offs, K); Qr_, Sr_, Er_, Mr_, P_ = cb.reduce(Qs, Svec, Ecoef, Phi, nflat, NV)
        nb = Phi.shape[1]; nc = nflat * nb; rec = {"label": label, "K": K}
        print(f"\n=== {label}   ({nb} basis functions per cell, {nc} coefficients) ===", flush=True)
        if not chsh_ok:
            Erc_ = np.hstack([Ecs[:, ts * NV:(ts + 1) * NV] @ Phi for ts in range(nflat)])
            rm = relaxed_miss_lp(np.vstack([Er_, Erc_]), Mr_, P_, np.concatenate([tgt, tgt_c]), nc, A.method, A.tlim)
            if rm.status == 0:
                c = rm.x[:nc].reshape(nflat, nb); rho = np.clip(np.concatenate([Phi @ c[ts] for ts in range(nflat)]), 0, None); rho /= rho.sum()
                miss = float(np.abs(np.concatenate([Ecoef @ rho, Ecs @ rho]) - np.concatenate([tgt, tgt_c])).max())
                rec.update(status="singlet unreachable (S_max < 2 sqrt 2); relaxed min miss", min_miss=miss, md_over_floor=float(cb.md_of(Qs, rho) / FLOOR))
                print(f"  relaxed min miss {miss:.6f}  (MD/floor of that census {rec['md_over_floor']:.4f})")
            else: rec.update(status="relaxed min-miss LP failed")
            out["results"].append(rec); json.dump(out, open(jpath, "w"), default=float); continue
        print("  price LP:", flush=True)
        res, _ = cb.price_lp(Qr_, Sr_, Er_, Mr_, P_, angs, A.eps, nflat, nb, A.method, A.tlim)
        if res.status == 0:
            c = res.x[:nc].reshape(nflat, nb); rho, v = cb.verify(c, Phi, nflat, Qs, Er_, angs, A.eps)
            rec.update(status="feasible" if v["ok"] else "feasible-unverified", price=v["price"], verify=v, solver=A.method)
        else:
            print("  min-miss LP (certificate by value):", flush=True)
            rmm = cb.miss_lp(Sr_, Er_, Mr_, P_, angs, A.eps, nc, A.method, A.tlim)
            if rmm.status not in (0, 2): rmm = cb.miss_lp(Sr_, Er_, Mr_, P_, angs, A.eps, nc, A.method2, 7200)
            if rmm.status != 0:
                # CHSH equality unattainable in this band (or solver failure): fall back to the relaxed miss, which is always feasible
                Erc_ = np.hstack([Ecs[:, ts * NV:(ts + 1) * NV] @ Phi for ts in range(nflat)])
                rr = relaxed_miss_lp(np.vstack([Er_, Erc_]), Mr_, P_, np.concatenate([tgt, tgt_c]), nc, A.method, A.tlim)
                if rr.status == 0:
                    c = rr.x[:nc].reshape(nflat, nb); rho = np.clip(np.concatenate([Phi @ c[ts] for ts in range(nflat)]), 0, None); rho /= rho.sum()
                    miss = float(np.abs(np.concatenate([Ecoef @ rho, Ecs @ rho]) - np.concatenate([tgt, tgt_c])).max())
                    rec.update(status="CHSH value unreachable in this band; relaxed min miss", min_miss=miss, md_over_floor=float(cb.md_of(Qs, rho) / FLOOR))
                    print(f"  CHSH equality infeasible in this band -> relaxed min miss {miss:.6f} (MD/floor of that census {rec['md_over_floor']:.4f})")
                else: rec.update(status="undetermined (min-miss and relaxed min-miss LPs both failed)")
            else:
                c = rmm.x[:nc].reshape(nflat, nb); rho, v = cb.verify(c, Phi, nflat, Qs, Er_, angs, 0.0, band_required=False)
                miss = v["band_miss"]; rec.update(sigma=float(rmm.x[nc]), min_miss=miss)
                if miss > A.eps + 1e-6: rec.update(status="infeasible (certified by value)"); print(f"  INFEASIBLE at eps={A.eps}: min miss {miss:.5f} > eps")
                else:
                    res2, _ = cb.price_lp(Qr_, Sr_, Er_, Mr_, P_, angs, A.eps, nflat, nb, A.method2, 7200)
                    if res2.status == 0:
                        c = res2.x[:nc].reshape(nflat, nb); rho, v = cb.verify(c, Phi, nflat, Qs, Er_, angs, A.eps)
                        rec.update(status="feasible" if v["ok"] else "feasible-unverified", price=v["price"], verify=v, solver=A.method2)
                    else: rec.update(status=f"feasible, price undetermined; upper bound {v['price']:.4f}", upper_bound_price=v["price"])
        if "price" in rec:
            print(f"  MD/floor = {rec['price']:.4f}   {rec['status']}   [{time.time()-t0:.0f}s total]")
            if K is None: rho_full = rho
            elif K <= 8:
                pix = cb.spectrum(rho, offs, nflat, NV, K); rec["leading"] = [(float(a), int(ka), int(kb), float(ph)) for a, ka, kb, ph in pix[:10]]
                print("  leading harmonics (|F|, k_A, k_B, phi):  " + "  ".join(f"{a:.3f}({ka:+d},{kb:+d}){ph:+.0f}" for a, ka, kb, ph in pix[:5]))
        out["results"].append(rec); json.dump(out, open(jpath, "w"), default=float)

    # ---------- marginals under the priced census, raw and symmetrized ----------
    # exact symmetry of the model: (xi_A, xi_B, lambda) -> each + 90 deg leaves the dynamics invariant and flips both detector
    # signs; so averaging a census with its (+90, +90)-shifted twin keeps every correlator and the MD and kills the marginals
    # -- when the twin offsets are in the sampling (the lattice).  Both versions are printed.
    rho_m = rho_full if rho_full is not None else rho_miss
    if rho_m is not None:
        src = "full price census" if rho_full is not None else "relaxed min-miss census"
        def stats(r):
            EA = mA @ r; EB = mB @ r
            return dict(A_max=float(np.abs(EA).max()), B_max=float(np.abs(EB).max()), A_spread=float(EA.max() - EA.min()), B_spread=float(EB.max() - EB.min()), A=EA.tolist(), B=EB.tolist())
        out["marginals"] = dict(source=src, **stats(rho_m))
        print(f"\nmarginals under the {src} (diagnostic; <A>(t) with Alice at 0 and Bob at t, {nang} angles):")
        print(f"  raw vertex:  max |<A>| {out['marginals']['A_max']:.4f}  spread over Bob's angle {out['marginals']['A_spread']:.4f}   max |<B>| {out['marginals']['B_max']:.4f}  spread {out['marginals']['B_spread']:.4f}")
        sh = np.array([[(x + np.pi / 2) % np.pi, (y + np.pi / 2) % np.pi] for x, y in offs])
        partner = np.full(NV, -1)
        for v in range(NV):
            dd = np.abs(((offs - sh[v]) + np.pi / 2) % np.pi - np.pi / 2).max(axis=1); j = int(np.argmin(dd))
            if dd[j] < 1e-9: partner[v] = j
        if np.all(partner >= 0):
            r2 = rho_m.reshape(nflat, NV); rho_sym = 0.5 * (r2 + r2[:, partner]).reshape(-1)
            s2 = stats(rho_sym)
            E_raw = np.concatenate([Ecoef @ rho_m, Ecs @ rho_m]); E_sym = np.concatenate([Ecoef @ rho_sym, Ecs @ rho_sym])
            out["marginals_sym"] = dict(**s2, md_over_floor=float(cb.md_of(Qs, rho_sym) / FLOOR), md_over_floor_raw=float(cb.md_of(Qs, rho_m) / FLOOR),
                                        max_dE=float(np.abs(E_raw - E_sym).max()))
            print(f"  symmetrized (census averaged with its +90/+90 twin):  max |<A>| {s2['A_max']:.2e}  max |<B>| {s2['B_max']:.2e}   "
                  f"MD/floor {out['marginals_sym']['md_over_floor']:.6f} (raw {out['marginals_sym']['md_over_floor_raw']:.6f})   max change in any correlator {out['marginals_sym']['max_dE']:.2e}")
            print("  (the symmetrized census is an equally good census with zero marginals: the raw marginals are the vertex's, not the model's)")
        else:
            print(f"  ({int(np.sum(partner < 0))} offsets have no +90/+90 twin in this sampling: symmetrization not available off-lattice; raw marginals are vertex-dependent)")
        print("  (an unpolarized source has <A> = <B> = 0; a spread of <A> over Bob's angle would be setting-dependence of a marginal)")
    json.dump(out, open(jpath, "w"), default=float)

    print(f"\nSUMMARY  profile {A.profile} [{tag}]  {desc}")
    print(f"  S_max/2rt2 {smax/S_QM:.6f}   relaxed miss {out.get('miss_relaxed', float('nan')):.6f}")
    for r in out["results"]:
        if "price" in r: print(f"  {r['label']:>6} {r['price']:9.4f}   {r['status']} [{r.get('solver','')}]")
        elif "min_miss" in r: print(f"  {r['label']:>6} {'-':>9}   {r['status']}; min miss {r['min_miss']:.5f}")
        else: print(f"  {r['label']:>6} {'-':>9}   {r['status']}")
    if "marginals" in out: print(f"  marginals: raw vertex max|<A>| {out['marginals']['A_max']:.4f} (spread {out['marginals']['A_spread']:.4f})  max|<B>| {out['marginals']['B_max']:.4f}" + (f";  symmetrized max|<A>| {out['marginals_sym']['A_max']:.1e}, MD/floor {out['marginals_sym']['md_over_floor']:.4f} vs raw {out['marginals_sym']['md_over_floor_raw']:.4f}, correlators moved {out['marginals_sym']['max_dE']:.1e}" if "marginals_sym" in out else "  (no twin offsets: not symmetrizable)"))
    print(f"  total {time.time()-t0:.0f}s   json {jpath}")

if __name__ == "__main__":
    main()
