"""pixel_null_v1.py -- per-pixel permutation null for |F(2,-2)| (and (2,2) for reference),
same measure, same seed, same independent re-pairing as census_structure_v2.py (NPERM=400, SEED=20261015)."""
import json, sys, numpy as np
TAUS=[0.5,1.0,2.0,4.0,8.0]; SCALES=[1.0,1.8,2.6]; NPERM=400; SEED=20261015
def F(xa,xb,w,ka,kb): return np.sum(w*np.exp(2j*(ka*xa+kb*xb)))
print(f"{'census':8s} {'|F(2,-2)|':>9s} {'phase':>7s} {'null95':>7s} {'null99':>7s} {'z':>6s} {'p':>7s}   {'|F(2,2)|':>8s} {'null95':>7s}   1/sqrt(n_eff)")
for name in ['random0','random1','strat0','strat1']:
    d=json.load(open(f'rho_nsh32_{name}_eps0.01_nang161.json'))
    rho=np.asarray(d['rho'],float); offs=np.asarray(d['offsets'],float)
    nv=len(offs); rho=rho/rho.sum(); R=rho.reshape(len(TAUS),len(SCALES),nv); wo=R.sum(axis=(0,1))
    xa,xb=offs[:,0],offs[:,1]
    f22m=F(xa,xb,wo,2,-2); f22=F(xa,xb,wo,2,2)
    rng=np.random.default_rng(SEED); n_m=[]; n_p=[]
    for _ in range(NPERM):
        perm=rng.permutation(nv); xa_p,xb_p=xa[perm],xb[rng.permutation(nv)]
        n_m.append(abs(F(xa_p,xb_p,wo,2,-2))); n_p.append(abs(F(xa_p,xb_p,wo,2,2)))
    n_m=np.array(n_m); n_p=np.array(n_p)
    z=(abs(f22m)-n_m.mean())/n_m.std(); p=(n_m>=abs(f22m)).mean()
    print(f"{name:8s} {abs(f22m):9.3f} {np.degrees(np.angle(f22m)):+7.1f} {np.quantile(n_m,0.95):7.3f} {np.quantile(n_m,0.99):7.3f} {z:6.2f} {p:7.4f}   {abs(f22):8.3f} {np.quantile(n_p,0.95):7.3f}   {1/np.sqrt(1/np.sum(wo**2)):.3f}")
