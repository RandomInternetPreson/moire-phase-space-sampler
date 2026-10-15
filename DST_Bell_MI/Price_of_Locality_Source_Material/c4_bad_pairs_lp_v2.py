"""c4_bad_pairs_lp_v2.py (v2: asserts LP value 1/3 within 1e-7 on every pair run) -- faithful strategy-level LP for the 12 padded-C4
scenarios on which the ground-state construction is not (feasible with F=1/3).
Usage: python c4_bad_pairs_lp_v2.py <start> <stride>   (for parallel runs)"""
import sys, functools, importlib.util, itertools, time
from fractions import Fraction
print = functools.partial(print, flush=True)
spec = importlib.util.spec_from_file_location('m', 'pol_scope_checks_v3.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
import numpy as np
bad = [('ZXZI','XZZX'),('ZXIX','XZZX'),('XZII','XIXZ'),('XZII','XIYY'),('XZZX','XIXZ'),('XZZX','XIYY'),
       ('XZZX','YYZI'),('XZZX','YYIX'),('XZZX','YXYZ'),('XZZX','YXXY'),('XIXZ','XIYY'),('YXYZ','YXXY')]
start, stride = int(sys.argv[1]), int(sys.argv[2])
# rebuild the pieces of c4_scenario we need
n=4; PAULI=m.PAULI
gens=[(1,(1,3,0,0)),(1,(3,1,3,0)),(1,(0,3,1,3)),(1,(0,0,3,1))]; grp=m.stabilizer_group(gens)
core=[('IZXZ',+1),('IZYY',+1),('ZYYZ',+1),('ZYXY',-1)]
strat=np.array(list(itertools.product(*[list(itertools.product((1,-1),repeat=3))]*n)))
def run(settings):
    S=len(settings); eqs=[]; full=[]
    for name,sg in settings:
        t=tuple(PAULI[ch] for ch in name); parties=[i for i in range(n) if t[i]]
        rows=[np.ones(len(strat))]; rhs=[1.0]
        for r in range(1,len(parties)+1):
            for T in itertools.combinations(parties,r):
                sub=tuple(t[i] if i in T else 0 for i in range(n)); qv=grp.get(sub,0)
                prod=np.ones(len(strat))
                for i in T: prod=prod*strat[:,i,t[i]-1]
                rows.append(prod); rhs.append(float(qv))
                if r==len(parties): full.append(prod==sg)
        eqs.append((np.array(rows),rhs))
    F,_=m.solve_minmax_l1(len(strat),eqs)
    H=(~np.array(full)).sum(axis=0); s=Fraction(S-int(H.min()),S)
    return F,s
for idx in range(start, len(bad), stride):
    a,b=bad[idx]; t0=time.time()
    sett=core+[(a,grp[tuple(PAULI[ch] for ch in a)]),(b,grp[tuple(PAULI[ch] for ch in b)])]
    F,s=run(sett)
    print(f"{a} {b}: LP F_min = {F:.9f} = {Fraction(F).limit_denominator(400)}, s = {s}  [{time.time()-t0:.0f}s]")
    assert abs(F - 1/3) < 1e-7 and s == Fraction(5, 6), (a, b, F, s)
print("ASSERTED: LP value 1/3 and s = 5/6 on every pair run")
