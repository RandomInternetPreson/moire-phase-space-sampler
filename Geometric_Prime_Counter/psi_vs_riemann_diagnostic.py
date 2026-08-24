#!/usr/bin/env python3
"""
Diagnostic: is the geometric prime counter's error set by zero count,
or by the conversion step  pi(x) ~ psi(x)/log x ?

Run next to geometric_prime_counter_v2.py.  Requires numpy, scipy, mpmath.
This is the script that exposed the v2 conversion error (August 2026).
"""
import sys, importlib.util
from math import log, isqrt
import numpy as np, mpmath as mp

spec = importlib.util.spec_from_file_location("g", "geometric_prime_counter_v2.py")
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)

def psi_exact(N):                      # arithmetic; diagnosis only
    s = bytearray(b'\x01')*(N+1); s[0]=s[1]=0
    for i in range(2, isqrt(N)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return sum(log(p)*(int(log(N)/log(p))) for p in range(2, N+1) if s[p])

# prime-free scan of Z(t) for ~100 true zeros (the "global pass" of Sec. 5.3)
zs=[]; t=10.0; Zp=g.rs_Z(t)
while t < 250:
    t2=t+0.05; Z2=g.rs_Z(t2)
    if Zp*Z2 < 0:
        r=g.refine_zero((t+t2)/2, search_radius=0.1, step=0.02)
        if r: zs.append(r)
    t, Zp = t2, Z2
zs = sorted(set(round(z,9) for z in zs))

surf = g.build_surface(); geo = surf['approx_zeros']; newt = g.refine_all(geo, verbose=False)
print(f"surface: {len(geo)} candidates, {len(newt)} refined; scan: {len(zs)} zeros to t=250\n")

def J(x, zeros):
    """Riemann (1859): J(x) = li(x) - sum_rho li(x^rho) - log 2 + int_x^inf dt/(t(t^2-1)log t)."""
    L=log(x); v = mp.li(x) - mp.log(2)
    for gam in zeros: v -= 2*mp.re(mp.ei(mp.mpc(0.5,gam)*L))
    v += mp.quad(lambda u: 1/(u*(u**2-1)*mp.log(u)), [x, mp.inf])
    return float(v)

def pi_riemann(x, zeros):
    """pi(x) = J(x) - pi(x^1/2)/2 - pi(x^1/3)/3 - ...  (recursive; no Mobius, no primes)"""
    if x < 2: return 0.0
    v=J(x,zeros); k=2
    while x**(1/k) >= 2: v -= pi_riemann(x**(1/k), zeros)/k; k+=1
    return v

hdr=f"{'x':>9} {'pi(x)':>7} | {'x/logx':>7} {'psi_ex/logx':>11} {'17geo/logx':>10} {'100z/logx':>9} | {'Riem,17':>8} {'Riem,100':>8} | {'li(x)':>7}"
print(hdr); print('-'*len(hdr))
for x in [100,1000,10000,100000,1000000]:
    pe=g.sieve_exact(x); e=lambda v: f"{abs(v-pe)/pe*100:6.2f}%"
    print(f"{x:>9} {pe:>7} | {e(x/log(x)):>7} {e(psi_exact(x)/log(x)):>11} "
          f"{e(g.count_primes(x,newt)):>10} {e(g.count_primes(x,zs[:100])):>9} | "
          f"{e(pi_riemann(x,newt)):>8} {e(pi_riemann(x,zs[:100])):>8} | {e(float(mp.li(x))):>7}")
