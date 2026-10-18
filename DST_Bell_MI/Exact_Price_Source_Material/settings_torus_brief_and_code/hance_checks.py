#!/usr/bin/env python3
"""
hance_checks.py - numpy port of the independent reviewer's Mathematica checks
(alai-brief-checks.wl, sections 3-5), October 2026. Reproduces: the 12.4% fringe's
rms deviation at Brida et al.'s eight angles; the (2,-2) log-Fourier pixel produced by
constant and angle-dependent two-channel efficiency imbalance (-delta^4/4 and
-d1 e1/4); the multi-pair identity; and the Larsson / symmetrised Gisin-Gisin
flat-rate efficiency. Each printed value is compared with the reviewer's figure.
The reviewer's own script is not redistributed here.
"""
import numpy as np
pi=np.pi
# ---- Sec 5: fringe model from the §5 profile (e1 = 0.062 amplitude, overtone e1/6) ----
e1=0.062; e2=e1/6; V=0.978
R4=lambda d: 1+e1*np.cos(4*d)+e2*np.cos(8*d)
print("R4(0)/R4(45) =",R4(0)/R4(pi/4),"  cross-ratio =",(R4(0)/R4(pi/4))**2,"  R4(22.5)/mean =",R4(pi/8))
ang=np.arange(1,9)*pi/8
r12=(1+V*np.cos(2*ang))*R4(ang)
m=r12.mean(); Vfit=2*(r12*np.cos(2*ang)).sum()/r12.sum()
print("Brida rms deviation, 12.4% model:",np.sqrt(np.mean((r12/m-1-Vfit*np.cos(2*ang))**2)))
# ---- Sec 3: two-channel imbalance -> (2,-2) pixel, numeric FFT on a torus ----
M=96; g=pi*np.arange(M)/M; A,B=np.meshgrid(g,g,indexing='ij')
def pix22(dA,dB,vis=1.0):
    R=0.25*(1+ -vis*dA(A)*dB(B)*np.cos(2*(A-B)))   # S_A=S_B=1, E=-vis cos2(a-b), D/S=delta
    F=np.fft.fft2(np.log(R)); return 2*np.real(F[2,M-2])/M**2, 2*np.imag(F[2,M-2])/M**2
for d0 in (0.1,0.2):
    c,_=pix22(lambda a:d0+0*a,lambda b:d0+0*b); print(f"const imbalance {d0}: cos4 amp={c:+.2e}  predicted -d0^4/4={-d0**4/4:+.2e}")
for d1 in (0.01,0.03):
    c,q=pix22(lambda a:0.05+d1*np.cos(2*a),lambda b:0.05+d1*np.cos(2*b)); print(f"ripple {d1} (p=q): cos4 amp={c:+.2e} quad={q:+.1e}  predicted -d1^2/4={-d1**2/4:+.2e}")
c,q=pix22(lambda a:0.05+0.03*np.cos(2*a),lambda b:0.05+0.03*np.cos(2*b-pi/2)); print(f"ripple 3% p-q=90deg: cos4 amp={c:+.2e} quad={q:+.2e}")
# ---- Sec 3b: two pairs in one window ----
d=np.linspace(0,pi,7); c=np.cos(2*d); ENN=9/4+c**2/4
print("E[N_A N_B] - (19/8 + cos4d/8) max dev:",np.max(np.abs(ENN-(19/8+np.cos(4*d)/8))))
# ---- Sec 4: symmetrised Gisin-Gisin efficiency ----
print("Gisin-Gisin eta = 4/(pi+2) =",4/(pi+2),"  Larsson 4/(2+pi) =",4/(2+pi))
# ---- 2x2 counting: 5 sigma at 0.1% peak-to-peak on the cross-ratio ----
# ratio r = R(0)^2/R(45)^2 ; fringe amp a=0.0005 -> r-1 ~ 4a ; var(log r) = 4/N per setting
a=0.0005; N=2.5e7; print("2x2: signal 4a =",4*a," sigma =",np.sqrt(4/N),"  -> sigma count =",4*a/np.sqrt(4/N))
