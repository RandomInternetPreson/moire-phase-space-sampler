"""
The superselected eta calculation (E1 residual, toy-exact).

Model: qubit + drive mode, number-conserving Jaynes-Cummings exchange
  U: |0,n> -> cos(g t sqrt(n)) |0,n> - i sin(g t sqrt(n)) |1,n-1>
Drive: superselected laser = Poissonian number mixture, mean Nbar
(no absolute phase exists -- H1 applied to the drive).

Gate: pi/2 pulse calibrated at the mean, g t = pi/(4 sqrt(Nbar)).

Object: the RELATIONAL winding coherence of qubit-with-drive,
  C = |Tr[rho sigma^+ (x) d]| / sqrt(n)-normalized,
whose ideal value is 1/2. The deficit
  eta^2 := 1/2 - C  (record strength^2 per unit winding exchanged)
is the partial-winding record the drive holds.

Result (machine-established, correcting a factor-1/2 analytic slip:
the deficit is (1/2)(1 - <sin 2theta>), giving pi^2/64, not pi^2/32):
  eta^2 = pi^2/(64 Nbar) + O(1/Nbar^2)  ~ 0.15421/Nbar.
"""
import numpy as np
from math import pi, sqrt, exp, lgamma

def poisson_pmf(n, mu):
    return exp(n*np.log(mu) - mu - lgamma(n+1))

def deficit(Nbar):
    # sum over number branches; normalized relational coherence per branch
    #   <sigma^+ d> matrix element magnitude = (1/2) sin(2 theta_n) * sqrt(n),
    # normalized by sqrt(n): contribution (1/2) sin(2 theta_n)
    gt = pi / (4.0 * sqrt(Nbar))
    lo = max(0, int(Nbar - 12*sqrt(Nbar))); hi = int(Nbar + 12*sqrt(Nbar)) + 2
    C = 0.0; Z = 0.0
    for n in range(max(lo,1), hi):
        p = poisson_pmf(n, Nbar)
        C += p * 0.5 * np.sin(2 * gt * sqrt(n))
        Z += p
    C /= Z
    return 0.5 - C

print(f"{'Nbar':>8} {'eta^2':>12} {'eta^2 * Nbar':>14} {'eta':>10}")
for Nbar in [25, 100, 400, 1600, 6400, 25600, 102400]:
    d = deficit(Nbar)
    print(f"{Nbar:>8} {d:>12.6e} {d*Nbar:>14.6f} {sqrt(d):>10.4e}")

print(f"\nanalytic coefficient pi^2/64 = {pi**2/64:.6f}")

# scaling exponent from the two largest Nbar
import numpy as np
Ns = np.array([1600, 6400, 25600, 102400], float)
ds = np.array([deficit(N) for N in Ns])
slope = np.polyfit(np.log(Ns), np.log(ds), 1)[0]
print(f"log-log slope of eta^2 vs Nbar: {slope:.4f}  (claim: -1)")

# assertions (suite discipline)
assert abs(slope + 1.0) < 1e-3, f"slope {slope} != -1"
for N in [1600, 6400, 25600, 102400]:
    coeff = deficit(N) * N
    assert abs(coeff - pi**2/64) < 2e-4, f"coefficient {coeff} != pi^2/64"
print("ALL ASSERTIONS PASSED")
