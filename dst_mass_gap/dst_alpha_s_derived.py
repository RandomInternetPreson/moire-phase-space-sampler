"""
DST-derived α_s running: UPDATED with geometric measure 2π/ln(2π⁵)

Previous version: SU(3) geometric measure ≈ 1 → α_s(m_Z) = 0.133 (13% error)
This version:     SU(3) geometric measure = 2π/ln(2π⁵) → α_s(m_Z) = 0.116 (1.6% error)

The structural parallel:
  EM:   α   × L = 2π / ((4/3) × Vol(S²))   = 3/8      [volume enters linearly]
  QCD:  α_s × L = 2π / ln(Vol(SU(3)))       = 0.979    [volume enters logarithmically]
  
  Linear vs logarithmic ← commutativity vs non-commutativity
"""
import numpy as np

L = np.log(2.176e-8 / 9.109e-31)  # ln(m_Pl/m_e) = 51.528
Vol_SU3 = 2 * np.pi**5
ln_Vol = np.log(Vol_SU3)
measure_QCD = 2 * np.pi / ln_Vol
alpha_s_planck = measure_QCD / L

print("DST Strong Coupling: α_s × L = 2π / ln(2π⁵)")
print("=" * 65)
print(f"  SU(3) measure = 2π/ln(2π⁵) = {measure_QCD:.6f}")
print(f"  α_s(m_Pl) = {alpha_s_planck:.6f}")
print()

inv_alpha = 1.0 / alpha_s_planck
thresholds = [
    (1.22e19, 173.0, 6, "m_t", None),
    (173.0, 91.2, 5, "m_Z", 0.1179),
    (91.2, 4.18, 5, "m_b", 0.22),
    (4.18, 1.27, 4, "m_c", None),
    (1.27, 1.0, 3, "1 GeV", None),
    (1.0, 0.5, 3, "500 MeV", None),
]

print(f"{'Scale':>12} {'α_s':>8} {'1/α_s':>8}  Notes")
print(f"{'m_Pl':>12} {alpha_s_planck:8.4f} {1/alpha_s_planck:8.2f}  DST: 2π/ln(2π⁵)")

for mu_hi, mu_lo, nf, label, obs in thresholds:
    b0 = 11 - 2*nf/3
    inv_alpha += (b0/(2*np.pi)) * np.log(mu_lo/mu_hi)
    if inv_alpha <= 0:
        print(f"{label:>12} {'div':>8} {'<0':>8}  Landau pole")
        break
    a = 1.0/inv_alpha
    note = f"obs: {obs} → err={100*(a-obs)/obs:+.1f}%" if obs else ""
    if "1 GeV" in label: note = "NUCLEAR SCALE"
    print(f"{label:>12} {a:8.4f} {inv_alpha:8.2f}  {note}")

print(f"\nα_s(1 GeV) ≈ 0.34 — use in mass gap bootstrap")
