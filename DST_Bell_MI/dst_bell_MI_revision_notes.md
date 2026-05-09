# DST Bell/MI Paper — Revision Notes (v6 FINAL)

**Status:** All items complete. Paper is submittable.

Progression: `dst_bell_MI_v4_FINAL.docx` (input) → `v5` (draft) → **`v6`** (final, Barrett-Gisin closed via direct computation).

---

## What changed in v6 (from v5)

**Section 8.3 rewritten** — previously flagged the Barrett-Gisin geometric derivation as an open calculation; now closes it directly by computing the mutual information of Hall's density (which DST derives in Section 5).

### Key computation added to Section 8.3

For the DST-produced density ρ_XY(λ) with uniform marginal over settings on S²:

- **I(X, Y; λ) = 0.028 bits** (symmetric case, both settings averaged over uniform sphere-pair measure)
- **I(X; λ) = 0.0023 bits** (asymmetric case with Y marginalized uniformly — the Barrett-Gisin configuration)

Both sit far below the Barrett-Gisin upper bound of 1 bit (by factors of ~36 and ~430 respectively). The condensate mechanism that produces Hall's density therefore satisfies both Hall's variational-distance minimum (by saturation, Section 5) and Barrett-Gisin's mutual-information bound (by substantial margin) through a single derivation.

### Why this reframe is stronger than the "open calculation" framing

The v5 text treated Barrett-Gisin's 1-bit figure as a load-bearing number that needed a separate geometric derivation. On reexamination, the 1-bit figure is a sufficient upper bound from Barrett-Gisin's proof technique, not a fundamental quantity of any physical model. The real quantities are the actual MI values of the specific density DST produces — and those we can compute directly because DST fixes the density in Section 5.

The structural insight is that Hall's bound and Barrett-Gisin's bound aren't independent constraints DST has to satisfy separately. They're two different information-theoretic measurements (L¹ norm vs. relative entropy) on the same condensate-produced object. DST's derivation of the density simultaneously determines both values. This is the stronger version of the claim that was previously hedged.

### Side observation kept out of the paper

The average L¹ distance between ρ_XY(λ) and uniform has an exact analytical form:

    avg_L1 = (4 - π) / (2π) ≈ 0.1366

This sits numerically between Hall's sup-bound (√2−1)/3 ≈ 0.1381 and g₀² = 9/64 = 0.1406. Three geometric quantities within ~3% of each other, from three different operations on the same structure. Left out of v6 at Aaron's discretion to avoid muddying the headline claim.

---

## All v5 → v6 changes (cumulative from v4)

### ✅ 1. Hall vs. Barrett-Gisin citation clarity

Three separate results distinguished in the literature:

| Result | Authors, year | Measure | Bound |
|---|---|---|---|
| (√2−1)/3 ≈ 13.81% | Hall 2010, PRL 105, 250404 | Variational distance | Used in paper ✓ |
| ≤ 1 bit mutual info | Barrett & Gisin 2011, PRL 106, 100406 | Mutual information | Addressed in §8.3 ✓ |
| Kochen-Specker extension | Hall 2011, PRA 84, 022102 | Various | Cited as [6] |

### ✅ 2. Barrett-Gisin directly addressed in Section 8.3

Upgraded from "flagged as open" (v5) to "closed by direct computation" (v6). Section now has five substantive paragraphs covering Hall's role, BG's bound, DST's actual MI values, the structural reason for the BG margin, and the summary claim.

### ✅ 3. Section 6 tightened — multiplicative-coupling hypothesis

**Section 6.1 "Structure of the Two-Observer Correction"** explicitly labels the multiplicative form as a *structural hypothesis*, justifies why multiplicative-on-E, addresses why L enters a low-energy experiment.

**Section 6.3 "The Insertion Count Is Sharply Selective"** foregrounds the n=1,2,3 selectivity as empirical verification.

### ✅ 4. Section 4.1 chain-independence foreword

Moved to prominent position so readers encounter it before Step 3 of the density derivation.

### ✅ 5. Section 5.5 "A Conditional Reading"

Renamed and reframed to lead with the conditional framing. Resolves tension with Section 8.1 disclaimers.

### ✅ 6. Predictions sharpened (Section 7)

- Prediction 1: concrete numerical deviation (~4×10⁻³ at φ=π/4)
- Prediction 2: tripartite scaling as sharpest distinguishing test
- Prediction 3: loophole-free tests don't close MI loophole
- Prediction 4: three-channel test of g₀ (now paired with the BG-MI prediction in Section 8.3)

### ✅ 7. Polish items

- Wigner reference [11] cited in Section 5.5
- Tsirelson framing clarified as "geometric reframing" not re-derivation
- Equation font changed from Cambria Math to Cambria (fixes √2 line-break)

---

## Verification log

All numerical claims in v6 verified:

**Bell/MI geometric identities (algebraic):**
- g₀² = 9/64 = 0.140625 ✓
- Hall min = (√2−1)/3 = 0.138071... ✓
- Excess = 1.85% relative ✓
- f_agree(π/4) = 3/4 = 2g₀ ✓
- 4cos(π(1−2g₀)) = 2√2 ✓

**Two-observer correction (closure at 0.014%):**
- ε = 2g₀²/L = 0.00546 with L=51.53 ✓
- MI_corrected = 0.14064, agreement with g₀² to 0.014% ✓
- n=1: 0.90% error; n=2: 0.014%; n=3: 0.93% ✓

**Section 5 density verification:**
- Hall density reproduces ⟨AB⟩ = −cos φ to machine precision ✓
- DST density deviates from Hall at φ=π/4 by M ≈ 4×10⁻³ ✓

**NEW — Section 8.3 mutual information (numerical integration):**
- ρ_X(λ) marginal normalizes to 1.00 ✓
- I(X, Y; λ) = 0.0280 bits (averaged over uniform pair measure) ✓
- I(X; λ) = 0.00231 bits (asymmetric BG case) ✓
- χ² divergence = 0.00329, gives I ≈ 0.00237 via small-deviation approximation ✓ (agrees with direct)
- Max peak-to-trough in ρ_X: ~20% at Bell-optimal angle ✓
- Safety margin vs BG bound: ~36× (symmetric), ~430× (asymmetric) ✓

---

## Strategic framing

DST-companion-paper framing retained (appropriate for viXra:2604.0009 companion). The Trojan-horse reframing (lead with Hall↔α coincidence, layer DST as explanation) remains available if submitting to a Bell-foundations venue in the future.

---

## Open items remaining after v6

These are honestly flagged in the paper rather than closed gaps:

1. **First-principles derivation of the two-observer hypothesis.** Section 6.1 labels this a hypothesis; n-selectivity is the empirical test. Deriving the multiplicative two-observer form from the DST Lagrangian remains future work.

2. **Explicit numerical prediction for tripartite GHZ measurement dependence.** Section 7.2 states n=3 scaling qualitatively; quantitative prediction requires the multipartite analog of Hall's construction.

Section 8.3's Barrett-Gisin computation is no longer on this list — it's now closed by direct computation.

---

## Files delivered

- `dst_bell_MI_v6.docx` — final revised paper (14 pages)
- `dst_bell_MI_v6.pdf` — PDF render
- `dst_bell_MI_revision_notes.md` — this document
