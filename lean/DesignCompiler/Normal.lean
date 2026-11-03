import DesignCompiler.Spec
import Mathlib.Probability.Distributions.Gaussian
import Mathlib.Probability.CDF
import Mathlib.Tactic

namespace DesignCompiler
open MeasureTheory ProbabilityTheory Set

noncomputable def Phi (x : ℝ) : ℝ := cdf (gaussianReal 0 1) x

/-- This rejection rule differs from `abs (x + lambda) > z` only at the
lower boundary. A normal distribution puts zero mass at that boundary. -/
noncomputable def rejectionProbability (z lambda : ℝ) : ℝ :=
  (gaussianReal 0 1).real {x | x + lambda ≤ -z ∨ z < x + lambda}

/-- Exact rejection probability for a shifted standard-normal statistic.
The iid raw-data-to-statistic distribution theorem is a separate obligation. -/
theorem normal_rejection_formula (z lambda : ℝ) (hz : 0 ≤ z) :
    rejectionProbability z lambda = Phi (-z - lambda) + 1 - Phi (z - lambda) := by
  have he : {x : ℝ | x + lambda ≤ -z ∨ z < x + lambda} =
      Iic (-z - lambda) ∪ Ioi (z - lambda) := by
    ext x
    simp only [mem_setOf_eq, mem_union, mem_Iic, mem_Ioi]
    constructor <;> intro h
    · rcases h with h | h
      · left; linarith
      · right; linarith
    · rcases h with h | h
      · left; linarith
      · right; linarith
  have hd : Disjoint (Iic (-z - lambda)) (Ioi (z - lambda)) := by
    rw [Set.disjoint_left]
    intro x hx hy
    simp only [mem_Iic, mem_Ioi] at hx hy
    linarith
  unfold rejectionProbability Phi
  rw [he, measureReal_union hd measurableSet_Ioi]
  rw [← compl_Iic, measureReal_compl measurableSet_Iic]
  simp only [measureReal_univ_eq_one, cdf_eq_real]
  ring

/-- The normal source formula and the variance-coefficient IR have identical
noncentrality. This is an exact-real lowering theorem, not a floating-point one. -/
theorem normal_ncp_lowering (delta sigma n : ℝ) (hs : 0 < sigma) (_hn : 0 < n) :
    delta * Real.sqrt n / (sigma * Real.sqrt 2) =
      delta / Real.sqrt ((2 * sigma ^ 2) / n) := by
  rw [Real.sqrt_div (by positivity), Real.sqrt_mul (by positivity), Real.sqrt_sq hs.le]
  rw [div_div_eq_mul_div]
  ring

/-- The finite-stratum ATE variance decomposition is nonnegative under positivity.
Its identification as the planned estimator's influence variance is not assumed
to hold for arbitrary nuisance fitting procedures. -/
theorem ate_variance_term_nonnegative (w g v0 v1 tau psi : ℝ)
    (hw : 0 ≤ w) (hg0 : 0 < g) (hg1 : g < 1)
    (hv0 : 0 ≤ v0) (hv1 : 0 ≤ v1) :
    0 ≤ w * (v1 / g + v0 / (1 - g) + (tau - psi) ^ 2) := by
  apply mul_nonneg hw
  exact add_nonneg (add_nonneg (div_nonneg hv1 hg0.le)
    (div_nonneg hv0 (by linarith))) (sq_nonneg _)

end DesignCompiler
