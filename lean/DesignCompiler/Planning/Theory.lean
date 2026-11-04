import DesignCompiler.Planning
import Mathlib.Probability.Variance

open MeasureTheory ProbabilityTheory
open scoped BigOperators

namespace DesignCompiler.Planning

/-- For centered random components, this is the covariance propagation identity.
    Integrability of cross-products makes exchanging finite sums and expectation valid. -/
theorem contrast_second_moment {Ω ι : Type*} [MeasurableSpace Ω] [Fintype ι]
    (μ : Measure Ω) (X : ι → Ω → ℝ) (c : ι → ℝ)
    (h : ∀ i j, Integrable (fun ω => X i ω * X j ω) μ) :
    (∫ ω, (∑ i, c i * X i ω)^2 ∂μ) =
      ∑ i, ∑ j, c i * c j * ∫ ω, X i ω * X j ω ∂μ := by
  have expansion : ∀ ω, (∑ i, c i * X i ω)^2 =
      ∑ i, ∑ j, (c i * c j) * (X i ω * X j ω) := by
    intro ω
    simp only [pow_two, Finset.sum_mul, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    apply Finset.sum_congr rfl
    intro j _
    ring
  simp_rw [expansion]
  rw [integral_finset_sum]
  · apply Finset.sum_congr rfl
    intro i _
    rw [integral_finset_sum]
    · simp only [integral_const_mul]
    · intro j _; exact (h i j).const_mul _
  · intro i _
    exact integrable_finset_sum _ (fun j _ => (h i j).const_mul _)

/-- Averaging n independent sampling blocks gives the V/n scaling used by the
    planning IR. The blocks need not have Gaussian outcomes for this variance identity. -/
theorem independent_average_variance {Ω : Type*} [MeasurableSpace Ω]
    (μ : Measure Ω) [IsProbabilityMeasure μ] (n : Nat) (hn : 0 < n)
    (X : Fin n → Ω → ℝ) (v : ℝ)
    (hX : ∀ i, MemLp (X i) 2 μ)
    (hind : Set.Pairwise (Set.univ : Set (Fin n)) fun i j => IndepFun (X i) (X j) μ)
    (hv : ∀ i, ProbabilityTheory.variance (X i) μ = v) :
    ProbabilityTheory.variance (fun ω => (1 / (n : ℝ)) * ∑ i, X i ω) μ = v / n := by
  rw [variance_mul]
  have hs : ProbabilityTheory.variance (∑ i, X i) μ = ∑ i : Fin n, ProbabilityTheory.variance (X i) μ := by
    apply IndepFun.variance_sum
    · intro i _; exact hX i
    · simpa using hind
  have hfun : (fun ω => ∑ i, X i ω) = ∑ i, X i := by funext ω; simp
  rw [hfun, hs]
  simp only [hv, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast Nat.ne_of_gt hn
  field_simp
  ring

end DesignCompiler.Planning
