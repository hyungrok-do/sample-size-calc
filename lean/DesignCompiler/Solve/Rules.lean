import DesignCompiler.Solve.Attribute
import Mathlib.Tactic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Real.Sqrt

namespace DesignCompiler.Equation

/- Each rule is an independently checked equivalence. The compiler composes
   these through Lean's rewrite machinery; it does not recognize source strings. -/
@[equation_rule] theorem add_right (u a b : ℝ) :
    u + a = b ↔ u = b - a := by constructor <;> intro h <;> linarith

@[equation_rule] theorem add_left (u a b : ℝ) :
    a + u = b ↔ u = b - a := by constructor <;> intro h <;> linarith

@[equation_rule] theorem sub_right (u a b : ℝ) :
    u - a = b ↔ u = b + a := by constructor <;> intro h <;> linarith

@[equation_rule] theorem sub_left (u a b : ℝ) :
    a - u = b ↔ u = a - b := by constructor <;> intro h <;> linarith

@[equation_rule] theorem mul_right (u a b : ℝ) (ha : a ≠ 0) :
    u * a = b ↔ u = b / a := by exact (eq_div_iff ha).symm

@[equation_rule] theorem mul_left (u a b : ℝ) (ha : a ≠ 0) :
    a * u = b ↔ u = b / a := by rw [mul_comm]; exact (eq_div_iff ha).symm

@[equation_rule] theorem div_right (u a b : ℝ) (ha : a ≠ 0) :
    u / a = b ↔ u = b * a := by exact div_eq_iff ha

@[equation_rule] theorem negative (u b : ℝ) :
    -u = b ↔ u = -b := by constructor <;> intro h <;> linarith

@[equation_rule] theorem exponential (u b : ℝ) (hb : 0 < b) :
    Real.exp u = b ↔ u = Real.log b := by
  constructor
  · intro h; rw [← h, Real.log_exp]
  · intro h; rw [h, Real.exp_log hb]

@[equation_rule] theorem logarithm (u b : ℝ) (hu : 0 < u) :
    Real.log u = b ↔ u = Real.exp b := by
  constructor
  · intro h; rw [← h, Real.exp_log hu]
  · intro h; rw [h, Real.log_exp]

@[equation_rule] theorem square (u b : ℝ) (hb : 0 ≤ b) :
    u ^ 2 = b ↔ u = Real.sqrt b ∨ u = -Real.sqrt b := by
  have hs := Real.sq_sqrt hb
  constructor
  · intro h
    have hf : (u - Real.sqrt b) * (u + Real.sqrt b) = 0 := by nlinarith
    rcases mul_eq_zero.mp hf with h | h
    · left; linarith
    · right; linarith
  · rintro (rfl | rfl) <;> nlinarith

@[equation_rule] theorem square_negative (u b : ℝ) (hb : b < 0) :
    u ^ 2 = b ↔ False := by
  constructor
  · intro h; nlinarith [sq_nonneg u]
  · intro h; exact False.elim h

@[equation_rule] theorem square_root (u b : ℝ) (hb : 0 < b) :
    Real.sqrt u = b ↔ u = b ^ 2 := by
  rw [Real.sqrt_eq_iff_mul_self_eq_of_pos hb, eq_comm, sq]

@[equation_rule] theorem absolute (u b : ℝ) (hb : 0 ≤ b) :
    |u| = b ↔ u = b ∨ u = -b := by
  constructor
  · intro h; rcases le_total 0 u with hu | hu
    · left; simpa [abs_of_nonneg hu] using h
    · right; rw [abs_of_nonpos hu] at h; linarith
  · rintro (rfl | rfl) <;> simp [abs_of_nonneg hb]

@[equation_rule] theorem product_zero (u v : ℝ) :
    u * v = 0 ↔ u = 0 ∨ v = 0 := mul_eq_zero

@[equation_rule] theorem exponential_cases (u b : ℝ) :
    Real.exp u = b ↔ 0 < b ∧ u = Real.log b := by
  constructor
  · intro h
    have hb : 0 < b := h ▸ Real.exp_pos u
    exact ⟨hb, (exponential u b hb).mp h⟩
  · rintro ⟨hb, h⟩; exact (exponential u b hb).mpr h

@[equation_rule] theorem square_cases (u b : ℝ) :
    u ^ 2 = b ↔ 0 ≤ b ∧ (u = Real.sqrt b ∨ u = -Real.sqrt b) := by
  constructor
  · intro h
    have hb : 0 ≤ b := h ▸ sq_nonneg u
    exact ⟨hb, (square u b hb).mp h⟩
  · rintro ⟨hb, h⟩; exact (square u b hb).mpr h

@[equation_rule] theorem mul_right_cases (u a b : ℝ) :
    u * a = b ↔ (a = 0 ∧ b = 0) ∨ (a ≠ 0 ∧ u = b / a) := by
  by_cases ha : a = 0
  · simp [ha, eq_comm]
  · simp [ha, mul_right u a b ha]

@[equation_rule] theorem mul_left_cases (u a b : ℝ) :
    a * u = b ↔ (a = 0 ∧ b = 0) ∨ (a ≠ 0 ∧ u = b / a) := by
  rw [mul_comm]; exact mul_right_cases u a b

@[equation_rule] theorem completing_square (u b c : ℝ) :
    u^2 + b*u = c ↔ (u + b/2)^2 = c + (b/2)^2 := by
  constructor <;> intro h <;> nlinarith

@[equation_rule] theorem completing_square_right (u b c : ℝ) :
    u^2 + u*b = c ↔ (u + b/2)^2 = c + (b/2)^2 := by
  rw [mul_comm u b]; exact completing_square u b c

@[equation_rule] theorem completing_scaled_square (u a b c : ℝ) (ha : a ≠ 0) :
    a*u^2 + b*u = c ↔ (2*a*u + b)^2 = b^2 + 4*a*c := by
  have identity : (2*a*u + b)^2 - (b^2 + 4*a*c) = 4*a*(a*u^2 + b*u - c) := by ring
  constructor
  · intro h
    have hz : a*u^2 + b*u - c = 0 := by rw [h]; ring
    rw [hz, mul_zero] at identity
    exact sub_eq_zero.mp identity
  · intro h
    rw [h, sub_self] at identity
    have hz := (mul_eq_zero.mp identity.symm).resolve_left (mul_ne_zero (by norm_num) ha)
    exact sub_eq_zero.mp hz

@[equation_rule] theorem move_rhs (u v : ℝ) :
    u = v ↔ u - v = 0 := sub_eq_zero.symm

end DesignCompiler.Equation
