-- Users can extend equation solving by proving a theorem in Lean.
noncomputable def cubeLink (x : ℝ) : ℝ := x^3 + x

@[equation_rule] theorem cubeLink_zero (x : ℝ) : cubeLink x = 0 ↔ x = 0 := by
  unfold cubeLink
  constructor
  · intro h
    have factor : x * (x^2 + 1) = 0 := by nlinarith [h]
    rcases mul_eq_zero.mp factor with hzero | himpossible
    · exact hzero
    · nlinarith [sq_nonneg x]
  · intro h; simp [h]

solve custom_solution for x : ℝ where cubeLink (2*x - 3) = 0

#check custom_solution
