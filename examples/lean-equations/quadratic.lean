solve quadratic_solution for x : ℝ where x^2 + 2*x = 8

-- Use the generated theorem in an ordinary Lean proof.
example (x : ℝ) (h : x^2 + 2*x = 8) : x = 2 ∨ x = -4 := by
  exact (quadratic_solution x).mp h
