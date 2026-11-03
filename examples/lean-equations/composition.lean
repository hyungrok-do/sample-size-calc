-- This is Lean code. The CLI supplies `import DesignCompiler.Solve`.
noncomputable def response (x : ℝ) : ℝ := Real.exp (2 * x + 1)

solve response_solution for x : ℝ where response x = 5

-- The generated theorem is available to subsequent Lean code.
#check response_solution
