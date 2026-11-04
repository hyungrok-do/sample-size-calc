-- No assumption that a is nonzero: preserve empty and all-real solution sets.
variable (a y : ℝ)

solve conditional_solution for x : ℝ where a * x = y

-- The same C++ executable handles all three cases:
-- ./equation a=2 y=6   -> x = 3
-- ./equation a=0 y=6   -> no roots
-- ./equation a=0 y=0   -> all real values
