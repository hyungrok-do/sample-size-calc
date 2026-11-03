section
variable (a b y : ℝ) (ha : a ≠ 0) (hy : 0 < y)

solve parameter_solution for x : ℝ where Real.exp (a * x + b) = y

#check parameter_solution
end
