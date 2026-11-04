-- Assume Var(estimator) = variance / n. Solve its continuous variance target.
-- This equation alone does not prove integer minimality or a power guarantee.
section
variable (variance targetVariance : ℝ) (hv : 0 < variance) (ht : 0 < targetVariance)

solve precision_size for n : ℝ where variance / n = targetVariance

#check precision_size
end

-- Compile once, then evaluate for different inputs:
-- python scripts/lean_equations.py examples/lean-equations/precision.lean \
--   --cpp equation.cpp --run --param variance=2 --param targetVariance=0.01
