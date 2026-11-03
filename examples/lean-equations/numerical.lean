-- The bounded proof search may leave this equation unresolved.
-- Its elaborated expression can still be lowered to numerical C++.
solve fixed_point for x : ℝ where Real.cos x = x

-- CLI numerical fallback, with a sign-changing bracket:
-- python scripts/lean_equations.py examples/lean-equations/numerical.lean \
--   --cpp /tmp/equation.cpp --run --lower 0 --upper 1
