import SSC

-- The model supplies a covariance and the coefficients of an estimator's
-- linearization. Its variance is derived algebraically from those equations.
-- Establishing a linearization for a particular estimator is a distinct theorem.
namespace Linearization

def quadratic (v11 v12 v21 v22 g1 g2 : Rat) : Rat :=
  g1 * (v11 * g1 + v12 * g2) + g2 * (v21 * g1 + v22 * g2)

theorem symmetric_covariance (v11 v12 v22 g1 g2 : Rat) :
    quadratic v11 v12 v12 v22 g1 g2 =
      g1 * g1 * v11 + 2 * g1 * g2 * v12 + g2 * g2 * v22 := by
  unfold quadratic
  grind

-- For a ratio a/b, these are the supplied linearization coefficients.
-- They are ordinary definitions and can be replaced with another estimand.
def ratioVariance (a b v11 v12 v22 : Rat) : Rat :=
  quadratic v11 v12 v12 v22 (1 / b) (-a / (b * b))

theorem ratio_variance_expansion (a b v11 v12 v22 : Rat) :
    ratioVariance a b v11 v12 v22 =
      (1 / b) * (1 / b) * v11 +
      2 * (1 / b) * (-a / (b * b)) * v12 +
      (-a / (b * b)) * (-a / (b * b)) * v22 := by
  exact symmetric_covariance v11 v12 v22 (1 / b) (-a / (b * b))

#eval ratioVariance 2 4 9 1 4

end Linearization
