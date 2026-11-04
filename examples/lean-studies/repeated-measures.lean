open DesignCompiler.Planning

-- Two arms, two visits. n counts participants PER ARM.
-- Sigma includes within-participant covariance; arms are independent.
def study : Study := {
  name := "Treatment-by-visit contrast"
  estimand := "Between-arm difference in mean change"
  analysis := .knownCovarianceZ
  sizeUnit := "participants per arm"
  groups := [{ label := "Control" }, { label := "Treatment" }]
  model := .gaussianContrast {
    labels := ["Control baseline", "Control follow-up", "Treatment baseline", "Treatment follow-up"]
    means := [0, 0, 0, 1/2]
    coefficients := [1, -1, -1, 1]
    covariance := [[1, 3/5, 0, 0], [3/5, 1, 0, 0],
                   [0, 0, 1, 3/5], [0, 0, 3/5, 1]]
  }
  criterion := { alpha := 1/20, targetPower := 4/5, minN := 2, maxN := 10000 }
}

solve sample_size repeated_measures using study
-- Other targets use the same model:
-- solve power at_fifty using study at 50
-- solve effect detectable_change using study at 50
