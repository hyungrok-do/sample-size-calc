open DesignCompiler.Planning

-- n allocation blocks: one control and two treatment participants per block.
-- The treatment block mean has variance sigma1^2/2.
def controlSD : ℚ := 1
def treatmentSD : ℚ := 6/5

def study : Study := {
  name := "Unequal-allocation randomized trial"
  estimand := "Treatment minus control mean"
  analysis := .knownCovarianceZ
  sizeUnit := "allocation blocks (1 control : 2 treatment)"
  groups := [{ label := "Control", unitsPerN := 1 }, { label := "Treatment", unitsPerN := 2 }]
  model := .gaussianContrast {
    labels := ["Control outcome", "Mean of two treatment outcomes"]
    means := [0, 2/5]
    coefficients := [-1, 1]
    covariance := [[controlSD^2, 0], [0, treatmentSD^2/2]]
  }
}

solve sample_size randomized_trial using study
