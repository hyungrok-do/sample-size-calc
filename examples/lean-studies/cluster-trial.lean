open DesignCompiler.Planning

-- Known Gaussian covariance model for equal-sized clusters.
-- If covariance must be estimated, use a method accounting for that uncertainty.
def clusterSize : Nat := 20
def individualVariance : ℚ := 1
def icc : ℚ := 1/20
def clusterVariance : ℚ := clusterMeanVariance individualVariance icc clusterSize
  (by decide +kernel) (by decide +kernel) (by decide +kernel)

def study : Study := {
  name := "Parallel cluster-randomized trial"
  estimand := "Difference in mean outcome between treatment and control"
  analysis := .knownCovarianceZ
  sizeUnit := "clusters per arm"
  groups := [
    { label := "Control clusters", participantsPerUnit := clusterSize },
    { label := "Treatment clusters", participantsPerUnit := clusterSize }
  ]
  model := .gaussianContrast {
    labels := ["Control cluster mean", "Treatment cluster mean"]
    means := [0, 3/10]
    coefficients := [-1, 1]
    covariance := [[clusterVariance, 0], [0, clusterVariance]]
  }
  criterion := { minN := 2, maxN := 1000 }
}

solve sample_size cluster_trial using study
