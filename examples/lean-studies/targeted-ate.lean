open DesignCompiler.Planning

-- The binary-outcome variances are computed by ordinary Lean functions.
def bernoulliVariance (p : ℚ) : ℚ := p*(1-p)

def study : Study := {
  name := "Covariate-adjusted ATE planning with TMLE"
  estimand := "Population-average risk difference"
  analysis := .tmle
  model := .efficientATE [
    { weight := 3/5, propensity := 1/2, mean0 := 1/10, mean1 := 1/5,
      variance0 := bernoulliVariance (1/10), variance1 := bernoulliVariance (1/5) },
    { weight := 2/5, propensity := 1/2, mean0 := 2/5, mean1 := 11/20,
      variance0 := bernoulliVariance (2/5), variance1 := bernoulliVariance (11/20) }
  ]
  criterion := { alpha := 1/20, targetPower := 9/10, minN := 20, maxN := 10000 }
}

solve sample_size targeted_ate using study
