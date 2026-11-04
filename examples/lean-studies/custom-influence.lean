open DesignCompiler.Planning

-- Define a variance from a pilot influence-function distribution.
-- This is general Wald planning, conditional on the representation you supply.
def pilotWeights : List ℚ := [1/4, 1/2, 1/4]
def pilotInfluences : List ℚ := [-2, 0, 2]
def influenceVariance : ℚ := dot pilotWeights (pilotInfluences.map (·^2))

def study : Study := {
  name := "User-specified influence-function design"
  estimand := "Mean contrast on the chosen effect scale"
  analysis := .wald
  model := .asymptoticWald {
    alternative := 1/4
    variance := influenceVariance
    varianceSource := "Second moment of the specified centered pilot influence distribution"
    assumptions := ["The supplied pilot influence distribution transports to the planning population.",
                    "The planned estimator has this centered influence representation."]
  }
  criterion := { minN := 20, maxN := 10000 }
}

solve sample_size custom_design using study
