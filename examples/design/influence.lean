import DesignCompiler.Input
open DesignCompiler

def study : InfluenceInput := {
  name := "ATE with influence-function variance"
  alternative := 0.5
  nullValue := 0
  variance := 4
  unit := .participant
  estimand := "ATE"
  estimator := "AIPW"
  varianceSource := "Balanced trial, residual variance 1, constant effect"
  alpha := 0.05
  targetPower := 0.80
  minN := 1
  maxN := 10000
}
