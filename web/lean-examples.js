export const examples = {
  normal: `import DesignCompiler.Input
open DesignCompiler

def study : NormalInput := {
  name := "Two-arm normal comparison"
  alternative := 0.5
  nullValue := 0
  sd := 1
  alpha := 0.05
  targetPower := 0.80
  minN := 1
  maxN := 10000
}
`,
  influence: `import DesignCompiler.Input
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
`,
  cluster: `import DesignCompiler.Input
open DesignCompiler

def study : InfluenceInput := {
  name := "Planning with independent clusters"
  alternative := 0.5
  nullValue := 0
  variance := 2
  unit := .cluster
  estimand := "Cluster-average contrast"
  estimator := "Asymptotically linear cluster estimator"
  varianceSource := "Illustrative variance per independent cluster; replace for your design"
  alpha := 0.05
  targetPower := 0.80
  minN := 10
  maxN := 10000
}
`,
};
