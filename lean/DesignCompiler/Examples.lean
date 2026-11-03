import DesignCompiler.Spec
import Mathlib.Tactic.NormNum

namespace DesignCompiler.Examples
open DesignCompiler Lean

def criterion : Criterion := {
  alpha := ⟨1/20, by norm_num, by norm_num⟩
  target := ⟨4/5, by norm_num, by norm_num⟩
  aboveAlpha := by norm_num
  minN := 1
  maxN := 10000
  nonempty := by decide
  positive := by decide
  bounded := by decide
}

/-- The motivating calculation: 63 per arm, 126 total. -/
def normal : NormalStudy := {
  name := "Two-arm normal, known SD"
  alternative := 1/2
  sd := ⟨1, by norm_num⟩
  criterion := criterion
}

/-- V = 4 per participant for an efficient ATE estimator in a balanced trial
with homoscedastic conditional variance 1 and constant treatment effect 1/2. -/
def influence : InfluenceStudy := {
  name := "ATE with supplied influence variance"
  alternative := 1/2
  variance := ⟨4, by norm_num⟩
  unit := .participant
  independent := by decide
  estimand := "ATE on the outcome scale"
  estimator := "asymptotically linear AIPW"
  varianceSource := "Planning model: balanced treatment, residual variance 1, constant effect 1/2"
  criterion := criterion
}

def stratum : Stratum := {
  weight := ⟨1, by norm_num⟩
  propensity := ⟨1/2, by norm_num, by norm_num⟩
  mean0 := 0
  mean1 := 1/2
  variance0 := 1
  variance1 := 1
  nonnegative0 := by norm_num
  nonnegative1 := by norm_num
}

def ate : AteStudy := {
  name := "ATE EIF from a planning stratum"
  strata := [stratum]
  unitMass := by norm_num [stratum]
  estimator := .tmle
  criterion := criterion
}

def getExample (name : String) : Option Study :=
  match name with
  | "normal" => some (.normal normal)
  | "influence" => some (.influence influence)
  | "ate" => some (.ate ate)
  | _ => none

end DesignCompiler.Examples
