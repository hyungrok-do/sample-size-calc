import DesignCompiler.Spec

namespace DesignCompiler
open Lean

/-- Editor input. Semantic constraints are validated by the numerical backend.
The proof-carrying `NormalStudy` type remains available for library development. -/
structure NormalInput where
  name : String := "Normal study"
  alternative : ℚ := 1/2
  nullValue : ℚ := 0
  sd : ℚ := 1
  alpha : ℚ := 1/20
  targetPower : ℚ := 4/5
  minN : Nat := 1
  maxN : Nat := 1000000

structure InfluenceInput where
  name : String := "Influence-function study"
  alternative : ℚ := 1/2
  nullValue : ℚ := 0
  variance : ℚ := 4
  unit : CountUnit := .participant
  estimand : String := "ATE"
  estimator : String := "AIPW"
  varianceSource : String := "Supplied per-unit planning variance"
  alpha : ℚ := 1/20
  targetPower : ℚ := 4/5
  minN : Nat := 1
  maxN : Nat := 1000000

class CompileInput (α : Type) where
  toPlan : α → Json

instance : CompileInput NormalInput where
  toPlan s := envelope s.name "normal_known"
    [("parameters", Json.mkObj [("alternative", rationalJson s.alternative),
      ("nullValue", rationalJson s.nullValue), ("sd", rationalJson s.sd),
      ("alpha", rationalJson s.alpha), ("targetPower", rationalJson s.targetPower)]),
     ("search", Json.mkObj [("minN", toJson s.minN), ("maxN", toJson s.maxN)])]

instance : CompileInput InfluenceInput where
  toPlan s := envelope s.name "influence_wald"
    [("parameters", Json.mkObj [("alternative", rationalJson s.alternative),
      ("nullValue", rationalJson s.nullValue), ("variance", rationalJson s.variance),
      ("independentUnit", toJson s.unit.wire), ("estimand", toJson s.estimand),
      ("estimator", toJson s.estimator), ("varianceSource", toJson s.varianceSource),
      ("alpha", rationalJson s.alpha), ("targetPower", rationalJson s.targetPower)]),
     ("search", Json.mkObj [("minN", toJson s.minN), ("maxN", toJson s.maxN)])]

def compileInput {α : Type} [CompileInput α] (s : α) : Json := CompileInput.toPlan s

end DesignCompiler
