import Lean
import Mathlib.Data.Rat.Defs

namespace DesignCompiler
open Lean

/-- Rational planning inputs retain their exact representation in the wire plan. -/
structure Positive where
  value : ℚ
  positive : 0 < value

structure Probability where
  value : ℚ
  positive : 0 < value
  belowOne : value < 1

structure Criterion where
  alpha : Probability
  target : Probability
  aboveAlpha : alpha.value < target.value
  minN : Nat
  maxN : Nat
  nonempty : minN ≤ maxN
  positive : 0 < minN
  bounded : maxN ≤ 10000000

inductive CountUnit where
  | perArm | participant | cluster
  deriving DecidableEq, Repr

def CountUnit.wire : CountUnit → String
  | .perArm => "participants_per_arm"
  | .participant => "participant"
  | .cluster => "cluster"

structure NormalStudy where
  name : String
  alternative : ℚ
  nullValue : ℚ := 0
  sd : Positive
  criterion : Criterion

structure InfluenceStudy where
  name : String
  alternative : ℚ
  nullValue : ℚ := 0
  variance : Positive
  unit : CountUnit
  independent : unit ≠ .perArm
  estimand : String
  estimator : String
  varianceSource : String
  criterion : Criterion

structure Stratum where
  weight : Positive
  propensity : Probability
  mean0 : ℚ
  mean1 : ℚ
  variance0 : ℚ
  variance1 : ℚ
  nonnegative0 : 0 ≤ variance0
  nonnegative1 : 0 ≤ variance1

inductive TargetedEstimator where
  | aipw | tmle
  deriving DecidableEq, Repr

def TargetedEstimator.wire : TargetedEstimator → String
  | .aipw => "aipw"
  | .tmle => "tmle"

structure AteStudy where
  name : String
  strata : List Stratum
  unitMass : (strata.map (fun s => s.weight.value)).sum = 1
  nullValue : ℚ := 0
  estimator : TargetedEstimator
  criterion : Criterion

inductive Study where
  | normal (study : NormalStudy)
  | influence (study : InfluenceStudy)
  | ate (study : AteStudy)

structure WaldKernel where
  effect : ℚ
  variance : ℚ
  unit : CountUnit

def lowerNormal (s : NormalStudy) : WaldKernel :=
  ⟨s.alternative - s.nullValue, 2 * s.sd.value ^ 2, .perArm⟩

def lowerInfluence (s : InfluenceStudy) : WaldKernel :=
  ⟨s.alternative - s.nullValue, s.variance.value, s.unit⟩

theorem lowerNormal_units (s : NormalStudy) : (lowerNormal s).unit = .perArm := rfl
theorem lowerInfluence_units (s : InfluenceStudy) : (lowerInfluence s).unit = s.unit := rfl
theorem lowerInfluence_variance (s : InfluenceStudy) :
    (lowerInfluence s).variance = s.variance.value := rfl

def rationalJson (x : ℚ) : Json := Json.str s!"{x.num}/{x.den}"
def searchJson (c : Criterion) : Json := Json.mkObj [("minN", toJson c.minN), ("maxN", toJson c.maxN)]
def criterionFields (c : Criterion) : List (String × Json) :=
  [("alpha", rationalJson c.alpha.value), ("targetPower", rationalJson c.target.value)]
def envelope (name method : String) (fields : List (String × Json)) : Json :=
  Json.mkObj ([("schema", toJson "ssc.design-plan"), ("version", toJson (1 : Nat)),
    ("name", toJson name), ("method", toJson method)] ++ fields)

def compile : Study → Json
  | .normal s => envelope s.name "normal_known"
      [("parameters", Json.mkObj (criterionFields s.criterion ++
        [("alternative", rationalJson s.alternative), ("nullValue", rationalJson s.nullValue),
         ("sd", rationalJson s.sd.value)])), ("search", searchJson s.criterion)]
  | .influence s => envelope s.name "influence_wald"
      [("parameters", Json.mkObj (criterionFields s.criterion ++
        [("alternative", rationalJson s.alternative), ("nullValue", rationalJson s.nullValue),
         ("variance", rationalJson (lowerInfluence s).variance),
         ("independentUnit", toJson (lowerInfluence s).unit.wire),
         ("estimand", toJson s.estimand), ("estimator", toJson s.estimator),
         ("varianceSource", toJson s.varianceSource)])), ("search", searchJson s.criterion)]
  | .ate s => envelope s.name "ate_eif"
      [("parameters", Json.mkObj (criterionFields s.criterion ++
        [("nullValue", rationalJson s.nullValue), ("estimator", toJson s.estimator.wire),
         ("strata", toJson (s.strata.map fun x => Json.mkObj
           [("weight", rationalJson x.weight.value), ("propensity", rationalJson x.propensity.value),
            ("mean0", rationalJson x.mean0), ("mean1", rationalJson x.mean1),
            ("variance0", rationalJson x.variance0), ("variance1", rationalJson x.variance1)]))])),
       ("search", searchJson s.criterion)]

end DesignCompiler
