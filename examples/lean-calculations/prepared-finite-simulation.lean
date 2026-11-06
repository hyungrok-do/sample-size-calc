import SSC
open Lean SSC.Numeric SSC.Probability SSC.Calculation

-- Runtime parameters choose an observation count, Bernoulli probability,
-- rejection threshold, seed and replicate count. Every draw is retained.
def sscEvaluate (input : Json) : Except String Json := do
  let n ← input.getObjValAs? Nat "n"
  let numerator ← input.getObjValAs? Nat "numerator"
  let denominator ← input.getObjValAs? Nat "denominator"
  let cutoff ← input.getObjValAs? Nat "cutoff"
  let replicates ← input.getObjValAs? Nat "replicates"
  let seed ← input.getObjValAs? Nat "seed"
  if n == 0 || n > 50 then throw "n must be between 1 and 50."
  if denominator == 0 || numerator > denominator then throw "Invalid probability."
  if replicates > 500 then throw "This example permits at most 500 replicates."
  let law := (Law.bernoulli ((numerator : Rat) / denominator)).iid n
  let model : Specification Nat (List Bool) := {
    null := fun _ => law
    alternative := fun _ => law
    reject := fun _ => .ofBool (fun xs => (xs.filter id).length >= cutoff) }
  let answer ← evaluate model n { route := .simulation, replicates, seed := seed.toUInt64 }
  return Json.mkObj [
    ("nullSuccesses", toJson answer.typeI.successes),
    ("alternativeSuccesses", toJson answer.power.successes),
    ("nullSeed", toJson answer.typeI.seed),
    ("alternativeSeed", toJson answer.power.seed),
    ("replicates", toJson answer.power.replicates),
    ("lower", answer.power.lower.map floatJson |>.getD .null),
    ("upper", answer.power.upper.map floatJson |>.getD .null),
    ("method", toJson "simulation"),
    ("concreteResultKernelChecked", toJson false)]
