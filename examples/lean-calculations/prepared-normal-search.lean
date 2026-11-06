import SSC
open Lean SSC.Numeric SSC.Probability SSC.Calculation

-- Ordinary Lean entry point for the optional prepare/evaluate browser API.
-- The JSON boundary validates data; it never interpolates parameters into source.
-- Inputs use integers to preserve their meaning across the JSON boundary.
def sscEvaluate (input : Json) : Except String Json := do
  let effectMilli ← input.getObjValAs? Nat "effectMilli"
  let alphaMillion ← input.getObjValAs? Nat "alphaMillion"
  let targetMilli ← input.getObjValAs? Nat "targetMilli"
  let maxN ← input.getObjValAs? Nat "maxN"
  if effectMilli == 0 || effectMilli > 10000 then
    throw "effectMilli must be between 1 and 10000."
  if alphaMillion == 0 || alphaMillion >= 1000000 then
    throw "alphaMillion must be between 1 and 999999."
  if targetMilli == 0 || targetMilli >= 1000 then
    throw "targetMilli must be between 1 and 999."
  if maxN < 2 || maxN > 1000 then
    throw "maxN must be between 2 and 1000."
  let critical ← normalQuantile (1 - nf alphaMillion / 2000000)
  let effect := nf effectMilli / 1000
  let model : Specification Nat Float := {
    null := fun _ => Law.normal 0 1
    alternative := fun n => Law.normal (effect * (nf n / 2).sqrt) 1
    reject := fun _ => .outside (-critical) critical }
  let requirements : Requirements := { power := (targetMilli : Rat) / 1000 }
  let answer ← search model ((List.range (maxN-1)).map (·+2)).toArray
    requirements { route := .numerical }
  let chosen : Option (Evaluation Nat) := answer.selected.bind
    (fun i => (answer.rows[i]?).map (fun row => row.evaluation))
  let selected := chosen.map (fun e => e.design)
  return Json.mkObj [
    ("selected", toJson selected),
    ("evaluatedCandidates", toJson answer.rows.size),
    ("power", chosen.map (fun e => floatJson e.power.value) |>.getD .null),
    ("typeI", chosen.map (fun e => floatJson e.typeI.value) |>.getD .null),
    ("earlierUnresolved", toJson answer.earlierUnresolved),
    ("method", toJson "numerical"),
    ("concreteResultKernelChecked", toJson false)]

-- Example input:
-- {"effectMilli":500,"alphaMillion":50000,"targetMilli":800,"maxN":200}
-- This searches n=2,...,200 per group, returning 63 for that input.
-- Probability calculations use the existing floating-point normal CDF.
-- The selected n is an executable result, not a kernel-certified minimum.
