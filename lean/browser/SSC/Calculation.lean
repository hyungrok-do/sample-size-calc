import SSC.Probability
import Lean.Data.Json

namespace SSC.Calculation
open Lean SSC.Numeric SSC.Probability

/-- One specification for both direct calculation and simulation. The design
type can contain n, allocations, clusters, effects, follow-up, or other parameters.
The outcome can be a scalar statistic or a complete user-defined dataset. -/
structure Specification (Design Outcome : Type) where
  null : Design → Law Outcome
  alternative : Design → Law Outcome
  reject : Design → Event Outcome

inductive Route where
  | automatic | exact | numerical | simulation
  deriving Repr, BEq

structure Options where
  route : Route := .automatic
  exactBudget : Nat := 10000
  replicates : Nat := 20000
  seed : UInt64 := 20261104
  confidence : Float := 0.95
  maxDraws : Nat := 10000000
  numericalGuard : Float := 1e-7
  deriving Repr

def Options.validate (o : Options) : Result Unit := do
  if o.replicates < 100 || o.replicates > 1000000 then throw "Use 100 through 1000000 simulation replicates."
  if !o.confidence.isFinite || o.confidence <= 0 || o.confidence >= 1 then throw "Confidence must be between zero and one."
  if !o.numericalGuard.isFinite || o.numericalGuard < 0 || o.numericalGuard > 0.01 then throw "Numerical decision guard must be between 0 and 0.01."
  if o.exactBudget > 1000000 || o.maxDraws == 0 || o.maxDraws > 100000000 then throw "Invalid calculation work budget."

structure Estimate where
  method : String
  value : Float
  exact : Option Rat := none
  lower : Option Float := none
  upper : Option Float := none
  successes : Option Nat := none
  replicates : Nat := 0
  seed : Option String := none
  confidence : Option Float := none
  reason : String := ""
  deriving Repr

private def simulation {α : Type} (law : Law α) (event : Event α)
    (o : Options) (failureProbability : Float) (reason : String) : Draw Estimate := do
  event.validate
  let initial ← get
  let mut successes := 0
  for _ in [:o.replicates] do
    let x ← sample law
    if ← event.accepts x then successes := successes + 1
  let value := nf successes / nf o.replicates
  -- Fixed simulation budget: no optional stopping. A union bound across all
  -- null/alternative evaluations in a search gives simultaneous coverage.
  let radius := ((2/failureProbability).log / (2*nf o.replicates)).sqrt
  return {
    method := "simulation"
    value := value
    lower := some (max 0 (value-radius))
    upper := some (min 1 (value+radius))
    successes := some successes
    replicates := o.replicates
    seed := some (toString initial.seed)
    confidence := some (1-failureProbability)
    reason := reason }

private def evaluateProbability {α : Type} (law : Law α) (event : Event α)
    (o : Options) (failureProbability : Float) : Draw Estimate := do
  let mut reason := "Simulation explicitly requested for this specification."
  if o.route == .automatic || o.route == .exact then
    match exactProbability law event o.exactBudget with
    | .ok p => return { method := "exact", value := probabilityFloat p, exact := some p }
    | .error (.invalid message) => throw message
    | .error (.unavailable message) =>
      if o.route == .exact then throw message
      reason := message
  if o.route == .automatic || o.route == .numerical then
    match numericalProbability law event with
    | .ok p => return { method := "numerical", value := p, reason := "Floating-point CDF/integration; no certified numerical error bound." }
    | .error (.invalid message) => throw message
    | .error (.unavailable message) =>
      if o.route == .numerical then throw message
      reason := reason ++ " " ++ message
  simulation law event o failureProbability reason

structure Evaluation (Design : Type) where
  design : Design
  typeI : Estimate
  power : Estimate
  deriving Repr

private def evaluateWith {D X : Type} (spec : Specification D X) (design : D)
    (o : Options) (failureProbability : Float) : Draw (Evaluation D) := do
  let event := spec.reject design
  event.validate
  let typeI ← evaluateProbability (spec.null design) event o failureProbability
  let power ← evaluateProbability (spec.alternative design) event o failureProbability
  return { design, typeI, power }

/-- Evaluate power and type-I error at any supplied design, under one shared spec. -/
def evaluate {D X : Type} (spec : Specification D X) (design : D)
    (o : Options := {}) : Result (Evaluation D) := do
  o.validate
  let (answer, _) ← (evaluateWith spec design o ((1-o.confidence)/2)).run ⟨o.seed, o.maxDraws⟩
  return answer

structure Requirements where
  power : Rat := 4/5
  /-- If omitted, type-I error is reported but not used as a search constraint. -/
  maxTypeI : Option Rat := none
  deriving Repr

inductive Decision where
  | meets | fails | unresolved
  deriving Repr, BEq

private def atLeast (e : Estimate) (target : Rat) (guard : Float) : Decision :=
  match e.exact with
  | some p => if p >= target then .meets else .fails
  | none =>
    let t := probabilityFloat target
    match e.lower, e.upper with
    | some lo, some hi => if lo >= t then .meets else if hi < t then .fails else .unresolved
    | _, _ => if e.value >= t+guard then .meets else if e.value < t-guard then .fails else .unresolved

private def atMost (e : Estimate) (target : Rat) (guard : Float) : Decision :=
  match e.exact with
  | some p => if p <= target then .meets else .fails
  | none =>
    let t := probabilityFloat target
    match e.lower, e.upper with
    | some lo, some hi => if hi <= t then .meets else if lo > t then .fails else .unresolved
    | _, _ => if e.value <= t-guard then .meets else if e.value > t+guard then .fails else .unresolved

def decideRequirements {D : Type} (e : Evaluation D) (r : Requirements) (guard : Float) : Decision :=
  let power := atLeast e.power r.power guard
  let size := r.maxTypeI.map (fun a => atMost e.typeI a guard) |>.getD .meets
  if power == .fails || size == .fails then .fails
  else if power == .meets && size == .meets then .meets else .unresolved

structure SearchRow (D : Type) where
  evaluation : Evaluation D
  decision : Decision
  deriving Repr

structure Search (D : Type) where
  rows : Array (SearchRow D)
  selected : Option Nat
  earlierUnresolved : Array Nat
  requirements : Requirements
  options : Options
  deriving Repr

/-- Search an ordered candidate set without assuming monotonic power. Search
does not claim a global minimum outside that set, nor resolve uncertain rows. -/
def search {D X : Type} (spec : Specification D X) (candidates : Array D)
    (requirements : Requirements := {}) (o : Options := {}) : Result (Search D) := do
  o.validate
  if candidates.isEmpty || candidates.size > 10000 then throw "Supply 1 through 10000 ordered candidates."
  if requirements.power <= 0 || requirements.power >= 1 then throw "Target power must lie strictly between zero and one."
  if requirements.maxTypeI.any (fun a => a <= 0 || a >= 1) then throw "Type-I error limit must lie strictly between zero and one."
  let failureProbability := (1-o.confidence) / nf (2*candidates.size)
  let calculation : Draw (Search D) := do
    let mut rows := #[]
    let mut uncertain := #[]
    let mut selected := none
    for design in candidates do
      let i := rows.size
      let evaluation ← evaluateWith spec design o failureProbability
      let decision := decideRequirements evaluation requirements o.numericalGuard
      rows := rows.push { evaluation, decision }
      if decision == .unresolved then uncertain := uncertain.push i
      if decision == .meets then
        selected := some i
        break
    return { rows, selected, earlierUnresolved := uncertain, requirements, options := o }
  let (answer, _) ← calculation.run ⟨o.seed, o.maxDraws⟩
  return answer

/-- Preserve every binary64 bit. The bundled Lean Float.toString / ToJson rounds
to six decimal places, which is unsuitable for probability reports and tests. -/
def floatJson (x : Float) : Json := Id.run do
  if !x.isFinite then return .null
  if x == 0 then return .num 0
  let bits := x.toBits
  let exponent := ((bits >>> 52) &&& 2047).toNat
  let fraction := (bits &&& 0x000fffffffffffff).toNat
  let mantissa := if exponent == 0 then fraction else fraction + 4503599627370496
  let power : Int := if exponent == 0 then -1074 else (exponent : Int)-1075
  let sign : Int := if bits >>> 63 == 0 then 1 else -1
  if power >= 0 then return .num ⟨sign * (mantissa * 2^power.toNat : Nat), 0⟩
  let places := power.natAbs
  return .num ⟨sign * (mantissa * 5^places : Nat), places⟩

private def optionalFloat (x : Option Float) : Json := x.map floatJson |>.getD .null
private def Estimate.json (e : Estimate) : Json := Json.mkObj [
  ("method", toJson e.method), ("value", floatJson e.value),
  ("exact", toJson (e.exact.map toString)), ("lower", optionalFloat e.lower),
  ("upper", optionalFloat e.upper), ("successes", toJson e.successes),
  ("replicates", toJson e.replicates), ("seed", toJson e.seed),
  ("confidence", optionalFloat e.confidence), ("reason", toJson e.reason)]

private def Evaluation.json {D : Type} [ToJson D] (e : Evaluation D) : Json := Json.mkObj [
  ("design", toJson e.design), ("typeI", e.typeI.json), ("power", e.power.json)]

private def scope := "Exact denotes rational evaluation of the supplied finite program. Numerical results use floating-point approximations. Simulation intervals use fixed-budget Hoeffding bounds, assuming independent draws from the specified law; the pseudorandom generator and floating-point samplers are not kernel-certified. Theorems are reported separately."

/-- Standard Lean IO output. No custom command syntax or compiler dispatcher. -/
def emitEvaluation {D : Type} [ToJson D] (answer : Result (Evaluation D)) : IO Unit := do
  match answer with
  | .error message => throw (IO.userError message)
  | .ok e => IO.println ("SSC_CALCULATION " ++ (Json.mkObj [
      ("schema", toJson "ssc.calculation.v1"), ("kind", toJson "evaluation"),
      ("evaluation", e.json), ("scope", toJson scope)]).compress)

def emitSearch {D : Type} [ToJson D] (answer : Result (Search D)) : IO Unit := do
  match answer with
  | .error message => throw (IO.userError message)
  | .ok s =>
    let rows := s.rows.map fun row => Json.mkObj [
      ("evaluation", row.evaluation.json), ("decision", toJson (match row.decision with
        | .meets => "meets" | .fails => "fails" | .unresolved => "unresolved"))]
    IO.println ("SSC_CALCULATION " ++ (Json.mkObj [
      ("schema", toJson "ssc.calculation.v1"), ("kind", toJson "search"),
      ("rows", toJson rows), ("selected", toJson s.selected),
      ("earlierUnresolved", toJson s.earlierUnresolved),
      ("targetPower", toJson (toString s.requirements.power)),
      ("maxTypeI", toJson (s.requirements.maxTypeI.map toString)),
      ("simultaneousConfidence", floatJson s.options.confidence),
      ("scope", toJson scope)]).compress)

end SSC.Calculation
