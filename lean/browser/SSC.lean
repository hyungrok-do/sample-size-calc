import Lean
import SSC.Calculation
import SSC.Derivation

/-!
Browser reporting for ordinary Lean modules. Imported SSC libraries provide
general probability programs, certified finite derivation, numerical mathematics and simulation. Statistical
models live in ordinary source modules. Proof search is performed by Lean tactics
in the user's source; evaluation output is distinct from kernel-checked theorems.

The observer records actual elaborator goals. The exporter reads declarations
from the checked environment, including their full transitive axiom dependencies.
Neither component creates a mathematical claim.
-/

open Lean Elab Command Meta

namespace SSC

/-- A bounded search returns a least witness with a proof, or a proof that every
candidate in the domain fails. It assumes neither continuity nor monotonicity. -/
inductive SearchResult (P : Nat → Prop) (bound : Nat) where
  | found (value : Nat) (within : value < bound) (holds : P value)
      (minimal : ∀ m, m < value → ¬ P m)
  | exhausted (proof : ∀ m, m < bound → ¬ P m)

/-- Generic proof-carrying finite search over any decidable predicate. -/
def leastBelow (P : Nat → Prop) [DecidablePred P] :
    (bound : Nat) → SearchResult P bound
  | 0 => .exhausted (by intro m h; omega)
  | b + 1 =>
    match leastBelow P b with
    | .found n hn hp hmin => .found n (by omega) hp hmin
    | .exhausted h =>
      if hp : P b then .found b (by omega) hp h
      else .exhausted (by
        intro m hm
        by_cases hb : m < b
        · exact h m hb
        · have heq : m = b := by omega
          simpa [heq] using hp)

def SearchResult.value? {P : Nat → Prop} {b : Nat} : SearchResult P b → Option Nat
  | .found n .. => some n
  | .exhausted _ => none

/-- The proposition attached to each computed answer, including exhausted bounds. -/
def SearchResult.Spec (P : Nat → Prop) (b : Nat) : Option Nat → Prop
  | some n => n < b ∧ P n ∧ ∀ m, m < n → ¬ P m
  | none => ∀ m, m < b → ¬ P m

theorem SearchResult.sound {P : Nat → Prop} {b : Nat} (r : SearchResult P b) :
    Spec P b r.value? := by
  cases r with
  | found n hn hp hmin => exact ⟨hn, hp, hmin⟩
  | exhausted h => exact h

end SSC

namespace SSC.Report

private def strings (xs : Array String) : Json := toJson xs

private def approvedAxiom (n : Name) : Bool :=
  n == ``propext || n == ``Classical.choice || n == ``Quot.sound

private def displayName (n : Name) : String :=
  (privateToUserName? n |>.getD n).toString

/-- Record source tactics and the goals they actually transformed. -/
private def observe (_ : Syntax) : CommandElabM Unit := do
  if (← get).messages.hasErrors then return
  let trees := (← getInfoState).substituteLazy.get.trees
  let count ← IO.mkRef (0 : Nat)
  for tree in trees do
    discard <| tree.visitM (m := CommandElabM) (α := Unit)
      (postNode := fun ctx info _ _ => do
        let .ofTacticInfo tactic := info | return
        let some decl := ctx.parentDecl? | return
        let some range := tactic.stx.getRange? | return
        -- Synthetic automation internals are available in the proof term. Only
        -- original source tactics are useful as a navigable derivation trace.
        unless tactic.stx.getHeadInfo matches .original .. do return
        if (← count.get) >= 256 then
          if (← count.get) == 256 then
            IO.println s!"SSC_PROOF_TRACE_LIMIT {displayName decl}"
            count.modify (· + 1)
          return
        let before ← ({ ctx with mctx := tactic.mctxBefore }).ppGoals tactic.goalsBefore
        let after ← ({ ctx with mctx := tactic.mctxAfter }).ppGoals tactic.goalsAfter
        if before.pretty == after.pretty then return
        count.modify (· + 1)
        let step := Json.mkObj [
          ("declaration", displayName decl),
          ("line", toJson (ctx.fileMap.toPosition range.start).line),
          ("column", toJson (ctx.fileMap.toPosition range.start).column),
          ("endLine", toJson (ctx.fileMap.toPosition range.stop).line),
          ("tactic", toJson (String.Pos.Raw.extract ctx.fileMap.source range.start range.stop)),
          ("before", toJson before.pretty), ("after", toJson after.pretty)]
        IO.println s!"SSC_PROOF_STEP {step.compress}")

initialize Lean.Elab.Command.addLinter { run := observe }

/-- Called by the browser using Lean's standard `run_cmd` command. -/
def exportReport (nonce : String) : CommandElabM Unit := do
  let env ← getEnv
  -- Force kernel checking, including any asynchronous declarations, before
  -- examining results. Never certify an elaborator-only environment entry.
  let checked := env.checked.get
  let declarations := checked.constants.foldStage2 (fun xs n ci => xs.push (n, ci)) #[]
  let mut results : Array Json := #[]
  for (n, ci) in declarations do
    let axioms ← collectAxioms n
    let unapproved := axioms.filter (!approvedAxiom ·)
    unless unapproved.isEmpty do
      throwError "Unproved dependency in {n}: {unapproved.toList}. State hypotheses as theorem parameters and discharge the goal. Axioms, sorry, and native_decide are not proof certificates here."
    let some ranges ← findDeclarationRangesCore? n | continue
    if (privateToUserName? n |>.getD n).isInternal then continue
    let kind := match ci with
      | .thmInfo _ => "theorem"
      | .defnInfo v => if v.safety == .safe then "definition" else "unchecked definition"
      | .opaqueInfo v => if v.isUnsafe then "unchecked definition" else "opaque definition"
      | .inductInfo _ => "type"
      | _ => "generated"
    if kind == "generated" then continue
    let (type, body, parameters) ← liftTermElabM do
      let type := (← ppExpr ci.type).pretty
      let body ← withOptions (fun o => o.setBool `pp.proofs true) do
        match ci.value? (allowOpaque := true) with
        | some v => return (← ppExpr v).pretty
        | none => return ""
      let parameters ← forallTelescope ci.type fun xs _ => xs.mapM fun x => do
        let d ← x.fvarId!.getDecl
        let t ← inferType x
        return Json.mkObj [("name", toJson d.userName.toString),
          ("type", toJson (← ppExpr t).pretty), ("hypothesis", toJson (← isProp t))]
      return (type, body, parameters)
    let used := ci.type.getUsedConstants ++ ((ci.value? (allowOpaque := true)).map Expr.getUsedConstants |>.getD #[])
    let used := used.toList.mergeSort Name.lt |>.eraseDups |>.toArray
    results := results.push <| Json.mkObj [
      ("name", toJson (displayName n)), ("kind", toJson kind),
      ("line", toJson ranges.selectionRange.pos.line),
      ("column", toJson ranges.selectionRange.pos.column),
      ("statement", toJson type), ("body", toJson body),
      ("parameters", toJson parameters),
      ("dependencies", strings (used.map displayName)),
      ("axioms", strings (axioms.map Name.toString))]
  IO.println s!"SSC_PROOF_REPORT {nonce} {(Json.mkObj [("declarations", toJson results)]).compress}"

end SSC.Report
