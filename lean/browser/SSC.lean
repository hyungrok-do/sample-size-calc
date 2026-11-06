import Lean
import Lean.Replay
import Lean.Elab.Frontend
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

/-- Read actual elaborator trees. User stdout is never a proof trace. -/
private def collectTrace : CommandElabM (Array Json × Array String) := do
  let trees := (← getInfoState).substituteLazy.get.trees
  let steps ← IO.mkRef (#[] : Array Json)
  let limits ← IO.mkRef (#[] : Array String)
  let counts ← IO.mkRef ({} : Std.HashMap Name Nat)
  for tree in trees do
    discard <| tree.visitM (m := CommandElabM) (α := Unit)
      (postNode := fun ctx info _ _ => do
        let .ofTacticInfo tactic := info | return
        let some decl := ctx.parentDecl? | return
        let some range := tactic.stx.getRange? | return
        -- Synthetic automation internals are available in the proof term. Only
        -- original source tactics are useful as a navigable derivation trace.
        unless tactic.stx.getHeadInfo matches .original .. do return
        let count := (← counts.get).getD decl 0
        if count >= 256 then
          if count == 256 then
            limits.modify (·.push (displayName decl))
            counts.modify (·.insert decl (count + 1))
          return
        let before ← ({ ctx with mctx := tactic.mctxBefore }).ppGoals tactic.goalsBefore
        let after ← ({ ctx with mctx := tactic.mctxAfter }).ppGoals tactic.goalsAfter
        if before.pretty == after.pretty then return
        counts.modify (·.insert decl (count + 1))
        let step := Json.mkObj [
          ("declaration", displayName decl),
          ("line", toJson (ctx.fileMap.toPosition range.start).line),
          ("column", toJson (ctx.fileMap.toPosition range.start).column),
          ("endLine", toJson (ctx.fileMap.toPosition range.stop).line),
          ("tactic", toJson (String.Pos.Raw.extract ctx.fileMap.source range.start range.stop)),
          ("before", toJson before.pretty), ("after", toJson after.pretty)]
        steps.modify (·.push step))
  return (← steps.get, ← limits.get)

/-- Audit the completed environment, independently of the input's command syntax. -/
private def checkedReport : CommandElabM Json := do
  let started ← IO.monoMsNow
  let env ← getEnv
  -- Force kernel checking, including any asynchronous declarations, before
  -- examining results. Never certify an elaborator-only environment entry.
  let checked := env.checked.get
  let declarations := checked.constants.foldStage2 (fun xs n ci => xs.push (n, ci)) #[]
  -- A source module can disable elaborator kernel checks, even temporarily
  -- inside a custom tactic. Start with only the imported constants and replay ALL
  -- local declarations with the kernel's unconditional addDeclCore. Replay
  -- orders dependencies and reconstructs inductives/recursors, so unchecked
  -- constants cannot validate themselves or each other.
  let some imported := env.importEnv?
    | throwError "Cannot independently check declarations without the import environment."
  let constants := declarations.foldl (fun m (n, ci) => m.insert n ci)
    ({} : Std.HashMap Name ConstantInfo)
  discard <| imported.replay constants
  let replayMs := (← IO.monoMsNow) - started
  let mut axiomMs := 0
  let mut renderMs := 0
  let mut results : Array Json := #[]
  for (n, ci) in declarations do
    let auditStarted ← IO.monoMsNow
    let axioms ← collectAxioms n
    let unapproved := axioms.filter (!approvedAxiom ·)
    unless unapproved.isEmpty do
      throwError "Unproved dependency in {n}: {unapproved.toList}. State hypotheses as theorem parameters and discharge the goal. Axioms, sorry, and native_decide are not proof certificates here."
    axiomMs := axiomMs + ((← IO.monoMsNow) - auditStarted)
    let some ranges ← findDeclarationRangesCore? n | continue
    if (privateToUserName? n |>.getD n).isInternal then continue
    let kind := match ci with
      | .thmInfo _ => "theorem"
      | .defnInfo v => if v.safety == .safe then "definition" else "unchecked definition"
      | .opaqueInfo v => if v.isUnsafe then "unchecked definition" else "opaque definition"
      | .inductInfo _ => "type"
      | _ => "generated"
    if kind == "generated" then continue
    let renderStarted ← IO.monoMsNow
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
    renderMs := renderMs + ((← IO.monoMsNow) - renderStarted)
  let traceStarted ← IO.monoMsNow
  let (steps, limits) ← collectTrace
  return Json.mkObj [("declarations", toJson results), ("steps", toJson steps),
    ("traceLimits", toJson limits),
    ("performance", Json.mkObj [("kernelReplayMs", toJson replayMs),
      ("axiomAuditMs", toJson axiomMs), ("declarationRenderMs", toJson renderMs),
      ("traceRenderMs", toJson ((← IO.monoMsNow) - traceStarted))])]

/-- Standalone reporting helper. The browser does not trust calls in user source. -/
def exportReport (nonce : String) : CommandElabM Unit := do
  IO.println s!"SSC_PROOF_REPORT {nonce} {(← checkedReport).compress}"

private structure CachedImport where
  imports : Array Import
  isModule : Bool
  env : Environment

/-- A runtime interface for ordinary parameterized Lean programs. Evaluating
this function is compiled execution, not a new kernel proof of its output. -/
abbrev PreparedFunction := Json → Except String Json

private unsafe def prepareFunction (entry : String) : CommandElabM PreparedFunction := do
  let name := (entry.splitOn ".").foldl Name.str .anonymous
  let info ← getConstInfo name
  match info with
  | .defnInfo value =>
    unless value.safety == .safe do
      throwError "Prepared entry must be a safe definition."
  | _ => throwError "Prepared entry must be a definition."
  unless info.levelParams.isEmpty do
    throwError "Prepared entry must not have universe parameters."
  liftTermElabM do
    unless ← isDefEq info.type (mkConst ``PreparedFunction) do
      throwError "Prepared entry must have type Lean.Json → Except String Lean.Json."
    -- evalConst's result type is a host assertion, so the exact function type
    -- above is mandatory. This runs only after independent replay and audit.
    match (← getEnv).evalConst PreparedFunction (← getOptions) name with
    | .ok fn => pure fn
    | .error message => throwError "{message}"

private def evaluatePrepared (prepared : IO.Ref (Option PreparedFunction))
    (input : String) : IO String := do
  -- Empty input is the host's reset operation; valid JSON is never empty.
  if input.isEmpty then
    prepared.set none
    return (Json.mkObj [("ok", toJson true), ("value", Json.null)]).compress
  let result := do
    let some fn ← prepared.get | return Except.error "No prepared Lean function."
    return (Json.parse input).bind fn
  return (match ← result with
    | .ok value => Json.mkObj [("ok", toJson true), ("value", value)]
    | .error message => Json.mkObj [("ok", toJson false), ("error", toJson message)]).compress

/-- Prepare the ordinary Lean frontend one command at a time. The host executes
returned native frontend actions directly, so the kernel does not run beneath
an interpreted driver or a nested command evaluation stack. Diagnostics and info
trees are retained before each top-level command resets them. -/
private unsafe def prepareSource (imports : IO.Ref (Array CachedImport))
    (prepared : IO.Ref (Option PreparedFunction)) (source nonce entry : String) :
    IO (IO (Option (IO Unit))) := do
  prepared.set none
  let started ← IO.monoMsNow
  let inputCtx := Parser.mkInputContext source "/workspace/source.lean"
  let (header, parserState, headerMessages) ← Parser.parseHeader inputCtx
  let opts := Elab.async.set ({} : Options) false
  let key := Elab.HeaderSyntax.imports header
  let isModule := Elab.HeaderSyntax.isModule header
  let (env, headerMessages) ←
    match (← imports.get).find? (fun entry => entry.imports == key && entry.isModule == isModule) with
    | some entry => pure (entry.env.setMainModule (Name.mkSimple "SSCInput"), headerMessages)
    | none => do
      unsafe enableInitializersExecution
      let (env, log) ← Elab.processHeader header opts headerMessages inputCtx
        (leakEnv := true) (mainModule := Name.mkSimple "SSCInput")
      -- Initializers retain references to imported regions. Keep each clean
      -- import environment once, rather than leaking a new copy per request.
      if log.toList.isEmpty then
        imports.modify (·.push { imports := key, isModule, env })
      pure (env, log)
  let headerMs := (← IO.monoMsNow) - started
  let context : Elab.Frontend.Context := { inputCtx }
  let state ← IO.mkRef ({
    commandState := Command.mkState env {} opts
    parserState, cmdPos := parserState.pos } : Elab.Frontend.State)
  let messages ← IO.mkRef headerMessages
  let trees ← IO.mkRef ({} : PersistentArray InfoTree)
  let pending ← IO.mkRef (none : Option (Nat × Nat × String))
  let commands ← IO.mkRef (#[] : Array Json)
  let parseMs ← IO.mkRef 0
  let done ← IO.mkRef false
  return do
    if let some (commandStarted, line, kind) ← pending.get then
      let commandMs := (← IO.monoMsNow) - commandStarted
      commands.modify (·.push (Json.mkObj [("line", toJson line),
        ("kind", toJson kind), ("ms", toJson commandMs)]))
      let current ← state.get
      messages.modify (· ++ current.commandState.messages)
      trees.modify (·.append current.commandState.infoState.substituteLazy.get.trees)
      pending.set none
    if ← done.get then
      let current ← state.get
      let allMessages ← messages.get
      let allTrees ← trees.get
      state.set { current with commandState := {
        current.commandState with messages := allMessages
                                  infoState.trees := allTrees } }
      let result ← (show IO Json from do
        try
          let log ← messages.get
          let diagnostics ← liftM (m := BaseIO) <| log.toList.toArray.mapM Message.toJson
          if log.hasErrors then
            return Json.mkObj [("ok", toJson false), ("diagnostics", toJson diagnostics)]
          let report ← Elab.Frontend.runCommandElabM checkedReport context state
          if !entry.isEmpty then
            let fn ← Elab.Frontend.runCommandElabM (prepareFunction entry) context state
            prepared.set (some fn)
          return Json.mkObj [("ok", toJson true), ("report", report),
            ("diagnostics", toJson diagnostics),
            ("performance", Json.mkObj [("headerMs", toJson headerMs),
              ("parseMs", toJson (← parseMs.get)), ("commands", toJson (← commands.get))])]
        catch error =>
          return Json.mkObj [("ok", toJson false), ("error", toJson error.toString)])
      IO.println s!"SSC_COMPILE_RESULT {nonce} {result.compress}"
      return none
    let current ← state.get
    let parseStarted ← IO.monoMsNow
    let scope := current.commandState.scopes.head!
    let parserContext : Parser.ParserModuleContext := {
      env := current.commandState.env, options := scope.opts
      currNamespace := scope.currNamespace, openDecls := scope.openDecls }
    let (command, parserState, parseMessages) := Parser.parseCommand inputCtx
      parserContext current.parserState {}
    let elapsed := (← IO.monoMsNow) - parseStarted
    parseMs.modify (· + elapsed)
    messages.modify (· ++ parseMessages)
    state.set { current with parserState, cmdPos := current.parserState.pos }
    done.set (Parser.isTerminalCommand command)
    let line := inputCtx.fileMap.toPosition (command.getPos?.getD current.parserState.pos)
    pending.set (some ((← IO.monoMsNow), line.line, command.getKind.toString))
    return some (Elab.Frontend.elabCommandAtFrontend command context state)

/-- Persistent host entry points share the clean import cache and one replaceable
runtime function. User functions are not marked persistent. The host invokes
the frontend directly without nesting compilation inside command evaluation. -/
syntax (name := browserDriver) "#ssc_browser_driver" : command

@[command_elab browserDriver] private unsafe def elabBrowserDriver : CommandElab
  | `(#ssc_browser_driver) => do
      let some imported := (← getEnv).importEnv?
        | throwError "Missing trusted compiler import environment."
      let entry : CachedImport := { imports := imported.header.imports, isModule := false, env := imported }
      let imports ← IO.mkRef #[entry]
      let prepared ← IO.mkRef (none : Option PreparedFunction)
      let driver ← liftM (m := BaseIO) <| Runtime.markPersistent (prepareSource imports prepared)
      let evaluator ← liftM (m := BaseIO) <| Runtime.markPersistent (evaluatePrepared prepared)
      IO.println s!"SSC_BROWSER_DRIVER {(ptrAddrUnsafe driver).toNat}"
      IO.println s!"SSC_BROWSER_EVALUATOR {(ptrAddrUnsafe evaluator).toNat}"
  | _ => throwUnsupportedSyntax

end SSC.Report
