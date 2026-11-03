import DesignCompiler.Solve.Rules
import DesignCompiler.Solve.Lower
import Lean.Util.CollectAxioms

open Lean Meta Elab Command

namespace DesignCompiler.Equation

structure DerivationStep where
  rule : String
  before : String
  after : String
  premises : Array String := #[]
  deriving ToJson

structure EquationResult where
  name : String
  status : String
  input : String
  solution : String
  theoremType : String
  proof : String
  axioms : Array String
  steps : Array DerivationStep
  explored : Nat
  lowering : Lowering
  deriving ToJson

initialize resultExtension : SimplePersistentEnvExtension EquationResult (Array EquationResult) ←
  registerSimplePersistentEnvExtension {
    addEntryFn := Array.push
    addImportedFn := fun _ => #[]
  }

private def render (e : Expr) : MetaM String := return (← ppExpr e).pretty

private def occurs (x e : Expr) : Bool := (e.find? (· == x)).isSome

/-- A solution describes the unknown explicitly, or gives parameter-only conditions. -/
private partial def isolated (x e : Expr) : Bool :=
  if !occurs x e then true
  else if e.isAppOfArity ``Eq 3 then
    let args := e.getAppArgs
    (args[1]! == x && !occurs x args[2]!) ||
      (args[2]! == x && !occurs x args[1]!)
  else if e.isAppOfArity ``Or 2 || e.isAppOfArity ``And 2 then
    e.getAppArgs.all (isolated x)
  else false

private partial def complexity (x e : Expr) : Nat :=
  if !occurs x e then 0
  else if e == x then 1
  else if e.isAppOfArity ``Or 2 || e.isAppOfArity ``And 2 then
    e.getAppArgs.foldl (fun n a => max n (complexity x a)) 0
  else 1 + e.getAppArgs.foldl (fun n a => max n (complexity x a)) 0

private partial def bothSides (x e : Expr) : Bool :=
  if e.isAppOfArity ``Eq 3 then
    let args := e.getAppArgs
    occurs x args[1]! && occurs x args[2]!
  else if e.isAppOfArity ``And 2 || e.isAppOfArity ``Or 2 then
    e.getAppArgs.any (bothSides x)
  else false

/-- Unfold current-module definitions through Lean's reduction machinery.
    Imported mathematical operations retain their identities for theorem matching. -/
private partial def unfoldUser (e : Expr) (fuel : Nat := 64) : MetaM Expr := do
  if fuel == 0 then return e
  let e := e.consumeMData
  if e.isLet then return ← unfoldUser (e.letBody!.instantiate1 e.letValue!) (fuel-1)
  if let .proj type index value := e then
    let projection := Expr.proj type index (← unfoldUser value (fuel-1))
    if let some reduced ← reduceProj? projection then return ← unfoldUser reduced (fuel-1)
    return projection
  let env ← getEnv
  if let .const name _ := e.getAppFn then
    if (env.getModuleIdxFor? name).isNone then
      if let some value ← unfoldDefinition? e (ignoreTransparency := true) then
        return ← unfoldUser value (fuel-1)
  if e.isApp then
    let args ← e.getAppArgs.mapM (fun a => unfoldUser a (fuel-1))
    return (mkAppN e.getAppFn args).headBeta
  return e

private def normalize (e : Expr) : MetaM Simp.Result := do
  let mut theorems := #[← getSimpTheorems]
  for declaration in ← getLCtx do
    if !declaration.isImplementationDetail && (← isProp declaration.type) then
      theorems ← SimpTheoremsArray.addTheorem theorems (.fvar declaration.fvarId) declaration.toExpr
  let ctx ← Simp.mkContext {} (simpTheorems := theorems)
    (congrTheorems := ← getSimpCongrTheorems)
  Mathlib.Meta.NormNum.deriveSimp ctx true e

private def discharge (goal : MVarId) : TermElabM String := do
  let statement ← render (← goal.getType)
  let remaining ← Tactic.run goal do
    Tactic.evalTactic (← `(tactic|
      first | assumption | positivity | (norm_num; done) | (linarith; done)))
  unless remaining.isEmpty do throwError "Undischarged rule premise: {statement}"
  return statement

private structure Node where
  expression : Expr
  equality : Expr
  steps : Array DerivationStep
  saved : Term.SavedState

private def canonical (e proof : Expr) (steps : Array DerivationStep) :
    TermElabM (Expr × Expr × Array DerivationStep) := do
  let r ← normalize e
  let eq ← mkEqTrans proof (← r.getProof)
  let steps ← if r.expr == e then pure steps else
    pure <| steps.push {
      rule := "Lean simplification and exact arithmetic"
      before := ← render e
      after := ← render r.expr }
  return (r.expr, eq, steps)

/-- Bounded best-first search over theorem rewrites. Every edge carries a proof.
    A user theorem registered with equation_rule participates in the same search. -/
private def derive (x original : Expr) : TermElabM (Node × Nat) := do
  let (first, proof, steps) ← canonical original (← mkEqRefl original) #[]
  let start : Node := ⟨first, proof, steps, ← Term.saveState⟩
  let mut queue := #[start]
  let mut visited := ({} : Std.HashSet Expr).insert first
  let mut best := start
  let mut explored := 0
  let rules := ruleExtension.getState (← getEnv)
  while !queue.isEmpty && explored < 96 do
    queue := queue.qsort fun a b => complexity x a.expression < complexity x b.expression
    let some node := queue[0]? | break
    queue := queue.extract 1 queue.size
    node.saved.restore
    explored := explored + 1
    if isolated x node.expression then
      let unfolded ← unfoldUser node.expression
      if unfolded == node.expression then return (node, explored)
      let (e, eq, steps) ← canonical unfolded node.equality <| node.steps.push {
        rule := "Lean unfolds user definitions in the solution"
        before := ← render node.expression
        after := ← render unfolded }
      return (⟨e, eq, steps, ← Term.saveState⟩, explored)
    if complexity x node.expression < complexity x best.expression then best := node
    -- Definition unfolding is an independently verified definitional equality.
    let unfolded ← unfoldUser node.expression
    if unfolded != node.expression then
      let (e, eq, steps) ← canonical unfolded node.equality <|
        node.steps.push {
          rule := "Lean unfolds user definitions"
          before := ← render node.expression
          after := ← render unfolded }
      if !visited.contains e then
        visited := visited.insert e
        queue := queue.push ⟨e, eq, steps, ← Term.saveState⟩
    for rule in rules do
      node.saved.restore
      if rule == ``move_rhs && !bothSides x node.expression then continue
      try
        let goal ← mkFreshExprMVar node.expression
        let result ← goal.mvarId!.rewrite node.expression (← mkConstWithFreshMVarLevels rule)
          (config := { occs := .pos [1] })
        let premises ← result.mvarIds.toArray.mapM discharge
        let rewritten ← instantiateMVars result.eNew
        let edge ← instantiateMVars result.eqProof
        if rewritten.hasMVar || edge.hasMVar then throwError "Unresolved rule arguments"
        let proof ← mkEqTrans node.equality edge
        let (e, eq, steps) ← canonical rewritten proof <| node.steps.push {
          rule := rule.toString
          before := ← render node.expression
          after := ← render rewritten
          premises }
        if !visited.contains e && complexity x e ≤ complexity x node.expression + 1 then
          visited := visited.insert e
          queue := queue.push ⟨e, eq, steps, ← Term.saveState⟩
      catch _ => pure ()
    -- Ring normalization supplies a general polynomial transformation, not a template.
    node.saved.restore
    try
      let normalized ← Mathlib.Tactic.RingNF.M.run (← IO.mkRef {}) {} <|
        Mathlib.Tactic.RingNF.rewrite node.expression
      if normalized.expr != node.expression then
        let (e, eq, steps) ← canonical normalized.expr
          (← mkEqTrans node.equality (← normalized.getProof)) <| node.steps.push {
            rule := "Lean ring normalization"
            before := ← render node.expression
            after := ← render normalized.expr }
        if !visited.contains e then
          visited := visited.insert e
          queue := queue.push ⟨e, eq, steps, ← Term.saveState⟩
    catch _ => pure ()
  best.saved.restore
  return (best, explored)

private def compileEquation (name : Name) (binderStx : TSyntax `ident)
    (type condition : TSyntax `term) : CommandElabM Unit := do
  let declarationName := (← getCurrNamespace) ++ name
  let record ← runTermElabM fun sectionVariables => do
    let term ← `(fun ($binderStx:ident : $type) => ($condition : Prop))
    let function ← Term.elabTerm term none
    Term.synthesizeSyntheticMVarsNoPostponing
    let function ← instantiateMVars function
    lambdaTelescope function fun xs original => do
      unless xs.size == 1 do throwError "Expected exactly one solve variable"
      let x := xs[0]!
      let (node, explored) ← derive x original
      node.saved.restore
      let proof ← instantiateMVars node.equality
      let equivalence ← mkAppM ``Iff.of_eq #[proof]
      let arguments := sectionVariables ++ xs
      let closedProof ← instantiateMVars (← mkLambdaFVars arguments equivalence)
      let closedType ← instantiateMVars (← inferType closedProof)
      unless !closedProof.hasMVar && !closedType.hasMVar do
        throwError "Unresolved metavariables in the compiled derivation"
      withOptions (fun o => o.setBool `Elab.async false |>.setBool `debug.skipKernelTC false) do
        addDecl <| .thmDecl {
          name := declarationName
          levelParams := ← Term.getLevelNames
          type := closedType
          value := closedProof }
      let axioms ← collectAxioms declarationName
      let unexpected := axioms.filter fun n => !([``propext, ``Classical.choice, ``Quot.sound].contains n)
      unless unexpected.isEmpty do
        throwError "The derivation depends on unapproved axioms: {unexpected}"
      return {
        name := declarationName.toString
        status := if isolated x node.expression then "solved" else "unresolved"
        input := ← render original
        solution := ← render node.expression
        theoremType := ← render closedType
        proof := ← render closedProof
        axioms := axioms.map Name.toString
        steps := node.steps
        explored
        lowering := ← lowerEquation x (← unfoldUser original) (← unfoldUser node.expression)
          (isolated x node.expression) : EquationResult }
  modifyEnv fun env => resultExtension.addEntry env record
  logInfo m!"{record.name}: {record.status}\n{record.solution}"

syntax (name := solveNamed) "solve " ident " for " ident " : " term " where " term : command
syntax (name := solveAnonymous) "solve " "for " ident " : " term " where " term : command

elab_rules : command
  | `(solve $name:ident for $binderStx:ident : $type where $condition) =>
    compileEquation name.getId binderStx type condition
  | `(solve for $binderStx:ident : $type where $condition) => do
    let name ← liftCoreM <| mkFreshUserName `solution
    compileEquation name binderStx type condition

/-- The command-line adapter appends this command after the user's Lean module. -/
elab "#export_equations" : command => do
  let results := resultExtension.getState (← getEnv)
  liftIO <| IO.println ("SSC_LEAN_REPORT " ++ (toJson results).compress)

end DesignCompiler.Equation
