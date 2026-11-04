import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Real.Sqrt
import Lean

open Lean Meta

namespace DesignCompiler.Equation

structure NumericParameter where
  name : String
  cppName : String
  deriving ToJson

structure NumericAssumption where
  statement : String
  expression : String
  deriving ToJson

/-- One disjunct of an explicit solution. No root means all real values when
    the parameter-only condition holds. Conjunctions intersect these sets. -/
structure NumericBranch where
  condition : String := "true"
  root : Option String := none
  deriving ToJson

/-- C++ expressions are lowered from elaborated, typed Lean expressions only. -/
structure Lowering where
  kind : String := "unavailable"
  expression : String := ""
  roots : Array String := #[]
  parameters : Array NumericParameter := #[]
  assumptions : Array NumericAssumption := #[]
  branches : Array NumericBranch := #[]
  reason : String := ""
  deriving ToJson

-- Fix the mathematical operations at compiler-build time. User notation and
-- local typeclass instances must not silently change the meaning of emitted C++.
private noncomputable def realAdd (a b : ℝ) := a + b
private noncomputable def realSub (a b : ℝ) := a - b
private noncomputable def realMul (a b : ℝ) := a * b
private noncomputable def realDiv (a b : ℝ) := a / b
private noncomputable def realPow (a : ℝ) (n : ℕ) := a ^ n
private noncomputable def realNeg (a : ℝ) := -a
private noncomputable def realAbs (a : ℝ) := |a|
private noncomputable def realNat (n : ℕ) : ℝ := n
private noncomputable def realZero : ℝ := 0
private noncomputable def realOne : ℝ := 1
private def realLe (a b : ℝ) := a ≤ b
private def realLt (a b : ℝ) := a < b

private def requireStandard (e : Expr) (constant : Name) (args : Array Expr) : MetaM Unit := do
  unless ← isDefEq e (← mkAppM constant args) do
    throwError "C++ lowering requires standard real arithmetic for {← ppExpr e}"

private partial def numeric (parameters : Array Expr) (x e : Expr)
    (fuel : Nat := 256) : MetaM String := do
  if fuel == 0 then throwError "Expression exceeds the C++ lowering limit"
  let e := e.consumeMData
  unless ← isDefEq (← inferType e) (mkConst ``Real) do
    throwError "C++ evaluation currently supports real-valued expressions"
  if e == x then return "x"
  for i in [:parameters.size] do
    if e == parameters[i]! then return s!"p{i}"
  let args := e.getAppArgs
  let .const name _ := e.getAppFn | throwError "No C++ lowering for {← ppExpr e}"
  if name == ``OfNat.ofNat && args.size == 3 then
    if let .lit (.natVal n) := args[1]! then
      if n == 0 then requireStandard e ``realZero #[]
      else if n == 1 then requireStandard e ``realOne #[]
      else requireStandard e ``realNat #[args[1]!]
      return s!"checked({n}.0)"
  let lower := fun a => numeric parameters x a (fuel-1)
  if [``HAdd.hAdd, ``HSub.hSub, ``HMul.hMul, ``HDiv.hDiv].contains name && args.size ≥ 2 then
    let constant := if name == ``HAdd.hAdd then ``realAdd else if name == ``HSub.hSub then ``realSub
      else if name == ``HMul.hMul then ``realMul else ``realDiv
    requireStandard e constant #[args[args.size-2]!, args[args.size-1]!]
    let a ← lower args[args.size-2]!
    let b ← lower args[args.size-1]!
    if name == ``HDiv.hDiv then return s!"checked(lean_div({a}, {b}))"
    let op := if name == ``HAdd.hAdd then "+" else if name == ``HSub.hSub then "-" else "*"
    return s!"checked({a} {op} {b})"
  if name == ``HPow.hPow && args.size ≥ 2 then
    let exponent := args[args.size-1]!
    unless ← isDefEq (← inferType exponent) (mkConst ``Nat) do
      throwError "Only natural-number powers are lowered to C++"
    let some n ← getNatValue? exponent | throwError "The power must be a constant natural number"
    if n > 10000 then throwError "The power exceeds the C++ evaluation limit"
    requireStandard e ``realPow #[args[args.size-2]!, exponent]
    return s!"checked(std::pow({← lower args[args.size-2]!}, {n}.0))"
  if name == ``Neg.neg && args.size > 0 then
    requireStandard e ``realNeg #[args.back!]
    return s!"checked(-{← lower args.back!})"
  if name == ``abs && args.size > 0 then requireStandard e ``realAbs #[args.back!]
  let functions := [( ``Real.exp, "std::exp"), (``Real.log, "lean_log"),
    (``Real.sqrt, "lean_sqrt"), (``Real.sin, "std::sin"), (``Real.cos, "std::cos"),
    (``abs, "std::abs")]
  for (constant, function) in functions do
    if name == constant && args.size > 0 then return s!"checked({function}({← lower args.back!}))"
  throwError "No C++ lowering for {← ppExpr e}"

private def occurs (x e : Expr) : Bool := (e.find? (· == x)).isSome

/-- Only explicitly supported real predicates become runtime guards. An
    unsupported hypothesis disables export instead of being erased. -/
private partial def predicate (parameters : Array Expr) (x e : Expr)
    (fuel : Nat := 256) : MetaM String := do
  if fuel == 0 then throwError "Predicate exceeds the C++ lowering limit"
  let e := e.consumeMData
  if occurs x e then throwError "A runtime assumption or branch condition still contains the unknown"
  if e.isConstOf ``True then return "true"
  if e.isConstOf ``False then return "false"
  let args := e.getAppArgs
  let lower := fun a => predicate parameters x a (fuel-1)
  if e.isAppOfArity ``And 2 || e.isAppOfArity ``Or 2 || e.isAppOfArity ``Iff 2 then
    let op := if e.isAppOf ``And then "&&" else if e.isAppOf ``Or then "||" else "=="
    return s!"({← lower args[0]!} {op} {← lower args[1]!})"
  if e.isAppOfArity ``Not 1 then return s!"(!{← lower args[0]!})"
  if e.isAppOfArity ``Eq 3 || e.isAppOfArity ``Ne 3 then
    let op := if e.isAppOf ``Eq then "==" else "!="
    return s!"({← numeric parameters x args[1]!} {op} {← numeric parameters x args[2]!})"
  if e.isAppOfArity ``LE.le 4 || e.isAppOfArity ``LT.lt 4 then
    let isLe := e.isAppOf ``LE.le
    requireStandard e (if isLe then ``realLe else ``realLt) #[args[2]!, args[3]!]
    let op := if isLe then "<=" else "<"
    return s!"({← numeric parameters x args[2]!} {op} {← numeric parameters x args[3]!})"
  throwError "No C++ guard for the hypothesis or condition {← ppExpr e}"

private partial def solutionBranches (parameters : Array Expr) (x e : Expr) :
    MetaM (Array NumericBranch) := do
  if e.isConstOf ``False then return #[]
  if !occurs x e then return #[{ condition := ← predicate parameters x e }]
  if e.isAppOfArity ``Or 2 then
    let left ← solutionBranches parameters x e.getAppArgs[0]!
    let right ← solutionBranches parameters x e.getAppArgs[1]!
    if left.size + right.size > 256 then throwError "Too many explicit solution branches for C++"
    return left ++ right
  if e.isAppOfArity ``And 2 then
    let left ← solutionBranches parameters x e.getAppArgs[0]!
    let right ← solutionBranches parameters x e.getAppArgs[1]!
    if left.size * right.size > 256 then throwError "Too many explicit solution branches for C++"
    let mut result := #[]
    for a in left do
      for b in right do
        let condition := s!"({a.condition} && {b.condition})"
        let branch : NumericBranch := match a.root, b.root with
          | some u, some v => { condition := s!"({condition} && ({u} == {v}))", root := some u }
          | some u, none => { condition, root := some u }
          | none, some v => { condition, root := some v }
          | none, none => { condition }
        result := result.push branch
    return result
  if e.isAppOfArity ``Eq 3 then
    let args := e.getAppArgs
    if args[1]! == x && !occurs x args[2]! then return #[{ root := some (← numeric parameters x args[2]!) }]
    if args[2]! == x && !occurs x args[1]! then return #[{ root := some (← numeric parameters x args[1]!) }]
  throwError "C++ evaluation requires explicit roots with parameter-only conditions"

def lowerEquation (x original solution : Expr) (arguments : Array Expr)
    (hypotheses : Array Expr) (solved : Bool) : MetaM Lowering := do
  try
    unless ← isDefEq (← inferType x) (mkConst ``Real) do
      throwError "C++ evaluation currently supports a real unknown"
    let mut parameterExprs := #[]
    let mut parameters : Array NumericParameter := #[]
    for argument in arguments do
      let type ← inferType argument
      if ← isProp type then continue
      unless ← isDefEq type (mkConst ``Real) do
        throwError "C++ parameters currently must have type ℝ: {← ppExpr argument}"
      let name := (← argument.fvarId!.getDecl).userName.toString
      if name.contains '=' || (parameters.any (·.name == name)) then
        throwError "C++ parameters require distinct names without '='"
      parameters := parameters.push ⟨name, s!"p{parameterExprs.size}"⟩
      parameterExprs := parameterExprs.push argument
    let assumptions ← hypotheses.mapM fun h => do
      return { statement := (← ppExpr h).pretty,
               expression := ← predicate parameterExprs x h : NumericAssumption }
    if solved then
      let branches ← solutionBranches parameterExprs x solution
      let roots := if branches.all (·.condition == "true") then branches.filterMap (·.root) else #[]
      return { kind := "analytic_roots", roots, branches, parameters, assumptions }
    unless original.isAppOfArity ``Eq 3 do
      throwError "Numerical fallback requires one real equality"
    let args := original.getAppArgs
    let lhs ← numeric parameterExprs x args[1]!
    let rhs ← numeric parameterExprs x args[2]!
    return { kind := "bracketed_bisection", expression := s!"checked({lhs} - {rhs})", parameters, assumptions }
  catch error =>
    return { reason := (← error.toMessageData.toString) }

end DesignCompiler.Equation
