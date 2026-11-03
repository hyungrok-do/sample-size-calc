import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Real.Sqrt
import Lean

open Lean Meta

namespace DesignCompiler.Equation

/-- C++ expressions are lowered from elaborated, typed Lean expressions only. -/
structure Lowering where
  kind : String := "unavailable"
  expression : String := ""
  roots : Array String := #[]
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

private def requireStandard (e : Expr) (constant : Name) (args : Array Expr) : MetaM Unit := do
  unless ← isDefEq e (← mkAppM constant args) do
    throwError "C++ lowering requires standard real arithmetic for {← ppExpr e}"

private partial def numeric (x e : Expr) (fuel : Nat := 256) : MetaM String := do
  if fuel == 0 then throwError "Expression exceeds the C++ lowering limit"
  if e == x then return "x"
  let e := e.consumeMData
  unless ← isDefEq (← inferType e) (mkConst ``Real) do
    throwError "C++ evaluation currently supports real-valued expressions"
  let args := e.getAppArgs
  let .const name _ := e.getAppFn | throwError "No C++ lowering for {← ppExpr e}"
  if name == ``OfNat.ofNat && args.size == 3 then
    if let .lit (.natVal n) := args[1]! then
      if n == 0 then requireStandard e ``realZero #[]
      else if n == 1 then requireStandard e ``realOne #[]
      else requireStandard e ``realNat #[args[1]!]
      return s!"{n}.0"
  let lower := fun a => numeric x a (fuel-1)
  if [``HAdd.hAdd, ``HSub.hSub, ``HMul.hMul, ``HDiv.hDiv].contains name && args.size ≥ 2 then
    let constant := if name == ``HAdd.hAdd then ``realAdd else if name == ``HSub.hSub then ``realSub
      else if name == ``HMul.hMul then ``realMul else ``realDiv
    requireStandard e constant #[args[args.size-2]!, args[args.size-1]!]
    let a ← lower args[args.size-2]!
    let b ← lower args[args.size-1]!
    if name == ``HDiv.hDiv then return s!"lean_div({a}, {b})"
    let op := if name == ``HAdd.hAdd then "+" else if name == ``HSub.hSub then "-" else "*"
    return s!"({a} {op} {b})"
  if name == ``HPow.hPow && args.size ≥ 2 then
    let exponent := args[args.size-1]!
    unless ← isDefEq (← inferType exponent) (mkConst ``Nat) do
      throwError "Only natural-number powers are lowered to C++"
    let some n ← getNatValue? exponent | throwError "The power must be a constant natural number"
    if n > 10000 then throwError "The power exceeds the C++ evaluation limit"
    requireStandard e ``realPow #[args[args.size-2]!, exponent]
    return s!"std::pow({← lower args[args.size-2]!}, {n}.0)"
  if name == ``Neg.neg && args.size > 0 then
    requireStandard e ``realNeg #[args.back!]
    return s!"(-{← lower args.back!})"
  if name == ``abs && args.size > 0 then requireStandard e ``realAbs #[args.back!]
  let functions := [( ``Real.exp, "std::exp"), (``Real.log, "lean_log"),
    (``Real.sqrt, "lean_sqrt"), (``Real.sin, "std::sin"), (``Real.cos, "std::cos"),
    (``abs, "std::abs")]
  for (constant, function) in functions do
    if name == constant && args.size > 0 then return s!"{function}({← lower args.back!})"
  throwError "No C++ lowering for {← ppExpr e}"

private partial def roots (x e : Expr) : MetaM (Array String) := do
  if e.isAppOfArity ``Or 2 then
    return (← roots x e.getAppArgs[0]!) ++ (← roots x e.getAppArgs[1]!)
  if e.isConstOf ``False then return #[]
  if e.isAppOfArity ``Eq 3 then
    let args := e.getAppArgs
    if args[1]! == x then return #[← numeric x args[2]!]
    if args[2]! == x then return #[← numeric x args[1]!]
  throwError "C++ root evaluation requires a finite set of closed real solutions"

def lowerEquation (x original solution : Expr) (solved : Bool) : MetaM Lowering := do
  try
    unless ← isDefEq (← inferType x) (mkConst ``Real) do
      throwError "C++ evaluation currently supports a real unknown"
    if (original.find? fun e => e.isFVar && e != x).isSome then
      throwError "Provide values for parameters before generating a C++ executable"
    if solved then
      let values ← roots x solution
      return { kind := "analytic_roots", roots := values }
    unless original.isAppOfArity ``Eq 3 do
      throwError "Numerical fallback requires one real equality"
    let args := original.getAppArgs
    let lhs ← numeric x args[1]!
    let rhs ← numeric x args[2]!
    return { kind := "bracketed_bisection", expression := s!"({lhs} - {rhs})" }
  catch error =>
    return { reason := (← error.toMessageData.toString) }

end DesignCompiler.Equation
