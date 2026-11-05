import SSC
open SSC.Derivation SSC.Finite SSC.Probability
set_option maxRecDepth 10000
set_option maxHeartbeats 8000000

def transition (s : Nat) : Program Nat :=
  (Program.draw (Distribution.bernoulli (1/3) (by decide +kernel))).map
    (fun x => if x then s + 2 else s / 2)
def process (n : Nat) := Program.iterate transition 0 n
def family : (n : Fin 9) → Certificate (process n.val) := by finite_derive
def emptyFamily : (n : Fin 0) → Certificate (process n.val) := by finite_derive

theorem family_correct (n : Fin 9) (event : Nat → Bool) :
    (family n).distribution.probability event = (process n.val).probability event :=
  (family n).probability_eq event

run_cmd do
  let env ← Lean.getEnv
  let count : Nat := env.checked.get.constants.foldStage2 (fun count name _ =>
    if (name.toString.splitOn "_finiteCheckpoint").length > 1 then count + 1 else count) 0
  Lean.logInfo m!"SSC_FINITE_CHECKPOINTS {count}"

-- Multiple goals must survive a tactic invocation. Nested finite families are
-- also checked from their actual closed programs, including repeated prefixes.
def twoGoals : Certificate (Program.pure (3 : Nat)) × Certificate (Program.pure (4 : Nat)) := by
  constructor
  · finite_derive
  · finite_derive
def nestedFamily : (i : Fin 3) → (j : Fin 2) → Certificate (Program.pure (i.val + j.val)) := by
  finite_derive

-- A lawful quotation may select another representative of the same quotient.
-- These values are propositionally equal, but kernel conversion alone cannot
-- identify them. The checked DecidableEq fallback must remain available.
def universal : Setoid Bool where
  r := fun _ _ => True
  iseqv := ⟨fun _ => trivial, fun _ => trivial, fun _ _ => trivial⟩
abbrev Label := Quotient universal
def labelFalse : Label := Quotient.mk universal false
def labelTrue : Label := Quotient.mk universal true
instance : DecidableEq Label := fun a b =>
  .isTrue (Quotient.inductionOn₂ a b (fun _ _ => Quotient.sound (by trivial)))
instance : Lean.ToExpr Label where
  toTypeExpr := Lean.mkConst ``Label
  toExpr _ := Lean.mkConst ``labelFalse
def quotientCertificate : Certificate (Program.pure labelTrue) := by finite_derive
theorem quotientCertificate_correct :
    quotientCertificate.distribution.probability (fun _ => true) = 1 := by
  rw [quotientCertificate.probability_eq]
  rfl

-- Independent copy of the original finite sampling interpreter. This checks
-- exact draw sequences and work charges, including adaptive weight changes,
-- cache eviction, reuse of weights with different values, and failures.
def oldTick : Draw Unit := do
  let s ← get
  if s.remaining == 0 then throw "Simulation work budget exhausted; no estimate was returned."
  set { s with remaining := s.remaining - 1 }

def oldSample : Nat → {α : Type} → Law α → Draw α
  | 0, _, _ => throw "Reference interpreter fuel exhausted"
  | fuel+1, _, .defer next => oldSample fuel (next ())
  | fuel+1, _, .samplingOnly _ model => oldSample fuel model
  | fuel+1, _, .bind first next => do
    oldTick
    let x ← oldSample fuel first
    oldSample fuel (next x)
  | _+1, _, .pure x => do oldTick; return x
  | _+1, _, .finite xs => do
    if xs.isEmpty then throw "A finite law must have positive total mass."
    if xs.any (fun x => x.2 < 0) then throw "Probability weights cannot be negative."
    if xs.foldl (fun s x => s+x.2) 0 != 1 then throw "Probability weights must sum exactly to one."
    let u ← uniform01
    let mut mass : Rat := 0
    for (x, p) in xs do
      mass := mass+p
      if u < probabilityFloat mass then return x
    throw "Finite sampler failed to accumulate unit mass."
  | _+1, _, .invalid message => throw message
  | _+1, _, _ => throw "Reference test covers finite laws"

def oldCount (law : Law Nat) (accept : Nat → SSC.Numeric.Result Bool) (replicates : Nat) : Draw Nat := do
  let mut count := 0
  for _ in [:replicates] do
    if ← accept (← oldSample 100 law) then count := count + 1
  return count

def oldToLaw : {α : Type} → Program α → Law α
  | _, .pure x => .pure x
  | _, .draw d => .finite d.atoms.toArray
  | _, @Program.bind _ _ _ p k => .defer (fun _ => .bind (oldToLaw p) (fun x => oldToLaw (k x)))

theorem memoized_law_is_original : {α : Type} → (p : Program α) → p.toLaw = oldToLaw p
  | _, .pure _ => rfl
  | _, .draw _ => rfl
  | _, @Program.bind _ _ _ p k =>
    congrArg Law.defer (funext (fun _ =>
      (congrArg (fun first => Law.bind first (fun x => (k x).toLaw))
        (memoized_law_is_original p)).trans
      (congrArg (Law.bind (oldToLaw p)) (funext (fun x => memoized_law_is_original (k x))))))

def sameProgramRun (p : Program Nat) (seed : UInt64) (budget : Nat) : Bool :=
  let accept := fun n => .ok (n % 3 == 0)
  match (countSamples p.toLaw accept 100).run ⟨seed, budget⟩,
      (oldCount (oldToLaw p) accept 100).run ⟨seed, budget⟩ with
  | .ok (a, sa), .ok (b, sb) => a == b && sa.seed == sb.seed && sa.remaining == sb.remaining
  | .error a, .error b => a == b
  | _, _ => false

def sameRun (law : Law Nat) (accept : Nat → SSC.Numeric.Result Bool)
    (seed : UInt64) (budget : Nat) (replicates : Nat := 1000) : Bool :=
  match (countSamples law accept replicates).run ⟨seed, budget⟩,
      (oldCount law accept replicates).run ⟨seed, budget⟩ with
  | .ok (a, sa), .ok (b, sb) => a == b && sa.seed == sb.seed && sa.remaining == sb.remaining
  | .error a, .error b => a == b
  | _, _ => false

def adaptive : Law Nat := do
  let branch ← Law.finite ((List.range 12).map (fun n => (n, (1/12 : Rat)))).toArray
  let p : Rat := ((branch + 1 : Nat) : Rat) / 13
  let x ← Law.finite #[(branch, p), (branch + 20, 1-p)]
  pure (x + branch)

#eval (show IO Unit from do
  let laws := #[Law.finite #[(0, 0), (1, 1/3), (2, 2/3)], adaptive,
    Law.finite #[], Law.finite #[(0, -1), (1, 2)], Law.finite #[(0, 1/2)]]
  for law in laws do
    for seed in [1, 41001, 99167] do
      for budget in [0, 25, 1000000] do
        unless sameRun law (fun n => .ok (n % 3 == 0)) seed budget do
          throw (IO.userError "Cached sampling changed a draw, seed, budget or error")
  unless sameRun adaptive (fun _ => .error "Estimator failed") 1 1000000 do
    throw (IO.userError "Estimator failures changed")
  let large : Law Nat := .finite ((List.range 1025).map
    (fun n => (n, if n == 1024 then (1 : Rat) else 0))).toArray
  unless sameRun large (fun n => .ok (n == 1024)) 123 10000 3 do
    throw (IO.userError "Large-support cache bypass changed sampling")
  for i in [:9] do
    if h : i < 9 then
      let event := fun n => n % 3 == 1
      unless (family ⟨i, h⟩).distribution.probability event ==
          (derive (process i)).distribution.probability event do
        throw (IO.userError "A finite-family certificate has the wrong result")
  for seed in [1, 41001, 99167] do
    for budget in [0, 2, 15, 100000] do
      unless sameProgramRun (process 8) seed budget do
        throw (IO.userError "Memoized program construction changed observations or work charges")
  IO.println "SSC_PERFORMANCE_SEMANTICS_OK")
