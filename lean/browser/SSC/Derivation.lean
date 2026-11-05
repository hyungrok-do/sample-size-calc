import SSC.Finite
import SSC.Calculation

/-!
Proof-producing compilation of finite observation programs. A certificate
equates expectations for EVERY function of the output, not merely one selected
test. Combining equal intermediate states is justified by finite additivity.
-/
namespace SSC.Derivation
open SSC.Finite

/-- A finite observation program. Bind may depend on every previous observation.
Decidable equality is needed only to merge identical intermediate states. -/
inductive Program : Type → Type 1 where
  | pure {α : Type} (value : α) : Program α
  | draw {α : Type} (distribution : Distribution α) : Program α
  | bind {α β : Type} [DecidableEq α] (first : Program α)
      (next : α → Program β) : Program β

namespace Program

def expect : {α : Type} → Program α → (α → Rat) → Rat
  | _, .pure x, f => f x
  | _, .draw d, f => d.expect f
  | _, @Program.bind _ _ _ p k, f => p.expect (fun x => (k x).expect f)

def probability (p : Program α) (event : α → Bool) : Rat :=
  p.expect (fun x => if event x then 1 else 0)

theorem expect_congr (p : Program α) (f g : α → Rat) (h : ∀ x, f x = g x) :
    p.expect f = p.expect g := congrArg p.expect (funext h)

theorem expect_const (p : Program α) (c : Rat) : p.expect (fun _ => c) = c := by
  induction p with
  | pure x => rfl
  | draw d => exact d.expect_const c
  | bind p k ih ikh =>
    change p.expect (fun x => (k x).expect (fun _ => c)) = c
    rw [expect_congr p _ (fun _ => c) (fun x => ikh x)]
    exact ih

/-- The unchanged observation program interpreted by the existing sampler. -/
def toLaw : {α : Type} → Program α → SSC.Probability.Law α
  | _, .pure x => .pure x
  | _, .draw d => .finite d.atoms.toArray
  | _, @Program.bind _ _ _ p k =>
    let suspended := Thunk.mk (fun _ => .bind p.toLaw (fun x => (k x).toLaw))
    .defer (fun _ => suspended.get)

/-- Memoization changes construction cost, not the law or its interpreter steps. -/
theorem toLaw_bind [DecidableEq α] (p : Program α) (k : α → Program β) :
    (p.bind k).toLaw = .defer (fun _ => .bind p.toLaw (fun x => (k x).toLaw)) := rfl

def map [DecidableEq α] (p : Program α) (f : α → β) : Program β :=
  p.bind (fun x => .pure (f x))

@[simp] theorem expect_map [DecidableEq α] (p : Program α) (f : α → β) (g : β → Rat) :
    (p.map f).expect g = p.expect (fun x => g (f x)) := rfl

/-- Independent observations, retaining the complete dataset. -/
def iid [DecidableEq α] (p : Program α) : Nat → Program (List α)
  | 0 => .pure []
  | n+1 => (iid p n).bind (fun xs => p.map (fun x => xs ++ [x]))

/-- The statistic is updated as observations arrive. The compiler may merge
equal states, but the sampler still draws and processes every observation. -/
def iidFold [DecidableEq α] [DecidableEq σ] (p : Program α)
    (step : σ → α → σ) (initial : σ) : Nat → Program σ
  | 0 => .pure initial
  | n+1 => (iidFold p step initial n).bind (fun s => p.map (step s))

/-- General finite-state observation processes, including dependent or adaptive
transitions. No independence assumption is inserted by the compiler. -/
def iterate [DecidableEq σ] (transition : σ → Program σ) (initial : σ) : Nat → Program σ
  | 0 => .pure initial
  | n+1 => (iterate transition initial n).bind transition

-- Keep the original defining equations available to user proofs. The compiler
-- uses tail-recursive construction only through these checked equalities.
private def iidRuntime [DecidableEq α] (p : Program α) (n : Nat) : Program (List α) :=
  Nat.repeat (fun previous => previous.bind (fun xs => p.map (fun x => xs ++ [x]))) n (.pure [])

@[csimp] theorem iid_eq_runtime : @iid = @iidRuntime := by
  funext α inst p n
  induction n with
  | zero => rfl
  | succ n ih =>
    change (iid p n).bind _ = (iidRuntime p n).bind _
    rw [ih]

private def iidFoldRuntime [DecidableEq α] [DecidableEq σ] (p : Program α)
    (step : σ → α → σ) (initial : σ) (n : Nat) : Program σ :=
  Nat.repeat (fun previous => previous.bind (fun s => p.map (step s))) n (.pure initial)

@[csimp] theorem iidFold_eq_runtime : @iidFold = @iidFoldRuntime := by
  funext α σ instα instσ p step initial n
  induction n with
  | zero => rfl
  | succ n ih =>
    change (iidFold p step initial n).bind _ = (iidFoldRuntime p step initial n).bind _
    rw [ih]

private def iterateRuntime [DecidableEq σ] (transition : σ → Program σ)
    (initial : σ) (n : Nat) : Program σ :=
  Nat.repeat (fun previous => previous.bind transition) n (.pure initial)

@[csimp] theorem iterate_eq_runtime : @iterate = @iterateRuntime := by
  funext σ inst transition initial n
  induction n with
  | zero => rfl
  | succ n ih =>
    change (iterate transition initial n).bind _ = (iterateRuntime transition initial n).bind _
    rw [ih]

/-- Generic observation-to-statistic derivation, with arbitrary state and update. -/
theorem iidFold_eq_iid [DecidableEq α] [DecidableEq σ] (p : Program α)
    (step : σ → α → σ) (initial : σ) (n : Nat) (f : σ → Rat) :
    (p.iidFold step initial n).expect f =
      (p.iid n).expect (fun xs => f (xs.foldl step initial)) := by
  induction n generalizing f with
  | zero => rfl
  | succ n ih =>
    simp only [iidFold, iid, expect, expect_map]
    rw [ih]
    apply expect_congr
    intro xs
    apply expect_congr
    intro x
    simp [List.foldl_append]

end Program

structure Certificate (p : Program α) where
  distribution : Distribution α
  correct : ∀ f, distribution.expect f = p.expect f

namespace Certificate

theorem probability_eq {p : Program α} (c : Certificate p) (event : α → Bool) :
    c.distribution.probability event = p.probability event :=
  (c.distribution.probability_eq_expect event).trans (c.correct _)

end Certificate

/-- Total reference compiler, useful in proofs. Production evaluation uses the
bounded compiler below so an expanding state space can fall back to simulation. -/
def derive : {α : Type} → [DecidableEq α] → (p : Program α) → Certificate p
  | _, _, .pure x => ⟨.pure x, fun _ => Distribution.expect_pure _ _⟩
  | _, _, .draw d => ⟨d.compact, fun _ => Distribution.expect_compact _ _⟩
  | _, _, @Program.bind _ _ inst p k =>
    letI := inst
    let left := derive p
    let right := fun x => derive (k x)
    { distribution := (left.distribution.bind (fun x => (right x).distribution)).compact
      correct := by
        intro f
        rw [Distribution.expect_compact, Distribution.expect_bind]
        change expectation left.distribution.atoms (fun x => (right x).distribution.expect f) = _
        rw [expectation_congr _ _ (fun x => (k x).expect f) (fun x => (right x).correct f)]
        exact left.correct _ }

/-- Compute only table data during conversion. Keeping recursive certificate
records out of this path avoids repeatedly expanding their proof fields. -/
def deriveAtoms : {α : Type} → [DecidableEq α] → Program α → Table α
  | _, _, .pure x => [(x, 1)]
  | _, _, .draw d => collect d.atoms
  | _, _, @Program.bind _ _ inst p k =>
    letI := inst
    collect (flatten (deriveAtoms p) (fun x => deriveAtoms (k x)))

theorem deriveAtoms_eq : {α : Type} → [DecidableEq α] → (p : Program α) →
    deriveAtoms p = (derive p).distribution.atoms
  | _, _, .pure _ => rfl
  | _, _, .draw _ => rfl
  | _, _, @Program.bind _ _ inst p k => by
    letI := inst
    exact congrArg collect
      ((congrArg (fun xs => flatten xs (fun x => deriveAtoms (k x))) (deriveAtoms_eq p)).trans
        (congrArg (flatten (derive p).distribution.atoms) (funext (fun x => deriveAtoms_eq (k x)))))

theorem expectation_deriveAtoms [DecidableEq α] (p : Program α) (f : α → Rat) :
    expectation (deriveAtoms p) f = p.expect f := by
  rw [deriveAtoms_eq]
  exact (derive p).correct f

theorem deriveAtoms_nonnegative [DecidableEq α] (p : Program α) :
    Nonnegative (deriveAtoms p) := by
  rw [deriveAtoms_eq]
  exact (derive p).distribution.nonnegative

namespace Certificate

/-- Compose already checked intermediate evidence with an arbitrary next program. -/
def bindDerived [DecidableEq α] [DecidableEq β] {p : Program α}
    (left : Certificate p) (next : α → Program β) : Certificate (p.bind next) where
  distribution := {
    atoms := flattenMerged left.distribution.atoms (fun x => deriveAtoms (next x))
    nonnegative := flattenMerged_nonnegative _ left.distribution.nonnegative _
      (fun x => deriveAtoms_nonnegative (next x))
    normalized := by
      rw [expectation_flattenMerged]
      rw [expectation_congr _ _ (fun _ => 1)
        (fun x => (expectation_deriveAtoms (next x) _).trans (Program.expect_const (next x) 1))]
      exact left.distribution.normalized }
  correct := by
    intro f
    change expectation (flattenMerged left.distribution.atoms (fun x => deriveAtoms (next x))) f =
      p.expect (fun x => (next x).expect f)
    rw [expectation_flattenMerged]
    rw [expectation_congr _ _ (fun x => (next x).expect f)
      (fun x => expectation_deriveAtoms (next x) f)]
    exact left.correct _

/-- A proposed table replaces a computed one only after an equality proof. -/
def repack {p : Program α} (certificate : Certificate p) (atoms : Table α)
    (checked : certificate.distribution.atoms = atoms) : Certificate p where
  distribution := {
    atoms := atoms
    nonnegative := by rw [← checked]; exact certificate.distribution.nonnegative
    normalized := by rw [← checked]; exact certificate.distribution.normalized }
  correct := by
    intro f
    change expectation atoms f = _
    rw [← checked]
    exact certificate.correct f

end Certificate

abbrev Work := StateT Nat Option

private def spend (amount : Nat) : Work Unit := fun remaining =>
  if amount > remaining then none else some ((), remaining - amount)

private structure Branches (xs : Table α) (k : α → Program β) where
  atoms : Table β
  nonnegative : Nonnegative atoms
  correct : ∀ f, expectation atoms f = expectation xs (fun x => (k x).expect f)

private def compileBranches [DecidableEq β] (k : α → Program β)
    (compileNext : (x : α) → Work (Certificate (k x))) :
    (xs : Table α) → Nonnegative xs → Work (Branches xs k)
  | [], _ => return ⟨[], by intro x hx; contradiction, fun _ => rfl⟩
  | (x, p) :: xs, hx => do
    let head ← compileNext x
    let tail ← compileBranches k compileNext xs (fun y hy => hx y (by simp [hy]))
    let scaled := scale p head.distribution.atoms
    -- Bound list construction and the worst-case equality comparisons in merge.
    let count := scaled.length + tail.atoms.length
    spend (scaled.length * count + 1)
    let atoms := merge scaled tail.atoms
    return {
      atoms := atoms
      nonnegative := merge_nonnegative _ _
        (scale_nonnegative p (hx (x, p) (by simp)) _ head.distribution.nonnegative) tail.nonnegative
      correct := by
        intro f
        simp only [atoms, expectation_merge, scaled,
          expectation_scale_weights, expectation]
        rw [show expectation head.distribution.atoms f = (k x).expect f from head.correct f,
          tail.correct f] }

/-- A bounded proof-producing compiler. Failure produces no certificate. All
successful branches carry the same all-expectations theorem as `derive`. -/
private def compileDepth : Nat → {α : Type} → [DecidableEq α] → (p : Program α) → Work (Certificate p)
  | 0, _, _, _ => fun _ => none
  | _+1, _, _, .pure x => do
    spend 1
    return ⟨.pure x, fun _ => Distribution.expect_pure _ _⟩
  | _+1, _, _, .draw d => do
    spend (d.atoms.length * d.atoms.length + 1)
    return ⟨d.compact, fun _ => Distribution.expect_compact _ _⟩
  | depth+1, _, _, @Program.bind _ _ inst p k => do
    letI := inst
    spend 1
    let left ← compileDepth depth p
    let branches ← compileBranches k (fun x => compileDepth depth (k x))
      left.distribution.atoms left.distribution.nonnegative
    return {
      distribution := {
        atoms := branches.atoms
        nonnegative := branches.nonnegative
        normalized := by
          rw [branches.correct]
          change left.distribution.expect (fun x => (k x).expect (fun _ => 1)) = 1
          rw [left.correct]
          exact Program.expect_const (.bind p k) 1 }
      correct := by
        intro f
        change expectation branches.atoms f = _
        rw [branches.correct]
        exact left.correct _ }

/-- Bound interpreter recursion as well as table work. Deep programs remain
available to the stack-safe sampler; no truncated certificate is returned. -/
def compile {α : Type} [DecidableEq α] (p : Program α) : Work (Certificate p) :=
  compileDepth 64 p

def compileWithin [DecidableEq α] (p : Program α) (budget : Nat := 1000000) :
    Option (Certificate p) := ((compile p).run budget).map Prod.fst

def derivedProbability [DecidableEq α] (p : Program α) (event : α → Bool)
    (budget : Nat := 1000000) : Option Rat :=
  (compileWithin p budget).map (fun c => c.distribution.probability event)

/-- A reported exact probability can be connected to the observation program by
checking its computation, then applying this generic soundness theorem. -/
theorem derivedProbability_sound [DecidableEq α] (p : Program α) (event : α → Bool)
    (budget : Nat) (q : Rat) (h : derivedProbability p event budget = some q) :
    p.probability event = q := by
  unfold derivedProbability at h
  cases hc : compileWithin p budget with
  | none => simp [hc] at h
  | some c =>
    simp only [hc, Option.map_some, Option.some.injEq] at h
    exact (c.probability_eq event).symm.trans h

/-- Preparing a rule at a design is separate from applying it to each dataset.
This also permits joint endpoint rules under the specified joint model. -/
inductive Rejection (α : Type) where
  | test (accept : α → Bool)
  | both (left right : Rejection α)
  | either (left right : Rejection α)

def Rejection.accepts (rule : Rejection α) (x : α) : Bool :=
  match rule with
  | .test accept => accept x
  | .both left right => left.accepts x && right.accepts x
  | .either left right => left.accepts x || right.accepts x

def Rejection.all (rules : List (Rejection α)) : Rejection α :=
  rules.foldr .both (.test (fun _ => true))

def Rejection.any (rules : List (Rejection α)) : Rejection α :=
  rules.foldr .either (.test (fun _ => false))

/-- The same source programs supply the formal probabilities, exact compilation,
and observation-by-observation simulation. The rejection rule is also shared. -/
structure Specification (Design Outcome : Type) where
  null : Design → Program Outcome
  alternative : Design → Program Outcome
  reject : Design → Rejection Outcome

namespace Specification

def typeI (s : Specification D X) (d : D) : Rat := (s.null d).probability (s.reject d).accepts
def power (s : Specification D X) (d : D) : Rat := (s.alternative d).probability (s.reject d).accepts

private def lower [DecidableEq X] (p : Program X) (budget : Nat) (simulate : Bool) :
    SSC.Probability.Law X :=
  if simulate then p.toLaw
  else match compileWithin p budget with
    | some certificate => .finite certificate.distribution.atoms.toArray
    | none => .samplingOnly "Finite distribution derivation exceeded its work budget; sampling the original observation program." p.toLaw

def toCalculation [DecidableEq X] (s : Specification D X) (budget : Nat := 1000000)
    (simulate : Bool := false) : SSC.Calculation.Specification D X where
  null := fun d => lower (s.null d) budget simulate
  alternative := fun d => lower (s.alternative d) budget simulate
  reject := fun d => .ofBool (s.reject d).accepts

def evaluate [DecidableEq X] (s : Specification D X) (design : D)
    (options : SSC.Calculation.Options := {}) (budget : Nat := 1000000) :=
  SSC.Calculation.evaluate (s.toCalculation budget (options.route == .simulation)) design options

def search [DecidableEq X] (s : Specification D X) (candidates : Array D)
    (requirements : SSC.Calculation.Requirements := {})
    (options : SSC.Calculation.Options := {}) (budget : Nat := 1000000) :=
  SSC.Calculation.search (s.toCalculation budget (options.route == .simulation)) candidates requirements options

end Specification

/-!
The optional tactic stages closed computations into small, independently checked
certificates. Evaluation only PROPOSES literal tables. Each proposal is checked
with an ordinary kernel reduction before it can be used in the next step.
-/
namespace Staging
open Lean Meta Elab Tactic

register_option ssc.profileFinite : Bool := {
  defValue := false
  descr := "Report finite certificate evaluation, quotation, kernel checking and registration times." }

private structure StagedCache where
  certificates : Std.HashMap Expr Expr := {}
  rationals : Std.HashMap Rat Expr := {}

private abbrev CheckedPrefixes := IO.Ref StagedCache

private structure QuotedAtom where
  value : Expr
  weight : Rat

private def quoteTable [ToExpr α] (atoms : Table α) : Array QuotedAtom :=
  atoms.toArray.map fun (value, weight) => ⟨toExpr value, weight⟩

private def checkedRational (memo : CheckedPrefixes) (q : Rat) : MetaM Expr := do
  if let some value := (← memo.get).rationals[q]? then return value
  let num := match q.num with
    | .ofNat n => mkApp (mkConst ``Int.ofNat) (mkNatLit n)
    | .negSucc n => mkApp (mkConst ``Int.negSucc) (mkNatLit n)
  let den := mkNatLit q.den
  let nonzero ← mkAppM ``Ne #[den, mkNatLit 0]
  let coprime ← mkAppM ``Nat.Coprime #[mkNatLit q.num.natAbs, den]
  let value := mkApp4 (mkConst ``Rat.mk') num den
    (← mkDecideProof nonzero) (← mkDecideProof coprime)
  let name ← mkAuxDeclName (kind := `_finiteRat)
  -- The evaluator supplies only integer data. The kernel checks both rational
  -- invariants once; later tables project data without repeating division/GCD.
  withOptions (fun o => Lean.debug.skipKernelTC.set (Lean.Elab.async.set o false) false) do
    addAndCompile (.defnDecl {
      name := name
      levelParams := []
      type := mkConst ``Rat
      value := value
      hints := .regular 0
      safety := .safe })
  let result := mkConst name
  memo.modify fun s => { s with rationals := s.rationals.insert q result }
  return result

private unsafe def checkpoint (memo : CheckedPrefixes) (certificate : Expr) : MetaM Expr := do
  let started ← IO.monoMsNow
  let distribution ← mkAppM ``Certificate.distribution #[certificate]
  let atoms ← mkAppM ``Distribution.atoms #[distribution]
  let quotation ← mkAppM ``quoteTable #[atoms]
  let quoted ← evalExpr (Array QuotedAtom)
    (mkApp (mkConst ``Array [.zero]) (mkConst ``QuotedAtom)) quotation
  let evaluated ← IO.monoMsNow
  let α := (← inferType distribution).getAppArgs[0]!
  let pairType := mkApp2 (mkConst ``Prod [.zero, .zero]) α (mkConst ``Rat)
  let mut proposed := mkApp (mkConst ``List.nil [.zero]) pairType
  for atom in quoted.reverse do
    let pair := mkApp4 (mkConst ``Prod.mk [.zero, .zero]) α (mkConst ``Rat)
      atom.value (← checkedRational memo atom.weight)
    proposed := mkApp3 (mkConst ``List.cons [.zero]) pairType pair proposed
  let proposedAt ← IO.monoMsNow
  let proposition ← mkEq atoms proposed
  -- addDecl checks this proof in the kernel. No evaluator-trust axiom is used.
  let checked ← withOptions (fun o => Lean.debug.skipKernelTC.set (Lean.Elab.async.set o false) false) do
    -- Compare the actual computed table directly with its quoted normal form.
    -- Reflecting through DecidableEq builds and reduces a second list/rational
    -- equality program. Reflexivity needs only the kernel's conversion check.
    -- Unusual ToExpr/DecidableEq instances may require propositional equality,
    -- so retain the original checked path after restoring failed speculation.
    let saved ← saveState
    try
      mkAuxLemma [] proposition (← mkEqRefl proposed) (kind? := some `_finiteEquality)
    catch _ =>
      saved.restore
      mkAuxLemma [] proposition (← mkDecideProof proposition) (kind? := some `_finiteEquality)
  let checkedAt ← IO.monoMsNow
  let value ← mkAppM ``Certificate.repack #[certificate, proposed, mkConst checked]
  let name ← mkAuxDeclName (kind := `_finiteCheckpoint)
  let type ← inferType value
  withOptions (fun o => Lean.debug.skipKernelTC.set (Lean.Elab.async.set o false) false) do
    addAndCompile (.defnDecl {
      name := name
      levelParams := []
      type := type
      value := value
      hints := .regular 0
      safety := .safe })
  if ssc.profileFinite.get (← getOptions) then
    logInfo m!"finite_derive: {quoted.size} states; evaluate={evaluated-started}ms quote={proposedAt-evaluated}ms kernel={checkedAt-proposedAt}ms register={(← IO.monoMsNow)-checkedAt}ms"
  return mkConst name

-- The memo is local to one tactic invocation. It contains only closed,
-- kernel-checked certificates in this environment, never evaluator answers.
-- Structural expression equality resolves hash collisions; no model names or
-- distribution families participate in lookup.
private unsafe def stage (memo : CheckedPrefixes) : Nat → Expr → MetaM Expr
  | 0, _ => throwError "finite_derive exceeded 256 intermediate checkpoints. Split the program into proved parts."
  | fuel+1, program => do
    let reduced ← whnf program
    if let some certificate := (← memo.get).certificates[reduced]? then
      return certificate
    let certificate ← if reduced.isAppOfArity ``Program.bind 5 then do
        let args := reduced.getAppArgs
        let left ← stage memo fuel args[3]!
        mkAppM ``Certificate.bindDerived #[left, args[4]!]
      else mkAppM ``derive #[program]
    let checked ← checkpoint memo certificate
    memo.modify fun s => { s with certificates := s.certificates.insert reduced checked }
    return checked

private unsafe def stageGoal (memo : CheckedPrefixes) (remaining : IO.Ref Nat) : Nat → TacticM Unit
  | 0 => throwError "finite_derive supports at most 256 finite-family cases. Split the family into proved parts."
  | fuel+1 => withMainContext do
    let goal ← getMainGoal
    let target ← whnf (← instantiateMVars (← goal.getType))
    if target.isForall then
      let domain ← whnf target.bindingDomain!
      unless domain.isAppOfArity ``Fin 1 do
        throwError "finite_derive expects a closed Certificate or a family indexed by Fin."
      let bound ← whnf domain.getAppArgs[0]!
      let .lit (.natVal n) := bound | throwError "The finite-family bound must be a closed natural number."
      if n > fuel then
        throwError "finite_derive supports at most 256 finite-family cases. Split the family into proved parts."
      if n == 0 then
        evalTactic (← `(tactic| exact fun i => Fin.elim0 i))
      else
        evalTactic (← `(tactic| refine Fin.cases ?_ ?_))
        let goals ← getGoals
        for next in goals do
          setGoals [next]
          stageGoal memo remaining fuel
      return
    unless target.isAppOfArity ``Certificate 2 do
      throwError "finite_derive expects a goal of type SSC.Derivation.Certificate program."
    let program := target.getAppArgs[1]!
    if program.hasFVar || program.hasMVar then
      throwError "finite_derive needs a closed program. Use symbolic derivation theorems for programs with parameters."
    if (← remaining.get) == 0 then
      throwError "finite_derive supports at most 256 finite-family cases. Split the family into proved parts."
    remaining.modify (· - 1)
    let certificate ← stage memo 256 program
    goal.assign certificate
    replaceMainGoal []

/-- Check intermediate finite distributions separately, then compose their
proofs. A `(i : Fin n) → Certificate (program i)` goal shares checked prefixes
across the finite family. Intermediate types need `Lean.ToExpr` instances. -/
syntax (name := finiteDerive) "finite_derive" : tactic

@[tactic finiteDerive] unsafe def evalFiniteDerive : Tactic := fun _ => do
  let goal ← getMainGoal
  let rest := (← getGoals).tail
  setGoals [goal]
  let memo ← IO.mkRef ({} : StagedCache)
  let remaining ← IO.mkRef 256
  stageGoal memo remaining 257
  setGoals rest

end Staging

end SSC.Derivation
