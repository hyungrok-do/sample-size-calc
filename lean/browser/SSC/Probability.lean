import SSC.Numeric

/-!
A compositional probability program has two interpreters: finite rational
enumeration and sampling. Direct scalar laws also have numerical CDFs. A
specification is never translated into another statistical model for simulation.

This is an executable probability library, not a formalized measure theory.
Finite probabilities use exact rational arithmetic. Continuous laws use floating
point samplers; their numerical correctness is validated separately from proofs.
-/
namespace SSC.Probability
open SSC.Numeric

inductive ScalarLaw where
  | normal (mean sd : Float)
  | uniform (lo hi : Float)
  | exponential (rate : Float)
  | gamma (shape scale : Float)
  | beta (a b : Float)
  | studentT (df noncentrality : Float)
  | chiSquare (df noncentrality : Float)
  | fisherF (df1 df2 noncentrality : Float)
  deriving Repr

def ScalarLaw.validate : ScalarLaw → Result Unit
  | .normal m s => if m.isFinite && s.isFinite && s > 0 then .ok () else .error "Normal law needs a finite mean and positive SD."
  | .uniform a b => if a.isFinite && b.isFinite && a < b && (b-a).isFinite then .ok () else .error "Uniform endpoints must be finite and increasing, with a finite span."
  | .exponential r => if r.isFinite && r > 0 then .ok () else .error "Exponential rate must be finite and positive."
  | .gamma a s => if a.isFinite && s.isFinite && a > 0 && s > 0 then .ok () else .error "Gamma shape and scale must be finite and positive."
  | .beta a b => if a.isFinite && b.isFinite && a > 0 && b > 0 then .ok () else .error "Beta shapes must be finite and positive."
  | .studentT df nc => if df.isFinite && nc.isFinite && df > 0 then .ok () else .error "Student t needs positive df and finite noncentrality."
  | .chiSquare df nc => if df.isFinite && nc.isFinite && df > 0 && nc >= 0 then .ok () else .error "Chi-square needs positive df and nonnegative noncentrality."
  | .fisherF d1 d2 nc => if d1.isFinite && d2.isFinite && nc.isFinite && d1 > 0 && d2 > 0 && nc >= 0 then .ok () else .error "F needs positive degrees of freedom and nonnegative noncentrality."

/-- Ordinary Lean constructors and bind allow arbitrary records, longitudinal
data, latent variables, missingness, transformations, and user-written estimators. -/
inductive Law : Type → Type 1 where
  /-- Suspend program construction until the interpreter reaches this node. -/
  | defer {α : Type} (next : Unit → Law α) : Law α
  | pure {α : Type} (value : α) : Law α
  | finite {α : Type} (weights : Array (α × Rat)) : Law α
  | binomial (n : Nat) (p : Rat) : Law Nat
  | poisson (rate : Float) : Law Nat
  | scalar (distribution : ScalarLaw) : Law Float
  | bind {α β : Type} (first : Law α) (next : α → Law β) : Law β
  /-- Preserve the model while requiring sampling after a bounded derivation fails. -/
  | samplingOnly {α : Type} (reason : String) (model : Law α) : Law α
  | invalid {α : Type} (message : String) : Law α

instance : Monad Law where
  pure := Law.pure
  bind := Law.bind

def Law.normal (mean sd : Float) : Law Float := .scalar (.normal mean sd)
def Law.uniform (lo hi : Float) : Law Float := .scalar (.uniform lo hi)
def Law.exponential (rate : Float) : Law Float := .scalar (.exponential rate)
def Law.gamma (shape scale : Float) : Law Float := .scalar (.gamma shape scale)
def Law.beta (a b : Float) : Law Float := .scalar (.beta a b)
def Law.studentT (df : Float) (noncentrality : Float := 0) : Law Float := .scalar (.studentT df noncentrality)
def Law.chiSquare (df : Float) (noncentrality : Float := 0) : Law Float := .scalar (.chiSquare df noncentrality)
def Law.fisherF (df1 df2 : Float) (noncentrality : Float := 0) : Law Float := .scalar (.fisherF df1 df2 noncentrality)
def Law.bernoulli (p : Rat) : Law Bool := .finite #[(false, 1-p), (true, p)]
def Law.iid {α : Type} (law : Law α) : Nat → Law (List α)
  | 0 => pure []
  | n+1 => do
    let x ← law
    let xs ← law.iid n
    pure (x :: xs)

/-- Scalar rejection regions can be integrated through CDFs. An arbitrary test
can return an error (e.g. failed model fitting); failures are never discarded. -/
inductive Event : Type → Type 1 where
  | lower (threshold : Float) : Event Float
  | upper (threshold : Float) : Event Float
  | inside (lo hi : Float) : Event Float
  | outside (lo hi : Float) : Event Float
  | test {α : Type} (decide : α → Result Bool) : Event α
  | invalid {α : Type} (message : String) : Event α

def Event.ofBool {α : Type} (decide : α → Bool) : Event α := .test (fun x => .ok (decide x))

def Event.validate {α : Type} : Event α → Result Unit
  | .lower c | .upper c => if c.isFinite then .ok () else .error "Rejection threshold must be finite."
  | .inside a b | .outside a b => if a.isFinite && b.isFinite && a < b then .ok () else .error "Rejection endpoints must be finite and increasing."
  | .test _ => .ok ()
  | .invalid message => .error message

def Event.accepts {α : Type} : Event α → α → Result Bool
  | .lower c, x => do let _ ← finite x; return x <= c
  | .upper c, x => do let _ ← finite x; return x >= c
  | .inside a b, x => do let _ ← finite x; return a <= x && x <= b
  | .outside a b, x => do let _ ← finite x; return x <= a || x >= b
  | .test f, x => f x
  | .invalid message, _ => .error message

private def checkWeights {α : Type} (xs : Array (α × Rat)) : Result Unit := do
  if xs.isEmpty then throw "A finite law must have positive total mass."
  if xs.any (fun x => x.2 < 0) then throw "Probability weights cannot be negative."
  if xs.foldl (fun s x => s+x.2) 0 != 1 then throw "Probability weights must sum exactly to one."

private def checkBinomial (n : Nat) (p : Rat) : Result Unit := do
  if p < 0 || p > 1 then throw "Binomial probability must be between zero and one."
  if n > 1000000 then throw "Binomial trial count exceeds the executable limit (1000000)."

private def checkPoisson (rate : Float) : Result Unit := do
  if !rate.isFinite || rate < 0 then throw "Poisson rate must be finite and nonnegative."
  if rate > 100000 then throw "Poisson rate exceeds the sampler limit (100000)."

/-- Exact binomial probabilities from the ratio of successive mass terms. -/
def binomialWeights (n : Nat) (p : Rat) : Array (Nat × Rat) := Id.run do
  if p == 0 then return #[(0, 1)]
  if p == 1 then return #[(n, 1)]
  let mut mass := (1-p)^n
  let mut out := #[(0, mass)]
  for k in List.range n do
    mass := mass * (n-k : Nat) * p / ((k+1 : Nat) * (1-p))
    out := out.push (k+1, mass)
  return out

inductive EvaluationError where
  | invalid (message : String)
  | unavailable (message : String)
  deriving Repr, DecidableEq

private def validated (r : Result α) : Except EvaluationError α := r.mapError .invalid
private def spend (n : Nat := 1) : StateT Nat (Except EvaluationError) Unit := do
  let remaining ← get
  if n > remaining then throw (.unavailable "Finite enumeration exceeded its work budget.")
  set (remaining-n)

private def enumerateAtom : {α : Type} → Law α → StateT Nat (Except EvaluationError) (Array (α × Rat))
  | _, .samplingOnly reason _ => throw (.unavailable reason)
  | _, .invalid message => throw (.invalid message)
  | _, .pure x => do
    spend
    return #[(x, 1)]
  | _, .finite xs => do
    validated (checkWeights xs)
    spend (xs.size+1)
    return xs.filter (fun x => x.2 > 0)
  | _, .binomial n p => do
    validated (checkBinomial n p)
    spend (n+2)
    return binomialWeights n p
  | _, .poisson rate => do
    validated (checkPoisson rate)
    spend
    if rate == 0 then return #[(0, 1)]
    throw (.unavailable "Poisson has unbounded support.")
  | _, .scalar d => do
    validated d.validate
    throw (.unavailable "Continuous law requires numerical integration or simulation.")
  | _, .bind .. | _, .defer .. =>
    throw (.invalid "Internal error: composite law passed to atomic enumeration.")

/-- Continuations retain the output type and the weight of each pending branch. -/
private inductive EnumerationStack : Type → Type → Type 1 where
  | done : EnumerationStack α α
  | first (next : α → Law β) (rest : EnumerationStack β γ) : EnumerationStack α γ
  | branch (next : α → Law β) (inputs : Array (α × Rat)) (index : Nat)
      (weight : Rat) (output : Array (β × Rat)) (rest : EnumerationStack β γ) :
      EnumerationStack β γ

private inductive EnumerationStep (α : Type) where
  | visit (law : Law α)
  | value (atoms : Array (α × Rat))

-- This loop mirrors the recursive finite interpreter, including its work
-- charges and branch order. Pending binds live on the heap, not the WASM stack.
private def enumerateLoop (fuel : Nat) {α β : Type} (step : EnumerationStep α)
    (stack : EnumerationStack α β) (budget : Nat) :
    Except EvaluationError (Array (β × Rat) × Nat) :=
  match fuel with
  | 0 => .error (.unavailable "Finite enumeration exceeded its work budget.")
  | fuel+1 => do
    match step with
    | .visit (.defer next) =>
      let (_, budget) ← spend.run budget
      enumerateLoop fuel (.visit (next ())) stack budget
    | .visit (.bind first next) =>
      let (_, budget) ← spend.run budget
      enumerateLoop fuel (.visit first) (.first next stack) budget
    | .visit atom =>
      let (xs, budget) ← (enumerateAtom atom).run budget
      enumerateLoop fuel (.value xs) stack budget
    | .value xs =>
      match (dependent := true) (generalizing := true) α, stack, xs with
      | _, .done, xs => return (xs, budget)
      | _, .first next rest, xs =>
        if h : 0 < xs.size then
          let (x, p) := xs[0]
          enumerateLoop fuel (.visit (next x)) (.branch next xs 1 p #[] rest) budget
        else
          enumerateLoop fuel (.value #[]) rest budget
      | _, .branch next inputs index weight output rest, ys =>
        let (_, budget) ← (spend ys.size).run budget
        let mut output := output
        for (y, q) in ys do output := output.push (y, weight*q)
        if h : index < inputs.size then
          let (x, p) := inputs[index]
          enumerateLoop fuel (.visit (next x)) (.branch next inputs (index+1) p output rest) budget
        else
          enumerateLoop fuel (.value output) rest budget

/-- Bounded, stack-safe finite interpreter. It never truncates or renormalizes a law. -/
def enumerate {α : Type} (law : Law α) : StateT Nat (Except EvaluationError) (Array (α × Rat)) :=
  -- Charged visits/branches and their administrative returns take at most
  -- four transitions per work unit. Fuel makes the interpreter total and
  -- kernel-reducible without introducing a fixed structural-depth limit.
  fun budget => enumerateLoop (4*budget+1) (.visit law) .done budget

def exactProbability {α : Type} (law : Law α) (event : Event α)
    (budget : Nat := 10000) : Except EvaluationError Rat := do
  validated event.validate
  let (xs, _) ← (enumerate law).run budget
  let mut probability : Rat := 0
  for (x, p) in xs do
    if ← validated (event.accepts x) then probability := probability + p
  return probability

/-- Convert a probability to binary64 without overflowing large rational numerators. -/
def probabilityFloat (p : Rat) : Float :=
  let numerator := p.num.natAbs
  let ns := numerator.log2 - 52
  let ds := p.den.log2 - 52
  let ratio := nf (numerator >>> ns) / nf (p.den >>> ds)
  let exponent := if ns >= ds then nf (ns-ds) else -(nf (ds-ns))
  ratio * exponent.exp2

structure RandomState where
  seed : UInt64
  remaining : Nat
  deriving Repr

abbrev Draw := StateT RandomState Result

private def tick : Draw Unit := do
  let s ← get
  if s.remaining == 0 then throw "Simulation work budget exhausted; no estimate was returned."
  set { s with remaining := s.remaining-1 }

/-- SplitMix64 with an open (0,1) grid of 52 random bits. Reproducible across WASM runs. -/
def uniform01 : Draw Float := do
  tick
  let s ← get
  let seed := s.seed + 0x9e3779b97f4a7c15
  let z := (seed ^^^ (seed >>> 30)) * 0xbf58476d1ce4e5b9
  let z := (z ^^^ (z >>> 27)) * 0x94d049bb133111eb
  let z := z ^^^ (z >>> 31)
  set { s with seed }
  return ((z >>> 12).toFloat + 0.5) / 4503599627370496

private def standardNormal : Draw Float := do
  let u ← uniform01
  let v ← uniform01
  return (-2*u.log).sqrt * (6.2831853071795864769*v).cos

private def gammaDraw (shape : Float) : Draw Float := do
  let a := if shape < 1 then shape+1 else shape
  let d := a - 1/3
  let c := 1 / (9*d).sqrt
  for _ in [:10000] do
    let x ← standardNormal
    let v0 := 1+c*x
    if v0 > 0 then
      let v := v0*v0*v0
      let u ← uniform01
      if u < 1-0.0331*x*x*x*x || u.log < x*x/2 + d*(1-v+v.log) then
        let value := d*v
        if shape < 1 then return value * ((← uniform01).log/shape).exp
        return value
  throw "Gamma sampler did not converge; no replicates were discarded."

private def poissonDraw (rate : Float) : Draw Nat := do
  checkPoisson rate
  if rate == 0 then return 0
  -- Poisson additivity keeps the product sampler away from underflow.
  let chunks := (rate/20).toUInt64.toNat + 1
  let limit := (-(rate/nf chunks)).exp
  let mut total := 0
  for _ in [:chunks] do
    let mut product := 1.0
    let mut finished := false
    for k in [:10000] do
      product := product * (← uniform01)
      if product <= limit then
        total := total + k
        finished := true
        break
    if !finished then throw "Poisson sampler did not converge; no replicates were discarded."
  return total

private def chiSquareDraw (df nc : Float) : Draw Float := do
  let k ← poissonDraw (nc/2)
  return 2 * (← gammaDraw (df/2 + nf k))

def ScalarLaw.sample (d : ScalarLaw) : Draw Float := do
  d.validate
  let value ← match d with
    | .normal m s => do pure (m + s * (← standardNormal))
    | .uniform a b => do pure (a + (b-a) * (← uniform01))
    | .exponential r => do pure (-(← uniform01).log / r)
    | .gamma a s => do pure (s * (← gammaDraw a))
    | .beta a b => do
      let x ← gammaDraw a
      let y ← gammaDraw b
      pure (x/(x+y))
    | .studentT df nc => do
      let z ← standardNormal
      let v ← gammaDraw (df/2)
      pure ((z+nc)/(2*v/df).sqrt)
    | .chiSquare df nc => chiSquareDraw df nc
    | .fisherF d1 d2 nc => do
      let u ← chiSquareDraw d1 nc
      let v ← chiSquareDraw d2 0
      pure ((u/d1)/(v/d2))
  finite value

private def sampleAtom : {α : Type} → Law α → Draw α
  | _, .invalid message => throw message
  | _, .pure x => do
    tick
    return x
  | _, .finite xs => do
    checkWeights xs
    let u ← uniform01
    let mut mass : Rat := 0
    for (x, p) in xs do
      mass := mass+p
      if u < probabilityFloat mass then return x
    throw "Finite sampler failed to accumulate unit mass."
  | _, .binomial n p => do
    checkBinomial n p
    tick
    if p == 0 then return 0
    if p == 1 then return n
    let threshold := probabilityFloat p
    let mut count := 0
    for _ in [:n] do if (← uniform01) < threshold then count := count+1
    return count
  | _, .poisson rate => do
    tick
    poissonDraw rate
  | _, .scalar d => do
    tick
    d.sample
  | _, .bind .. | _, .defer .. | _, .samplingOnly .. =>
    throw "Internal error: composite law passed to the atomic sampler."

/-- A typed continuation stack keeps nested binds off the browser call stack. -/
private inductive SampleStack : Type → Type → Type 1 where
  | done : SampleStack α α
  | next (continuation : α → Law β) (rest : SampleStack β γ) : SampleStack α γ

private partial def sampleLoop {α β : Type} (law : Law α) (stack : SampleStack α β)
    (state : RandomState) : Result (β × RandomState) := do
  match law with
  | .defer next => sampleLoop (next ()) stack state
  | .samplingOnly _ model => sampleLoop model stack state
  | .bind first next =>
    let (_, state) ← tick.run state
    sampleLoop first (.next next stack) state
  | atom =>
    let (value, state) ← (sampleAtom atom).run state
    match (dependent := true) (generalizing := true) α, stack, value with
    | _, .done, value => return (value, state)
    | _, .next next rest, value => sampleLoop (next value) rest state

def sample {α : Type} (law : Law α) : Draw α :=
  fun state => sampleLoop law .done state

private def poissonMixture (lambda : Float) (cdf : Nat → Result Float) : Result Float := do
  if lambda == 0 then return (← cdf 0)
  if lambda > 300 then throw "Noncentral mixture exceeds the numerical budget."
  let mut weight := (-lambda).exp
  let mut mass := weight
  let mut value := weight * (← cdf 0)
  for k in [1:2001] do
    weight := weight * lambda / nf k
    mass := mass+weight
    value := value + weight * (← cdf k)
    if nf k > lambda && 1-mass < 1e-12 then return clamp01 (← finite value)
  throw "Noncentral mixture did not converge."

private def standardize (x mean sd : Float) : Float :=
  let difference := x-mean
  if difference.isFinite then difference/sd else x/sd-mean/sd

def ScalarLaw.cdf (d : ScalarLaw) (x : Float) : Result Float := do
  d.validate
  let _ ← finite x
  match d with
  | .normal m s => finite (normalCDF (standardize x m s))
  | .uniform a b => finite (clamp01 ((x-a)/(b-a)))
  | .exponential r => pure (if x <= 0 then 0 else 1-(-r*x).exp)
  | .gamma a s => gammaCDF a (x/s)
  | .beta a b => betaCDF a b x
  | .studentT df nc => noncentralTCDF df nc x
  | .chiSquare df nc => poissonMixture (nc/2) (fun k => gammaCDF (df/2+nf k) (x/2))
  | .fisherF d1 d2 nc =>
    if x <= 0 then return 0
    let z := 1 / (1 + (d2.log-d1.log-x.log).exp)
    poissonMixture (nc/2) (fun k => betaCDF (d1/2+nf k) (d2/2) z)

def ScalarLaw.survival (d : ScalarLaw) (x : Float) : Result Float := do
  d.validate
  let _ ← finite x
  match d with
  | .normal m s => finite (normalCDF (-standardize x m s))
  | .studentT df nc => noncentralTCDF df (-nc) (-x)
  | _ => return clamp01 (1 - (← d.cdf x))

def numericalProbability {α : Type} (law : Law α) (event : Event α) : Except EvaluationError Float := do
  validated event.validate
  match (dependent := true) (generalizing := true) α, law, event with
  | _, .scalar d, .lower c =>
    validated d.validate
    (d.cdf c).mapError .unavailable
  | _, .scalar d, .upper c =>
    validated d.validate
    (d.survival c).mapError .unavailable
  | _, .scalar d, .inside a b =>
    validated d.validate
    let lo ← (d.cdf a).mapError .unavailable
    let hi ← (d.cdf b).mapError .unavailable
    return clamp01 (hi-lo)
  | _, .scalar d, .outside a b =>
    validated d.validate
    let lo ← (d.cdf a).mapError .unavailable
    let hi ← (d.survival b).mapError .unavailable
    return clamp01 (lo+hi)
  | _, _, _ => throw (.unavailable "No direct integral is available for this composed law and event.")

end SSC.Probability
