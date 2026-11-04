import SSC
open SSC.Numeric SSC.Probability SSC.Calculation
set_option maxRecDepth 10000

-- Correlated lognormal paired observations. The actual differences are tested.
-- The distribution of their studentized mean is not declared to be Student t.
def difference (shift rho : Float) : Law Float := do
  let shared ← Law.normal 0 1
  let residual ← Law.normal 0 1
  let before := shared.exp
  let after := (shift + rho*shared + (1-rho*rho).sqrt*residual).exp
  pure (after-before)

def pairedTest (critical : Float) (xs : List Float) : Result Bool := do
  let n := nf xs.length
  if xs.length < 2 then throw "At least two pairs are required."
  let mean := xs.foldl (· + ·) 0 / n
  let variance := xs.foldl (fun s x => s + (x-mean)*(x-mean)) 0 / (n-1)
  if !variance.isFinite || variance <= 0 then throw "Paired variance is invalid; no replicate is discarded."
  return (mean / (variance/n).sqrt).abs >= critical

def paired : Specification Nat (List Float) where
  null := fun n => (difference 0 0.65).iid n
  alternative := fun n => (difference 0.5 0.65).iid n
  reject := fun n =>
    match tQuantile (nf n-1) 0.975 with
    | .ok c => .test (pairedTest c)
    | .error message => .invalid message

-- Automatic fallback samples these same paired outcomes and applies pairedTest.
-- Fixed-budget intervals can leave candidates unresolved near 80% power.
-- The null rejection rate is reported; it is not assumed to equal 5%.
#eval emitSearch (search paired #[12, 24, 36, 48, 60]
  { power := 4/5 } { replicates := 2000, seed := 4203, maxDraws := 10000000 })
