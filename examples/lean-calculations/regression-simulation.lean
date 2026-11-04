import SSC
open SSC.Numeric SSC.Probability SSC.Calculation
set_option maxRecDepth 10000

-- Random covariates and heteroscedastic errors. The estimator below uses the
-- HC1 sandwich variance. It is applied identically under H0 and H1.
def observation (slope : Float) : Law (Float × Float) := do
  let x ← Law.uniform (-1) 1
  let e ← Law.normal 0 1
  pure (x, slope*x + (0.5 + x.abs)*e)

def robustSlope (critical : Float) (xs : List (Float × Float)) : Result Bool := do
  let n := nf xs.length
  if xs.length < 3 then throw "Regression needs at least three observations."
  let mx := xs.foldl (fun s xy => s+xy.1) 0 / n
  let my := xs.foldl (fun s xy => s+xy.2) 0 / n
  let sxx := xs.foldl (fun s xy => s+(xy.1-mx)*(xy.1-mx)) 0
  if sxx <= 0 then throw "Singular regression design."
  let b := xs.foldl (fun s xy => s+(xy.1-mx)*(xy.2-my)) 0 / sxx
  let meat := xs.foldl (fun s xy =>
    let residual := xy.2-my-b*(xy.1-mx)
    s+(xy.1-mx)*(xy.1-mx)*residual*residual) 0
  let variance := n/(n-2)*meat/(sxx*sxx)
  if !variance.isFinite || variance <= 0 then throw "Invalid sandwich variance."
  return (b/variance.sqrt).abs >= critical

def regression : Specification Nat (List (Float × Float)) where
  null := fun n => (observation 0).iid n
  alternative := fun n => (observation 0.8).iid n
  reject := fun _ => match normalQuantile 0.975 with
    | .ok critical => .test (robustSlope critical)
    | .error message => .invalid message

#eval emitSearch (search regression #[30, 60, 90, 120]
  { power := 4/5 } { replicates := 1500, seed := 4204, maxDraws := 10000000 })
