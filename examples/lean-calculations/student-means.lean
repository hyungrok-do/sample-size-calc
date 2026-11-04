import SSC
open SSC.Numeric SSC.Probability SSC.Calculation

-- Equal-sized independent normal groups, common unknown variance.
-- The pooled t statistic has df=2n-2 and noncentrality delta*sqrt(n/2)/SD.
-- Changing delta or SD changes these equations, without changing the compiler.
def studentMeans (delta sd : Float) : Specification Nat Float where
  null := fun n => Law.studentT (2 * nf n - 2)
  alternative := fun n =>
    if !sd.isFinite || sd <= 0 || !delta.isFinite then .invalid "SD must be positive and the effect finite."
    else Law.studentT (2 * nf n - 2) (delta * (nf n/2).sqrt / sd)
  reject := fun n =>
    match tQuantile (2 * nf n - 2) 0.975 with
    | .ok c => .outside (-c) c
    | .error message => .invalid message

#eval emitSearch (search (studentMeans 0.5 1)
  ((List.range 90).map (· + 2)).toArray)

-- A design can include a detectable effect as well as sample size.
-- Search is over the supplied effect grid, with its order and resolution explicit.
def effects : Specification Float Float where
  null := fun _ => (studentMeans 0 1).null 64
  alternative := fun delta => (studentMeans delta 1).alternative 64
  reject := fun _ => (studentMeans 0 1).reject 64

#eval emitSearch (search effects
  ((List.range 31).map (fun i => 0.30 + nf i/100)).toArray)
