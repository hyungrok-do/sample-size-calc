import SSC
open SSC.Numeric SSC.Probability SSC.Calculation

-- n is the number per group. For independent normal observations with known
-- SD 1, the standardized mean difference has this sampling distribution.
-- The model equations are ordinary Lean source, not a compiler case.
def normalMeans : Specification Nat Float where
  null := fun _ => Law.normal 0 1
  alternative := fun n => Law.normal (0.5 * (nf n / 2).sqrt) 1
  reject := fun _ => match normalQuantile 0.975 with
    | .ok critical => .outside (-critical) critical
    | .error message => .invalid message

-- Ordered candidates are searched without assuming monotone power.
#eval emitSearch (search normalMeans ((List.range 100).map (· + 2)).toArray)

-- The same distributions and rejection event are used in the simulation run.
#eval emitEvaluation (evaluate normalMeans 63
  { route := .simulation, replicates := 10000, seed := 4201 })
