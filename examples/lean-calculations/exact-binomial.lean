import SSC
open SSC.Numeric SSC.Probability SSC.Calculation
set_option maxRecDepth 4096
set_option maxHeartbeats 8000000

-- One-sided rejection region is calibrated from the supplied null distribution.
def criticalCount (n : Nat) (p0 alpha : Rat) : Nat := Id.run do
  let weights := binomialWeights n p0
  let mut tail : Rat := 1
  for k in [:n+1] do
    if tail <= alpha then return k
    tail := tail - (weights.find? (fun item => item.1 == k) |>.map Prod.snd |>.getD 0)
  return n+1

def binomialTest (p0 p1 alpha : Rat) : Specification Nat Nat where
  null := fun n => .binomial n p0
  alternative := fun n => .binomial n p1
  reject := fun n =>
    let cutoff := criticalCount n p0 alpha
    .ofBool (fun successes => successes >= cutoff)

def trial := binomialTest (1/10) (3/10) (1/20)
#eval emitSearch (search trial ((List.range 40).map (·+1)).toArray
  { power := 4/5, maxTypeI := some (1/20) })

-- This is a kernel-checked claim about an exact evaluation of the program.
-- The general probability library is not itself a formalized probability space.
theorem coin_probability :
    (exactProbability (.binomial 3 (1/2)) (.ofBool (fun k => k == 3))).toOption = some (1/8) := by
  decide +kernel

-- Force the second interpreter on precisely the same trial specification.
#eval emitEvaluation (evaluate trial 25 { route := .simulation, replicates := 10000, seed := 4202 })
