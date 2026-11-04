import SSC
open SSC.Finite SSC.Derivation
open SSC.Calculation (emitEvaluation emitSearch)
set_option maxRecDepth 10000
set_option maxHeartbeats 8000000

namespace PairedCategorical

-- Each participant supplies a dependent pair of binary outcomes.
-- The joint law, rather than independence of the two responses, is specified.
def nullPair : Distribution (Bool × Bool) where
  atoms := [((false, false), 7/20), ((false, true), 3/20),
    ((true, false), 3/20), ((true, true), 7/20)]
  nonnegative := by
    intro x hx
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
    rcases hx with rfl | rfl | rfl | rfl <;> decide +kernel
  normalized := by decide +kernel

def alternativePair : Distribution (Bool × Bool) where
  atoms := [((false, false), 1/4), ((false, true), 1/20),
    ((true, false), 7/20), ((true, true), 7/20)]
  nonnegative := by
    intro x hx
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
    rcases hx with rfl | rfl | rfl | rfl <;> decide +kernel
  normalized := by decide +kernel

def difference (pair : Bool × Bool) : Int :=
  (if pair.1 then 1 else 0) - (if pair.2 then 1 else 0)

def update (score : Int) (pair : Bool × Bool) := score + difference pair
def score (law : Distribution (Bool × Bool)) (n : Nat) : Program Int :=
  (Program.draw law).iidFold update 0 n

theorem score_distribution_from_pairs (law : Distribution (Bool × Bool))
    (n : Nat) (f : Int → Rat) :
    (derive (score law n)).distribution.expect f =
      ((Program.draw law).iid n).expect (fun pairs => f (pairs.foldl update 0)) := by
  calc
    _ = (score law n).expect f := (derive (score law n)).correct f
    _ = _ := Program.iidFold_eq_iid _ _ _ _ _

def upper (cut : Nat) (s : Int) : Bool := (cut : Int) ≤ s
def calibrate (n : Nat) (tailProbability : Nat → Rat) : Nat :=
  ((List.range (n+2)).find? (fun c => tailProbability c ≤ 1/20)).getD (n+1)

def cutoff (n : Nat) : Nat :=
  calibrate n (fun c => (derive (score nullPair n)).distribution.probability (upper c))

def specification : Specification Nat Int where
  null := score nullPair
  alternative := score alternativePair
  reject := fun n => .test (upper (cutoff n))

-- Intermediate tables are proposed by evaluation and checked one step at a time.
-- The resulting certificate still proves the complete observation-model claim.
def checkedNull : Certificate (score nullPair 23) := by finite_derive
def checkedAlternative : Certificate (score alternativePair 23) := by finite_derive

theorem calibration_from_certificate (n : Nat) (certificate : Certificate (score nullPair n)) :
    cutoff n = calibrate n (fun c => certificate.distribution.probability (upper c)) := by
  have tails :
      (fun c => (derive (score nullPair n)).distribution.probability (upper c)) =
      (fun c => certificate.distribution.probability (upper c)) := by
    funext c
    exact ((derive (score nullPair n)).probability_eq (upper c)).trans
      (certificate.probability_eq (upper c)).symm
  exact congrArg (calibrate n) tails

theorem selected_cutoff : cutoff 23 = 5 :=
  (calibration_from_certificate 23 checkedNull).trans (by decide +kernel)

theorem target_from_joint_observations :
    ((Program.draw nullPair).iid 23).probability
      (fun pairs => upper 5 (pairs.foldl update 0)) ≤ 1/20 ∧
    4/5 ≤ ((Program.draw alternativePair).iid 23).probability
      (fun pairs => upper 5 (pairs.foldl update 0)) := by
  have h0 : checkedNull.distribution.probability (upper 5) =
      ((Program.draw nullPair).iid 23).probability
        (fun pairs => upper 5 (pairs.foldl update 0)) := by
    calc
      _ = (score nullPair 23).probability (upper 5) := checkedNull.probability_eq _
      _ = _ := Program.iidFold_eq_iid _ _ _ _ _
  have h1 : checkedAlternative.distribution.probability (upper 5) =
      ((Program.draw alternativePair).iid 23).probability
        (fun pairs => upper 5 (pairs.foldl update 0)) := by
    calc
      _ = (score alternativePair 23).probability (upper 5) := checkedAlternative.probability_eq _
      _ = _ := Program.iidFold_eq_iid _ _ _ _ _
  rw [← h0, ← h1]
  decide +kernel

#eval emitSearch (specification.search ((List.range 35).map (· + 1)).toArray
  { power := 4/5, maxTypeI := some (1/20) })

-- Deliberately small compilation budget: simulate the original pairs and score.
#eval emitEvaluation (specification.evaluate 23
  { replicates := 2000, seed := 41002 } (budget := 1))

end PairedCategorical
