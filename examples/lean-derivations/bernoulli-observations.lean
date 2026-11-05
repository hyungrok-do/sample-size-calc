import SSC
open SSC.Finite SSC.Derivation
open SSC.Calculation (emitEvaluation emitSearch)

set_option maxRecDepth 10000
set_option maxHeartbeats 8000000

namespace ObservedTrial

-- The input is one participant's outcome, not a supplied binomial statistic.
def participant (p : Rat) (hp : 0 ≤ p ∧ p ≤ 1) : Program Bool :=
  .draw (.bernoulli p hp)

def update (successes : Nat) (outcome : Bool) : Nat :=
  successes + if outcome then 1 else 0

def count (p : Rat) (hp : 0 ≤ p ∧ p ≤ 1) (n : Nat) : Program Nat :=
  (participant p hp).iidFold update 0 n

-- A symbolic proof for every sample size, probability and statistic function.
theorem distribution_from_observations (p : Rat) (hp : 0 ≤ p ∧ p ≤ 1)
    (n : Nat) (f : Nat → Rat) :
    (derive (count p hp n)).distribution.expect f =
      ((participant p hp).iid n).expect (fun xs => f (xs.foldl update 0)) := by
  calc
    _ = (count p hp n).expect f := (derive (count p hp n)).correct f
    _ = _ := Program.iidFold_eq_iid _ _ _ _ _

-- Independence and the observation probabilities imply this recurrence.
-- It is a theorem in the input, not a statistical case in the compiler.
theorem count_recurrence (p : Rat) (hp : 0 ≤ p ∧ p ≤ 1) (n : Nat) (f : Nat → Rat) :
    (count p hp (n+1)).expect f =
      (count p hp n).expect (fun k => (1-p)*f k + p*f (k+1)) := by
  simp only [count, Program.iidFold, Program.expect]
  apply Program.expect_congr
  intro k
  simp only [Program.expect, Program.map, participant, Distribution.bernoulli,
    Distribution.expect, expectation, update, Bool.false_eq_true, ↓reduceIte,
    Nat.add_zero, Nat.add_one]
  grind

def nullCount (n : Nat) := count (1/10) (by decide +kernel) n
def alternativeCount (n : Nat) := count (3/10) (by decide +kernel) n
def upper (c : Nat) (k : Nat) : Bool := c ≤ k

-- Calibrate the actual rejection rule from the derived null distribution.
def cutoff (n : Nat) : Nat :=
  let d := (derive (nullCount n)).distribution
  ((List.range (n+2)).find? (fun c => d.probability (upper c) ≤ 1/20)).getD (n+1)

def specification : Specification Nat Nat where
  null := nullCount
  alternative := alternativeCount
  reject := fun n => .test (upper (cutoff n))

-- The generic compiler supplies a proof for this exact probability calculation.
theorem typeI_from_observations (n : Nat) :
    (derive (nullCount n)).distribution.probability (upper (cutoff n)) =
      ((participant (1/10) (by decide +kernel)).iid n).probability
        (fun xs => upper (cutoff n) (xs.foldl update 0)) := by
  rw [Distribution.probability_eq_expect]
  exact distribution_from_observations _ _ _ _

theorem power_from_observations (n : Nat) :
    (derive (alternativeCount n)).distribution.probability (upper (cutoff n)) =
      ((participant (3/10) (by decide +kernel)).iid n).probability
        (fun xs => upper (cutoff n) (xs.foldl update 0)) := by
  rw [Distribution.probability_eq_expect]
  exact distribution_from_observations _ _ _ _

-- These are concrete kernel checks connected to the participant model by the
-- preceding derivation, not merely theorems about a supplied binomial formula.
theorem selected_cutoff : cutoff 25 = 6 := by decide +kernel

theorem target_from_participant_model :
    ((participant (1/10) (by decide +kernel)).iid 25).probability
      (fun xs => upper 6 (xs.foldl update 0)) ≤ 1/20 ∧
    4/5 ≤ ((participant (3/10) (by decide +kernel)).iid 25).probability
      (fun xs => upper 6 (xs.foldl update 0)) := by
  have h0 := typeI_from_observations 25
  have h1 := power_from_observations 25
  rw [selected_cutoff] at h0 h1
  rw [← h0, ← h1]
  decide +kernel

def meetsTarget (n : Nat) : Prop :=
  specification.typeI n ≤ 1/20 ∧ 4/5 ≤ specification.power n

def compiledTarget (n : Nat) : Prop :=
  (derive (nullCount n)).distribution.probability (upper (cutoff n)) ≤ 1/20 ∧
  4/5 ≤ (derive (alternativeCount n)).distribution.probability (upper (cutoff n))

instance (n : Nat) : Decidable (compiledTarget n) := inferInstanceAs (Decidable (_ ∧ _))

theorem compiled_iff_observation_target (n : Nat) : compiledTarget n ↔ meetsTarget n := by
  unfold compiledTarget meetsTarget Specification.typeI Specification.power
  rw [(derive (nullCount n)).probability_eq, (derive (alternativeCount n)).probability_eq]
  rfl

-- A minimum for this specified sequence of calibrated rejection rules.
theorem minimum_sample_size : meetsTarget 25 ∧ ∀ n, n < 25 → ¬ meetsTarget n := by
  constructor
  · apply (compiled_iff_observation_target 25).mp
    decide +kernel
  · have earlier : ∀ n : Fin 25, ¬ compiledTarget n.val := by decide +kernel
    intro n hn h
    exact earlier ⟨n, hn⟩ ((compiled_iff_observation_target n).mpr h)

#eval emitSearch (specification.search ((List.range 40).map (· + 1)).toArray
  { power := 4/5, maxTypeI := some (1/20) })

-- Forced simulation draws all 25 participants from the very same program.
#eval emitEvaluation (specification.evaluate 25
  { route := .simulation, replicates := 10000, seed := 41001 })

end ObservedTrial
