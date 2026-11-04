import SSC

-- Independent Bernoulli outcomes: H₀ p = 1/10, H₁ p = 3/10.
-- Reject when at least c of n outcomes are successes. The equations below
-- specify the model and the rejection rule; no analysis selector is involved.
set_option maxRecDepth 4096
set_option maxHeartbeats 8000000

namespace Trial

-- Integer probability weights, with common denominator b^n.
-- Appending a failure or success gives the next row by convolution.
def weights (a b : Nat) : Nat → List Nat
  | 0 => [1]
  | n + 1 =>
    let previous := weights a b n
    List.zipWith (fun failure success => (b - a) * failure + a * success)
      (previous ++ [0]) (0 :: previous)

def rejectionProbability (a b n c : Nat) : Rat :=
  ((weights a b n).drop c).sum / (b ^ n : Nat)

-- Exact integer comparisons for alpha ≤ 1/20 and power ≥ 4/5.
-- All n+2 cutoffs, including always/never rejecting, are considered.
def adequate (n : Nat) : Prop :=
  let nullWeights := weights 1 10 n
  let alternativeWeights := weights 3 10 n
  0 < n ∧ (List.range (n + 2)).any (fun c =>
    20 * (nullWeights.drop c).sum ≤ 10 ^ n &&
    4 * 10 ^ n ≤ 5 * (alternativeWeights.drop c).sum) = true

instance (n : Nat) : Decidable (adequate n) := inferInstanceAs (Decidable (_ ∧ _))

-- This generic search constructs the witness AND a minimality proof.
def search := SSC.leastBelow adequate 40
#eval search.value?

-- A concrete rejection rule attaining the criterion at the minimum.
theorem rejection_rule_meets_target :
    rejectionProbability 1 10 25 6 ≤ 1/20 ∧
    4/5 ≤ rejectionProbability 3 10 25 6 := by decide +kernel

theorem selected_value : search.value? = some 25 := by decide +kernel

-- The proof is assembled from the computed witness and the general theorem
-- about proof-carrying search. There is no binomial solver in the compiler.
theorem minimum_sample_size :
    adequate 25 ∧ ∀ m, m < 25 → ¬ adequate m := by
  have certificate := search.sound
  rw [selected_value] at certificate
  exact ⟨certificate.2.1, certificate.2.2⟩

#eval rejectionProbability 1 10 25 6
#eval rejectionProbability 3 10 25 6

end Trial
