import SSC

-- A mathematical derivation with symbolic parameters. The power function,
-- admissible domain and approximation bounds remain arbitrary.
-- Supplying validated numerical bounds is a separate obligation.
namespace Bounds

def IsMinimum (domain : Nat → Prop) (power : Nat → Rat) (target : Rat)
    (n : Nat) : Prop :=
  domain n ∧ target ≤ power n ∧
    ∀ m, domain m → m < n → power m < target

theorem minimum_from_bounds
    (domain : Nat → Prop) (power lower upper : Nat → Rat)
    (target : Rat) (n : Nat)
    (hLower : ∀ m, domain m → lower m ≤ power m)
    (hUpper : ∀ m, domain m → power m ≤ upper m)
    (hDomain : domain n)
    (hTarget : target ≤ lower n)
    (hEarlier : ∀ m, domain m → m < n → upper m < target) :
    IsMinimum domain power target n := by
  have reaches : target ≤ power n := by
    have bound := hLower n hDomain
    grind
  have excludesEarlier : ∀ m, domain m → m < n → power m < target := by
    intro m admissible smaller
    have bound := hUpper m admissible
    have below := hEarlier m admissible smaller
    grind
  exact ⟨hDomain, reaches, excludesEarlier⟩

-- Bounds that overlap the target leave a real proof obligation. Neither
-- monotonicity nor a numerical error tolerance has been silently assumed.
theorem target_separated (lo hi target : Rat)
    (hOrder : lo ≤ hi) (hMeets : target ≤ lo) : ¬ hi < target := by
  grind

end Bounds
