import SSC
universe u v
namespace LanguageStress
def polymorphic {α : Type u} (x : α) : α := x
def higherOrder {α : Type u} {β : Type v} (f : α → β) (xs : List α) : List β := xs.map f
inductive Outcome where
  | control
  | treated (response : Int)
  deriving DecidableEq
class Score (α : Type u) where
  score : α → Int
instance : Score Outcome where
  score
    | .control => 0
    | .treated x => x
structure Record (α : Type u) where
  value : α
  count : Nat
def dependent (n : Nat) : Fin (n+1) := ⟨n, Nat.lt_succ_self n⟩
private def hidden := (7 : Nat)
def θ₀ : Rat := 1/3
mutual
  def evenCount : Nat → Bool
    | 0 => true
    | n+1 => oddCount n
  def oddCount : Nat → Bool
    | 0 => false
    | n+1 => evenCount n
end
syntax "duplicate!" term:arg : term
macro_rules
  | `(duplicate! $x) => `(($x, $x))
section
variable {α : Type u}
theorem symbolic (f : α → α) (h : ∀ x, f x = x) : f = id := by
  funext x
  exact h x
end
theorem concrete :
    Score.score (Outcome.treated (-3)) = -3 ∧
    (dependent 12).val = 12 ∧ evenCount 40 = true := by
  decide +kernel
theorem quoted_macro : duplicate! hidden = (7, 7) := rfl
#eval IO.println ("SSC_STRESS_LANGUAGE " ++
  (Lean.toJson [polymorphic 13, (higherOrder (· + 1) [1, 2, 3]).sum,
    (dependent 12).val]).compress)
end LanguageStress
