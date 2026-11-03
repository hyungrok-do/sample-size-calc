import Mathlib.Tactic

namespace DesignCompiler

inductive Verdict where
  | adequate | inadequate | unresolved
  deriving DecidableEq, Repr
inductive SearchResult where
  | found (n : Nat) | unresolved (n : Nat) | exhausted
  deriving DecidableEq, Repr

/-- Never skips an unresolved count. Exhaustive search does not require monotonicity. -/
def scan (classify : Nat → Verdict) (start : Nat) : Nat → SearchResult
  | 0 => .exhausted
  | fuel + 1 => match classify start with
    | .adequate => .found start
    | .unresolved => .unresolved start
    | .inadequate => scan classify (start + 1) fuel

theorem scan_found (classify : Nat → Verdict) (fuel start n : Nat)
    (h : scan classify start fuel = .found n) :
    start ≤ n ∧ n < start + fuel ∧ classify n = .adequate ∧
      ∀ k, start ≤ k → k < n → classify k = .inadequate := by
  induction fuel generalizing start with
  | zero => simp [scan] at h
  | succ fuel ih =>
    cases hc : classify start with
    | adequate =>
      simp [scan, hc] at h
      subst n
      exact ⟨le_rfl, by omega, hc, by intro k hk hkn; omega⟩
    | unresolved => simp [scan, hc] at h
    | inadequate =>
      have h' : scan classify (start + 1) fuel = .found n := by simpa [scan, hc] using h
      obtain ⟨hlo, hhi, ha, hp⟩ := ih (start + 1) h'
      refine ⟨by omega, by omega, ha, ?_⟩
      intro k hk hkn
      by_cases he : k = start
      · simpa [he] using hc
      · exact hp k (by omega) hkn

theorem scan_exhausted (classify : Nat → Verdict) (fuel start : Nat)
    (h : scan classify start fuel = .exhausted) :
    ∀ k, start ≤ k → k < start + fuel → classify k = .inadequate := by
  induction fuel generalizing start with
  | zero => intro k h1 h2; omega
  | succ fuel ih =>
    cases hc : classify start with
    | adequate => simp [scan, hc] at h
    | unresolved => simp [scan, hc] at h
    | inadequate =>
      have h' : scan classify (start + 1) fuel = .exhausted := by simpa [scan, hc] using h
      intro k hk hkn
      by_cases he : k = start
      · simpa [he] using hc
      · exact ih (start + 1) h' k (by omega) (by omega)

/-- A conditional certificate theorem. Numerical bounds are hypotheses, not axioms. -/
theorem adequate_of_lower_bound (power lower target : ℝ)
    (bound : lower ≤ power) (adequate : target ≤ lower) : target ≤ power :=
  le_trans adequate bound

theorem inadequate_of_upper_bound (power upper target : ℝ)
    (bound : power ≤ upper) (inadequate : upper < target) : power < target :=
  lt_of_le_of_lt bound inadequate

end DesignCompiler
