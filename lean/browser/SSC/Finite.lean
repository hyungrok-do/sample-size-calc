import Lean

namespace SSC.Finite

abbrev Table (α : Type) := List (α × Rat)

def expectation : Table α → (α → Rat) → Rat
  | [], _ => 0
  | (x, p) :: xs, f => p * f x + expectation xs f

def Nonnegative (xs : Table α) : Prop := ∀ x ∈ xs, 0 ≤ x.2

theorem expectation_congr (xs : Table α) (f g : α → Rat)
    (h : ∀ x, f x = g x) : expectation xs f = expectation xs g := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp only [expectation, h, ih]

theorem expectation_add (xs : Table α) (f g : α → Rat) :
    expectation xs (fun x => f x + g x) = expectation xs f + expectation xs g := by
  induction xs with
  | nil => simp only [expectation]; grind
  | cons x xs ih => simp only [expectation, ih]; grind

theorem expectation_scale (xs : Table α) (c : Rat) (f : α → Rat) :
    expectation xs (fun x => c * f x) = c * expectation xs f := by
  induction xs with
  | nil => simp [expectation]
  | cons x xs ih => simp only [expectation, ih]; grind

theorem expectation_nonnegative (xs : Table α) (hx : Nonnegative xs)
    (f : α → Rat) (hf : ∀ x, 0 ≤ f x) : 0 ≤ expectation xs f := by
  induction xs with
  | nil => change 0 ≤ (0 : Rat); decide
  | cons x xs ih =>
    exact Rat.add_nonneg (Rat.mul_nonneg (hx x (by simp)) (hf x.1))
      (ih (fun y hy => hx y (by simp [hy])))

theorem expectation_append (xs ys : Table α) (f : α → Rat) :
    expectation (xs ++ ys) f = expectation xs f + expectation ys f := by
  induction xs with
  | nil => simp only [List.nil_append, expectation]; grind
  | cons x xs ih => simp only [List.cons_append, expectation, ih]; grind

/-- Balanced exact summation avoids retaining one kernel reduction frame per
atom while normalizing large rational weights. Zero fuel remains exact. -/
def expectationBalanced : Nat → Table α → (α → Rat) → Rat
  | 0, xs, f => expectation xs f
  | fuel+1, xs, f =>
    if xs.length ≤ 8 then expectation xs f
    else
      let mid := xs.length / 2
      expectationBalanced fuel (xs.take mid) f +
        expectationBalanced fuel (xs.drop mid) f

theorem expectationBalanced_eq (fuel : Nat) (xs : Table α) (f : α → Rat) :
    expectationBalanced fuel xs f = expectation xs f := by
  induction fuel generalizing xs with
  | zero => rfl
  | succ fuel ih =>
    simp only [expectationBalanced, ih]
    split
    · rfl
    · rw [← expectation_append, List.take_append_drop]

/-- Boolean events need only selected masses, not multiplication by 0 and 1. -/
def eventMass : Table α → (α → Bool) → Rat
  | [], _ => 0
  | (x, p) :: xs, event =>
    if event x then p + eventMass xs event else eventMass xs event

theorem eventMass_eq (xs : Table α) (event : α → Bool) :
    eventMass xs event = expectation xs (fun x => if event x then 1 else 0) := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [eventMass, expectation, ih]
    cases h : event x.1 <;> simp [h, Rat.mul_one, Rat.mul_zero, Rat.zero_add]

def eventMassBalanced : Nat → Table α → (α → Bool) → Rat
  | 0, xs, event => eventMass xs event
  | fuel+1, xs, event =>
    if xs.length ≤ 8 then eventMass xs event
    else
      let mid := xs.length / 2
      eventMassBalanced fuel (xs.take mid) event +
        eventMassBalanced fuel (xs.drop mid) event

theorem eventMassBalanced_eq (fuel : Nat) (xs : Table α) (event : α → Bool) :
    eventMassBalanced fuel xs event = expectation xs (fun x => if event x then 1 else 0) := by
  induction fuel generalizing xs with
  | zero => exact eventMass_eq _ _
  | succ fuel ih =>
    simp only [eventMassBalanced, ih, eventMass_eq]
    split
    · rfl
    · rw [← expectation_append, List.take_append_drop]

/-- Generic algebraic identities avoid rational GCD/division work at pure
program nodes and zero-weight branches. Each branch has an equality proof. -/
def multiplyWeight (p q : Rat) : Rat :=
  if q = 1 then p else if p = 1 then q else if q = 0 then 0 else if p = 0 then 0 else p*q

@[simp] theorem multiplyWeight_eq (p q : Rat) : multiplyWeight p q = p*q := by
  unfold multiplyWeight
  split
  · rename_i h; subst q; exact (Rat.mul_one p).symm
  · split
    · rename_i h; subst p; exact (Rat.one_mul q).symm
    · split
      · rename_i h; subst q; exact (Rat.mul_zero p).symm
      · split
        · rename_i h; subst p; exact (Rat.zero_mul q).symm
        · rfl

def scale (p : Rat) (xs : Table α) : Table α := xs.map (fun (x, q) => (x, multiplyWeight p q))

theorem expectation_scale_weights (p : Rat) (xs : Table α) (f : α → Rat) :
    expectation (scale p xs) f = p * expectation xs f := by
  induction xs with
  | nil => simp [scale, expectation]
  | cons x xs ih => simp only [scale, List.map_cons, expectation, multiplyWeight_eq] at *; grind

theorem scale_nonnegative (p : Rat) (hp : 0 ≤ p) (xs : Table α)
    (hx : Nonnegative xs) : Nonnegative (scale p xs) := by
  intro y hy
  obtain ⟨x, hmem, rfl⟩ := List.mem_map.mp hy
  simpa only [multiplyWeight_eq] using Rat.mul_nonneg hp (hx x hmem)

def flatten (xs : Table α) (k : α → Table β) : Table β :=
  xs.flatMap (fun (x, p) => scale p (k x))

theorem expectation_flatten (xs : Table α) (k : α → Table β) (f : β → Rat) :
    expectation (flatten xs k) f = expectation xs (fun x => expectation (k x) f) := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [flatten, List.flatMap_cons, expectation_append,
      expectation_scale_weights, expectation] at *
    rw [ih]

theorem flatten_nonnegative (xs : Table α) (hx : Nonnegative xs)
    (k : α → Table β) (hk : ∀ x, Nonnegative (k x)) : Nonnegative (flatten xs k) := by
  intro y hy
  obtain ⟨x, hmem, hy⟩ := List.mem_flatMap.mp hy
  exact scale_nonnegative x.2 (hx x hmem) (k x.1) (hk x.1) y hy

/-- Insert a mass at its exact outcome, combining equal outcomes. -/
def insert [DecidableEq α] (x : α) (p : Rat) : Table α → Table α
  | [] => [(x, p)]
  | (y, q) :: ys => if x = y then (y, p+q) :: ys else (y, q) :: insert x p ys

theorem expectation_insert [DecidableEq α] (x : α) (p : Rat)
    (xs : Table α) (f : α → Rat) :
    expectation (insert x p xs) f = p * f x + expectation xs f := by
  induction xs with
  | nil => simp [insert, expectation]
  | cons y xs ih =>
    simp only [insert]
    split
    · rename_i h; simp only [expectation, h]; grind
    · simp only [expectation, ih]; grind

theorem insert_nonnegative [DecidableEq α] (x : α) (p : Rat) (hp : 0 ≤ p)
    (xs : Table α) (hx : Nonnegative xs) : Nonnegative (insert x p xs) := by
  induction xs with
  | nil => intro y hy; simp only [insert, List.mem_singleton] at hy; subst y; exact hp
  | cons y xs ih =>
    have hy := hx y (by simp)
    have ht : Nonnegative xs := fun z hz => hx z (by simp [hz])
    simp only [insert]
    split
    · intro z hz
      simp only [List.mem_cons] at hz
      rcases hz with rfl | hz
      · exact Rat.add_nonneg hp hy
      · exact ht z hz
    · intro z hz
      simp only [List.mem_cons] at hz
      rcases hz with rfl | hz
      · exact hy
      · exact ih ht z hz

def collect [DecidableEq α] : Table α → Table α
  | [] => []
  | (x, p) :: xs => insert x p (collect xs)

theorem expectation_collect [DecidableEq α] (xs : Table α) (f : α → Rat) :
    expectation (collect xs) f = expectation xs f := by
  induction xs with
  | nil => rfl
  | cons x xs ih => simp only [collect, expectation_insert, expectation, ih]

theorem collect_nonnegative [DecidableEq α] (xs : Table α) (hx : Nonnegative xs) :
    Nonnegative (collect xs) := by
  induction xs with
  | nil => intro x hx; contradiction
  | cons x xs ih =>
    exact insert_nonnegative x.1 x.2 (hx x (by simp)) _
      (ih (fun y hy => hx y (by simp [hy])))

/-- Add a new table into an already collected table, preserving its distinct states. -/
def merge [DecidableEq α] : Table α → Table α → Table α
  | [], ys => ys
  | (x, p) :: xs, ys => insert x p (merge xs ys)

theorem expectation_merge [DecidableEq α] (xs ys : Table α) (f : α → Rat) :
    expectation (merge xs ys) f = expectation xs f + expectation ys f := by
  induction xs with
  | nil => simp only [merge, expectation]; grind
  | cons x xs ih => simp only [merge, expectation_insert, expectation, ih]; grind

theorem merge_nonnegative [DecidableEq α] (xs ys : Table α)
    (hx : Nonnegative xs) (hy : Nonnegative ys) : Nonnegative (merge xs ys) := by
  induction xs with
  | nil => exact hy
  | cons x xs ih =>
    exact insert_nonnegative x.1 x.2 (hx x (by simp)) _
      (ih (fun y hmem => hx y (by simp [hmem])))

/-- Combine tables in a balanced tree. The fallback at zero fuel is exact too.
Balancing limits the nesting of table operations during kernel reduction. -/
def mergeMany [DecidableEq α] : Nat → List (Table α) → Table α
  | 0, tables => tables.foldr merge []
  | _+1, [] => []
  | _+1, [table] => table
  | fuel+1, tables =>
    let mid := tables.length / 2
    merge (mergeMany fuel (tables.take mid)) (mergeMany fuel (tables.drop mid))

theorem expectation_mergeMany [DecidableEq α] (fuel : Nat)
    (tables : List (Table α)) (f : α → Rat) :
    expectation (mergeMany fuel tables) f = (tables.map (fun t => expectation t f)).sum := by
  induction fuel generalizing tables with
  | zero =>
    induction tables with
    | nil => rfl
    | cons t ts ih =>
      simp only [mergeMany, List.foldr_cons, expectation_merge, List.map_cons, List.sum_cons] at *
      rw [ih]
  | succ fuel ih =>
    cases tables with
    | nil => rfl
    | cons t ts =>
      cases ts with
      | nil => simp only [mergeMany, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]; grind
      | cons u us =>
        simp only [mergeMany, expectation_merge, ih]
        rw [← List.sum_append, ← List.map_append, List.take_append_drop]

theorem mergeMany_nonnegative [DecidableEq α] (fuel : Nat)
    (tables : List (Table α)) (h : ∀ t ∈ tables, Nonnegative t) :
    Nonnegative (mergeMany fuel tables) := by
  induction fuel generalizing tables with
  | zero =>
    induction tables with
    | nil => intro x hx; contradiction
    | cons t ts ih =>
      exact merge_nonnegative _ _ (h t (by simp))
        (ih (fun x hx => h x (by simp [hx])))
  | succ fuel ih =>
    cases tables with
    | nil => intro x hx; contradiction
    | cons t ts =>
      cases ts with
      | nil => exact h t (by simp)
      | cons u us =>
        exact merge_nonnegative _ _
          (ih _ (fun x hx => h x (List.mem_of_mem_take hx)))
          (ih _ (fun x hx => h x (List.mem_of_mem_drop hx)))

/-- Compose branch masses using balanced exact merging. -/
def flattenMerged [DecidableEq β] (xs : Table α) (k : α → Table β) : Table β :=
  mergeMany xs.length (xs.map (fun (x, p) => scale p (k x)))

theorem expectation_flattenMerged [DecidableEq β] (xs : Table α)
    (k : α → Table β) (f : β → Rat) :
    expectation (flattenMerged xs k) f = expectation xs (fun x => expectation (k x) f) := by
  simp only [flattenMerged, expectation_mergeMany, List.map_map]
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [List.map_cons, List.sum_cons, Function.comp_apply,
      expectation_scale_weights, expectation, ih]

theorem flattenMerged_nonnegative [DecidableEq β] (xs : Table α)
    (hx : Nonnegative xs) (k : α → Table β) (hk : ∀ x, Nonnegative (k x)) :
    Nonnegative (flattenMerged xs k) := by
  apply mergeMany_nonnegative
  intro t ht
  obtain ⟨x, hmem, rfl⟩ := List.mem_map.mp ht
  exact scale_nonnegative x.2 (hx x hmem) _ (hk x.1)

namespace Ordered
variable {α : Type}

/-- Ordered stream merging avoids scanning the whole accumulated table for
 each incoming state. Equality is checked independently of `compare`, so even an
 arbitrary user ordering preserves every expectation and nonnegative weight. -/
def mergeHead [DecidableEq α] [Ord α] (x : α × Rat) (xs : Table α)
    (rest : Table α → Table α) : Table α → Table α
  | [] => x :: xs
  | y :: ys =>
    if x.1 = y.1 then (x.1, x.2 + y.2) :: rest ys
    else if compare x.1 y.1 == .lt then x :: rest (y :: ys)
    else y :: mergeHead x xs rest ys

def merge [DecidableEq α] [Ord α] : Table α → Table α → Table α
  | [], ys => ys
  | x :: xs, ys => mergeHead x xs (merge xs) ys

theorem expectation_mergeHead [DecidableEq α] [Ord α] (x : α × Rat) (xs : Table α)
    (rest : Table α → Table α) (f : α → Rat)
    (hrest : ∀ ys, expectation (rest ys) f = expectation xs f + expectation ys f)
    (ys : Table α) : expectation (mergeHead x xs rest ys) f =
      expectation (x :: xs) f + expectation ys f := by
  induction ys with
  | nil => simp [mergeHead, expectation, Rat.add_zero]
  | cons y ys ih =>
    simp only [mergeHead]
    split
    · rename_i h
      simp only [expectation, hrest, h]
      grind
    · split
      · simp only [expectation, hrest]
        grind
      · simp only [expectation, ih]
        grind

theorem expectation_merge [DecidableEq α] [Ord α] (xs ys : Table α) (f : α → Rat) :
    expectation (merge xs ys) f = expectation xs f + expectation ys f := by
  induction xs generalizing ys with
  | nil => simp [merge, expectation, Rat.zero_add]
  | cons x xs ih => exact expectation_mergeHead x xs (merge xs) f ih ys

theorem nonnegative_cons (x : α × Rat) (xs : Table α)
    (hx : 0 ≤ x.2) (hxs : Nonnegative xs) : Nonnegative (x :: xs) := by
  intro y hy
  rcases List.mem_cons.mp hy with rfl | hy
  · exact hx
  · exact hxs y hy

theorem mergeHead_nonnegative [DecidableEq α] [Ord α] (x : α × Rat) (xs : Table α)
    (rest : Table α → Table α) (hx : 0 ≤ x.2) (hxs : Nonnegative xs)
    (hrest : ∀ ys, Nonnegative ys → Nonnegative (rest ys))
    (ys : Table α) (hys : Nonnegative ys) : Nonnegative (mergeHead x xs rest ys) := by
  induction ys with
  | nil => exact nonnegative_cons x xs hx hxs
  | cons y ys ih =>
    have hy := hys y (by simp)
    have ht : Nonnegative ys := fun z hz => hys z (by simp [hz])
    simp only [mergeHead]
    split
    · exact nonnegative_cons _ _ (Rat.add_nonneg hx hy) (hrest ys ht)
    · split
      · exact nonnegative_cons x _ hx (hrest (y :: ys) hys)
      · exact nonnegative_cons y _ hy (ih ht)

theorem merge_nonnegative [DecidableEq α] [Ord α] (xs ys : Table α)
    (hx : Nonnegative xs) (hy : Nonnegative ys) : Nonnegative (merge xs ys) := by
  induction xs generalizing ys with
  | nil => exact hy
  | cons x xs ih =>
    exact mergeHead_nonnegative x xs (merge xs) (hx x (by simp))
      (fun z hz => hx z (by simp [hz]))
      (fun ys hy => ih ys (fun z hz => hx z (by simp [hz])) hy) ys hy

def mergeMany [DecidableEq α] [Ord α] : Nat → List (Table α) → Table α
  | 0, tables => tables.foldr merge []
  | _+1, [] => []
  | _+1, [table] => table
  | fuel+1, tables =>
    let mid := tables.length / 2
    merge (mergeMany fuel (tables.take mid)) (mergeMany fuel (tables.drop mid))

theorem expectation_mergeMany [DecidableEq α] [Ord α] (fuel : Nat)
    (tables : List (Table α)) (f : α → Rat) :
    expectation (mergeMany fuel tables) f = (tables.map (fun t => expectation t f)).sum := by
  induction fuel generalizing tables with
  | zero =>
    induction tables with
    | nil => rfl
    | cons t ts ih =>
      simp only [mergeMany, List.foldr_cons, expectation_merge, List.map_cons, List.sum_cons] at *
      rw [ih]
  | succ fuel ih =>
    cases tables with
    | nil => rfl
    | cons t ts =>
      cases ts with
      | nil => simp only [mergeMany, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]; grind
      | cons u us =>
        simp only [mergeMany, expectation_merge, ih]
        rw [← List.sum_append, ← List.map_append, List.take_append_drop]

theorem mergeMany_nonnegative [DecidableEq α] [Ord α] (fuel : Nat)
    (tables : List (Table α)) (h : ∀ t ∈ tables, Nonnegative t) :
    Nonnegative (mergeMany fuel tables) := by
  induction fuel generalizing tables with
  | zero =>
    induction tables with
    | nil => intro x hx; contradiction
    | cons t ts ih =>
      exact merge_nonnegative _ _ (h t (by simp))
        (ih (fun x hx => h x (by simp [hx])))
  | succ fuel ih =>
    cases tables with
    | nil => intro x hx; contradiction
    | cons t ts =>
      cases ts with
      | nil => exact h t (by simp)
      | cons u us =>
        exact merge_nonnegative _ _
          (ih _ (fun x hx => h x (List.mem_of_mem_take hx)))
          (ih _ (fun x hx => h x (List.mem_of_mem_drop hx)))

/-- Arrange and combine outcomes using balanced stream merges. -/
def sort [DecidableEq α] [Ord α] (xs : Table α) : Table α :=
  mergeMany 64 (xs.map (fun x => [x]))

theorem expectation_sort [DecidableEq α] [Ord α] (xs : Table α) (f : α → Rat) :
    expectation (sort xs) f = expectation xs f := by
  rw [sort, expectation_mergeMany]
  induction xs with
  | nil => rfl
  | cons x xs ih => simp only [List.map_cons, List.sum_cons, expectation, Rat.add_zero, ih]

theorem sort_nonnegative [DecidableEq α] [Ord α] (xs : Table α)
    (hx : Nonnegative xs) : Nonnegative (sort xs) := by
  apply mergeMany_nonnegative
  intro table ht
  obtain ⟨x, hmem, rfl⟩ := List.mem_map.mp ht
  exact nonnegative_cons x [] (hx x hmem) (by intro y hy; contradiction)

def flattenMerged [DecidableEq β] [Ord β] (xs : Table α) (k : α → Table β) : Table β :=
  mergeMany xs.length (xs.map (fun (x, p) => scale p (sort (k x))))

theorem expectation_flattenMerged [DecidableEq β] [Ord β] (xs : Table α)
    (k : α → Table β) (f : β → Rat) :
    expectation (flattenMerged xs k) f = expectation xs (fun x => expectation (k x) f) := by
  simp only [flattenMerged, expectation_mergeMany, List.map_map]
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [List.map_cons, List.sum_cons, Function.comp_apply,
      expectation_scale_weights, expectation_sort, expectation, ih]

theorem flattenMerged_nonnegative [DecidableEq β] [Ord β] (xs : Table α)
    (hx : Nonnegative xs) (k : α → Table β) (hk : ∀ x, Nonnegative (k x)) :
    Nonnegative (flattenMerged xs k) := by
  apply mergeMany_nonnegative
  intro t ht
  obtain ⟨x, hmem, rfl⟩ := List.mem_map.mp ht
  exact scale_nonnegative x.2 (hx x hmem) _ (sort_nonnegative _ (hk x.1))

end Ordered

/-- A finite probability distribution, including proofs of positivity and mass one. -/
structure Distribution (α : Type) where
  atoms : Table α
  nonnegative : Nonnegative atoms
  normalized : expectation atoms (fun _ => 1) = 1

namespace Distribution

def expect (d : Distribution α) (f : α → Rat) := expectation d.atoms f
def probability (d : Distribution α) (event : α → Bool) : Rat :=
  -- The depth bound controls summation shape only; zero fuel remains exact.
  eventMassBalanced 64 d.atoms event

theorem probability_eq_expect (d : Distribution α) (event : α → Bool) :
    d.probability event = d.expect (fun x => if event x then 1 else 0) :=
  eventMassBalanced_eq _ _ _

theorem expect_const (d : Distribution α) (c : Rat) : d.expect (fun _ => c) = c := by
  have h := expectation_scale d.atoms c (fun _ => 1)
  simp only [Rat.mul_one, d.normalized] at h
  exact h

theorem expect_nonnegative (d : Distribution α) (f : α → Rat)
    (hf : ∀ x, 0 ≤ f x) : 0 ≤ d.expect f := expectation_nonnegative _ d.nonnegative _ hf

theorem probability_bounds (d : Distribution α) (event : α → Bool) :
    0 ≤ d.probability event ∧ d.probability event ≤ 1 := by
  rw [probability_eq_expect]
  have h0 := d.expect_nonnegative (fun x => if event x then 1 else 0) (by intro x; split <;> decide)
  have h1 := d.expect_nonnegative (fun x => if event x then 0 else 1) (by intro x; split <;> decide)
  have hs := expectation_add d.atoms (fun x => if event x then 1 else 0)
    (fun x => if event x then 0 else 1)
  have hc : expectation d.atoms (fun x => (if event x then 1 else 0) +
      (if event x then 0 else 1)) = 1 := by
    calc
      _ = expectation d.atoms (fun _ => 1) := by apply expectation_congr; intro x; split <;> grind
      _ = 1 := d.normalized
  unfold expect at *
  grind

def pure (x : α) : Distribution α where
  atoms := [(x, 1)]
  nonnegative := by intro y hy; simp only [List.mem_singleton] at hy; subst y; change 0 ≤ (1 : Rat); decide
  normalized := by simp only [expectation]; grind

@[simp] theorem expect_pure (x : α) (f : α → Rat) : (pure x).expect f = f x := by
  simp only [pure, expect, expectation]; grind

def compact [DecidableEq α] (d : Distribution α) : Distribution α where
  atoms := collect d.atoms
  nonnegative := collect_nonnegative _ d.nonnegative
  normalized := by rw [expectation_collect]; exact d.normalized

@[simp] theorem expect_compact [DecidableEq α] (d : Distribution α) (f : α → Rat) :
    d.compact.expect f = d.expect f := expectation_collect _ _

def bind (d : Distribution α) (k : α → Distribution β) : Distribution β where
  atoms := flatten d.atoms (fun x => (k x).atoms)
  nonnegative := flatten_nonnegative _ d.nonnegative _ (fun x => (k x).nonnegative)
  normalized := by
    rw [expectation_flatten]
    have h := expectation_congr d.atoms (fun x => expectation (k x).atoms (fun _ => 1))
      (fun _ => 1) (fun x => (k x).normalized)
    rw [h]; exact d.normalized

@[simp] theorem expect_bind (d : Distribution α) (k : α → Distribution β) (f : β → Rat) :
    (d.bind k).expect f = d.expect (fun x => (k x).expect f) := by
  exact expectation_flatten d.atoms (fun x => (k x).atoms) f

/-- Compose using balanced merging to bound individual reduction depth. -/
def bindBalanced [DecidableEq β] (d : Distribution α) (k : α → Distribution β) : Distribution β where
  atoms := flattenMerged d.atoms (fun x => (k x).atoms)
  nonnegative := flattenMerged_nonnegative _ d.nonnegative _ (fun x => (k x).nonnegative)
  normalized := by
    rw [expectation_flattenMerged]
    rw [expectation_congr _ _ (fun _ => 1) (fun x => (k x).normalized)]
    exact d.normalized

@[simp] theorem expect_bindBalanced [DecidableEq β] (d : Distribution α)
    (k : α → Distribution β) (f : β → Rat) :
    (d.bindBalanced k).expect f = d.expect (fun x => (k x).expect f) := by
  exact expectation_flattenMerged d.atoms (fun x => (k x).atoms) f

def bernoulli (p : Rat) (hp : 0 ≤ p ∧ p ≤ 1) : Distribution Bool where
  atoms := [(false, 1-p), (true, p)]
  nonnegative := by
    intro x hx
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
    rcases hx with rfl | rfl <;> simp only <;> grind
  normalized := by simp only [expectation]; grind

end Distribution
end SSC.Finite
