import DesignCompiler.Normal
import DesignCompiler.Search

namespace DesignCompiler.Planning

structure Group where
  label : String
  unitsPerN : Nat := 1
  participantsPerUnit : Nat := 1
  deriving Repr, BEq

structure Criterion where
  alpha : ℚ := 1/20
  targetPower : ℚ := 4/5
  minN : Nat := 2
  maxN : Nat := 100000
  deriving Repr

/-- Covariance of the vector averaged over one independent sampling block.
    An arbitrary linear contrast is allowed, including correlated measurements. -/
structure Contrast where
  labels : List String
  means : List ℚ
  covariance : List (List ℚ)
  coefficients : List ℚ
  deriving Repr

structure Stratum where
  weight : ℚ
  propensity : ℚ
  mean0 : ℚ
  mean1 : ℚ
  variance0 : ℚ
  variance1 : ℚ
  deriving Repr

structure Influence where
  alternative : ℚ
  variance : ℚ
  varianceSource : String
  assumptions : List String
  deriving Repr

inductive Model where
  | gaussianContrast (design : Contrast)
  | efficientATE (strata : List Stratum)
  | asymptoticWald (design : Influence)
  deriving Repr

inductive Analysis where
  | knownCovarianceZ | aipw | tmle | wald
  deriving Repr, BEq

def Analysis.label : Analysis → String
  | .knownCovarianceZ => "Linear contrast with known covariance, two-sided z test"
  | .aipw => "AIPW with efficient influence-function standard error"
  | .tmle => "TMLE with efficient influence-function standard error"
  | .wald => "Asymptotically linear estimator with Wald inference"

structure Study where
  name : String
  estimand : String
  analysis : Analysis
  model : Model
  groups : List Group := [{ label := "participants" }]
  sizeUnit : String := "participants"
  nullValue : ℚ := 0
  criterion : Criterion := {}
  deriving Repr

def dot (a b : List ℚ) : ℚ := (a.zipWith (· * ·) b).sum

def quadraticForm (c : List ℚ) (covariance : List (List ℚ)) : ℚ :=
  dot c (covariance.map (dot c))

def diagonalVariance (d : Contrast) : ℚ :=
  ((List.range d.coefficients.length).map fun i =>
    (d.coefficients[i]!)^2 * (d.covariance[i]!)[i]!).sum

def ate (strata : List Stratum) : ℚ :=
  (strata.map fun s => s.weight * (s.mean1 - s.mean0)).sum

def ateResidualVariance (strata : List Stratum) : ℚ :=
  (strata.map fun s => s.weight * (s.variance1 / s.propensity + s.variance0 / (1-s.propensity))).sum

def ateHeterogeneityVariance (strata : List Stratum) : ℚ :=
  (strata.map fun s => s.weight * (s.mean1-s.mean0-ate strata)^2).sum

def alternative (s : Study) : ℚ := match s.model with
  | .gaussianContrast d => dot d.coefficients d.means
  | .efficientATE strata => ate strata
  | .asymptoticWald d => d.alternative

def variance (s : Study) : ℚ := match s.model with
  | .gaussianContrast d => quadraticForm d.coefficients d.covariance
  | .efficientATE strata => ateResidualVariance strata + ateHeterogeneityVariance strata
  | .asymptoticWald d => d.variance

def effect (s : Study) : ℚ := alternative s - s.nullValue

/-- Exact rational Schur-complement validation, including singular PSD matrices. -/
def psd : Nat → List (List ℚ) → Bool
  | 0, rows => rows.isEmpty
  | fuel+1, rows => match rows with
    | [] => true
    | row :: rest => match row with
      | [] => false
      | pivot :: edge =>
        if pivot < 0 then false
        else if pivot == 0 then edge.all (· == 0) && psd fuel (rest.map List.tail)
        else psd fuel (rest.map fun r =>
          r.tail.zipWith (fun a b => a - r.head! * b / pivot) edge)

def validate (s : Study) : Except String Unit := do
  if s.name.isEmpty || s.estimand.isEmpty || s.sizeUnit.isEmpty then
    throw "Describe the study, estimand, analysis and independent counting unit."
  let c := s.criterion
  if !(0 < c.alpha && c.alpha < c.targetPower && c.targetPower < 1) then
    throw "Require 0 < alpha < targetPower < 1."
  if c.minN == 0 || c.minN > c.maxN || c.maxN > 10000000 then
    throw "Use a nonempty positive integer search domain up to 10,000,000."
  if s.groups.isEmpty || s.groups.any (fun g => g.label.isEmpty || g.unitsPerN == 0 || g.participantsPerUnit == 0) then
    throw "Supply positive allocation counts and participants per independent unit."
  if (s.groups.map (fun g => g.unitsPerN * g.participantsPerUnit)).sum > 1000000 then
    throw "The sampling-block allocation is too large (at most 1,000,000 participants per n)."
  match s.model with
  | .gaussianContrast d =>
    if s.analysis != .knownCovarianceZ then throw "The Gaussian contrast model requires knownCovarianceZ analysis."
    let k := d.coefficients.length
    if k == 0 || k > 32 || d.labels.length != k || d.means.length != k ||
        d.covariance.length != k || d.covariance.any (·.length != k) then
      throw "The contrast, means, labels and square covariance matrix must agree (1–32 dimensions)."
    for i in List.range k do
      for j in List.range k do
        if (d.covariance[i]!)[j]! != (d.covariance[j]!)[i]! then throw "Covariance matrix must be symmetric."
    unless psd k d.covariance do throw "Covariance matrix must be positive semidefinite."
  | .efficientATE strata =>
    if s.analysis != .aipw && s.analysis != .tmle then throw "The efficient ATE variance requires AIPW or TMLE analysis."
    if strata.isEmpty || strata.length > 1000 then throw "Supply 1–1000 baseline strata."
    if (strata.map (·.weight)).sum != 1 then throw "Stratum probabilities must sum exactly to one."
    for r in strata do
      if !(0 < r.weight && 0 < r.propensity && r.propensity < 1 && 0 ≤ r.variance0 && 0 ≤ r.variance1) then
        throw "Require positive stratum probabilities, treatment positivity and nonnegative variances."
    if s.groups != [{ label := "participants" }] then
      -- A participant EIF cannot be reinterpreted as a cluster EIF by relabeling.
      throw "The finite-stratum ATE model counts independent participants; keep the default groups."
    if s.sizeUnit != "participants" then throw "The ATE influence function is per independent participant."
  | .asymptoticWald d =>
    if s.analysis != .wald then throw "A supplied influence variance requires wald analysis."
    if d.varianceSource.isEmpty || d.assumptions.isEmpty then
      throw "State the influence-variance source and asymptotic conditions."
  if !(0 < variance s) then throw "The contrast or influence variance must be strictly positive."

/-- Covariance of a cluster average from an exchangeable individual covariance.
    m is the fixed number of individuals observed per cluster. -/
def clusterMeanVariance (individualVariance icc : ℚ) (m : Nat)
    (_positiveSize : 0 < m := by decide)
    (_nonnegativeVariance : 0 ≤ individualVariance := by decide)
    (_iccBounds : 0 ≤ icc ∧ icc ≤ 1 := by decide) : ℚ :=
  individualVariance * (icc + (1-icc) / m)

/-- The off-diagonal covariances are included, rather than discarded by treating
    repeated measurements as independent. -/
theorem two_measurement_contrast_variance (a b v0 v1 covariance : ℚ) :
    quadraticForm [a,b] [[v0,covariance],[covariance,v1]] =
      a^2*v0 + b^2*v1 + 2*a*b*covariance := by
  simp [quadraticForm, dot]
  ring

theorem exchangeable_average_variance (v rho m : ℚ) (hm : m ≠ 0) :
    (m*v + m*(m-1)*(v*rho)) / m^2 = v*(rho + (1-rho)/m) := by
  field_simp
  ring

end DesignCompiler.Planning
