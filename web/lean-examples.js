export const examples = {
  normal: `import DesignCompiler.Input
open DesignCompiler

def study : NormalInput := {
  name := "Two-arm normal comparison"
  alternative := 0.5
  nullValue := 0
  sd := 1
  alpha := 0.05
  targetPower := 0.80
  minN := 1
  maxN := 10000
}
`,
  influence: `import DesignCompiler.Input
open DesignCompiler

def study : InfluenceInput := {
  name := "ATE with influence-function variance"
  alternative := 0.5
  nullValue := 0
  variance := 4
  unit := .participant
  estimand := "ATE"
  estimator := "AIPW"
  varianceSource := "Balanced trial, residual variance 1, constant effect"
  alpha := 0.05
  targetPower := 0.80
  minN := 1
  maxN := 10000
}
`,
  cluster: `import DesignCompiler.Input
open DesignCompiler

def study : InfluenceInput := {
  name := "Planning with independent clusters"
  alternative := 0.5
  nullValue := 0
  variance := 2
  unit := .cluster
  estimand := "Cluster-average contrast"
  estimator := "Asymptotically linear cluster estimator"
  varianceSource := "Illustrative variance per independent cluster; replace for your design"
  alpha := 0.05
  targetPower := 0.80
  minN := 10
  maxN := 10000
}
`,
};

export const programs = {
  "normalPower": "-- The compiler derives the Gaussian rejection probability from this event.\ndef effect : Real := 0.5\ndef alpha : Real := 0.05\ndef critical : Real := normalQuantile(1 - alpha / 2)\n\ndef power (n : Real) : Real :=\n  let z := normal(effect * sqrt(n / 2), 1)\n  probability(abs(z) > critical)\n\nsolve for n : Nat in [2, 10000]\nwhere power(n) >= 0.80\n",
  "effectEquation": "-- n is fixed. Solve the user-defined power equation for a treatment effect.\ndef n : Real := 63\ndef critical : Real := normalQuantile(0.975)\n\ndef power (effect : Real) : Real :=\n  probability(abs(normal(effect * sqrt(n / 2), 1)) > critical)\n\nsolve for effect : Real in [0.01, 2] tolerance 0.00000001\nwhere power(effect) = 0.80\n",
  "equation": "-- An arbitrary equation, with no predefined statistical procedure.\ndef objective (x : Real) : Real := x^2 + 2*x\n\nsolve for x : Real in [0, 10]\nwhere objective(x) = 8\n",
  "simulation": "-- Skewed outcomes and an estimated standard error: simulate the actual test.\n-- This is the rejection probability of the specified z-threshold rule.\ndef effect : Real := 0.8\ndef critical : Real := normalQuantile(0.975)\n\ndef outcome (shift : Real) : Real :=\n  let x := normal(0, 1)\n  shift + 0.3*x + exp(normal(0, 0.6))\n\ndef reject (n : Nat) : Bool :=\n  let control := sample(n, fun i => outcome(0))\n  let treated := sample(n, fun i => outcome(effect))\n  let se := sqrt(variance(control) / n + variance(treated) / n)\n  abs(mean(treated) - mean(control)) / se > critical\n\ndef power (n : Nat) : Real := probability(reject(n))\n\nsimulation { replicates := 5000, seed := 2026, confidence := 0.95 }\nsolve for n : Nat in [10, 60] step 10\nwhere power(n) >= 0.80\n",
  "integration": "-- Composite Simpson quadrature, with an explicit even panel count.\ndef area (upper : Real) : Real :=\n  integrate(fun x => exp(-x*x), 0, upper, 200)\n\nsolve for upper : Real in [0, 3]\nwhere area(upper) = 0.7\n",
  "estimate": "-- Events with nonlinear random transformations automatically use simulation.\ndef event (shift : Real) : Bool := exp(normal(0, 1)) + shift > 2\nsimulation { replicates := 20000, seed := 19, confidence := 0.95 }\nevaluate probability(event(0.2))\n"
};
