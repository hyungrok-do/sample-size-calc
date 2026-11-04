import DesignCompiler.Planning.Theory
import DesignCompiler.Solve
import Lean.Replay

open Lean Meta Elab Command

namespace DesignCompiler.Planning

initialize planningResults : SimplePersistentEnvExtension Json (Array Json) ←
  registerSimplePersistentEnvExtension { addEntryFn := Array.push, addImportedFn := fun _ => #[] }

private def q := DesignCompiler.rationalJson
private def text (s : String) := toJson s
private def rats (xs : List ℚ) := toJson (xs.map q)
private def showRat (r : ℚ) : String := if r.den == 1 then s!"{r.num}" else s!"{r.num}/{r.den}"
private def step (title formula detail : String) : Json :=
  Json.mkObj [("title",text title),("formula",text formula),("detail",text detail)]

private def modelReport (s : Study) : Json × List String × Array Json := Id.run do
  match s.model with
  | .gaussianContrast d =>
    let model := Json.mkObj [("kind",text "gaussian_contrast"),("labels",toJson d.labels),
      ("means",rats d.means),("contrast",rats d.coefficients),
      ("covariance",toJson (d.covariance.map rats)),
      ("diagonalVariance",q (diagonalVariance d)),("crossCovariance",q (variance s-diagonalVariance d))]
    let assumptions := ["Independent sampling blocks with the displayed allocation and covariance.",
      "The block outcome vector is jointly Gaussian with the specified population covariance known.",
      "The analysis is a linear contrast of block means; estimating covariance instead requires a different reference distribution or an asymptotic approximation."]
    return (model, assumptions, #[
      step "Estimand and planned estimator" s!"psi = c^T mu = {showRat (alternative s)}"
        s!"{s.estimand}. Estimate with c^T times the vector of sample means.",
      step "Propagate the full covariance" s!"V = c^T Sigma c = {showRat (variance s)}; Var(psi_hat) = V/n"
        s!"Diagonal contribution: {showRat (diagonalVariance d)}. Cross-covariance contribution: {showRat (variance s-diagonalVariance d)}. Both enter V; n counts independent sampling blocks."])
  | .efficientATE strata =>
    let model := Json.mkObj [("kind",text "efficient_ate"),
      ("strata",toJson (strata.map fun r => Json.mkObj [("weight",q r.weight),("propensity",q r.propensity),
        ("mean0",q r.mean0),("mean1",q r.mean1),("variance0",q r.variance0),("variance1",q r.variance1)])),
      ("residualVariance",q (ateResidualVariance strata)),("heterogeneityVariance",q (ateHeterogeneityVariance strata))]
    let assumptions := ["Independent participants and finite baseline strata with treatment positivity.",
      "Consistency and conditional exchangeability for the causal ATE interpretation.",
      "The planned AIPW/TMLE estimator attains the stated efficient influence function with a negligible remainder and consistent standard error.",
      "Nuisance estimation must meet its convergence conditions (with suitable sample splitting or empirical-process conditions). Double robustness alone does not establish efficiency.",
      "Power uses an asymptotic normal approximation, without a finite-sample guarantee."]
    return (model, assumptions, #[
      step "Identify the target ATE" s!"psi = sum_w p(w) [mu1(w)-mu0(w)] = {showRat (alternative s)}"
        s.estimand,
      step "Construct the estimator's influence function"
        "D = A/g(W) [Y-mu1(W)] - (1-A)/(1-g(W)) [Y-mu0(W)] + mu1(W)-mu0(W)-psi"
        "Use the efficient influence function for the planned AIPW/TMLE analysis under the listed conditions.",
      step "Derive the influence variance"
        s!"V = E[v1(W)/g(W) + v0(W)/(1-g(W))] + Var[tau(W)] = {showRat (ateResidualVariance strata)} + {showRat (ateHeterogeneityVariance strata)} = {showRat (variance s)}"
        "The first component is outcome noise weighted by treatment probabilities. The second is treatment-effect heterogeneity."])
  | .asymptoticWald d =>
    return (Json.mkObj [("kind",text "asymptotic_wald"),("varianceSource",text d.varianceSource)],
      d.assumptions ++ ["Asymptotic linearity, a negligible remainder, and a consistent standard error on the specified effect scale.",
        "Power is an asymptotic normal approximation, without a finite-sample guarantee."], #[
      step "Estimand and estimator" s!"psi = {showRat (alternative s)}" s.estimand,
      step "Influence representation" s!"sqrt(n) (psi_hat-psi) converges to Normal(0,V), V = {showRat (variance s)}" d.varianceSource])

private def checkedValue (name function : Name) (source : Expr) (value : ℚ) : TermElabM (String × Array Name) := do
  let literal ← Term.elabTermAndSynthesize (← `(($(quote value.num) : ℚ) / ($(quote value.den) : ℚ))) (some (mkConst ``Rat))
  let proposition ← instantiateMVars (← mkEq (← mkAppM function #[source]) literal)
  let proof ← instantiateMVars (← mkDecideProof proposition)
  withOptions (fun o => o.setBool `Elab.async false |>.setBool `debug.skipKernelTC false) do
    addDecl <| .thmDecl { name, levelParams := [], type := proposition, value := proof }
  let axioms ← collectAxioms name
  unless axioms.all (fun a => [``propext, ``Classical.choice, ``Quot.sound].contains a) do
    throwError "Statistical arithmetic certificate depends on an unapproved axiom"
  return ((← ppExpr proposition).pretty, axioms)

private def compileStudy (name : Name) (source : TSyntax `term) (mode : String) (count : Nat) : CommandElabM Unit := do
  let fullName := (← getCurrNamespace) ++ name
  let report ← runTermElabM fun _ => do
    let expr ← Term.elabTermAndSynthesize source (some (mkConst ``Study))
    if expr.hasFVar || expr.hasMVar then throwError "Supply a closed planning study with numerical planning assumptions."
    let study ← unsafe evalExpr Study (mkConst ``Study) expr
    match validate study with
    | .error e => throwError e
    | .ok _ => pure ()
    if mode != "sample_size" && (count == 0 || count > 10000000) then
      throwError "Power and detectable-effect queries require a count from 1 to 10,000,000."
    if mode == "effect" then
      if let .efficientATE _ := study.model then
        throwError "For ATE power, change the planning strata explicitly: its influence variance can change with the effect."
    let (varianceType, ax1) ← checkedValue (fullName.appendAfter "_variance") ``variance expr (variance study)
    let (effectType, ax2) ← checkedValue (fullName.appendAfter "_effect") ``effect expr (effect study)
    let (model, assumptions, initialSteps) := modelReport study
    let c := study.criterion
    let steps := initialSteps ++ #[
      step "Specify the null and rejection rule"
        s!"H0: psi = {showRat study.nullValue}; reject if |Z| > z_(1-alpha/2), alpha = {showRat c.alpha}"
        "Z = (psi_hat - psi0) / sqrt(V/n). This version uses a two-sided Wald rejection rule.",
      step "Obtain the sampling distribution under the alternative"
        s!"Z ~ Normal(lambda(n),1), lambda(n) = ({showRat (effect study)}) sqrt(n/({showRat (variance study)}))"
        (match study.model with | .gaussianContrast _ => "Exact under the declared Gaussian known-covariance model." | _ => "Asymptotic normal approximation under the stated estimator conditions."),
      step "Integrate both rejection tails"
        "power(n) = Phi(-z - lambda(n)) + 1 - Phi(z - lambda(n))"
        "Uses the checked shifted-normal rejection-probability identity. The two tails are retained in the numerical calculation.",
      step "Solve the design criterion"
        (if mode == "sample_size" then s!"min n in [{c.minN}, {c.maxN}] with power(n) >= {showRat c.targetPower}"
         else if mode == "power" then s!"Evaluate power at n = {count}"
         else s!"Find the smallest |psi-psi0| attaining power {showRat c.targetPower} at n = {count}")
        (if mode == "sample_size" then "Compile the criterion into an ascending integer search. Compare the selected count with its predecessor."
         else if mode == "power" then "Keep the specified estimand, covariance and allocation fixed."
         else "Hold V fixed when inverting the normal noncentrality; this is not valid if the planning variance changes with the effect.")]
    return Json.mkObj [("kind",text "statistical_design"),("name",text fullName.toString),("title",text study.name),
      ("solveFor",text mode),("fixedN",toJson count),("estimand",text study.estimand),("analysis",text study.analysis.label),
      ("sizeUnit",text study.sizeUnit),("model",model),("assumptions",toJson assumptions),("steps",toJson steps),
      ("groups",toJson (study.groups.map fun g => Json.mkObj [("label",text g.label),("unitsPerN",toJson g.unitsPerN),("participantsPerUnit",toJson g.participantsPerUnit)])),
      ("statistics",Json.mkObj [("alternative",q (alternative study)),("effect",q (effect study)),("nullValue",q study.nullValue),("variance",q (variance study))]),
      ("criterion",Json.mkObj [("alpha",q c.alpha),("targetPower",q c.targetPower),("minN",toJson c.minN),("maxN",toJson c.maxN)]),
      ("theoremType",text s!"{fullName}_variance : {varianceType}\n{fullName}_effect : {effectType}"),
      ("axioms",toJson ((ax1 ++ ax2).toList.eraseDups.map Name.toString)),
      ("verification",text "Lean checks exact rational effect and variance certificates and the generic normal-tail identity. The sampling-model assumptions, CLT conditions and floating-point executable are not certified by those arithmetic proofs.")]
  modifyEnv fun env => planningResults.addEntry env report
  logInfo m!"{fullName}: statistical {mode} plan compiled"

-- The query targets are identifiers, so ordinary Lean definitions named
-- `power` or `effect` remain valid elsewhere in a user module.
syntax (name := planCountQuery) "solve " ident ident " using " term : command
syntax (name := planFixedQuery) "solve " ident ident " using " term " at " num : command

elab_rules : command
  | `(solve $mode:ident $name:ident using $source at $count:num) => do
    let target := mode.getId.toString
    unless target == "power" || target == "effect" do
      throwError "A fixed-count query must solve power or effect."
    compileStudy name.getId source target count.getNat
  | `(solve $mode:ident $name:ident using $source) => do
    unless mode.getId == `sample_size do
      throwError "Use solve sample_size NAME using STUDY, or supply 'at N' for power/effect."
    compileStudy name.getId source "sample_size" 0

elab "#export_compiler" : command => do
  let env ← getEnv
  -- Audit the whole submitted module, including declarations unused by solve.
  -- A sorry warning can be silenced, so diagnostics alone are not proof evidence.
  let constants := env.constants.foldStage2 (fun m n ci => m.insert n ci)
    ({} : Std.HashMap Name ConstantInfo)
  let some imported := env.importEnv?
    | throwError "Cannot independently check declarations without the import environment."
  discard <| liftIO <| imported.replay constants
  let declarations := env.constants.foldStage2 (fun ns n _ => ns.push n) #[]
  for name in declarations do
    let axioms ← collectAxioms name
    unless axioms.all (fun a => [``propext, ``Classical.choice, ``Quot.sound].contains a) do
      throwError "Unproved dependency in {name}: {axioms.toList}. Unfinished proofs and user axioms are not accepted."
  let equations := (DesignCompiler.Equation.resultExtension.getState env).map toJson
  let studies := planningResults.getState env
  liftIO <| IO.println ("SSC_LEAN_REPORT " ++ (toJson (studies ++ equations)).compress)

end DesignCompiler.Planning
