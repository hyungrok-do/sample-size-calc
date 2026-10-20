// Group-sequential procedures for two means, two proportions and the logrank
// test, built on the boundary/exit-probability engine in groupseq.cpp.
#include <map>
#include <memory>

#include "procedures/common.hpp"
#include "procedures/groupseq.hpp"
#include "util/memo.hpp"
#include "util/numtext.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

const char* G_GS = "Group-Sequential Design";
const char* G_FUT = "Futility Boundary";

const std::vector<std::string> SPENDING = {"obf", "pocock", "power", "hsd"};

// ------------------------------------------------------------------ shared parameters
std::vector<ParamDef> gs_params(const std::string& upper_label, const std::string& lower_label) {
  const Clause one_sided = is_not("alternative", {"two"});
  ParamDef power = power_param();
  power.hint(
      "Probability of rejecting H0 when H1 is true (1 − β). With β-spending futility boundaries it also sets the β "
      "that the futility boundaries spend, also when solving for power.");
  power.show(when_any({Conj{not_solving({"power"})}, Conj{one_sided, is("futility", {"beta"})}}));
  return {
      select("alternative", "Alternative Hypothesis",
             {{"two", "Two-sided (symmetric boundaries)", {}}, {"upper", upper_label, {}}, {"lower", lower_label, {}}},
             "two")
          .grp(G_TEST)
          .hint("Two-sided designs reject for large |Z| at any analysis; futility boundaries need a one-sided design."),
      num("looks", "K (Number of Analyses, Including Final)", "3")
          .grp(G_GS)
          .sym("K")
          .range(1, 20)
          .whole()
          .single()
          .sugg({"2", "3", "4", "5"})
          .hint("Total number of analyses, counting the final one; 3 means two interim analyses and a final "
                "analysis."),
      select("spacing", "Analysis Times", {{"equal", "Equally spaced information", {}}, {"custom", "Custom information fractions", {}}}, "equal")
          .grp(G_GS)
          .hint("Hold the analyses at equal steps of information (subjects or events), or give the fraction of the "
                "total information at each analysis."),
      vec("tfrac", "Information Fractions (e.g. 0.5 0.75 1)", "0.5 0.75 1")
          .grp(G_GS)
          .sym("Fractions")
          .above(0)
          .hint("Relative information (sample size or events) at each analysis; rescaled so the last is 1.")
          .show(when(is("spacing", {"custom"}))),
      select("boundary", "Efficacy Boundary",
             {{"obf", "O'Brien–Fleming type spending (Lan–DeMets)", {}},
              {"pocock", "Pocock type spending (Lan–DeMets)", {}},
              {"power", "Power family spending α·t^ρ (Kim–DeMets)", {}},
              {"hsd", "Hwang–Shih–DeCani spending (γ)", {}},
              {"wt_obf", "O'Brien–Fleming (classical)", {}},
              {"wt_pocock", "Pocock (classical)", {}},
              {"wt", "Wang–Tsiatis (Δ)", {}}},
             "obf")
          .grp(G_GS)
          .hint("How the chance of stopping for efficacy is spread over the analyses. O'Brien–Fleming-type "
                "boundaries are strict early, leaving the final test near alpha; Pocock-type ones are alike at "
                "every look."),
      num("rho", "ρ (Power Family Exponent)", "3")
          .grp(G_GS)
          .sym("ρ")
          .above(0)
          .single()
          .sugg({"1", "2", "3"})
          .hint("Exponent of the spending function α·t^ρ: larger values spend less alpha early (3 is close to "
                "O'Brien–Fleming, 1 to Pocock).")
          .show(when(is("boundary", {"power"}))),
      num("gamma", "γ (Hwang–Shih–DeCani Parameter)", "-4")
          .grp(G_GS)
          .sym("γ")
          .range(-40, 40)
          .single()
          .sugg({"-4", "-2", "1"})
          .hint("Shape of the Hwang–Shih–DeCani spending function: negative values spend less alpha early (−4 is "
                "close to O'Brien–Fleming, 1 to Pocock).")
          .show(when(is("boundary", {"hsd"}))),
      num("wtdelta", "Δ (Wang–Tsiatis Shape)", "0.25")
          .grp(G_GS)
          .sym("Δ")
          .range(-0.5, 0.5)
          .single()
          .sugg({"0", "0.25", "0.5"})
          .hint("Shape of the Wang–Tsiatis boundaries: 0 gives O'Brien–Fleming and 0.5 Pocock boundaries; values "
                "in between are intermediate.")
          .show(when(is("boundary", {"wt"}))),
      select("futility", "Futility Boundary",
             {{"none", "None", {}},
              {"beta", "β-spending (Lan–DeMets type)", when(Conj{one_sided, is("boundary", SPENDING)})},
              {"z", "Fixed z-values", when(one_sided)}},
             "none")
          .grp(G_FUT)
          .hint("Boundaries for stopping early when the data favor H0. β-spending boundaries spend the type II error "
                "β = 1 − power over the analyses, like alpha-spending does for the efficacy boundaries.")
          .show(when(one_sided)),
      select("bspend", "β-Spending Function",
             {{"obf", "O'Brien–Fleming type", {}}, {"pocock", "Pocock type", {}}, {"power", "Power family β·t^ρ", {}},
              {"hsd", "Hwang–Shih–DeCani (γ)", {}}},
             "obf")
          .grp(G_FUT)
          .hint("How the type II error β is spread over the analyses by the futility boundaries; the choices match "
                "the alpha-spending functions of the efficacy boundary.")
          .show(when(Conj{one_sided, is("futility", {"beta"})})),
      num("brho", "ρ (β-Spending Exponent)", "2")
          .grp(G_FUT)
          .sym("ρ (β)")
          .above(0)
          .single()
          .sugg({"1", "2", "3"})
          .hint("Exponent of the β-spending function β·t^ρ: larger values make stopping for futility early less "
                "likely.")
          .show(when(Conj{one_sided, is("futility", {"beta"}), is("bspend", {"power"})})),
      num("bgamma", "γ (β-Spending Parameter)", "-2")
          .grp(G_FUT)
          .sym("γ (β)")
          .range(-40, 40)
          .single()
          .sugg({"-4", "-2", "1"})
          .hint("Shape of the Hwang–Shih–DeCani β-spending function: negative values make stopping for futility "
                "early less likely.")
          .show(when(Conj{one_sided, is("futility", {"beta"}), is("bspend", {"hsd"})})),
      vec("fz", "Futility Z Values (Interim Analyses)", "0")
          .grp(G_FUT)
          .sym("Futility Z")
          .hint("Stop for futility when Z ≤ the value: one value per interim analysis, or a single value used at every "
                "interim analysis.")
          .show(when(Conj{one_sided, is("futility", {"z"})})),
      select("binding", "Futility Boundaries Are",
             {{"nonbinding", "Non-binding (alpha computed without them)", {}},
              {"binding", "Binding (alpha computed with them)", {}}},
             "nonbinding")
          .grp(G_FUT)
          .hint("Non-binding boundaries keep the type I error at alpha even if the trial continues past them.")
          .show(when(Conj{one_sided, is_not("futility", {"none"})})),
      power,
      alpha_param("0.05", "Alpha (Total; Split Between Tails if Two-Sided)"),
  };
}

std::vector<double> gs_fractions(const Scenario& s) {
  const int K = static_cast<int>(s("looks"));
  std::vector<double> t;
  if (s.is("spacing", "custom")) {
    t = s.v("tfrac");
    if (static_cast<int>(t.size()) != K)
      throw CalcError("enter " + std::to_string(K) + " information fractions (one per analysis)");
    const double last = t.back();
    for (double& v : t) v /= last;
  } else {
    for (int k = 1; k <= K; ++k) t.push_back(static_cast<double>(k) / K);
  }
  return t;
}

bool two_sided(const Scenario& s) { return s.is("alternative", "two"); }

std::string futility_kind(const Scenario& s) {
  if (two_sided(s) || !s.has("futility") || s("looks") < 2) return "none";
  return s.s("futility");
}

double spending_param(const std::string& name, const Scenario& s, const char* rho, const char* gamma) {
  return name == "power" ? s(rho) : (name == "hsd" ? s(gamma) : 0.0);
}

struct DesignInfo {
  gs::Design d;
  gs::ExitProbs h0;  // exit probabilities under H0 (futility rule followed)
};
using DesignPtr = std::shared_ptr<const DesignInfo>;

// Designs depend only on a few inputs but are needed for every scenario and
// every point of the power curve, so they are built once and shared between
// threads (a waiting thread blocks until the first one has finished).
DesignPtr cached_design(const std::string& key, const std::function<DesignInfo()>& build, const RunContext& ctx) {
  static OnceCache<std::string, DesignInfo> cache;
  // another run cancelled while building it: built again unless this one is cancelled too
  return cache.get(key, build, [&] { ctx.check(); });
}

DesignPtr gs_design(const Scenario& s, const RunContext& ctx) {
  const auto t = gs_fractions(s);
  const double alpha = s("alpha");
  const bool two = two_sided(s);
  const std::string& b = s.s("boundary");
  const std::string fk = futility_kind(s);
  std::vector<double> fz;
  if (fk == "z") {
    fz = s.v("fz");
    const std::size_t m = t.size() - 1;
    if (fz.size() == 1 && m > 1) fz.assign(m, fz[0]);
    if (fz.size() != m)
      throw CalcError("enter " + std::to_string(m) + " futility values (one per interim analysis) or a single value");
  }
  const double beta = fk == "beta" ? 1 - s("power") : 0.0;
  if (fk == "beta" && !(1 - beta > alpha)) throw CalcError("the power must exceed alpha");

  std::string key;
  for (double v : t) key += text::format_roundtrip(v) + ",";
  key += "|" + text::format_roundtrip(alpha) + (two ? "|2|" : "|1|") + b;
  if (b == "power" || b == "hsd" || b == "wt")
    key += ":" + text::format_roundtrip(b == "power" ? s("rho") : (b == "hsd" ? s("gamma") : s("wtdelta")));
  key += "|" + fk;
  if (fk != "none") key += s.is("binding", "binding") ? ":binding" : ":nonbinding";
  if (fk == "beta") {
    const std::string& bs = s.s("bspend");
    key += ":" + bs + ":" + text::format_roundtrip(spending_param(bs, s, "brho", "bgamma")) + ":" +
           text::format_roundtrip(beta);
  }
  for (double v : fz) key += ":" + text::format_roundtrip(v);

  return cached_design(
      key,
      [&] {
        gs::Efficacy eff;
        if (b == "wt_obf" || b == "wt_pocock" || b == "wt") {
          eff.kind = gs::Efficacy::Kind::WangTsiatis;
          eff.delta = b == "wt_obf" ? 0.0 : (b == "wt_pocock" ? 0.5 : s("wtdelta"));
        } else {
          eff.spend = gs::alpha_spending(b, alpha, two, spending_param(b, s, "rho", "gamma"));
        }
        gs::Futility fut;
        if (fk == "beta") {
          const std::string& bs = s.s("bspend");
          fut.kind = gs::Futility::Kind::BetaSpending;
          fut.beta = beta;
          fut.spend = gs::spending(bs, beta, spending_param(bs, s, "brho", "bgamma"));
        } else if (fk == "z") {
          fut.kind = gs::Futility::Kind::Fixed;
          fut.z = fz;
        }
        fut.binding = s.is("binding", "binding");
        DesignInfo info;
        info.d = gs::make_design(t, alpha, two, eff, fut, [&] { ctx.check(); });
        info.h0 = gs::exit_probs(info.d, 0.0);
        return info;
      },
      ctx);
}

std::string boundary_name(const Scenario& s) {
  static const std::map<std::string, std::string> names = {
      {"obf", "O'Brien–Fleming-type spending"}, {"pocock", "Pocock-type spending"}, {"power", "power-family spending"},
      {"hsd", "Hwang–Shih–DeCani spending"},     {"wt_obf", "classical O'Brien–Fleming"}, {"wt_pocock", "classical Pocock"},
      {"wt", "Wang–Tsiatis"}};
  return names.at(s.s("boundary"));
}

std::string futility_name(const Scenario& s) {
  const std::string fk = futility_kind(s);
  if (fk == "none") return "none";
  static const std::map<std::string, std::string> fn = {
      {"obf", "O'Brien–Fleming type"}, {"pocock", "Pocock type"}, {"power", "power family"}, {"hsd", "Hwang–Shih–DeCani"}};
  const std::string kind = fk == "beta" ? "β-spending (" + fn.at(s.s("bspend")) + ")" : "fixed z-values";
  return kind + (s.is("binding", "binding") ? ", binding" : ", non-binding");
}

// signed drift in the direction of the alternative (|·| for two-sided designs)
double directed(const Scenario& s, double eta) {
  if (two_sided(s)) return std::fabs(eta);
  return s.is("alternative", "lower") ? -eta : eta;
}

// +1 when large values of the effect estimate favor H1, −1 otherwise
double direction(const Scenario& s, double effect) {
  if (two_sided(s)) return effect < 0 ? -1.0 : 1.0;
  return s.is("alternative", "lower") ? -1.0 : 1.0;
}

// drift of the fixed-sample test with the same alpha and power
// Drift of the fixed-sample test with the same α and power (NaN if none).  As in
// the group-sequential power, two-sided power counts both rejection regions.
double fixed_drift(const Scenario& s, double power, double scale) {
  const double a = s("alpha");
  const double one_tail = scale * stats::qnorm(two_sided(s) ? a / 2 : a, false) + stats::qnorm(power);
  if (!two_sided(s)) return one_tail;
  const double c = scale * stats::qnorm(a / 2, false);
  auto excess = [&](double eta) { return stats::pnorm(eta - c) + stats::pnorm(-eta - c) - power; };
  if (excess(0) >= 0) return NAN;
  return numerics::brent(excess, 0, one_tail, 1e-13);  // the other region only adds power
}

// Drift needed for the target power.
double target_drift(const Scenario& s, const DesignInfo& info, double scale) {
  if (scale == 1 && std::isfinite(info.d.eta_design)) return info.d.eta_design;
  return gs::drift_for_power(info.d, s("power"), scale);
}

void check_power_target(const Scenario& s) {
  if (!(s("power") > s("alpha"))) throw CalcError("the power must exceed alpha");
}

// Common outputs: power, expected sizes, inflation and the boundary tables.
// `effect(k, z)` converts a z boundary of analysis k to the effect scale.
// `eta_target` is the drift solved for the target power (NaN when solving for power).
void gs_report(const Scenario& s, Row& r, const DesignInfo& info, double eta, double eta_target, double scale,
               double nmax, const std::function<double(std::size_t, double)>& effect, const std::string& unit,
               const std::string& effect_label) {
  const gs::Design& d = info.d;
  const auto e1 = gs::exit_probs(d, eta, scale);
  const auto& e0 = info.h0;
  const std::size_t K = d.t.size();
  const bool fut = d.has_futility();
  const bool beta_spending = !d.beta_spent.empty();
  auto zcell = [](double z) { return std::isfinite(z) ? Cell(z) : Cell(std::string(z > 0 ? "∞" : "−∞")); };
  auto ecell = [&](std::size_t k, double z) { return std::isfinite(z) ? Cell(effect(k, z)) : Cell(std::string("—")); };

  Table tb;
  tb.title = "Boundaries";
  tb.columns = {{"k", "Analysis"},
                {"t", "Info Fraction"},
                {"n", unit},
                {"z", d.two_sided ? "Efficacy Z (±)" : "Efficacy Z"},
                {"p", "Nominal p"},
                {"eff", "Efficacy " + effect_label},
                {"spent", "Cum. α"}};
  tb.formats = {"int", "num", "num", "num", "prob", "num", "prob"};
  if (fut) {
    tb.columns.push_back({"fz", "Futility Z"});
    tb.columns.push_back({"feff", "Futility " + effect_label});
    tb.formats.push_back("num");
    tb.formats.push_back("num");
    if (beta_spending) {
      tb.columns.push_back({"bspent", "Cum. β"});
      tb.formats.push_back("prob");
    }
  }
  for (std::size_t k = 0; k < K; ++k) {
    const double z = d.b[k];
    std::vector<Cell> row = {Cell(static_cast<double>(k + 1)), Cell(d.t[k]), Cell(nmax * d.t[k]), zcell(z),
                             Cell(std::isfinite(z) ? stats::pnorm(z, false) : 0.0), ecell(k, z), Cell(d.alpha_spent[k])};
    if (fut) {
      const double a = d.a[k];
      const bool last = k + 1 == K;
      row.push_back(last ? Cell(std::string("—")) : zcell(a));
      row.push_back(last || !std::isfinite(a) ? Cell(std::string("—")) : ecell(k, a));
      if (beta_spending) row.push_back(Cell(d.beta_spent[k]));
    }
    tb.rows.push_back(std::move(row));
  }
  r.tables.push_back(std::move(tb));

  Table sp;
  sp.title = "Stopping Probabilities";
  sp.columns = {{"k", "Analysis"}, {"e1", "P(Efficacy | H1)"}};
  sp.formats = {"int", "prob"};
  if (fut) {
    sp.columns.push_back({"f1", "P(Futility | H1)"});
    sp.formats.push_back("prob");
  }
  sp.columns.push_back({"cum", "Cumulative Power"});
  sp.columns.push_back({"e0", "P(Efficacy | H0)"});
  sp.formats.push_back("prob");
  sp.formats.push_back("prob");
  if (fut) {
    sp.columns.push_back({"f0", "P(Futility | H0)"});
    sp.formats.push_back("prob");
  }
  double cum = 0;
  for (std::size_t k = 0; k < K; ++k) {
    const double p1 = e1.up[k] + e1.lo[k];
    cum += p1;
    std::vector<Cell> row = {Cell(static_cast<double>(k + 1)), Cell(p1)};
    if (fut) row.push_back(Cell(e1.fut[k]));
    row.push_back(Cell(cum));
    row.push_back(Cell(e0.up[k] + e0.lo[k]));
    if (fut) row.push_back(Cell(e0.fut[k]));
    sp.rows.push_back(std::move(row));
  }
  r.tables.push_back(std::move(sp));

  put_power(s, r, e1.reject);
  r.set("looks", static_cast<double>(K));
  r.set("en1", nmax * e1.info);
  r.set("en0", nmax * e0.info);
  const bool solved = std::isfinite(eta_target);
  const double p = solved ? s("power") : e1.reject;
  if (p > s("alpha") && p < 1) {
    const double eta_p = solved ? eta_target : eta;
    const double f = fixed_drift(s, p, scale);
    if (f > 0) r.set("inflation", eta_p * eta_p / (f * f));
  }
  r.set("bound", boundary_name(s));
  r.set("futility", futility_name(s));
  r.set("alpha", s("alpha"));
  if (fut && !d.binding) {
    double a0 = 0;
    for (std::size_t k = 0; k < K; ++k) a0 += e0.up[k];
    r.notes.push_back("With the futility boundaries followed, the type I error is " + fnum(a0, 5) +
                      " (non-binding boundaries are ignored when spending alpha).");
  }
}

std::vector<Column> gs_columns(const Scenario& s, const std::vector<Column>& effect_cols, const char* e_label = "E(N|H1)",
                               const char* e0_label = "E(N|H0)") {
  auto c = power_columns(s);
  for (const auto& x : two_group_columns()) c.push_back(x);
  c.push_back(col("looks", "K", "int", "Number of analyses"));
  for (const auto& x : effect_cols) c.push_back(x);
  for (auto x : {col("en1", e_label, "num", "Expected size when H1 is true"),
                 col("en0", e0_label, "num", "Expected size when H0 is true"),
                 col("inflation", "Inflation", "num", "Maximum information relative to a fixed design"),
                 col("alpha", "Alpha", "alpha"), col("bound", "Boundary", "text")})
    c.push_back(x);
  if (futility_kind(s) != "none") c.push_back(col("futility", "Futility", "text"));
  return c;
}

std::string text_cell(const Row& r, const std::string& key) {
  for (const auto& kv : r.cells)
    if (kv.first == key)
      if (const std::string* p = std::get_if<std::string>(&kv.second)) return *p;
  return {};
}

std::string design_words(const Scenario& s, const Row& r) {
  std::string w = fint(r.num("looks")) + (r.num("looks") == 1 ? " analysis (" : " analyses (") + text_cell(r, "bound") + ", " +
                  (two_sided(s) ? "two-sided" : "one-sided") + " overall α = " + fnum(r.num("alpha"));
  const std::string f = text_cell(r, "futility");
  if (!f.empty() && f != "none") w += "; futility: " + f;
  return w + ")";
}

const std::vector<std::string> GS_REFS = {
    "Jennison, C. & Turnbull, B. W. (2000). Group Sequential Methods with Applications to Clinical Trials. Chapman & "
    "Hall/CRC.",
    "Lan, K. K. G. & DeMets, D. L. (1983). Discrete sequential boundaries for clinical trials. Biometrika 70, 659–663.",
    "O'Brien, P. C. & Fleming, T. R. (1979). A multiple testing procedure for clinical trials. Biometrics 35, 549–556.",
    "Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. Biometrika 64, "
    "191–199.",
    "Kim, K. & DeMets, D. L. (1987). Design and analysis of group sequential tests based on the type I error spending "
    "rate function. Biometrika 74, 149–154.",
    "Hwang, I. K., Shih, W. J. & DeCani, J. S. (1990). Group sequential designs using a family of type I error "
    "probability spending functions. Statistics in Medicine 9, 1439–1445.",
    "Wang, S. K. & Tsiatis, A. A. (1987). Approximately optimal one-parameter boundaries for group sequential trials. "
    "Biometrics 43, 193–199.",
    "Pampallona, S., Tsiatis, A. A. & Kim, K. (2001). Interim monitoring of group sequential trials using spending "
    "functions for the type I and type II error probabilities. Drug Information Journal 35, 1113–1121."};

const std::vector<std::pair<std::string, std::string>> GS_DEFS = {
    {"K", "Number of analyses (interim analyses plus the final analysis)."},
    {"Info Fraction", "Proportion of the maximum information (sample size or events) at an analysis."},
    {"Efficacy Z", "Critical value of the standardized statistic for rejecting H0 (±Z for two-sided designs)."},
    {"Nominal p", "One-sided p-value threshold equivalent to the efficacy Z boundary."},
    {"Cum. α", "Type I error spent up to and including the analysis."},
    {"Futility Z", "The trial stops for futility when the statistic is at or below this value."},
    {"Cum. β", "Type II error spent by the futility boundaries up to and including the analysis."},
    {"E(N|H1), E(N|H0)", "Expected sample size (or events) when H1 / H0 is true, with all stopping rules followed."},
    {"Inflation", "Maximum information relative to the fixed-sample design with the same α and power."}};

// ======================================================================= two means
Procedure gs_means() {
  Procedure p;
  p.id = "gs.means";
  p.name = "Group-Sequential Tests for Two Means";
  p.category = "Group-Sequential";
  p.subcategory = "Two Means";
  p.order = 1;
  p.summary =
      "Maximum sample size, boundaries and expected sample size for a group-sequential comparison of two means "
      "(z-test, known σ), with alpha-spending (O'Brien–Fleming, Pocock, power family, Hwang–Shih–DeCani) or "
      "Wang–Tsiatis efficacy boundaries and optional β-spending or fixed futility boundaries.";
  p.keywords = {"group sequential", "interim analysis", "alpha spending", "o'brien-fleming", "pocock", "lan-demets",
                "early stopping", "futility", "beta spending", "two means", "data monitoring", "stopping rules", "wang-tsiatis", "hwang-shih-decani"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Maximum Total)"}, {"power", "Power", "power", "Power"}};
  p.params = gs_params("One-sided (H1: μ1 > μ2)", "One-sided (H1: μ1 < μ2)");
  for (auto& q : two_group_params(2, "50")) p.params.push_back(q);
  p.params.push_back(num("mu1", "μ1 (Mean of Group 1)", "1")
                         .grp(G_EFFECT)
                         .sym("μ1")
                         .sugg({"1", "0.5 1"})
                         .hint("Mean of group 1 (e.g. treatment) under the alternative."));
  p.params.push_back(num("mu2", "μ2 (Mean of Group 2)", "0")
                         .grp(G_EFFECT)
                         .sym("μ2")
                         .hint("Mean of group 2 (e.g. control); μ1 − μ2 is the difference to detect."));
  p.params.push_back(num("sigma", "σ (Common Standard Deviation)", "2")
                         .grp(G_EFFECT)
                         .sym("σ")
                         .above(0)
                         .sugg({"1", "2"})
                         .hint("Standard deviation of the outcome within each group, assumed the same in both "
                               "groups."));
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double diff = s("mu1") - s("mu2"), sigma = s("sigma");
    if (diff == 0) throw CalcError("μ1 equals μ2");
    auto eta_of = [&](double n1, double n2) { return directed(s, diff / (sigma * std::sqrt(1 / n1 + 1 / n2))); };
    if (eta_of(10, 10) < 0) throw CalcError("the difference has the sign opposite to the one-sided alternative");
    const DesignPtr info = gs_design(s, ctx);
    TwoN n;
    double eta = NAN;
    if (s.solve == "n") {
      check_power_target(s);
      eta = target_drift(s, *info, 1.0);
      n = solve_two_group_n(s, eta_of, eta, 2, 2 * sigma * sigma * eta * eta / (diff * diff), true, ctx);
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    const double sg = direction(s, diff);
    gs_report(s, r, *info, eta_of(n.n1, n.n2), eta, 1.0, n.total(),
              [&](std::size_t k, double z) { return sg * z * sigma * std::sqrt((1 / n.n1 + 1 / n.n2) / info->d.t[k]); },
              "N at Analysis", "(μ1 − μ2)");
    r.set("mu1", s("mu1"));
    r.set("mu2", s("mu2"));
    r.set("sigma", sigma);
  };
  p.columns = [](const Scenario& s) {
    return gs_columns(s, {col("mu1", "μ1", "num"), col("mu2", "μ2", "num"), col("sigma", "σ", "num")});
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A group-sequential design with " + design_words(s, r) + " and a maximum of " + fint(r.num("n1")) + " and " +
           fint(r.num("n2")) + " subjects per group (" + fint(r.num("n")) + " in total) " + achieves_singular(s, r) +
           " to detect a difference of " + fnum(r.num("mu1") - r.num("mu2")) + " between the means when σ = " +
           fnum(r.num("sigma")) + ". The expected total sample size is " + fnum(r.num("en1"), 1) + " under H1 and " +
           fnum(r.num("en0"), 1) + " under H0.";
  };
  p.definitions = GS_DEFS;
  p.definitions.push_back({"Efficacy / Futility (μ1 − μ2)", "Observed difference in means at a boundary."});
  p.references = GS_REFS;
  p.notes = {"The test statistic is the z-test with known σ; with σ estimated, a t-test needs a few more subjects "
             "when the sample size is small."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Maximum, Group 1)";
  p.example_solve = "n";
  p.example = {{"looks", "3"}, {"power", "0.9"}, {"mu1", "1"}, {"mu2", "0"}, {"sigma", "2"}};
  return p;
}
SSC_REGISTER(gs_means);

// ======================================================================= two proportions
Procedure gs_props() {
  Procedure p;
  p.id = "gs.props";
  p.name = "Group-Sequential Tests for Two Proportions";
  p.category = "Group-Sequential";
  p.subcategory = "Two Proportions";
  p.order = 1;
  p.summary =
      "Group-sequential comparison of two proportions with the z-test of the difference (pooled or unpooled "
      "variance), alpha-spending or Wang–Tsiatis efficacy boundaries and optional futility boundaries; reports "
      "maximum and expected sample sizes.";
  p.keywords = {"group sequential", "interim analysis", "two proportions", "alpha spending", "futility",
                "binary outcome", "data monitoring", "stopping rules"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Maximum Total)"}, {"power", "Power", "power", "Power"}};
  p.params = gs_params("One-sided (H1: P1 > P2)", "One-sided (H1: P1 < P2)");
  p.params.push_back(select("stat", "Test Statistic",
                            {{"unpooled", "Z-test, unpooled variance", {}}, {"pooled", "Z-test, pooled variance", {}}},
                            "unpooled")
                         .grp(G_TEST)
                         .hint("Pooled: variance estimated under H0 (equivalent to the chi-square test); unpooled: "
                               "separate variances (Wald)."));
  for (auto& q : two_group_params(2, "100")) p.params.push_back(q);
  p.params.push_back(num("p1", "P1 (Group 1 Proportion)", "0.5")
                         .grp(G_EFFECT)
                         .sym("P1")
                         .range(0, 1, true, true)
                         .sugg({"0.5", "0.4 0.5"})
                         .hint("The share of subjects in group 1 (e.g. treatment) expected to have the outcome under "
                               "the alternative, between 0 and 1."));
  p.params.push_back(num("p2", "P2 (Group 2 Proportion)", "0.35")
                         .grp(G_EFFECT)
                         .sym("P2")
                         .range(0, 1, true, true)
                         .sugg({"0.35"})
                         .hint("The share of subjects in group 2 (e.g. control) expected to have the outcome, between "
                               "0 and 1."));
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double p1 = s("p1"), p2 = s("p2"), diff = p1 - p2;
    if (diff == 0) throw CalcError("P1 equals P2");
    const bool pooled = s.is("stat", "pooled");
    // standard errors under H1 and (pooled) under H0 for group sizes n1, n2
    auto se1 = [&](double n1, double n2) { return std::sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2); };
    auto se0 = [&](double n1, double n2) {
      const double pb = (n1 * p1 + n2 * p2) / (n1 + n2);
      return std::sqrt(pb * (1 - pb) * (1 / n1 + 1 / n2));
    };
    // the statistic is (P̂1 − P̂2)/SE; dividing by SE1/SE gives a unit-variance statistic with drift diff/SE1
    auto scale_of = [&](double n1, double n2) { return pooled ? se0(n1, n2) / se1(n1, n2) : 1.0; };
    auto eta_of = [&](double n1, double n2) { return directed(s, diff / se1(n1, n2)); };
    if (eta_of(10, 10) < 0) throw CalcError("the difference has the sign opposite to the one-sided alternative");
    const DesignPtr info = gs_design(s, ctx);
    TwoN n;
    double eta = NAN;
    if (s.solve == "n") {
      check_power_target(s);
      const double v = p1 * (1 - p1) + p2 * (1 - p2);
      if (pooled) {
        eta = target_drift(s, *info, scale_of(100, 100));
        n = solve_two_group_n(
            s, [&](double n1, double n2) { return gs::exit_probs(info->d, eta_of(n1, n2), scale_of(n1, n2)).reject; },
            s("power"), 2, v * eta * eta / (diff * diff), true, ctx);
      } else {
        eta = target_drift(s, *info, 1.0);
        n = solve_two_group_n(s, eta_of, eta, 2, v * eta * eta / (diff * diff), true, ctx);
      }
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    const double sg = direction(s, diff);
    if (pooled && std::isfinite(eta)) eta = target_drift(s, *info, scale_of(n.n1, n.n2));
    gs_report(s, r, *info, eta_of(n.n1, n.n2), eta, scale_of(n.n1, n.n2), n.total(),
              [&](std::size_t k, double z) {
                const double t = info->d.t[k];
                return sg * z * (pooled ? se0(n.n1 * t, n.n2 * t) : se1(n.n1 * t, n.n2 * t));
              },
              "N at Analysis", "(P1 − P2)");
    r.set("p1", p1);
    r.set("p2", p2);
  };
  p.columns = [](const Scenario& s) { return gs_columns(s, {col("p1", "P1", "num"), col("p2", "P2", "num")}); };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A group-sequential design with " + design_words(s, r) + " and a maximum of " + fint(r.num("n1")) + " and " +
           fint(r.num("n2")) + " subjects per group " + achieves_singular(s, r) + " to detect P1 = " + fnum(r.num("p1")) +
           " versus P2 = " + fnum(r.num("p2")) + " with the " +
           (s.is("stat", "pooled") ? "pooled" : "unpooled") + " z-test. The expected total sample size is " +
           fnum(r.num("en1"), 1) + " under H1 and " + fnum(r.num("en0"), 1) + " under H0.";
  };
  p.definitions = GS_DEFS;
  p.definitions.push_back({"Efficacy / Futility (P1 − P2)", "Observed difference P̂1 − P̂2 at a boundary (at the "
                                                            "anticipated variance)."});
  p.references = GS_REFS;
  p.notes = {"The z-test uses the normal approximation. With small proportions or small groups its actual type I "
             "error can exceed α (P1 = 0.05, P2 = 0.15, three Pocock analyses at one-sided α = 0.025: 0.030 by "
             "simulation)."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Maximum, Group 1)";
  p.example_solve = "n";
  p.example = {{"looks", "3"}, {"power", "0.8"}, {"p1", "0.5"}, {"p2", "0.35"}};
  return p;
}
SSC_REGISTER(gs_props);

// ======================================================================= logrank
double gs_event_prob(double lam, double A, double F) {
  if (A <= 0) return -std::expm1(-lam * F);
  // without cancellation for small hazards, as event_prob of the logrank procedures
  const double x = lam * A;
  const double tail = x < 1e-3 ? x / 2 - x * x / 6 + x * x * x / 24 - x * x * x * x / 120 : (x + std::expm1(-x)) / x;
  return -std::expm1(-lam * F) + std::exp(-lam * F) * tail;
}

Procedure gs_logrank() {
  Procedure p;
  p.id = "gs.logrank";
  p.name = "Group-Sequential Logrank Tests (Hazard Ratio)";
  p.category = "Group-Sequential";
  p.subcategory = "Survival";
  p.order = 1;
  p.summary =
      "Group-sequential logrank test planned on the number of events (Schoenfeld), with analyses at fractions of the "
      "maximum events and optional futility boundaries; subjects follow from exponential survival with uniform "
      "accrual and additional follow-up.";
  p.keywords = {"group sequential", "logrank", "survival", "hazard ratio", "events", "interim analysis", "futility", "interim survival analysis", "data monitoring"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Maximum Total)"}, {"power", "Power", "power", "Power"}};
  p.params = gs_params("One-sided (H1: HR > 1)", "One-sided (H1: HR < 1)");
  for (auto& q : two_group_params(1, "200")) p.params.push_back(q);
  for (auto q : {num("hr", "HR (Hazard Ratio, Group 1/Group 2)", "0.7")
                     .grp(G_EFFECT)
                     .sym("HR")
                     .above(0)
                     .sugg({"0.7", "0.6 0.7 0.8"})
                     .hint("Group 1's hazard (event rate) divided by group 2's; 0.7 means a 30% lower risk of the "
                           "event at any time in group 1."),
                 num("med2", "Median Survival in Group 2", "12")
                     .grp(G_EFFECT)
                     .sym("Median2")
                     .above(0)
                     .sugg({"12"})
                     .hint("Median survival time of group 2 (e.g. control), assuming exponential survival, in the "
                           "time unit of accrual and follow-up."),
                 num("accrual", "Accrual Time", "12")
                     .grp("Study Duration")
                     .sym("Accrual")
                     .at_least(0)
                     .sugg({"12"})
                     .hint("Length of the period over which subjects enter the study at an even rate; 0 if all "
                           "start together."),
                 num("follow", "Additional Follow-Up", "12")
                     .grp("Study Duration")
                     .sym("Follow-Up")
                     .at_least(0)
                     .sugg({"12"})
                     .hint("Time from the end of accrual to the final analysis; the last subject to enter is "
                           "followed for this long.")})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double hr = s("hr"), lam2 = std::log(2.0) / s("med2"), lam1 = hr * lam2;
    if (hr == 1) throw CalcError("HR = 1: there is no effect to detect");
    if (s("accrual") + s("follow") <= 0) throw CalcError("the study duration (accrual + follow-up) must be positive");
    const double pe1 = gs_event_prob(lam1, s("accrual"), s("follow")), pe2 = gs_event_prob(lam2, s("accrual"), s("follow"));
    auto events = [&](double n1, double n2) { return n1 * pe1 + n2 * pe2; };
    // Schoenfeld: the logrank statistic has drift √(D w1 w2)·ln HR (sign: HR > 1 gives positive values)
    auto eta_of = [&](double n1, double n2) {
      const double w1 = n1 / (n1 + n2);
      return directed(s, std::sqrt(events(n1, n2) * w1 * (1 - w1)) * std::log(hr));
    };
    if (eta_of(100, 100) < 0) throw CalcError("the hazard ratio is on the side opposite to the one-sided alternative");
    const DesignPtr info = gs_design(s, ctx);
    TwoN n;
    double eta = NAN;
    if (s.solve == "n") {
      check_power_target(s);
      eta = target_drift(s, *info, 1.0);
      const double ev = eta * eta / (0.25 * std::log(hr) * std::log(hr));
      n = solve_two_group_n(s, eta_of, eta, 1, ev / (pe1 + pe2), true, ctx);
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    const double dmax = events(n.n1, n.n2);
    const double w1 = n.n1 / n.total();
    const double sg = direction(s, std::log(hr));
    gs_report(s, r, *info, eta_of(n.n1, n.n2), eta, 1.0, dmax,
              [&](std::size_t k, double z) { return std::exp(sg * z / std::sqrt(dmax * info->d.t[k] * w1 * (1 - w1))); },
              "Events at Analysis", "(HR)");
    r.set("events", dmax);
    r.set("hr", hr);
    r.set("med2", s("med2"));
  };
  p.columns = [](const Scenario& s) {
    return gs_columns(s, {col("events", "Max Events", "num"), col("hr", "HR", "num"), col("med2", "Median2", "num")},
                      "E(Events|H1)", "E(Events|H0)");
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A group-sequential logrank design with " + design_words(s, r) + " requires up to " +
           fint(std::ceil(r.num("events") - 1e-9)) + " events from " + fint(r.num("n")) + " subjects (" +
           fint(r.num("n1")) + " and " + fint(r.num("n2")) + " per group) and " + achieves_singular(s, r) +
           " to detect a hazard ratio of " + fnum(r.num("hr")) + "; the expected number of events is " +
           fnum(r.num("en1"), 1) + " under H1 and " + fnum(r.num("en0"), 1) + " under H0.";
  };
  p.definitions = GS_DEFS;
  p.definitions.push_back({"Efficacy / Futility (HR)", "Observed hazard ratio at a boundary (Schoenfeld's "
                                                       "approximation)."});
  p.references = GS_REFS;
  p.references.push_back("Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing "
                         "survival distributions. Biometrika 68, 316–319.");
  p.notes = {"The events come from Schoenfeld's formula. With equal groups it agrees with simulation; with unequal "
             "groups and a hazard ratio far from 1 it overstates the power when the larger group has the higher "
             "hazard (HR 0.65 with 1:2 allocation: 0.80, against 0.78 by simulation) and understates it in the "
             "opposite case."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Maximum, Group 1)";
  p.curve_min = 1;
  p.example_solve = "n";
  p.example = {{"looks", "3"}, {"power", "0.8 0.9"}, {"hr", "0.7"}, {"med2", "12"}};
  return p;
}
SSC_REGISTER(gs_logrank);

}  // namespace
}  // namespace ssc::procs
