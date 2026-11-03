// Group-sequential procedures for two means, two proportions and the logrank
// test, built on the boundary/exit-probability engine in groupseq.cpp.
#include <map>
#include <memory>
#include <set>

#include "procedures/apart.hpp"
#include "procedures/common.hpp"
#include "procedures/derive_misc.hpp"
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

// The futility values of a fixed futility rule, one per interim analysis.
std::vector<double> gs_futility_z(const Scenario& s, std::size_t interims) {
  std::vector<double> fz = s.v("fz");
  if (fz.size() == 1 && interims > 1) fz.assign(interims, fz[0]);
  if (fz.size() != interims)
    throw CalcError("enter " + std::to_string(interims) +
                    " futility values (one per interim analysis) or a single value");
  return fz;
}

// What a design depends on, as the key of the design cache (and of the cost of a request).
std::string design_key(const Scenario& s) {
  const auto t = gs_fractions(s);
  const double alpha = s("alpha");
  const bool two = two_sided(s);
  const std::string& b = s.s("boundary");
  const std::string fk = futility_kind(s);
  const std::vector<double> fz = fk == "z" ? gs_futility_z(s, t.size() - 1) : std::vector<double>{};
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
  return key;
}

gs::Efficacy gs_efficacy(const Scenario& s) {
  const std::string& b = s.s("boundary");
  gs::Efficacy eff;
  if (b == "wt_obf" || b == "wt_pocock" || b == "wt") {
    eff.kind = gs::Efficacy::Kind::WangTsiatis;
    eff.delta = b == "wt_obf" ? 0.0 : (b == "wt_pocock" ? 0.5 : s("wtdelta"));
  } else {
    eff.spend = gs::alpha_spending(b, s("alpha"), two_sided(s), spending_param(b, s, "rho", "gamma"));
  }
  return eff;
}

// β-spending futility boundaries of a scenario that has them (β = 1 − power).
gs::Futility beta_futility(const Scenario& s) {
  const std::string& bs = s.s("bspend");
  gs::Futility fut;
  fut.kind = gs::Futility::Kind::BetaSpending;
  fut.beta = 1 - s("power");
  fut.spend = gs::spending(bs, fut.beta, spending_param(bs, s, "brho", "bgamma"));
  return fut;
}

DesignPtr gs_design(const Scenario& s, const RunContext& ctx) {
  const std::string key = design_key(s);  // which checks the inputs
  const auto t = gs_fractions(s);
  const double alpha = s("alpha");
  const bool two = two_sided(s);
  const std::string fk = futility_kind(s);
  const std::vector<double> fz = fk == "z" ? gs_futility_z(s, t.size() - 1) : std::vector<double>{};
  return cached_design(
      key,
      [&] {
        const gs::Efficacy eff = gs_efficacy(s);
        gs::Futility fut;
        if (fk == "beta") {
          fut = beta_futility(s);
        } else if (fk == "z") {
          fut.kind = gs::Futility::Kind::Fixed;
          fut.z = fz;
        }
        fut.binding = s.is("binding", "binding");
        DesignInfo info;
        info.d = gs::make_design(t, alpha, two, eff, fut, [&] { ctx.check(); });
        info.h0 = gs::exit_probs(info.d, 0.0, 1, 0, [&] { ctx.check(); });
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
  const double one_tail = scale * stats::qnorm(two_sided(s) ? half_alpha(a) : a, false) + stats::qnorm(power);
  if (!two_sided(s)) return one_tail;
  const double c = scale * stats::qnorm(half_alpha(a), false);
  auto excess = [&](double eta) { return stats::pnorm(eta - c) + stats::pnorm(-eta - c) - power; };
  const double at_zero = excess(0);
  if (at_zero >= 0) return NAN;
  // the other region only adds power: at a tiny alpha, less than the rounding error of the one region's
  const double at_one = excess(one_tail);
  if (!(at_one >= 0)) return one_tail;
  return numerics::brent(excess, 0, one_tail, at_zero, at_one, 1e-13);
}

// Drift needed for the target power.  Stop is checked at every analysis of every design it tries: with a tiny alpha
// or a flat spending function, the finest grid takes some 10 ms an analysis, a second a solve (main checked only
// between solves, and took up to 10 s to stop a list of such designs).  Kept by design, power and scale (like the
// designs), so that a row's derivation shows the drift without searching for it again.
double target_drift(const Scenario& s, const DesignInfo& info, double scale, const RunContext& ctx) {
  if (scale == 1 && std::isfinite(info.d.eta_design)) return info.d.eta_design;
  static OnceCache<std::string, double> cache(256);
  const std::string key =
      design_key(s) + "|" + text::format_roundtrip(s("power")) + "|" + text::format_roundtrip(scale);
  return *cache.get(
      key, [&] { return gs::drift_for_power(info.d, s("power"), scale, 0, [&] { ctx.check(); }); },
      [&] { ctx.check(); });
}

void check_power_target(const Scenario& s) {
  if (!(s("power") > s("alpha"))) throw CalcError("the power must exceed alpha");
}

// Common outputs: power, expected sizes, inflation and the boundary tables.
// `effect(k, z)` converts a z boundary of analysis k to the effect scale.
// `eta_target` is the drift solved for the target power (NaN when solving for power).
void gs_report(const Scenario& s, Row& r, const DesignInfo& info, double eta, double eta_target, double scale,
               double nmax, const std::function<double(std::size_t, double)>& effect, const std::string& unit,
               const std::string& effect_label, const RunContext& ctx) {
  const gs::Design& d = info.d;
  const auto e1 = gs::exit_probs(d, eta, scale, 0, [&] { ctx.check(); });
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
    row.push_back(Cell(stats::clamp01(cum)));  // the quadrature's error can take it past 1, as the power (exit_probs)
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

// ------------------------------------------------------------------ derivations
// The spending function `name` of the level `lvl` (a symbol), as gs::spending computes it, with its parameter (an
// input) in place.
std::string spending_text(const std::string& name, const std::string& lvl, double param) {
  using namespace derive;
  if (name == "obf") return lvl + "(t) = 2 − 2Φ(Φ⁻¹(1 − " + lvl + "/2) / √t), O'Brien–Fleming type";
  if (name == "pocock") return lvl + "(t) = " + lvl + " ln(1 + (e − 1)t), Pocock type";
  if (name == "power") return lvl + "(t) = " + lvl + " t^ρ (ρ = " + gv(param) + ")";
  if (std::fabs(param) < 1e-10) return lvl + "(t) = " + lvl + " t (Hwang–Shih–DeCani, γ = " + gv(param) + ")";
  return lvl + "(t) = " + lvl + " (1 − e^(−γt)) / (1 − e^(−γ)) (Hwang–Shih–DeCani, γ = " + gv(param) + ")";
}

// How a design's boundaries came about: the efficacy boundaries from their spending function (or the Wang–Tsiatis
// family), and the futility boundaries.
void gs_boundary_steps(derive::Steps& out, const Scenario& s, const gs::Design& d) {
  using namespace derive;
  const std::string& b = s.s("boundary");
  const std::string cross = d.two_sided ? "|Z_k| ≥ b_k" : "Z_k ≥ b_k";
  // the fractions are k/K, or those given rescaled so that the last is 1
  const std::string tb = "t = " + dvlist(d.t) + "; b = " + dvlist(d.b) + "; cumulative α = " +
                         dvlist(d.alpha_spent) + " (in all, α = " + gv(s("alpha")) + ")";
  if (b == "wt_obf" || b == "wt_pocock" || b == "wt") {
    // Δ is an input of the Wang–Tsiatis family, and 0 or ½ for its O'Brien–Fleming and Pocock members
    const std::string delta = b == "wt" ? gv(s("wtdelta")) : b == "wt_pocock" ? "½" : "0";
    out.push_back({"Efficacy boundaries",
                   "reject at analysis k when " + cross + ", b_k = C t_k^(Δ − ½) (Wang–Tsiatis, Δ = " + delta +
                       "), C such that P_H0(reject at some analysis) = α",
                   "C = b_K = " + dv(d.b.back()) + "; " + tb});
  } else {
    // a two-sided design spends α/2 in each tail with the one-sided function (gs::alpha_spending)
    const std::string lvl = d.two_sided ? "α′" : "α";
    const std::string inc = lvl + "(t_k) − " + lvl + "(t_{k−1})";
    out.push_back({"Efficacy boundaries",
                   "reject at analysis k when " + cross + ", b_k such that P_H0(no earlier stop, " + cross + ") = " +
                       (d.two_sided ? "2(" + inc + ")" : inc) + ", " +
                       spending_text(b, lvl, spending_param(b, s, "rho", "gamma")) +
                       (d.two_sided ? ", α′ = α/2 for each tail" : ""),
                   tb});
  }
  if (!d.has_futility()) return;
  const std::vector<double> a(d.a.begin(), d.a.end() - 1);  // the last is b_K
  const std::string bind = d.binding ? "; binding: the efficacy boundaries spend α among the paths above them"
                                     : "; non-binding: α is spent as if there were none";
  if (!d.beta_spent.empty()) {
    const std::string& bs = s.s("bspend");
    const std::vector<double> spent(d.beta_spent.begin(), d.beta_spent.end() - 1);
    out.push_back({"Futility boundaries",
                   "stop for futility at interim analysis k when Z_k ≤ a_k, a_k such that P_η(no earlier stop, Z_k ≤ "
                   "a_k) = β(t_k) − β(t_{k−1}) under the drift η at which the boundaries spend β = 1 − power in all, " +
                       spending_text(bs, "β", spending_param(bs, s, "brho", "bgamma")) + bind,
                   "β = 1 − " + gv(s("power")) + " = " + dv(1 - s("power")) + ", η = " + dv(d.eta_design) +
                       "; a = " + dvlist(a) + "; cumulative β = " + dvlist(spent)});
  } else {
    std::string z;
    for (double v : s.v("fz")) z += (z.empty() ? "" : ", ") + gv(v);
    out.push_back({"Futility boundaries",
                   "stop for futility at interim analysis k when Z_k ≤ a_k = min(z_k, b_k), z_k the values given (one "
                   "for all when one is given)" + bind,
                   "z = " + z + "; a = " + dvlist(a)});
  }
}

// The drift for the target power, when solving for N (with the boundaries × the scale s of a pooled statistic).
Step gs_target_drift_step(const Scenario& s, const gs::Design& d, double eta, double scale) {
  using namespace derive;
  if (scale == 1 && std::isfinite(d.eta_design))
    return {"Drift for the target power",
            "η = the drift (the mean of Z_K) at which the β-spending boundaries spend β = 1 − power in all, where the "
            "design's power is 1 − β",
            "1 − β = " + gv(s("power")) + " at η = " + dv(eta)};
  return {"Drift for the target power",
          "η = the drift (the mean of Z_K) at which Σₖ P_η(first crossing of the efficacy boundary at analysis k) = "
          "1 − β, the target power" +
              std::string(scale == 1 ? "" : ", with the boundaries × s of the sizes found (by their power) below"),
          "1 − β = " + gv(s("power")) + " at η = " + dv(eta)};
}

// What a drift stands for, the size that the inflation factor multiplies: the information I = (η/δ)², or the
// events D = η²/(w1 (1 − w1) θ²); and, where the allocation sets one from it, the number of subjects.
struct GsSize {
  std::string formula;                          // the size at the drift η_fixed, in symbols
  std::function<double(double)> value;          // the size at a drift
  std::function<std::string(double)> numbers;   // the formula with a drift's numbers
  std::function<std::string(double)> subjects;  // the subjects a size takes, in words and numbers ("" for none)
  std::string sym = "I";                        // the size's symbol
};

// The fixed design's drift for the power p and the inflation factor IF = (η/η_fixed)², as gs_report computes them
// (when the row has one): the fixed design's size, and the maximum size IF times it. Solving for N, p is the target
// 1 − β (given); solving for power, the power found (derived).
void gs_inflation_step(derive::Steps& out, const Scenario& s, const Row& r, double eta, double p, bool target,
                       double scale, const GsSize& size) {
  using namespace derive;
  if (!r.has("inflation")) return;
  const double f = fixed_drift(s, p, scale), inflation = eta * eta / (f * f), fixed = size.value(f);
  const std::string sz = scale == 1 ? "" : "s ", szn = scale == 1 ? "" : dv(scale) + " × ";
  const std::string p_sym = target ? "1 − β" : "P", p_text = target ? gv(p) : dvp(p);
  std::string law, numbers;
  if (two_sided(s)) {
    const double z = stats::qnorm(half_alpha(s("alpha")), false), c = scale * z;
    law = "η_fixed solves Φ(η − z) + Φ(−η − z) = " + p_sym + ", z = " + sz + "Φ⁻¹(1 − α/2)";
    numbers = (scale == 1 ? "" : "z = " + szn + dv(z) + " = " + dv(c) + "; ") + "Φ(" + dv(f) + " − " + dv(c) +
              ") + Φ(" + dv(-f) + " − " + dv(c) + ") = " + p_text + " at η_fixed = " + dv(f);
  } else {
    law = "η_fixed = " + sz + "Φ⁻¹(1 − α) + Φ⁻¹(" + p_sym + ")";
    numbers = "η_fixed = " + szn + "Φ⁻¹(1 − " + gv(s("alpha")) + ") + Φ⁻¹(" + p_text + ") = " + szn +
              dv(stats::qnorm(s("alpha"), false)) + " + " + dvs(stats::qnorm(p)) + " = " + dv(f);
  }
  const std::string fixed_subjects = size.subjects ? size.subjects(fixed) : "";
  out.push_back({"Fixed design and inflation",
                 "IF = (η / η_fixed)², " + law + " the drift of a fixed design (one analysis) with the same power" +
                     (target ? "" : " P as this design") + "; the fixed design needs " + size.formula +
                     ", this design IF times as much",
                 numbers + "; IF = (" + dv(eta) + " / " + dv(f) + ")² = " + dv(inflation) + "; fixed design: " +
                     size.numbers(f) + " = " + dv(fixed) + (fixed_subjects.empty() ? "" : ", " + fixed_subjects) +
                     "; this design: " + size.sym + " = " + dv(inflation) + " × " + dv(fixed) + " = " +
                     dv(inflation * fixed)});
}

// The power: the chance of crossing an efficacy boundary at some analysis, under the drift eta (as gs_report).
Step gs_power_step(const gs::Design& d, double eta, double scale, double power) {
  using namespace derive;
  const gs::ExitProbs e = gs::exit_probs(d, eta, scale);
  std::string parts;
  for (std::size_t k = 0; k < d.t.size(); ++k) parts += (k ? " + " : "") + dvp(e.up[k] + e.lo[k]);
  return {"Power",
          "Σₖ P_η(first crossing of the efficacy boundary at analysis k)" +
              std::string(scale == 1 ? "" : ", with the boundaries × s") +
              (d.has_futility() ? ", a trial stopped for futility not rejecting" : ""),
          (d.t.size() > 1 ? parts + " = " : "") + dvp(power)};
}

// "δ = μ1 − μ2" toward the alternative: the difference itself, its negative or its size; and with the numbers a
// and b (given) in place, ending in its value.
std::string toward(const Scenario& s, const std::string& sym, const std::string& a, const std::string& b) {
  if (two_sided(s)) return sym + " = |" + a + " − " + b + "|";
  return s.is("alternative", "lower") ? sym + " = " + b + " − " + a : sym + " = " + a + " − " + b;
}
std::string toward_numbers(const Scenario& s, const std::string& sym, double a, double b, double value) {
  using namespace derive;
  if (two_sided(s)) return sym + " = |" + gv(a) + " − " + gvs(b) + "| = " + dv(value);
  return s.is("alternative", "lower") ? sym + " = " + gv(b) + " − " + gvs(a) + " = " + dv(value)
                                      : sym + " = " + gv(a) + " − " + gvs(b) + " = " + dv(value);
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

// The time of a request in units of a typical scenario (a few ms), for the limit of scenarios. Each design is
// built once: about 3.5 units with 20 analyses, 11 with β-spending futility (whose design depends on the power
// too); a power from it is quick, and a search for the sample size takes about 3.5 ms per analysis (80 ms with
// 20 analyses).  The work grows with the square of the grid's r, which is finer for a tiny alpha and for nearly flat
// spending functions (see gs::grid_resolution).
double gs_cost(const std::vector<Scenario>& scenarios) {
  std::set<std::string> designs;
  double cost = 0;
  for (const auto& s : scenarios) {
    const double looks = s("looks");
    double w = 1;
    try {
      const gs::Futility fut = futility_kind(s) == "beta" ? beta_futility(s) : gs::Futility{};
      const double grid = gs::grid_resolution(s("alpha"), gs_fractions(s), gs_efficacy(s), fut) / 32.0;
      w = grid * grid;
    } catch (const CalcError&) {
      // the scenario stops at once with its message
    }
    cost += w * (s.solve == "power" ? 1 : 1 + looks / 2);
    try {
      if (designs.insert(design_key(s)).second)
        cost += w * (futility_kind(s) == "beta" ? 1 + looks / 2 : 1 + looks / 8);
    } catch (const CalcError&) {
      // the scenario stops at once with its message
    }
  }
  return cost;
}

// ======================================================================= two means
// The drift of the comparison of two means with n1 and n2 subjects: the mean of its z statistic.
double means_drift(const Scenario& s, double diff, double sigma, double n1, double n2) {
  return directed(s, diff / (sigma * std::sqrt(1 / n1 + 1 / n2)));
}

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
    // with no effect the power is the design's type I error, which no sample size raises
    if (diff == 0 && s.solve == "n") throw CalcError("μ1 equals μ2");
    auto eta_of = [&](double n1, double n2) { return means_drift(s, diff, sigma, n1, n2); };
    if (eta_of(10, 10) < 0) throw CalcError("the difference has the sign opposite to the one-sided alternative");
    const DesignPtr info = gs_design(s, ctx);
    TwoN n;
    double eta = NAN;
    if (s.solve == "n") {
      check_power_target(s);
      eta = target_drift(s, *info, 1.0, ctx);
      n = solve_two_group_n(s, eta_of, eta, 2, 2 * sigma * sigma * eta * eta / (diff * diff), true, ctx);
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    const double sg = direction(s, diff);
    gs_report(s, r, *info, eta_of(n.n1, n.n2), eta, 1.0, n.total(),
              [&](std::size_t k, double z) { return sg * z * sigma * std::sqrt((1 / n.n1 + 1 / n.n2) / info->d.t[k]); },
              "N at Analysis", "(μ1 − μ2)", ctx);
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
  // Solving for N: the drift the design needs for the target power, the fixed design's size and the inflation, and
  // the drift and power at the sizes found and one design below; solving for power, the power at the sizes given.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const RunContext ctx;
    const DesignPtr info = gs_design(s, ctx);
    const double mu1 = r.num("mu1"), mu2 = r.num("mu2"), diff = mu1 - mu2, sigma = r.num("sigma");
    const blocks::TwoN n{r.num("n1"), r.num("n2")};
    const double delta = directed(s, diff), eta = means_drift(s, diff, sigma, n.n1, n.n2), pw = r.num("power");
    const double inf = 1 / (sigma * sigma * (1 / n.n1 + 1 / n.n2));
    const bool solve_n = s.solve == "n";
    const double target = solve_n ? target_drift(s, *info, 1.0, ctx) : NAN;
    const auto [N1, N2] = group_sizes(s, n, solve_n);
    const std::string sg = gv(sigma);
    // the fixed design's information, and the subjects it takes where group 2 follows group 1 (I = N1 / (σ²(1 + 1/r)))
    const std::string& alloc = s.s("alloc");
    const double k = alloc == "ratio" ? s("r") : 1;
    const GsSize size{"I = (η_fixed / δ)² of information",
                      [&](double e) { return (e / delta) * (e / delta); },
                      [&](double e) { return "I = (" + dv(e) + " / " + dv(delta) + ")²"; },
                      [&](double i) -> std::string {
                        if (alloc == "equal")
                          return "N1 = N2 = 2σ² I = 2 × " + sg + "² × " + dv(i) + " = " + dv(2 * sigma * sigma * i);
                        if (alloc == "ratio")
                          return "N1 = σ²(1 + 1/r) I = " + sg + "² × (1 + 1/" + gv(k) + ") × " + dv(i) + " = " +
                                 dv(sigma * sigma * (1 + 1 / k) * i);
                        return "";
                      }};
    Derivation out;
    out.find = solve_n ? find_size(two_group_unknown(s), s("power"))
                       : find_power("with " + group_sizes_text(s, n, false));
    gs_boundary_steps(out.steps, s, info->d);
    if (solve_n) {
      out.steps.push_back(gs_target_drift_step(s, info->d, target, 1.0));
      gs_inflation_step(out.steps, s, r, target, s("power"), true, 1.0, size);
    }
    out.steps.push_back({"Drift", "η = δ √I, I = 1 / (σ²(1/N1 + 1/N2)) the information, " + toward(s, "δ", "μ1", "μ2"),
                         toward_numbers(s, "δ", mu1, mu2, delta) + "; I = 1 / (" + sg + "² × (1/" + N1 + " + 1/" +
                             N2 + ")) = " + dv(inf) + "; η = " + dv(delta) + " × √" + dv(inf) + " = " + dv(eta)});
    out.steps.push_back(gs_power_step(info->d, eta, 1.0, pw));
    if (solve_n) {
      out.steps.push_back(two_group_search_of(
          s, [&](double x1, double x2) { return means_drift(s, diff, sigma, x1, x2); }, n, target, 2, "drift",
          "η"));
      out.answer =
          group_sizes_text(s, n, true) + " (N = " + dvi(n.total()) + "), with power " + dvp(pw, s("power"));
    } else {
      gs_inflation_step(out.steps, s, r, eta, pw, false, 1.0, size);
      out.answer = "power = " + dvp(pw);
    }
    return out;
  };
  p.definitions = GS_DEFS;
  p.definitions.push_back({"Efficacy / Futility (μ1 − μ2)", "Observed difference in means at a boundary."});
  p.references = GS_REFS;
  p.notes = {"The test statistic is the z-test with known σ; with σ estimated, a t-test needs a few more subjects "
             "when the sample size is small."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Maximum, Group 1)";
  p.example_solve = "n";
  p.cost = gs_cost;
  p.example = {{"looks", "3"}, {"power", "0.9"}, {"mu1", "1"}, {"mu2", "0"}, {"sigma", "2"}};
  return p;
}
SSC_REGISTER(gs_means);

// ======================================================================= two proportions
// Standard errors of P̂1 − P̂2 under H1 and (pooled) under H0 for group sizes n1, n2.
double props_se1(double p1, double p2, double n1, double n2) {
  return std::sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2);
}
double props_se0(double p1, double p2, double n1, double n2) {
  const double pb = (n1 * p1 + n2 * p2) / (n1 + n2);
  return std::sqrt(pb * (1 - pb) * (1 / n1 + 1 / n2));
}

// The statistic is (P̂1 − P̂2)/SE; dividing by SE1/SE gives a unit-variance statistic with drift diff/SE1, whose
// boundaries are multiplied by SE/SE1 (1 for the unpooled statistic).
double props_scale(bool pooled, double p1, double p2, double n1, double n2) {
  return pooled ? props_se0(p1, p2, n1, n2) / props_se1(p1, p2, n1, n2) : 1.0;
}
double props_drift(const Scenario& s, double p1, double p2, double n1, double n2) {
  return directed(s, (p1 - p2) / props_se1(p1, p2, n1, n2));
}

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
    // with no effect the power is the design's type I error, which no sample size raises
    if (diff == 0 && s.solve == "n") throw CalcError("P1 equals P2");
    const bool pooled = s.is("stat", "pooled");
    // standard errors under H1 and (pooled) under H0 for group sizes n1, n2
    auto se1 = [&](double n1, double n2) { return props_se1(p1, p2, n1, n2); };
    auto se0 = [&](double n1, double n2) { return props_se0(p1, p2, n1, n2); };
    auto scale_of = [&](double n1, double n2) { return props_scale(pooled, p1, p2, n1, n2); };
    auto eta_of = [&](double n1, double n2) { return props_drift(s, p1, p2, n1, n2); };
    if (eta_of(10, 10) < 0) throw CalcError("the difference has the sign opposite to the one-sided alternative");
    const DesignPtr info = gs_design(s, ctx);
    TwoN n;
    double eta = NAN;
    if (s.solve == "n") {
      check_power_target(s);
      const double v = p1 * (1 - p1) + p2 * (1 - p2);
      if (pooled) {
        eta = target_drift(s, *info, scale_of(100, 100), ctx);
        n = solve_two_group_n(
            s,
            [&](double n1, double n2) {
              return gs::exit_probs(info->d, eta_of(n1, n2), scale_of(n1, n2), 0, [&] { ctx.check(); }).reject;
            },
            s("power"), 2, v * eta * eta / (diff * diff), true, ctx);
      } else {
        eta = target_drift(s, *info, 1.0, ctx);
        n = solve_two_group_n(s, eta_of, eta, 2, v * eta * eta / (diff * diff), true, ctx);
      }
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    const double sg = direction(s, diff);
    if (pooled && std::isfinite(eta)) eta = target_drift(s, *info, scale_of(n.n1, n.n2), ctx);
    gs_report(s, r, *info, eta_of(n.n1, n.n2), eta, scale_of(n.n1, n.n2), n.total(),
              [&](std::size_t k, double z) {
                const double t = info->d.t[k];
                return sg * z * (pooled ? se0(n.n1 * t, n.n2 * t) : se1(n.n1 * t, n.n2 * t));
              },
              "N at Analysis", "(P1 − P2)", ctx);
    r.set("p1", p1);
    r.set("p2", p2);
  };
  p.columns = [](const Scenario& s) { return gs_columns(s, {col("p1", "P1", "num"), col("p2", "P2", "num")}); };
  p.statement = [](const Scenario& s, const Row& r) {
    const auto [p1, p2] = fnum_apart(r.num("p1"), r.num("p2"));  // 0.30001 versus 0.3, not 0.3 versus 0.3
    return "A group-sequential design with " + design_words(s, r) + " and a maximum of " + fint(r.num("n1")) + " and " +
           fint(r.num("n2")) + " subjects per group " + achieves_singular(s, r) + " to detect P1 = " + p1 +
           " versus P2 = " + p2 + " with the " +
           (s.is("stat", "pooled") ? "pooled" : "unpooled") + " z-test. The expected total sample size is " +
           fnum(r.num("en1"), 1) + " under H1 and " + fnum(r.num("en0"), 1) + " under H0.";
  };
  // Solving for N: the drift the design needs for the target power, the fixed design's size and the inflation, and
  // the drift and power at the sizes found and one design below; solving for power, the power at the sizes given.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const RunContext ctx;
    const DesignPtr info = gs_design(s, ctx);
    const double p1 = r.num("p1"), p2 = r.num("p2"), pw = r.num("power");
    const blocks::TwoN n{r.num("n1"), r.num("n2")};
    const double n1 = n.n1, n2 = n.n2;
    const bool pooled = s.is("stat", "pooled"), solve_n = s.solve == "n";
    const double delta = directed(s, p1 - p2), eta = props_drift(s, p1, p2, n1, n2);
    const double se1 = props_se1(p1, p2, n1, n2), scale = props_scale(pooled, p1, p2, n1, n2), inf = 1 / (se1 * se1);
    const double target = solve_n ? target_drift(s, *info, scale, ctx) : NAN;
    const auto [N1, N2] = group_sizes(s, n, solve_n);
    const std::string P1 = gv(p1), P2 = gv(p2);
    // the fixed design's information, and the subjects it takes where group 2 follows group 1
    // (I = N1 / (P1(1 − P1) + P2(1 − P2)/r))
    const std::string& alloc = s.s("alloc");
    const double k = alloc == "ratio" ? s("r") : 1, per = p1 * (1 - p1) + p2 * (1 - p2) / k;
    const std::string per_text = P1 + " × (1 − " + P1 + ") + " + P2 + " × (1 − " + P2 + ")" +
                                 (alloc == "ratio" ? " / " + gv(k) : std::string());
    const GsSize size{"I = (η_fixed / δ)² of information",
                      [&](double e) { return (e / delta) * (e / delta); },
                      [&](double e) { return "I = (" + dv(e) + " / " + dv(delta) + ")²"; },
                      [&](double i) -> std::string {
                        if (alloc != "equal" && alloc != "ratio") return "";
                        return std::string(alloc == "equal" ? "N1 = N2 = (P1(1 − P1) + P2(1 − P2)) I = ("
                                                            : "N1 = (P1(1 − P1) + P2(1 − P2)/r) I = (") +
                               per_text + ") × " + dv(i) + " = " + dv(per * i);
                      }};
    Derivation out;
    out.find = solve_n ? find_size(two_group_unknown(s), s("power"))
                       : find_power("with " + group_sizes_text(s, n, false));
    gs_boundary_steps(out.steps, s, info->d);
    if (pooled) {
      const double pb = (n1 * p1 + n2 * p2) / (n1 + n2), se0 = props_se0(p1, p2, n1, n2);
      out.steps.push_back({"Pooled statistic",
                           "s = SE0 / SE1, SE0 = √(P̄(1 − P̄)(1/N1 + 1/N2)) with P̄ = (N1 P1 + N2 P2) / (N1 + N2), SE1 "
                           "= √(P1(1 − P1)/N1 + P2(1 − P2)/N2): (P̂1 − P̂2)/SE0 has SD SE1/SE0 under H1, so the "
                           "boundaries are multiplied by s",
                           "P̄ = (" + N1 + " × " + P1 + " + " + N2 + " × " + P2 + ") / " + dvi(n1 + n2) + " = " +
                               dv(pb) + "; SE0 = " + dv(se0) + ", SE1 = " + dv(se1) + "; s = " + dv(scale)});
    }
    if (solve_n) {
      out.steps.push_back(gs_target_drift_step(s, info->d, target, scale));
      gs_inflation_step(out.steps, s, r, target, s("power"), true, scale, size);
    }
    out.steps.push_back(
        {"Drift",
         "η = δ √I, I = 1 / (P1(1 − P1)/N1 + P2(1 − P2)/N2) the information, " + toward(s, "δ", "P1", "P2"),
         toward_numbers(s, "δ", p1, p2, delta) + "; I = 1 / (" + P1 + " × (1 − " + P1 + ") / " + N1 + " + " + P2 +
             " × (1 − " + P2 + ") / " + N2 + ") = " + dv(inf) + "; η = " + dv(delta) + " × √" + dv(inf) + " = " +
             dv(eta)});
    out.steps.push_back(gs_power_step(info->d, eta, scale, pw));
    if (solve_n) {
      if (pooled)
        // with the scale changing with the sizes, the search reaches the target power itself
        out.steps.push_back(two_group_search(
            s,
            [&](double x1, double x2) {
              if (x1 == n1 && x2 == n2) return pw;
              return gs::exit_probs(info->d, props_drift(s, p1, p2, x1, x2), props_scale(true, p1, p2, x1, x2))
                  .reject;
            },
            n, s("power"), 2));
      else
        out.steps.push_back(two_group_search_of(
            s, [&](double x1, double x2) { return props_drift(s, p1, p2, x1, x2); }, n, target, 2, "drift", "η"));
      out.answer =
          group_sizes_text(s, n, true) + " (N = " + dvi(n.total()) + "), with power " + dvp(pw, s("power"));
    } else {
      gs_inflation_step(out.steps, s, r, eta, pw, false, scale, size);
      out.answer = "power = " + dvp(pw);
    }
    return out;
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
  p.cost = gs_cost;
  p.example = {{"looks", "3"}, {"power", "0.8"}, {"p1", "0.5"}, {"p2", "0.35"}};
  return p;
}
SSC_REGISTER(gs_props);

// ======================================================================= logrank
double gs_event_prob(double lam, double A, double F) {
  if (std::isinf(lam)) return 1;  // a median of about 0: each subject has the event at once
  if (A <= 0) return -std::expm1(-lam * F);
  // without cancellation for small hazards, as event_prob of the logrank procedures
  const double x = lam * A;
  const double tail = x < 1e-3   ? x / 2 - x * x / 6 + x * x * x / 24 - x * x * x * x / 120
                      : x < 1e300 ? (x + std::expm1(-x)) / x
                                  : 1;
  return -std::expm1(-lam * F) + std::exp(-lam * F) * tail;
}

// The hazards of an exponential survival in each group, and each group's chance of an event by the final analysis.
struct LogrankRates {
  double lam1, lam2, pe1, pe2;
};
LogrankRates logrank_rates(const Scenario& s) {
  const double hr = s("hr"), lam2 = std::log(2.0) / s("med2"), lam1 = hr * lam2;
  return {lam1, lam2, gs_event_prob(lam1, s("accrual"), s("follow")), gs_event_prob(lam2, s("accrual"), s("follow"))};
}

// The expected events of n1 and n2 subjects, and the logrank statistic's drift √(D w1 w2)·ln HR (Schoenfeld; HR > 1
// gives positive values).
double logrank_events(double pe1, double pe2, double n1, double n2) { return n1 * pe1 + n2 * pe2; }
double logrank_drift(const Scenario& s, double hr, double pe1, double pe2, double n1, double n2) {
  const double w1 = n1 / (n1 + n2);
  return directed(s, std::sqrt(logrank_events(pe1, pe2, n1, n2) * w1 * (1 - w1)) * std::log(hr));
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
    const double hr = s("hr");
    // with no effect the power is the design's type I error, which no sample size raises
    if (hr == 1 && s.solve == "n") throw CalcError("HR = 1: there is no effect to detect");
    if (s("accrual") + s("follow") <= 0) throw CalcError("the study duration (accrual + follow-up) must be positive");
    const LogrankRates rates = logrank_rates(s);
    const double pe1 = rates.pe1, pe2 = rates.pe2;
    auto events = [&](double n1, double n2) { return logrank_events(pe1, pe2, n1, n2); };
    // Schoenfeld: the logrank statistic has drift √(D w1 w2)·ln HR (sign: HR > 1 gives positive values)
    auto eta_of = [&](double n1, double n2) { return logrank_drift(s, hr, pe1, pe2, n1, n2); };
    if (eta_of(100, 100) < 0) throw CalcError("the hazard ratio is on the side opposite to the one-sided alternative");
    const DesignPtr info = gs_design(s, ctx);
    TwoN n;
    double eta = NAN;
    if (s.solve == "n") {
      check_power_target(s);
      eta = target_drift(s, *info, 1.0, ctx);
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
              "Events at Analysis", "(HR)", ctx);
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
  // Solving for N: the drift the design needs for the target power, the fixed design's events and the inflation,
  // the events the sizes found bring, and the drift and power there and one design below; solving for power, the
  // power at the sizes given.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const RunContext ctx;
    const DesignPtr info = gs_design(s, ctx);
    const LogrankRates rates = logrank_rates(s);
    const blocks::TwoN n{r.num("n1"), r.num("n2")};
    const double hr = r.num("hr"), n1 = n.n1, n2 = n.n2, A = s("accrual"), F = s("follow"), pw = r.num("power");
    const double D = logrank_events(rates.pe1, rates.pe2, n1, n2), w1 = n1 / (n1 + n2);
    const double theta = directed(s, std::log(hr)), eta = logrank_drift(s, hr, rates.pe1, rates.pe2, n1, n2);
    const bool solve_n = s.solve == "n";
    const double target = solve_n ? target_drift(s, *info, 1.0, ctx) : NAN;
    const auto [N1, N2] = group_sizes(s, n, solve_n);
    // the fixed design's events at this design's w1, and the subjects that bring them (D = N (w1 p1 + (1 − w1) p2))
    const double per = w1 * rates.pe1 + (1 - w1) * rates.pe2;
    const GsSize size{"D = η_fixed² / (w1 (1 − w1) θ²) events at this design's w1",
                      [&](double e) { return e * e / (w1 * (1 - w1) * theta * theta); },
                      [&](double e) {
                        return "D = " + dv(e) + "² / (" + dv(w1) + " × " + dv(1 - w1) + " × " + dv(theta) + "²)";
                      },
                      [&](double d) {
                        return "N = D / (w1 p1 + (1 − w1) p2) = " + dv(d) + " / (" + dv(w1) + " × " + dv(rates.pe1) +
                               " + " + dv(1 - w1) + " × " + dv(rates.pe2) + ") = " + dv(d / per);
                      },
                      "D"};
    Derivation out;
    out.find = solve_n ? find_size(two_group_unknown(s), s("power"))
                       : find_power("with " + group_sizes_text(s, n, false));
    gs_boundary_steps(out.steps, s, info->d);
    const std::string prob_f = A > 0 ? "p = 1 − e^(−λF)(1 − e^(−λA)) / (λA), the chance of an event with entry spread "
                                       "evenly over the accrual time A and follow-up F after it"
                                     : "p = 1 − e^(−λF), the chance of an event within the follow-up F";
    // each group's chance of an event, which turns events into subjects (and is needed before the fixed design's N)
    out.steps.push_back({"Chance of an event", prob_f + ", λ2 = ln 2 / median2, λ1 = HR λ2",
                         "λ2 = ln 2 / " + gv(s("med2")) + " = " + dv(rates.lam2) + ", λ1 = " + gv(hr) + " × " +
                             dv(rates.lam2) + " = " + dv(rates.lam1) + "; A = " + gv(A) + ", F = " + gv(F) +
                             ": p1 = " + dv(rates.pe1) + ", p2 = " + dv(rates.pe2)});
    if (solve_n) {
      out.steps.push_back(gs_target_drift_step(s, info->d, target, 1.0));
      gs_inflation_step(out.steps, s, r, target, s("power"), true, 1.0, size);
    }
    out.steps.push_back({"Expected events", "D = N1 p1 + N2 p2",
                         "D = " + N1 + " × " + dv(rates.pe1) + " + " + N2 + " × " + dv(rates.pe2) + " = " + dv(D)});
    const std::string th = two_sided(s) ? "θ = |ln HR|" : s.is("alternative", "lower") ? "θ = −ln HR" : "θ = ln HR";
    const std::string th_numbers = two_sided(s) ? "θ = |ln " + gv(hr) + "|"
                                   : s.is("alternative", "lower") ? "θ = −ln " + gv(hr)
                                                                  : "θ = ln " + gv(hr);
    out.steps.push_back({"Drift", "η = θ √(D w1 (1 − w1)), w1 = N1 / N, " + th + " (Schoenfeld)",
                         th_numbers + " = " + dv(theta) + "; w1 = " + N1 + " / " + dvi(n1 + n2) + " = " + dv(w1) +
                             "; η = " + dv(theta) + " × √(" + dv(D) + " × " + dv(w1) + " × " + dv(1 - w1) + ") = " +
                             dv(eta)});
    out.steps.push_back(gs_power_step(info->d, eta, 1.0, pw));
    if (solve_n) {
      out.steps.push_back(two_group_search_of(
          s, [&](double x1, double x2) { return logrank_drift(s, hr, rates.pe1, rates.pe2, x1, x2); }, n, target, 1,
          "drift", "η"));
      out.answer = group_sizes_text(s, n, true) + " (N = " + dvi(n.total()) + ", " + dv(D, 1) +
                   " events expected by the final analysis), with power " + dvp(pw, s("power"));
    } else {
      gs_inflation_step(out.steps, s, r, eta, pw, false, 1.0, size);
      out.answer = "power = " + dvp(pw);
    }
    return out;
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
  p.cost = gs_cost;
  p.example = {{"looks", "3"}, {"power", "0.8 0.9"}, {"hr", "0.7"}, {"med2", "12"}};
  return p;
}
SSC_REGISTER(gs_logrank);

}  // namespace
}  // namespace ssc::procs
