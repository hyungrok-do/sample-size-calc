// Reusable building blocks for procedure definitions: standard parameters
// (alternative hypothesis, power, alpha, dropout, sample sizes, two-group
// allocation), the matching sample-size solvers and text formatting helpers
// for summary statements.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "engine/procedure.hpp"
#include "engine/solve.hpp"

namespace ssc::blocks {

// ------------------------------------------------------------------ groups
inline const char* G_TEST = "Test";
inline const char* G_POWER = "Power and Alpha";
inline const char* G_N = "Sample Size";
inline const char* G_EFFECT = "Effect Size";
inline const char* G_OPTIONS = "Options";

// ------------------------------------------------------------------ hypotheses
enum class Alt { Two, Upper, Lower };

/// Select "alternative" with two-sided / upper / lower choices.
ParamDef alternative(const std::string& two, const std::string& upper, const std::string& lower,
                     const std::string& def = "two");
Alt alt(const Scenario& s);
/// "two-sided" / "one-sided" wording.
std::string sidedness(Alt a);

// ------------------------------------------------------------------ standard params
ParamDef power_param(const std::string& def = "0.80");
ParamDef alpha_param(const std::string& def = "0.05", const std::string& label = "Alpha (Significance Level)");
ParamDef dropout_param();
/// A sample-size input that is hidden when solving for N.
ParamDef n_param(const std::string& id, const std::string& label, const std::string& def, double min_n);

/// Two-group allocation: "alloc", "n1", "n2", "r", "ntotal", "pct1".
std::vector<ParamDef> two_group_params(double min_n, const std::string& n1_default = "20");

struct TwoN {
  double n1 = 0;
  double n2 = 0;
  double total() const { return n1 + n2; }
};

/// Group sizes when they are inputs (solving for power or an effect).
TwoN two_group_n(const Scenario& s);

/// Smallest group sizes (per the allocation rule) reaching the target power.
/// `guess_n1` may be NaN.  `power(n1, n2)` must be safe for any n ≥ min_n.
TwoN solve_two_group_n(const Scenario& s, const std::function<double(double, double)>& power, double target,
                       double min_n, double guess_n1, bool monotone, const RunContext& ctx,
                       double max_n = 1e8);

void put_two_group(Row& r, const TwoN& n);
std::vector<Column> two_group_columns();

// ------------------------------------------------------------------ columns
Column col(const std::string& key, const std::string& label, const std::string& fmt, const std::string& tip = "");
/// Power columns: "Target Power" + "Actual Power" when solving for N, else "Power".
std::vector<Column> power_columns(const Scenario& s);

// ------------------------------------------------------------------ text
std::string fnum(double x, int max_decimals = 4);
std::string fpct(double p, int decimals = 1);
std::string fint(double x);
/// "power of 80.7%" style.
std::string power_text(double p);

}  // namespace ssc::blocks
