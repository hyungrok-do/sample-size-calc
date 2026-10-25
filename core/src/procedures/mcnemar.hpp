// McNemar's test for paired yes/no outcomes, shared by the paired-proportions procedure and the comparison of two
// diagnostic tests on the same subjects.
#pragma once

#include <cmath>
#include <string>
#include <unordered_map>

#include "engine/procedure.hpp"

namespace ssc::procs {

/// Rejection probabilities of the exact conditional test by number of discordant pairs, which do not depend on the
/// number of pairs: a sample-size search computes each once.
using McNemarCache = std::unordered_map<int, double>;

/// Power of McNemar's test with n pairs whose discordant cells have probabilities p10 and p01, for the scenario's
/// "alternative" (upper: p10 > p01), "alpha" and "method" ("exact" conditional test, or Connor's "normal"
/// approximation).  q10 and q01, when given, are 1 − p10 and 1 − p01 (for a p10 that rounds to 1).  spread, when
/// given, receives how far the power can be from the formula's for the rounding of the normal approximation's
/// numerator: of its critical value, of √ψ and of √n·(P10 − P01) (0 for the exact test).
double mcnemar_power(const Scenario& s, int n, double p10, double p01, McNemarCache& cache, const RunContext& ctx,
                     double q10 = NAN, double q01 = NAN, double* spread = nullptr);

/// Why the normal approximation's power cannot be shown, for the cells p10 and p01 (with q10 = 1 − p10 and
/// q01 = 1 − p01) and the spread that mcnemar_power gave: a pair's difference with next to no variance, and how far
/// the rounding of the statistic can move the power.
std::string mcnemar_rounding(double p10, double p01, double q10, double q01, double spread);

}  // namespace ssc::procs
