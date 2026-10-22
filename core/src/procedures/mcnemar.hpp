// McNemar's test for paired yes/no outcomes, shared by the paired-proportions procedure and the comparison of two
// diagnostic tests on the same subjects.
#pragma once

#include <unordered_map>

#include "engine/procedure.hpp"

namespace ssc::procs {

/// Rejection probabilities of the exact conditional test by number of discordant pairs, which do not depend on the
/// number of pairs: a sample-size search computes each once.
using McNemarCache = std::unordered_map<int, double>;

/// Power of McNemar's test with n pairs whose discordant cells have probabilities p10 and p01, for the scenario's
/// "alternative" (upper: p10 > p01), "alpha" and "method" ("exact" conditional test, or Connor's "normal"
/// approximation).
double mcnemar_power(const Scenario& s, int n, double p10, double p01, McNemarCache& cache, const RunContext& ctx);

}  // namespace ssc::procs
