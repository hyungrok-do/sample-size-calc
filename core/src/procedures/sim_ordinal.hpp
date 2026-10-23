// The proportional odds fit of the ordinal simulation (sim_ordinal.cpp), which the tests check directly.
#pragma once

namespace ssc::procs::ordinal {

constexpr int CMAX = 20;  // categories

struct PoFit {
  bool ok = false;
  bool separated = false;  // complete separation: the estimate is infinite
  bool single = false;     // every outcome in one category (or a group without subjects): nothing to fit; a fit
                           // that is not ok otherwise did not converge
  double beta = 0, se = 0, loglik = 0, loglik0 = 0;
};

/// Proportional odds fit of a 2 × k table (counts of group 1 and group 2 by category, lowest first):
/// logit P(Y ≤ j | x) = θ_j − βx with x = 1 in group 1, so β > 0 when group 1 tends to higher categories.
/// Empty categories are merged away; Newton's method with step halving on the concave log-likelihood.
/// loglik0 is the log-likelihood without the group effect.
PoFit po_fit(const double* group1, const double* group2, int k);

}  // namespace ssc::procs::ordinal
