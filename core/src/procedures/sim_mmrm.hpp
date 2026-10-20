// The mixed model for repeated measures (MMRM) of the simulation sim.mmrm, which the tests check directly:
// two groups, a mean for every group and visit, an unstructured covariance matrix common to the groups,
// and monotone missing data (a patient observed up to some visit, then not).
#pragma once

#include <array>

namespace ssc::procs::mmrm {

constexpr int TMAX = 12;  // visits

/// Sufficient statistics: for each group and last observed visit k (1 … T), the number of patients and
/// the sums of their responses and of the products of their responses.
struct Stats {
  int visits = 0;
  struct Pattern {
    double count = 0;
    std::array<double, TMAX> sum{};
    std::array<double, TMAX * TMAX> cross{};
  };
  std::array<std::array<Pattern, TMAX>, 2> pattern{};  // [group][k − 1]
  /// Add a patient of group g (0 or 1) observed at visits 1 … k (k ≥ 1) with responses y[0 … k − 1].
  void add(int g, int k, const double* y);
};

struct Fit {
  bool ok = false;
  double delta = 0;  // estimated mean of group 1 minus group 2 at the last visit
  double se = 0;     // its model-based standard error
  double observed_last = 0;  // patients observed at the last visit
  std::array<std::array<double, TMAX>, 2> mu{};  // estimated means by group and visit
  std::array<double, TMAX * TMAX> sigma{};      // estimated covariance matrix (row-major, visits × visits)
};

/// Maximum likelihood (reml = false) or restricted maximum likelihood (true) fit.  With monotone missing
/// data the likelihood factors into regressions of each visit on the earlier ones (Anderson 1957), which
/// give the covariance matrix in closed form; the means are then the generalized least squares estimates.
/// REML divides each regression's residual sum of squares by its patients minus 2 (the group means), which
/// is exact REML with complete data.
Fit fit(const Stats& s, bool reml = true);

}  // namespace ssc::procs::mmrm
