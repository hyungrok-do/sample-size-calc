// Parts of the non-proportional hazards simulation (sim_survival.cpp) that the
// tests check directly: the test statistics of one data set and the
// multivariate normal probability behind the MaxCombo p-value.
#pragma once

#include <vector>

namespace ssc::procs::nph {

/// One subject's observation: the time from entry to the event or to censoring.
struct Obs {
  double time;
  int group;  // 0: group 1 (treatment), 1: group 2 (control)
  int event;  // 1: event, 0: censored
};

/// Weighted logrank statistics for k ≤ 4 Fleming–Harrington weights Ŝ(t−)^ρ (1 − Ŝ(t−))^γ, Ŝ the pooled
/// Kaplan–Meier estimate: the scores U = Σ w (d1 − d·n1/n), observed minus expected events in group 1,
/// and their covariance matrix V (hypergeometric, row-major k × k).
struct WeightedLogrank {
  int k = 0;
  double u[4] = {0, 0, 0, 0};
  double v[16] = {};
};

/// `sorted`: the observations in increasing order of time; n1 and n2 subjects in groups 1 and 2.
WeightedLogrank weighted_logrank(const std::vector<Obs>& sorted, long n1, long n2, int k, const double* rho,
                                 const double* gamma);

/// Restricted mean survival time up to tau of each group (area under the Kaplan–Meier curve) and its
/// variance Σ A_j² d_j / (n_j (n_j − d_j)), A_j the area from the j-th event time to tau.
struct Rmst {
  double rmst[2] = {0, 0};
  double var[2] = {0, 0};
};

Rmst restricted_mean(const std::vector<Obs>& sorted, long n1, long n2, double tau);

/// P(a < X < b) for X ~ N(0, R), R an m × m correlation matrix (row-major), 1 ≤ m ≤ 4; the limits may be
/// infinite.  Deterministic: Genz's separation of variables with a fixed lattice rule (error about 1e-5),
/// whose points are shared between up to `threads` threads (the same result on any number).
double mvn_probability(int m, const double* a, const double* b, const double* r, int threads = 1);

/// P(X ∉ (a, b)) = 1 − P(a < X < b), also where it is far below the lattice's error: the sum over j of the
/// probability that X_j is the first variable outside its limits, each computed with X_j first, so that its tail
/// is a factor of full relative precision (and the rest is to about 1e-5 of the term).
double mvn_outside(int m, const double* a, const double* b, const double* r, int threads = 1);

}  // namespace ssc::procs::nph
