// Logistic random-intercept model for clustered yes/no outcomes:
//   s_c | b_c ~ Binomial(m_c, expit(β0 + β1·x_c + b_c)),  b_c ~ N(0, σ²),
// fitted by maximum likelihood with adaptive Gauss–Hermite quadrature (as
// glmer with nAGQ > 1, here 21 points by default), for the simulation procedures.
//
// With only cluster-level covariates a cluster enters the likelihood through
// (x, m, s) alone, so each one is a single one-dimensional integral over its
// random effect, written as ∫ exp(f(z)) φ(z) dz with b = σz.  The quadrature is
// centered at the mode of each cluster's integrand and scaled by its curvature.
// The likelihood is maximized by damped Newton steps on (β0, β1, σ), with the
// exact gradient of the quadrature (including the movement of its nodes with
// the parameters) and the Hessian of the sum with the nodes held fixed.  The
// likelihood is even in σ, so σ needs no constraint.
#pragma once

#include <vector>

namespace ssc::sim::glmm {

struct Cluster {
  double x = 0;  // group indicator (or another cluster-level covariate)
  double m = 0;  // size
  double s = 0;  // successes
  double count = 1;  // clusters with these values
};

/// Gauss–Hermite nodes and weights for the weight exp(−x²), n points (Newton's method on Hermite polynomials).
void gauss_hermite(int n, std::vector<double>& nodes, std::vector<double>& weights);

/// Log-likelihood of the clusters (without the binomial coefficients), with `nodes` quadrature points.
double loglik(const std::vector<Cluster>& data, double b0, double b1, double sigma, int nodes = 21);

/// Log of ∫ exp(s(η + σz) − m log(1 + e^(η+σz))) φ(z) dz for one cluster, by adaptive quadrature.
double cluster_loglik(double eta, double m, double s, double sigma, int nodes = 21);

/// Why a fit failed.
enum class Failure {
  None,
  RunOff,       // σ stayed past SIGMA_MAX or ended there: the subjects of every cluster agree, or see SIGMA_MAX
  Arithmetic,   // the likelihood could not be computed, or Newton's method stopped short of the maximum
  Information,  // the observed information at the maximum is not positive definite: no standard error
};

struct Fit {
  bool ok = false;
  Failure failure = Failure::None;
  double b0 = 0, b1 = 0, sigma = 0;  // σ ≥ 0
  double se_b1 = 0;                  // from the observed information
  double loglik = 0;
  int iterations = 0;
  bool boundary = false;  // σ estimated as 0, or held there
};

/// Called with `arg`, when set, as each pass over the clusters begins and every 4096 clusters within it: a caller whose
/// passes are long (a million clusters of different sizes, up to about a second each) checks there whether to stop,
/// and throws to stop the fit.
struct Poll {
  void (*call)(const void*) = nullptr;
  const void* arg = nullptr;
};

/// The largest σ a fit may end at.  The likelihood peaks at a finite σ unless the subjects of every cluster agree,
/// but where most of them do the quadrature can move the peak far past it: its error for a cluster whose subjects
/// agree grows with σ (for 1 to 5 subjects up to 0.01 of the cluster's log-likelihood at 6, the simulated designs'
/// largest σ, and 0.12 at 30 with 21 points; 0.0008 and 0.06 with 64), and summed over a few hundred such clusters
/// it moved the peak of one study's likelihood from 5.7 to about 32 with 21 points.
constexpr double SIGMA_MAX = 30;

/// Maximum likelihood fit from the starting values (b0, b1, sigma).  With one subject per cluster and two groups, σ
/// is not identified and is held at 0 (ordinary logistic regression).
Fit fit(const std::vector<Cluster>& data, double b0, double b1, double sigma, int nodes = 21, Poll poll = {});

/// E[expit(β + σZ)] for Z ~ N(0, 1): the population-averaged proportion of a cluster-specific β.
double marginal_mean(double beta, double sigma);
/// Var[expit(β + σZ)]: the between-cluster variance of the cluster proportions.
double between_variance(double beta, double sigma);
/// The intracluster correlation on the proportion scale, Var[expit(β + σZ)] / (p̄(1 − p̄)), which β and −β share:
/// taken where p̄ ≤ 1/2, so that 1 − p̄ keeps its digits next to 1.
double icc_of(double beta, double sigma);
/// β with E[expit(β + σZ)] = p.
double beta_for_marginal(double p, double sigma);
/// σ giving an intracluster correlation (on the proportion scale) of `icc` when the cluster-specific
/// intercept is `beta(σ)` (icc_of).  Throws std::domain_error when the largest power of 2 up to `limit` gives less.
template <class BetaOf>
double sigma_for_icc(double icc, BetaOf&& beta_of, double limit = 1e3);

}  // namespace ssc::sim::glmm

#include <cmath>
#include <stdexcept>

template <class BetaOf>
double ssc::sim::glmm::sigma_for_icc(double icc, BetaOf&& beta_of, double limit) {
  if (icc <= 0) return 0;
  auto f = [&](double sigma) { return icc_of(beta_of(sigma), sigma) - icc; };
  // the proportion-scale ICC grows with σ from 0 towards 1
  double lo = 0, hi = 1;
  while (f(hi) < 0) {
    lo = hi;
    hi *= 2;
    if (hi > limit) throw std::domain_error("intracluster correlation too large for the proportions");
  }
  for (int i = 0; i < 200 && hi - lo > 1e-13 * hi; ++i) {
    const double mid = 0.5 * (lo + hi);
    (f(mid) < 0 ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}
