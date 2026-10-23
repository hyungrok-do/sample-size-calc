// Linear mixed models y_i = X_i β + Z_i b_i + e_i, b_i ~ N(0, D), e_i ~ N(0, σ² I),
// fitted by restricted maximum likelihood (REML), for the simulation procedures.
//
// Simulated studies repeat a few designs many times (every subject with the
// same visits and arm, every cluster of the same size and arm), so the data
// are kept as *patterns*: one design (X, Z) with the sums of the responses
// over the clusters that share it.  Each cluster's y is split into its
// projection on the columns of Z, given by the least-squares coefficients
// b = (Z'Z)⁺ Z'y, and the rest, (I − P)y.  The REML criterion only needs, per
// design, Z'Z, Z'X and the within-cluster part X'(I − P)X and, per pattern,
// Σ b, Σ b b' and the within-cluster sums Σ X'(I − P)y and Σ y'(I − P)y, so
// one evaluation costs O(patterns), whatever the number of subjects.
//
// The two parts are kept apart, and the between-cluster part is taken in the
// coordinates b, because anything else loses digits when the clusters differ
// far more than their members do (an intercept SD 10⁵ times the residual SD,
// slopes that spread the subjects apart over a long follow-up): y'y, or Z'y,
// is then almost all between-cluster variation, and the residual sum of
// squares REML needs would be the small difference of large numbers.  With
// G = D/σ² and K = I + Z'Z G (q × q), V/σ² = I + Z G Z' has the determinant
// |K| and the inverse (I − P) + Z (Z'Z)⁺ W (Z'Z)⁺ Z', so that y'V⁻¹y σ² =
// y'(I − P)y + b'Wb, with W = K⁻¹ Z'Z (= ((Z'Z)⁻¹ + G)⁻¹ for a Z of full
// rank) computed entry by entry in closed form (q ≤ 2); both stay valid at
// the boundary (singular D).  σ² is profiled out; the relative
// covariance factor Λ (G = ΛΛ') is optimized: Brent's method for a random
// intercept, Nelder–Mead for a random intercept and slope.  Tests of c'β use
// a t statistic with Satterthwaite's degrees of freedom (as lmerTest does),
// from the numerical Hessian of the REML criterion in (σ², D).
#pragma once

#include <cstdint>
#include <vector>

namespace ssc::sim::lmm {

constexpr int PMAX = 12;  // fixed effects (a stepped wedge with 10 steps: intercept, 10 periods, treatment)
constexpr int QMAX = 2;  // random effects per cluster

/// All clusters that share the design matrices X (n × p) and Z (n × q).  P = Z (Z'Z)⁺ Z' is the projection on
/// the columns of one cluster's Z, b = (Z'Z)⁺ Z'y a cluster's coefficients on them (for a random intercept, the
/// cluster's mean), and X'(I − P)X and y'(I − P)y are what is left of X'X and y'y within the cluster.
struct Pattern {
  double count = 0;  // clusters with this design (an expected count in analytic approximations)
  std::int64_t n = 0;  // observations per cluster (a long has 32 bits in WebAssembly and on Windows)
  double xx[PMAX * PMAX] = {};  // X'(I − P)X, Z'Z, Z'X of one cluster (row-major)
  double zz[QMAX * QMAX] = {};
  double zx[QMAX * PMAX] = {};
  double zi[QMAX * QMAX] = {};  // (Z'Z)⁺ and X's coefficients on Z, C = (Z'Z)⁺ Z'X (q × p)
  double cx[QMAX * PMAX] = {};
  double xy[PMAX] = {};  // sums over the clusters: Σ X'(I − P)y, Σ b, Σ b b', Σ y'(I − P)y
  double b[QMAX] = {};
  double bb[QMAX * QMAX] = {};
  double yy = 0;

  /// Set the design from X (n × p) and Z (n × q), row-major.
  void set_design(long rows, int p, int q, const double* x, const double* z);
  /// Set (Z'Z)⁺ and C from Z'Z and Z'X: set_design does, and a design whose products are written directly must.
  void set_coefficients(int p, int q);
  /// Add one cluster's responses (n values); the design must be set.
  void add_cluster(int p, int q, const double* x, const double* z, const double* y);
  /// Add one cluster whose responses are Z u + e (u: q values, e: n values): the sums come from e alone and u,
  /// so that large random effects u cost no digits.
  void add_cluster(int p, int q, const double* x, const double* z, const double* u, const double* e);
  /// Add one cluster whose design gives the sufficient statistics directly: X'(I − P)y, b and y'(I − P)y.
  void add_sums(int p, int q, const double* xy_, const double* b_, double yy_);
  void clear_responses();
};

struct Data {
  int p = 0;
  int q = 0;
  std::vector<Pattern> patterns;
  double observations() const;
  double clusters() const;
};

struct Fit {
  bool ok = false;
  double theta[3] = {};  // Λ: q = 1 → θ0; q = 2 → [[θ0, 0], [θ1, θ2]]
  double g[QMAX * QMAX] = {};  // G = D/σ²
  double sigma2 = 0;
  double beta[PMAX] = {};
  double cov[PMAX * PMAX] = {};  // covariance of β̂
  double criterion = 0;          // −2 × REML log-likelihood
  int evaluations = 0;
  bool boundary = false;  // a variance component (or the correlation) at the edge of its range
};

/// REML fit; `theta0` (may be null) is the starting value of Λ.
Fit fit(const Data& d, const double* theta0 = nullptr);

struct Wald {
  double estimate = 0;
  double se = 0;
  double df = 0;
  double t = 0;
  bool fallback = false;  // Satterthwaite's df unavailable (boundary), `fallback_df` used
};

/// t test of c'β = 0; df_method: 0 Satterthwaite, 1 the given df, 2 normal (df = ∞).
Wald wald(const Data& d, const Fit& f, const double* c, int df_method, double fixed_df);

/// −2 × REML log-likelihood at (σ², G) with β profiled out; NaN if V is not positive definite.
double criterion(const Data& d, double sigma2, const double* g);

/// Covariance of β̂ at (σ², G) (the generalized least squares covariance), into `cov` (p × p).
bool beta_covariance(const Data& d, double sigma2, const double* g, double* cov);

}  // namespace ssc::sim::lmm
