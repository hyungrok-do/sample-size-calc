// Linear mixed models y_i = X_i β + Z_i b_i + e_i, b_i ~ N(0, D), e_i ~ N(0, σ² I),
// fitted by restricted maximum likelihood (REML), for the simulation procedures.
//
// Simulated studies repeat a few designs many times (every subject with the
// same visits and arm, every cluster of the same size and arm), so the data
// are kept as *patterns*: one design (X, Z) with the sums of the responses
// over the clusters that share it.  The REML criterion only needs X'X, Z'Z,
// Z'X per design and Σ X'y, Σ Z'y, Σ (Z'y)(Z'y)', Σ y'y per pattern, so one
// evaluation costs O(patterns), whatever the number of subjects.
//
// With G = D/σ² and K = I + Z'Z G (q × q), the inverse and determinant of
// V/σ² = I + Z G Z' follow from V⁻¹σ² = I − Z G K⁻¹ Z' and |V/σ²| = |K|, which
// stay valid at the boundary (singular D).  σ² is profiled out; the relative
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

/// All clusters that share the design matrices X (n × p) and Z (n × q).
struct Pattern {
  double count = 0;  // clusters with this design (an expected count in analytic approximations)
  std::int64_t n = 0;  // observations per cluster (a long has 32 bits in WebAssembly and on Windows)
  double xx[PMAX * PMAX] = {};  // X'X, Z'Z, Z'X of one cluster (row-major)
  double zz[QMAX * QMAX] = {};
  double zx[QMAX * PMAX] = {};
  double xy[PMAX] = {};  // sums over the clusters: Σ X'y, Σ Z'y, Σ (Z'y)(Z'y)', Σ y'y
  double zy[QMAX] = {};
  double zyzy[QMAX * QMAX] = {};
  double yy = 0;

  /// Set the design from X (n × p) and Z (n × q), row-major.
  void set_design(long rows, int p, int q, const double* x, const double* z);
  /// Add one cluster's responses (n values).
  void add_cluster(int p, int q, const double* x, const double* z, const double* y);
  /// Add one cluster whose design gives the sufficient statistics directly.
  void add_sums(int p, int q, const double* xy_, const double* zy_, double yy_);
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
