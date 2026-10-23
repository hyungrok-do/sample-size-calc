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
// one evaluation costs O(patterns), whatever the number of subjects.  Clusters
// whose X does not vary within the cluster, with a random intercept, need
// less, and are kept by arm and size (Arm).
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
#include <limits>
#include <vector>

namespace ssc::sim::lmm {

constexpr int PMAX = 12;  // fixed effects (a stepped wedge with 10 steps: intercept, 10 periods, treatment)
constexpr int QMAX = 2;  // random effects per cluster

/// All clusters that share the design matrices X (n × p) and Z (n × q).  P = Z (Z'Z)⁺ Z' is the projection on
/// the columns of one cluster's Z, b = (Z'Z)⁺ Z'y a cluster's coefficients on them (for a random intercept, the
/// cluster's mean), and X'(I − P)X and y'(I − P)y are what is left of X'X and y'y within the cluster.
///
/// Z's rows may be given in other coordinates than the random effects' (set_spanned): the random effects'
/// covariates are then Z S, for a q × q S with |S| = 1, where Z itself is far from singular and Z S need not be
/// (a random intercept given at a time far from two visits close together).  Z'Z, (Z'Z)⁺, C, b and the sums are
/// in the random effects' coordinates; the rows of Z, and the errors' coefficients on them, in Z's.
struct Pattern {
  double count = 0;  // clusters with this design (an expected count in analytic approximations)
  std::int64_t n = 0;  // observations per cluster (a long has 32 bits in WebAssembly and on Windows)
  double xx[PMAX * PMAX] = {};  // X'(I − P)X, Z'Z, Z'X of one cluster (row-major)
  double zz[QMAX * QMAX] = {};
  double zx[QMAX * PMAX] = {};
  double zi[QMAX * QMAX] = {};  // (Z'Z)⁺ and X's coefficients on Z, C = (Z'Z)⁺ Z'X (q × p)
  double cx[QMAX * PMAX] = {};
  int rank_z = 0;               // the rank and determinant of Z'Z: what is left within a cluster has n − rank_z
  double det_z = 0;             // dimensions
  double zl[QMAX * QMAX] = {};  // (Z'Z)⁺ in Z's coordinates, and S⁻¹, which takes coefficients from them
  double si[QMAX * QMAX] = {1, 0, 0, 1};
  double xy[PMAX] = {};  // sums over the clusters: Σ X'(I − P)y, Σ b, Σ b b', Σ y'(I − P)y
  double b[QMAX] = {};
  double bb[QMAX * QMAX] = {};
  double yy = 0;

  /// Set the design from X (n × p) and Z (n × q), row-major.
  void set_design(long rows, int p, int q, const double* x, const double* z);
  /// Set the design of an X that lies in the span of Z from Z (n × q, in its own coordinates), S (q × q, |S| = 1;
  /// null for the identity) and X's coefficients on the random effects' covariates Z S, X = Z S C (C: q × p,
  /// exact, so that nothing of X is left within a cluster).  g, when given, is the G the clusters are expected
  /// to have: a direction in which Z S is nearly singular is kept only if the random effects show in it.
  void set_spanned(long rows, int p, int q, const double* z, const double* c, const double* s = nullptr,
                   const double* g = nullptr);
  /// Set (Z'Z)⁺, C, the rank and the determinant from Z'Z and Z'X: set_design does, and a design whose products
  /// are written directly must.
  void set_coefficients(int p, int q);
  /// The rows of (I − P)X (n × p) from X and Z: what is left of X within a cluster of this design.  They depend
  /// only on the design, so a caller adding many clusters forms them once and passes them to add_effects.
  void within_rows(int p, int q, const double* x, const double* z, double* rx) const;
  /// Add one cluster's responses (n values); the design must be set.
  void add_cluster(int p, int q, const double* x, const double* z, const double* y);
  /// Add one cluster whose responses are Z S u + e (u: q values, e: n values), given X's rows within the cluster
  /// (rx, from within_rows): the sums come from e alone and u, so that large random effects u cost no digits.
  void add_effects(int p, int q, const double* rx, const double* z, const double* u, const double* e);
  /// add_effects with q = Q, and p = P unless P is 0, fixed at compile time, so that the loops over Z's and X's
  /// columns unroll where many clusters are added (a longitudinal simulation adds every subject this way).
  template <int P, int Q>
  void add_effects(int p, const double* rx, const double* z, const double* u, const double* e);
  /// Add one cluster whose design gives the sufficient statistics directly: X'(I − P)y, b and y'(I − P)y.
  void add_sums(int p, int q, const double* xy_, const double* b_, double yy_);
  void clear_responses();
};

template <int P, int Q>
void Pattern::add_effects(int p_, const double* rx, const double* z, const double* u, const double* e) {
  const int p = P > 0 ? P : p_;
  // e's coefficients on Z (and S⁻¹ times them on Z S) give b = u + (Z'Z)⁺ Z'e and e's residuals; X'(I − P)e
  // from X's residuals
  double ze[Q] = {}, coef_e[Q] = {};
  for (std::int64_t r = 0; r < n; ++r)
    for (int i = 0; i < Q; ++i) ze[i] += z[r * Q + i] * e[r];
  for (int i = 0; i < Q; ++i)
    for (int j = 0; j < Q; ++j) coef_e[i] += zl[i * Q + j] * ze[j];
  double sx[P > 0 ? P : PMAX] = {}, sb[Q], syy = 0;
  if (n > rank_z) {
    for (std::int64_t r = 0; r < n; ++r) {
      double v = e[r];
      for (int i = 0; i < Q; ++i) v -= z[r * Q + i] * coef_e[i];
      for (int a = 0; a < p; ++a) sx[a] += rx[r * p + a] * v;
      syy += v * v;
    }
  }
  for (int i = 0; i < Q; ++i) {
    sb[i] = u[i];
    for (int j = 0; j < Q; ++j) sb[i] += si[i * Q + j] * coef_e[j];
  }
  // add_sums
  count += 1;
  for (int a = 0; a < p; ++a) xy[a] += sx[a];
  for (int i = 0; i < Q; ++i) {
    b[i] += sb[i];
    for (int j = 0; j < Q; ++j) bb[i * Q + j] += sb[i] * sb[j];
  }
  yy += syy;
}

/// All clusters of one arm of a design whose X is the same row x on every observation of a cluster (a
/// cluster-randomized comparison), with a random intercept (q = 1), by size.  Nothing of X is left within such a
/// cluster, and its coefficients on Z are x itself, so that of a Pattern a size needs only its count and the sums
/// Σ b, Σ b² and Σ y'(I − P)y (Σ X'(I − P)y is 0): 40 bytes in place of 1,544, for sizes that vary so widely that
/// almost every cluster has one of its own (a million clusters took 1.9 GB as patterns).
struct Arm {
  struct Size {
    double n = 0;      // observations per cluster
    double count = 0;  // clusters of this size
    double b = 0;      // sums over them: Σ b (b: the cluster's mean), Σ b², Σ y'(I − P)y
    double bb = 0;
    double yy = 0;
  };
  double x[PMAX] = {};
  std::vector<Size> sizes;
};

struct Data {
  int p = 0;
  int q = 0;
  std::vector<Pattern> patterns;
  std::vector<Arm> arms;  // with q = 1 only; taken after the patterns
  /// Random effects (q = 2) whose covariates take two values only (a random slope with two visits): σ² and D are
  /// not told apart then, only the covariance of the two measurements is, and a D of rank 1 (θ2 = 0) reaches each
  /// one once, so the fit holds D there.
  bool two_rows = false;
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

/// −2 × REML log-likelihood at (σ², G) with β profiled out; NaN if V is not positive definite.  det_g as for
/// beta_covariance.
double criterion(const Data& d, double sigma2, const double* g,
                 double det_g = std::numeric_limits<double>::quiet_NaN());

/// Covariance of β̂ at (σ², G) (the generalized least squares covariance), into `cov` (p × p).  det_g: |G| from a
/// factor of G, when G is known as ΛΛ' (it is used if G is so nearly singular that |G| from its entries is not).
bool beta_covariance(const Data& d, double sigma2, const double* g, double* cov,
                     double det_g = std::numeric_limits<double>::quiet_NaN());

}  // namespace ssc::sim::lmm
