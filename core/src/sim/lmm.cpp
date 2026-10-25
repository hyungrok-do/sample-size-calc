#include "sim/lmm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace ssc::sim::lmm {

namespace {

/// Whether a 2 × 2 Z'Z with determinant det has full rank: one whose determinant is negligible next to the
/// product of its diagonal (one visit, or visits at the same time) has rank 1.
bool full_rank(const double* a, double det) { return det > 1e-12 * a[0] * a[3]; }
bool full_rank(const double* a) { return full_rank(a, a[0] * a[3] - a[1] * a[2]); }

/// Rank of Z'Z (q ≤ 2).
int rank(const double* a, int q) {
  if (q == 2 && full_rank(a)) return 2;
  return a[0] + (q == 2 ? a[3] : 0) > 0 ? 1 : 0;
}

/// Row r of (I − P)X: x_r minus z_r times X's coefficients on Z.
void within_row(const double* xr, const double* zr, const double* coef, int p, int q, double* res) {
  for (int b = 0; b < p; ++b) {
    res[b] = xr[b];
    for (int i = 0; i < q; ++i) res[b] -= zr[i] * coef[i * p + b];
  }
}

/// The Moore–Penrose inverse of a 2 × 2 A of rank 1: A / tr(A)².
void rank1_inverse(const double* a, double* ai) {
  const double tr = a[0] + a[3];
  for (int i = 0; i < 4; ++i) ai[i] = tr > 0 ? a[i] / (tr * tr) : 0;
}

/// M' A M and M A M' for 2 × 2 matrices.
void congruent(const double* m, const double* a, double* out, bool transposed) {
  double t[4];
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      t[i * 2 + j] = transposed ? m[i] * a[j] + m[2 + i] * a[2 + j] : m[i * 2] * a[j] + m[i * 2 + 1] * a[2 + j];
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      out[i * 2 + j] = transposed ? t[i * 2] * m[j] + t[i * 2 + 1] * m[2 + j]
                                  : t[i * 2] * m[j * 2] + t[i * 2 + 1] * m[j * 2 + 1];
}

}  // namespace

void Pattern::set_coefficients(int p, int q) {
  // the Moore–Penrose inverse (of a 2 × 2 Z'Z of rank 1, A⁺ = A / tr(A)²), and C = A⁺ Z'X with the adjugate
  // applied before dividing by the determinant, so that a column of X that is a column of Z, or a multiple of
  // one (the intercept, the time, the group), gets coefficients of exactly 0 and 1 and leaves exactly nothing
  // within the cluster: a rounding error there would be weighed against the between-cluster information, which
  // a large variance component makes tiny
  const double* a = zz;
  if (q == 1) {
    zi[0] = a[0] > 0 ? 1 / a[0] : 0;
    for (int b = 0; b < p; ++b) cx[b] = a[0] > 0 ? zx[b] / a[0] : 0;
  } else if (full_rank(a)) {
    const double det = a[0] * a[3] - a[1] * a[2];
    zi[0] = a[3] / det, zi[1] = -a[1] / det, zi[2] = -a[2] / det, zi[3] = a[0] / det;
    for (int b = 0; b < p; ++b) {
      cx[b] = (a[3] * zx[b] - a[1] * zx[p + b]) / det;
      cx[p + b] = (a[0] * zx[p + b] - a[2] * zx[b]) / det;
    }
  } else {
    rank1_inverse(a, zi);
    const double tr = a[0] + a[3];
    for (int b = 0; b < p; ++b) {
      cx[b] = tr > 0 ? (a[0] * zx[b] + a[1] * zx[p + b]) / (tr * tr) : 0;
      cx[p + b] = tr > 0 ? (a[2] * zx[b] + a[3] * zx[p + b]) / (tr * tr) : 0;
    }
  }
  rank_z = rank(a, q);
  det_z = q == 2 ? a[0] * a[3] - a[1] * a[2] : a[0];
  std::copy(zi, zi + q * q, zl);
  si[0] = 1, si[1] = 0, si[2] = 0, si[3] = 1;
}

void Pattern::set_spanned(long rows, int p, int q, const double* z, const double* c, const double* s,
                          const double* g) {
  n = rows;
  std::fill(std::begin(xx), std::end(xx), 0.0);
  std::fill(std::begin(zz), std::end(zz), 0.0);
  for (long r = 0; r < rows; ++r) {
    const double* zr = z + r * q;
    for (int i = 0; i < q; ++i)
      for (int j = 0; j < q; ++j) zz[i * q + j] += zr[i] * zr[j];
  }
  const double one[4] = {1, 0, 0, 1};
  if (!s || q == 1 || (s[0] == 1 && s[1] == 0 && s[2] == 0 && s[3] == 1)) s = one;
  // Z'X = Z'Z S C, then (Z'Z)⁺, the rank and |Z'Z| in Z's coordinates
  for (int i = 0; i < q; ++i)
    for (int b = 0; b < p; ++b) {
      double t = 0;
      for (int j = 0; j < q; ++j)
        for (int k = 0; k < q; ++k) t += zz[i * q + j] * s[j * q + k] * c[k * p + b];
      zx[i * p + b] = t;
    }
  set_coefficients(p, q);
  std::copy(c, c + q * p, cx);
  if (s == one) return;
  // in the random effects' coordinates: Z'Z → S'Z'Z S (whose determinant is the same), Z'X → S'Z'X and
  // coefficients S⁻¹ b; (Z'Z)⁻¹ → S⁻¹ (Z'Z)⁻¹ S⁻ᵀ from the inverse in Z's coordinates, where Z'Z is far from
  // singular
  const double a[4] = {zz[0], zz[1], zz[2], zz[3]};
  congruent(s, a, zz, true);
  for (int b = 0; b < p; ++b) {
    const double x0 = zx[b], x1 = zx[p + b];
    zx[b] = s[0] * x0 + s[2] * x1, zx[p + b] = s[1] * x0 + s[3] * x1;
  }
  si[0] = s[3], si[1] = -s[1], si[2] = -s[2], si[3] = s[0];
  bool full = rank_z == 2;
  if (full && !full_rank(zz, det_z) && g) {
    // a direction in which Z S is nearly singular (two visits close together, far from where the random
    // effects are given) is kept as long as the random effects show in it at more than 10⁻⁶ of the residual
    // variance: λ (v'Gv), with λ the small eigenvalue of S'Z'Z S and v its eigenvector.  Below, it is left within
    // the cluster (rank 1), which errs by as little, while keeping it would make the clusters' coefficients
    // swing along it far more than the random effects do, and their sums of squares lose their digits
    const double tr = zz[0] + zz[3];
    const double lambda = det_z / (0.5 * tr + std::sqrt(std::fmax(0.25 * tr * tr - det_z, 0.0)));
    const double v0 = zz[0] >= zz[3] ? -zz[1] : -zz[3], v1 = zz[0] >= zz[3] ? zz[0] : zz[2];
    const double vgv = (v0 * (g[0] * v0 + g[1] * v1) + v1 * (g[2] * v0 + g[3] * v1)) / (v0 * v0 + v1 * v1);
    full = lambda * vgv > 1e-6;
  }
  if (full) {
    congruent(si, zl, zi, false);
  } else {
    rank_z = std::min(rank_z, 1);
    rank1_inverse(zz, zi);
    congruent(s, zi, zl, false);
  }
}

void Pattern::set_design(long rows, int p, int q, const double* x, const double* z) {
  n = rows;
  std::fill(std::begin(xx), std::end(xx), 0.0);
  std::fill(std::begin(zz), std::end(zz), 0.0);
  std::fill(std::begin(zx), std::end(zx), 0.0);
  for (long r = 0; r < rows; ++r) {
    const double* xr = x + r * p;
    const double* zr = z + r * q;
    for (int i = 0; i < q; ++i) {
      for (int j = 0; j < q; ++j) zz[i * q + j] += zr[i] * zr[j];
      for (int b = 0; b < p; ++b) zx[i * p + b] += zr[i] * xr[b];
    }
  }
  set_coefficients(p, q);
  // a column of X that is a multiple s of a column of Z (the intercept, the group, the time, the group's time) gets
  // the coefficients s and 0 exactly: set_coefficients gives them where each product is rounded on its own, but a
  // compiler that fuses a multiply and an add (as on ARM) leaves a rounding error in the one that should vanish, which
  // is weighed within the cluster against the between-cluster information (with slopes 10⁵⁹ residual SDs apart, the
  // approximate power came out 1 for 0.97)
  if (rank_z == q) {
    for (int b = 0; b < p; ++b)
      for (int i = 0; i < q; ++i) {
        double s = 0;
        bool found = false, multiple = true;
        for (long r = 0; r < rows && multiple; ++r) {
          const double xv = x[r * p + b], zv = z[r * q + i];
          if (!found && zv != 0) s = xv / zv, found = true;
          multiple = xv == s * zv;
        }
        if (!multiple) continue;
        for (int j = 0; j < q; ++j) cx[j * p + b] = j == i ? s : 0;
        break;
      }
  }
  // X'(I − P)X from the rows' residuals on Z (none when the cluster has no more rows than Z has rank)
  if (rows <= rank_z) return;
  for (long r = 0; r < rows; ++r) {
    double res[PMAX];
    within_row(x + r * p, z + r * q, cx, p, q, res);
    for (int a = 0; a < p; ++a)
      for (int b = 0; b < p; ++b) xx[a * p + b] += res[a] * res[b];
  }
}

void Pattern::within_rows(int p, int q, const double* x, const double* z, double* rx) const {
  for (long r = 0; r < n; ++r) within_row(x + r * p, z + r * q, cx, p, q, rx + r * p);
}

void Pattern::add_cluster(int p, int q, const double* x, const double* z, const double* y) {
  std::vector<double> rx(static_cast<std::size_t>(n * p));
  within_rows(p, q, x, z, rx.data());
  const double none[QMAX] = {};
  add_effects(p, q, rx.data(), z, none, y);
}

void Pattern::add_effects(int p, int q, const double* rx, const double* z, const double* u, const double* e) {
  if (q == 1)
    add_effects<0, 1>(p, rx, z, u, e);
  else
    add_effects<0, 2>(p, rx, z, u, e);
}

void Pattern::add_sums(int p, int q, const double* xy_, const double* b_, double yy_) {
  count += 1;
  for (int a = 0; a < p; ++a) xy[a] += xy_[a];
  for (int i = 0; i < q; ++i) {
    b[i] += b_[i];
    for (int j = 0; j < q; ++j) bb[i * q + j] += b_[i] * b_[j];
  }
  yy += yy_;
}

void Pattern::clear_responses() {
  count = 0;
  std::fill(std::begin(xy), std::end(xy), 0.0);
  std::fill(std::begin(b), std::end(b), 0.0);
  std::fill(std::begin(bb), std::end(bb), 0.0);
  yy = 0;
}

double Data::observations() const {
  double s = 0;
  for (const auto& pt : patterns) s += pt.count * static_cast<double>(pt.n);
  for (const auto& arm : arms)
    for (const auto& sz : arm.sizes) s += sz.count * sz.n;
  return s;
}

double Data::clusters() const {
  double s = 0;
  for (const auto& pt : patterns) s += pt.count;
  for (const auto& arm : arms)
    for (const auto& sz : arm.sizes) s += sz.count;
  return s;
}

namespace {

constexpr double LOG_2PI = 1.8378770664093454836;

/// In-place lower Cholesky factor of a symmetric n × n matrix; false unless positive definite.
bool cholesky(double* a, int n) {
  for (int j = 0; j < n; ++j) {
    double s = a[j * n + j];
    for (int k = 0; k < j; ++k) s -= a[j * n + k] * a[j * n + k];
    if (!(s > 0)) return false;
    const double d = std::sqrt(s);
    a[j * n + j] = d;
    for (int i = j + 1; i < n; ++i) {
      double t = a[i * n + j];
      for (int k = 0; k < j; ++k) t -= a[i * n + k] * a[j * n + k];
      a[i * n + j] = t / d;
    }
  }
  return true;
}

/// Whether β's covariance from X'V*⁻¹X keeps at least three digits, from its Cholesky factor L: each pivot L_jj² is
/// what is left of the diagonal entry (Σ_k L_jk²) once the earlier columns are taken out, and the inverse loses
/// about as many digits as that share has zeros.  Nothing is left when the design cannot tell the fixed effects
/// apart (every subject of a group seen once), and next to nothing with random slopes 10⁸ times as variable as the
/// intercepts at time 0 when time is measured from a visit 1 time unit away, where the level all but follows from
/// the slope (an approximate power of 100% for 8.6% once; the longitudinal design now measures time from where the
/// two are uncorrelated).
bool keeps_digits(const double* l, int n) {
  for (int j = 0; j < n; ++j) {
    double big = 0, s = 0;
    for (int k = 0; k <= j; ++k) big = std::fmax(big, std::fabs(l[j * n + k]));
    for (int k = 0; k <= j; ++k) s += (l[j * n + k] / big) * (l[j * n + k] / big);
    const double pivot = l[j * n + j] / big;
    if (!(pivot * pivot > 1e-13 * s)) return false;
  }
  return true;
}

/// Solve L L' x = b in place, L from cholesky().
void chol_solve(const double* l, int n, double* b) {
  for (int i = 0; i < n; ++i) {
    double s = b[i];
    for (int k = 0; k < i; ++k) s -= l[i * n + k] * b[k];
    b[i] = s / l[i * n + i];
  }
  for (int i = n - 1; i >= 0; --i) {
    double s = b[i];
    for (int k = i + 1; k < n; ++k) s -= l[k * n + i] * b[k];
    b[i] = s / l[i * n + i];
  }
}

void chol_inverse(const double* l, int n, double* inv) {
  for (int j = 0; j < n; ++j) {
    double e[PMAX * PMAX] = {};
    e[j] = 1;
    chol_solve(l, n, e);
    for (int i = 0; i < n; ++i) inv[i * n + j] = e[i];
  }
}

/// The generalized least squares pieces at G = D/σ².
struct Gls {
  double m[PMAX * PMAX];  // X'V*⁻¹X, then its Cholesky factor
  double v[PMAX];         // X'V*⁻¹y
  double beta[PMAX];
  double logdet_k = 0;  // Σ log |K| = log |V*|
  double logdet_m = 0;  // log |X'V*⁻¹X|
  double rss = 0;       // (y − Xβ̂)'V*⁻¹(y − Xβ̂)
  double within = 0;    // y'V*⁻¹y within the clusters, Σ y'(I − P)y, which does not depend on G, and between them,
  double between = 0;   // Σ b'Wb
};

/// For a cluster with Z'Z = A (pt.zz, with A⁺, rank and determinant): W = K⁻¹ A, the between-cluster part of V*⁻¹
/// in the coordinates b, and |K| = |I + A G| (det_g = |G|); false unless V* is positive definite.  Computed so
/// that nothing cancels: for a 2 × 2 A of full rank, |K| = 1 + tr(AG) + |A||G| and W = (A⁻¹ + G)⁻¹, whose
/// determinant is |K|/|A|; for a 1 × 1 A, or one of rank 1, W = A / (1 + tr(AG)).
bool between(const Pattern& pt, const double* g, double det_g, int q, double* w, double& detk) {
  const double *a = pt.zz, *ai = pt.zi;
  const double tr_ag = q == 1 ? a[0] * g[0] : a[0] * g[0] + a[1] * g[2] + a[2] * g[1] + a[3] * g[3];
  if (q == 1 || pt.rank_z < 2) {
    detk = 1 + tr_ag;
    if (!(detk > 0)) return false;
    for (int i = 0; i < q * q; ++i) w[i] = a[i] / detk;
    return true;
  }
  const double det_a = pt.det_z;
  detk = 1 + tr_ag + det_a * det_g;
  const double e00 = ai[0] + g[0], e11 = ai[3] + g[3], e01 = 0.5 * (ai[1] + g[1] + ai[2] + g[2]);
  if (!(detk > 0) || !(e00 > 0)) return false;
  const double det_e = detk / det_a;
  w[0] = e11 / det_e, w[1] = w[2] = -e01 / det_e, w[3] = e00 / det_e;
  return true;
}

/// gls's sums over the patterns, with q = Q fixed at compile time: the loops over Z's columns unroll, which took a
/// tenth off a stepped wedge's time (q = 1, 12 fixed effects).
template <int Q>
bool pattern_sums(const Data& d, const double* g, double det_g, Gls& w, double& s) {
  const int p = d.p;
  for (const Pattern& pt : d.patterns) {
    if (pt.count == 0) continue;
    double wm[Q * Q], detk;
    if (!between(pt, g, det_g, Q, wm, detk)) return false;
    w.logdet_k += pt.count * std::log(detk);
    // W C (q × p), C = (Z'Z)⁺ Z'X the coefficients of X on Z
    double wc[Q * PMAX];
    for (int i = 0; i < Q; ++i)
      for (int b = 0; b < p; ++b) {
        double t = 0;
        for (int j = 0; j < Q; ++j) t += wm[i * Q + j] * pt.cx[j * p + b];
        wc[i * p + b] = t;
      }
    const double c = pt.count;
    for (int a = 0; a < p; ++a) {
      for (int b = a; b < p; ++b) {
        double t = 0;
        for (int i = 0; i < Q; ++i) t += pt.cx[i * p + a] * wc[i * p + b];
        w.m[a * p + b] += c * (pt.xx[a * p + b] + t);
      }
      double t = 0;
      for (int i = 0; i < Q; ++i) t += wc[i * p + a] * pt.b[i];
      w.v[a] += pt.xy[a] + t;
    }
    double tr = 0;
    for (int i = 0; i < Q; ++i)
      for (int j = 0; j < Q; ++j) tr += wm[i * Q + j] * pt.bb[j * Q + i];
    s += pt.yy + tr;
    w.within += pt.yy;
    w.between += tr;
  }
  return true;
}

/// The generalized least squares pieces at G; det_g, when given, is |G| from a factor of G, which takes the place of
/// |G| from its entries when that is mostly rounding (G nearly singular: a random intercept given at visits so far
/// from time 0 that it all but follows from the slope).
bool gls(const Data& d, const double* g, Gls& w, double det_g = std::numeric_limits<double>::quiet_NaN()) {
  const int p = d.p, q = std::min(d.q, QMAX);
  if (d.poll) d.poll(d.poll_arg);
  std::fill(w.m, w.m + p * p, 0.0);
  std::fill(w.v, w.v + p, 0.0);
  double s = 0;
  w.logdet_k = 0;
  w.within = w.between = 0;
  if (q == 2) {
    const double det_entries = g[0] * g[3] - g[1] * g[2];
    if (!std::isfinite(det_g) || std::fabs(det_entries) > 1e-8 * std::fabs(g[0] * g[3])) det_g = det_entries;
  }
  if (!(q == 1 ? pattern_sums<1>(d, g, det_g, w, s) : pattern_sums<2>(d, g, det_g, w, s))) return false;
  // the same sums for clusters by arm and size: Z'Z = n, C = x and nothing of X within the cluster.  The products are
  // formed as pattern_sums forms them (in statements of their own, into sums that start at 0), so that a compiler
  // that fuses a multiply and an add within a statement (clang's default, as on Apple silicon) fuses the same ones in
  // both, and the fits stay the same to the bit; one that fuses across statements (GCC with FMA) can round the two
  // apart (in 104 of 400 random studies, by up to 6 · 10⁻⁹ in the standard error and 9 · 10⁻⁵ in the degrees of freedom).
  // A caller with a million sizes is polled every 65,536 of them as well
  std::size_t polled = 0;
  for (const Arm& arm : d.arms) {
    const double* x = arm.x;
    for (const Arm::Size& sz : arm.sizes) {
      if (d.poll && ++polled % 65536 == 0) d.poll(d.poll_arg);
      if (sz.count == 0) continue;
      const double tr_ag = sz.n * g[0];
      const double detk = 1 + tr_ag;
      if (!(detk > 0)) return false;
      const double wm = sz.n / detk;
      w.logdet_k += sz.count * std::log(detk);
      double wc[PMAX];
      for (int b = 0; b < p; ++b) wc[b] = wm * x[b];
      for (int a = 0; a < p; ++a) {
        for (int b = a; b < p; ++b) {
          const double t = x[a] * wc[b];
          w.m[a * p + b] += sz.count * t;
        }
        double t = 0;
        t += wc[a] * sz.b;
        w.v[a] += t;
      }
      double tr = 0;
      tr += wm * sz.bb;
      s += sz.yy + tr;
      w.within += sz.yy;
      w.between += tr;
    }
  }
  for (int a = 0; a < p; ++a)
    for (int b = 0; b < a; ++b) w.m[a * p + b] = w.m[b * p + a];
  if (!cholesky(w.m, p)) return false;
  w.logdet_m = 0;
  for (int a = 0; a < p; ++a) w.logdet_m += 2 * std::log(w.m[a * p + a]);
  std::copy(w.v, w.v + p, w.beta);
  chol_solve(w.m, p, w.beta);
  double vb = 0;
  for (int a = 0; a < p; ++a) vb += w.v[a] * w.beta[a];
  w.rss = s - vb;  // callers that need a likelihood check that it is positive
  return true;
}

void lambda_to_g(int q, const double* theta, double* g) {
  if (q == 1) {
    g[0] = theta[0] * theta[0];
  } else {
    const double a = theta[0], b = theta[1], c = theta[2];
    g[0] = a * a;
    g[1] = g[2] = a * b;
    g[3] = b * b + c * c;
  }
}

/// REML criterion with σ² profiled out, as a function of Λ.
double profiled(const Data& d, const double* theta, Gls& w) {
  double g[QMAX * QMAX];
  lambda_to_g(d.q, theta, g);
  if (!gls(d, g, w) || !(w.rss > 0)) return std::numeric_limits<double>::infinity();
  const double r = d.observations() - d.p;
  return r * (1 + LOG_2PI + std::log(w.rss / r)) + w.logdet_k + w.logdet_m;
}

/// The REML criterion at (σ², G) in two parts: the one that depends on σ² alone, r (log 2π + log σ²) +
/// Σ y'(I − P)y / σ², which is nearly all of it when the clusters have many observations, and the rest.  False where
/// criterion() is NaN.
bool criterion_parts(const Data& d, double sigma2, const double* g, double det_g, double* parts) {
  Gls w;
  if (!(sigma2 > 0) || !gls(d, g, w, det_g) || !(w.rss > 0)) return false;
  double vb = 0;
  for (int a = 0; a < d.p; ++a) vb += w.v[a] * w.beta[a];
  parts[0] = (d.observations() - d.p) * (LOG_2PI + std::log(sigma2)) + w.within / sigma2;
  parts[1] = w.logdet_k + w.logdet_m + (w.between - vb) / sigma2;
  return true;
}

/// Brent's minimization of f on [a, b] (Brent 1973, "localmin").
template <class F>
double brent_min(F&& f, double a, double b, double tol, int& evals) {
  const double golden = 0.3819660112501051;
  double x = a + golden * (b - a), w = x, v = x;
  double fx = f(x), fw = fx, fv = fx;
  ++evals;
  double dist = 0, e = 0;
  for (int it = 0; it < 200; ++it) {
    const double m = 0.5 * (a + b);
    const double t1 = tol * std::fabs(x) + 1e-12, t2 = 2 * t1;
    if (std::fabs(x - m) <= t2 - 0.5 * (b - a)) break;
    bool use_golden = true;
    if (std::fabs(e) > t1) {
      double r = (x - w) * (fx - fv), qq = (x - v) * (fx - fw), pp = (x - v) * qq - (x - w) * r;
      qq = 2 * (qq - r);
      if (qq > 0) pp = -pp;
      qq = std::fabs(qq);
      const double etemp = e;
      e = dist;
      if (std::fabs(pp) < std::fabs(0.5 * qq * etemp) && pp > qq * (a - x) && pp < qq * (b - x)) {
        dist = pp / qq;
        const double u = x + dist;
        if (u - a < t2 || b - u < t2) dist = x < m ? t1 : -t1;
        use_golden = false;
      }
    }
    if (use_golden) {
      e = (x < m ? b : a) - x;
      dist = golden * e;
    }
    const double u = x + (std::fabs(dist) >= t1 ? dist : (dist > 0 ? t1 : -t1));
    const double fu = f(u);
    ++evals;
    if (fu <= fx) {
      (u < x ? b : a) = x;
      v = w, fv = fw, w = x, fw = fx, x = u, fx = fu;
    } else {
      (u < x ? a : b) = u;
      if (fu <= fw || w == x) {
        v = w, fv = fw, w = u, fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u, fv = fu;
      }
    }
  }
  return x;
}

/// Nelder–Mead minimization in `n` dimensions (n ≤ 3), restarted once from its result; the first simplex steps
/// `step` times scale(x, j) along each coordinate j.
template <class F, class Scale>
void nelder_mead(F&& f, int n, double* x, double step, int max_evals, int& evals, Scale&& scale) {
  for (int round = 0; round < 2; ++round) {
    double pts[4][3], val[4];
    for (int i = 0; i <= n; ++i) {
      for (int j = 0; j < n; ++j) pts[i][j] = x[j];
      if (i > 0) pts[i][i - 1] += step * scale(x, i - 1);
      val[i] = f(pts[i]);
      ++evals;
    }
    while (evals < max_evals) {
      int order[4] = {0, 1, 2, 3};
      std::sort(order, order + n + 1, [&](int a, int b) { return val[a] < val[b]; });
      const int best = order[0], worst = order[n], second = order[n - 1];
      double size = 0;
      for (int i = 0; i <= n; ++i)
        for (int j = 0; j < n; ++j) size = std::fmax(size, std::fabs(pts[i][j] - pts[best][j]));
      if (std::fabs(val[worst] - val[best]) <= 1e-10 * (1 + std::fabs(val[best])) && size < 1e-7) break;
      double centroid[3] = {};
      for (int i = 0; i <= n; ++i)
        if (i != worst)
          for (int j = 0; j < n; ++j) centroid[j] += pts[i][j] / n;
      auto along = [&](double t, double* out) {
        for (int j = 0; j < n; ++j) out[j] = centroid[j] + t * (pts[worst][j] - centroid[j]);
        ++evals;
        return f(out);
      };
      double xr[3], xe[3], xc[3];
      const double fr = along(-1, xr);
      if (fr < val[best]) {
        const double fe = along(-2, xe);
        if (fe < fr) {
          std::copy(xe, xe + n, pts[worst]);
          val[worst] = fe;
        } else {
          std::copy(xr, xr + n, pts[worst]);
          val[worst] = fr;
        }
      } else if (fr < val[second]) {
        std::copy(xr, xr + n, pts[worst]);
        val[worst] = fr;
      } else {
        const double fc = fr < val[worst] ? along(-0.5, xc) : along(0.5, xc);
        if (fc < std::fmin(fr, val[worst])) {
          std::copy(xc, xc + n, pts[worst]);
          val[worst] = fc;
        } else {  // shrink towards the best point
          for (int i = 0; i <= n; ++i) {
            if (i == best) continue;
            for (int j = 0; j < n; ++j) pts[i][j] = pts[best][j] + 0.5 * (pts[i][j] - pts[best][j]);
            val[i] = f(pts[i]);
            ++evals;
          }
        }
      }
    }
    int best = 0;
    for (int i = 1; i <= n; ++i)
      if (val[i] < val[best]) best = i;
    std::copy(pts[best], pts[best] + n, x);
    step *= 0.1;
  }
}

/// Each cluster's random-effect coefficients as far as its own measurements tell them, the mean of (Z'Z)⁺ over the
/// clusters (relative to σ²): the scale of D's entries below which the criterion no longer sees them.
void noise_scale(const Data& d, double* nu) {
  double n0 = 0, n1 = 0, count = 0;
  for (const Pattern& pt : d.patterns) {
    n0 += pt.count * pt.zi[0];
    if (d.q == 2) n1 += pt.count * pt.zi[3];
    count += pt.count;
  }
  for (const Arm& arm : d.arms)
    for (const Arm::Size& sz : arm.sizes) n0 += sz.count * (1 / sz.n), count += sz.count;  // as a Pattern's zi
  nu[0] = count > 0 ? n0 / count : 1;
  nu[1] = count > 0 ? n1 / count : 1;
}

/// Whether some cluster has more observations than its random effects have dimensions, which is what tells the
/// residual variance apart from theirs: with random slopes, a subject measured three times or more.
bool residual_identified(const Data& d) {
  for (const Pattern& pt : d.patterns)
    if (pt.count > 0 && pt.n > pt.rank_z) return true;
  for (const Arm& arm : d.arms)
    for (const Arm::Size& sz : arm.sizes)
      if (sz.count > 0 && sz.n > 1) return true;
  return false;
}

/// Whether D is held at rank 1 (θ2 = 0): with two visits (Data::two_rows), and in a study none of whose subjects
/// was measured three times (heavy dropout), where the same ridge of equal criteria runs to σ² = 0.
bool rank_one(const Data& d) { return d.q == 2 && (d.two_rows || !residual_identified(d)); }

/// The slope at G = 0 of the criterion with σ² profiled out, for a random intercept: Σ count n − Σ count n² C M⁻¹ C'
/// − Σ n² (b − Cβ̂)²/σ̂² over the clusters, with each record's residuals taken about their own mean first, (bb −
/// b²/count) + count (b/count − Cβ̂)², where a record of one cluster has no spread about its mean (bb is b² rounded:
/// with a cluster of 5.3 · 10⁸ subjects, that rounding, which a compiler that fuses a multiply and an add sees,
/// weighs 0.09 in the slope).  Expanded as exact_derivatives expands them, the squares lose every digit of it where
/// one cluster has 10⁸ subjects (−1.5 for a slope of +1.9).
///
/// Clusters by arm whose arms each have a mean of their own (as many fixed effects as arms) take every term in closed
/// form: at G = 0 the fit is each arm's mean weighted by size, C M⁻¹ C' = 1/N for an arm of N observations, so that
/// count n − count n² C M⁻¹ C' = count n (N − n)/N, and σ̂² comes from the same residuals.  Through M⁻¹ and β̂, they
/// lost digits where one cluster is nearly all of its arm: with clusters of 1, 5.3 · 10⁸ and 2 subjects in one arm and
/// 3 of 1 in the other, the slope came out +9.09, or −6.79 where the compiler fuses a multiply and an add, for +4.15,
/// and there θ stayed at 1.3 · 10⁻⁵, with 1.4 · 10⁸ degrees of freedom for 5.3 · 10⁸.
bool slope_at_zero(const Data& d, double& slope) {
  const int p = d.p;
  const double zero[1] = {0};
  Gls w;
  if (d.q != 1 || !gls(d, zero, w) || !(w.rss > 0)) return false;
  auto spread = [](double count, double b, double bb) {
    const double mean_part = b * (b / count);
    return count > 1 ? bb - mean_part : 0;
  };
  bool own = d.arms.size() == static_cast<std::size_t>(p);
  for (const Pattern& pt : d.patterns) own = own && pt.count == 0;
  if (own) {
    double rss = 0, lev = 0, sq = 0;
    for (const Arm& arm : d.arms) {
      double nn = 0, sum = 0;
      for (const Arm::Size& sz : arm.sizes) nn += sz.count * sz.n, sum += sz.n * sz.b;
      const double mean = sum / nn;
      for (const Arm::Size& sz : arm.sizes) {
        if (sz.count == 0) continue;
        const double dev = sz.b / sz.count - mean, ss = spread(sz.count, sz.b, sz.bb) + sz.count * dev * dev;
        rss += sz.yy + sz.n * ss;
        lev += sz.count * sz.n * ((nn - sz.n) / nn);
        sq += sz.n * sz.n * ss;
      }
    }
    slope = lev - sq / (rss / (d.observations() - p));
    return std::isfinite(slope);
  }
  const double s2 = w.rss / (d.observations() - p);
  double minv[PMAX * PMAX];
  chol_inverse(w.m, p, minv);
  slope = 0;
  auto add = [&](double count, double n, const double* c, double b, double bb) {
    double cb = 0, cmc = 0;
    for (int a = 0; a < p; ++a) {
      cb += c[a] * w.beta[a];
      for (int k = 0; k < p; ++k) cmc += c[a] * minv[a * p + k] * c[k];
    }
    const double mean = b / count - cb;
    slope += count * n - count * n * n * cmc - n * n * (spread(count, b, bb) + count * mean * mean) / s2;
  };
  for (const Pattern& pt : d.patterns)
    if (pt.count > 0) add(pt.count, pt.zz[0], pt.cx, pt.b[0], pt.bb[0]);
  for (const Arm& arm : d.arms)
    for (const Arm::Size& sz : arm.sizes)
      if (sz.count > 0) add(sz.count, sz.n, arm.x, sz.b, sz.bb);
  return std::isfinite(slope);
}

/// Exact derivatives for Satterthwaite's degrees of freedom at (σ², D) = s0 (1, G), with respect to x = (σ²/s0, δ),
/// D = s0 (G + Σ_a scale_a δ_a E_a) over D's entries a (E_a: the symmetric unit matrix of D00, D01 or D11; scale_a
/// only keeps the numbers in range): the gradient and Hessian of the REML criterion, and Var(c'β̂)/s0 and its gradient.
/// An s0 of 0 is σ² profiled out at G, rss/(N − p).
struct Exact {
  int k = 0;  // 1 + D's entries
  double grad[4] = {}, hess[16] = {}, var = 0, dvar[4] = {};
  double fg[3] = {}, fgg[9] = {}, rss_g[3] = {};  // in G = G0 + Σ_a scale_a δ_a E_a at fixed σ²: the criterion's
                                                   // gradient and Hessian, and the gradient of rss/s0
};

/// The criterion is r log σ² + Σ log |K| + log |M| + rss/σ² with M = X'V*⁻¹X, M, rss and K functions of G alone
/// through each cluster's W (gls), and dW = −W dG W, d log |K| = tr(W dG): every derivative in G is a sum over the
/// patterns of 2 × 2 products (U_a = W Ẽ_a W, and T_ab = W Ẽ_a W Ẽ_b W + W Ẽ_b W Ẽ_a W for the second), and those in
/// (σ², D) follow from G = D/σ².  Nothing is differenced, so a direction in which the criterion barely moves (the
/// variance along a singular D's missing direction, or of a level the visits are far from) keeps its digits.
bool exact_derivatives(const Data& d, double s0, const double* g, double det_g, const double* scale, const double* c,
                       Exact& e) {
  const int p = d.p, q = std::min(d.q, QMAX), nd = q == 2 ? 3 : 1, k = 1 + nd;
  Gls w;
  if (!gls(d, g, w, det_g) || !(w.rss > 0)) return false;
  if (s0 == 0) s0 = w.rss / (d.observations() - p);
  double minv[PMAX * PMAX];
  chol_inverse(w.m, p, minv);
  // the directions Ẽ_a = scale_a E_a, and G = Σ γ_a Ẽ_a
  double dir[3][QMAX * QMAX] = {}, gam[3];
  if (q == 1) {
    dir[0][0] = scale[0];
    gam[0] = g[0] / scale[0];
  } else {
    dir[0][0] = scale[0], dir[1][1] = dir[1][2] = scale[1], dir[2][3] = scale[2];
    gam[0] = g[0] / scale[0], gam[1] = g[1] / scale[1], gam[2] = g[3] / scale[2];
  }
  // sums over the clusters of the derivatives of Σ log |K| (L), M, X'V*⁻¹y (v) and y'V*⁻¹y (s), first (a) and second
  // (a, b ≤ a)
  double La[3] = {}, Lab[9] = {}, sa[3] = {}, sab[9] = {}, Ma[3][PMAX * PMAX] = {}, Mab[9][PMAX * PMAX] = {};
  double va[3][PMAX] = {}, vab[9][PMAX] = {};
  auto mul = [q](const double* a, const double* b, double* out) {
    for (int i = 0; i < q; ++i)
      for (int j = 0; j < q; ++j) {
        double t = 0;
        for (int l = 0; l < q; ++l) t += a[i * q + l] * b[l * q + j];
        out[i * q + j] = t;
      }
  };
  // M += sign · count · C'mC, v += sign · C'm Σb, s += sign · tr(m Σbb'), for a symmetric q × q m
  auto add = [&](double count, const double* cx, const double* sb, const double* sbb, const double* m, double sign,
                 double* M, double* v, double& s) {
    double mc[QMAX * PMAX];
    for (int i = 0; i < q; ++i)
      for (int b = 0; b < p; ++b) {
        double t = 0;
        for (int j = 0; j < q; ++j) t += m[i * q + j] * cx[j * p + b];
        mc[i * p + b] = t;
      }
    for (int a = 0; a < p; ++a) {
      for (int b = 0; b < p; ++b) {
        double t = 0;
        for (int i = 0; i < q; ++i) t += cx[i * p + a] * mc[i * p + b];
        M[a * p + b] += sign * count * t;
      }
      double t = 0;
      for (int i = 0; i < q; ++i) t += mc[i * p + a] * sb[i];
      v[a] += sign * t;
    }
    double tr = 0;
    for (int i = 0; i < q; ++i)
      for (int j = 0; j < q; ++j) tr += m[i * q + j] * sbb[j * q + i];
    s += sign * tr;
  };
  auto cluster = [&](double count, const double* wm, const double* cx, const double* sb, const double* sbb) {
    double u[3][QMAX * QMAX], t[QMAX * QMAX];
    for (int a = 0; a < nd; ++a) {
      mul(wm, dir[a], t);  // W Ẽ_a
      for (int i = 0; i < q; ++i) La[a] += count * t[i * q + i];
      mul(t, wm, u[a]);
      add(count, cx, sb, sbb, u[a], -1, Ma[a], va[a], sa[a]);
    }
    for (int a = 0; a < nd; ++a)
      for (int b = 0; b <= a; ++b) {
        double ue[QMAX * QMAX], t1[QMAX * QMAX], t2[QMAX * QMAX], tab[QMAX * QMAX];
        mul(u[a], dir[b], ue);
        for (int i = 0; i < q; ++i) Lab[a * 3 + b] -= count * ue[i * q + i];
        mul(ue, wm, t1);
        mul(u[b], dir[a], ue);
        mul(ue, wm, t2);
        for (int i = 0; i < q; ++i)
          for (int j = 0; j < q; ++j)
            tab[i * q + j] = 0.5 * (t1[i * q + j] + t2[i * q + j] + t1[j * q + i] + t2[j * q + i]);
        add(count, cx, sb, sbb, tab, 1, Mab[a * 3 + b], vab[a * 3 + b], sab[a * 3 + b]);
      }
  };
  if (q == 2) {
    const double det_entries = g[0] * g[3] - g[1] * g[2];
    if (!std::isfinite(det_g) || std::fabs(det_entries) > 1e-8 * std::fabs(g[0] * g[3])) det_g = det_entries;
  }
  for (const Pattern& pt : d.patterns) {
    if (pt.count == 0) continue;
    double wm[QMAX * QMAX], detk;
    if (!between(pt, g, det_g, q, wm, detk)) return false;
    cluster(pt.count, wm, pt.cx, pt.b, pt.bb);
  }
  std::size_t polled = 0;
  for (const Arm& arm : d.arms)
    for (const Arm::Size& sz : arm.sizes) {
      if (d.poll && ++polled % 65536 == 0) d.poll(d.poll_arg);
      if (sz.count == 0) continue;
      const double tr_ag = sz.n * g[0];  // as between() forms it, for compilers that fuse a multiply and an add
      const double detk = 1 + tr_ag;
      const double wm = sz.n / detk;
      cluster(sz.count, &wm, arm.x, &sz.b, &sz.bb);
    }
  // with β = M⁻¹v and u_a = v_a − M_a β: rss_a = s_a − 2β'v_a + β'M_aβ, rss_ab = s_ab − 2β'v_ab + β'M_abβ −
  // 2 u_a'M⁻¹u_b; (log |M|)_a = tr(M⁻¹M_a), (log |M|)_ab = tr(M⁻¹M_ab) − tr(M⁻¹M_a M⁻¹M_b); with κ = c'M⁻¹c,
  // κ_a = −(M⁻¹c)'M_a(M⁻¹c)
  const double* beta = w.beta;
  auto quad = [p](const double* M, const double* x, const double* y) {
    double t = 0;
    for (int a = 0; a < p; ++a)
      for (int b = 0; b < p; ++b) t += x[a] * M[a * p + b] * y[b];
    return t;
  };
  auto dot = [p](const double* x, const double* y) {
    double t = 0;
    for (int a = 0; a < p; ++a) t += x[a] * y[a];
    return t;
  };
  auto product = [p](const double* A, const double* B, double* out) {
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j) {
        double t = 0;
        for (int l = 0; l < p; ++l) t += A[i * p + l] * B[l * p + j];
        out[i * p + j] = t;
      }
  };
  double mc[PMAX] = {};
  if (c)
    for (int a = 0; a < p; ++a)
      for (int b = 0; b < p; ++b) mc[a] += minv[a * p + b] * c[b];
  const double kappa = c ? dot(c, mc) : 0, rho = w.rss / s0, r = d.observations() - p;
  double Fa[3], Fab[9], rhoa[3], ua[3][PMAX], mua[3][PMAX], mma[3][PMAX * PMAX];
  for (int a = 0; a < nd; ++a) {
    product(minv, Ma[a], mma[a]);
    for (int i = 0; i < p; ++i) {
      double t = va[a][i];
      for (int j = 0; j < p; ++j) t -= Ma[a][i * p + j] * beta[j];
      ua[a][i] = t;
    }
    for (int i = 0; i < p; ++i) mua[a][i] = 0;
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j) mua[a][i] += minv[i * p + j] * ua[a][j];
    double tr = 0;
    for (int i = 0; i < p; ++i) tr += mma[a][i * p + i];
    rhoa[a] = (sa[a] - 2 * dot(beta, va[a]) + quad(Ma[a], beta, beta)) / s0;
    Fa[a] = La[a] + tr + rhoa[a];
    e.dvar[1 + a] = -quad(Ma[a], mc, mc);
  }
  for (int a = 0; a < nd; ++a)
    for (int b = 0; b <= a; ++b) {
      double mmab[PMAX * PMAX], tr = 0, trmm = 0;
      product(minv, Mab[a * 3 + b], mmab);
      for (int i = 0; i < p; ++i) {
        tr += mmab[i * p + i];
        for (int j = 0; j < p; ++j) trmm += mma[a][i * p + j] * mma[b][j * p + i];
      }
      const double rhoab =
          (sab[a * 3 + b] - 2 * dot(beta, vab[a * 3 + b]) + quad(Mab[a * 3 + b], beta, beta) - 2 * dot(ua[a], mua[b])) /
          s0;
      Fab[a * 3 + b] = Fab[b * 3 + a] = Lab[a * 3 + b] + tr - trmm + rhoab;
    }
  // x = (σ²/s0, δ) through G = (G0 + Σ_a scale_a δ_a E_a)/(σ²/s0), at σ² = s0 and δ = 0
  double fg = 0, rg = 0, ggf = 0;
  for (int a = 0; a < nd; ++a) {
    fg += Fa[a] * gam[a];
    rg += rhoa[a] * gam[a];
    for (int b = 0; b < nd; ++b) ggf += Fab[a * 3 + b] * gam[a] * gam[b];
  }
  e.k = k;
  for (int a = 0; a < nd; ++a) {
    e.fg[a] = Fa[a], e.rss_g[a] = rhoa[a];
    for (int b = 0; b < nd; ++b) e.fgg[a * 3 + b] = Fab[a * 3 + b];
  }
  e.grad[0] = r - rho - fg;
  e.hess[0] = 2 * rho - r + 2 * rg + ggf + 2 * fg;
  for (int a = 0; a < nd; ++a) {
    double t = 0;
    for (int b = 0; b < nd; ++b) t += Fab[a * 3 + b] * gam[b];
    e.grad[1 + a] = Fa[a];
    e.hess[(1 + a) * k] = e.hess[1 + a] = -rhoa[a] - t - Fa[a];
    for (int b = 0; b < nd; ++b) e.hess[(1 + a) * k + 1 + b] = Fab[a * 3 + b];
  }
  e.var = kappa;
  double t = 0;
  for (int a = 0; a < nd; ++a) t += e.dvar[1 + a] * gam[a];
  e.dvar[0] = kappa - t;
  return true;
}

/// Whether Nelder–Mead's result is a minimum of the criterion with σ² profiled out, from its exact derivatives in Λ
/// (in G: the gradient F_a and the Hessian F_ab − rss_a rss_b / r), and if it is not (or `force`), Newton's method from
/// there.  Where Λ's entries are 10⁹ times the residual SD, the simplex's steps of 10⁻⁷ are below their rounding, and
/// in a direction the criterion barely sees (the level at a time origin far from the visits) it collapsed at a saddle
/// 0.05 to 0.9 above the optimum, or took all 4,000 evaluations, which put Satterthwaite's degrees of freedom, taken
/// off the optimum, up to 9% off (9.13 for exactly 10) or left them to the given ones.  A minimum is a point whose
/// scaled Hessian has no eigenvalue below −10⁻⁶ and whose Newton step promises to lower the criterion by at most 10⁻⁹
/// of it (ten times what Nelder–Mead stops at); it stays where it is, to the bit.  Each step is taken only if the
/// criterion falls (halved up to 30 times; a Hessian that is not positive definite is damped towards its diagonal).
/// `held`: θ2 stays 0 (D at rank 1).  `beat`: the criterion of a fit elsewhere, which this one gives up on once it is
/// not going to reach it (the criterion, less twice what a step from a positive definite Hessian promises, above it);
/// with none, it never gives up.
bool polish(const Data& d, double* theta, bool held, int& evals, bool force,
            double beat = std::numeric_limits<double>::infinity()) {
  const int nt = held ? 2 : 3;
  const double r = d.observations() - d.p;
  double nu[2];
  noise_scale(d, nu);
  Gls w;
  double f = profiled(d, theta, w);
  ++evals;
  bool settled = true;
  for (int it = 0; it < 50 && std::isfinite(f); ++it) {
    double g[QMAX * QMAX];
    lambda_to_g(2, theta, g);
    const double s00 = g[0] + nu[0], s11 = g[3] + nu[1], scale[3] = {s00, std::sqrt(s00 * s11), s11};
    const double det_g = theta[0] * theta[2] * (theta[0] * theta[2]);
    Exact e;
    if (!exact_derivatives(d, w.rss / r, g, det_g, scale, nullptr, e)) return settled;
    // δ_a = (G(θ)_a − G0_a)/scale_a: J_aj = ∂δ_a/∂θ_j, and the second derivatives of δ_a below
    const double J[3][3] = {{2 * theta[0] / s00, 0, 0},
                            {theta[1] / scale[1], theta[0] / scale[1], 0},
                            {0, 2 * theta[1] / s11, 2 * theta[2] / s11}};
    double hg[9], gt[3];
    for (int i = 0; i < nt; ++i) {
      gt[i] = 0;
      for (int a = 0; a < 3; ++a) gt[i] += e.fg[a] * J[a][i];
      for (int j = 0; j < nt; ++j) {
        double t = 0;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) t += J[a][i] * (e.fgg[a * 3 + b] - e.rss_g[a] * e.rss_g[b] / r) * J[b][j];
        hg[i * nt + j] = t;
      }
    }
    hg[0] += 2 * e.fg[0] / s00;
    hg[1] += e.fg[1] / scale[1], hg[nt] += e.fg[1] / scale[1];
    hg[nt + 1] += 2 * e.fg[2] / s11;
    if (nt == 3) hg[8] += 2 * e.fg[2] / s11;
    // the step, scaled by the Hessian's diagonal and damped until positive definite
    double step[3], sc[3];
    bool solved = false;
    for (int i = 0; i < nt; ++i) sc[i] = hg[i * nt + i] != 0 ? 1 / std::sqrt(std::fabs(hg[i * nt + i])) : 1;
    if (it == 0) {
      // a minimum: no eigenvalue of the scaled Hessian below −10⁻⁶ (a saddle, where the simplex collapsed, has one
      // of −10³)
      double h[9];
      for (int i = 0; i < nt; ++i)
        for (int j = 0; j < nt; ++j) h[i * nt + j] = hg[i * nt + j] * sc[i] * sc[j] + (i == j ? 1e-6 : 0);
      settled = cholesky(h, nt);
    }
    bool undamped = true;
    for (int i = 0; i < nt; ++i) undamped = undamped && hg[i * nt + i] > 0;
    for (double damp = 0; !solved && damp < 1e8; damp = damp == 0 ? 1e-8 : 10 * damp) {
      double h[9];
      for (int i = 0; i < nt; ++i) {
        for (int j = 0; j < nt; ++j) h[i * nt + j] = hg[i * nt + j] * sc[i] * sc[j];
        h[i * nt + i] = std::fabs(h[i * nt + i]) * (1 + damp) + damp;
        step[i] = -gt[i] * sc[i];
      }
      if (cholesky(h, nt)) {
        chol_solve(h, nt, step);
        solved = true;
      } else {
        undamped = false;
      }
    }
    if (!solved) return false;
    double moved = 0, promise = 0;
    for (int i = 0; i < nt; ++i) {
      step[i] *= sc[i];
      promise -= 0.5 * gt[i] * step[i];
      moved = std::fmax(moved, std::fabs(step[i]) / std::fmax(std::fabs(theta[i]), 1e-300));
    }
    if (it == 0) {
      settled = settled && !(promise > 1e-9 * (1 + std::fabs(f)));
      if (settled && !force) return true;
    }
    if (undamped && f - 2 * promise > beat) return settled;
    bool better = false;
    for (int half = 0; half < 30 && !better; ++half) {
      double next[3] = {theta[0], theta[1], nt == 3 ? theta[2] : 0};
      for (int i = 0; i < nt; ++i) next[i] = theta[i] + step[i];
      Gls wn;
      const double fn = profiled(d, next, wn);
      ++evals;
      if (fn < f) {
        std::copy(next, next + 3, theta);
        f = fn;
        w = wn;
        better = true;
      } else {
        for (int i = 0; i < nt; ++i) step[i] *= 0.5;
      }
    }
    if (!better || moved < 1e-12) return settled;
  }
  return settled;
}

/// Satterthwaite's degrees of freedom V² / (∇V' H⁻¹ ∇V) from the Hessian H of the criterion and the gradient of V in
/// k coordinates (scaled by H's diagonal first, and ∇V by V, whose square passes the largest double when the residual
/// variance is 10⁻¹⁵⁰ of the random effects'); NaN unless H is positive definite.
double exact_df(int k, const double* hess, const double* dvar, double var) {
  double h[16], y[4], gs[4], sc[4];
  for (int i = 0; i < k; ++i) {
    if (!(hess[i * k + i] > 0)) return std::numeric_limits<double>::quiet_NaN();
    sc[i] = 1 / std::sqrt(hess[i * k + i]);
  }
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < k; ++j) h[i * k + j] = hess[i * k + j] * sc[i] * sc[j];
    gs[i] = y[i] = dvar[i] / var * sc[i];
  }
  if (!cholesky(h, k)) return std::numeric_limits<double>::quiet_NaN();
  chol_solve(h, k, y);
  double quad = 0;
  for (int i = 0; i < k; ++i) quad += gs[i] * y[i];
  return 1 / quad;
}

/// The same data with the random effects' coordinates turned so that the unit vector (c0, c1) is the first of them:
/// with R = [[c0, −c1], [c1, c0]], b → R'b, Z'Z → R'Z'Z R, (Z'Z)⁺ → R'(Z'Z)⁺R, Z'X → R'Z'X, C → R'C and S⁻¹ → R'S⁻¹.
/// The criterion, β̂ and its covariance are the same at R'GR as before at G.
Data turned(const Data& d, double c0, double c1) {
  Data out = d;
  const double r[4] = {c0, -c1, c1, c0};
  // v[0] and v[stride], one vector's two coordinates
  auto turn = [c0, c1](double* v, int stride) {
    const double v0 = v[0], v1 = v[stride];
    v[0] = c0 * v0 + c1 * v1, v[stride] = c0 * v1 - c1 * v0;
  };
  for (Pattern& pt : out.patterns) {
    double a[4];
    std::copy(pt.zz, pt.zz + 4, a);
    congruent(r, a, pt.zz, true);
    std::copy(pt.zi, pt.zi + 4, a);
    congruent(r, a, pt.zi, true);
    std::copy(pt.bb, pt.bb + 4, a);
    congruent(r, a, pt.bb, true);
    turn(pt.b, 1);
    for (int b = 0; b < d.p; ++b) turn(pt.zx + b, d.p), turn(pt.cx + b, d.p);
    for (int j = 0; j < 2; ++j) turn(pt.si + j, 2);
  }
  return out;
}


}  // namespace

double criterion(const Data& d, double sigma2, const double* g, double det_g) {
  Gls w;
  if (!(sigma2 > 0) || !gls(d, g, w, det_g) || !(w.rss > 0)) return std::numeric_limits<double>::quiet_NaN();
  const double r = d.observations() - d.p;
  return r * (LOG_2PI + std::log(sigma2)) + w.logdet_k + w.logdet_m + w.rss / sigma2;
}

bool beta_covariance(const Data& d, double sigma2, const double* g, double* cov, double det_g) {
  Gls w;
  if (!gls(d, g, w, det_g) || (d.q == 2 && !keeps_digits(w.m, d.p))) return false;
  chol_inverse(w.m, d.p, cov);
  for (int i = 0; i < d.p * d.p; ++i) cov[i] *= sigma2;
  return true;
}

Fit fit(const Data& d, const double* theta0) {
  Fit out;
  const int p = d.p, q = d.q;
  if (p < 1 || p > PMAX || q < 1 || q > QMAX || (q != 1 && !d.arms.empty()) || d.observations() <= p) return out;
  Gls w;
  if (q == 1 && !residual_identified(d)) {
    // every cluster a single observation: the criterion is the same at every θ (V = σ²(1 + θ²) I), and θ = 0 gives
    // the residual variance all of it, with N − p degrees of freedom; the fits stopped anywhere, and the Hessian was
    // singular
    out.theta[0] = 0;
  } else if (q == 1) {
    auto f = [&](double t) {
      const double th[1] = {t};
      return profiled(d, th, w);
    };
    // bracket the minimum: grow the upper end while the criterion keeps falling
    double hi = theta0 && theta0[0] > 0 ? 3 * theta0[0] : 1, fhi = f(hi);
    out.evaluations = 1;
    for (double next = 3 * hi; next < 1e6; next *= 3) {
      const double fn = f(next);
      ++out.evaluations;
      if (!(fn < fhi)) {
        hi = next;
        break;
      }
      hi = next;
      fhi = fn;
    }
    double t = brent_min(f, 0, hi, 1e-9, out.evaluations);
    const double f0 = f(0);
    ++out.evaluations;
    if (f0 <= f(t)) t = 0;
    // a θ² between 10⁻¹⁰ (below which the fit is on the boundary anyway) and 10⁻⁶ of the clusters' own noise (the
    // mean of 1/n) is told from 0 by the criterion's slope at 0, not by its values, whose differences there are
    // rounding with clusters of 10⁶ subjects: fits ended at θ of 10⁻⁵ to 2 · 10⁻⁴ where the optimum was 0 (with 3
    // clusters a group of 10⁷ subjects on average, CV 5), where Satterthwaite's degrees of freedom, taken as at a
    // minimum, came out of finite differences not at all, and of exact derivatives as 10⁻⁸ to 0.0007, or 77,347 for
    // 3.6 · 10⁸
    double nu[2], slope;
    noise_scale(d, nu);
    if (t * t >= 1e-10 && t * t < 1e-6 * nu[0] && slope_at_zero(d, slope) && !(slope < 0)) t = 0;
    // with clusters of more than 100 observations on average the criterion's rounding, which grows with the
    // observations, hides its minimum from Brent's method (θ² up to 6·10⁻⁴ off with 10⁷ subjects a cluster, which
    // took Satterthwaite's degrees of freedom 0.2% off): Newton's method in θ², from the exact slope and curvature of
    // the criterion with σ² profiled out (F_a and F_aa − rss_a²/r), takes it there, a step at a time while each is
    // under half of θ², until one under 10⁻⁶ of it (after which the next would be below rounding)
    const double r = d.observations() - p;
    for (int it = 0; it < 20 && t > 0 && d.observations() > 100 * d.clusters(); ++it) {
      const double g[1] = {t * t}, scale[3] = {t * t + nu[0], 0, 0};
      Exact e;
      ++out.evaluations;
      if (!exact_derivatives(d, 0, g, std::numeric_limits<double>::quiet_NaN(), scale, nullptr, e)) break;
      const double curvature = e.fgg[0] - e.rss_g[0] * e.rss_g[0] / r;
      const double step = -e.fg[0] / curvature * scale[0];
      if (!(curvature > 0) || !(std::fabs(step) < 0.5 * g[0])) break;
      t = std::sqrt(g[0] + step);
      if (std::fabs(step) < 1e-6 * g[0]) break;
    }
    out.theta[0] = t;
  } else {
    double x[3] = {1, 0, 1};
    if (theta0) std::copy(theta0, theta0 + 3, x);
    // each coordinate steps by its own size (at least 1), except that Λ's off-diagonal θ1 sets the correlation
    // θ1/√(θ1² + θ2²): where θ2 is over 2²⁰ times that step (slopes that spread the subjects over a million residual
    // SDs, from a correlation near 0), it steps by θ2, as steps of its own size would leave the correlation where it
    // started (with an intercept SD 10⁴⁷ times below the slopes', the fits stopped where the Hessian was not positive
    // definite).  Only there: the fits that converge with steps of their own size keep them
    auto scale = [](const double* th, int j) {
      const double own = std::fmax(1.0, std::fabs(th[j]));
      return j == 1 && std::fabs(th[2]) > 1048576 * own ? std::fabs(th[2]) : own;
    };
    // every point of a ridge from a singular D to σ² = 0 has the same V and criterion when no subject was measured
    // three times: the fits stopped anywhere on it, and the Hessian in (σ², D) was singular, which left
    // Satterthwaite's degrees of freedom to rounding or to the given ones (in 69 of 200 studies with visits at 0 and 1
    // and 60 subjects a group, and 21 of 200 with five visits and all but 12% of the subjects seen once); held at
    // rank 1, they are those of the covariance of the two measurements
    const bool held = rank_one(d);
    if (held) x[2] = 0;
    auto minimize = [&](double* th, int& evals) {
      if (held) {
        nelder_mead(
            [&](const double* y) {
              const double t[3] = {y[0], y[1], 0};
              return profiled(d, t, w);
            },
            2, th, 0.5, 4000, evals, scale);
      } else {
        nelder_mead([&](const double* y) { return profiled(d, y, w); }, 3, th, 0.5, 4000, evals, scale);
      }
    };
    minimize(x, out.evaluations);
    if (!polish(d, x, held, out.evaluations, out.evaluations >= 4000)) {
      // short of a minimum: a fresh simplex from there, which the criterion's flat valleys need (the simplex
      // collapsed 0.05 to 0.9 above the optimum in 10 of 200 studies with Λ's entries 10⁹ times the residual SD),
      // kept if it ends lower
      double y[3] = {x[0], x[1], x[2]};
      int evals = 0;
      minimize(y, evals);
      polish(d, y, held, evals, true);
      out.evaluations += evals;
      Gls wx, wy;
      if (profiled(d, y, wy) < profiled(d, x, wx)) std::copy(y, y + 3, x);
    }
    // near the boundary (the intercept's variance at the time where the subjects' lines cross with the least spread,
    // |G|/G22, below a tenth of its variance at time 0 or below 10⁻⁶ σ²) the fit can stop short of a boundary optimum
    // where the criterion is flat across its valley (221 degrees of freedom for 148 in a study whose fit ended
    // 3·10⁻¹⁰ of the criterion above it), or well above it (0.04 above, with intercepts and slopes correlated 0.9992,
    // 22.78 for 22.86): the best fit with D at rank 1, from the one that keeps G's diagonal, is taken where it is no
    // worse; so is a fit that ended on the boundary (θ2 = 0), which was left as it was
    const double cross = x[0] * x[2] * (x[0] * x[2]), g22 = x[1] * x[1] + x[2] * x[2];
    if (!held && x[0] != 0 && cross < std::fmax(0.1 * x[0] * x[0], 1e-6) * g22) {
      double y[3] = {x[0], std::copysign(std::sqrt(g22), x[1]), 0};
      int evals = 0;
      Gls wx, wy;
      const double fx = profiled(d, x, wx);
      polish(d, y, true, evals, true, fx);
      out.evaluations += evals + 1;
      if (!(profiled(d, y, wy) > fx)) std::copy(y, y + 3, x);
    }
    std::copy(x, x + 3, out.theta);
  }
  // a fit whose X'V⁻¹X has lost its digits has no estimate of β: with every subject of a group seen once, nothing
  // tells that group's slope apart, and its estimate came out exactly 0 or 1 with a standard error of 10⁷; with two
  // visits and one subject of a group measured twice and two of the other, the criterion has no minimum (it falls
  // without bound as the covariance of the two measurements turns singular), and the fits ended with σ² 10⁻²² of
  // the random effects' variances (the longitudinal design fails such studies before the fit).  Nor has one whose
  // estimate the arithmetic cannot reach: with slopes 5 · 10⁷ times the residual SD and a single subject beyond one a
  // group to tell their variance, the sums of squares keep too few digits for the fit to see the REML optimum, and
  // the standard errors of 11 of 12 such studies were 2% to 18% off it (the twelfth had 14 degrees of freedom for 8)
  out.criterion = profiled(d, out.theta, w);
  if (!std::isfinite(out.criterion) || !keeps_digits(w.m, p)) return out;
  lambda_to_g(q, out.theta, out.g);
  const double r = d.observations() - p;
  out.sigma2 = w.rss / r;
  std::copy(w.beta, w.beta + p, out.beta);
  chol_inverse(w.m, p, out.cov);
  for (int i = 0; i < p * p; ++i) out.cov[i] *= out.sigma2;
  if (q == 1) {
    out.boundary = out.g[0] < 1e-10;
  } else {
    const double g00 = out.g[0], g11 = out.g[3], scale = std::fmax(1.0, std::fmax(g00, g11));
    out.boundary = g00 < 1e-10 * scale || g11 < 1e-10 * scale ||
                   std::fabs(out.g[1]) > (1 - 1e-6) * std::sqrt(g00 * g11);
  }
  out.ok = true;
  return out;
}

Wald wald(const Data& d, const Fit& f, const double* c, int df_method, double fixed_df) {
  Wald out;
  const int p = d.p, q = d.q;
  double var = 0;
  for (int a = 0; a < p; ++a) {
    out.estimate += c[a] * f.beta[a];
    for (int b = 0; b < p; ++b) var += c[a] * f.cov[a * p + b] * c[b];
  }
  out.se = std::sqrt(std::fmax(var, 0.0));
  out.t = out.se > 0 ? out.estimate / out.se : 0;
  if (df_method == 2) {
    out.df = std::numeric_limits<double>::infinity();
    return out;
  }
  if (df_method == 1) {
    out.df = fixed_df;
    return out;
  }
  // Satterthwaite: df = 2 V² / (∇V' Cov(φ̂) ∇V), φ = (σ², D), Cov(φ̂) = 2 H⁻¹ (H: Hessian of the criterion).
  const double s2 = f.sigma2;
  auto good = [](double df) { return std::isfinite(df) && df > 0; };
  double nu[2];
  noise_scale(d, nu);
  // from exact derivatives (exact_derivatives), in units of the fit's σ², at the fitted D in full: G = g in the data
  // dd, whose noise scales are nn (f.g in d, or both turned near the boundary, below)
  auto interior = [&](const Data& dd, const double* g, const double* nn, const double* det_g) {
    const double s00 = std::fabs(g[0]) + nn[0], s11 = q == 2 ? std::fabs(g[3]) + nn[1] : 0;
    const double scale[3] = {s00, std::sqrt(s00 * s11), s11};
    Exact e;
    if (!exact_derivatives(dd, s2, g, det_g ? *det_g : std::numeric_limits<double>::quiet_NaN(), scale, c, e))
      return std::numeric_limits<double>::quiet_NaN();
    return exact_df(e.k, e.hess, e.dvar, e.var);
  };
  double df = std::numeric_limits<double>::quiet_NaN();
  if (q == 1) {
    // a random intercept: from finite differences of the criterion and of V, and from exact derivatives where they
    // fail (a Hessian that is not positive definite); one estimated at zero is held there (φ = σ² only, which gives
    // N − p).  Clusters of more than 100 observations on average (n̄) take exact derivatives first, and where those
    // fail the criterion in two parts with the intercept's step floor from n̄ (below); smaller ones keep the steps and
    // differences they had, to the bit
    const double nbar = d.observations() / d.clusters();
    const bool large = nbar > 100;
    auto satterthwaite = [&](int k, const double* phi) {
      auto unpack = [&](const double* x, double& sigma2, double* g) {
        sigma2 = x[0];
        g[0] = k > 1 ? x[1] / sigma2 : 0;
      };
      auto crit = [&](const double* x, double* parts) {
        double sigma2, g[1];
        unpack(x, sigma2, g);
        if (!large)
          parts[0] = criterion(d, sigma2, g), parts[1] = 0;
        else if (!criterion_parts(d, sigma2, g, std::numeric_limits<double>::quiet_NaN(), parts))
          parts[0] = parts[1] = std::numeric_limits<double>::quiet_NaN();
      };
      auto cvar = [&](const double* x) {
        double sigma2, g[1], cov[PMAX * PMAX];
        unpack(x, sigma2, g);
        if (!beta_covariance(d, sigma2, g, cov)) return std::numeric_limits<double>::quiet_NaN();
        double v = 0;
        for (int a = 0; a < p; ++a)
          for (int b = 0; b < p; ++b) v += c[a] * cov[a * p + b] * c[b];
        return v;
      };
      // a variance's step is relative to itself, with a floor of 10⁻² σ² near 0.  A random intercept shows in
      // V = σ² I + σb² J from σb² of about σ²/n, so with clusters of more than 100 observations on average its floor
      // is σ²/n̄: with 10⁷ subjects a cluster and an ICC of 10⁻⁶, steps of 10⁻⁶ σ² took σb² to 0 and back, and the
      // degrees of freedom came out 3 to 40 times too large where the Hessian was positive definite at all
      double h[2];
      const double least = large ? s2 / nbar : 1e-2 * s2;
      for (int j = 0; j < k; ++j) h[j] = j == 0 ? 1e-4 * s2 : 1e-4 * (2 * std::fabs(phi[j]) + least);
      // the derivatives are taken in units of the steps and of V, each a power of 2 (hs = h 2^e, V 2^m), which keeps
      // every product and quotient below exact to the last bit and away from under- and overflow
      double hess[4], grad[2], x[2], hs[2];
      for (int j = 0; j < k; ++j) hs[j] = std::ldexp(h[j], -std::ilogb(h[j]));
      const int m = var > 0 ? -std::ilogb(var) : 0;
      // with large clusters the criterion's two parts (criterion_parts) are differenced each on its own: in the
      // direction of σb², and across, the part that depends on σ² alone is the same at every point and drops out
      // exactly, where in the whole criterion its rounding swamped the differences (with 10⁷ subjects a cluster,
      // degrees of freedom 10 to 15% off and Hessians that were not positive definite; 2 · 10⁻⁵ off with 1,000)
      double f0[2], fp[2], fm[2], fpp[2], fpm[2], fmp[2], fmm[2];
      crit(phi, f0);
      bool ok = std::isfinite(f0[0] + f0[1]) && var > 0;
      for (int a = 0; a < k && ok; ++a) {
        std::copy(phi, phi + k, x);
        x[a] = phi[a] + h[a];
        crit(x, fp);
        const double vp = cvar(x);
        x[a] = phi[a] - h[a];
        crit(x, fm);
        const double vm = cvar(x);
        hess[a * k + a] = ((fp[0] - 2 * f0[0] + fm[0]) + (fp[1] - 2 * f0[1] + fm[1])) / (hs[a] * hs[a]);
        grad[a] = std::ldexp(vp - vm, m) / (2 * hs[a]);
        for (int b = 0; b < a; ++b) {
          std::copy(phi, phi + k, x);
          x[a] = phi[a] + h[a], x[b] = phi[b] + h[b];
          crit(x, fpp);
          x[b] = phi[b] - h[b];
          crit(x, fpm);
          x[a] = phi[a] - h[a];
          crit(x, fmm);
          x[b] = phi[b] + h[b];
          crit(x, fmp);
          const double cross = large ? ((fpp[0] - fpm[0]) - (fmp[0] - fmm[0])) + ((fpp[1] - fpm[1]) - (fmp[1] - fmm[1]))
                                     : fpp[0] - fpm[0] - fmp[0] + fmm[0];
          hess[a * k + b] = hess[b * k + a] = cross / (4 * hs[a] * hs[b]);
        }
        ok = std::isfinite(hess[a * k + a]) && std::isfinite(grad[a]);
      }
      if (!ok || !cholesky(hess, k)) return std::numeric_limits<double>::quiet_NaN();
      double y[2];
      std::copy(grad, grad + k, y);
      chol_solve(hess, k, y);  // H⁻¹ ∇V
      double quad = 0;
      for (int a = 0; a < k; ++a) quad += grad[a] * y[a];
      const double v = std::ldexp(var, m);
      return v * v / quad;  // 2 V² / (∇V' · 2 H⁻¹ · ∇V)
    };
    const double phi[2] = {s2, s2 * f.g[0]};
    if (large && !f.boundary) {
      // from exact derivatives first: the two parts' differences still left up to 6% in the degrees of freedom of
      // studies with 3 clusters a group of 10⁷ subjects on average (CV 5)
      df = interior(d, f.g, nu, nullptr);
      if (!good(df)) df = satterthwaite(2, phi);
    } else {
      df = satterthwaite(f.boundary ? 1 : 2, phi);
      if (!good(df) && !f.boundary) df = interior(d, f.g, nu, nullptr);
    }
  } else {
    // with random slopes, from exact derivatives: finite differences in D left up to 1% in the degrees of freedom of
    // fits whose intercepts and slopes are correlated 0.9 (0.18% where theory gives exactly N − 2), and near a
    // singular D, where the steps reach matrices that are not covariances, they left the Hessian not positive definite
    // or the degrees of freedom to rounding (55 for a 60-digit 16.8 with visits before time 0; 9.05 to 10.7 for exactly
    // 10 with calendar years; 56 of 100 studies on N − 2 where the subjects' lines cross 10¹⁴ time units before the
    // visits).  A fit on the boundary, with D = w w' singular (the subjects' lines all cross at one time, or their
    // slopes do not vary), is taken on that boundary, in (σ², w): the variance in the missing direction is held at 0,
    // as a random intercept estimated at zero is, and as lmerTest's derivatives in the factor Λ hold it.  Whether a
    // fit near the boundary is on it is told by the criterion's slope off the boundary, at w from Λ's first column:
    // the fits stop short of θ2 = 0 where the criterion is flat in it.  Only where neither way gives a positive
    // definite Hessian at an optimum are the given degrees of freedom used.
    //
    // Near the boundary when the intercept's variance at the time where the subjects' lines cross with the least
    // spread (|D|/D22, with |G| from Λ: from G's entries it is the small difference of large numbers) is below a tenth
    // of its variance at time 0 (a correlation beyond 0.95); flat when that variance, or the slope's (in units of the
    // squared root mean square τ² of the slope's covariate), is below 10⁻⁶ σ²
    double tau = 1, tt = 0, nn = 0;
    for (const Pattern& pt : d.patterns) tt += pt.count * pt.zz[3], nn += pt.count * pt.zz[0];
    if (tt > 0 && nn > 0) tau = std::sqrt(tt / nn);
    const double det_g = f.theta[0] * f.theta[2] * (f.theta[0] * f.theta[2]);
    const bool flat_slope = !(f.g[3] * tau * tau > 1e-6), flat = flat_slope || !(det_g > 1e-6 * f.g[3]);
    if (!flat && !(det_g < 0.1 * (f.g[0] * f.g[3]))) {
      df = interior(d, f.g, nu, &det_g);
    } else {
      // w from Λ's first column (G = w w' when θ2 = 0), or the slope's alone when the intercept does not vary
      // (σ² alone when neither does); held at rank 1 by the fit when no subject was measured three times, where it is
      // no boundary
      const bool none = flat_slope && !(f.g[0] > 1e-6);
      const double w0 = none || f.theta[0] == 0 ? 0 : f.theta[0];
      const double w1 = none ? 0 : f.theta[0] != 0 ? f.theta[1] : std::hypot(f.theta[1], f.theta[2]);
      // the derivatives are taken with the random effects' coordinates turned so that w is the first of them, where
      // G = diag(|w|², 0).  In their own coordinates, with w long next to σ, W = (A⁻¹ + G)⁻¹ has an eigenvalue near
      // 1/|w|² beside ones near 1, which its entries there keep to a share of their digits only, and each of D's
      // entries reaches into the direction off w, whose information swamps what tells w's length apart: with two
      // visits at 809.796 and 819.952 (σe 0.197, Λ some 7,000 times σe), where theory gives exactly N − 2 = 18 with
      // complete data, 15 of 200 studies were more than 0.1% off (up to 4.2%), and with dropout 0.5 one of 200 fell
      // back to the given degrees of freedom and two had 13.55 and 19.27 where 60-digit fits give 18.0004 and 18.0001
      const double len = std::hypot(w0, w1), c0 = len > 0 ? w0 / len : 1, c1 = len > 0 ? w1 / len : 0;
      const Data dw = turned(d, c0, c1);
      double nw[2];
      noise_scale(dw, nw);
      const double gw[4] = {len * len, 0, 0, 0};
      // D in full at the fit is taken in the same turned coordinates, where its derivatives lost digits the same way at
      // a fit near the boundary whose optimum is inside (with visits at 30 to 40, an intercept SD of 5 and σe 0.0005,
      // three of 400 studies had 35.5, 21.6 and 19.2 degrees of freedom where 60-digit fits give 8.000; with three
      // visits near 875 a study fell back to N − 2 for a 60-digit 13.9997): with w from Λ's first column, G = w w' +
      // θ2² e e' (e = (0, 1)) turns to diag(|w|², 0) + θ2² u u' with u = (c1, c0), formed from Λ without cancelling;
      // G = w w' when w is the slope's alone, and stays as it is when w = 0
      const double t2 = w0 != 0 ? f.theta[2] * f.theta[2] : 0;
      double gt[4] = {len * len + t2 * c1 * c1, t2 * c0 * c1, t2 * c0 * c1, t2 * c0 * c0};
      if (!(len > 0)) std::copy(f.g, f.g + 4, gt);
      const double s00 = gw[0] + nw[0], s11 = nw[1], scale[3] = {s00, std::sqrt(s00 * s11), s11};
      Exact e;
      const bool held = exact_derivatives(dw, s2, gw, 0, scale, c, e);
      const bool rank1 = rank_one(d);
      bool boundary = flat || rank1 || none;
      // the criterion's derivative along v v' (v ⟂ w, the second coordinate): rising off the boundary, the fit is on it
      if (held && !rank1 && !none) boundary = !(e.grad[3] < 0);
      auto on_boundary = [&]() {
        if (!held) return std::numeric_limits<double>::quiet_NaN();
        if (none) return exact_df(1, e.hess, e.dvar, e.var);
        // in x = (σ²/s0, η), w = (|w|, 0) + (√scale_00, √scale_11) η: J_aj = ∂δ_a/∂η_j and the Hessians of δ_a in η
        // (the criterion's gradient in δ does not vanish on the boundary)
        const double r0 = std::sqrt(s00), r1 = std::sqrt(s11);
        const double J[3][2] = {{r0 * 2 * len / scale[0], 0}, {0, r1 * len / scale[1]}, {0, 0}};
        const double D2[3][2][2] = {{{r0 * r0 * 2 / scale[0], 0}, {0, 0}},
                                    {{0, r0 * r1 / scale[1]}, {r0 * r1 / scale[1], 0}},
                                    {{0, 0}, {0, r1 * r1 * 2 / scale[2]}}};
        double hx[9], gx[3], fx[3];
        hx[0] = e.hess[0];
        gx[0] = e.dvar[0];
        fx[0] = e.grad[0];
        for (int j = 0; j < 2; ++j) {
          double hj = 0, gj = 0, fj = 0;
          for (int a = 0; a < 3; ++a)
            hj += e.hess[1 + a] * J[a][j], gj += e.dvar[1 + a] * J[a][j], fj += e.grad[1 + a] * J[a][j];
          hx[1 + j] = hx[(1 + j) * 3] = hj;
          gx[1 + j] = gj;
          fx[1 + j] = fj;
          for (int i = 0; i < 2; ++i) {
            double t = 0;
            for (int a = 0; a < 3; ++a) {
              for (int b = 0; b < 3; ++b) t += J[a][i] * e.hess[(1 + a) * 4 + 1 + b] * J[b][j];
              t += e.grad[1 + a] * D2[a][i][j];
            }
            hx[(1 + i) * 3 + 1 + j] = t;
          }
        }
        // Satterthwaite's degrees of freedom belong to an optimum: where a Newton step along the boundary would take
        // the criterion lower by more than 10⁻⁸ of it (½ ∇F'H⁻¹∇F, exact_df's quadratic form with the criterion's
        // gradient), w w' is not the boundary optimum, and the boundary gives none.  At the boundary optima of the
        // designs above and of others with heavy dropout such a step promises 1.3 · 10⁻⁹ or less.  It promised 0.005
        // from a fit near the boundary whose optimum is inside (Λ's θ2 1.07 beside θ1 3.67), which took the boundary's
        // 48.9 degrees of freedom for the interior's 28.02 (60-digit fits); and 5 · 10⁻⁵ to 0.5 in 11 of 200 studies
        // with a residual SD 10⁻⁶ of the slopes' (times 0, 0.59, 0.82, 1.25 and 2.07, dropout 0.899, 5 subjects a
        // group), whose fits stop where the criterion in doubles cannot see the boundary optimum, and which had 5.3 to
        // 9.4 degrees of freedom there for its 8
        if (0.5 / exact_df(3, hx, fx, 1) > 1e-8 * (1 + std::fabs(f.criterion)))
          return std::numeric_limits<double>::quiet_NaN();
        return exact_df(3, hx, gx, e.var);
      };
      // the other way where the first fails, except that D held at rank 1 by the fit stays there
      if (boundary) {
        df = on_boundary();
        if (!good(df) && !rank1) df = interior(dw, gt, nw, &det_g);
      } else {
        df = interior(dw, gt, nw, &det_g);
        if (!good(df)) df = on_boundary();
      }
    }
  }
  if (good(df)) {
    out.df = std::fmin(df, d.observations() - p);
    return out;
  }
  out.df = fixed_df;
  out.fallback = true;
  return out;
}

}  // namespace ssc::sim::lmm
