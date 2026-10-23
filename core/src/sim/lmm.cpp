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

/// The generalized least squares pieces at G; det_g, when given, is |G| from a factor of G, which takes the place of
/// |G| from its entries when that is mostly rounding (G nearly singular: a random intercept given at visits so far
/// from time 0 that it all but follows from the slope).
bool gls(const Data& d, const double* g, Gls& w, double det_g = std::numeric_limits<double>::quiet_NaN()) {
  const int p = d.p, q = std::min(d.q, QMAX);  // bounded, so that the compiler unrolls the loops over Z's columns
  std::fill(w.m, w.m + p * p, 0.0);
  std::fill(w.v, w.v + p, 0.0);
  double s = 0;
  w.logdet_k = 0;
  if (q == 2) {
    const double det_entries = g[0] * g[3] - g[1] * g[2];
    if (!std::isfinite(det_g) || std::fabs(det_entries) > 1e-8 * std::fabs(g[0] * g[3])) det_g = det_entries;
  }
  for (const Pattern& pt : d.patterns) {
    if (pt.count == 0) continue;
    double wm[QMAX * QMAX], detk;
    if (!between(pt, g, det_g, q, wm, detk)) return false;
    w.logdet_k += pt.count * std::log(detk);
    // W C (q × p), C = (Z'Z)⁺ Z'X the coefficients of X on Z
    double wc[QMAX * PMAX];
    for (int i = 0; i < q; ++i)
      for (int b = 0; b < p; ++b) {
        double t = 0;
        for (int j = 0; j < q; ++j) t += wm[i * q + j] * pt.cx[j * p + b];
        wc[i * p + b] = t;
      }
    const double c = pt.count;
    for (int a = 0; a < p; ++a) {
      for (int b = a; b < p; ++b) {
        double t = 0;
        for (int i = 0; i < q; ++i) t += pt.cx[i * p + a] * wc[i * p + b];
        w.m[a * p + b] += c * (pt.xx[a * p + b] + t);
      }
      double t = 0;
      for (int i = 0; i < q; ++i) t += wc[i * p + a] * pt.b[i];
      w.v[a] += pt.xy[a] + t;
    }
    double tr = 0;
    for (int i = 0; i < q; ++i)
      for (int j = 0; j < q; ++j) tr += wm[i * q + j] * pt.bb[j * q + i];
    s += pt.yy + tr;
  }
  // the same sums for clusters by arm and size: Z'Z = n, C = x and nothing of X within the cluster
  for (const Arm& arm : d.arms) {
    const double* x = arm.x;
    for (const Arm::Size& sz : arm.sizes) {
      if (sz.count == 0) continue;
      const double detk = 1 + sz.n * g[0];
      if (!(detk > 0)) return false;
      const double wm = sz.n / detk;
      w.logdet_k += sz.count * std::log(detk);
      double wc[PMAX];
      for (int b = 0; b < p; ++b) wc[b] = wm * x[b];
      for (int a = 0; a < p; ++a) {
        for (int b = a; b < p; ++b) w.m[a * p + b] += sz.count * (x[a] * wc[b]);
        w.v[a] += wc[a] * sz.b;
      }
      s += sz.yy + wm * sz.bb;
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
  if (q == 1) {
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
    if (d.two_rows) {
      // every point of a ridge from a singular D to σ² = 0 has the same V and criterion: the fits stopped anywhere
      // on it, and the Hessian in (σ², D) was singular, which left Satterthwaite's degrees of freedom to rounding
      // or to the given ones (in 69 of 200 studies with visits at 0 and 1 and 60 subjects a group); held at rank 1,
      // they are those of the covariance of the two measurements
      x[2] = 0;
      nelder_mead(
          [&](const double* th) {
            const double t[3] = {th[0], th[1], 0};
            return profiled(d, t, w);
          },
          2, x, 0.5, 4000, out.evaluations, scale);
    } else {
      nelder_mead([&](const double* th) { return profiled(d, th, w); }, 3, x, 0.5, 4000, out.evaluations, scale);
    }
    std::copy(x, x + 3, out.theta);
  }
  // a fit whose X'V⁻¹X has lost its digits has no estimate of β: with every subject of a group seen once, nothing
  // tells that group's slope apart, and its estimate came out exactly 0 or 1 with a standard error of 10⁷; with two
  // visits and one subject of a group measured twice and two of the other, the criterion has no minimum (it falls
  // without bound as the covariance of the two measurements turns singular), and the fits ended with σ² 10⁻²² of
  // the random effects' variances
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
  // A random intercept estimated at zero is held there (φ = σ² only, which gives N − p).
  //
  // With random slopes, D is the covariance of the intercept and slope at a time t0 + δ (δ in the design's time
  // unit, D' = T D T' with T = [[1, δ], [0, 1]]): at the design's own origin when the fitted D is far from singular
  // (δ = 0), and otherwise at the time where the fitted intercept and slope are uncorrelated (δ = −D12/D22), which
  // changes nothing in exact arithmetic (the criterion and V are the same functions of D'), but keeps the
  // finite-difference steps away from where D stops being a covariance matrix: where the subjects' lines all but
  // cross away from the design's origin (fitted correlation 0.99), a step of 10⁻⁴ is a large part of the distance to
  // singularity, and the Hessian came out not positive definite in every study (at a correlation of 0.9, the degrees
  // of freedom were 0.1% off).  A fit on the boundary, with a singular D = w w' (the subjects' lines all cross at one
  // time, or their slopes do not vary: a variance in some direction below the steps), is taken on that boundary, in
  // (σ², w), where V and the criterion are even in w: the variance in the missing direction is held at 0, as a
  // random intercept estimated at zero is, and as lmerTest's derivatives in the factor Λ do (its degrees of freedom
  // are these to four digits or more on the same fits).  The steps in D went through the boundary, into matrices
  // that are not covariances, which left the result to rounding (2.3 or 13.3 for lmerTest's 149.5 or 90.4; a
  // 40-digit limit can be negative), or to the given degrees of freedom when the Hessian was not positive definite.
  // Only where both ways fail are the given degrees of freedom used.
  enum class Vary { Full, Rank1 };  // D' in full, or D = w w' (φ = σ², w0, w1)
  const double s2 = f.sigma2;
  double tau = 1;
  if (q == 2) {
    double tt = 0, nn = 0;
    for (const Pattern& pt : d.patterns) tt += pt.count * pt.zz[3], nn += pt.count * pt.zz[0];
    if (tt > 0 && nn > 0) tau = std::sqrt(tt / nn);
  }
  auto satterthwaite = [&](int k, const double* phi, Vary vary, double shift) {
    const bool full = q == 2 && vary == Vary::Full;
    auto unpack = [&](const double* x, double& sigma2, double* g, double& det_g) {
      sigma2 = x[0];
      det_g = std::numeric_limits<double>::quiet_NaN();
      if (q == 1) {
        g[0] = k > 1 ? x[1] / sigma2 : 0;
        return;
      }
      // D' and then D = T⁻¹ D' T⁻ᵀ, |D| = |D'|
      const double w0 = vary == Vary::Rank1 && k > 1 ? x[1] : 0, w1 = vary == Vary::Rank1 && k > 1 ? x[2] : 0;
      const double d00 = full ? x[1] : w0 * w0, d01 = full ? x[2] : w0 * w1, d11 = full ? x[3] : w1 * w1;
      if (shift == 0) {
        g[0] = d00 / sigma2;
        g[1] = g[2] = d01 / sigma2;
        g[3] = d11 / sigma2;
        if (!full) det_g = 0;
        return;
      }
      g[0] = (d00 - shift * (2 * d01 - shift * d11)) / sigma2;
      g[1] = g[2] = (d01 - shift * d11) / sigma2;
      g[3] = d11 / sigma2;
      det_g = (d00 * d11 - d01 * d01) / sigma2 / sigma2;
    };
    auto crit = [&](const double* x) {
      double sigma2, g[QMAX * QMAX], det_g;
      unpack(x, sigma2, g, det_g);
      return criterion(d, sigma2, g, det_g);
    };
    auto cvar = [&](const double* x) {
      double sigma2, g[QMAX * QMAX], det_g, cov[PMAX * PMAX];
      unpack(x, sigma2, g, det_g);
      if (!beta_covariance(d, sigma2, g, cov, det_g)) return std::numeric_limits<double>::quiet_NaN();
      double v = 0;
      for (int a = 0; a < p; ++a)
        for (int b = 0; b < p; ++b) v += c[a] * cov[a * p + b] * c[b];
      return v;
    };
    // finite-difference steps on the variance components standardized by the random slope's covariate (its root
    // mean square τ: τ·D12 and τ²·D22 have the outcome's squared units), so that the degrees of freedom do not
    // depend on the time unit
    const double unit[4] = {1, 1, q == 2 ? tau : 1, tau * tau};
    // a variance's step is relative to itself (a step relative to the other variance as well would be most of it
    // when that one is 10⁴ times larger), the covariance's to the geometric mean of the two variances, and w's to
    // its length
    double h[4];
    const double dscale = full ? std::sqrt(std::fabs(phi[1] * phi[3] * unit[3])) : 0;
    const double length = vary == Vary::Rank1 && k > 1 ? std::hypot(phi[1], phi[2] * tau) : 0;
    for (int j = 0; j < k; ++j) {
      const double size = std::fabs(phi[j]) * unit[j];
      h[j] = j == 0                    ? 1e-4 * s2
             : vary == Vary::Rank1 ? 1e-4 * (length + 1e-1 * std::sqrt(s2)) / unit[j]
                                   : 1e-4 * ((full && j == 2 ? size + dscale : 2 * size) + 1e-2 * s2) / unit[j];
    }
    // the derivatives are taken in units of the steps and of V, each a power of 2 (hs = h 2^e, V 2^m), which keeps
    // every product and quotient below exact to the last bit and away from under- and overflow: h², and V², pass
    // below the smallest double once the residual variance is 10⁻¹⁵⁰ of the random effects' in the analysis units
    double hess[16], grad[4], x[4], hs[4];
    for (int j = 0; j < k; ++j) hs[j] = std::ldexp(h[j], -std::ilogb(h[j]));
    const int m = var > 0 ? -std::ilogb(var) : 0;
    const double f0 = crit(phi);
    bool ok = std::isfinite(f0) && var > 0;
    for (int a = 0; a < k && ok; ++a) {
      std::copy(phi, phi + k, x);
      x[a] = phi[a] + h[a];
      const double fp = crit(x), vp = cvar(x);
      x[a] = phi[a] - h[a];
      const double fm = crit(x), vm = cvar(x);
      hess[a * k + a] = (fp - 2 * f0 + fm) / (hs[a] * hs[a]);
      grad[a] = std::ldexp(vp - vm, m) / (2 * hs[a]);
      for (int b = 0; b < a; ++b) {
        double fpp, fpm, fmp, fmm;
        std::copy(phi, phi + k, x);
        x[a] = phi[a] + h[a], x[b] = phi[b] + h[b];
        fpp = crit(x);
        x[b] = phi[b] - h[b];
        fpm = crit(x);
        x[a] = phi[a] - h[a];
        fmm = crit(x);
        x[b] = phi[b] + h[b];
        fmp = crit(x);
        hess[a * k + b] = hess[b * k + a] = (fpp - fpm - fmp + fmm) / (4 * hs[a] * hs[b]);
      }
      ok = std::isfinite(hess[a * k + a]) && std::isfinite(grad[a]);
    }
    if (!ok || !cholesky(hess, k)) return std::numeric_limits<double>::quiet_NaN();
    double y[4];
    std::copy(grad, grad + k, y);
    chol_solve(hess, k, y);  // H⁻¹ ∇V
    double quad = 0;
    for (int a = 0; a < k; ++a) quad += grad[a] * y[a];
    const double v = std::ldexp(var, m);
    return v * v / quad;  // 2 V² / (∇V' · 2 H⁻¹ · ∇V)
  };
  double df = std::numeric_limits<double>::quiet_NaN();
  if (q == 1) {
    const double phi[2] = {s2, s2 * f.g[0]};
    df = satterthwaite(f.boundary ? 1 : 2, phi, Vary::Full, 0);
  } else {
    // on the boundary when the slope's variance (in units of τ²), or else the intercept's at the crossing time
    // (|D|/D22, with |G| from Λ: from G's entries it is the small difference of large numbers), is below the steps'
    // floor of 10⁻⁶ σ²
    const double det_g = f.theta[0] * f.theta[2] * (f.theta[0] * f.theta[2]);
    const bool flat_slope = !(f.g[3] * tau * tau > 1e-6), flat = flat_slope || !(det_g > 1e-6 * f.g[3]);
    if (flat) {
      // w from Λ's first column (G = w w' when θ2 = 0), or the slope's alone when the intercept does not vary; σ²
      // alone when neither does
      const double sd = std::sqrt(s2), w0 = f.theta[0] != 0 ? sd * f.theta[0] : 0;
      const double w1 = f.theta[0] != 0 ? sd * f.theta[1] : sd * std::hypot(f.theta[1], f.theta[2]);
      const double held[3] = {s2, w0, w1};
      df = satterthwaite(flat_slope && !(f.g[0] > 1e-6) ? 1 : 3, held, Vary::Rank1, 0);
    }
    if (!(std::isfinite(df) && df > 0)) {
      const double shift = !flat && det_g < 0.1 * (f.g[0] * f.g[3]) ? -f.g[1] / f.g[3] : 0;
      const double phi[4] = {s2, shift != 0 ? s2 * (det_g / f.g[3]) : s2 * f.g[0], shift != 0 ? 0 : s2 * f.g[1],
                             s2 * f.g[3]};
      df = satterthwaite(4, phi, Vary::Full, shift);
    }
  }
  if (std::isfinite(df) && df > 0) {
    out.df = std::fmin(df, d.observations() - p);
    return out;
  }
  out.df = fixed_df;
  out.fallback = true;
  return out;
}

}  // namespace ssc::sim::lmm
