#include "sim/lmm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace ssc::sim::lmm {

namespace {

/// Whether a 2 × 2 Z'Z has full rank: one whose determinant is negligible next to the product of its diagonal
/// (one visit, or visits at the same time) has rank 1.
bool full_rank(const double* a) { return a[0] * a[3] - a[1] * a[2] > 1e-12 * a[0] * a[3]; }

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
    const double tr = a[0] + a[3];
    for (int i = 0; i < 4; ++i) zi[i] = tr > 0 ? a[i] / (tr * tr) : 0;
    for (int b = 0; b < p; ++b) {
      cx[b] = tr > 0 ? (a[0] * zx[b] + a[1] * zx[p + b]) / (tr * tr) : 0;
      cx[p + b] = tr > 0 ? (a[2] * zx[b] + a[3] * zx[p + b]) / (tr * tr) : 0;
    }
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
  if (rows <= rank(zz, q)) return;
  for (long r = 0; r < rows; ++r) {
    double res[PMAX];
    within_row(x + r * p, z + r * q, cx, p, q, res);
    for (int a = 0; a < p; ++a)
      for (int b = 0; b < p; ++b) xx[a * p + b] += res[a] * res[b];
  }
}

void Pattern::add_cluster(int p, int q, const double* x, const double* z, const double* y) {
  const double none[QMAX] = {};
  add_cluster(p, q, x, z, none, y);
}

void Pattern::add_cluster(int p, int q, const double* x, const double* z, const double* u, const double* e) {
  // e's coefficients on Z give b = u + (Z'Z)⁺ Z'e and e's residuals; X'(I − P)e from X's residuals
  double ze[QMAX] = {}, coef_e[QMAX] = {};
  for (long r = 0; r < n; ++r)
    for (int i = 0; i < q; ++i) ze[i] += z[r * q + i] * e[r];
  for (int i = 0; i < q; ++i)
    for (int j = 0; j < q; ++j) coef_e[i] += zi[i * q + j] * ze[j];
  double sx[PMAX] = {}, sb[QMAX], syy = 0;
  if (n > rank(zz, q)) {
    for (long r = 0; r < n; ++r) {
      double v = e[r], res[PMAX];
      for (int i = 0; i < q; ++i) v -= z[r * q + i] * coef_e[i];
      within_row(x + r * p, z + r * q, cx, p, q, res);
      for (int a = 0; a < p; ++a) sx[a] += res[a] * v;
      syy += v * v;
    }
  }
  for (int i = 0; i < q; ++i) sb[i] = u[i] + coef_e[i];
  add_sums(p, q, sx, sb, syy);
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
  return s;
}

double Data::clusters() const {
  double s = 0;
  for (const auto& pt : patterns) s += pt.count;
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

/// For a cluster with Z'Z = A (and A⁺ = ai): W = K⁻¹ A, the between-cluster part of V*⁻¹ in the coordinates b,
/// and |K| = |I + A G|; false unless V* is positive definite.  Computed so that nothing cancels: for a 2 × 2 A
/// of full rank, |K| = 1 + tr(AG) + |A||G| and W = (A⁻¹ + G)⁻¹, whose determinant is |K|/|A|; for a
/// 1 × 1 A, or one of rank 1, W = A / (1 + tr(AG)).
bool between(const double* a, const double* ai, const double* g, int q, double* w, double& detk) {
  const double tr_ag = q == 1 ? a[0] * g[0] : a[0] * g[0] + a[1] * g[2] + a[2] * g[1] + a[3] * g[3];
  if (q == 1 || !full_rank(a)) {
    detk = 1 + tr_ag;
    if (!(detk > 0)) return false;
    for (int i = 0; i < q * q; ++i) w[i] = a[i] / detk;
    return true;
  }
  const double det_a = a[0] * a[3] - a[1] * a[2];
  detk = 1 + tr_ag + det_a * (g[0] * g[3] - g[1] * g[2]);
  const double e00 = ai[0] + g[0], e11 = ai[3] + g[3], e01 = 0.5 * (ai[1] + g[1] + ai[2] + g[2]);
  if (!(detk > 0) || !(e00 > 0)) return false;
  const double det_e = detk / det_a;
  w[0] = e11 / det_e, w[1] = w[2] = -e01 / det_e, w[3] = e00 / det_e;
  return true;
}

bool gls(const Data& d, const double* g, Gls& w) {
  const int p = d.p, q = std::min(d.q, QMAX);  // bounded, so that the compiler unrolls the loops over Z's columns
  std::fill(w.m, w.m + p * p, 0.0);
  std::fill(w.v, w.v + p, 0.0);
  double s = 0;
  w.logdet_k = 0;
  for (const Pattern& pt : d.patterns) {
    if (pt.count == 0) continue;
    double wm[QMAX * QMAX], detk;
    if (!between(pt.zz, pt.zi, g, q, wm, detk)) return false;
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

/// Nelder–Mead minimization in `n` dimensions (n ≤ 3), restarted once from its result.
template <class F>
void nelder_mead(F&& f, int n, double* x, double step, int max_evals, int& evals) {
  for (int round = 0; round < 2; ++round) {
    double pts[4][3], val[4];
    for (int i = 0; i <= n; ++i) {
      for (int j = 0; j < n; ++j) pts[i][j] = x[j];
      if (i > 0) pts[i][i - 1] += step * std::fmax(1.0, std::fabs(x[i - 1]));
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

double criterion(const Data& d, double sigma2, const double* g) {
  Gls w;
  if (!(sigma2 > 0) || !gls(d, g, w) || !(w.rss > 0)) return std::numeric_limits<double>::quiet_NaN();
  const double r = d.observations() - d.p;
  return r * (LOG_2PI + std::log(sigma2)) + w.logdet_k + w.logdet_m + w.rss / sigma2;
}

bool beta_covariance(const Data& d, double sigma2, const double* g, double* cov) {
  Gls w;
  if (!gls(d, g, w)) return false;
  chol_inverse(w.m, d.p, cov);
  for (int i = 0; i < d.p * d.p; ++i) cov[i] *= sigma2;
  return true;
}

Fit fit(const Data& d, const double* theta0) {
  Fit out;
  const int p = d.p, q = d.q;
  if (p < 1 || p > PMAX || q < 1 || q > QMAX || d.observations() <= p) return out;
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
    nelder_mead([&](const double* th) { return profiled(d, th, w); }, 3, x, 0.5, 4000, out.evaluations);
    std::copy(x, x + 3, out.theta);
  }
  out.criterion = profiled(d, out.theta, w);
  if (!std::isfinite(out.criterion)) return out;
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
  const double s2 = f.sigma2;
  double phi[4];
  int k = 1;
  phi[0] = s2;
  if (q == 1) {
    if (!f.boundary) phi[k++] = s2 * f.g[0];
  } else {
    phi[k++] = s2 * f.g[0];
    phi[k++] = s2 * f.g[1];
    phi[k++] = s2 * f.g[3];
  }
  auto unpack = [&](const double* x, double& sigma2, double* g) {
    sigma2 = x[0];
    if (q == 1) {
      g[0] = k > 1 ? x[1] / sigma2 : 0;
    } else {
      g[0] = x[1] / sigma2;
      g[1] = g[2] = x[2] / sigma2;
      g[3] = x[3] / sigma2;
    }
  };
  auto crit = [&](const double* x) {
    double sigma2, g[QMAX * QMAX];
    unpack(x, sigma2, g);
    return criterion(d, sigma2, g);
  };
  auto cvar = [&](const double* x) {
    double sigma2, g[QMAX * QMAX], cov[PMAX * PMAX];
    unpack(x, sigma2, g);
    if (!beta_covariance(d, sigma2, g, cov)) return std::numeric_limits<double>::quiet_NaN();
    double v = 0;
    for (int a = 0; a < p; ++a)
      for (int b = 0; b < p; ++b) v += c[a] * cov[a * p + b] * c[b];
    return v;
  };
  // finite-difference steps on the variance components standardized by the random slope's covariate (its root
  // mean square τ: τ·D12 and τ²·D22 have the outcome's squared units), so that the degrees of freedom do not
  // depend on the time unit
  double tau = 1;
  if (q == 2) {
    double tt = 0, nn = 0;
    for (const Pattern& pt : d.patterns) tt += pt.count * pt.zz[3], nn += pt.count * pt.zz[0];
    if (tt > 0 && nn > 0) tau = std::sqrt(tt / nn);
  }
  const double unit[4] = {1, 1, q == 2 ? tau : 1, tau * tau};
  // a variance's step is relative to itself (a step relative to the other variance as well would be most of it
  // when that one is 10⁴ times larger), the covariance's to the geometric mean of the two variances
  double h[4];
  const double dscale = q == 2 ? std::sqrt(std::fabs(phi[1] * phi[3] * unit[3])) : 0;
  for (int j = 0; j < k; ++j) {
    const double size = std::fabs(phi[j]) * unit[j];
    h[j] = j == 0 ? 1e-4 * s2 : 1e-4 * ((j == 2 ? size + dscale : 2 * size) + 1e-2 * s2) / unit[j];
  }
  double hess[16], grad[4], x[4];
  const double f0 = crit(phi);
  bool ok = std::isfinite(f0);
  for (int a = 0; a < k && ok; ++a) {
    std::copy(phi, phi + k, x);
    x[a] = phi[a] + h[a];
    const double fp = crit(x), vp = cvar(x);
    x[a] = phi[a] - h[a];
    const double fm = crit(x), vm = cvar(x);
    hess[a * k + a] = (fp - 2 * f0 + fm) / (h[a] * h[a]);
    grad[a] = (vp - vm) / (2 * h[a]);
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
      hess[a * k + b] = hess[b * k + a] = (fpp - fpm - fmp + fmm) / (4 * h[a] * h[b]);
    }
    ok = std::isfinite(hess[a * k + a]) && std::isfinite(grad[a]);
  }
  if (ok && cholesky(hess, k)) {
    double y[4];
    std::copy(grad, grad + k, y);
    chol_solve(hess, k, y);  // H⁻¹ ∇V
    double quad = 0;
    for (int a = 0; a < k; ++a) quad += grad[a] * y[a];
    const double df = var * var / quad;  // 2 V² / (∇V' · 2 H⁻¹ · ∇V)
    if (std::isfinite(df) && df > 0) {
      out.df = std::fmin(df, d.observations() - p);
      return out;
    }
  }
  out.df = fixed_df;
  out.fallback = true;
  return out;
}

}  // namespace ssc::sim::lmm
