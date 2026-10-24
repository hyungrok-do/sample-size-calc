#include "sim/glmm.hpp"

#include <algorithm>
#include <cmath>

namespace ssc::sim::glmm {

namespace {

constexpr double LOG_SQRT_PI = 0.57236494292470008707;  // log √π
constexpr double SQRT2 = 1.41421356237309504880;

double expit(double t) { return t >= 0 ? 1 / (1 + std::exp(-t)) : std::exp(t) / (1 + std::exp(t)); }

struct Rule {
  std::vector<double> x, w, wexp, logw;  // nodes, weights, weights × exp(x²) and its log
};

Rule make_rule(int n) {
  Rule r;
  gauss_hermite(n, r.x, r.w);
  r.wexp.resize(r.w.size());
  r.logw.resize(r.w.size());
  for (std::size_t k = 0; k < r.w.size(); ++k) {
    r.wexp[k] = r.w[k] * std::exp(r.x[k] * r.x[k]);
    r.logw[k] = std::log(r.wexp[k]);
  }
  return r;
}

const Rule& rule(int n) {
  static const Rule r21 = make_rule(21), r40 = make_rule(40);
  if (n == 21) return r21;
  if (n == 40) return r40;
  thread_local Rule other;
  thread_local int other_n = 0;
  if (other_n != n) {
    other = make_rule(n);
    other_n = n;
  }
  return other;
}

/// Mode ẑ of g(z) = s·t − m·log(1 + e^t) − z²/2 with t = η + σz (strictly concave), and −g''(ẑ).
/// Newton's method safeguarded by bisection (plain Newton can cycle when σ²m is large): g'(z) =
/// σ(s − m·expit(t)) − z decreases, and is positive below min(σs, σ(s − m)) and negative above max(σs, σ(s − m)).
double mode(double eta, double m, double s, double sigma, double start, double& curvature) {
  double lo = std::min(sigma * s, sigma * (s - m)), hi = std::max(sigma * s, sigma * (s - m));
  double z = std::max(lo, std::min(hi, start));
  for (int it = 0; it < 200 && lo < hi; ++it) {
    const double p = expit(eta + sigma * z);
    const double g1 = sigma * (s - m * p) - z, g2 = -sigma * sigma * m * p * (1 - p) - 1;
    if (g1 > 0)
      lo = z;
    else if (g1 < 0)
      hi = z;
    else
      break;
    double next = z - g1 / g2;
    if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
    const bool done = std::fabs(next - z) < 1e-12 * (1 + std::fabs(z));
    z = next;
    if (done) break;
  }
  const double p = expit(eta + sigma * z);
  curvature = sigma * sigma * m * p * (1 - p) + 1;
  return z;
}

/// g(z) and expit(t) at once, sharing the exponential.
double g_and_p(double eta, double m, double s, double sigma, double z, double& p) {
  const double t = eta + sigma * z, e = std::exp(-std::fabs(t));
  p = t >= 0 ? 1 / (1 + e) : e / (1 + e);
  return s * t - m * (std::fmax(t, 0.0) + std::log1p(e)) - 0.5 * z * z;
}

/// One cluster's log-likelihood and its derivatives in (β0, β1, σ), with the quadrature nodes centered
/// at the mode `z0` and scaled by `tau` = (−g''(z0))^(−1/2).  The gradient is exact: it includes the terms
/// from the nodes moving with the parameters (which vanish for an exact integral but not for the
/// quadrature).  The Hessian holds the nodes fixed, which is accurate up to the quadrature error.
double cluster_terms(const Cluster& c, double b0, double b1, double sigma, double z0, double tau, const Rule& r,
                     double* grad, double* hess) {
  const double eta = b0 + b1 * c.x;
  const std::size_t q = std::min<std::size_t>(r.x.size(), 64);
  double gz[64], zk[64], pz[64];
  double gmax = -INFINITY;
  for (std::size_t k = 0; k < q; ++k) {
    zk[k] = z0 + SQRT2 * tau * r.x[k];
    gz[k] = g_and_p(eta, c.m, c.s, sigma, zk[k], pz[k]);
    gmax = std::max(gmax, gz[k] + r.logw[k]);
  }
  // ∫ e^g dz / √(2π) ≈ (τ/√π) Σ w_k e^(x_k²) e^(g(z_k)); u_k is a term of the sum divided by e^gmax
  double sum = 0;
  double u[64];
  for (std::size_t k = 0; k < q; ++k) sum += (u[k] = r.wexp[k] * std::exp(gz[k] - gmax));
  const double ll = gmax + std::log(sum) + std::log(tau) - LOG_SQRT_PI;
  if (grad) {
    double gm[3] = {0, 0, 0}, h[9] = {0};
    for (std::size_t k = 0; k < q; ++k) {
      const double pk = u[k] / sum;
      const double p = pz[k];
      const double a = c.s - c.m * p, cc = -c.m * p * (1 - p), z = zk[k];
      const double gk[3] = {a, a * c.x, a * z};
      const double hk[9] = {cc, cc * c.x, cc * z, cc * c.x, cc * c.x * c.x, cc * c.x * z, cc * z, cc * c.x * z,
                            cc * z * z};
      for (int i = 0; i < 3; ++i) {
        gm[i] += pk * gk[i];
        for (int j = 0; j < 3; ++j) h[3 * i + j] += pk * (hk[3 * i + j] + gk[i] * gk[j]);
      }
    }
    // the node terms: ∂ℓ/∂z0 · dz0/dθ + ∂ℓ/∂τ · dτ/dθ, with dz0/dθ from g'(z0) = 0 and τ = c^(−1/2), where
    // c = σ²m·p(1 − p) + 1 at the mode (g' and g'' are derivatives in z)
    double dl_dz0 = 0, dl_dtau = 1 / tau;
    for (std::size_t k = 0; k < q; ++k) {
      const double g1 = sigma * (c.s - c.m * pz[k]) - zk[k], pk = u[k] / sum;
      dl_dz0 += pk * g1;
      dl_dtau += pk * g1 * SQRT2 * r.x[k];
    }
    const double p0 = expit(eta + sigma * z0), v0 = p0 * (1 - p0), curv = 1 / (tau * tau);
    const double w = sigma * sigma * c.m * v0 * (1 - 2 * p0);  // ∂c/∂t
    const double dg1[3] = {-sigma * c.m * v0, -sigma * c.m * v0 * c.x, c.s - c.m * p0 - sigma * c.m * v0 * z0};
    const double dc[3] = {w, w * c.x, 2 * sigma * c.m * v0 + w * z0};  // ∂c/∂θ with z0 held fixed
    for (int i = 0; i < 3; ++i) {
      const double dz0 = dg1[i] / curv, dtau = -0.5 * tau * (dc[i] + w * sigma * dz0) / curv;
      grad[i] += c.count * (gm[i] + dl_dz0 * dz0 + dl_dtau * dtau);
      for (int j = 0; j < 3; ++j) hess[3 * i + j] += c.count * (h[3 * i + j] - gm[i] * gm[j]);
    }
  }
  return ll;
}

/// Total log-likelihood (and derivatives), re-centering every cluster's nodes; `modes` keeps the modes.
double evaluate(const std::vector<Cluster>& data, double b0, double b1, double sigma, const Rule& r,
                std::vector<double>& modes, double* grad, double* hess, Poll poll = {}) {
  if (grad) {
    std::fill(grad, grad + 3, 0.0);
    std::fill(hess, hess + 9, 0.0);
  }
  double total = 0;
  for (std::size_t i = 0; i < data.size(); ++i) {
    if (poll.call && i % 4096 == 0) poll.call(poll.arg);
    const Cluster& c = data[i];
    double curvature;
    modes[i] = mode(b0 + b1 * c.x, c.m, c.s, sigma, modes[i], curvature);
    total += c.count * cluster_terms(c, b0, b1, sigma, modes[i], 1 / std::sqrt(curvature), r, grad, hess);
  }
  return total;
}

/// Solve the 3 × 3 (or the leading 2 × 2) system A d = b by Cholesky; false unless A is positive definite.
bool solve_spd(const double* a_in, const double* b, double* d, int n) {
  double a[9];
  std::copy(a_in, a_in + 9, a);
  double l[9] = {0};
  for (int j = 0; j < n; ++j) {
    double s = a[3 * j + j];
    for (int k = 0; k < j; ++k) s -= l[3 * j + k] * l[3 * j + k];
    if (!(s > 0)) return false;
    l[3 * j + j] = std::sqrt(s);
    for (int i = j + 1; i < n; ++i) {
      double t = a[3 * i + j];
      for (int k = 0; k < j; ++k) t -= l[3 * i + k] * l[3 * j + k];
      l[3 * i + j] = t / l[3 * j + j];
    }
  }
  double y[3];
  for (int i = 0; i < n; ++i) {
    double s = b[i];
    for (int k = 0; k < i; ++k) s -= l[3 * i + k] * y[k];
    y[i] = s / l[3 * i + i];
  }
  for (int i = n - 1; i >= 0; --i) {
    double s = y[i];
    for (int k = i + 1; k < n; ++k) s -= l[3 * k + i] * d[k];
    d[i] = s / l[3 * i + i];
  }
  return true;
}

}  // namespace

void gauss_hermite(int n, std::vector<double>& x, std::vector<double>& w) {
  // Newton's method on the orthonormal Hermite polynomials with the classical starting values
  // (Press et al., Numerical Recipes, gauher); nodes in decreasing order.
  x.assign(static_cast<std::size_t>(n), 0.0);
  w.assign(static_cast<std::size_t>(n), 0.0);
  const double pim4 = 0.75112554446494248286;  // π^(−1/4)
  double z = 0, pp = 0;
  for (int i = 0; i < (n + 1) / 2; ++i) {
    if (i == 0)
      z = std::sqrt(2.0 * n + 1) - 1.85575 * std::pow(2.0 * n + 1, -0.16667);
    else if (i == 1)
      z -= 1.14 * std::pow(static_cast<double>(n), 0.426) / z;
    else if (i == 2)
      z = 1.86 * z - 0.86 * x[0];
    else if (i == 3)
      z = 1.91 * z - 0.91 * x[1];
    else
      z = 2 * z - x[static_cast<std::size_t>(i - 2)];
    for (int it = 0; it < 100; ++it) {
      double p1 = pim4, p2 = 0;
      for (int j = 1; j <= n; ++j) {
        const double p3 = p2;
        p2 = p1;
        p1 = z * std::sqrt(2.0 / j) * p2 - std::sqrt((j - 1.0) / j) * p3;
      }
      pp = std::sqrt(2.0 * n) * p2;
      const double z1 = z;
      z = z1 - p1 / pp;
      if (std::fabs(z - z1) <= 1e-15 * std::max(1.0, std::fabs(z))) break;
    }
    x[static_cast<std::size_t>(i)] = z;
    x[static_cast<std::size_t>(n - 1 - i)] = -z;
    w[static_cast<std::size_t>(i)] = w[static_cast<std::size_t>(n - 1 - i)] = 2 / (pp * pp);
  }
}

double cluster_loglik(double eta, double m, double s, double sigma, int nodes) {
  double curvature;
  const double z0 = mode(eta, m, s, sigma, 0, curvature);
  Cluster c;
  c.x = 0, c.m = m, c.s = s;
  return cluster_terms(c, eta, 0, sigma, z0, 1 / std::sqrt(curvature), rule(nodes), nullptr, nullptr);
}

double loglik(const std::vector<Cluster>& data, double b0, double b1, double sigma, int nodes) {
  std::vector<double> modes(data.size(), 0.0);
  return evaluate(data, b0, b1, sigma, rule(nodes), modes, nullptr, nullptr);
}

Fit fit(const std::vector<Cluster>& data, double b0, double b1, double sigma, int nodes, Poll poll) {
  // Newton's method on the adaptive-quadrature log-likelihood, whose nodes are re-centered at the modes
  // for every parameter value tried.  The derivatives are those of the quadrature sum with its nodes held
  // where they are, which equal the log-likelihood's up to the quadrature error.  (Holding the nodes
  // fixed over several steps does not work for large clusters: each cluster's integrand is then so
  // narrow that the estimates cannot move away from the nodes.)
  Fit out;
  auto fail = [&](Failure why) {
    out.ok = false;
    out.failure = why;
    return out;
  };
  const Rule& r = rule(nodes);
  const std::size_t n = data.size();
  // σ is not identified when no cluster has more than one subject and the covariate takes at most two values (two
  // groups): each cluster then enters the likelihood only through its group's population-averaged proportion,
  // which β and σ trade off along a ridge of maxima, and the standard error of β1 is that of the ridge.  σ is held
  // at 0 there, which is ordinary logistic regression.
  bool hold = true;
  std::vector<double> xs;
  for (std::size_t i = 0; i < n && hold; ++i) {
    if (std::find(xs.begin(), xs.end(), data[i].x) == xs.end()) xs.push_back(data[i].x);
    hold = data[i].m <= 1 && xs.size() <= 2;
  }
  // with σ held, its derivatives are left out (its Newton step is then 0)
  auto pin = [&](double* g, double* h) {
    if (!hold) return;
    g[2] = 0;
    h[2] = h[5] = h[6] = h[7] = 0;
    h[8] = -1;
  };
  std::vector<double> modes(n, 0.0), trial_modes(n);
  double th[3] = {b0, b1, hold ? 0 : (sigma <= 0 ? 0.3 : sigma)};
  double grad[3], hess[9];
  double ll = evaluate(data, th[0], th[1], th[2], r, modes, grad, hess, poll);
  pin(grad, hess);
  if (!std::isfinite(ll)) return fail(Failure::Arithmetic);
  bool converged = false;
  for (int it = 0; it < 100 && !converged; ++it) {
    ++out.iterations;
    // ascent direction d = (−H + λI)⁻¹ ∇ℓ, with λ raised until −H + λI is positive definite
    double neg[9], d[3];
    bool solved = false;
    for (double lambda = 0; !solved && lambda < 1e8; lambda = lambda == 0 ? 1e-8 : lambda * 10) {
      for (int i = 0; i < 9; ++i) neg[i] = -hess[i];
      for (int i = 0; i < 3; ++i) neg[4 * i] += lambda;
      solved = solve_spd(neg, grad, d, 3);
    }
    if (!solved) return fail(Failure::Arithmetic);
    const double decrement = grad[0] * d[0] + grad[1] * d[1] + grad[2] * d[2];
    if (decrement < 1e-12 * (1 + std::fabs(ll))) {
      // a last Newton step: too small for the log-likelihood to show its gain, it still squares the
      // remaining error
      for (int i = 0; i < 3; ++i) th[i] += d[i];
      ll = evaluate(data, th[0], th[1], th[2], r, modes, grad, hess, poll);
      pin(grad, hess);
      converged = std::isfinite(ll);
      break;
    }
    // the trial points are evaluated with their derivatives, which the accepted one then supplies
    double step = 1, trial[3], next = -INFINITY, tgrad[3], thess[9];
    for (int half = 0; half < 50; ++half, step /= 2) {
      for (int i = 0; i < 3; ++i) trial[i] = th[i] + step * d[i];
      std::copy(modes.begin(), modes.end(), trial_modes.begin());
      next = evaluate(data, trial[0], trial[1], trial[2], r, trial_modes, tgrad, thess, poll);
      pin(tgrad, thess);
      if (std::isfinite(next) && next >= ll) break;
    }
    if (!(std::isfinite(next) && next >= ll)) {
      // no progress possible: converged up to rounding (and the quadrature error)
      converged = decrement < 1e-6 * (1 + std::fabs(ll));
      break;
    }
    std::copy(trial, trial + 3, th);
    std::copy(tgrad, tgrad + 3, grad);
    std::copy(thess, thess + 9, hess);
    modes.swap(trial_modes);
    ll = next;
    // β is left free: separated data, whose β has no finite estimate, are the caller's to catch, and otherwise it
    // has one, however large (a log odds ratio of 30 between groups with 3e-7 and 1 − 3e-7 of 10⁷ subjects a cluster)
    if (std::fabs(th[2]) > SIGMA_MAX) return fail(Failure::RunOff);
  }
  if (!converged) return fail(Failure::Arithmetic);
  out.ok = true;
  out.b0 = th[0];
  out.b1 = th[1];
  out.sigma = std::fabs(th[2]);
  out.loglik = ll;
  out.boundary = out.sigma < 1e-5;
  // standard error of β1 from the observed information, the Hessian of the log-likelihood by forward
  // differences of its exact gradient: the Hessian with the nodes held fixed misses their movement, which
  // understates the standard error when the clusters' integrands are far from Gaussian (by a quarter at a
  // latent ICC of 0.9).  With σ at (or near) zero, σ is held fixed there, where the integrands do not
  // depend on the random effect and the fixed-node Hessian is exact.
  double info[9];
  for (int i = 0; i < 9; ++i) info[i] = -hess[i];
  for (int j = 0; j < 3 && !out.boundary; ++j) {
    const double h = 1e-5 * (1 + std::fabs(th[j]));
    double t[3] = {th[0], th[1], th[2]}, gj[3], hj[9];
    t[j] += h;
    std::copy(modes.begin(), modes.end(), trial_modes.begin());
    if (!std::isfinite(evaluate(data, t[0], t[1], t[2], r, trial_modes, gj, hj, poll)))
      return fail(Failure::Information);
    for (int i = 0; i < 3; ++i) info[3 * i + j] = -(gj[i] - grad[i]) / h;
  }
  for (int i = 0; i < 3; ++i)  // symmetrize
    for (int j = 0; j < i; ++j) info[3 * i + j] = info[3 * j + i] = 0.5 * (info[3 * i + j] + info[3 * j + i]);
  double e1[3] = {0, 1, 0}, col[3];
  const bool full = !out.boundary && solve_spd(info, e1, col, 3);
  if (!full && !solve_spd(info, e1, col, 2)) return fail(Failure::Information);
  out.se_b1 = std::sqrt(std::max(col[1], 0.0));
  return out;
}

double marginal_mean(double beta, double sigma) {
  const Rule& r = rule(40);
  double s = 0;
  for (std::size_t k = 0; k < r.x.size(); ++k) s += r.w[k] * expit(beta + sigma * SQRT2 * r.x[k]);
  return s / std::exp(LOG_SQRT_PI);
}

double between_variance(double beta, double sigma) {
  if (sigma == 0) return 0;  // which the quadrature gives only up to rounding
  const Rule& r = rule(40);
  double s = 0, s2 = 0;
  for (std::size_t k = 0; k < r.x.size(); ++k) {
    const double p = expit(beta + sigma * SQRT2 * r.x[k]);
    s += r.w[k] * p;
    s2 += r.w[k] * p * p;
  }
  const double norm = std::exp(LOG_SQRT_PI);
  return std::max(0.0, s2 / norm - (s / norm) * (s / norm));
}

double beta_for_marginal(double p, double sigma) {
  // [−60, 60] holds β for every p below 1 down to 8.8e-27 at σ = 0 and 5.7e-19 at σ = 6; smaller ones widen it
  double lo = -60, hi = 60;
  while (marginal_mean(lo, sigma) > p && lo > -1e3) lo *= 2;
  for (int i = 0; i < 200 && hi - lo > 1e-13; ++i) {
    const double mid = 0.5 * (lo + hi);
    (marginal_mean(mid, sigma) < p ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

}  // namespace ssc::sim::glmm
