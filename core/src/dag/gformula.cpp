#include "gformula.hpp"
#include "sim/random.hpp"
#include "statistics.hpp"
#include "stats/stats.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace ssc::dag {
namespace {
// Common random numbers for the chain's draws, independent of the study data.
constexpr std::uint64_t DRAWS = 0x47464f524d554c41ULL;
// Rows times draws per row: enough draws to keep the Monte Carlo error of the
// standardized means small next to the sampling error.
constexpr int DRAWS_PER_STUDY = 1000, MAX_DRAWS_PER_ROW = 20;

struct Fitted {
  int node = 0, p = 0;
  std::string family;
  std::vector<Term> terms;
  Vec means, scales, theta, covariance;
  double shape = 1;  // Gamma shape, held at its estimate
};
bool survival(const std::string &family) {
  return family == "weibull" || family == "exponential";
}
bool scaled(const std::string &family) {
  return family == "gaussian" || family == "uniform" || family == "lognormal";
}
double expit(double eta) {
  return eta >= 0 ? 1 / (1 + std::exp(-eta)) : std::exp(eta) / (1 + std::exp(eta));
}
double linear(const Fitted &f, const Vec &theta, const double *row) {
  double eta = theta[0];
  for (std::size_t j = 0; j < f.terms.size(); ++j) {
    const auto &t = f.terms[j];
    const double x = row[t.left] * (t.right < 0 ? 1 : row[t.right]);
    eta += theta[j + 1] * (x - f.means[j + 1]) / f.scales[j + 1];
  }
  return eta;
}
// The extra parameter after the coefficients: log SD, or log Weibull shape.
double extra(const Fitted &f, const Vec &theta) {
  return f.family == "exponential" ? 0 : theta[f.p];
}
Fitted fit(const Graph &g, const Data &d, const ChainModel &c) {
  const auto &node = g.nodes[c.node];
  Fitted f;
  f.node = c.node;
  f.family = node.family;
  f.terms = c.terms;
  const auto m = model_matrix(d, c.terms, true);
  f.p = m.p;
  f.means = m.means;
  f.scales = m.scales;
  if (survival(node.family)) {
    Vec times = d.times[c.node], events = d.events[c.node];
    if (times.empty()) {  // data without an observation rule: every time seen
      times = d.values[c.node];
      events.assign(d.n, 1);
    }
    const auto model = fit_weibull(m, times, events, node.family == "exponential");
    f.theta = model.beta;
    f.covariance = model.covariance;
    return f;
  }
  if (scaled(node.family)) {
    Vec y = d.values[c.node];
    if (node.family == "lognormal")
      for (double &v : y)
        v = std::log(v);
    const auto model = fit_model(m, y, "gaussian");
    // log SD appended, with the large-sample variance 1 / (2 (n - p))
    f.theta = model.beta;
    f.theta.push_back(.5 * std::log(model.dispersion));
    const int q = m.p + 1;
    f.covariance.assign(static_cast<std::size_t>(q) * q, 0);
    for (int i = 0; i < m.p; ++i)
      for (int j = 0; j < m.p; ++j)
        f.covariance[i * q + j] = model.covariance[i * m.p + j];
    f.covariance[q * q - 1] = 1 / (2.0 * (m.n - m.p));
    return f;
  }
  const auto model = fit_model(m, d.values[c.node], node.family);
  f.theta = model.beta;
  f.covariance = model.covariance;
  if (node.family == "gamma")
    f.shape = 1 / model.dispersion;
  return f;
}
// Standardized noise for one draw, fixed for the study: the parameters only
// shift and scale it, so the g-formula is smooth in them.
double noise(const Fitted &f, double u) {
  if (f.family == "gaussian" || f.family == "lognormal")
    return stats::qnorm(u);
  if (f.family == "gamma")
    return stats::qgamma(u, f.shape) / f.shape;
  if (survival(f.family))
    return -std::log1p(-u);
  return 0;
}
double draw(const Fitted &f, const Vec &theta, const double *row, double u,
            double z) {
  const double eta = linear(f, theta, row);
  if (f.family == "gaussian")
    return eta + std::exp(extra(f, theta)) * z;
  if (f.family == "uniform")
    return eta + std::exp(extra(f, theta)) * std::sqrt(12.0) * (u - .5);
  if (f.family == "lognormal")
    return std::exp(eta + std::exp(extra(f, theta)) * z);
  if (f.family == "gamma")
    return z * std::exp(eta);
  if (f.family == "bernoulli")
    return u < expit(eta) ? 1 : 0;
  if (f.family == "poisson")
    return stats::qpois(u, std::exp(eta));
  return std::exp((std::log(z) - eta) / std::exp(extra(f, theta)));
}
// The outcome's mean, P(Y = 1), S(horizon) or RMST to the horizon given its
// parents.
double functional(const Fitted &f, const Vec &theta, const double *row,
                  const std::string &measure, double horizon) {
  const double eta = linear(f, theta, row);
  if (survival(f.family)) {
    const double k = std::exp(extra(f, theta)),
                 x = std::exp(k * std::log(horizon) + eta);
    if (measure == "survival_difference")
      return std::exp(-x);
    // RMST = integral of exp(-(t / lambda)^k) to the horizon
    //      = lambda Gamma(1 + 1/k) P(1/k, (horizon / lambda)^k)
    return std::exp(-eta / k + stats::lgammafn(1 + 1 / k)) *
           stats::pgamma_raw(x, 1 / k, true);
  }
  if (f.family == "bernoulli")
    return expit(eta);
  if (f.family == "gaussian" || f.family == "uniform")
    return eta;
  if (f.family == "lognormal") {
    const double sd = std::exp(extra(f, theta));
    return std::exp(eta + sd * sd / 2);
  }
  return std::exp(eta);  // Gamma and Poisson means
}
bool ratio(const std::string &measure) {
  return measure == "risk_ratio" || measure == "odds_ratio";
}
// Lower Cholesky factor of a covariance; a direction without variance is left
// at zero.
Vec root(const Vec &a, int q) {
  Vec l(static_cast<std::size_t>(q) * q);
  for (int i = 0; i < q; ++i)
    for (int j = 0; j <= i; ++j) {
      double v = a[i * q + j];
      for (int k = 0; k < j; ++k)
        v -= l[i * q + k] * l[j * q + k];
      if (i == j)
        l[i * q + i] = v > 0 ? std::sqrt(v) : 0;
      else
        l[i * q + j] = l[j * q + j] > 0 ? v / l[j * q + j] : 0;
    }
  return l;
}
} // namespace

GFormula gformula(const Graph &g, const Data &d,
                  const std::vector<ChainModel> &chain, int exposure,
                  double reference, double active, const std::string &measure,
                  double horizon) {
  const int n = d.n, width = static_cast<int>(g.nodes.size()),
            links = static_cast<int>(chain.size());
  std::vector<Fitted> fits;
  for (const auto &c : chain)
    fits.push_back(fit(g, d, c));
  const int draws = links > 1 ? std::clamp((DRAWS_PER_STUDY + n - 1) / n, 1,
                                           MAX_DRAWS_PER_ROW)
                              : 1,
            cells = n * draws;
  std::vector<Vec> uniforms(links - 1), noises(links - 1);
  for (int j = 0; j + 1 < links; ++j) {
    uniforms[j].resize(cells);
    noises[j].resize(cells);
    for (int c = 0; c < cells; ++c) {
      sim::Rng rng(DRAWS, static_cast<std::uint64_t>(fits[j].node),
                   static_cast<std::uint64_t>(c));
      uniforms[j][c] = rng.uniform();
      noises[j][c] = noise(fits[j], uniforms[j][c]);
    }
  }
  // Every cell's variables under each level: observed values, the comparison
  // variable set, and the chain drawn in order. Kept so that a change in one
  // model redraws only that model and the ones after it.
  std::vector<Vec> base(2, Vec(static_cast<std::size_t>(cells) * width));
  for (int arm = 0; arm < 2; ++arm)
    for (int c = 0; c < cells; ++c) {
      double *row = &base[arm][static_cast<std::size_t>(c) * width];
      for (int v = 0; v < width; ++v)
        row[v] = d.values[v][c / draws];
      row[exposure] = arm ? active : reference;
    }
  std::vector<Vec> thetas;
  for (const auto &f : fits)
    thetas.push_back(f.theta);
  // Row means of the outcome functional for both levels; models before `from`
  // keep the draws stored in `base`, which the unperturbed pass fills in.
  const auto evaluate = [&](const std::vector<Vec> &th, int from, bool store,
                            Vec rows[2]) {
    Vec scratch(width);
    double mean[2] = {0, 0};
    for (int arm = 0; arm < 2; ++arm) {
      if (rows)
        rows[arm].assign(n, 0);
      for (int c = 0; c < cells; ++c) {
        double *row = &base[arm][static_cast<std::size_t>(c) * width];
        double *use = store ? row : scratch.data();
        if (!store)
          std::copy(row, row + width, use);
        for (int j = from; j + 1 < links; ++j) {
          const double v = draw(fits[j], th[j], use, uniforms[j][c], noises[j][c]);
          if (!std::isfinite(v))
            throw Error("numerical", "A g-formula draw overflowed.");
          use[fits[j].node] = v;
        }
        const double y = functional(fits.back(), th.back(), use, measure, horizon);
        if (!std::isfinite(y))
          throw Error("numerical", "A g-formula outcome mean overflowed.");
        mean[arm] += y;
        if (rows)
          rows[arm][c / draws] += y / draws;
      }
      mean[arm] /= cells;
    }
    return std::array<double, 2>{mean[0], mean[1]};
  };
  const auto scale = [&](const std::array<double, 2> &m) {
    if (!ratio(measure))
      return m[1] - m[0];
    if (!(m[0] > 0 && m[1] > 0) ||
        (measure == "odds_ratio" && !(m[0] < 1 && m[1] < 1)))
      throw Error("boundary_response",
                  "A standardized risk lies on the boundary of a probability.");
    return measure == "risk_ratio"
               ? std::log(m[1] / m[0])
               : std::log(m[1] / (1 - m[1])) - std::log(m[0] / (1 - m[0]));
  };
  Vec rows[2];
  const auto psi = evaluate(thetas, 0, true, rows);
  const double contrast = scale(psi);
  // Rows: the contrast's influence through the covariate distribution.
  Vec influence(n);
  for (int i = 0; i < n; ++i)
    influence[i] = measure == "risk_ratio"
                       ? rows[1][i] / psi[1] - rows[0][i] / psi[0]
                   : measure == "odds_ratio"
                       ? rows[1][i] / (psi[1] * (1 - psi[1])) -
                             rows[0][i] / (psi[0] * (1 - psi[0]))
                       : rows[1][i] - rows[0][i];
  double variance = moments(influence).variance / n;
  // Fitted parameters: each model's covariance, by symmetric differences along
  // the columns of its Cholesky factor.
  for (int j = 0; j < links; ++j) {
    const int q = static_cast<int>(fits[j].theta.size());
    const auto l = root(fits[j].covariance, q);
    for (int col = 0; col < q; ++col) {
      auto up = thetas, down = thetas;
      bool moved = false;
      for (int i = 0; i < q; ++i) {
        up[j][i] += l[i * q + col];
        down[j][i] -= l[i * q + col];
        moved |= l[i * q + col] != 0;
      }
      if (!moved)
        continue;
      const double half = (scale(evaluate(up, j, false, nullptr)) -
                           scale(evaluate(down, j, false, nullptr))) / 2;
      variance += half * half;
    }
  }
  GFormula out;
  out.testEstimate = contrast;
  out.estimate = ratio(measure) ? std::exp(contrast) : contrast;
  out.se = std::sqrt(variance);
  return out;
}
} // namespace ssc::dag
