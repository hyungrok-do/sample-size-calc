// Fits use each replicate's observed data, never the DGP's true variances or
// coefficients.
#include "statistics.hpp"
#include "stats/stats.hpp"
#include <algorithm>
#include <functional>
#include <numeric>

namespace ssc::dag {
static void bad(const char *code, const char *message) {
  throw Error(code, message);
}
Moments moments(const Vec &v) {
  Moments out;
  double m2 = 0;
  for (double x : v) {
    const double d = x - out.mean;
    out.mean += d / ++out.n;
    m2 += d * (x - out.mean);
  }
  out.variance = out.n > 1 ? m2 / (out.n - 1) : 0;
  return out;
}
double p_value(double z, double df, const std::string &side) {
  if (std::isnan(z) || std::isnan(df) || df <= 0)
    bad("numerical", "Invalid test statistic or degrees of freedom.");
  // Both tails are evaluated directly by the calculator core, avoiding
  // cancellation for small p-values.
  if (side == "greater")
    return std::isfinite(df) ? stats::pt(z, df, false) : stats::pnorm(z, false);
  if (side == "less")
    return std::isfinite(df) ? stats::pt(z, df, true) : stats::pnorm(z, true);
  return std::min(1.0,
                  2 * (std::isfinite(df) ? stats::pt(std::fabs(z), df, false)
                                         : stats::pnorm(std::fabs(z), false)));
}
static Vec cholesky(const Vec &a, int p) {
  Vec l(p * p);
  double largest = 0;
  for (int i = 0; i < p; ++i)
    largest = std::max(largest, std::fabs(a[i * p + i]));
  if (!(largest > 0))
    bad("singular_model", "The analysis information matrix is singular.");
  for (int i = 0; i < p; ++i)
    for (int j = 0; j <= i; ++j) {
      double v = a[i * p + j];
      for (int k = 0; k < j; ++k)
        v -= l[i * p + k] * l[j * p + k];
      if (i == j) {
        if (!std::isfinite(v) || v <= largest * 1e-11)
          bad("singular_model",
              "The analysis is singular or numerically unidentified.");
        l[i * p + j] = std::sqrt(v);
      } else
        l[i * p + j] = v / l[j * p + j];
    }
  return l;
}
static Vec solve(const Vec &l, Vec out, int p) {
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < i; ++j)
      out[i] -= l[i * p + j] * out[j];
    out[i] /= l[i * p + i];
  }
  for (int i = p - 1; i >= 0; --i) {
    for (int j = i + 1; j < p; ++j)
      out[i] -= l[j * p + i] * out[j];
    out[i] /= l[i * p + i];
  }
  return out;
}
static double covariance_diagonal(const Vec &l, int p, int target) {
  Vec unit(p);
  unit.at(target) = 1;
  return solve(l, unit, p).at(target);
}
Matrix model_matrix(const Data &d, const std::vector<Term> &terms,
                    bool intercept) {
  Matrix m;
  m.n = d.n;
  const int offset = intercept ? 1 : 0;
  m.p = static_cast<int>(terms.size()) + offset;
  if (m.n <= m.p + 1)
    bad("too_few_observations", "Too few observations for the planned model.");
  m.x.resize(static_cast<std::size_t>(m.n) * m.p);
  m.scales.resize(m.p, 1);
  if (intercept)
    for (int i = 0; i < m.n; ++i)
      m.x[i * m.p] = 1;
  for (int j = 0; j < static_cast<int>(terms.size()); ++j) {
    const auto &term = terms[j];
    Vec values(m.n);
    for (int i = 0; i < m.n; ++i) {
      values[i] = d.values[term.left][i] *
                  (term.right < 0 ? 1 : d.values[term.right][i]);
      if (!std::isfinite(values[i]))
        bad("numerical", "A model term overflowed.");
    }
    const auto s = moments(values);
    const double scale = std::sqrt(s.variance);
    if (!(scale > 0) || !std::isfinite(scale))
      bad("constant_predictor", "A planned predictor has no variation or "
                                "finite scale in this replicate.");
    m.scales[j + offset] = scale;
    for (int i = 0; i < m.n; ++i)
      m.x[i * m.p + j + offset] = (values[i] - s.mean) / scale;
  }
  return m;
}
static double predictor(const Matrix &m, const Vec &beta, int i) {
  double x = 0;
  for (int j = 0; j < m.p; ++j)
    x += m.x[i * m.p + j] * beta[j];
  return x;
}
static void symmetrize(Vec &a, int p) {
  for (int j = 0; j < p; ++j)
    for (int k = 0; k < j; ++k)
      a[k * p + j] = a[j * p + k];
}
static Fit finish_fit(const Matrix &m, const Vec &beta, const Vec &info,
                      int target, double dispersion = 1, double df = INFINITY) {
  const auto factor = cholesky(info, m.p);
  const double se =
      std::sqrt(dispersion * covariance_diagonal(factor, m.p, target)) /
      m.scales[target];
  if (!(se > 0) || !std::isfinite(se))
    bad("zero_variance", "No finite positive regression standard error.");
  return {beta[target] / m.scales[target], se, df};
}
static Vec least_squares(const Matrix &m, const Vec &y, Vec &info) {
  info.assign(m.p * m.p, 0);
  Vec score(m.p);
  for (int i = 0; i < m.n; ++i)
    for (int j = 0; j < m.p; ++j) {
      const double x = m.x[i * m.p + j];
      score[j] += x * y[i];
      for (int k = 0; k <= j; ++k)
        info[j * m.p + k] += x * m.x[i * m.p + k];
    }
  symmetrize(info, m.p);
  return solve(cholesky(info, m.p), score, m.p);
}
Fit fit_linear(const Matrix &m, const Vec &y, int target) {
  Vec info;
  const auto beta = least_squares(m, y, info);
  double sse = 0;
  for (int i = 0; i < m.n; ++i) {
    const double r = y[i] - predictor(m, beta, i);
    sse += r * r;
  }
  return finish_fit(m, beta, info, target, sse / (m.n - m.p), m.n - m.p);
}
struct Evaluation {
  double ll = 0;
  Vec score, info, beta;
};
static double max_abs(const Vec &v) {
  double out = 0;
  for (double x : v)
    out = std::max(out, std::fabs(x));
  return out;
}
static Evaluation
maximize(Vec beta, const std::function<Evaluation(const Vec &, bool)> &evaluate,
         int size) {
  for (int iteration = 0; iteration < 60; ++iteration) {
    auto current = evaluate(beta, true);
    const int p = static_cast<int>(beta.size());
    if (!std::isfinite(current.ll))
      bad("numerical", "The model likelihood is not finite.");
    const auto step = solve(cholesky(current.info, p), current.score, p);
    if (max_abs(step) < 1e-8 * (1 + max_abs(beta)) &&
        max_abs(current.score) / size < 1e-7) {
      current.beta = beta;
      return current;
    }
    bool accepted = false;
    for (double scale = 1; scale >= 1.0 / 1048576; scale /= 2) {
      Vec trial = beta;
      for (int j = 0; j < p; ++j)
        trial[j] += scale * step[j];
      const double ll = evaluate(trial, false).ll;
      if (std::isfinite(ll) &&
          ll >= current.ll - 1e-11 * (1 + std::fabs(current.ll))) {
        beta = std::move(trial);
        accepted = true;
        break;
      }
    }
    if (!accepted)
      bad("nonconvergence", "The model line search failed.");
  }
  bad("nonconvergence", "The analysis did not converge; separation or an "
                        "infinite coefficient may be present.");
  return {};
}
static Evaluation maximize_glm(const Matrix &m, const Vec &y,
                               const std::string &family) {
  const double average = moments(y).mean;
  if (!(average > 0) || (family == "bernoulli" && average >= 1))
    bad("boundary_response", "All responses lie on a distribution boundary.");
  const auto evaluate = [&](const Vec &beta, bool derivatives) {
    Evaluation e;
    e.score.resize(m.p);
    e.info.resize(m.p * m.p);
    for (int i = 0; i < m.n; ++i) {
      const double eta = predictor(m, beta, i),
                   mu = family == "bernoulli"
                            ? (eta >= 0 ? 1 / (1 + std::exp(-eta))
                                        : std::exp(eta) / (1 + std::exp(eta)))
                            : std::exp(eta);
      if (!std::isfinite(mu) || (mu <= 0 && family != "bernoulli")) {
        e.ll = -INFINITY;
        return e;
      }
      e.ll += family == "bernoulli" ? y[i] * eta - std::max(0.0, eta) -
                                          std::log1p(std::exp(-std::fabs(eta)))
              : family == "poisson" ? y[i] * eta - mu
                                    : -y[i] / mu - eta;
      if (!derivatives)
        continue;
      const double residual = family == "gamma" ? y[i] / mu - 1 : y[i] - mu,
                   weight = family == "bernoulli" ? mu * (1 - mu)
                            : family == "poisson" ? mu
                                                  : 1;
      for (int j = 0; j < m.p; ++j) {
        const double x = m.x[i * m.p + j];
        e.score[j] += x * residual;
        for (int k = 0; k <= j; ++k)
          e.info[j * m.p + k] += weight * x * m.x[i * m.p + k];
      }
    }
    symmetrize(e.info, m.p);
    return e;
  };
  Vec initial(m.p);
  initial[0] = family == "bernoulli" ? std::log(average / (1 - average))
                                     : std::log(average);
  return maximize(initial, evaluate, m.n);
}
Fit fit_glm(const Matrix &m, const Vec &y, int target,
            const std::string &family) {
  const auto fit = maximize_glm(m, y, family);
  double dispersion = 1;
  if (family == "gamma") {
    double pearson = 0;
    for (int i = 0; i < m.n; ++i) {
      const double r = y[i] / std::exp(predictor(m, fit.beta, i)) - 1;
      pearson += r * r;
    }
    dispersion = pearson / (m.n - m.p);
  }
  return finish_fit(m, fit.beta, fit.info, target, dispersion);
}
Vec fitted_means(const Matrix &m, const Vec &y, const std::vector<int> &rows,
                 const std::string &family) {
  Matrix part;
  part.n = static_cast<int>(rows.size());
  part.p = m.p;
  part.scales = m.scales;
  if (part.n <= part.p + 1)
    bad("too_few_observations",
        "Too few observations in an arm for the outcome model.");
  part.x.reserve(static_cast<std::size_t>(part.n) * part.p);
  Vec response;
  response.reserve(part.n);
  for (int i : rows) {
    part.x.insert(part.x.end(), m.x.begin() + i * m.p,
                  m.x.begin() + (i + 1) * m.p);
    response.push_back(y[i]);
  }
  Vec info;
  const auto beta = family == "gaussian"
                        ? least_squares(part, response, info)
                        : maximize_glm(part, response, family).beta;
  Vec out(m.n);
  for (int i = 0; i < m.n; ++i) {
    const double eta = predictor(m, beta, i);
    out[i] = family == "gaussian"    ? eta
             : family == "bernoulli" ? (eta >= 0 ? 1 / (1 + std::exp(-eta))
                                                 : std::exp(eta) /
                                                       (1 + std::exp(eta)))
                                     : std::exp(eta);
    if (!std::isfinite(out[i]))
      bad("numerical", "A fitted mean overflowed.");
  }
  return out;
}
Fit fit_cox(const Matrix &m, const Vec &times, const Vec &events, int target) {
  const int count =
      static_cast<int>(std::accumulate(events.begin(), events.end(), 0.0));
  if (count <= m.p + 1)
    bad("too_few_events", "Too few observed events for the planned Cox model.");
  std::vector<int> order(m.n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
                   [&](int a, int b) { return times[a] > times[b]; });
  const auto evaluate = [&](const Vec &beta, bool derivatives) {
    Evaluation e;
    e.score.resize(m.p);
    e.info.resize(m.p * m.p);
    Vec eta(m.n), s1(m.p), s2(m.p * m.p);
    for (int i = 0; i < m.n; ++i)
      eta[i] = predictor(m, beta, i);
    double shift = -INFINITY, s0 = 0;
    for (int at = 0; at < m.n;) {
      int end = at;
      while (end < m.n && times[order[end]] == times[order[at]])
        ++end;
      int deaths = 0;
      double eventEta = 0;
      Vec eventX(m.p);
      for (int k = at; k < end; ++k) {
        const int i = order[k];
        if (eta[i] > shift) {
          const double f = std::exp(shift - eta[i]);
          s0 *= f;
          for (double &v : s1)
            v *= f;
          for (double &v : s2)
            v *= f;
          shift = eta[i];
        }
        const double w = std::exp(eta[i] - shift);
        s0 += w;
        if (derivatives)
          for (int j = 0; j < m.p; ++j) {
            const double x = m.x[i * m.p + j];
            s1[j] += w * x;
            for (int l = 0; l <= j; ++l)
              s2[j * m.p + l] += w * x * m.x[i * m.p + l];
          }
        if (events[i]) {
          ++deaths;
          eventEta += eta[i];
          if (derivatives)
            for (int j = 0; j < m.p; ++j)
              eventX[j] += m.x[i * m.p + j];
        }
      }
      if (deaths) {
        if (!(s0 > 0)) {
          e.ll = -INFINITY;
          return e;
        }
        e.ll += eventEta - deaths * (shift + std::log(s0));
        if (derivatives)
          for (int j = 0; j < m.p; ++j) {
            e.score[j] += eventX[j] - deaths * s1[j] / s0;
            for (int l = 0; l <= j; ++l)
              e.info[j * m.p + l] +=
                  deaths * (s2[j * m.p + l] / s0 - s1[j] * s1[l] / (s0 * s0));
          }
      }
      at = end;
    }
    symmetrize(e.info, m.p);
    return e;
  };
  const auto fit = maximize(Vec(m.p), evaluate, count);
  return finish_fit(m, fit.beta, fit.info, target);
}
KM kaplan_meier(const Vec &times, const Vec &events, double horizon) {
  const int n = static_cast<int>(times.size());
  if (n < 2)
    bad("empty_arm", "An analysis arm has fewer than two observations.");
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
                   [&](int a, int b) { return times[a] < times[b]; });
  KM k;
  double risk = n, previous = 0, greenwood = 0;
  std::vector<std::pair<double, double>> increments;
  for (int at = 0; at < n && times[order[at]] <= horizon;) {
    const double time = times[order[at]];
    int end = at;
    double deaths = 0;
    while (end < n && times[order[end]] == time) {
      deaths += events[order[end]];
      ++end;
    }
    k.rmst += k.survival * (time - previous);
    previous = time;
    if (deaths) {
      const double q = risk > deaths ? deaths / (risk * (risk - deaths)) : 0;
      increments.push_back({k.rmst, q});
      greenwood += q;
      k.survival *= 1 - deaths / risk;
    }
    risk -= end - at;
    at = end;
  }
  if (times[order.back()] < horizon && k.survival > 0)
    bad("insufficient_followup", "Follow-up ends before the estimand horizon "
                                 "while the KM curve is positive.");
  k.rmst += k.survival * (horizon - previous);
  k.survivalVariance = k.survival * k.survival * greenwood;
  for (const auto &t : increments) {
    const double tail = k.rmst - t.first;
    k.rmstVariance += tail * tail * t.second;
  }
  return k;
}
} // namespace ssc::dag
