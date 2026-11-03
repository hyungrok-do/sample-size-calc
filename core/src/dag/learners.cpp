#include "learners.hpp"
#include "statistics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <set>

namespace ssc::dag {
const std::vector<std::string> LEARNERS = {"mean", "glm", "glm.interaction", "gam", "gbm", "kernel"};
namespace {
constexpr int FOLDS = 5;
constexpr double TINY = 1e-6;
double expit(double x) { return x >= 0 ? 1 / (1 + std::exp(-x)) : std::exp(x) / (1 + std::exp(x)); }
bool binary(const std::string &family) { return family == "bernoulli"; }

// A learner predicts every row from a fit to `rows`.
using Learner = std::function<Vec(const std::vector<int> &rows)>;

Vec mean_learner(const Vec &y, const std::vector<int> &rows, int n, bool logistic) {
  double total = 0;
  for (int i : rows) total += y[i];
  double m = total / std::max<std::size_t>(1, rows.size());
  if (logistic) m = std::clamp(m, TINY, 1 - TINY);
  return Vec(n, m);
}
// A GLM on fixed terms: its model matrix is built once (covariates only).
Learner glm_learner(const Data &d, const std::vector<Term> &terms, const Vec &y, bool logistic) {
  if (terms.empty()) return [&y, n = d.n, logistic](const std::vector<int> &rows) { return mean_learner(y, rows, n, logistic); };
  auto matrix = std::make_shared<Matrix>(model_matrix(d, terms, true));
  return [matrix, &y, logistic](const std::vector<int> &rows) {
    return fitted_means(*matrix, y, rows, logistic ? "bernoulli" : "gaussian");
  };
}
// Distinct values of a column, up to `limit` + 1.
std::size_t distinct(const Vec &values, std::size_t limit) {
  std::set<double> seen;
  for (double v : values) { seen.insert(v); if (seen.size() > limit) break; }
  return seen.size();
}
double quantile(const Vec &sorted, double q) {
  const double at = q * (sorted.size() - 1);
  const auto lo = static_cast<std::size_t>(std::floor(at));
  const auto hi = std::min(lo + 1, sorted.size() - 1);
  return sorted[lo] + (at - lo) * (sorted[hi] - sorted[lo]);
}

// Gradient boosting with depth-2 regression trees on binned covariates:
// squared error, or the logistic loss with Newton leaf values.
class Booster {
public:
  Booster(const Data &d, const std::vector<int> &variables, const Vec &y, bool logistic)
      : n_(d.n), p_(static_cast<int>(variables.size())), y_(y), logistic_(logistic), bins_(p_), count_(p_) {
    for (int j = 0; j < p_; ++j) {
      const Vec &x = d.values[variables[j]];
      Vec sorted(x);
      std::sort(sorted.begin(), sorted.end());
      Vec edges;
      for (int b = 1; b < BINS; ++b) {
        const double e = quantile(sorted, static_cast<double>(b) / BINS);
        if (edges.empty() || e > edges.back()) edges.push_back(e);
      }
      bins_[j].resize(n_);
      for (int i = 0; i < n_; ++i)
        bins_[j][i] = static_cast<std::uint8_t>(std::lower_bound(edges.begin(), edges.end(), x[i]) - edges.begin());
      count_[j] = static_cast<int>(edges.size()) + 1;
    }
  }
  // The number of rounds is chosen by early stopping on every fifth training
  // row, then the model is refitted on all of them.
  Vec operator()(const std::vector<int> &rows) const {
    std::vector<int> fit, held;
    for (std::size_t k = 0; k < rows.size(); ++k) (k % 5 == 4 ? held : fit).push_back(rows[k]);
    int best = ROUNDS;
    if (held.size() >= 10 && fit.size() >= 20) boost(fit, ROUNDS, &held, &best);
    return boost(rows, std::max(best, 1), nullptr, nullptr);
  }
private:
  Vec boost(const std::vector<int> &rows, int rounds, const std::vector<int> *held, int *bestRounds) const {
    double total = 0;
    for (int i : rows) total += y_[i];
    const double average = total / std::max<std::size_t>(1, rows.size());
    const double start = logistic_ ? std::log(std::clamp(average, 1e-3, 1 - 1e-3) / (1 - std::clamp(average, 1e-3, 1 - 1e-3))) : average;
    Vec f(n_, start), g(n_), h(n_);
    std::vector<int> left, right, ll, lr, rl, rr;
    double bestLoss = INFINITY;
    const auto loss = [&] {
      double sum = 0;
      for (int i : *held) {
        if (logistic_) { const double q = std::clamp(expit(f[i]), 1e-12, 1 - 1e-12); sum -= y_[i] * std::log(q) + (1 - y_[i]) * std::log1p(-q); }
        else sum += (y_[i] - f[i]) * (y_[i] - f[i]);
      }
      return sum;
    };
    if (held) { bestLoss = loss(); *bestRounds = 0; }
    for (int round = 0; round < rounds && p_; ++round) {
      for (int i : rows) {
        if (logistic_) { const double q = expit(f[i]); g[i] = y_[i] - q; h[i] = std::max(q * (1 - q), 1e-6); }
        else { g[i] = y_[i] - f[i]; h[i] = 1; }
      }
      const Split root = split(rows, g, h);
      if (root.feature < 0) break;
      partition(rows, root, left, right);
      const Split a = split(left, g, h), b = split(right, g, h);
      // leaves: the Newton step, shrunk; rows reached in the training data
      // decide each leaf's value, and every row of the data is routed
      double values[4];
      const auto leaf = [&](const std::vector<int> &part) {
        double gs = 0, hs = 0;
        for (int i : part) { gs += g[i]; hs += h[i]; }
        return RATE * gs / (hs + LAMBDA);
      };
      if (a.feature >= 0) { partition(left, a, ll, lr); values[0] = leaf(ll); values[1] = leaf(lr); }
      else values[0] = values[1] = leaf(left);
      if (b.feature >= 0) { partition(right, b, rl, rr); values[2] = leaf(rl); values[3] = leaf(rr); }
      else values[2] = values[3] = leaf(right);
      for (int i = 0; i < n_; ++i) {
        const bool goLeft = bins_[root.feature][i] <= root.threshold;
        const Split &next = goLeft ? a : b;
        const bool second = next.feature >= 0 && bins_[next.feature][i] > next.threshold;
        f[i] += values[(goLeft ? 0 : 2) + (second ? 1 : 0)];
      }
      if (held) {
        const double now = loss();
        if (now < bestLoss - 1e-12) { bestLoss = now; *bestRounds = round + 1; }
        else if (round + 1 - *bestRounds >= PATIENCE) break;
      }
    }
    Vec out(n_);
    for (int i = 0; i < n_; ++i) out[i] = logistic_ ? std::clamp(expit(f[i]), TINY, 1 - TINY) : f[i];
    return out;
  }
  static constexpr int BINS = 32, ROUNDS = 1000, MIN_LEAF = 5, PATIENCE = 50;
  static constexpr double RATE = .1, LAMBDA = 1;
  struct Split { int feature = -1, threshold = 0; };
  Split split(const std::vector<int> &rows, const Vec &g, const Vec &h) const {
    Split best;
    if (static_cast<int>(rows.size()) < 2 * MIN_LEAF) return best;
    double G = 0, H = 0;
    for (int i : rows) { G += g[i]; H += h[i]; }
    const double parent = G * G / (H + LAMBDA);
    double bestGain = 1e-12;
    double gb[BINS], hb[BINS];
    int nb[BINS];
    for (int j = 0; j < p_; ++j) {
      std::fill(gb, gb + count_[j], 0.0); std::fill(hb, hb + count_[j], 0.0); std::fill(nb, nb + count_[j], 0);
      for (int i : rows) { const int b = bins_[j][i]; gb[b] += g[i]; hb[b] += h[i]; ++nb[b]; }
      double gl = 0, hl = 0; int nl = 0;
      for (int b = 0; b + 1 < count_[j]; ++b) {
        gl += gb[b]; hl += hb[b]; nl += nb[b];
        const int nr = static_cast<int>(rows.size()) - nl;
        if (nl < MIN_LEAF || nr < MIN_LEAF) continue;
        const double gain = gl * gl / (hl + LAMBDA) + (G - gl) * (G - gl) / (H - hl + LAMBDA) - parent;
        if (gain > bestGain) { bestGain = gain; best.feature = j; best.threshold = b; }
      }
    }
    return best;
  }
  void partition(const std::vector<int> &rows, const Split &s, std::vector<int> &left, std::vector<int> &right) const {
    left.clear(); right.clear();
    for (int i : rows) (bins_[s.feature][i] <= s.threshold ? left : right).push_back(i);
  }
  int n_, p_;
  const Vec &y_;
  bool logistic_;
  std::vector<std::vector<std::uint8_t>> bins_;
  std::vector<int> count_;
};

// Kernel ridge regression approximated by random Fourier features of a
// Gaussian kernel on the standardized covariates: bandwidth by the median
// heuristic, penalty by generalized cross-validation. A binary response is
// fitted as a probability and clipped.
class KernelRidge {
public:
  KernelRidge(const Data &d, const std::vector<int> &variables, const Vec &y, bool logistic)
      : n_(d.n), y_(y), logistic_(logistic), phi_(static_cast<std::size_t>(d.n) * FEATURES) {
    const int p = static_cast<int>(variables.size());
    std::vector<Vec> z(p, Vec(n_));
    for (int j = 0; j < p; ++j) {
      const Vec &x = d.values[variables[j]];
      double m = 0, m2 = 0;
      for (int i = 0; i < n_; ++i) { const double delta = x[i] - m; m += delta / (i + 1); m2 += delta * (x[i] - m); }
      const double sd = n_ > 1 ? std::sqrt(m2 / (n_ - 1)) : 1;
      for (int i = 0; i < n_; ++i) z[j][i] = (x[i] - m) / (sd > 0 ? sd : 1);
    }
    // median distance between pairs of the first rows
    Vec distances;
    const int sample = std::min(n_, 150);
    for (int a = 0; a < sample; ++a) for (int b = a + 1; b < sample; ++b) {
      double s = 0; for (int j = 0; j < p; ++j) s += (z[j][a] - z[j][b]) * (z[j][a] - z[j][b]);
      distances.push_back(std::sqrt(s));
    }
    double bandwidth = 1;
    if (!distances.empty()) {
      std::nth_element(distances.begin(), distances.begin() + distances.size() / 2, distances.end());
      bandwidth = distances[distances.size() / 2] > 0 ? distances[distances.size() / 2] : 1;
    }
    std::uint64_t state = 0x6b65726e656cULL;  // a fixed generator: the learner is deterministic
    const auto uniform = [&state] { state = state * 6364136223846793005ULL + 1442695040888963407ULL; return ((state >> 11) + .5) / 9007199254740992.0; };
    const auto normal = [&] { return std::sqrt(-2 * std::log(uniform())) * std::cos(6.283185307179586 * uniform()); };
    std::vector<double> omega(static_cast<std::size_t>(FEATURES) * p), phase(FEATURES);
    for (int m = 0; m < FEATURES; ++m) {
      for (int j = 0; j < p; ++j) omega[m * p + j] = normal() / bandwidth;
      phase[m] = 6.283185307179586 * uniform();
    }
    for (int i = 0; i < n_; ++i) for (int m = 0; m < FEATURES; ++m) {
      double dot = phase[m];
      for (int j = 0; j < p; ++j) dot += omega[m * p + j] * z[j][i];
      phi_[static_cast<std::size_t>(i) * FEATURES + m] = std::sqrt(2.0 / FEATURES) * std::cos(dot);
    }
  }
  Vec operator()(const std::vector<int> &rows) const {
    const int r = static_cast<int>(rows.size()), M = FEATURES;
    double ybar = 0;
    for (int i : rows) ybar += y_[i];
    ybar /= std::max(1, r);
    // centred features and response: an unpenalized intercept
    Vec mean(M);
    for (int i : rows) for (int m = 0; m < M; ++m) mean[m] += phi_[static_cast<std::size_t>(i) * M + m] / r;
    Vec gram(M * M), cross(M);
    double yy = 0;
    for (int i : rows) {
      const double *f = &phi_[static_cast<std::size_t>(i) * M];
      const double yc = y_[i] - ybar;
      yy += yc * yc;
      for (int a = 0; a < M; ++a) {
        const double fa = f[a] - mean[a];
        cross[a] += fa * yc;
        for (int b = 0; b <= a; ++b) gram[a * M + b] += fa * (f[b] - mean[b]);
      }
    }
    for (int a = 0; a < M; ++a) for (int b = 0; b < a; ++b) gram[b * M + a] = gram[a * M + b];
    Vec values, vectors;
    eigen(gram, M, values, vectors);
    Vec projected(M);
    for (int k = 0; k < M; ++k) for (int a = 0; a < M; ++a) projected[k] += vectors[a * M + k] * cross[a];
    double trace = 0;
    for (double v : values) trace += std::max(0.0, v);
    double bestScore = INFINITY, lambda = trace / M;
    for (int step = -8; step <= 4; ++step) {
      const double l = trace / M * std::pow(10.0, step / 2.0);
      if (!(l > 0)) continue;
      double rss = yy, df = 0;
      for (int k = 0; k < M; ++k) {
        const double d = std::max(0.0, values[k]), c2 = projected[k] * projected[k];
        rss += -2 * c2 / (d + l) + d * c2 / ((d + l) * (d + l));
        df += d / (d + l);
      }
      if (!(r - df > 1)) continue;
      const double score = r * std::max(rss, 0.0) / ((r - df) * (r - df));
      if (score < bestScore) { bestScore = score; lambda = l; }
    }
    Vec beta(M);
    for (int k = 0; k < M; ++k) {
      const double scaled = projected[k] / (std::max(0.0, values[k]) + lambda);
      for (int a = 0; a < M; ++a) beta[a] += vectors[a * M + k] * scaled;
    }
    double offset = ybar;
    for (int a = 0; a < M; ++a) offset -= beta[a] * mean[a];
    Vec out(n_);
    for (int i = 0; i < n_; ++i) {
      double v = offset;
      for (int a = 0; a < M; ++a) v += beta[a] * phi_[static_cast<std::size_t>(i) * M + a];
      out[i] = logistic_ ? std::clamp(v, .005, .995) : v;
    }
    return out;
  }
private:
  static constexpr int FEATURES = 48;
  // Symmetric eigendecomposition by cyclic Jacobi rotations; vectors[a * M + k]
  // holds component a of eigenvector k.
  static void eigen(Vec a, int M, Vec &values, Vec &vectors) {
    vectors.assign(M * M, 0);
    for (int k = 0; k < M; ++k) vectors[k * M + k] = 1;
    for (int sweep = 0; sweep < 30; ++sweep) {
      double off = 0;
      for (int p = 0; p < M; ++p) for (int q = p + 1; q < M; ++q) off += a[p * M + q] * a[p * M + q];
      if (off < 1e-22) break;
      for (int p = 0; p < M; ++p) for (int q = p + 1; q < M; ++q) {
        const double apq = a[p * M + q];
        if (std::fabs(apq) < 1e-300) continue;
        const double theta = (a[q * M + q] - a[p * M + p]) / (2 * apq);
        const double t = (theta >= 0 ? 1 : -1) / (std::fabs(theta) + std::sqrt(theta * theta + 1));
        const double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < M; ++k) {
          const double akp = a[k * M + p], akq = a[k * M + q];
          a[k * M + p] = c * akp - s * akq; a[k * M + q] = s * akp + c * akq;
        }
        for (int k = 0; k < M; ++k) {
          const double apk = a[p * M + k], aqk = a[q * M + k];
          a[p * M + k] = c * apk - s * aqk; a[q * M + k] = s * apk + c * aqk;
        }
        for (int k = 0; k < M; ++k) {
          const double vkp = vectors[k * M + p], vkq = vectors[k * M + q];
          vectors[k * M + p] = c * vkp - s * vkq; vectors[k * M + q] = s * vkp + c * vkq;
        }
      }
    }
    values.resize(M);
    for (int k = 0; k < M; ++k) values[k] = a[k * M + k];
  }
  int n_;
  const Vec &y_;
  bool logistic_;
  Vec phi_;
};

// Non-negative least squares by cyclic coordinate descent on the normal
// equations (a handful of learners).
Vec nnls(const std::vector<Vec> &columns, const Vec &y, const std::vector<int> &rows, const std::vector<bool> &use) {
  const int k = static_cast<int>(columns.size());
  std::vector<double> gram(k * k), cross(k);
  for (int a = 0; a < k; ++a) {
    if (!use[a]) continue;
    for (int i : rows) cross[a] += columns[a][i] * y[i];
    for (int b = 0; b < k; ++b) if (use[b]) for (int i : rows) gram[a * k + b] += columns[a][i] * columns[b][i];
  }
  Vec w(k);
  for (int sweep = 0; sweep < 1000; ++sweep) {
    double change = 0;
    for (int a = 0; a < k; ++a) {
      if (!use[a] || !(gram[a * k + a] > 0)) continue;
      double r = cross[a];
      for (int b = 0; b < k; ++b) r -= gram[a * k + b] * w[b];
      const double next = std::max(0.0, w[a] + r / gram[a * k + a]);
      change = std::max(change, std::fabs(next - w[a]));
      w[a] = next;
    }
    if (change < 1e-12) break;
  }
  return w;
}
} // namespace

Learned super_learner(const Data &d, const std::vector<int> &variables, const Vec &y,
                      const std::vector<int> &rows, const std::string &family) {
  const bool logistic = binary(family);
  const int n = d.n, k = static_cast<int>(LEARNERS.size());
  if (rows.size() < 2 * FOLDS) throw Error("too_few_observations", "Too few observations for the Super Learner's cross-validation.");
  // The learners' fixed parts depend on the covariates only.
  std::vector<Term> linear, interaction, additive;
  for (int v : variables) linear.push_back(Term{v});
  interaction = linear;
  for (std::size_t a = 0; a < variables.size(); ++a)
    for (std::size_t b = a + 1; b < variables.size() && interaction.size() < 40; ++b)
      interaction.push_back(Term{variables[a], variables[b]});
  additive = linear;
  for (int v : variables) {
    if (distinct(d.values[v], 4) <= 4) continue;  // binary and few-valued variables stay linear
    Vec sorted(d.values[v]);
    std::sort(sorted.begin(), sorted.end());
    Vec knots;
    for (double q : {.05, .35, .65, .95}) { const double x = quantile(sorted, q); if (knots.empty() || x > knots.back()) knots.push_back(x); }
    if (knots.size() < 3) continue;
    for (int index = 0; index + 2 < static_cast<int>(knots.size()); ++index) {
      Term t{v};
      t.transform.kind = "spline"; t.transform.knots = knots; t.transform.index = index;
      additive.push_back(t);
    }
  }
  std::vector<Learner> library(k);
  std::vector<bool> use(k, true);
  library[0] = [&y, n, logistic](const std::vector<int> &r) { return mean_learner(y, r, n, logistic); };
  const std::vector<Term> *bases[3] = {&linear, &interaction, &additive};
  for (int l = 1; l <= 3; ++l) {
    try { library[l] = glm_learner(d, *bases[l - 1], y, logistic); }
    catch (const Error &) { use[l] = false; }
  }
  if (variables.empty()) use[1] = use[2] = use[3] = use[4] = use[5] = false;
  else if (bases[1]->size() == linear.size()) use[2] = false;   // one variable: no products
  if (additive.size() == linear.size()) use[3] = false;          // nothing to smooth
  auto booster = std::make_shared<Booster>(d, variables, y, logistic);
  library[4] = [booster](const std::vector<int> &r) { return (*booster)(r); };
  if (variables.empty()) use[5] = false;
  else {
    auto kernel = std::make_shared<KernelRidge>(d, variables, y, logistic);
    library[5] = [kernel](const std::vector<int> &r) { return (*kernel)(r); };
  }
  // folds: by position, separately for each outcome value of a binary response
  std::vector<int> fold(n, -1);
  int counter[2] = {0, 0};
  for (int i : rows) { const int c = logistic && y[i] > .5 ? 1 : 0; fold[i] = counter[c]++ % FOLDS; }
  std::vector<Vec> cv(k, Vec(n));
  for (int v = 0; v < FOLDS; ++v) {
    std::vector<int> train;
    for (int i : rows) if (fold[i] != v) train.push_back(i);
    for (int l = 0; l < k; ++l) {
      if (!use[l]) continue;
      try {
        const Vec pred = library[l](train);
        for (int i : rows) if (fold[i] == v) {
          if (!std::isfinite(pred[i])) throw Error("numerical", "A learner's prediction overflowed.");
          cv[l][i] = pred[i];
        }
      } catch (const Error &) { use[l] = false; }
    }
  }
  Vec w = nnls(cv, y, rows, use);
  double total = 0;
  for (double x : w) total += x;
  if (!(total > 0)) {
    // no learner improves on zero: keep the one with the smallest CV risk
    int best = 0; double risk = INFINITY;
    for (int l = 0; l < k; ++l) if (use[l]) {
      double r = 0; for (int i : rows) r += (y[i] - cv[l][i]) * (y[i] - cv[l][i]);
      if (r < risk) { risk = r; best = l; }
    }
    w.assign(k, 0); w[best] = 1; total = 1;
  }
  for (double &x : w) x /= total;
  Learned out;
  out.prediction.assign(n, 0);
  for (int l = 0; l < k; ++l) {
    if (!(w[l] > 1e-10)) { w[l] = 0; continue; }
    Vec pred;
    try { pred = library[l](rows); }
    catch (const Error &) {
      // a learner that cross-validated but cannot be refitted on all rows
      w[l] = 0; continue;
    }
    for (int i = 0; i < n; ++i) out.prediction[i] += w[l] * pred[i];
  }
  double kept = 0;
  for (double x : w) kept += x;
  if (!(kept > 0)) throw Error("nonconvergence", "No Super Learner candidate could be fitted.");
  for (int i = 0; i < n; ++i) {
    out.prediction[i] /= kept;
    if (logistic) out.prediction[i] = std::clamp(out.prediction[i], TINY, 1 - TINY);
  }
  for (double &x : w) x /= kept;
  out.weights = w;
  return out;
}
} // namespace ssc::dag
