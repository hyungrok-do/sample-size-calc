#include "learners.hpp"
#include "statistics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
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
  Booster(const Data &d, const std::vector<int> &variables, const Vec &y, bool logistic, const std::vector<int> *groups)
      : n_(d.n), p_(static_cast<int>(variables.size())), y_(y), logistic_(logistic), groups_(groups), bins_(p_), count_(p_) {
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
  // row, or every fifth cluster when rows are clustered (rows of a held-out
  // cluster seen in the fit would reward memorizing clusters), then the model
  // is refitted on all of them.
  Vec operator()(const std::vector<int> &rows) const {
    std::vector<int> fit, held;
    if (groups_) {
      std::map<int, int> order;
      for (int i : rows) order.emplace((*groups_)[i], static_cast<int>(order.size()));
      const int count = static_cast<int>(order.size()), every = std::min(5, count);
      for (int i : rows) (count >= 2 && order[(*groups_)[i]] % every == every - 1 ? held : fit).push_back(i);
    } else
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
        if (now < bestLoss * (1 - 1e-12)) { bestLoss = now; *bestRounds = round + 1; }
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
    double G = 0, H = 0, scale = 0;
    for (int i : rows) { G += g[i]; H += h[i]; scale += g[i] * g[i] / h[i]; }
    const double parent = G * G / (H + LAMBDA);
    // a split must gain a relative amount of the loss: outcomes on any scale split alike
    double bestGain = 1e-12 * scale;
    if (!(bestGain > 0)) return best;
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
  const std::vector<int> *groups_;
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

// The convex combination of the learners' cross-validated predictions with
// the smallest squared error: min ||sum_l w_l (cv_l - y)||^2 over weights
// w >= 0 summing to one. Working with the residual columns makes the weights
// independent of the outcome's location (an outcome near 1000 with SD 1 gets
// the same weights as one near 0), unlike non-negative least squares on the
// raw predictions followed by normalization. The minimum lies in the relative
// interior of some face of the simplex, where it solves the face's equality-
// constrained problem, so every support is tried (at most 63 for six learners)
// and the best feasible solution kept; a single learner is always feasible.
Vec simplex_weights(const std::vector<Vec> &columns, const Vec &y, const std::vector<int> &rows, const std::vector<bool> &use) {
  const int k = static_cast<int>(columns.size());
  std::vector<double> gram(static_cast<std::size_t>(k) * k);
  for (int a = 0; a < k; ++a) if (use[a])
    for (int b = a; b < k; ++b) if (use[b]) {
      double s = 0;
      for (int i : rows) s += (columns[a][i] - y[i]) * (columns[b][i] - y[i]);
      gram[a * k + b] = gram[b * k + a] = s;
    }
  Vec best(k);
  double bestRisk = INFINITY;
  for (int mask = 1; mask < (1 << k); ++mask) {
    std::vector<int> s;
    bool usable = true;
    for (int l = 0; l < k; ++l) if (mask >> l & 1) { if (!use[l]) usable = false; s.push_back(l); }
    if (!usable) continue;
    const int m = static_cast<int>(s.size());
    if (m == 1) {
      if (gram[s[0] * k + s[0]] < bestRisk) { bestRisk = gram[s[0] * k + s[0]]; best.assign(k, 0); best[s[0]] = 1; }
      continue;
    }
    // solve G_SS v = 1 by Gaussian elimination with partial pivoting; a ridge
    // far below the scale of G separates learners with identical predictions
    std::vector<double> a(static_cast<std::size_t>(m) * (m + 1));
    double trace = 0;
    for (int r = 0; r < m; ++r) trace += gram[s[r] * k + s[r]];
    const double ridge = 1e-12 * std::max(trace / m, 1e-300);
    for (int r = 0; r < m; ++r) {
      for (int c = 0; c < m; ++c) a[r * (m + 1) + c] = gram[s[r] * k + s[c]] + (r == c ? ridge : 0);
      a[r * (m + 1) + m] = 1;
    }
    bool singular = false;
    for (int c = 0; c < m && !singular; ++c) {
      int pivot = c;
      for (int r = c + 1; r < m; ++r) if (std::fabs(a[r * (m + 1) + c]) > std::fabs(a[pivot * (m + 1) + c])) pivot = r;
      if (!(std::fabs(a[pivot * (m + 1) + c]) > 0)) { singular = true; break; }
      if (pivot != c) for (int j = 0; j <= m; ++j) std::swap(a[c * (m + 1) + j], a[pivot * (m + 1) + j]);
      for (int r = c + 1; r < m; ++r) {
        const double f = a[r * (m + 1) + c] / a[c * (m + 1) + c];
        for (int j = c; j <= m; ++j) a[r * (m + 1) + j] -= f * a[c * (m + 1) + j];
      }
    }
    if (singular) continue;
    Vec v(m);
    for (int r = m - 1; r >= 0; --r) {
      double x = a[r * (m + 1) + m];
      for (int j = r + 1; j < m; ++j) x -= a[r * (m + 1) + j] * v[j];
      v[r] = x / a[r * (m + 1) + r];
    }
    double total = 0;
    for (double x : v) total += x;
    if (!(total > 0) || !std::isfinite(total)) continue;
    bool feasible = true;
    for (double &x : v) { x /= total; if (x < -1e-10) feasible = false; x = std::max(x, 0.0); }
    if (!feasible) continue;
    double risk = 0;
    for (int r = 0; r < m; ++r) for (int c = 0; c < m; ++c) risk += v[r] * v[c] * gram[s[r] * k + s[c]];
    if (risk < bestRisk) { bestRisk = risk; best.assign(k, 0); for (int r = 0; r < m; ++r) best[s[r]] = v[r]; }
  }
  return best;
}
} // namespace

Learned super_learner(const Data &d, const std::vector<int> &variables, const Vec &y,
                      const std::vector<int> &rows, const std::string &family, const std::vector<int> *groups) {
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
  auto booster = std::make_shared<Booster>(d, variables, y, logistic, groups);
  library[4] = [booster](const std::vector<int> &r) { return (*booster)(r); };
  // kernel ridge chooses its penalty by GCV, which treats rows as independent;
  // with clusters it can overfit them, and the cluster folds below judge it
  if (variables.empty()) use[5] = false;
  else {
    auto kernel = std::make_shared<KernelRidge>(d, variables, y, logistic);
    library[5] = [kernel](const std::vector<int> &r) { return (*kernel)(r); };
  }
  // folds: by position, separately for each outcome value of a binary
  // response; whole clusters when the rows are clustered
  std::vector<int> fold(n, -1);
  if (groups) {
    std::map<int, int> clusterFold;
    for (int i : rows) {
      const auto found = clusterFold.find((*groups)[i]);
      fold[i] = found != clusterFold.end() ? found->second : (clusterFold[(*groups)[i]] = static_cast<int>(clusterFold.size()) % FOLDS);
    }
    if (clusterFold.size() < FOLDS) throw Error("few_units", "Too few independent clusters for the Super Learner's cross-validation.");
  } else {
    int counter[2] = {0, 0};
    for (int i : rows) { const int c = logistic && y[i] > .5 ? 1 : 0; fold[i] = counter[c]++ % FOLDS; }
  }
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
  if (std::none_of(use.begin(), use.end(), [](bool u) { return u; }))
    throw Error("nonconvergence", "No Super Learner candidate could be fitted.");
  Vec w = simplex_weights(cv, y, rows, use);
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
