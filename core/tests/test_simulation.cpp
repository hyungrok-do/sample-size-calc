// Simulation procedures: random streams and variates, the REML fitter against
// closed forms and a brute-force likelihood, Satterthwaite's degrees of
// freedom, simulated power against exact power, type I error, reproducibility
// and the sample-size search.
#include <atomic>
#include <chrono>
#include <thread>
#include <vector>

#include "engine/runner.hpp"
#include "sim/lmm.hpp"
#include "sim/random.hpp"
#include "stats/stats.hpp"
#include "testing.hpp"
#include "util/json.hpp"
#include "util/numtext.hpp"

using namespace ssc;
namespace lmm = ssc::sim::lmm;

namespace {

json::Value run(const std::string& request, Token* token = nullptr) {
  Token local;
  return run_request(json::parse(request), token ? token : &local);
}

json::Value first_row(const std::string& request) {
  const json::Value res = run(request);
  CHECK(res.get("ok").as_bool(false));
  return res.get("rows").elements().at(0);
}

double num(const json::Value& row, const char* key) { return row.get(key).as_number(NAN); }

// Moments of n draws.
template <class F>
void moments(F&& draw, int n, double& mean, double& var) {
  double s = 0, s2 = 0;
  for (int i = 0; i < n; ++i) {
    const double x = draw();
    s += x;
    s2 += x * x;
  }
  mean = s / n;
  var = (s2 - s * s / n) / (n - 1);
}

// Dense REML criterion (−2 log-likelihood) of clusters (X_i, Z_i, y_i) at (σ², G), for checking lmm::criterion.
struct DenseCluster {
  int n;
  std::vector<double> x, z, y;
};

bool dense_chol(std::vector<double>& a, int n) {
  for (int j = 0; j < n; ++j) {
    double s = a[j * n + j];
    for (int k = 0; k < j; ++k) s -= a[j * n + k] * a[j * n + k];
    if (!(s > 0)) return false;
    a[j * n + j] = std::sqrt(s);
    for (int i = j + 1; i < n; ++i) {
      double t = a[i * n + j];
      for (int k = 0; k < j; ++k) t -= a[i * n + k] * a[j * n + k];
      a[i * n + j] = t / a[j * n + j];
    }
  }
  return true;
}

void dense_solve(const std::vector<double>& l, int n, double* b) {
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

double dense_criterion(const std::vector<DenseCluster>& cl, int p, int q, double sigma2, const double* g) {
  std::vector<double> xvx(p * p, 0.0), xvy(p, 0.0);
  double yvy = 0, logdet = 0;
  int total = 0;
  for (const auto& c : cl) {
    const int n = c.n;
    total += n;
    std::vector<double> v(n * n);
    for (int a = 0; a < n; ++a)
      for (int b = 0; b < n; ++b) {
        double zgz = 0;
        for (int i = 0; i < q; ++i)
          for (int j = 0; j < q; ++j) zgz += c.z[a * q + i] * g[i * q + j] * c.z[b * q + j];
        v[a * n + b] = sigma2 * ((a == b ? 1 : 0) + zgz);
      }
    if (!dense_chol(v, n)) return NAN;
    for (int a = 0; a < n; ++a) logdet += 2 * std::log(v[a * n + a]);
    std::vector<double> vy(c.y);
    dense_solve(v, n, vy.data());
    for (int a = 0; a < n; ++a) yvy += c.y[a] * vy[a];
    for (int col = 0; col < p; ++col) {
      std::vector<double> vx(n);
      for (int a = 0; a < n; ++a) vx[a] = c.x[a * p + col];
      dense_solve(v, n, vx.data());
      for (int a = 0; a < n; ++a) {
        xvy[col] += vx[a] * c.y[a];
        for (int row = 0; row < p; ++row) xvx[row * p + col] += c.x[a * p + row] * vx[a];
      }
    }
  }
  std::vector<double> l(xvx);
  if (!dense_chol(l, p)) return NAN;
  double logdet_m = 0;
  for (int a = 0; a < p; ++a) logdet_m += 2 * std::log(l[a * p + a]);
  std::vector<double> beta(xvy);
  dense_solve(l, p, beta.data());
  double fitted = 0;
  for (int a = 0; a < p; ++a) fitted += xvy[a] * beta[a];
  return logdet + logdet_m + (yvy - fitted) + (total - p) * 1.8378770664093454836;  // log 2π
}

}  // namespace

// ------------------------------------------------------------------ random numbers

TEST_CASE("sim: random streams are reproducible and distinct") {
  sim::Rng a(7, 1, 2, 3), b(7, 1, 2, 3), c(7, 1, 2, 4), d(8, 1, 2, 3);
  bool same = true, differ_c = false, differ_d = false;
  for (int i = 0; i < 100; ++i) {
    const auto x = a.next(), y = b.next(), u = c.next(), w = d.next();
    same = same && x == y;
    differ_c = differ_c || x != u;
    differ_d = differ_d || x != w;
  }
  CHECK(same);
  CHECK(differ_c);
  CHECK(differ_d);
  CHECK(sim::stream_seed(1, 2, 3) != sim::stream_seed(1, 3, 2));  // order matters
}

TEST_CASE("sim: variates have the right moments") {
  sim::Rng rng(2026);
  double m, v;
  moments([&] { return rng.uniform(); }, 400000, m, v);
  CHECK_CLOSE(m, 0.5, 0, 0.002);
  CHECK_CLOSE(v, 1.0 / 12, 0, 0.001);
  moments([&] { return rng.normal(); }, 400000, m, v);
  CHECK_CLOSE(m, 0, 0, 0.008);
  CHECK_CLOSE(v, 1, 0, 0.01);
  for (double shape : {0.4, 2.5, 30.0}) {
    moments([&] { return rng.gamma(shape); }, 200000, m, v);
    CHECK_CLOSE(m, shape, 0.01, 0);
    CHECK_CLOSE(v, shape, 0.03, 0);
  }
  moments([&] { return rng.chisq(7); }, 200000, m, v);
  CHECK_CLOSE(m, 7, 0.01, 0);
  CHECK_CLOSE(v, 14, 0.03, 0);
  for (auto np : {std::make_pair(1200L, 0.3), std::make_pair(50L, 0.9), std::make_pair(20L, 0.05)}) {
    const double n = static_cast<double>(np.first), p = np.second;
    moments([&] { return static_cast<double>(rng.binomial(np.first, p)); }, 200000, m, v);
    CHECK_CLOSE(m, n * p, 0.005, 0.01);
    CHECK_CLOSE(v, n * p * (1 - p), 0.03, 0.01);
  }
}

TEST_CASE("sim: binomial variates follow the binomial distribution") {
  // inversion (n·min(p, 1 − p) < 10) and transformed rejection (BTRS) against the exact probabilities,
  // by a chi-square statistic over cells with at least 20 expected draws (the tails pooled)
  sim::Rng rng(31);
  const long draws = 200000;
  for (auto np : {std::make_pair(30L, 0.2), std::make_pair(25L, 0.4), std::make_pair(200L, 0.85),
                  std::make_pair(5000L, 0.03), std::make_pair(1000000L, 0.5)}) {
    const long n = np.first;
    const double p = np.second;
    std::vector<long> count(static_cast<std::size_t>(n + 1), 0);
    bool in_range = true;
    for (long i = 0; i < draws; ++i) {
      const long k = rng.binomial(n, p);
      in_range = in_range && k >= 0 && k <= n;
      if (in_range) ++count[static_cast<std::size_t>(k)];
    }
    CHECK(in_range);
    double chi2 = 0, e = 0, o = 0, seen = 0;  // expected and observed draws of the open cell; expected so far
    int cells = 0;
    for (long k = 0; k <= n; ++k) {
      const double ek = draws * std::exp(std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0) +
                                         static_cast<double>(k) * std::log(p) +
                                         static_cast<double>(n - k) * std::log1p(-p));
      e += ek;
      seen += ek;
      o += static_cast<double>(count[static_cast<std::size_t>(k)]);
      if (k == n || (e >= 20 && draws - seen >= 20)) {  // the last cell takes the upper tail
        chi2 += (o - e) * (o - e) / e;
        ++cells;
        e = o = 0;
      }
    }
    const double df = cells - 1;
    CHECK(chi2 < df + 4 * std::sqrt(2 * df));
  }
}

// ------------------------------------------------------------------ REML

TEST_CASE("lmm: REML of a balanced random-intercept design equals the ANOVA estimates") {
  // k clusters per group of m subjects; y = μ + δ·group + b + e
  const int k = 6, m = 5;
  sim::Rng rng(11);
  lmm::Data d;
  d.p = 2;
  d.q = 1;
  d.patterns.resize(2);
  double x[2 * m * 2], z[m];
  std::vector<double> means[2];
  double ssw = 0;
  for (int g = 0; g < 2; ++g) {
    for (int j = 0; j < m; ++j) {
      x[2 * j] = 1;
      x[2 * j + 1] = g == 0 ? 1 : 0;
      z[j] = 1;
    }
    d.patterns[g].set_design(m, 2, 1, x, z);
    for (int c = 0; c < k; ++c) {
      double y[m], mean = 0;
      const double b = 0.8 * rng.normal();
      for (int j = 0; j < m; ++j) mean += (y[j] = 1 + (g == 0 ? 0.5 : 0) + b + rng.normal()) / m;
      for (int j = 0; j < m; ++j) ssw += (y[j] - mean) * (y[j] - mean);
      means[g].push_back(mean);
      d.patterns[g].add_cluster(2, 1, x, z, y);
    }
  }
  double ssb = 0, grand[2] = {0, 0};
  for (int g = 0; g < 2; ++g) {
    for (double v : means[g]) grand[g] += v / k;
    for (double v : means[g]) ssb += m * (v - grand[g]) * (v - grand[g]);
  }
  const double msw = ssw / (2 * k * (m - 1)), msb = ssb / (2 * k - 2);
  const lmm::Fit f = lmm::fit(d);
  CHECK(f.ok);
  CHECK(msb > msw);  // an interior estimate for this seed
  CHECK_CLOSE(f.sigma2, msw, 1e-7, 0);
  CHECK_CLOSE(f.sigma2 * f.g[0], (msb - msw) / m, 1e-6, 0);
  CHECK_CLOSE(f.beta[1], grand[0] - grand[1], 1e-10, 1e-12);
  // the Wald test of the group effect is the two-sample t-test on the cluster means, df = 2k − 2
  const double c[2] = {0, 1};
  const lmm::Wald w = lmm::wald(d, f, c, 0, 0);
  CHECK(!w.fallback);
  CHECK_CLOSE(w.se, std::sqrt(msb / m * 2.0 / k), 1e-6, 0);
  CHECK_CLOSE(w.df, 2 * k - 2, 1e-4, 0);
}

TEST_CASE("lmm: random slopes on complete balanced data reproduce the two-stage analysis") {
  // GLS = OLS on each subject's line; the REML variance of the slope difference is the pooled
  // variance of the subjects' slopes, so the test is the two-sample t-test with N − 2 df.
  const int n = 25, K = 5;
  const double times[K] = {0, 1, 2, 3, 4};
  lmm::Data d;
  d.p = 4;
  d.q = 2;
  d.patterns.resize(2);
  double x[K * 4], z[K * 2];
  std::vector<double> slopes[2];
  sim::Rng rng(5);
  for (int g = 0; g < 2; ++g) {
    const double ind = g == 0 ? 1 : 0;
    for (int j = 0; j < K; ++j) {
      x[4 * j] = 1, x[4 * j + 1] = ind, x[4 * j + 2] = times[j], x[4 * j + 3] = ind * times[j];
      z[2 * j] = 1, z[2 * j + 1] = times[j];
    }
    d.patterns[g].set_design(K, 4, 2, x, z);
    for (int u = 0; u < n; ++u) {
      double y[K];
      const double b0 = rng.normal(), b1 = 0.5 * rng.normal();
      double st = 0, sy = 0, sty = 0, stt = 0;
      for (int j = 0; j < K; ++j) {
        y[j] = 2 + (0.4 * ind + b1) * times[j] + b0 + 0.5 * rng.normal();
        st += times[j], sy += y[j], sty += times[j] * y[j], stt += times[j] * times[j];
      }
      slopes[g].push_back((sty - st * sy / K) / (stt - st * st / K));
      d.patterns[g].add_cluster(4, 2, x, z, y);
    }
  }
  double mean[2] = {0, 0}, ss = 0;
  for (int g = 0; g < 2; ++g) {
    for (double v : slopes[g]) mean[g] += v / n;
    for (double v : slopes[g]) ss += (v - mean[g]) * (v - mean[g]);
  }
  const lmm::Fit f = lmm::fit(d);
  CHECK(f.ok);
  CHECK(!f.boundary);
  const double c[4] = {0, 0, 0, 1};
  const lmm::Wald w = lmm::wald(d, f, c, 0, 0);
  CHECK_CLOSE(w.estimate, mean[0] - mean[1], 1e-9, 1e-12);
  CHECK_CLOSE(w.se, std::sqrt(ss / (2 * n - 2) * 2.0 / n), 1e-5, 0);
  CHECK_CLOSE(w.df, 2 * n - 2, 1e-3, 0);
}

TEST_CASE("lmm: random effects that dwarf the residual SD leave the two-stage analysis intact") {
  // the same design with the subjects' intercepts 10⁸ or 10¹⁰⁰ residual SDs apart, or their slopes spreading
  // them 10⁵ residual SDs apart over a follow-up of 10⁶ time units: y'y, and Z'y, are then almost all
  // between-subject variation, which once lost the digits of the residual sum of squares (the slopes' standard
  // error came out 1.6 and 48 times too large, and any number with intercepts 10¹⁰⁰ apart)
  const int n = 25, K = 5;
  struct Scale {
    double intercept, time, slope;
  };
  for (const Scale sc : {Scale{1, 1, 1}, Scale{1e8, 1, 1}, Scale{1e100, 1, 1}, Scale{1, 2.5e5, 2.5e5}}) {
    lmm::Data d;
    d.p = 4;
    d.q = 2;
    d.patterns.resize(2);
    double times[K], x[2][K * 4], z[K * 2];
    for (int j = 0; j < K; ++j) times[j] = j * sc.time;
    std::vector<double> slopes[2];
    sim::Rng rng(5);
    for (int g = 0; g < 2; ++g) {
      const double ind = g == 0 ? 1 : 0;
      for (int j = 0; j < K; ++j) {
        x[g][4 * j] = 1, x[g][4 * j + 1] = ind, x[g][4 * j + 2] = times[j], x[g][4 * j + 3] = ind * times[j];
        z[2 * j] = 1, z[2 * j + 1] = times[j];
      }
      d.patterns[g].set_design(K, 4, 2, x[g], z);
      for (int u = 0; u < n; ++u) {
        // y = Z u + e: the subject's line and the errors around it
        const double line[2] = {2 + sc.intercept * rng.normal(), (0.4 * ind + 0.5 * rng.normal()) * sc.slope / sc.time};
        double e[K], st = 0, se = 0, stt = 0, ste = 0;
        for (int j = 0; j < K; ++j) {
          e[j] = 0.5 * rng.normal();
          st += times[j], se += e[j], stt += times[j] * times[j], ste += times[j] * e[j];
        }
        slopes[g].push_back(line[1] + (ste - st * se / K) / (stt - st * st / K));
        d.patterns[g].add_cluster(4, 2, x[g], z, line, e);
      }
    }
    double mean[2] = {0, 0}, ss = 0;
    for (int g = 0; g < 2; ++g) {
      for (double v : slopes[g]) mean[g] += v / n;
      for (double v : slopes[g]) ss += (v - mean[g]) * (v - mean[g]);
    }
    const double theta0[3] = {2 * sc.intercept, 0, sc.slope / sc.time};
    const lmm::Fit f = lmm::fit(d, theta0);
    CHECK(f.ok);
    const double c[4] = {0, 0, 0, 1};
    const lmm::Wald w = lmm::wald(d, f, c, 0, 0);
    CHECK_CLOSE(w.estimate, mean[0] - mean[1], 1e-7, 0);
    CHECK_CLOSE(w.se, std::sqrt(ss / (2 * n - 2) * 2.0 / n), 1e-5, 0);
    if (sc.intercept < 1e100) CHECK_CLOSE(w.df, 2 * n - 2, 1e-3, 0);  // beyond, the Hessian over/underflows
  }
}

TEST_CASE("lmm: the pattern formulas equal a dense REML likelihood with dropout") {
  // unbalanced data (monotone dropout), random intercepts and slopes, compared with V built in full
  const double times[4] = {0, 1, 3, 6};
  sim::Rng rng(99);
  lmm::Data d;
  d.p = 4;
  d.q = 2;
  d.patterns.resize(8);
  std::vector<DenseCluster> dense;
  for (int g = 0; g < 2; ++g)
    for (int k = 1; k <= 4; ++k) {
      double x[16], z[8];
      for (int j = 0; j < k; ++j) {
        x[4 * j] = 1, x[4 * j + 1] = g == 0, x[4 * j + 2] = times[j], x[4 * j + 3] = (g == 0) * times[j];
        z[2 * j] = 1, z[2 * j + 1] = times[j];
      }
      d.patterns[g * 4 + k - 1].set_design(k, 4, 2, x, z);
    }
  for (int g = 0; g < 2; ++g)
    for (int u = 0; u < 14; ++u) {
      const int k = 1 + static_cast<int>(rng.uniform() * 4);
      DenseCluster c{k, {}, {}, {}};
      const double b0 = rng.normal(), b1 = 0.3 * rng.normal();
      for (int j = 0; j < k; ++j) {
        const double t = times[j];
        c.x.insert(c.x.end(), {1.0, g == 0 ? 1.0 : 0.0, t, (g == 0) * t});
        c.z.insert(c.z.end(), {1.0, t});
        c.y.push_back(0.2 * (g == 0) * t + b0 + b1 * t + 0.7 * rng.normal());
      }
      d.patterns[g * 4 + k - 1].add_cluster(4, 2, c.x.data(), c.z.data(), c.y.data());
      dense.push_back(c);
    }
  const double gs[3][4] = {{1.0, 0.2, 0.2, 0.3}, {0.5, -0.1, -0.1, 0.05}, {0.0, 0.0, 0.0, 0.2}};
  for (const auto& g : gs)
    for (double s2 : {0.4, 1.3})
      CHECK_CLOSE(lmm::criterion(d, s2, g), dense_criterion(dense, 4, 2, s2, g), 1e-10, 1e-9);
  // the fit minimizes the criterion: no nearby admissible (σ², G) does better
  const lmm::Fit f = lmm::fit(d);
  CHECK(f.ok);
  const double best = lmm::criterion(d, f.sigma2, f.g);
  CHECK_CLOSE(best, f.criterion, 1e-9, 1e-9);
  for (double h : {-1e-3, 1e-3}) {
    CHECK(lmm::criterion(d, f.sigma2 * (1 + h), f.g) >= best - 1e-9);
    for (int i : {0, 1, 3}) {
      double g[4] = {f.g[0], f.g[1], f.g[2], f.g[3]};
      g[i] += h * (0.1 + std::fabs(g[i]));
      g[2] = g[1];
      if (g[0] < 0 || g[3] < 0 || g[1] * g[1] > g[0] * g[3]) continue;  // not a covariance matrix
      CHECK(lmm::criterion(d, f.sigma2, g) >= best - 1e-7);
    }
  }
}

// ------------------------------------------------------------------ simulated power

TEST_CASE("sim: quick approximations equal the closed-form procedures") {
  const json::Value a = first_row(
      R"({"procedure":"sim.lmm.cluster","solve_for":"power","options":{"approximate":true},)"
      R"("inputs":{"k1":"8","m":"12","icc":"0.08","mu1":"0.45","df":"fixed"}})");
  const json::Value b = first_row(
      R"({"procedure":"means.two.cluster","solve_for":"power",)"
      R"("inputs":{"k1":"8","m":"12","icc":"0.08","mu1":"0.45"}})");
  CHECK_CLOSE(num(a, "power"), num(b, "power"), 1e-12, 0);
  CHECK(a.get("method").to_text() == "approximation");
  const json::Value c = first_row(
      R"({"procedure":"sim.lmm.longitudinal","solve_for":"n","options":{"approximate":true},)"
      R"("inputs":{"power":"0.9","delta":"0.3","sd1":"0.4","sde":"1.2","m":"6","duration":"5","df":"fixed"}})");
  const json::Value e = first_row(
      R"({"procedure":"means.two.slopes","solve_for":"n",)"
      R"("inputs":{"power":"0.9","delta":"0.3","sdb":"0.4","sde":"1.2","m":"6","duration":"5"}})");
  CHECK_EQ(num(c, "n1"), num(e, "n1"));
  CHECK_CLOSE(num(c, "power"), num(e, "power"), 1e-12, 0);
}

TEST_CASE("sim: simulated power agrees with exact power for balanced designs") {
  // equal cluster sizes and complete balanced longitudinal data: the mixed-model test is the exact
  // t-test on cluster means (or on subjects' slopes), so the simulation must agree with the closed form
  // within its Monte Carlo error (4 standard errors).
  const json::Value a = first_row(
      R"({"procedure":"sim.lmm.cluster","solve_for":"power",)"
      R"("inputs":{"k1":"8","m":"10","icc":"0.2","mu1":"0.6","df":"fixed","reps":"20000","seed":"1"}})");
  const json::Value b = first_row(
      R"({"procedure":"means.two.cluster","solve_for":"power","inputs":{"k1":"8","m":"10","icc":"0.2","mu1":"0.6"}})");
  CHECK_CLOSE(num(a, "power"), num(b, "power"), 0, 4 * num(a, "power_se"));
  CHECK(num(a, "power_lo") < num(a, "power") && num(a, "power") < num(a, "power_hi"));
  const json::Value c = first_row(
      R"({"procedure":"sim.lmm.longitudinal","solve_for":"power",)"
      R"("inputs":{"n1":"30","delta":"0.4","sd1":"0.5","sde":"1","df":"fixed","reps":"20000","seed":"2"}})");
  const json::Value e = first_row(
      R"({"procedure":"means.two.slopes","solve_for":"power",)"
      R"("inputs":{"n1":"30","delta":"0.4","sdb":"0.5","sde":"1"}})");
  CHECK_CLOSE(num(c, "power"), num(e, "power"), 0, 4 * num(c, "power_se"));
  CHECK_CLOSE(num(c, "mean_delta"), 0.4, 0, 0.01);
  CHECK_EQ(num(c, "failures"), 0.0);
}

TEST_CASE("sim: the type I error is close to alpha") {
  for (const char* req :
       {R"({"procedure":"sim.lmm.cluster","solve_for":"power",)"
        R"("inputs":{"k1":"6","m":"15","cv":"0.7","icc":"0.05","mu1":"0","reps":"20000","seed":"3"}})",
        R"({"procedure":"sim.lmm.longitudinal","solve_for":"power",)"
        R"("inputs":{"n1":"20","delta":"0","sd1":"0.5","rho":"0.3","dropout":"0.3","reps":"20000","seed":"4"}})"}) {
    const json::Value r = first_row(req);
    CHECK_CLOSE(num(r, "power"), 0.05, 0, 0.0065);  // about 4 Monte Carlo standard errors
  }
}

TEST_CASE("sim: results depend on the seed only, not on the number of threads") {
  const std::string body =
      R"("procedure":"sim.lmm.longitudinal","solve_for":"power",)"
      R"("inputs":{"n1":"12","delta":"0.5","dropout":"0.2","reps":"400","seed":"77"})";
  const json::Value one = first_row("{" + body + R"(,"options":{"threads":1}})");
  const json::Value four = first_row("{" + body + R"(,"options":{"threads":4}})");
  CHECK_EQ(num(one, "power"), num(four, "power"));
  CHECK_EQ(num(one, "mean_delta"), num(four, "mean_delta"));
  CHECK_EQ(num(one, "mean_df"), num(four, "mean_df"));
}

TEST_CASE("sim: the sample-size search returns the smallest size reaching the target") {
  const json::Value r = first_row(
      R"({"procedure":"sim.lmm.cluster","solve_for":"k",)"
      R"("inputs":{"power":"0.8","m":"15","cv":"0.5","icc":"0.1","mu1":"0.5","reps":"1000","seed":"9"}})");
  const double k1 = num(r, "k1");
  CHECK(num(r, "power") >= 0.8);
  // with the same seed, one cluster fewer falls short (common random numbers make this reproducible)
  const json::Value below = first_row(
      R"({"procedure":"sim.lmm.cluster","solve_for":"power","inputs":{"k1":")" + text::format_g(k1 - 1, 10) +
      R"(","m":"15","cv":"0.5","icc":"0.1","mu1":"0.5","reps":"1000","seed":"9"}})");
  CHECK(num(below, "power") < 0.8);
  const json::Value at = first_row(
      R"({"procedure":"sim.lmm.cluster","solve_for":"power","inputs":{"k1":")" + text::format_g(k1, 10) +
      R"(","m":"15","cv":"0.5","icc":"0.1","mu1":"0.5","reps":"1000","seed":"9"}})");
  CHECK_EQ(num(at, "power"), num(r, "power"));
}

TEST_CASE("sim: progress is reported and a simulation can be cancelled") {
  Token token;
  run(R"({"procedure":"sim.lmm.longitudinal","solve_for":"power","inputs":{"reps":"300"}})", &token);
  CHECK_CLOSE(token.progress(), 1.0, 0, 1e-12);

  Token stop;
  std::atomic<bool> finished{false};
  json::Value res;
  std::thread worker([&] {
    res = run(R"({"procedure":"sim.lmm.longitudinal","solve_for":"power","options":{"power_curve":false},)"
              R"("inputs":{"n1":"400","m":"8","reps":"200000"}})",
              &stop);
    finished = true;
  });
  std::this_thread::sleep_for(std::chrono::milliseconds(300));
  const double partial = stop.progress();
  stop.cancel();
  worker.join();
  CHECK(finished.load());
  CHECK(partial > 0 && partial < 1);
  CHECK(res.get("cancelled").as_bool(false));
}
