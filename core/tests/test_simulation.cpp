// Simulation procedures: random streams and variates, the REML fitter against
// closed forms and a brute-force likelihood, Satterthwaite's degrees of
// freedom, simulated power against exact power, type I error, reproducibility
// and the sample-size search.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <thread>
#include <vector>

#include "engine/runner.hpp"
#include "procedures/sim_common.hpp"
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

TEST_CASE("lmm: clusters kept by arm and size fit as their patterns do, to the bit") {
  // a cluster-randomized study with sizes that vary (a few shared, most of their own), as one Pattern per arm and
  // size in the order the sizes first appear, and as the same sizes in two Arms
  for (const double icc : {0.0, 0.05, 0.6}) {
    sim::Rng rng(29);
    lmm::Data by_pattern, by_arm;
    by_pattern.p = by_arm.p = 2;
    by_pattern.q = by_arm.q = 1;
    by_arm.arms.resize(2);
    for (int g = 0; g < 2; ++g) {
      std::vector<long> seen;
      lmm::Arm& arm = by_arm.arms[g];
      arm.x[0] = 1, arm.x[1] = g == 0 ? 1 : 0;
      for (int c = 0; c < 40; ++c) {
        const long m = c % 3 == 0 ? 5 : 1 + static_cast<long>(rng.uniform() * 300);
        const double dm = static_cast<double>(m);
        const double mean = std::sqrt(icc) * rng.normal() + std::sqrt((1 - icc) / dm) * rng.normal();
        const double ss = m > 1 ? (1 - icc) * rng.chisq(dm - 1) : 0;
        std::size_t i = std::find(seen.begin(), seen.end(), m) - seen.begin();
        if (i == seen.size()) {
          seen.push_back(m);
          lmm::Pattern pt;
          pt.n = m;
          pt.zz[0] = dm;
          pt.zx[0] = dm, pt.zx[1] = dm * arm.x[1];
          pt.set_coefficients(2, 1);
          by_pattern.patterns.push_back(pt);
          arm.sizes.push_back({dm});
        }
        const double xy[2] = {0, 0};
        by_pattern.patterns[by_pattern.patterns.size() - seen.size() + i].add_sums(2, 1, xy, &mean, ss);
        lmm::Arm::Size& sz = arm.sizes[i];
        sz.count += 1, sz.b += mean, sz.bb += mean * mean, sz.yy += ss;
      }
    }
    CHECK(by_arm.arms[0].sizes.size() + by_arm.arms[1].sizes.size() == by_pattern.patterns.size());
    CHECK(by_arm.observations() == by_pattern.observations());
    CHECK(by_arm.clusters() == 80);
    const double theta0 = std::sqrt(icc / (1 - icc));
    const lmm::Fit a = lmm::fit(by_pattern, &theta0), b = lmm::fit(by_arm, &theta0);
    CHECK(a.ok && b.ok);
    CHECK(a.theta[0] == b.theta[0] && a.sigma2 == b.sigma2 && a.criterion == b.criterion);
    for (int i = 0; i < 2; ++i) CHECK(a.beta[i] == b.beta[i]);
    for (int i = 0; i < 4; ++i) CHECK(a.cov[i] == b.cov[i]);
    const double c[2] = {0, 1};
    for (int method = 0; method < 3; ++method) {
      const lmm::Wald wa = lmm::wald(by_pattern, a, c, method, 78), wb = lmm::wald(by_arm, b, c, method, 78);
      CHECK(wa.estimate == wb.estimate && wa.se == wb.se && wa.df == wb.df && wa.fallback == wb.fallback);
    }
  }
  // arms are clusters with a random intercept alone
  lmm::Data d;
  d.p = 2;
  d.q = 2;
  d.arms.resize(1);
  d.arms[0].x[0] = 1;
  d.arms[0].sizes.push_back({10, 5, 0.1, 0.2, 40});
  CHECK(!lmm::fit(d).ok);
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
      double rx[K * 4];
      d.patterns[g].within_rows(4, 2, x[g], z, rx);
      for (int u = 0; u < n; ++u) {
        // y = Z u + e: the subject's line and the errors around it
        const double line[2] = {2 + sc.intercept * rng.normal(), (0.4 * ind + 0.5 * rng.normal()) * sc.slope / sc.time};
        double e[K], st = 0, se = 0, stt = 0, ste = 0;
        for (int j = 0; j < K; ++j) {
          e[j] = 0.5 * rng.normal();
          st += times[j], se += e[j], stt += times[j] * times[j], ste += times[j] * e[j];
        }
        slopes[g].push_back(line[1] + (ste - st * se / K) / (stt - st * st / K));
        d.patterns[g].add_effects(4, 2, rx, z, line, e);
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
    CHECK(!w.fallback);  // the Hessian in units of its steps: their squares once overflowed at 10¹⁰⁰
    CHECK_CLOSE(w.df, 2 * n - 2, 1e-3, 0);
  }
}

TEST_CASE("lmm: two visits close together and far from where the random effects are given keep what they tell") {
  // visits at −1000001, −1000000 and 0 (in units of τ = 1000001) and the random effects given at time 0: from
  // there, the first two visits' Z'Z is singular to within 10⁻¹², and was once taken as of rank 1, which dropped
  // the random slope's spread in their difference (a power of 96% for an exact 70.5%).  With Z from the first
  // visit and G at time 0 through S, the GLS variance of the slope difference is that of a 50-digit calculation;
  // with a slope SD that the two visits cannot see, the direction is left within the subject (rank 1), and the
  // variance is as close
  const double tau = 1000001, s[4] = {1, -1, 0, 1};  // the first visit is τ before time 0
  struct Case {
    double sd1, exact;  // the exact variance per τ², 10 subjects of each group with two visits and 40 with three
  };
  for (const Case cs : {Case{10, 4000008000004.08}, Case{1e-9, 0.0733334511111115}}) {
    const double g[4] = {1, 0, 0, cs.sd1 * tau * (cs.sd1 * tau)};
    lmm::Data d;
    d.p = 4;
    d.q = 2;
    for (int grp = 0; grp < 2; ++grp)
      for (int k = 2; k <= 3; ++k) {
        const double z[6] = {1, 0, 1, 1 / tau, 1, 1};  // time from the first visit
        const double ind = grp == 0 ? 1 : 0, c[8] = {1, ind, 0, 0, 0, 0, 1, ind};
        lmm::Pattern pt;
        pt.set_spanned(k, 4, 2, z, c, s, g);
        CHECK(pt.rank_z == (k == 2 && cs.sd1 < 1 ? 1 : 2));
        pt.count = k == 2 ? 10 : 40;
        d.patterns.push_back(pt);
      }
    double cov[lmm::PMAX * lmm::PMAX];
    CHECK(lmm::beta_covariance(d, 1, g, cov, g[0] * g[3]));
    CHECK_CLOSE(cov[15], cs.exact, 1e-10, 0);
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

TEST_CASE("lmm: Satterthwaite's degrees of freedom do not depend on where time starts") {
  // subjects whose lines all but cross ten time units before the first visit (an intercept SD of 100 at time 0,
  // slopes of SD 10 correlated 0.99 with it), fitted with time from 0 and from the crossing: at time 0 the fitted
  // correlation is about 0.99, and steps of 10⁻⁴ there once left the Hessian not positive definite in every study
  // (the test then used N − 2 without a word); taken where the fitted D is diagonal, the degrees of freedom are
  // the same from either origin
  const double times[4] = {0, 1, 2, 3}, l00 = 100, l10 = 9.9, l11 = 10 * std::sqrt(1 - 0.99 * 0.99);
  lmm::Wald w[2];
  for (int o = 0; o < 2; ++o) {
    const double origin = o == 0 ? 0 : -9.9;
    lmm::Data d;
    d.p = 4;
    d.q = 2;
    d.patterns.resize(8);
    std::vector<double> x[8], z[8];
    for (int g = 0; g < 2; ++g)
      for (int k = 1; k <= 4; ++k) {
        const int i = g * 4 + k - 1;
        for (int j = 0; j < k; ++j) {
          const double t = times[j] - origin;
          x[i].insert(x[i].end(), {1.0, g == 0 ? 1.0 : 0.0, t, g == 0 ? t : 0.0});
          z[i].insert(z[i].end(), {1.0, t});
        }
        d.patterns[i].set_design(k, 4, 2, x[i].data(), z[i].data());
      }
    sim::Rng rng(11);
    for (int g = 0; g < 2; ++g)
      for (int u = 0; u < 100; ++u) {
        const double z0 = rng.normal(), z1 = rng.normal(), b0 = l00 * z0, b1 = l10 * z0 + l11 * z1;
        const int k = 1 + static_cast<int>(rng.uniform() * 4), i = g * 4 + k - 1;
        double y[4];
        for (int j = 0; j < k; ++j) y[j] = b0 + b1 * times[j] + rng.normal();
        d.patterns[i].add_cluster(4, 2, x[i].data(), z[i].data(), y);
      }
    // Λ at the assumed variances from this origin
    const double a = l00 + origin * l10, b = origin * l11;
    const double g00 = a * a + b * b, g01 = a * l10 + b * l11, g11 = l10 * l10 + l11 * l11;
    const double theta0[3] = {std::sqrt(g00), g01 / std::sqrt(g00), std::sqrt(g11 - g01 * g01 / g00)};
    const lmm::Fit f = lmm::fit(d, theta0);
    CHECK(f.ok);
    const double c[4] = {0, 0, 0, 1};
    w[o] = lmm::wald(d, f, c, 0, 198);
    CHECK(!w[o].fallback);
  }
  CHECK_CLOSE(w[1].estimate, w[0].estimate, 1e-7, 0);
  CHECK_CLOSE(w[1].se, w[0].se, 1e-6, 0);
  CHECK_CLOSE(w[1].df, w[0].df, 1e-5, 0);
}

TEST_CASE("lmm: a fit that cannot tell a group's slope apart has no estimate") {
  // two visits far from the time origin, and every subject of group 1 seen only at the first: nothing tells that
  // group's slope from its level, and X'V⁻¹X is singular but for rounding; the estimate came out as whatever the
  // rounding made it (exactly 0 or 1 in the simulations), with a standard error of 10⁷
  const double times[2] = {1000.3, 1000.7};
  lmm::Data d;
  d.p = 4;
  d.q = 2;
  d.patterns.resize(4);
  std::vector<double> x[4], z[4];
  for (int g = 0; g < 2; ++g)
    for (int k = 1; k <= 2; ++k) {
      const int i = g * 2 + k - 1;
      for (int j = 0; j < k; ++j) {
        x[i].insert(x[i].end(), {1.0, g == 0 ? 1.0 : 0.0, times[j], g == 0 ? times[j] : 0.0});
        z[i].insert(z[i].end(), {1.0, times[j]});
      }
      d.patterns[i].set_design(k, 4, 2, x[i].data(), z[i].data());
    }
  sim::Rng rng(3);
  for (int g = 0; g < 2; ++g)
    for (int u = 0; u < 12; ++u) {
      const int k = g == 0 ? 1 : 1 + u % 2, i = g * 2 + k - 1;
      const double b0 = rng.normal(), b1 = 0.5 * rng.normal();
      const double y[2] = {b0 + b1 * times[0] + rng.normal(), b0 + b1 * times[1] + rng.normal()};
      d.patterns[i].add_cluster(4, 2, x[i].data(), z[i].data(), y);
    }
  CHECK(!lmm::fit(d).ok);
  d.two_rows = true;
  CHECK(!lmm::fit(d).ok);
}

TEST_CASE("lmm: with two visits, Satterthwaite's degrees of freedom are those of the two measurements' covariance") {
  // with random intercepts and slopes and two visits, the data tell only the 2 × 2 covariance of a subject's two
  // measurements, not σ² and D apart: the fit holds D at rank 1, where (σ², w) reach each covariance once, and the
  // estimate, its standard error and the degrees of freedom are those of a 40-digit REML fit in that covariance
  // (Satterthwaite's df by finite differences in its three entries: 39.0415755857863).  Fitted in full, the fit
  // stops anywhere on a ridge of equal criteria, the Hessian in (σ², D) is singular, and this study fell back to
  // the given degrees of freedom
  const double times[2] = {0, 1};
  lmm::Data d;
  d.p = 4;
  d.q = 2;
  d.patterns.resize(4);
  std::vector<double> x[4], z[4];
  for (int g = 0; g < 2; ++g)
    for (int k = 1; k <= 2; ++k) {
      const int i = g * 2 + k - 1;
      for (int j = 0; j < k; ++j) {
        x[i].insert(x[i].end(), {1.0, g == 0 ? 1.0 : 0.0, times[j], g == 0 ? times[j] : 0.0});
        z[i].insert(z[i].end(), {1.0, times[j]});
      }
      d.patterns[i].set_design(k, 4, 2, x[i].data(), z[i].data());
    }
  for (int g = 0; g < 2; ++g)
    for (int u = 0; u < 30; ++u) {
      const int k = u % 3 == 0 ? 1 : 2, i = g * 2 + k - 1;
      const double b0 = std::cos(1.7 * u + g + 0.5), b1 = 0.5 * std::sin(2.3 * u + g);
      double y[2];
      for (int j = 0; j < k; ++j) y[j] = b0 + b1 * times[j] + 0.7 * std::sin(3.1 * u + 1.3 * j + g);
      d.patterns[i].add_cluster(4, 2, x[i].data(), z[i].data(), y);
    }
  d.two_rows = true;
  const lmm::Fit f = lmm::fit(d);
  CHECK(f.ok);
  CHECK(f.theta[2] == 0);
  const double c[4] = {0, 0, 0, 1};
  const lmm::Wald w = lmm::wald(d, f, c, 0, 58);
  CHECK(!w.fallback);
  CHECK_CLOSE(w.estimate, 0.00421217577174542, 0, 1e-8);
  CHECK_CLOSE(w.se, 0.242236572215089, 1e-7, 0);
  CHECK_CLOSE(w.df, 39.0415755857863, 1e-6, 0);
}

// ------------------------------------------------------------------ simulated power

TEST_CASE("sim: a power is refused only when no study was analyzed, or the arithmetic broke down in most") {
  // studies that fail by design (no wins or no losses, a model that cannot be fitted to the few subjects left)
  // count as not rejecting, a fair account of a small design; refusing once more than half failed turned such
  // designs into errors that said to check for an extreme value
  namespace simc = procs::simc;
  Scenario s;
  s.solve = "power";
  s.num["seed"] = 1;
  auto outcome = [&](long failures, simc::Failures kind) {
    sim::Summary sm;
    sm.reps = 400;
    sm.failures = failures;
    sm = sim::with_rejections(sm, 0);
    Row r;
    try {
      simc::put_summary(s, r, sm, kind, "each had no wins or no losses");
    } catch (const CalcError& e) {
      return std::string(e.what());
    }
    CHECK(r.num("power") == 0 && r.num("failures") == static_cast<double>(failures));
    return std::string();
  };
  CHECK(outcome(399, simc::Failures::Design).empty());
  CHECK(outcome(400, simc::Failures::Design) ==
        "the simulated power cannot be computed: none of the 400 simulated studies could be analyzed (each had no "
        "wins or no losses)");
  CHECK(outcome(200, simc::Failures::Numerical).empty());
  CHECK(outcome(201, simc::Failures::Numerical) ==
        "the simulated power cannot be computed: 201 of the 400 simulated studies could not be analyzed; check for "
        "an extreme value");
}

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

TEST_CASE("sim: the longitudinal approximation is the exact power with slopes that dwarf the intercepts") {
  // slopes 10⁶ to 10¹⁵ times as variable as the intercepts at time 0, with visits away from it: measured from the
  // first visit, X'V⁻¹X was singular to within 10⁻¹³ (an answer refused, or three digits); measured from where the
  // estimates of the mean level and slope are uncorrelated, the GLS power is that of a 300-digit calculation
  struct Case {
    const char *times, *sd1, *rho, *dropout, *n1, *delta;
    double exact;
  };
  for (const Case cs : {Case{"1 2 4", "1e8", "0.5", "0.3", "10", "2.5e7", 0.0864997134942078},
                        Case{"1 2 4", "1e15", "0.5", "0.3", "10", "2.5e14", 0.0864997134942078},
                        Case{"1 2 4", "5e6", "0", "0", "100", "1414213.5623731", 0.516005273976174},
                        Case{"2020 2021 2022 2023 2024", "1e6", "0", "0", "100", "282842.712474633", 0.516005273976174}}) {
    const json::Value r = first_row(
        std::string(R"({"procedure":"sim.lmm.longitudinal","solve_for":"power","options":{"approximate":true},)") +
        R"("inputs":{"timing":"custom","df":"z","times":")" + cs.times + R"(","sd1":")" + cs.sd1 + R"(","rho":")" +
        cs.rho + R"(","dropout":")" + cs.dropout + R"(","n1":")" + cs.n1 + R"(","delta":")" + cs.delta + R"("}})");
    CHECK_CLOSE(num(r, "power"), cs.exact, 0, 1e-9);
  }
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
