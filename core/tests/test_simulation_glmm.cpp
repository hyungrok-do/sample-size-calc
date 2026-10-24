// Logistic mixed model for clustered yes/no outcomes: quadrature, fit, scale
// conversions and the simulation procedure sim.glmm.cluster.
#include <cmath>
#include <vector>

#include "engine/runner.hpp"
#include "procedures/sim_common.hpp"
#include "sim/glmm.hpp"
#include "sim/random.hpp"
#include "testing.hpp"
#include "util/json.hpp"
#include "util/numtext.hpp"

using namespace ssc;
namespace glmm = ssc::sim::glmm;

namespace {

json::Value first_row(const std::string& request) {
  Token token;
  const json::Value res = run_request(json::parse(request), &token);
  CHECK(res.get("ok").as_bool(false));
  return res.get("rows").elements().at(0);
}

double num(const json::Value& row, const char* key) { return row.get(key).as_number(NAN); }

/// log ∫ exp(s(η + σz) − m log(1 + e^(η+σz))) φ(z) dz by Simpson's rule on a fine grid.
double brute_force(double eta, double m, double s, double sigma) {
  const int n = 40000;
  const double a = -12, b = 12, h = (b - a) / n;
  std::vector<double> g(n + 1);
  double gmax = -INFINITY;
  for (int i = 0; i <= n; ++i) {
    const double z = a + i * h, t = eta + sigma * z;
    const double sp = t > 0 ? t + std::log1p(std::exp(-t)) : std::log1p(std::exp(t));
    g[i] = s * t - m * sp - 0.5 * z * z;
    gmax = std::max(gmax, g[i]);
  }
  double sum = 0;
  for (int i = 0; i <= n; ++i) sum += (i == 0 || i == n ? 1 : (i % 2 ? 4 : 2)) * std::exp(g[i] - gmax);
  return gmax + std::log(sum * h / 3) - 0.5 * std::log(2 * 3.14159265358979323846);
}

double logit(double p) { return std::log(p / (1 - p)); }

}  // namespace

TEST_CASE("glmm: Gauss-Hermite rules integrate polynomials exactly") {
  for (int n : {3, 9, 20, 40}) {
    std::vector<double> x, w;
    glmm::gauss_hermite(n, x, w);
    double m0 = 0, m2 = 0, m4 = 0, odd = 0;
    for (int k = 0; k < n; ++k) {
      m0 += w[k];
      m2 += w[k] * x[k] * x[k];
      m4 += w[k] * std::pow(x[k], 4);
      odd += w[k] * std::pow(x[k], 3);
    }
    const double sqrt_pi = 1.7724538509055160273;
    CHECK_CLOSE(m0, sqrt_pi, 1e-13, 0);
    CHECK_CLOSE(m2, sqrt_pi / 2, 1e-13, 0);
    CHECK_CLOSE(m4, 3 * sqrt_pi / 4, 1e-12, 0);
    CHECK_CLOSE(odd, 0, 0, 1e-13);
  }
}

TEST_CASE("glmm: adaptive quadrature matches the integral of each cluster") {
  // the last three have a large σ²m, where plain Newton's method for the mode cycles
  const double cases[][4] = {{-1.2, 30, 5, 0.6},    {0.4, 8, 8, 1.5},       {-3, 50, 0, 0.3},
                             {2, 100, 97, 2.5},     {-0.5, 1, 1, 1.0},      {0.0, 20, 10, 0.0},
                             {-3.5, 100, 50, 1.0},  {-3.5, 100, 50, 8.0},   {-3.46, 80, 41, 3.7}};
  for (const auto& c : cases)
    CHECK_CLOSE(glmm::cluster_loglik(c[0], c[1], c[2], c[3]), brute_force(c[0], c[1], c[2], c[3]), 0, 5e-7);
}

TEST_CASE("glmm: a large between-cluster SD is estimated without bias") {
  // ICC 0.575 on the proportion scale (σ ≈ 3.7): with many clusters the fit recovers the parameters
  const double sigma = glmm::sigma_for_icc(0.575, [](double s) { return glmm::beta_for_marginal(0.2, s); });
  const double b0 = glmm::beta_for_marginal(0.2, sigma), b1 = glmm::beta_for_marginal(0.3, sigma) - b0;
  double sum_b1 = 0, sum_sigma = 0;
  const int reps = 20;
  for (int rep = 0; rep < reps; ++rep) {
    std::vector<glmm::Cluster> data;
    for (int g = 0; g < 2; ++g)
      for (int c = 0; c < 300; ++c) {
        sim::Rng rng(77, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g),
                     static_cast<std::uint64_t>(c));
        const double m = 10 + std::floor(50 * rng.uniform());
        const double eta = b0 + (g == 0 ? b1 : 0) + sigma * rng.normal();
        const double s = static_cast<double>(rng.binomial(static_cast<long>(m), 1 / (1 + std::exp(-eta))));
        data.push_back({g == 0 ? 1.0 : 0.0, m, s, 1});
      }
    const glmm::Fit f = glmm::fit(data, b0, b1, sigma);
    CHECK(f.ok);
    sum_b1 += f.b1;
    sum_sigma += f.sigma;
  }
  // the averages have standard errors of about 0.075 (β1) and 0.04 (σ)
  CHECK_CLOSE(sum_b1 / reps, b1, 0, 0.3);
  CHECK_CLOSE(sum_sigma / reps, sigma, 0, 0.16);
}

TEST_CASE("glmm: without overdispersion the fit is the logistic regression") {
  // every cluster of a group has the same proportion: σ̂ = 0 and β̂ are the groups' logits
  std::vector<glmm::Cluster> data = {{1, 20, 10, 5}, {0, 20, 6, 5}};
  const glmm::Fit f = glmm::fit(data, 0, 0, 0.5);
  CHECK(f.ok);
  CHECK(f.boundary);
  CHECK_CLOSE(f.b0, logit(0.3), 1e-7, 1e-9);
  CHECK_CLOSE(f.b1, logit(0.5) - logit(0.3), 1e-7, 1e-9);
  CHECK_CLOSE(f.se_b1, std::sqrt(1 / (100 * 0.25) + 1 / (100 * 0.21)), 1e-6, 0);
}

TEST_CASE("glmm: clusters of one subject are fitted by logistic regression") {
  // σ is not identified: β and σ trade off along a ridge of maxima, and the Wald test of β1 took its standard error
  // from that ridge (a simulated power of 1.5% for 96%).  σ is held at 0, from any starting value.
  const std::vector<glmm::Cluster> data = {{1, 1, 1, 30}, {1, 1, 0, 70}, {0, 1, 1, 20}, {0, 1, 0, 80}};
  for (double start : {0.0, 0.8, 3.0}) {
    const glmm::Fit f = glmm::fit(data, -1, 0.3, start);
    CHECK(f.ok);
    CHECK(f.boundary);
    CHECK_EQ(f.sigma, 0.0);
    CHECK_CLOSE(f.b0, logit(0.2), 1e-9, 1e-12);
    CHECK_CLOSE(f.b1, logit(0.3) - logit(0.2), 1e-9, 1e-12);
    CHECK_CLOSE(f.se_b1, std::sqrt(1 / 30.0 + 1 / 70.0 + 1 / 20.0 + 1 / 80.0), 1e-9, 0);
  }
  const json::Value r = first_row(
      R"({"procedure":"sim.glmm.cluster","solve_for":"power","options":{"power_curve":false},)"
      R"("inputs":{"p1":"0.3","p2":"0.2","m":"1","k1":"500","icc":"0","test":"z","reps":"400"}})");
  CHECK_CLOSE(num(r, "power"), num(r, "approx_power"), 0, 4 * std::sqrt(0.96 * 0.04 / 400));
  CHECK_EQ(num(r, "icc2"), 0.0);  // which the quadrature left at 7e-16
  CHECK(r.get("_statement").as_string().find("500 and 500 clusters of 1 subject achieve") == 0);
}

TEST_CASE("glmm: the fit is a stationary point of the adaptive-quadrature likelihood") {
  std::vector<glmm::Cluster> data;
  sim::Rng rng(8);
  for (int g = 0; g < 2; ++g)
    for (int c = 0; c < 18; ++c) {
      const double m = 10 + std::floor(30 * rng.uniform());
      const double eta = -1 + 0.5 * (g == 0) + 0.7 * rng.normal();
      const double s = static_cast<double>(rng.binomial(static_cast<long>(m), 1 / (1 + std::exp(-eta))));
      data.push_back({g == 0 ? 1.0 : 0.0, m, s, 1});
    }
  const glmm::Fit f = glmm::fit(data, -1, 0.5, 0.7);
  CHECK(f.ok);
  CHECK(!f.boundary);
  const double th[3] = {f.b0, f.b1, f.sigma}, h = 1e-5;
  for (int i = 0; i < 3; ++i) {
    double up[3] = {th[0], th[1], th[2]}, down[3] = {th[0], th[1], th[2]};
    up[i] += h;
    down[i] -= h;
    const double grad =
        (glmm::loglik(data, up[0], up[1], up[2]) - glmm::loglik(data, down[0], down[1], down[2])) / (2 * h);
    CHECK_CLOSE(grad, 0, 0, 2e-6);
  }
  CHECK_CLOSE(glmm::loglik(data, f.b0, f.b1, f.sigma), f.loglik, 1e-11, 1e-9);
}

TEST_CASE("glmm: fits converge for very large clusters, single-subject clusters and a large SD") {
  // each cluster's integrand is then extremely narrow, flat in σ, or far from Gaussian; the fit must
  // still reach a stationary point of the quadrature likelihood
  struct Design {
    double m, cv, sigma;
    int k;
  };
  for (const Design& dsg : {Design{1e6, 0.5, 0.6, 20}, Design{1, 0, 0.8, 40}, Design{30, 0.4, 3.7, 60}}) {
    const procs::simc::ClusterSizes sizes(dsg.m, dsg.cv);
    for (int rep = 0; rep < 10; ++rep) {
      std::vector<glmm::Cluster> data;
      for (int g = 0; g < 2; ++g)
        for (int c = 0; c < dsg.k; ++c) {
          sim::Rng rng(91, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g),
                       static_cast<std::uint64_t>(c));
          const long m = sizes.draw(rng);
          const double eta = -1.4 + (g == 0 ? 0.5 : 0) + dsg.sigma * rng.normal();
          data.push_back({g == 0 ? 1.0 : 0.0, static_cast<double>(m),
                          static_cast<double>(rng.binomial(m, 1 / (1 + std::exp(-eta)))), 1});
        }
      const glmm::Fit f = glmm::fit(data, -1.4, 0.5, dsg.sigma);
      CHECK(f.ok);
      if (!f.ok) continue;
      CHECK(f.iterations < 20);
      const double th[3] = {f.b0, f.b1, f.sigma};
      for (int i = 0; i < 3; ++i) {  // central differences; the step suits the size of the log-likelihood
        const double h = dsg.m > 1e3 ? 1e-4 : 1e-5;
        double up[3] = {th[0], th[1], th[2]}, down[3] = {th[0], th[1], th[2]};
        up[i] += h;
        down[i] -= h;
        const double grad =
            (glmm::loglik(data, up[0], up[1], up[2]) - glmm::loglik(data, down[0], down[1], down[2])) / (2 * h);
        CHECK_CLOSE(grad, 0, 0, 1e-9 * std::fabs(f.loglik) / h + 1e-5);
      }
    }
  }
}

TEST_CASE("glmm: a large log odds ratio is estimated; only the between-cluster SD is held to a bound") {
  // 50 + 50 clusters of 10⁷ with 2 to 4 subjects without the outcome (group 1) or with it (group 2): neither
  // separated nor all alike, with a log odds ratio of 30.04, the logistic regression's (σ̂ = 0).  A guard on the
  // estimates (|β1| > 30 or |β0| > 50) failed such fits as run off.
  std::vector<glmm::Cluster> data;
  const double m = 1e7, few1[5] = {2, 3, 4, 3, 3}, few2[5] = {3, 2, 3, 4, 3};
  double s1 = 0, s2 = 0;
  for (int i = 0; i < 5; ++i) {
    data.push_back({1, m, m - few1[i], 10});
    data.push_back({0, m, few2[i], 10});
    s1 += 10 * (m - few1[i]);
    s2 += 10 * few2[i];
  }
  const double n = 50 * m, b0 = logit(3e-7);
  const glmm::Fit f = glmm::fit(data, b0, -2 * b0, 0);
  CHECK(f.ok);
  CHECK_CLOSE(f.b1, std::log(s1 / (n - s1)) - std::log(s2 / (n - s2)), 1e-9, 0);
  CHECK_CLOSE(f.se_b1, std::sqrt(1 / s1 + 1 / (n - s1) + 1 / s2 + 1 / (n - s2)), 1e-6, 0);
  CHECK(f.b1 > 30);
  // where the subjects of every cluster agree, σ has no finite estimate: a fit stops below SIGMA_MAX or fails
  const std::vector<glmm::Cluster> agree = {{1, 1, 1, 60}, {1, 1, 0, 40}, {1, 2, 2, 3},
                                            {0, 1, 1, 30}, {0, 1, 0, 70}, {0, 2, 0, 4}};
  for (double start : {0.0, 0.5, 3.0}) {
    const glmm::Fit g = glmm::fit(agree, -1, 1, start);
    CHECK(g.ok ? g.sigma <= glmm::SIGMA_MAX : g.failure != glmm::Failure::None);
  }
  CHECK(glmm::fit(agree, -1, 1, 0.5).failure == glmm::Failure::RunOff);
  // in the procedure: P1 0.9999997, P2 3e-7 with clusters of 10⁷ and ICC 0 gave 30% (70 of 100 studies "grew
  // without bound"), P1 0.9999999, P2 1e-7 no power at all
  const json::Value r = first_row(
      R"({"procedure":"sim.glmm.cluster","solve_for":"power","options":{"power_curve":false},)"
      R"("inputs":{"p1":"0.9999997","p2":"3e-7","m":"1e7","k1":"10","icc":"0","reps":"100"}})");
  CHECK_EQ(num(r, "failures"), 0.0);
  CHECK_EQ(num(r, "power"), 1.0);
}

TEST_CASE("glmm: the standard error comes from the Hessian of the quadrature likelihood") {
  // at a latent ICC of 0.9 the integrands are far from Gaussian, and holding the nodes fixed for the
  // Hessian understated the standard error by about a fifth; compare with a central-difference Hessian
  const double sigma = 5.44, b0 = glmm::beta_for_marginal(0.2, sigma);
  const procs::simc::ClusterSizes sizes(30, 0.4);
  std::vector<glmm::Cluster> data;
  for (int g = 0; g < 2; ++g)
    for (int c = 0; c < 155; ++c) {
      sim::Rng rng(4242, 3, static_cast<std::uint64_t>(g + 1), static_cast<std::uint64_t>(c));
      const long m = sizes.draw(rng);
      const double eta = b0 + sigma * rng.normal();
      data.push_back({g == 0 ? 1.0 : 0.0, static_cast<double>(m),
                      static_cast<double>(rng.binomial(m, 1 / (1 + std::exp(-eta)))), 1});
    }
  const glmm::Fit f = glmm::fit(data, b0, 0, sigma);
  CHECK(f.ok);
  const double th[3] = {f.b0, f.b1, f.sigma}, h = 1e-3;
  double info[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      auto at = [&](double di, double dj) {
        double t[3] = {th[0], th[1], th[2]};
        t[i] += di;
        t[j] += dj;
        return glmm::loglik(data, t[0], t[1], t[2]);
      };
      info[i][j] = -(at(h, h) - at(h, -h) - at(-h, h) + at(-h, -h)) / (4 * h * h);
    }
  // (info⁻¹)₁₁ by cofactors
  const double det = info[0][0] * (info[1][1] * info[2][2] - info[1][2] * info[2][1]) -
                     info[0][1] * (info[1][0] * info[2][2] - info[1][2] * info[2][0]) +
                     info[0][2] * (info[1][0] * info[2][1] - info[1][1] * info[2][0]);
  const double v11 = (info[0][0] * info[2][2] - info[0][2] * info[2][0]) / det;
  CHECK_CLOSE(f.se_b1, std::sqrt(v11), 2e-3, 0);
}

TEST_CASE("glmm: population-averaged proportions and ICC scales convert exactly") {
  for (double p : {0.03, 0.3, 0.8})
    for (double sigma : {0.0, 0.5, 2.0})
      CHECK_CLOSE(glmm::marginal_mean(glmm::beta_for_marginal(p, sigma), sigma), p, 1e-11, 1e-13);
  // below what an intercept of −60 gives (8.8e-27 at σ = 0, 5.7e-19 at σ = 6) β was held at −60
  for (double p : {1e-30, 1e-27, 1e-20})
    for (double sigma : {0.0, 1.0, 6.0})
      CHECK_CLOSE(glmm::marginal_mean(glmm::beta_for_marginal(p, sigma), sigma), p, 1e-11, 0);
  const double sigma = glmm::sigma_for_icc(0.05, [](double s) { return glmm::beta_for_marginal(0.2, s); });
  const double b = glmm::beta_for_marginal(0.2, sigma);
  CHECK_CLOSE(glmm::between_variance(b, sigma) / (0.2 * 0.8), 0.05, 1e-9, 0);
  // the procedure reproduces the inputs on the chosen scales
  const json::Value r = first_row(
      R"({"procedure":"sim.glmm.cluster","solve_for":"power","options":{"approximate":true},)"
      R"("inputs":{"k1":"10","m":"25","p1":"0.35","p2":"0.25","icc":"0.04"}})");
  CHECK_CLOSE(num(r, "pm1"), 0.35, 1e-10, 0);
  CHECK_CLOSE(num(r, "pm2"), 0.25, 1e-10, 0);
  CHECK_CLOSE(num(r, "icc2"), 0.04, 1e-8, 0);
  const json::Value c = first_row(
      R"({"procedure":"sim.glmm.cluster","solve_for":"power","options":{"approximate":true},)"
      R"("inputs":{"k1":"10","m":"25","p1":"0.35","p2":"0.25","icc":"0.1","pscale":"conditional",)"
      R"("iccscale":"latent"}})");
  CHECK_CLOSE(num(c, "or"), (0.35 / 0.65) / (0.25 / 0.75), 1e-12, 0);
  CHECK_CLOSE(num(c, "icc_latent"), 0.1, 1e-12, 0);
  CHECK(num(c, "pm1") < 0.35 + 0.05 && num(c, "pm1") != 0.35);  // averaging over clusters moves the proportion
}

TEST_CASE("glmm: simulated type I error and power") {
  // with 30 clusters per arm the maximum-likelihood Wald t-test is close to nominal
  // (with few clusters it rejects somewhat too often, as maximum likelihood underestimates σ)
  const json::Value h0 = first_row(
      R"({"procedure":"sim.glmm.cluster","solve_for":"power",)"
      R"("inputs":{"k1":"30","m":"30","p1":"0.2","p2":"0.2","icc":"0.05","reps":"10000","seed":"5"}})");
  CHECK_CLOSE(num(h0, "power"), 0.05, 0, 0.0075);
  CHECK_EQ(num(h0, "failures"), 0.0);
  // with many clusters the simulation agrees with the analytic approximation
  const json::Value h1 = first_row(
      R"({"procedure":"sim.glmm.cluster","solve_for":"power",)"
      R"("inputs":{"k1":"30","m":"40","cv":"0.3","p1":"0.27","p2":"0.2","icc":"0.02","reps":"4000","seed":"6"}})");
  CHECK_CLOSE(num(h1, "power"), num(h1, "approx_power"), 0, 0.035);
  CHECK(num(h1, "mean_logor") > 0);
}

TEST_CASE("glmm: reproducible across threads, and the search is minimal") {
  const std::string body =
      R"("procedure":"sim.glmm.cluster","solve_for":"power",)"
      R"("inputs":{"k1":"12","m":"20","cv":"0.6","p1":"0.35","p2":"0.2","reps":"300","seed":"41"})";
  const json::Value one = first_row("{" + body + R"(,"options":{"threads":1}})");
  const json::Value four = first_row("{" + body + R"(,"options":{"threads":4}})");
  CHECK_EQ(num(one, "power"), num(four, "power"));
  CHECK_EQ(num(one, "mean_logor"), num(four, "mean_logor"));
  const json::Value r = first_row(
      R"({"procedure":"sim.glmm.cluster","solve_for":"k",)"
      R"("inputs":{"power":"0.8","m":"20","p1":"0.35","p2":"0.2","reps":"600","seed":"3"}})");
  const double k1 = num(r, "k1");
  CHECK(num(r, "power") >= 0.8);
  const json::Value below = first_row(
      R"({"procedure":"sim.glmm.cluster","solve_for":"power","inputs":{"k1":")" + text::format_g(k1 - 1, 10) +
      R"(","m":"20","p1":"0.35","p2":"0.2","reps":"600","seed":"3"}})");
  CHECK(num(below, "power") < 0.8);
}
