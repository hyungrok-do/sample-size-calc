// Distribution functions vs. reference values computed independently of the engine
// (each fixture's first line says with what).
// Fixture format (CSV, '#' comments):  fn,arg1 arg2 ...,expected,rtol,atol
#include <fstream>
#include <map>

#include "stats/dispatch.hpp"
#include "stats/numerics.hpp"
#include "stats/stats.hpp"
#include "testing.hpp"
#include "util/numtext.hpp"

namespace {

struct ErrStat {
  int n = 0;
  int bad = 0;
  double max_rel = 0;
};

void run_fixture(const std::string& file) {
  std::ifstream in(std::string(SSC_FIXTURE_DIR) + "/" + file);
  CHECK(in.good());
  if (!in.good()) return;
  std::map<std::string, ErrStat> stats;
  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> f;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) f.push_back(tok);
    if (f.size() < 5) continue;
    std::vector<double> args;
    std::stringstream as(f[1]);
    double v;
    while (as >> v) args.push_back(v);
    double expected = 0, rtol = 0, atol = 0;
    // (std::stod throws on subnormal values; parse_double does not)
    if (!ssc::text::parse_double(f[2], expected) || !ssc::text::parse_double(f[3], rtol) ||
        !ssc::text::parse_double(f[4], atol)) {
      ssc::test::report_failure("unparsable line: " + line, file.c_str(), lineno);
      continue;
    }
    double got = 0;
    ssc::test::count_check();
    auto& st = stats[f[0]];
    ++st.n;
    if (!ssc::stats::eval(f[0], args.data(), static_cast<int>(args.size()), got)) {
      ssc::test::report_failure("unknown function or arity: " + f[0], file.c_str(), lineno);
      continue;
    }
    if (std::isfinite(expected) && expected != 0) {
      const double rel = std::fabs(got - expected) / std::fabs(expected);
      if (rel > st.max_rel && std::fabs(got - expected) > atol) st.max_rel = rel;
    }
    if (!ssc::test::close(got, expected, rtol, atol)) {
      ++st.bad;
      std::ostringstream os;
      os.precision(17);
      os << f[0] << "(" << f[1] << ") = " << got << ", expected " << expected;
      ssc::test::report_failure(os.str(), file.c_str(), lineno);
    }
  }
  for (const auto& kv : stats)
    std::cout << "    " << kv.first << ": " << kv.second.n << " values, max rel err " << kv.second.max_rel
              << (kv.second.bad ? "  <-- failures: " + std::to_string(kv.second.bad) : std::string()) << "\n";
}

}  // namespace

TEST_CASE("stats: reference values") { run_fixture("stats_reference.csv"); }

TEST_CASE("stats: high-precision reference values") { run_fixture("stats_high_precision.csv"); }

TEST_CASE("stats: quantile round trips") {
  using namespace ssc::stats;
  const double ps[] = {1e-300, 1e-100, 1e-20, 1e-8, 1e-3, 0.025, 0.2, 0.5, 0.8, 0.975, 0.999};
  for (double p : ps) {
    // the round trip amplifies the quantile's rounding error by ~|log p|
    const double tol = p < 1e-10 ? 2e-11 : 1e-12;
    for (double df : {1.0, 2.0, 3.0, 4.5, 10.0, 47.3, 1e3, 1e6}) {
      if (p < 1e-200 && df < 3) continue;
      const double t = qt(p, df);
      CHECK_CLOSE(pt(t, df), p, tol, 0);
      const double tu = qt(p, df, false);
      CHECK_CLOSE(pt(tu, df, false), p, tol, 0);
    }
    for (double a : {0.05, 0.5, 1.0, 3.7, 50.0, 1e4}) {
      const double x = qgamma(p, a);
      if (x > 0) CHECK_CLOSE(pgamma_raw(x, a), p, tol, 0);
      const double xu = qgamma(p, a, false);
      if (std::isfinite(xu)) CHECK_CLOSE(pgamma_raw(xu, a, false), p, tol, 0);
    }
    for (auto ab : {std::make_pair(0.5, 0.5), std::make_pair(2.0, 3.0), std::make_pair(0.1, 30.0),
                    std::make_pair(400.0, 2.5), std::make_pair(1e4, 3e4)}) {
      const double x = qbeta(p, ab.first, ab.second);
      if (x > 0 && x < 1) CHECK_CLOSE(pbeta(x, ab.first, ab.second), p, 1e-11, 0);
      const auto xy = qbeta_pair(p, ab.first, ab.second, false);
      if (xy.first > 0 && xy.second > 0)
        CHECK_CLOSE(pbeta_raw(xy.first, ab.first, ab.second, false, xy.second), p, 1e-11, 0);
    }
    const double z = qnorm(p);
    CHECK_CLOSE(pnorm(z), p, tol, 0);
  }
}

TEST_CASE("stats: quadrature and root finding") {
  using namespace ssc;
  // ∫_0^π sin = 2 ; ∫_0^1 x^20 = 1/21 ; ∫ normal density over (−8, 8) ≈ 1
  CHECK_CLOSE(numerics::integrate([](double x) { return std::sin(x); }, 0, stats::PI), 2.0, 1e-13, 0);
  CHECK_CLOSE(numerics::integrate([](double x) { return std::pow(x, 20); }, 0, 1), 1.0 / 21, 1e-13, 0);
  CHECK_CLOSE(numerics::gauss_legendre_integrate([](double x) { return std::pow(x, 39); }, 0, 1, 20), 1.0 / 40,
              1e-13, 0);
  CHECK_CLOSE(numerics::integrate([](double x) { return stats::dnorm(x); }, -8, 8), 1 - 2 * stats::pnorm(-8), 1e-13, 0);
  const double r = numerics::brent([](double x) { return x * x * x - 2; }, 0, 2, 1e-15);
  CHECK_CLOSE(r, std::cbrt(2.0), 1e-14, 0);
  // Gauss–Legendre weights sum to 2 for many orders
  for (int n : {2, 5, 16, 64, 128}) {
    double s = 0;
    for (double w : stats::gauss_legendre(n).w) s += w;
    CHECK_CLOSE(s, 2.0, 1e-13, 0);
  }
}

TEST_CASE("stats: noncentral consistency identities") {
  using namespace ssc::stats;
  // lower + upper = 1 and reflection symmetry for the noncentral t
  for (double df : {1.0, 3.0, 9.5, 40.0, 500.0, 1e5}) {
    for (double ncp : {-12.0, -3.0, -0.5, 0.7, 2.0, 8.0, 45.0}) {
      for (double t : {-30.0, -2.0, -0.3, 0.0, 0.4, 1.9, 6.0, 50.0}) {
        const double lo = pnt(t, df, ncp, true), up = pnt(t, df, ncp, false);
        CHECK_CLOSE(lo + up, 1.0, 0, 2e-14);
        CHECK_CLOSE(pnt(-t, df, -ncp, false), lo, 1e-12, 1e-300);
      }
    }
  }
  for (double df : {1.0, 4.0, 30.0}) {
    for (double ncp : {0.5, 5.0, 60.0, 400.0}) {
      for (double x : {0.1, 3.0, 20.0, 90.0, 600.0}) {
        CHECK_CLOSE(pnchisq(x, df, ncp, true) + pnchisq(x, df, ncp, false), 1.0, 0, 2e-14);
      }
    }
  }
}
