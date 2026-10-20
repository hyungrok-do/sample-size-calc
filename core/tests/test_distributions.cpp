// Cross-checks of distributions computed by two unrelated methods, and of the
// large-parameter paths (Temme's expansion, Poisson-mixture integrals).
#include "procedures/integrals.hpp"
#include "stats/stats.hpp"
#include "testing.hpp"

using namespace ssc;

namespace {

// Sample correlation through the slope t-statistic: conditionally on the X
// values, t = r√(n−2)/√(1−r²) is noncentral t with noncentrality ρ√W/√(1−ρ²),
// W ~ χ²_{n−1} — a χ²-mixture of noncentral t distributions.
double pcorr_mixture(double r0, double n, double rho, bool lower) {
  const double df = n - 2;
  const double t0 = r0 * std::sqrt(df) / std::sqrt((1 - r0) * (1 + r0));
  if (rho == 0) return stats::pt(t0, df, lower);
  const double k = rho * std::sqrt(n - 1) / std::sqrt((1 - rho) * (1 + rho));
  return procs::expect_sd_ratio([&](double u) { return stats::pnt(t0, df, k * u, lower); }, n - 1, INFINITY, 1e-15,
                                1e-12);
}

}  // namespace

TEST_CASE("distributions: correlation density method vs noncentral-t mixture") {
  for (double n : {4.0, 7.0, 20.0, 60.0, 250.0}) {
    for (double rho : {-0.85, -0.3, 0.0, 0.25, 0.6, 0.93}) {
      const double sd = (1 - rho * rho) / std::sqrt(n - 1);
      for (double k : {-4.0, -2.0, -0.5, 0.0, 1.0, 2.5, 4.0}) {
        const double x = std::max(-0.999, std::min(0.999, rho + k * sd));
        for (bool lower : {true, false}) {
          const double ref = pcorr_mixture(x, n, rho, lower);
          if (ref < 1e-250) continue;
          CHECK_CLOSE(stats::pcorr(x, n, rho, lower), ref, 2e-9, 0);
        }
        CHECK_CLOSE(stats::pcorr(x, n, rho, true) + stats::pcorr(x, n, rho, false), 1.0, 0, 1e-14);
      }
    }
  }
}

TEST_CASE("distributions: correlation under rho = 0 is a t-distribution") {
  for (double n : {3.0, 5.0, 12.0, 40.0, 500.0, 1e5}) {
    for (double x : {-0.9, -0.4, -0.05, 0.02, 0.3, 0.75}) {
      const double t = x * std::sqrt(n - 2) / std::sqrt((1 - x) * (1 + x));
      for (bool lower : {true, false}) {
        const double ref = stats::pt(t, n - 2, lower);
        if (ref < 1e-250) continue;
        CHECK_CLOSE(stats::pcorr(x, n, 0.0, lower), ref, 1e-12, 0);
      }
    }
  }
}

TEST_CASE("distributions: incomplete gamma for large shapes (Temme) vs series") {
  // reference: the defining series P(a, x) = x^a e^{−x}/Γ(a + 1) Σ x^k/((a+1)…(a+k)),
  // with compensated (Kahan) summation so that it is accurate on every compiler
  for (double a : {120.0, 900.0, 25000.0}) {
    for (double k : {-6.0, -2.0, -0.3, 0.0, 0.4, 3.0}) {
      const double x = a + k * std::sqrt(a);
      double s = 1, c = 0, t = 1;
      for (int i = 1; i < 2000000; ++i) {
        t *= x / (a + i);
        const double y = t - c, u = s + y;
        c = (u - s) - y;
        s = u;
        if (t < s * 1e-19) break;
      }
      // x^a e^{−x}/Γ(a + 1) from the saddle-point form (accurate for large a)
      const double ref = s * stats::dpois_raw(a, x);
      CHECK_CLOSE(stats::pgamma_raw(x, a, true), ref, 3e-13, 0);
      CHECK_CLOSE(stats::pgamma_raw(x, a, true) + stats::pgamma_raw(x, a, false), 1.0, 0, 4e-15);
    }
  }
}

TEST_CASE("distributions: huge noncentrality stays accurate and fast") {
  // P(χ'² ≤ df + λ) → 1/2 + skewness term; the upper tail of the F test is ≈ 1.
  for (double lam : {1e6, 1e9, 1e12, 1e15}) {
    const double p = stats::pnchisq(3 + lam, 3, lam, true);
    CHECK(p > 0.5 && p < 0.5 + 0.3 / std::sqrt(lam));
    CHECK_CLOSE(stats::pnf(3.2, 3, 96, lam, false), 1.0, 0, 1e-14);
    const double d = std::sqrt(2 * lam);
    CHECK_CLOSE(stats::pnt(d, 50, d, true), stats::pchisq(50, 50, false), 2e-3, 0);
  }
}
