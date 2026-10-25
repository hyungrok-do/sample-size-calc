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

TEST_CASE("distributions: a small beta tail beside one large shape, from TOMS 708's expansion") {
  // Where the series' complement is too close to 1 (x above ½ and the lower tail below 0.01, or the mirror): BGRAT
  // for b ≤ 1 and the series for x ≤ 0.7, raised to b by positive terms.  References: the series of the other tail in
  // 60 to 400 digits, 1 less it (mpmath); R's pbeta agrees to 2e-15 and less.
  struct Case {
    double x, a, b;
    bool lower;
    double ref;
  };
  for (const Case& c : {Case{0.9983406602323119, 39900.33457501, 28.405605073134232, true, 5.3941505331920082e-8},
                        Case{0.9925807368451423, 6500.026766827606, 0.3842282849380524, true, 3.7330702638813437e-23},
                        Case{0.6, 50.5, 15, true, 1.4817443899652601e-3},
                        Case{0.999998, 4500000.5, 1, true, 1.2340857002292161e-4},
                        Case{0.999998, 4500000.5, 1.5, true, 4.3984526648444234e-4},
                        Case{2e-6, 1.5, 4500000.5, false, 4.3984526638384715e-4},
                        Case{0.25, 2.5, 120, false, 1.3249309816841355e-13}})
    CHECK_CLOSE(stats::pbeta(c.x, c.a, c.b, c.lower), c.ref, 1e-14, 0);
  // the t distribution's tail, ½ I_x(ν/2, ½) with x = ν/(ν + t²)
  CHECK_CLOSE(stats::pt(5, 30, false), 1.1648342733503898e-5, 1e-14, 0);
  CHECK_CLOSE(stats::pt(3, 1000, false), 1.3833545221190962e-3, 1e-14, 0);
  CHECK_CLOSE(stats::pt(4.5, 1e6, false), 3.3980552524769493e-6, 1e-14, 0);
}

TEST_CASE("distributions: the upper incomplete gamma at a small shape, and the beta tail beside one") {
  // Q(a, x) for a < ½ and x < a + 1 was 1 − P, about a·E1(x) known only to 1e-16 (2.3e-5 relative at a = 1e-10);
  // BGRAT takes it for a tiny second shape b, and pbeta(0.99, 30, 1e-15) was 3.3e-16 for 9.1e-16.  References:
  // mpmath's gammainc, and the beta tail integrated in 60 digits.
  struct Case {
    double x, a, ref;
  };
  for (const Case& c : {Case{1.0, 1e-10, 2.1938393441796778575e-11}, Case{0.5, 1e-3, 5.6006665647074988868e-4},
                        Case{1.2, 0.3, 0.063394719569259104781}, Case{1e-300, 1e-30, 6.9019831223331222984e-28},
                        Case{0.5615, 1e-6, 4.9149342226249250718e-7}, Case{1.4, 0.49, 0.091717618248318038954}})
    CHECK_CLOSE(stats::pgamma_raw(c.x, c.a, false), c.ref, 1e-14, 0);
  CHECK_CLOSE(stats::pbeta(0.99, 30, 1e-15, true), 9.1437730900816014257e-16, 1e-14, 0);
  CHECK_CLOSE(stats::pbeta(0.01, 1e-15, 30, false), 9.1437730900816014257e-16, 1e-14, 0);
  CHECK_CLOSE(stats::pf(5, 2e-12, 100, false), 2.5454400870275017447e-11, 1e-14, 0);
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
