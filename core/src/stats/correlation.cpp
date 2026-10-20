// Exact distribution of the Pearson correlation coefficient r of n pairs from a
// bivariate normal distribution with correlation ρ (Fisher 1915; Hotelling
// 1953):
//
//   f(r) = (n − 2) Γ(n − 1) (1 − ρ²)^{(n−1)/2} (1 − r²)^{(n−4)/2}
//          / (√(2π) Γ(n − ½) (1 − ρr)^{n−3/2}) · ₂F₁(½, ½; n − ½; (1 + ρr)/2).
//
// Probabilities integrate the density of Fisher's z = atanh(r),
// f(tanh z)·sech²z, which is smooth, free of the endpoint singularities of f
// and close to a normal density with standard deviation 1/√(n − 3).  The
// smaller tail is always integrated directly, so tail probabilities keep
// their relative accuracy.
#include <cmath>

#include "numerics.hpp"
#include "stats.hpp"

namespace ssc::stats {

namespace {

// ₂F₁(½, ½; c; x) as a power series (|x| < 1, c not zero or a negative integer).
double f21_half_series(double c, double x) {
  double s = 1, t = 1;
  for (int k = 0; k < 1000000; ++k) {
    t *= x * (k + 0.5) * (k + 0.5) / ((c + k) * (k + 1));
    s += t;
    if (std::fabs(t) <= 1e-17 * std::fabs(s) && c + k > 0) break;
  }
  return s;
}

// ₂F₁(½, ½; c; w) for 0 ≤ w < 1 and c > 2.  Near w = 1 the terms of the power
// series fall off only like k^−c, which takes tens of thousands of terms for
// small c (few pairs), so there the connection formula to 1 − w is used
// (Abramowitz & Stegun 15.3.6, with Euler's transformation of its second
// term), both of whose series converge geometrically:
//   F(½,½;c;w) = Γ(c)Γ(c−1)/Γ(c−½)² · F(½,½;2−c;1−w) + ((1−w)/w)^{c−1} / sin(πc) · F(½,½;c;1−w).
// It needs c − 1 not to be an integer; c = n − ½ for whole n.
double f21_half(double c, double w) {
  const double frac = c - std::floor(c);
  if (w <= 0.75 || c >= 30 || frac < 1e-6 || frac > 1 - 1e-6) return f21_half_series(c, w);
  const double a1 = std::exp(lgamma_diff(c - 0.5, 0.5) - lgamma_diff(c - 1, 0.5));  // Γ(c)Γ(c−1)/Γ(c−½)²
  return a1 * f21_half_series(2 - c, 1 - w) +
         std::pow((1 - w) / w, c - 1) / std::sin(PI * c) * f21_half_series(c, 1 - w);
}

// ln cosh x, accurate for small and large |x|
double log_cosh(double x) {
  const double a = std::fabs(x);
  if (a < 1) {
    const double h = std::sinh(0.5 * a);
    return std::log1p(2 * h * h);
  }
  return a + std::log1p(std::exp(-2 * a)) - 0.693147180559945309417232121458;
}

// Density of z = atanh(r).  With ρ = tanh ζ, 1 − ρr = cosh(z − ζ)/(cosh ζ cosh z),
// so that
//   ln g(z) = K − ½ ln cosh ζ + ½ ln cosh z − (n − 3/2) ln cosh(z − ζ) + ln ₂F₁,
//   K = ln(n − 2) + ln Γ(n − 1) − ln Γ(n − ½) − ln √(2π),
// where every term is O(1) near the bulk: no cancellation of large terms even
// for very large n.
struct ZDensity {
  double n, zeta, lcz, k0;
  ZDensity(double n_, double rho) : n(n_), zeta(std::atanh(rho)) {
    lcz = log_cosh(zeta);
    k0 = std::log(n - 2) - lgamma_diff(n - 1, 0.5) - LN_SQRT_2PI - 0.5 * lcz;
  }
  double operator()(double z) const {
    const double lcd = log_cosh(z - zeta), lc = log_cosh(z);
    const double omrr = std::exp(lcd - lcz - lc);  // 1 − ρr
    return std::exp(k0 + 0.5 * lc - (n - 1.5) * lcd + std::log(f21_half(n - 0.5, 1 - 0.5 * omrr)));
  }
};

}  // namespace

double dcorr(double r, double n, double rho) {
  if (std::isnan(r) || !(n > 2) || !(std::fabs(rho) < 1)) return NaN;
  if (r <= -1 || r >= 1) return 0;
  return ZDensity(n, rho)(std::atanh(r)) / ((1 - r) * (1 + r));
}

double pcorr(double x, double n, double rho, bool lower) {
  if (std::isnan(x) || std::isnan(rho) || !(n > 2) || !(std::fabs(rho) <= 1)) return NaN;
  if (x <= -1) return lower ? 0.0 : 1.0;
  if (x >= 1) return lower ? 1.0 : 0.0;
  if (std::fabs(rho) == 1) return (x >= rho) == lower ? 1.0 : 0.0;  // r = ρ almost surely
  const ZDensity g(n, rho);
  const double zx = std::atanh(x), zm = std::atanh(rho);
  const double s = 1 / std::sqrt(std::max(n - 3, 1.0));
  // integrate from z(x) outwards, through the bulk if it lies on that side,
  // until the density is negligible
  const double end = lower ? std::min(zx, zm) - 45 * s : std::max(zx, zm) + 45 * s;
  const double width = 1.5 * s;
  double total = 0, a = zx;
  while (lower ? a > end : a < end) {
    const double b = lower ? std::max(end, a - width) : std::min(end, a + width);
    total += numerics::integrate(g, std::min(a, b), std::max(a, b), 1e-300, 1e-14, 20);
    a = b;
  }
  return clamp01(total);
}

}  // namespace ssc::stats
