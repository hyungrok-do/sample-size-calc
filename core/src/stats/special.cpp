// Elementary special functions.
//
// Many probabilities in power analysis reduce to x^a(1−x)^b/B(a,b) or
// λ^x e^{−λ}/Γ(x+1).  Computing these via differences of lnΓ values loses
// accuracy for large arguments, so (following C. Loader, "Fast and accurate
// computation of binomial probabilities", 2000) we write them through the
// Stirling remainder δ(x) and the deviance term bd0(x, m), both of which are
// available to full relative precision.
#include "stats.hpp"

namespace ssc::stats {

namespace {
// B_{2k} / (2k(2k−1)), k = 1..8
constexpr double S1 = 1.0 / 12.0;
constexpr double S2 = -1.0 / 360.0;
constexpr double S3 = 1.0 / 1260.0;
constexpr double S4 = -1.0 / 1680.0;
constexpr double S5 = 1.0 / 1188.0;
constexpr double S6 = -691.0 / 360360.0;
constexpr double S7 = 1.0 / 156.0;
constexpr double S8 = -3617.0 / 122400.0;

inline double stirling_series(double x) {
  const double t = 1.0 / x;
  const double t2 = t * t;
  return t * (S1 + t2 * (S2 + t2 * (S3 + t2 * (S4 + t2 * (S5 + t2 * (S6 + t2 * (S7 + t2 * S8)))))));
}
}  // namespace

double stirlerr(double x) {
  if (x >= 10.0) return stirling_series(x);
  // δ(x) = δ(x+1) + (x + ½) ln(1 + 1/x) − 1
  double s = 0.0;
  double y = x;
  while (y < 10.0) {
    s += (y + 0.5) * std::log1p(1.0 / y) - 1.0;
    y += 1.0;
  }
  return s + stirling_series(y);
}

double lgammafn(double x) {
  if (std::isnan(x)) return NaN;
  if (x <= 0) {
    if (x == std::floor(x)) return INF;
    return std::log(PI / std::fabs(std::sin(PI * x))) - lgammafn(1.0 - x);
  }
  if (std::isinf(x)) return INF;
  if (x < 1e-300) return -std::log(x);
  return (x - 0.5) * std::log(x) - x + LN_SQRT_2PI + stirlerr(x);
}

double lgamma_diff(double b, double a) {
  if (b >= 8.0) {
    return (b - 0.5) * std::log1p(a / b) + a * std::log(a + b) - a + stirlerr(a + b) - stirlerr(b);
  }
  return lgammafn(a + b) - lgammafn(b);
}

double lbeta(double a, double b) {
  const double p = std::fmin(a, b);
  const double q = std::fmax(a, b);
  if (p < 0) return NaN;
  if (p == 0) return INF;
  if (p >= 10.0) {
    const double s = p + q;
    return LN_SQRT_2PI - 0.5 * std::log(q) + (p - 0.5) * std::log(p / s) + q * std::log1p(-p / s) +
           stirlerr(p) + stirlerr(q) - stirlerr(s);
  }
  if (q >= 10.0) return lgammafn(p) - lgamma_diff(q, p);
  return lgammafn(p) + lgammafn(q) - lgammafn(p + q);
}

double bd0(double x, double m) {
  if (x == 0) return m;
  if (!std::isfinite(x) || !std::isfinite(m) || m == 0) return INF;
  const double d = x - m;
  if (std::fabs(d) < 0.1 * (x + m)) {
    const double v = d / (x + m);
    double s = d * v;
    if (std::fabs(s) < DBL_MIN_) return s;
    double ej = 2.0 * x * v;
    const double v2 = v * v;
    for (int j = 1; j < 1000; ++j) {
      ej *= v2;
      const double s1 = s + ej / (2 * j + 1);
      if (s1 == s) return s1;
      s = s1;
    }
    return s;
  }
  return x * std::log(x / m) + m - x;
}

double dpois_raw(double x, double lambda) {
  if (lambda == 0) return x == 0 ? 1.0 : 0.0;
  if (!std::isfinite(lambda) || x < 0) return 0.0;
  if (x == 0) return std::exp(-lambda);
  if (x <= lambda * DBL_MIN_) return std::exp(-lambda);
  if (lambda < x * DBL_MIN_) return std::exp(-lambda + x * std::log(lambda) - lgammafn(x + 1));
  return std::exp(-stirlerr(x) - bd0(x, lambda)) / std::sqrt(2.0 * PI * x);
}

double dbinom_raw(double x, double n, double p, double q) {
  if (p == 0) return x == 0 ? 1.0 : 0.0;
  if (q == 0) return x == n ? 1.0 : 0.0;
  if (x == 0) {
    if (n == 0) return 1.0;
    return std::exp(p < 0.1 ? -bd0(n, n * q) - n * p : n * std::log(q));
  }
  if (x == n) return std::exp(q < 0.1 ? -bd0(n, n * p) - n * q : n * std::log(p));
  if (x < 0 || x > n) return 0.0;
  const double lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);
  const double lf = LN_2PI + std::log(x) + std::log1p(-x / n);
  return std::exp(lc - 0.5 * lf);
}

double beta_prefactor(double a, double b, double x, double y) {
  if (x <= 0 || y <= 0) return 0.0;
  // x^a y^b / B(a,b) = √(ab/(2πn)) exp(δ(n) − δ(a) − δ(b) − bd0(a, nx) − bd0(b, ny))
  const double n = a + b;
  const double e = stirlerr(n) - stirlerr(a) - stirlerr(b) - bd0(a, n * x) - bd0(b, n * y);
  return std::sqrt(a * b / (2.0 * PI * n)) * std::exp(e);
}

double lchoose(double n, double k) { return -std::log(n + 1.0) - lbeta(n - k + 1.0, k + 1.0); }

}  // namespace ssc::stats
