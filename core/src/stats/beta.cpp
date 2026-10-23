// Regularized incomplete beta function and the beta, Student t, F and
// binomial distributions.
//
// I_x(a, b) is evaluated from the classical continued fraction (DLMF 8.17.22)
// with the modified Lentz algorithm, applied to whichever of I_x(a,b) or
// I_{1−x}(b,a) converges fastest, so the tail that is computed directly is the
// small one and neither tail suffers from 1 − (tiny) cancellation.  The
// prefactor uses beta_prefactor(), which keeps full relative precision for
// large shape parameters (large df / sample sizes).
#include "stats.hpp"
#include "numerics.hpp"

namespace ssc::stats {

namespace {
constexpr double TINY = 1e-300;

double beta_cf(double a, double b, double x) {
  const double qab = a + b, qap = a + 1, qam = a - 1;
  double c = 1;
  double d = 1 - qab * x / qap;
  if (std::fabs(d) < TINY) d = TINY;
  d = 1 / d;
  double h = d;
  // bounded before conversion: shapes beyond 1e16 (a noncentral F with an enormous noncentrality) pass any int
  const int maxIt = 20000 + static_cast<int>(std::fmin(20 * std::sqrt(std::fmax(a, b)), 1e8));
  for (int m = 1; m <= maxIt; ++m) {
    const double m2 = 2.0 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (std::fabs(d) < TINY) d = TINY;
    c = 1 + aa / c;
    if (std::fabs(c) < TINY) c = TINY;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (std::fabs(d) < TINY) d = TINY;
    c = 1 + aa / c;
    if (std::fabs(c) < TINY) c = TINY;
    d = 1 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1) <= 4e-16) break;
  }
  return h;
}

// Positive-term hypergeometric series (DLMF 8.17.8):
//   I_x(a,b) = x^a y^b / (a B(a,b)) · Σ_n (a+b)_n / (a+1)_n · x^n.
// Efficient when (a+b)·x is moderate and x ≤ ½; all terms are positive so the
// sum is accurate to a few ulps (the continued fraction needs O(1/(1−ρ))
// iterations near its switch point when the other shape parameter is small).
double beta_series(double a, double b, double x, double y) {
  const double pre = beta_prefactor(a, b, x, y) / a;
  if (pre == 0) return 0;
  double term = 1, sum = 1;
  for (int n = 0; n < 100000; ++n) {
    term *= (a + b + n) * x / (a + 1 + n);
    sum += term;
    if (term < 1e-17 * sum && (a + b + n) * x < (a + 1 + n)) break;
  }
  return pre * sum;
}

// x − ln(1 + x), accurate for small x.
double rlog1(double x) {
  if (std::fabs(x) > 0.1) return x - std::log1p(x);
  double term = x * x, sum = 0;  // x²/2 − x³/3 + x⁴/4 − …
  for (int k = 2; k < 60; ++k) {
    sum += (k % 2 ? -term : term) / k;
    term *= x;
    if (std::fabs(term) < 1e-17 * std::fabs(sum)) break;
  }
  return sum;
}

// exp(x²)·erfc(x) for x ≥ 0, from its continued fraction where erfc alone would underflow.
double erfcx(double x) {
  if (x < 2) return std::exp(x * x) * std::erfc(x);
  double t = x;  // 1/(x + ½/(x + 1/(x + 3/2/(x + …))))
  for (int k = 200; k >= 1; --k) t = x + 0.5 * k / t;
  return 1 / (1.77245385090551602730 * t);  // √π
}

// I_x(a, b) for large a and b near the mean, with λ = a − (a + b)·x ≥ 0: the uniform asymptotic expansion of
// Didonato & Morris (1992, ACM TOMS 708, BASYM), in terms of erfc.  The series and the continued fraction need
// about √a terms there and lose accuracy beyond a ≈ 1e8.
double beta_asymptotic(double a, double b, double lambda) {
  constexpr int N = 20;  // terms (even)
  constexpr double e0 = 1.12837916709551257;  // 2/√π
  constexpr double e1 = 0.353553390593273762;  // 2^(−3/2)
  constexpr double eps = 1e-15;
  double a0[N + 1] = {}, b0[N + 1] = {}, c[N + 1] = {}, d[N + 1] = {};
  const double f = a * rlog1(-lambda / a) + b * rlog1(lambda / b);
  const double t = std::exp(-f);
  if (t == 0) return 0;
  const double z0 = std::sqrt(f), z = z0 / e1 * 0.5, z2 = f + f;
  double h, r0, r1, w0;
  if (a < b) {
    h = a / b;
    r0 = 1 / (h + 1);
    r1 = (b - a) / b;
    w0 = 1 / std::sqrt(a * (h + 1));
  } else {
    h = b / a;
    r0 = 1 / (h + 1);
    r1 = (b - a) / a;
    w0 = 1 / std::sqrt(b * (h + 1));
  }
  a0[0] = r1 * 2.0 / 3;
  c[0] = -0.5 * a0[0];
  d[0] = -c[0];
  double j0 = 0.5 / e0 * erfcx(z0), j1 = e1;
  double sum = j0 + d[0] * w0 * j1;
  double sh = 1, h2 = h * h, hn = 1, w = w0, znm1 = z, zn = z2;
  for (int n = 2; n <= N; n += 2) {
    hn *= h2;
    a0[n - 1] = 2 * r0 * (h * hn + 1) / (n + 2.0);
    const int np1 = n + 1;
    sh += hn;
    a0[np1 - 1] = 2 * r1 * sh / (n + 3.0);
    for (int i = n; i <= np1; ++i) {
      const double r = -0.5 * (i + 1.0);
      b0[0] = r * a0[0];
      for (int m = 2; m <= i; ++m) {
        double bsum = 0;
        for (int j = 1; j <= m - 1; ++j) bsum += (j * r - (m - j)) * a0[j - 1] * b0[m - j - 1];
        b0[m - 1] = r * a0[m - 1] + bsum / m;
      }
      c[i - 1] = b0[i - 1] / (i + 1.0);
      double dsum = 0;
      for (int j = 1; j <= i - 1; ++j) dsum += d[i - j - 1] * c[j - 1];
      d[i - 1] = -(dsum + c[i - 1]);
    }
    j0 = e1 * znm1 + (n - 1.0) * j0;
    j1 = e1 * zn + n * j1;
    znm1 *= z2;
    zn *= z2;
    w *= w0;
    const double t0 = d[n - 1] * w * j0;
    w *= w0;
    const double t1 = d[np1 - 1] * w * j1;
    sum += t0 + t1;
    if (std::fabs(t0) + std::fabs(t1) <= eps * sum) break;
  }
  const double bcorr = stirlerr(a) + stirlerr(b) - stirlerr(a + b);
  return e0 * t * std::exp(-bcorr) * sum;
}

// P(Z ≥ z) for Z ~ Beta(c, d) by adaptive Gauss–Kronrod integration of the
// density, used when z lies beyond the bulk of a distribution concentrated
// near c/(c+d) (one shape small, the other large) and the complement computed
// by the series is too close to 1 to be subtracted safely.
double beta_upper_tail_integral(double c, double d, double z) {
  const double s = c + d;
  const double lb = lbeta(c, d);
  auto g = [&](double w) {
    const double t = z + w / s;
    if (t >= 1) return 0.0;
    return std::exp((c - 1) * std::log(t) + (d - 1) * std::log1p(-t) - lb) / s;
  };
  double wmax = std::max(80.0, 4 * (c + 20));
  if (z + wmax / s > 1) wmax = (1 - z) * s;
  const double pts[] = {0.0, 0.5, 2.0, 6.0, 15.0, 35.0, 80.0, 160.0, 320.0};
  double total = 0, prev = 0;
  for (double q : pts) {
    if (q <= prev) continue;
    const double hi = std::min(q, wmax);
    if (hi > prev) total += numerics::integrate(g, prev, hi, 0.0, 1e-15, 30);
    prev = hi;
    if (prev >= wmax) break;
  }
  if (prev < wmax) total += numerics::integrate(g, prev, wmax, 0.0, 1e-15, 30);
  return total;
}

double beta_guess(double p, double a, double b, bool lower) {
  const double pLow = lower ? p : 1 - p;
  const double pUp = lower ? 1 - p : p;
  const double lb = lbeta(a, b);
  // power-law tails: I_x ≈ x^a/(aB) near 0; 1 − I_x ≈ (1−x)^b/(bB) near 1
  const double xLowTail = std::exp((std::log(pLow) + std::log(a) + lb) / a);
  const double yUpTail = std::exp((std::log(pUp) + std::log(b) + lb) / b);
  if (a > 1 && b > 1) {
    // normal-based approximation (cf. Abramowitz & Stegun 26.5.22)
    const double z = qnorm(p, !lower);
    const double r = (z * z - 3) / 6;
    const double s = 1 / (2 * a - 1);
    const double t = 1 / (2 * b - 1);
    const double h = 2 / (s + t);
    const double w = z * std::sqrt(h + r) / h - (t - s) * (r + 5.0 / 6 - 2 / (3 * h));
    double x = a / (a + b * std::exp(2 * w));
    if (!(x > 0 && x < 1)) x = a / (a + b);
    if (pLow < 1e-8 && xLowTail < x) x = xLowTail;
    if (pUp < 1e-8 && 1 - yUpTail > x) x = 1 - yUpTail;
    return x;
  }
  const double mean = a / (a + b);
  if (xLowTail < mean) return std::fmin(std::fmax(xLowTail, 1e-300), 1 - 1e-16);
  if (1 - yUpTail > mean) return std::fmax(std::fmin(1 - yUpTail, 1 - 1e-16), 1e-300);
  return mean;
}

// Safeguarded Halley iteration on f(x) = I_x(a,b) − target (bracket [0, 1]).
double solve_beta(double p, double a, double b, bool lower, double x0) {
  double lo = 0, hi = 1;
  double x = std::fmin(std::fmax(x0, 1e-300), 1 - 1e-16);
  for (int it = 0; it < 500; ++it) {
    const double y = 1 - x;
    const double F = pbeta_raw(x, a, b, lower, y);
    const double diff = lower ? F - p : p - F;
    if (diff == 0) return x;
    if (diff > 0)
      hi = x;
    else
      lo = x;
    const double dens = dbeta(x, a, b, y);
    double xNew = NaN;
    if (dens > 0 && std::isfinite(dens)) {
      const double t = diff / dens;
      const double corr = 1 - 0.5 * t * ((a - 1) / x - (b - 1) / y);
      const double step = (corr > 0.2 && corr < 5) ? t / corr : t;
      xNew = x - step;
    }
    if (!(xNew > lo && xNew < hi)) {
      if (lo > 0 && hi / lo > 1e3)
        xNew = std::sqrt(lo * hi);  // geometric bisection across magnitudes
      else if (lo == 0 && x < 1e-3)
        xNew = x * 1e-3;
      else
        xNew = 0.5 * (lo + hi);
      if (!(xNew > lo && xNew < hi)) xNew = 0.5 * (lo + hi);
    }
    if (std::fabs(xNew - x) <= 2e-16 * xNew) return xNew;
    if (hi - lo <= 4e-16 * hi) return 0.5 * (lo + hi);
    x = xNew;
  }
  return x;
}
}  // namespace

double pbeta_raw(double x, double a, double b, bool lower, double y) {
  if (std::isnan(x) || std::isnan(a) || std::isnan(b)) return NaN;
  if (x <= 0) return lower ? 0.0 : 1.0;
  if (y <= 0) return lower ? 1.0 : 0.0;
  if (a == 0 || b == 0 || std::isinf(a) || std::isinf(b)) {
    if (a == 0 && b == 0) return NaN;
    if (a == 0 || a / b == 0) return lower ? 1.0 : 0.0;  // mass at 0
    return lower ? 0.0 : 1.0;                              // mass at 1
  }
  // Both shape parameters large, within 3% of the smaller one from the mean (beyond lies a tail far below any
  // double's resolution of the other): the asymptotic expansion, for the tail below the mean.
  if (std::fmin(a, b) >= 1e6) {
    // a − (a + b)·x; from y when a > b (as TOMS 708 writes it), which keeps its precision there
    const double lambda = a > b ? (a + b) * y - b : a - (a + b) * x;
    if (std::fabs(lambda) <= 0.03 * std::fmin(a, b)) {
      if (lambda >= 0) {
        const double I = beta_asymptotic(a, b, lambda);  // I_x(a, b)
        return lower ? I : 1 - I;
      }
      const double J = beta_asymptotic(b, a, -lambda);  // I_y(b, a) = 1 − I_x(a, b)
      return lower ? 1 - J : J;
    }
    // beyond, the tail is below e^(−4000) from 1e7 on (λ is over 95 of its standard deviations, √(ab/(a+b+1)),
    // from the mean): it rounds to 0
    if (std::fmin(a, b) >= 1e7) return (lambda > 0) == lower ? 0.0 : 1.0;
  }
  // One shape parameter astronomically larger than the other (the Poisson terms of a noncentral F with an
  // enormous noncentrality): the series and the continued fraction cannot resolve a + m there.  As b → ∞,
  // −(b + (a − 1)/2)·ln(1 − x) has the Gamma(a) distribution, with a relative error of about 16·(a/b)², under
  // 2e-17 here.
  if (std::fmax(a, b) >= 1e15 && std::fmin(a, b) <= 1e-9 * std::fmax(a, b)) {
    if (b > a) return pgamma_raw((b + (a - 1) / 2) * (x < 0.5 ? -std::log1p(-x) : -std::log(y)), a, lower);
    // I_x(a, b) = 1 − I_y(b, a)
    return pgamma_raw((a + (b - 1) / 2) * (y < 0.5 ? -std::log1p(-y) : -std::log(x)), b, !lower);
  }
  // Series in whichever variable makes (a+b)·var small (one shape parameter
  // small, the other large, near the mean).  Accept the result when the
  // complement does not lose more than ~3 digits; otherwise use the CF.
  const double s = a + b;
  if (y <= 0.5 && s * y <= b + 60) {
    const double J = std::fmin(beta_series(b, a, y, x), 1.0);  // = 1 − I_x(a,b)
    if (!lower) return J;
    if (J <= 0.99) return 1 - J;
    return beta_upper_tail_integral(b, a, y);  // I_x(a,b) = P(Y ≥ y), Y = 1 − X
  }
  if (x <= 0.5 && s * x <= a + 60) {
    const double I = std::fmin(beta_series(a, b, x, y), 1.0);
    if (lower) return I;
    if (I <= 0.99) return 1 - I;
    return beta_upper_tail_integral(a, b, x);  // 1 − I_x(a,b) = P(X ≥ x)
  }
  if (x > (a + 1) / (a + b + 2)) {
    double r = beta_prefactor(b, a, y, x) / b * beta_cf(b, a, y);  // = 1 − I_x(a,b)
    if (r > 1) r = 1;
    return lower ? 1 - r : r;
  }
  double r = beta_prefactor(a, b, x, y) / a * beta_cf(a, b, x);  // = I_x(a,b)
  if (r > 1) r = 1;
  return lower ? r : 1 - r;
}

double pbeta(double x, double a, double b, bool lower) { return pbeta_raw(x, a, b, lower, 1 - x); }

double dbeta(double x, double a, double b, double y) {
  if (x < 0 || y < 0) return 0;
  if (x == 0) return a < 1 ? INF : (a == 1 ? b : 0.0);
  if (y == 0) return b < 1 ? INF : (b == 1 ? a : 0.0);
  return beta_prefactor(a, b, x, y) / (x * y);
}

std::pair<double, double> qbeta_pair(double p, double a, double b, bool lower) {
  if (std::isnan(p) || p < 0 || p > 1 || !(a > 0) || !(b > 0)) return {NaN, NaN};
  if (lower ? p == 0 : p == 1) return {0.0, 1.0};
  if (lower ? p == 1 : p == 0) return {1.0, 0.0};
  // quantiles that underflow (the power-law tail is exact there)
  const double lb = lbeta(a, b);
  const double pLow = lower ? p : 1 - p, pUp = lower ? 1 - p : p;
  if (pLow <= 0.5 && (std::log(pLow) + std::log(a) + lb) / a < -744) return {0.0, 1.0};
  if (pUp <= 0.5 && (std::log(pUp) + std::log(b) + lb) / b < -744) return {1.0, 0.0};
  const double x0 = beta_guess(p, a, b, lower);
  if (x0 > 0.5) {
    const double y = solve_beta(p, b, a, !lower, 1 - x0);
    return {1 - y, y};
  }
  const double x = solve_beta(p, a, b, lower, x0);
  return {x, 1 - x};
}

double qbeta(double p, double a, double b, bool lower) { return qbeta_pair(p, a, b, lower).first; }

// ------------------------------------------------------------------ Student t

double pt(double t, double df, bool lower) {
  if (std::isnan(t) || std::isnan(df) || df <= 0) return NaN;
  if (std::isinf(t)) return ((t < 0) == lower) ? 0.0 : 1.0;
  if (std::isinf(df) || df > 1e20) return pnorm(t, lower);
  const double t2 = t * t;
  double x, y;
  if (t2 < df) {
    const double r = t2 / df;
    x = 1 / (1 + r);
    y = r / (1 + r);
  } else {
    const double r = df / t2;
    x = r / (1 + r);
    y = 1 / (1 + r);
  }
  const double tail = 0.5 * pbeta_raw(x, df / 2, 0.5, true, y);  // P(T > |t|)
  if (t > 0) return lower ? 1 - tail : tail;
  return lower ? tail : 1 - tail;
}

double dt(double t, double df) {
  if (std::isinf(df) || df > 1e20) return dnorm(t);
  const double lc = lgamma_diff(df / 2, 0.5) - 0.5 * std::log(df * PI);
  return std::exp(lc - (df + 1) / 2 * std::log1p(t * t / df));
}

double qt(double p, double df, bool lower) {
  if (std::isnan(p) || p < 0 || p > 1 || !(df > 0)) return NaN;
  if (std::isinf(df) || df > 1e20) return qnorm(p, lower);
  if (p == 0) return lower ? -INF : INF;
  if (p == 1) return lower ? INF : -INF;
  if (p == 0.5) return 0;
  double s, sign;
  if (lower) {
    if (p < 0.5) { s = p; sign = -1; } else { s = 1 - p; sign = 1; }
  } else {
    if (p < 0.5) { s = p; sign = 1; } else { s = 1 - p; sign = -1; }
  }
  double q;
  if (df == 1) {
    q = 1 / std::tan(PI * s);
  } else if (df == 2) {
    q = (1 - 2 * s) / std::sqrt(2 * s * (1 - s));
  } else {
    const auto xy = qbeta_pair(2 * s, df / 2, 0.5, true);
    q = std::sqrt(df * xy.second / xy.first);
    for (int i = 0; i < 3; ++i) {  // Newton polish on P(T > q) = s
      const double f = pt(q, df, false) - s;
      const double dens = dt(q, df);
      if (!(dens > 0)) break;
      const double step = f / dens;
      const double qn = q + step;
      if (!std::isfinite(qn) || qn <= 0) break;
      q = qn;
      if (std::fabs(step) <= 1e-15 * q) break;
    }
  }
  return sign * q;
}

// ------------------------------------------------------------------ F

double pf(double f, double df1, double df2, bool lower) {
  if (std::isnan(f) || !(df1 > 0) || !(df2 > 0)) return NaN;
  if (f <= 0) return lower ? 0.0 : 1.0;
  if (std::isinf(f)) return lower ? 1.0 : 0.0;
  if (std::isinf(df2) || df2 > 1e15) return pgamma_raw(df1 * f / 2, df1 / 2, lower);
  const double den = df1 * f + df2;
  return pbeta_raw(df1 * f / den, df1 / 2, df2 / 2, lower, df2 / den);
}

double qf(double p, double df1, double df2, bool lower) {
  if (std::isnan(p) || p < 0 || p > 1 || !(df1 > 0) || !(df2 > 0)) return NaN;
  if (std::isinf(df2) || df2 > 1e15) return 2 * qgamma(p, df1 / 2, lower) / df1;
  const auto xy = qbeta_pair(p, df1 / 2, df2 / 2, lower);
  if (xy.second == 0) return INF;
  return df2 / df1 * xy.first / xy.second;
}

// ------------------------------------------------------------------ binomial

double dbinom(double k, double n, double p) {
  if (k < 0 || k > n || k != std::floor(k)) return 0;
  return dbinom_raw(k, n, p, 1 - p);
}

double pbinom(double k, double n, double p, bool lower) {
  const double kk = std::floor(k + 1e-7);
  if (kk < 0) return lower ? 0.0 : 1.0;
  if (kk >= n) return lower ? 1.0 : 0.0;
  if (p <= 0) return lower ? 1.0 : 0.0;
  if (p >= 1) return lower ? 0.0 : 1.0;
  return pbeta_raw(1 - p, n - kk, kk + 1, lower, p);  // P(X ≤ k) = I_{1−p}(n − k, k + 1)
}

double qbinom(double p, double n, double prob) {
  if (std::isnan(p) || p < 0 || p > 1) return NaN;
  if (p <= 0) return 0;
  if (p >= 1) return n;
  p *= 1 - 64 * DBL_EPSILON_;  // left-continuity fuzz (exact CDF values such as 1/2)
  double lo = -1, hi = n;  // invariant: P(X ≤ lo) < p ≤ P(X ≤ hi)
  const double mu = n * prob, sd = std::sqrt(n * prob * (1 - prob));
  double g = std::floor(mu + sd * qnorm(p));
  if (g > lo && g < hi) {
    if (pbinom(g, n, prob) >= p)
      hi = g;
    else
      lo = g;
  }
  while (hi - lo > 1) {
    const double mid = std::floor(0.5 * (lo + hi));
    if (pbinom(mid, n, prob) >= p)
      hi = mid;
    else
      lo = mid;
  }
  return hi;
}

}  // namespace ssc::stats
