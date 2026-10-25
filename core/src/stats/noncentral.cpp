// Noncentral beta / F, chi-square and t distributions.
//
// Each is a Poisson mixture of central distributions.  Instead of summing from
// j = 0 (O(λ) terms, underflow for large λ, the reason classical codes switch
// to crude normal approximations), we follow Benton & Krishnamoorthy (2003,
// CSDA 43:249–267): evaluate the central CDF once at the mode of the Poisson
// weights and walk outwards in both directions with the exact recurrences
//
//   I_x(c+1, b) = I_x(c, b) − g(c),   g(c+1) = g(c)·x(c+b)/(c+1)
//   P(s+1, y)   = P(s, y)   − h(s),   h(s+1) = h(s)·y/(s+1)
//
// so every additional term costs a few flops.  Summation stops when a rigorous
// geometric bound on the neglected Poisson mass is below 1e-16 of the sum (of
// 2^-54 at least where the sum is taken from 1, see complement()), or when the
// weights fall below the smallest normal double.
// Cost is O(√λ).  For very large noncentrality (Poisson mean above 1e5) the
// sum is replaced by an integral that costs O(1) — see poisson_mixture().
#include <cmath>

#include "numerics.hpp"
#include "stats.hpp"

namespace ssc::stats {

namespace {
constexpr double REL_EPS = 1e-16;
constexpr long MAX_TERMS = 5000000;
constexpr double LARGE_MEAN = 1e5;
// What a sum taken from 1 needs: beside 1, a tail below 2^-54 is 0 (see complement()).  Summed to 1e-16 of itself, a
// tail of 1e-300 or 0 ran on until its weights reached 0, and they stay at the smallest denormal while the ratio from
// one to the next is above ½: for half the Poisson mean (49,601 terms for a power of 1 at λ = 198,400).  Integrated
// to 1e-16 of itself (Poisson means above 1e5), a tail of 5e-109 took 13,701 terms where 321 do.
constexpr double LEAST_FROM_ONE = 0x1p-54;

// Poisson weight e^{−μ}μ^j/Γ(j + 1) at j = μ + δ, from the offset δ so that it
// stays accurate when j itself cannot resolve δ (μ beyond 2^53 / 1e3):
//   p = exp(−δ(j) − μ·g(δ/μ)) / √(2πj),  g(t) = (1 + t)·ln(1 + t) − t.
double poisson_weight(double mu, double delta) {
  const double j = mu + delta, t = delta / mu;
  double g;
  if (std::fabs(t) < 0.1) {
    g = 0;  // Σ_{n≥2} (−t)^n / (n(n − 1))
    double tn = t * t;
    for (int n = 2; n < 60; ++n) {
      const double term = tn / (n * (n - 1.0));
      g += (n % 2 == 0) ? term : -term;
      if (std::fabs(term) <= 1e-17 * std::fabs(g)) break;
      tn *= t;
    }
  } else {
    g = (1 + t) * std::log1p(t) - t;
  }
  return std::exp(-stirlerr(j) - mu * g) / std::sqrt(2 * PI * j);
}

// Σ_j p_j(μ) T(μ + δ) for a large Poisson mean μ, with the term written as a
// function of the offset δ = j − μ.  When T varies slowly on the
// unit scale, the sum equals the integral of the continuous extension
// p(j) = e^{−μ}μ^j/Γ(j + 1) times T(j) up to terms of order exp(−2π²μ)
// (Poisson summation formula).  With j = μ + √μ·u the integrand decays like a
// Gaussian in u, so the trapezoidal rule converges geometrically (error about
// 2·exp(−2π²/h²) for step h); halving the step from ½ to ¼ checks the result,
// and adaptive Gauss–Kronrod takes over if T has a sharp feature.  Both to
// their tolerance of the sum, or of `least` if larger (see LEAST_FROM_ONE).
template <class T>
double poisson_mixture_offset(double mu, T&& term, double least = 0) {
  const double sd = std::sqrt(mu);
  auto f = [&](double u) {
    const double delta = sd * u;
    return mu + delta <= 0 ? 0.0 : poisson_weight(mu, delta) * term(delta);
  };
  constexpr int U = 40;  // the Poisson mass beyond ±40 sd is below 1e-300
  double coarse = 0;
  for (int i = -2 * U; i <= 2 * U; ++i) coarse += f(0.5 * i);
  double fine = coarse;
  for (int i = -2 * U; i < 2 * U; ++i) fine += f(0.5 * i + 0.25);
  coarse *= 0.5 * sd;
  fine *= 0.25 * sd;
  const double scale = std::fmax(std::fabs(fine), least);
  if (std::fabs(fine - coarse) <= 1e-14 * scale) return fine;
  const double tol = std::max(1e-300, 1e-16 * scale);
  double s = 0;
  for (int i = -U; i < U; i += 4) s += sd * numerics::integrate(f, i, i + 4, tol / 20, 1e-14);
  return s;
}

// The same with the term written as a function of j.
template <class T>
double poisson_mixture(double mu, T&& term, double least = 0) {
  return poisson_mixture_offset(mu, [&](double delta) { return term(mu + delta); }, least);
}
}  // namespace

// A tail near 1 is 1 less the other one, which is summed instead.  Summed directly, the terms of a tail near 1 add up
// to it only to a few units in the last place, up or down from one argument to the next: a power next to 1 did not
// grow steadily with n or the effect, and searches stopped wherever it first crossed the target (a sample size of
// 797 for 208 at a target of 0.9999999999999999, of 199 for 200 at 0.9999999999999995).  The smaller tail is told
// by the central distribution at the Poisson weights' mode (or mean), the term of most weight: the sum is near 1
// only where every term of much weight is.  1 − q rounded to the nearest double can reach a target that 1 − q falls
// short of (1 − 1.2e-16 rounds to 0.9999999999999999, as 1 − 1.07e-16 does): complement() rounds it down.
double complement(double q) {
  const double p = 1 - q;  // 1 − p is exact (p ≥ ½)
  return 1 - p >= q || q < 0x1p-54 ? p : std::nextafter(p, 0.0);
}

namespace {
double pnbeta_tail(double x, double a, double b, double ncp, bool lower, double y) {
  const double mu = ncp / 2;
  if (mu > LARGE_MEAN) {
    const bool flip = pbeta_raw(x, a + mu, b, lower, y) > 0.5;  // the other tail is summed (see complement())
    const bool low = lower != flip;
    const double sum = clamp01(poisson_mixture(
        mu, [&](double j) { return pbeta_raw(x, a + j, b, low, y); }, flip ? LEAST_FROM_ONE : 0.0));
    return flip ? complement(sum) : sum;
  }
  const double k = std::floor(mu);
  const double pk = dpois_raw(k, mu);
  const double ak = a + k;
  double Ik = pbeta_raw(x, ak, b, lower, y);  // lower: I, upper: J = 1 − I
  const bool flip = Ik > 0.5;  // the other tail is summed (see complement())
  if (flip) {
    lower = !lower;
    Ik = pbeta_raw(x, ak, b, lower, y);
  }
  const double least = flip ? LEAST_FROM_ONE : 0.0;
  const double gk = beta_prefactor(ak, b, x, y) / ak;
  double sum = pk * Ik;

  // backward j = k−1, …, 0
  double p = pk, I = Ik, g = gk;
  for (double j = k - 1; j >= 0; j -= 1) {
    const double c = a + j;
    g = (g == 0) ? beta_prefactor(c, b, x, y) / c : g * (c + 1) / (x * (c + b));
    I = lower ? I + g : I - g;
    if (I < 0) I = 0;
    p = p * (j + 1) / mu;
    sum += p * I;
    const double r = j / mu;  // p_{j−1}/p_j
    const double bound = p * r / (1 - r) * (lower ? 1.0 : I);
    if (bound <= REL_EPS * std::fmax(sum, least) || p < DBL_MIN_) break;
  }
  // forward j = k+1, …
  p = pk;
  I = Ik;
  g = gk;
  for (long jj = 1; jj < MAX_TERMS; ++jj) {
    const double j = k + jj;
    const double cPrev = a + j - 1;
    if (g == 0) g = beta_prefactor(cPrev, b, x, y) / cPrev;
    I = lower ? I - g : I + g;
    if (I < 0) I = 0;
    g = g * x * (cPrev + b) / (cPrev + 1);
    p = p * mu / j;
    sum += p * I;
    const double r = mu / (j + 1);
    if (r < 1) {
      const double bound = p * r / (1 - r) * (lower ? I : 1.0);
      if (bound <= REL_EPS * std::fmax(sum, least) || p < DBL_MIN_) break;
    }
  }
  return flip ? complement(clamp01(sum)) : clamp01(sum);
}
}  // namespace

double pnbeta(double x, double a, double b, double ncp, bool lower, double y) {
  if (std::isnan(x) || !(a > 0) || !(b > 0) || !(ncp >= 0)) return NaN;
  if (ncp == 0) return pbeta_raw(x, a, b, lower, y);
  if (x <= 0) return lower ? 0.0 : 1.0;
  if (y <= 0) return lower ? 1.0 : 0.0;
  if (std::isinf(ncp)) return lower ? 0.0 : 1.0;  // all the mass at 1 (the series would never end)
  return pnbeta_tail(x, a, b, ncp, lower, y);
}

double pnf(double f, double df1, double df2, double ncp, bool lower) {
  if (std::isnan(f) || !(df1 > 0) || !(df2 > 0) || !(ncp >= 0)) return NaN;
  if (f <= 0) return lower ? 0.0 : 1.0;
  if (std::isinf(f)) return lower ? 1.0 : 0.0;
  if (std::isinf(df2) || df2 > 1e15) return pnchisq(df1 * f, df1, ncp, lower);
  const double den = df1 * f + df2;
  return pnbeta(df1 * f / den, df1 / 2, df2 / 2, ncp, lower, df2 / den);
}

namespace {
double pnchisq_tail(double x, double df, double ncp, bool lower) {
  const double mu = ncp / 2;
  const double yv = x / 2;
  if (mu > LARGE_MEAN) {
    // shape s_j = df/2 + j and the offset yv − s_j, kept exact as j = mu + δ
    const double d0 = (yv - mu) - df / 2;
    const bool flip = pgamma_offset(df / 2 + mu, d0, lower) > 0.5;  // the other tail is summed (see complement())
    const bool low = lower != flip;
    const double sum = clamp01(poisson_mixture_offset(
        mu, [&](double delta) { return pgamma_offset(df / 2 + mu + delta, d0 - delta, low); },
        flip ? LEAST_FROM_ONE : 0.0));
    return flip ? complement(sum) : sum;
  }
  const double k = std::floor(mu);
  const double pk = dpois_raw(k, mu);
  const double sk = df / 2 + k;
  auto central = [&] { return (sk == 0) ? (lower ? 1.0 : 0.0) : pgamma_raw(yv, sk, lower); };
  double Pk = central();
  const bool flip = Pk > 0.5;  // the other tail is summed (see complement())
  if (flip) {
    lower = !lower;
    Pk = central();
  }
  const double least = flip ? LEAST_FROM_ONE : 0.0;
  const double hk = dpois_raw(sk, yv);  // P(s) − P(s+1)
  double sum = pk * Pk;

  double p = pk, P = Pk, h = hk;
  for (double j = k - 1; j >= 0; j -= 1) {
    const double s = df / 2 + j;
    h = (h == 0) ? dpois_raw(s, yv) : h * (s + 1) / yv;
    P = lower ? P + h : P - h;
    if (P < 0) P = 0;
    p = p * (j + 1) / mu;
    sum += p * P;
    const double r = j / mu;
    const double bound = p * r / (1 - r) * (lower ? 1.0 : P);
    if (bound <= REL_EPS * std::fmax(sum, least) || p < DBL_MIN_) break;
  }
  p = pk;
  P = Pk;
  h = hk;
  for (long jj = 1; jj < MAX_TERMS; ++jj) {
    const double j = k + jj;
    const double sPrev = df / 2 + j - 1;
    if (h == 0) h = dpois_raw(sPrev, yv);
    P = lower ? P - h : P + h;
    if (P < 0) P = 0;
    h = h * yv / (sPrev + 1);
    p = p * mu / j;
    sum += p * P;
    const double r = mu / (j + 1);
    if (r < 1) {
      const double bound = p * r / (1 - r) * (lower ? P : 1.0);
      if (bound <= REL_EPS * std::fmax(sum, least) || p < DBL_MIN_) break;
    }
  }
  return flip ? complement(clamp01(sum)) : clamp01(sum);
}
}  // namespace

double pnchisq(double x, double df, double ncp, bool lower) {
  if (std::isnan(x) || !(df >= 0) || !(ncp >= 0)) return NaN;
  if (ncp == 0) return pchisq(x, df, lower);
  if (x <= 0) return lower ? 0.0 : 1.0;
  if (std::isinf(x)) return lower ? 1.0 : 0.0;
  if (std::isinf(ncp)) return lower ? 0.0 : 1.0;  // all the mass at infinity (the series would never end)
  return pnchisq_tail(x, df, ncp, lower);
}

namespace {
// Noncentral t for t ≥ 0 (Guenther 1978; Lenth 1989):
//   F(t) = Φ(−δ) + ½ Σ_j [ p_j I_x(j+½, ν/2) + (δ/√2) q_j I_x(j+1, ν/2) ],
//   x = t²/(t²+ν), p_j = e^{−λ}λ^j/j!, q_j = e^{−λ}λ^j/Γ(j+3/2), λ = δ²/2.
// Upper tail: ½ Σ_j [ p_j (1 − I) + (δ/√2) q_j (1 − I) ], summed directly.  phi_lower = Φ(−δ).  For δ ≥ 0 the
// tail near 1 is 1 less the other one (see complement()); for δ < 0 the upper tail is the smaller one (see pnt_pos).
double pnt_series(double t, double df, double del, bool lower, double phi_lower) {
  const double t2 = t * t;
  double x, y;
  if (t2 < df) {
    const double r = t2 / df;
    x = r / (1 + r);
    y = 1 / (1 + r);
  } else {
    const double r = df / t2;
    x = 1 / (1 + r);
    y = r / (1 + r);
  }
  const double b = df / 2;
  const double lam = del * del / 2;
  const double d = del / std::sqrt(2.0);
  if (lam > LARGE_MEAN) {
    const bool flip = del >= 0 && pbeta_raw(x, lam + 0.5, b, lower, y) > 0.5;
    const bool low = lower != flip;
    // q_j / p_j = Γ(j + 1)/Γ(j + 3/2)
    const double s = poisson_mixture(
        lam,
        [&](double j) {
          return pbeta_raw(x, j + 0.5, b, low, y) +
                 d * std::exp(-lgamma_diff(j + 1, 0.5)) * pbeta_raw(x, j + 1, b, low, y);
        },
        flip ? LEAST_FROM_ONE : 0.0);
    const double tail = clamp01(low ? phi_lower + 0.5 * s : 0.5 * s);
    return flip ? complement(tail) : tail;
  }
  const double k = std::floor(lam);
  const double pk = dpois_raw(k, lam);
  const double qk = pk * std::exp(-lgamma_diff(k + 1, 0.5));  // e^{−λ}λ^k/Γ(k+3/2)
  const double a1 = k + 0.5, a2 = k + 1;
  double I1k = pbeta_raw(x, a1, b, lower, y);
  const bool flip = del >= 0 && I1k > 0.5;
  if (flip) {
    lower = !lower;
    I1k = pbeta_raw(x, a1, b, lower, y);
  }
  const double least = flip ? LEAST_FROM_ONE : 0.0;
  const double I2k = pbeta_raw(x, a2, b, lower, y);
  const double g1k = beta_prefactor(a1, b, x, y) / a1;
  const double g2k = beta_prefactor(a2, b, x, y) / a2;
  double sum = pk * I1k + d * qk * I2k;

  double p = pk, q = qk, I1 = I1k, I2 = I2k, g1 = g1k, g2 = g2k;
  for (double j = k - 1; j >= 0; j -= 1) {
    const double c1 = j + 0.5, c2 = j + 1;
    g1 = (g1 == 0) ? beta_prefactor(c1, b, x, y) / c1 : g1 * (c1 + 1) / (x * (c1 + b));
    g2 = (g2 == 0) ? beta_prefactor(c2, b, x, y) / c2 : g2 * (c2 + 1) / (x * (c2 + b));
    if (lower) {
      I1 += g1;
      I2 += g2;
    } else {
      I1 = std::fmax(I1 - g1, 0.0);
      I2 = std::fmax(I2 - g2, 0.0);
    }
    p = p * (j + 1) / lam;
    q = q * (j + 1.5) / lam;
    sum += p * I1 + d * q * I2;
    const double r = (j + 0.5) / lam;  // bounds both weight ratios
    if (r < 1) {
      const double w = p + std::fabs(d) * q;
      const double bound = w * r / (1 - r) * (lower ? 1.0 : std::fmax(I1, I2));
      if (bound <= REL_EPS * std::fmax(std::fabs(sum), least) || w < DBL_MIN_) break;
    }
  }
  p = pk;
  q = qk;
  I1 = I1k;
  I2 = I2k;
  g1 = g1k;
  g2 = g2k;
  for (long jj = 1; jj < MAX_TERMS; ++jj) {
    const double j = k + jj;
    const double c1 = j - 0.5, c2 = j;  // previous shapes
    if (g1 == 0) g1 = beta_prefactor(c1, b, x, y) / c1;
    if (g2 == 0) g2 = beta_prefactor(c2, b, x, y) / c2;
    if (lower) {
      I1 = std::fmax(I1 - g1, 0.0);
      I2 = std::fmax(I2 - g2, 0.0);
    } else {
      I1 += g1;
      I2 += g2;
    }
    g1 = g1 * x * (c1 + b) / (c1 + 1);
    g2 = g2 * x * (c2 + b) / (c2 + 1);
    p = p * lam / j;
    q = q * lam / (j + 0.5);
    sum += p * I1 + d * q * I2;
    const double r = lam / (j + 1);
    if (r < 1) {
      const double w = p + std::fabs(d) * q;
      const double bound = w * r / (1 - r) * (lower ? std::fmax(I1, I2) : 1.0);
      if (bound <= REL_EPS * std::fmax(std::fabs(sum), least) || w < DBL_MIN_) break;
    }
  }
  const double tail = clamp01(lower ? phi_lower + 0.5 * sum : 0.5 * sum);
  return flip ? complement(tail) : tail;
}

// ln P(Z > y), also beyond y = 37, where P(Z > y) nears the smallest double: there from its asymptotic series.
double log_upper_normal(double y) {
  if (y < 37) return std::log(pnorm(y, false));
  const double r = 1 / (y * y);
  return -y * y / 2 - std::log(y) - LN_SQRT_2PI + std::log1p(-r * (1 - 3 * r * (1 - 5 * r * (1 - 7 * r))));
}

// e^x − 1 − x, also next to 0, where expm1(x) − x cancels (to 1e-10 of it at x = 1e-6)
double expm1_less_x(double x) {
  if (std::fabs(x) > 1e-3) return std::expm1(x) - x;
  double term = x * x / 2, sum = term;
  for (int n = 3; std::fabs(term) > 1e-17 * sum; ++n) {
    term *= x / n;
    sum += term;
  }
  return sum;
}

// P(T > t) for t > 0 and δ = −a < 0, which is P(Z > a + tS) with S = √(χ²_ν/ν), integrated over u = log S,
// whose density is exp(g(u)) / ∫exp(g), g(u) = νu − ν(e^{2u} − 1)/2 = −ν(e^{2u} − 1 − 2u)/2, times P(Z > a + t·e^u):
// the logarithm h(u) of the product is concave, and Gauss–Legendre takes it in panels of at most 1.5 over the range
// where it is within 40 of its largest value (as the SD ratio's nodes of sim.multiarm).
double pnt_far(double t, double df, double a) {
  // P(T > t) ≤ P(S < ½) + P(Z > a + t/2): where both are below the doubles, so is the tail (P(S < ½) is, from
  // ν = 2400 on).  With ν and t both large (ν = 1e19, t = 1e10) h is of the order of −1e19, its rounding far above 1,
  // and the peak and the range below lost their meaning: NaN or ±∞.
  if (df > 1e3 && pgamma_raw(df / 8, df / 2, true) + pnorm(a + t / 2, false) == 0) return 0;
  auto g = [&](double u) { return -df / 2 * expm1_less_x(2 * u); };
  auto h = [&](double u) { return g(u) + log_upper_normal(a + t * std::exp(u)); };
  // the largest value: below u = 0, where h falls already (g'(0) = 0), bracketed in steps that double going down,
  // then by golden section to a twentieth of the width of the peak of g, 1/√(2ν).  Beyond t = 1e150 from where
  // t·e^u is 1e150: the square in h overflowed at u = 0 and next to it, h was −∞ at both, and the search threw.
  constexpr double G = 0.6180339887498949;
  const double step = std::fmin(1 / std::sqrt(2 * df), 1.0), tol = 0.05 * step;
  double hi = std::fmin(0.0, std::log(1e150 / t)), mid = hi - step, h_mid = h(mid), lo = mid;
  if (h_mid > h(hi)) {  // further down: until h falls again, so that the largest value lies between lo and hi
    for (lo = mid - 2 * step;; lo = mid - 2 * (hi - mid)) {
      const double h_lo = h(lo);
      if (h_lo <= h_mid || lo < -1e3) break;
      hi = mid;
      mid = lo;
      h_mid = h_lo;
    }
  }
  double m1 = hi - G * (hi - lo), m2 = lo + G * (hi - lo), h1 = h(m1), h2 = h(m2);
  while (hi - lo > tol) {
    if (h1 < h2) {
      lo = m1;
      m1 = m2;
      h1 = h2;
      m2 = lo + G * (hi - lo);
      h2 = h(m2);
    } else {
      hi = m2;
      m2 = m1;
      h2 = h1;
      m1 = hi - G * (hi - lo);
      h1 = h(m1);
    }
  }
  const double top = (lo + hi) / 2, peak = h(top), level = peak - 40;
  // the ends of the range, bracketed in steps that double from 8 widths out (about where e^g falls by e^−32)
  auto end = [&](double dir) {
    double near = top, d = 8 * step, far = top + dir * d;
    while (h(far) > level) {
      near = far;
      far += dir * (d *= 2);
    }
    return numerics::brent([&](double v) { return h(v) - level; }, near, far, tol);
  };
  // Below uc, where t·e^u·(a + 1) < 2^-56, P(Z > a + t·e^u) is P(Z > a) to 1e-16 of itself (its relative slope is
  // below a + 1), and that part of the integral is P(Z > a) times P(S < e^uc), a chi-square tail.  The range within
  // 40 of the peak reaches some 40/ν below it, where g rises at the rate ν, and was taken in panels all the same: 7.5 s
  // at ν = 1e-5, ten times that for each factor of 10 below it, and below ν = 1.2e-8 a count of panels past the
  // largest int, which gave 0.  Above uc the range is at most 43 wide (ν from 5e-324 to 1e20, t from 1e-300 to 1e300).
  const double ub = end(1), uc = std::fmin(std::log(0x1p-56 / (t * (a + 1))), ub);
  const bool cut = uc >= top || h(uc) > level;
  const double ua = cut ? uc : end(-1);
  double below = 0;
  if (cut) {  // P(S < e^uc) = P(χ²_ν < ν·e^{2uc}) = P(k, x), k = ν/2, from x^k/Γ(k + 1) where x is below the doubles
    const double lx = std::log(df) - std::log(2.0) + 2 * uc;
    below = pnorm(a, false) * (lx > -700 ? pgamma_raw(std::exp(lx), df / 2, true)
                                         : std::exp(df / 2 * lx - lgammafn(df / 2 + 1)));
  }
  if (!(ub > ua)) return below;
  const GLRule& rule = gauss_legendre(48);
  const int panels = static_cast<int>(std::ceil((ub - ua) / 1.5));
  const double half = (ub - ua) / panels / 2;
  double sum = 0;
  for (int q = 0; q < panels; ++q)
    for (std::size_t i = 0; i < rule.x.size(); ++i) {
      const double u = ua + half * (2 * q + 1 + rule.x[i]), y = a + t * std::exp(u);
      sum += rule.w[i] *
             (y < 37 ? std::exp(g(u) - peak) * pnorm(y, false) : std::exp(g(u) + log_upper_normal(y) - peak));
    }
  // ln ∫exp(g) = lnΓ(k) + k − k·ln k − ln 2 with k = ν/2, by Stirling's series without cancellation
  const double k = df / 2, log_total = LN_SQRT_2PI - 0.5 * std::log(k) + stirlerr(k) - std::log(2.0);
  return below + std::exp(peak + std::log(half * sum) - log_total);
}

double pnt_pos(double t, double df, double del, bool lower) {
  // An enormous noncentrality (t > 0, δ > 0): Z is nothing beside δ, so T = (Z + δ)/S behaves as δ/S
  // with S = √(χ²_ν/ν), and P(T ≤ t) = P(χ²_ν ≥ ν(δ/t)²).  (Beyond δ ≈ 1.3e154, δ²/2 overflows and the
  // series cannot be summed at all.)
  if (del >= 1e8) {
    const double q = del / t;
    return pchisq(df * q * q, df, !lower);
  }
  const auto phi = pnorm_both(-del);  // {Φ(−δ), 1 − Φ(−δ)}
  if (t == 0) return lower ? phi.first : phi.second;
  if (del < 0) {
    // The upper tail, at most Φ(δ), is the smaller one.  Its series cancels down to it from terms of about ½ and
    // leaves it only to about 1e-16 (1.5e-16 for 1.2e-32 at t = 1.97, ν = 392 and δ = −9.9): as much as a power next
    // to 1 falls short of 1, and more, so that the power of a two-sided t-test, which adds this tail, moved up and
    // down by that much from one n to the next.  Below 1e-8 it is integrated instead.
    double q = phi.second < 1e-8 ? 0.0 : pnt_series(t, df, del, false, phi.first);
    if (q < 1e-8) q = pnt_far(t, df, -del);
    return lower ? complement(q) : q;
  }
  return pnt_series(t, df, del, lower, phi.first);
}
}  // namespace

double pnt(double t, double df, double ncp, bool lower) {
  if (std::isnan(t) || std::isnan(df) || std::isnan(ncp) || df <= 0) return NaN;
  if (ncp == 0) return pt(t, df, lower);
  if (std::isinf(t)) return ((t < 0) == lower) ? 0.0 : 1.0;
  if (std::isinf(df) || df > 1e20) {
    const auto pr = pnorm_both(t - ncp);
    return lower ? pr.first : pr.second;
  }
  // Far tails: T ≤ t ≤ 0 needs Z + δ ≤ 0, so P(T ≤ t) ≤ Φ(−δ), which is below the
  // smallest double once δ ≥ 39 (likewise for t ≥ 0 and δ ≤ −39).  The series
  // would only add up rounding noise there, and slowly.
  if (t <= 0 && ncp >= 39) return lower ? 0.0 : 1.0;
  if (t >= 0 && ncp <= -39) return lower ? 1.0 : 0.0;
  if (t < 0) return pnt_pos(-t, df, -ncp, !lower);
  return pnt_pos(t, df, ncp, lower);
}

}  // namespace ssc::stats
