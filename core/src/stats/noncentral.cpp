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
// geometric bound on the neglected Poisson mass is below 1e-16 of the sum.
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
// and adaptive Gauss–Kronrod takes over if T has a sharp feature.
template <class T>
double poisson_mixture_offset(double mu, T&& term) {
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
  if (std::fabs(fine - coarse) <= 1e-14 * std::fabs(fine)) return fine;
  const double tol = std::max(1e-300, 1e-16 * std::fabs(fine));
  double s = 0;
  for (int i = -U; i < U; i += 4) s += sd * numerics::integrate(f, i, i + 4, tol / 20, 1e-14);
  return s;
}

// The same with the term written as a function of j.
template <class T>
double poisson_mixture(double mu, T&& term) {
  return poisson_mixture_offset(mu, [&](double delta) { return term(mu + delta); });
}
}  // namespace

double pnbeta(double x, double a, double b, double ncp, bool lower, double y) {
  if (std::isnan(x) || !(a > 0) || !(b > 0) || !(ncp >= 0)) return NaN;
  if (ncp == 0) return pbeta_raw(x, a, b, lower, y);
  if (x <= 0) return lower ? 0.0 : 1.0;
  if (y <= 0) return lower ? 1.0 : 0.0;
  if (std::isinf(ncp)) return lower ? 0.0 : 1.0;  // all the mass at 1 (the series would never end)
  const double mu = ncp / 2;
  if (mu > LARGE_MEAN)
    return clamp01(poisson_mixture(mu, [&](double j) { return pbeta_raw(x, a + j, b, lower, y); }));
  const double k = std::floor(mu);
  const double pk = dpois_raw(k, mu);
  const double ak = a + k;
  const double Ik = pbeta_raw(x, ak, b, lower, y);  // lower: I, upper: J = 1 − I
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
    if (bound <= REL_EPS * sum || p == 0) break;
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
      if (bound <= REL_EPS * sum || p == 0) break;
    }
  }
  return clamp01(sum);
}

double pnf(double f, double df1, double df2, double ncp, bool lower) {
  if (std::isnan(f) || !(df1 > 0) || !(df2 > 0) || !(ncp >= 0)) return NaN;
  if (f <= 0) return lower ? 0.0 : 1.0;
  if (std::isinf(f)) return lower ? 1.0 : 0.0;
  if (std::isinf(df2) || df2 > 1e15) return pnchisq(df1 * f, df1, ncp, lower);
  const double den = df1 * f + df2;
  return pnbeta(df1 * f / den, df1 / 2, df2 / 2, ncp, lower, df2 / den);
}

double pnchisq(double x, double df, double ncp, bool lower) {
  if (std::isnan(x) || !(df >= 0) || !(ncp >= 0)) return NaN;
  if (ncp == 0) return pchisq(x, df, lower);
  if (x <= 0) return lower ? 0.0 : 1.0;
  if (std::isinf(x)) return lower ? 1.0 : 0.0;
  if (std::isinf(ncp)) return lower ? 0.0 : 1.0;  // all the mass at infinity (the series would never end)
  const double mu = ncp / 2;
  const double yv = x / 2;
  if (mu > LARGE_MEAN) {
    // shape s_j = df/2 + j and the offset yv − s_j, kept exact as j = mu + δ
    const double d0 = (yv - mu) - df / 2;
    return clamp01(poisson_mixture_offset(mu, [&](double delta) {
      return pgamma_offset(df / 2 + mu + delta, d0 - delta, lower);
    }));
  }
  const double k = std::floor(mu);
  const double pk = dpois_raw(k, mu);
  const double sk = df / 2 + k;
  const double Pk = (sk == 0) ? (lower ? 1.0 : 0.0) : pgamma_raw(yv, sk, lower);
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
    if (bound <= REL_EPS * sum || p == 0) break;
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
      if (bound <= REL_EPS * sum || p == 0) break;
    }
  }
  return clamp01(sum);
}

namespace {
// Noncentral t for t ≥ 0 (Guenther 1978; Lenth 1989):
//   F(t) = Φ(−δ) + ½ Σ_j [ p_j I_x(j+½, ν/2) + (δ/√2) q_j I_x(j+1, ν/2) ],
//   x = t²/(t²+ν), p_j = e^{−λ}λ^j/j!, q_j = e^{−λ}λ^j/Γ(j+3/2), λ = δ²/2.
// Upper tail: ½ Σ_j [ p_j (1 − I) + (δ/√2) q_j (1 − I) ], summed directly.
double pnt_pos(double t, double df, double del, bool lower) {
  // An enormous noncentrality (t > 0, δ > 0): Z is nothing beside δ, so T = (Z + δ)/S behaves as δ/S
  // with S = √(χ²_ν/ν), and P(T ≤ t) = P(χ²_ν ≥ ν(δ/t)²).  (Beyond δ ≈ 1.3e154, δ²/2 overflows and the
  // series below cannot be summed at all.)
  if (del >= 1e8) {
    const double q = del / t;
    return pchisq(df * q * q, df, !lower);
  }
  const auto phi = pnorm_both(-del);  // {Φ(−δ), 1 − Φ(−δ)}
  if (t == 0) return lower ? phi.first : phi.second;
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
    // q_j / p_j = Γ(j + 1)/Γ(j + 3/2)
    const double s = poisson_mixture(lam, [&](double j) {
      return pbeta_raw(x, j + 0.5, b, lower, y) + d * std::exp(-lgamma_diff(j + 1, 0.5)) * pbeta_raw(x, j + 1, b, lower, y);
    });
    return clamp01(lower ? phi.first + 0.5 * s : 0.5 * s);
  }
  const double k = std::floor(lam);
  const double pk = dpois_raw(k, lam);
  const double qk = pk * std::exp(-lgamma_diff(k + 1, 0.5));  // e^{−λ}λ^k/Γ(k+3/2)
  const double a1 = k + 0.5, a2 = k + 1;
  const double I1k = pbeta_raw(x, a1, b, lower, y);
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
      if (bound <= REL_EPS * std::fabs(sum) || w == 0) break;
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
      if (bound <= REL_EPS * std::fabs(sum) || w == 0) break;
    }
  }
  return clamp01(lower ? phi.first + 0.5 * sum : 0.5 * sum);
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
