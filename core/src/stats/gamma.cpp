// Regularized incomplete gamma function and the gamma / chi-square / Poisson
// distributions.  P(a, x) uses the power series for x < a + 1 and Q(a, x) a
// modified-Lentz continued fraction otherwise; both use the cancellation-free
// prefactor x^a e^{−x}/Γ(a + 1) = dpois_raw(a, x).
#include "stats.hpp"

namespace ssc::stats {

namespace {
constexpr double TINY = 1e-300;

#include "temme_coefficients.inc"

// m − ln(1 + m) without cancellation for small m.
double m_minus_log1p(double m) {
  if (std::fabs(m) >= 0.5) return m - std::log1p(m);
  double s = 0, mk = m * m;  // Σ_{k≥2} (−m)^k / k
  for (int k = 2; k < 200; ++k) {
    const double term = mk / k;
    s += (k % 2 == 0) ? term : -term;
    if (std::fabs(term) <= 1e-17 * s) break;
    mk *= m;
  }
  return s;
}

// Temme's uniform asymptotic expansion for large a (Temme 1979; DiDonato &
// Morris 1986), with λ = x/a = 1 + m and η = sign(m)·√(2(m − ln(1 + m))):
//   Q(a, x) = Φ̄(η√a) + e^{−aη²/2}/√(2πa)·Σ_n C_n(η) a^{−n},
//   P(a, x) = Φ(η√a)  − e^{−aη²/2}/√(2πa)·Σ_n C_n(η) a^{−n}.
// The Taylor series of C_n converge for |η| < 2√π; with |η| ≤ 1.5 and eight
// terms in 1/a for a ≥ 100 the result is accurate to about 1e-15 relative, in
// O(1) work where the series and the continued fraction need O(√a) terms.
constexpr double TEMME_MIN_A = 100;
constexpr double TEMME_MAX_ETA = 1.5;

double temme_eta(double m) {
  const double e = std::sqrt(2 * m_minus_log1p(m));
  return m < 0 ? -e : e;
}

double temme(double a, double m, double eta, bool lower) {
  double S = 0, an = 1;
  for (int n = 0; n < 8; ++n) {
    const double* c = TEMME_C[n];
    double v = 0;
    for (int k = TEMME_LEN[n] - 1; k >= 0; --k) v = v * eta + c[k];
    S += v * an;
    an /= a;
    if (std::fabs(an) < 1e-18) break;
  }
  const double R = std::exp(-a * 0.5 * eta * eta) / std::sqrt(2 * PI * a) * S;
  const double w = eta * std::sqrt(a);
  return clamp01(lower ? pnorm(w) - R : pnorm(w, false) + R);
}

// 1/Γ(1 + a) − 1 for 0 ≤ a ≤ ½, from the Taylor series of 1/Γ (Abramowitz & Stegun 6.1.34): to full relative
// precision next to a = 0, where 1/Γ(1 + a) itself is 1 to the last place.
double gam1(double a) {
  static constexpr double C[] = {
      0.57721566490153286,     -0.65587807152025388,    -0.042002635034095236,   0.16653861138229149,
      -0.042197734555544337,   -0.0096219715278769736,  0.0072189432466630995,   -0.0011651675918590651,
      -0.00021524167411495097, 0.00012805028238811619,  -2.0134854780788239e-05, -1.2504934821426707e-06,
      1.1330272319816959e-06,  -2.0563384169776071e-07, 6.1160951044814158e-09,  5.0020076444692229e-09,
      -1.1812745704870201e-09, 1.0434267116911005e-10,  7.7822634399050713e-12,  -3.6968056186422057e-12,
      5.1003702874544760e-13};
  double s = 0;
  for (int k = sizeof(C) / sizeof(C[0]) - 1; k >= 0; --k) s = s * a + C[k];
  return a * s;
}

// Q(a, x) for a < ½ and x < a + 1, where 1 − P(a, x) loses Q's digits: Q is about a·E1(x) there, and 1 − P was off
// by about 1e-16/a of it (2.3e-5 at a = 1e-10, x = 1).  As TOMS 708 (GRAT1) takes it, with w = x^a, l = w − 1,
// g = 1/Γ(1 + a), h = g − 1 and J = a·(x/(a + 1) − x²/(2!(a + 2)) + x³/(3!(a + 3)) − …), so that P = w·g·(1 − J):
//   Q = (w·J − l)·g − h,
// or 1 − P itself where w is below 0.875 (Q is then above 0.1).
double small_shape_upper(double x, double a) {
  double c = x, s = x / (a + 1);
  for (int n = 2; n < 100; ++n) {
    c *= -x / n;
    const double term = c / (a + n);
    s += term;
    if (std::fabs(term) <= 1e-17 * std::fabs(s)) break;
  }
  const double J = a * s, z = a * std::log(x), h = gam1(a), g = 1 + h;
  if (z <= -0.13394) return clamp01(1 - std::exp(z) * g * (1 - J));
  const double l = std::expm1(z);
  return clamp01(((l + 1) * J - l) * g - h);
}
}  // namespace

double pgamma_offset(double a, double d, bool lower) {
  if (a >= TEMME_MIN_A) {
    const double m = d / a;
    if (m > -1) {
      const double eta = temme_eta(m);
      if (std::fabs(eta) <= TEMME_MAX_ETA) return temme(a, m, eta, lower);
    }
  }
  return pgamma_raw(a + d, a, lower);
}

double pgamma_raw(double x, double a, bool lower) {
  if (std::isnan(x) || std::isnan(a) || a < 0) return NaN;
  if (x <= 0) return lower ? 0.0 : 1.0;
  if (std::isinf(x)) return lower ? 1.0 : 0.0;
  if (a == 0) return lower ? 1.0 : 0.0;
  if (a >= TEMME_MIN_A) {
    const double m = (x - a) / a;
    const double eta = temme_eta(m);
    if (std::fabs(eta) <= TEMME_MAX_ETA) return temme(a, m, eta, lower);
  }
  if (x < a + 1) {
    if (!lower && a < 0.5) return small_shape_upper(x, a);
    const double pre = dpois_raw(a, x);
    if (pre == 0) return lower ? 0.0 : 1.0;
    double sum = 1, term = 1, ap = a;
    for (int n = 1; n < 1000000; ++n) {
      ap += 1;
      term *= x / ap;
      sum += term;
      if (term < sum * 1e-17) break;
    }
    const double pl = pre * sum;
    return lower ? clamp01(pl) : clamp01(1.0 - pl);  // a NaN stays NaN (std::fmin(x, 1.0) gives 1)
  }
  const double pre = a * dpois_raw(a, x);  // x^a e^{−x} / Γ(a)
  if (pre == 0) return lower ? 1.0 : 0.0;
  double b = x + 1 - a;
  double c = 1.0 / TINY;
  double d = 1.0 / b;
  double h = d;
  for (int i = 1; i < 1000000; ++i) {
    const double an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (std::fabs(d) < TINY) d = TINY;
    c = b + an / c;
    if (std::fabs(c) < TINY) c = TINY;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < 1e-16) break;
  }
  const double qu = pre * h;
  return lower ? clamp01(1.0 - qu) : clamp01(qu);
}

double dgamma_std(double x, double a) {
  if (x < 0) return 0;
  if (x == 0) return a < 1 ? INF : (a == 1 ? 1.0 : 0.0);
  return (a / x) * dpois_raw(a, x);
}

double qgamma(double p, double a, bool lower) {
  if (std::isnan(p) || std::isnan(a) || p < 0 || p > 1 || a <= 0) return NaN;
  if (lower ? p == 0 : p == 1) return 0;
  if (lower ? p == 1 : p == 0) return INF;
  // Solve in the tail whose target probability is small, on a log scale:
  //  * lower tail: g(u) = ln P(a, e^u) − ln p is concave in u = ln x (ln X has a
  //    log-concave density), so Newton converges from any start;
  //  * upper tail: h(x) = ln Q(a, x) − ln q is nearly linear for large x.
  const bool useLower = lower ? (p <= 0.5) : (p > 0.5);
  const double tgt = useLower ? (lower ? p : 1 - p) : (lower ? 1 - p : p);
  const double ltgt = std::log(tgt);
  const double z = qnorm(tgt, useLower);  // normal deviate of the target
  const double c = 1.0 / (9.0 * a);
  const double w = 1.0 - c + z * std::sqrt(c);
  const double xWH = (w > 0) ? a * w * w * w : NaN;  // Wilson–Hilferty

  if (useLower) {
    // P(a,x) ≤ x^a/Γ(a+1) ⇒ the root is ≥ x_small (a rigorous lower bound)
    const double lxs = (ltgt + lgammafn(a + 1)) / a;
    if (lxs < -744) return 0;  // underflows
    double u = lxs;
    if (xWH > 0 && std::log(xWH) > u) u = std::log(xWH);
    double uLo = lxs, uHi = INF;
    for (int it = 0; it < 200; ++it) {
      const double x = std::exp(u);
      const double P = pgamma_raw(x, a, true);
      if (!(P > 0)) {  // far below: move right
        uLo = u;
        u = std::isfinite(uHi) ? 0.5 * (uLo + uHi) : u + 1;
        continue;
      }
      const double g = std::log(P) - ltgt;
      if (g == 0) return x;
      if (g > 0) uHi = std::fmin(uHi, u); else uLo = std::fmax(uLo, u);
      const double gp = x * dgamma_std(x, a) / P;
      double uNew = (gp > 0 && std::isfinite(gp)) ? u - g / gp : NaN;
      if (!(uNew >= uLo && uNew <= uHi)) uNew = std::isfinite(uHi) ? 0.5 * (uLo + uHi) : u + 1;
      if (std::fabs(uNew - u) <= 1e-15 * std::fmax(1.0, std::fabs(u))) return std::exp(uNew);
      u = uNew;
    }
    return std::exp(u);
  }

  // upper tail: initial guess from the asymptotic Q(a,x) ≈ x^{a−1}e^{−x}/Γ(a)
  double x = xWH;
  if (tgt < 1e-3 || !(x > 0)) {
    const double lg = lgammafn(a);
    double xa = std::fmax(a, -ltgt);
    for (int k = 0; k < 6; ++k) xa = std::fmax(-ltgt - lg + (a - 1) * std::log(xa), 1e-3);
    if (!(x > 0) || xa > a + 1) x = xa;
  }
  if (!(x > 0) || !std::isfinite(x)) x = std::fmax(a, 1.0);
  double lo = 0, hi = INF;
  for (int it = 0; it < 200; ++it) {
    const double Q = pgamma_raw(x, a, false);
    if (!(Q > 0)) {  // far above: move left
      hi = x;
      x = 0.5 * (lo + hi);
      continue;
    }
    const double h = std::log(Q) - ltgt;  // decreasing in x
    if (h == 0) return x;
    if (h > 0) lo = std::fmax(lo, x); else hi = std::fmin(hi, x);
    const double hp = -dgamma_std(x, a) / Q;
    double xNew = (hp < 0 && std::isfinite(hp)) ? x - h / hp : NaN;
    if (!(xNew > lo && xNew < hi)) xNew = std::isfinite(hi) ? 0.5 * (lo + hi) : 2 * x + 1;
    if (std::fabs(xNew - x) <= 1e-15 * xNew) return xNew;
    x = xNew;
  }
  return x;
}

double pchisq(double x, double df, bool lower) { return pgamma_raw(x / 2, df / 2, lower); }
double qchisq(double p, double df, bool lower) { return 2 * qgamma(p, df / 2, lower); }
double dchisq(double x, double df) { return 0.5 * dgamma_std(x / 2, df / 2); }

double dpois(double k, double lambda) {
  if (k < 0 || k != std::floor(k)) return 0;
  return dpois_raw(k, lambda);
}

double ppois(double k, double lambda, bool lower) {
  if (k < 0) return lower ? 0.0 : 1.0;
  if (lambda == 0) return lower ? 1.0 : 0.0;
  const double kk = std::floor(k + 1e-7);
  return pgamma_raw(lambda, kk + 1, !lower);  // P(X ≤ k) = Q(k + 1, λ)
}

double qpois(double p, double lambda, bool lower) {
  if (std::isnan(p) || p < 0 || p > 1 || lambda < 0) return NaN;
  if (lambda == 0) return 0;
  const double pl = (lower ? p : 1 - p) * (1 - 64 * DBL_EPSILON_);
  if (pl <= 0) return 0;
  if (pl >= 1) return INF;
  // normal start then local search
  double k = std::floor(lambda + std::sqrt(lambda) * qnorm(pl) + 0.5);
  if (k < 0) k = 0;
  if (lambda >= 1e15) return k;  // a step of one no longer changes k in floating point
  while (k > 0 && ppois(k - 1, lambda) >= pl) k -= 1;
  while (ppois(k, lambda) < pl) k += 1;
  return k;
}

}  // namespace ssc::stats
