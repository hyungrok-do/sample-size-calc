// Hypergeometric distribution, bivariate normal probabilities and
// Gauss–Legendre quadrature rules.
#include <map>
#include <memory>
#include <mutex>

#include "stats.hpp"

namespace ssc::stats {

// ------------------------------------------------------------ hypergeometric

double dhyper(double x, double m, double n, double k) {
  // X = number of white balls among k drawn from m white + n black
  if (x < 0 || x != std::floor(x)) return 0;
  if (x < k - n || x > k || x > m) return 0;
  if (k == 0) return x == 0 ? 1.0 : 0.0;
  const double N = m + n;
  const double p = k / N, q = (N - k) / N;
  // Loader's form: b(x; m, p) b(k − x; n, p) / b(k; N, p)
  const double p1 = dbinom_raw(x, m, p, q);
  const double p2 = dbinom_raw(k - x, n, p, q);
  const double p3 = dbinom_raw(k, N, p, q);
  return p1 * p2 / p3;
}

double phyper(double x, double m, double n, double k, bool lower) {
  x = std::floor(x + 1e-7);
  const double lo = std::fmax(0.0, k - n), hi = std::fmin(k, m);
  if (x < lo) return lower ? 0.0 : 1.0;
  if (x >= hi) return lower ? 1.0 : 0.0;
  // sum the smaller tail, walking the pmf with the ratio recurrence
  const double mean = k * m / (m + n);
  bool sumLower = x < mean;
  double total = 0;
  if (sumLower) {
    double d = dhyper(x, m, n, k);
    double t = x;
    while (t >= lo && d > 0) {
      total += d;
      // f(t−1)/f(t) = t (n − k + t) / ((m − t + 1)(k − t + 1))
      d *= t * (n - k + t) / ((m - t + 1) * (k - t + 1));
      t -= 1;
      if (d < total * 1e-17) break;
    }
    return lower ? clamp01(total) : clamp01(1.0 - total);  // a NaN stays NaN (std::fmin(x, 1.0) gives 1)
  }
  double t = x + 1;
  double d = dhyper(t, m, n, k);
  while (t <= hi && d > 0) {
    total += d;
    // f(t+1)/f(t) = (m − t)(k − t) / ((t + 1)(n − k + t + 1))
    d *= (m - t) * (k - t) / ((t + 1) * (n - k + t + 1));
    t += 1;
    if (d < total * 1e-17) break;
  }
  return lower ? clamp01(1.0 - total) : clamp01(total);
}

// ------------------------------------------------------------ Gauss–Legendre

const GLRule& gauss_legendre(int n) {
  static std::mutex mtx;
  static std::map<int, std::unique_ptr<GLRule>> cache;
  std::lock_guard<std::mutex> lock(mtx);
  auto it = cache.find(n);
  if (it != cache.end()) return *it->second;
  auto rule = std::make_unique<GLRule>();
  rule->x.assign(n, 0.0);
  rule->w.assign(n, 0.0);
  const int m = (n + 1) / 2;
  for (int i = 0; i < m; ++i) {
    double z = std::cos(PI * (i + 0.75) / (n + 0.5));
    double pp = 0;
    for (int it2 = 0; it2 < 100; ++it2) {
      double p1 = 1, p2 = 0;
      for (int j = 0; j < n; ++j) {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j + 1) * z * p2 - j * p3) / (j + 1);
      }
      pp = n * (z * p1 - p2) / (z * z - 1);
      const double z1 = z;
      z = z1 - p1 / pp;
      if (std::fabs(z - z1) < 1e-15) {
        // final derivative at converged node
        p1 = 1;
        p2 = 0;
        for (int j = 0; j < n; ++j) {
          const double p3 = p2;
          p2 = p1;
          p1 = ((2.0 * j + 1) * z * p2 - j * p3) / (j + 1);
        }
        pp = n * (z * p1 - p2) / (z * z - 1);
        break;
      }
    }
    rule->x[i] = -z;
    rule->x[n - 1 - i] = z;
    rule->w[i] = rule->w[n - 1 - i] = 2.0 / ((1 - z * z) * pp * pp);
  }
  const GLRule& ref = *rule;
  cache.emplace(n, std::move(rule));
  return ref;
}

// ------------------------------------------------------------ bivariate normal
// A. Genz (2004), "Numerical computation of rectangular bivariate and
// trivariate normal and t probabilities", Statistics and Computing 14:251–260,
// building on Drezner & Wesolowsky (1990).

double bvn_upper(double dh, double dk, double r) {
  if (std::isnan(dh) || std::isnan(dk) || std::isnan(r)) return NAN;
  if (dh == INFINITY || dk == INFINITY) return 0;
  if (dh == -INFINITY) return pnorm(-dk);  // the formulas below meet 0·∞ or ∞ − ∞ at an infinite limit
  if (dk == -INFINITY) return pnorm(-dh);
  int ng;
  if (std::fabs(r) < 0.3)
    ng = 6;
  else if (std::fabs(r) < 0.75)
    ng = 12;
  else
    ng = 20;
  const GLRule& gl = gauss_legendre(ng);
  const int lg = ng / 2;  // use the negative half of the nodes
  double h = dh, k = dk, hk = h * k, bvn = 0;
  if (std::fabs(r) < 0.925) {
    const double hs = (h * h + k * k) / 2;
    const double asr = std::asin(r);
    for (int i = 0; i < lg; ++i) {
      const double xi = gl.x[i], wi = gl.w[i];
      double sn = std::sin(asr * (xi + 1) / 2);
      bvn += wi * std::exp((sn * hk - hs) / (1 - sn * sn));
      sn = std::sin(asr * (-xi + 1) / 2);
      bvn += wi * std::exp((sn * hk - hs) / (1 - sn * sn));
    }
    bvn = bvn * asr / (4 * PI) + pnorm(-h) * pnorm(-k);
  } else {
    if (r < 0) {
      k = -k;
      hk = -hk;
    }
    if (std::fabs(r) < 1) {
      const double as = (1 - r) * (1 + r);
      double a = std::sqrt(as);
      const double bs = (h - k) * (h - k);
      const double c = (4 - hk) / 8;
      const double d = (12 - hk) / 16;
      double asr = -(bs / as + hk) / 2;
      if (asr > -100) bvn = a * std::exp(asr) * (1 - c * (bs - as) * (1 - d * bs / 5) / 3 + c * d * as * as / 5);
      if (hk > -100) {
        const double b = std::sqrt(bs);
        const double sp = SQRT_2PI * pnorm(-b / a);
        bvn -= std::exp(-hk / 2) * sp * b * (1 - c * bs * (1 - d * bs / 5) / 3);
      }
      a /= 2;
      for (int i = 0; i < lg; ++i) {
        for (int sgn = -1; sgn <= 1; sgn += 2) {
          const double xs0 = a * (sgn * gl.x[i] + 1);
          const double xs = xs0 * xs0;
          const double rs = std::sqrt(1 - xs);
          asr = -(bs / xs + hk) / 2;
          if (asr > -100) {
            const double sp = 1 + c * xs * (1 + d * xs);
            const double ep = std::exp(-hk * (1 - rs) / (2 * (1 + rs))) / rs;
            bvn += a * gl.w[i] * std::exp(asr) * (ep - sp);
          }
        }
      }
      bvn = -bvn / (2 * PI);
    }
    if (r > 0) {
      bvn += pnorm(-std::fmax(h, k));
    } else {
      bvn = -bvn;
      if (k > h) bvn += pnorm(k) - pnorm(h);
    }
  }
  return clamp01(bvn);
}

double pbvnorm(double h, double k, double r) { return bvn_upper(-h, -k, r); }

}  // namespace ssc::stats
