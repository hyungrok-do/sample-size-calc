// Random numbers for the simulation procedures.
//
// Every simulated quantity comes from its own stream, identified by the seed,
// the replicate and a unit (cluster or subject) index, so that
//  * results do not depend on the number of threads or the order of work, and
//  * a larger study contains the units of a smaller one: with the same seed,
//    the first k clusters of a replicate are identical whatever the total
//    number, which makes simulated power smooth in the sample size (common
//    random numbers) and the sample-size search stable.
//
// The generator is xoshiro256++ (Blackman & Vigna 2021), seeded through
// SplitMix64 (Steele, Lea & Flood 2014); the variates use exact methods
// (Marsaglia's polar method, Marsaglia–Tsang gamma, inversion or Hörmann's BTRS for binomial).
#pragma once

#include <cmath>
#include <cstdint>

namespace ssc::sim {

/// SplitMix64 step: advances `x` and returns a well-mixed 64-bit value.
inline std::uint64_t splitmix64(std::uint64_t& x) {
  std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

/// Seed of the stream (seed, a, b, c): distinct arguments give unrelated streams.
inline std::uint64_t stream_seed(std::uint64_t seed, std::uint64_t a, std::uint64_t b = 0, std::uint64_t c = 0) {
  std::uint64_t x = seed;
  std::uint64_t h = splitmix64(x);
  for (std::uint64_t v : {a, b, c}) {
    x = h ^ (v * 0xd1342543de82ef95ULL + 0x632be59bd9b4e019ULL);
    h = splitmix64(x);
  }
  return h;
}

class Rng {
 public:
  explicit Rng(std::uint64_t seed) {
    std::uint64_t x = seed;
    for (auto& w : s_) w = splitmix64(x);
  }
  Rng(std::uint64_t seed, std::uint64_t a, std::uint64_t b = 0, std::uint64_t c = 0)
      : Rng(stream_seed(seed, a, b, c)) {}

  std::uint64_t next() {
    const std::uint64_t result = rotl(s_[0] + s_[3], 23) + s_[0];
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }

  /// Uniform on the open interval (0, 1), with 53 random bits.
  double uniform() { return (static_cast<double>(next() >> 11) + 0.5) * 0x1.0p-53; }

  /// Standard normal (Marsaglia's polar method; the second value is kept).
  double normal() {
    if (has_spare_) {
      has_spare_ = false;
      return spare_;
    }
    double u, v, s;
    do {
      u = 2 * uniform() - 1;
      v = 2 * uniform() - 1;
      s = u * u + v * v;
    } while (s >= 1 || s == 0);
    const double f = std::sqrt(-2 * std::log(s) / s);
    spare_ = v * f;
    has_spare_ = true;
    return u * f;
  }

  double exponential() { return -std::log(uniform()); }

  /// Gamma(shape, 1) by Marsaglia & Tsang (2000); shape < 1 uses the boost U^(1/shape).
  double gamma(double shape) {
    if (shape < 1) return gamma(shape + 1) * std::pow(uniform(), 1 / shape);
    const double d = shape - 1.0 / 3, c = 1 / std::sqrt(9 * d);
    while (true) {
      double x, v;
      do {
        x = normal();
        v = 1 + c * x;
      } while (v <= 0);
      v = v * v * v;
      const double u = uniform();
      if (u < 1 - 0.0331 * x * x * x * x || std::log(u) < 0.5 * x * x + d * (1 - v + std::log(v))) return d * v;
    }
  }

  double chisq(double df) { return 2 * gamma(df / 2); }

  /// Binomial(n, p): inversion when n·min(p, 1 − p) < 10, otherwise Hörmann's transformed rejection
  /// with squeeze (BTRS), whose time does not grow with n.
  long binomial(long n, double p) {
    if (n <= 0 || p <= 0) return 0;
    if (p >= 1) return n;
    const bool flip = p > 0.5;
    const double pp = flip ? 1 - p : p;
    const long k = static_cast<double>(n) * pp < 10 ? binomial_inversion(n, pp) : binomial_btrs(n, pp);
    return flip ? n - k : k;
  }

 private:
  /// Inversion, in chunks small enough that q^n cannot underflow (p ≤ 1/2).
  long binomial_inversion(long n, double p) {
    const double ratio = p / (1 - p);
    long total = 0;
    for (long left = n; left > 0;) {
      const long m = left < 500 ? left : 500;  // (1/2)^500 is still a normal double
      left -= m;
      double f = std::pow(1 - p, static_cast<double>(m)), u = uniform();
      long k = 0;
      while (u > f && k < m) {
        u -= f;
        f *= ratio * static_cast<double>(m - k) / static_cast<double>(k + 1);
        ++k;
      }
      total += k;
    }
    return total;
  }

  /// log k! − ((k + ½) log(k + 1) − (k + 1) + ½ log 2π), the error of Stirling's formula.
  static double stirling_tail(double k) {
    static const double small[10] = {0.0810614667953272,  0.0413406959554092,  0.0276779256849983,
                                     0.02079067210376509, 0.0166446911898211,  0.0138761288230707,
                                     0.0118967099458917,  0.0104112652619720,  0.00925546218271273,
                                     0.00833056343336287};
    if (k <= 9) return small[static_cast<int>(k)];
    const double k1 = k + 1, k1sq = k1 * k1;
    return (1.0 / 12 - (1.0 / 360 - 1.0 / 1260 / k1sq) / k1sq) / k1;
  }

  /// BTRS (Hörmann 1993) for p ≤ 1/2 and np ≥ 10.
  long binomial_btrs(long n, double p) {
    const double nd = static_cast<double>(n), spq = std::sqrt(nd * p * (1 - p));
    const double b = 1.15 + 2.53 * spq, a = -0.0873 + 0.0248 * b + 0.01 * p, c = nd * p + 0.5;
    const double vr = 0.92 - 4.2 / b, r = p / (1 - p), alpha = (2.83 + 5.1 / b) * spq;
    const double m = std::floor((nd + 1) * p);
    const double hm = (m + 0.5) * std::log((m + 1) / (r * (nd - m + 1))) + stirling_tail(m) + stirling_tail(nd - m);
    while (true) {
      const double u = uniform() - 0.5, us = 0.5 - std::fabs(u);
      double v = uniform();
      const double k = std::floor((2 * a / us + b) * u + c);
      if (us >= 0.07 && v <= vr) return static_cast<long>(k);  // inside the box: always accepted
      if (k < 0 || k > nd) continue;
      v = std::log(v * alpha / (a / (us * us) + b));
      const double bound = hm + (nd + 1) * std::log((nd - m + 1) / (nd - k + 1)) +
                           (k + 0.5) * std::log(r * (nd - k + 1) / (k + 1)) - stirling_tail(k) -
                           stirling_tail(nd - k);
      if (v <= bound) return static_cast<long>(k);
    }
  }

  static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  std::uint64_t s_[4];
  double spare_ = 0;
  bool has_spare_ = false;
};

}  // namespace ssc::sim
