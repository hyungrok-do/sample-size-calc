// Exact, deliberately bounded rational arithmetic for derivation/proof
// emission. This is separate from binary64 execution. Unsupported expressions
// stay numerical.
#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>

namespace ssc::design::algebra {
struct Unsupported {};
struct Rational {
  std::int64_t n = 0, d = 1;
  Rational(std::int64_t numerator = 0, std::int64_t denominator = 1) {
    if (!denominator)
      throw Unsupported{};
    if (denominator < 0) {
      numerator = -numerator;
      denominator = -denominator;
    }
    auto g = std::gcd(numerator, denominator);
    n = numerator / g;
    d = denominator / g;
    // Products and cross-products of two admitted values fit int64_t.
    if (std::abs(n) > 1000000000 || d > 1000000000)
      throw Unsupported{};
  }
  std::string text() const {
    return d == 1 ? std::to_string(n)
                  : "(" + std::to_string(n) + " / " + std::to_string(d) + ")";
  }
  bool zero() const { return n == 0; }
};
inline Rational operator+(Rational a, Rational b) {
  return {a.n * b.d + b.n * a.d, a.d * b.d};
}
inline Rational operator-(Rational a, Rational b) {
  return {a.n * b.d - b.n * a.d, a.d * b.d};
}
inline Rational operator*(Rational a, Rational b) {
  return {a.n * b.n, a.d * b.d};
}
inline Rational operator/(Rational a, Rational b) {
  return {a.n * b.d, a.d * b.n};
}
inline bool operator<(Rational a, Rational b) { return a.n * b.d < b.n * a.d; }
inline Rational decimal(std::string s) {
  if (s.size() > 30)
    throw Unsupported{};
  int exponent = 0;
  auto e = s.find_first_of("eE");
  if (e != std::string::npos) {
    try {
      exponent = std::stoi(s.substr(e + 1));
    } catch (...) {
      throw Unsupported{};
    }
    s.resize(e);
  }
  if (exponent < -30 || exponent > 30)
    throw Unsupported{};
  auto dot = s.find('.');
  if (dot != std::string::npos) {
    exponent -= int(s.size() - dot - 1);
    s.erase(dot, 1);
  }
  if (s.empty() || s.size() > 9 || exponent < -9 || exponent > 9)
    throw Unsupported{};
  std::int64_t n = 0;
  for (char c : s) {
    if (c < '0' || c > '9')
      throw Unsupported{};
    n = n * 10 + c - '0';
  }
  Rational r(n);
  for (int i = 0; i < exponent; ++i)
    r = r * Rational(10);
  for (int i = 0; i > exponent; --i)
    r = r / Rational(10);
  return r;
}
using Polynomial = std::array<Rational, 3>;
struct Expression {
  Polynomial p{};
  std::string lean;
};
inline Expression constant(Rational r) { return {{r, {}, {}}, r.text()}; }
inline Expression binary(const std::string &op, const Expression &a,
                         const Expression &b) {
  Expression out;
  out.lean = "(" + a.lean + " " + op + " " + b.lean + ")";
  if (op == "+" || op == "-")
    for (int i = 0; i < 3; ++i)
      out.p[i] = op == "+" ? a.p[i] + b.p[i] : a.p[i] - b.p[i];
  else if (op == "*")
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        if (i + j > 2) {
          if (!a.p[i].zero() && !b.p[j].zero())
            throw Unsupported{};
        } else
          out.p[i + j] = out.p[i + j] + a.p[i] * b.p[j];
      }
  else if (op == "/") {
    if (!b.p[1].zero() || !b.p[2].zero() || b.p[0].zero())
      throw Unsupported{};
    for (int i = 0; i < 3; ++i)
      out.p[i] = a.p[i] / b.p[0];
  } else
    throw Unsupported{};
  return out;
}
inline Rational squareRoot(Rational r) {
  if (r.n < 0)
    throw Unsupported{};
  auto n = std::int64_t(std::sqrt(double(r.n))),
       d = std::int64_t(std::sqrt(double(r.d)));
  if (n * n != r.n || d * d != r.d)
    throw Unsupported{};
  return {n, d};
}
} // namespace ssc::design::algebra
