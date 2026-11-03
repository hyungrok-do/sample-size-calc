"""Emit C++17 from the numeric expressions lowered by the Lean metaprogram.

These evaluations use binary64 arithmetic. The Lean theorem certifies the
symbolic equivalence, not this numerical evaluation or the C++ toolchain.
"""

PRELUDE = r'''// Generated from elaborated Lean expressions. Requires C++17.
// Floating-point evaluation is approximate, not a Lean proof certificate.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
// Preserve Lean Real's totalized operations, including exceptional arguments.
double lean_div(double a, double b) { return b == 0 ? 0 : a / b; }
double lean_log(double a) { return a == 0 ? 0 : std::log(std::abs(a)); }
double lean_sqrt(double a) { return a <= 0 ? 0 : std::sqrt(a); }
'''


def emit_cpp(lowering):
    kind = lowering['kind']
    if kind == 'analytic_roots':
        values = ', '.join(lowering['roots'])
        return PRELUDE + '\nint main() {\n  std::vector<double> roots = {' + values + r'''};
  for (double root : roots) if (!std::isfinite(root)) {
    std::cerr << "Non-finite numerical evaluation\n"; return 2;
  }
  std::sort(roots.begin(), roots.end());
  roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
  std::cout << std::setprecision(17) << "{\"method\":\"analytic_roots\",\"kernelCertified\":false,\"roots\":[";
  for (size_t i = 0; i < roots.size(); ++i) std::cout << (i ? "," : "") << roots[i];
  std::cout << "]}\n";
}
'''
    if kind == 'bracketed_bisection':
        return PRELUDE + '\ndouble residual(double x) { return ' + lowering['expression'] + r'''; }
double bound(const char* text) {
  char* end = nullptr;
  double value = std::strtod(text, &end);
  if (end == text || *end || !std::isfinite(value)) throw std::runtime_error("Invalid finite bound");
  return value;
}
int main(int argc, char** argv) {
  try {
    if (argc != 3) throw std::runtime_error("Usage: equation LOWER UPPER (a continuous, sign-changing bracket)");
    double lo = bound(argv[1]), hi = bound(argv[2]);
    if (!(lo < hi)) throw std::runtime_error("LOWER must be less than UPPER");
    double fl = residual(lo), fh = residual(hi), x = lo, fx = fl;
    if (!std::isfinite(fl) || !std::isfinite(fh)) throw std::runtime_error("Non-finite endpoint residual");
    if (fl == 0) { x = lo; fx = fl; }
    else if (fh == 0) { x = hi; fx = fh; }
    else {
      if (std::signbit(fl) == std::signbit(fh)) throw std::runtime_error("No sign-changing bracket; this does not prove no solution exists");
      for (int iteration = 0; iteration < 200; ++iteration) {
        x = lo / 2 + hi / 2;
        fx = residual(x);
        if (!std::isfinite(fx)) throw std::runtime_error("Non-finite residual inside bracket");
        if (std::abs(fx) <= 1e-10) break;
        if (x == lo || x == hi) throw std::runtime_error("Precision exhausted before residual tolerance");
        if (std::signbit(fl) == std::signbit(fx)) { lo = x; fl = fx; }
        else hi = x;
      }
    }
    if (std::abs(fx) > 1e-10) throw std::runtime_error("Residual tolerance not reached");
    std::cout << std::setprecision(17) << "{\"method\":\"bracketed_bisection\",\"solution\":" << x
              << ",\"residual\":" << fx << ",\"kernelCertified\":false}\n";
  } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
}
'''
    return ''
