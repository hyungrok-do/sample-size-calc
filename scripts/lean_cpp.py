"""Package C++17 expressions and guarded solution branches lowered by Lean.

Binary64 values and guard comparisons are approximate. Only the symbolic
Lean equivalence is kernel checked, not this exporter or the C++ execution.
"""
import json

PRELUDE = r'''// Generated from elaborated Lean expressions. Requires C++17.
// Floating-point values and guard comparisons are not Lean certificates.
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
double checked(double value) {
  if (!std::isfinite(value)) throw std::runtime_error("Non-finite numerical evaluation");
  return value;
}
// Preserve Lean Real's totalized operations, including exceptional arguments.
double lean_div(double a, double b) { return b == 0 ? 0 : a / b; }
double lean_log(double a) { return a == 0 ? 0 : std::log(std::abs(a)); }
double lean_sqrt(double a) { return a <= 0 ? 0 : std::sqrt(a); }
double finite_input(const std::string& text) {
  char* end = nullptr;
  errno = 0;
  double value = std::strtod(text.c_str(), &end);
  if (end == text.c_str() || *end || !std::isfinite(value) || errno == ERANGE)
    throw std::runtime_error("Invalid finite bound or parameter: " + text);
  return value;
}
std::vector<double> read_parameters(int argc, char** argv, int start,
                                    const std::vector<std::string>& names) {
  std::vector<double> values(names.size());
  std::vector<bool> supplied(names.size(), false);
  for (int i = start; i < argc; ++i) {
    std::string argument(argv[i]);
    auto equal = argument.find('=');
    if (equal == std::string::npos) throw std::runtime_error("Expected parameter NAME=VALUE: " + argument);
    auto name = argument.substr(0, equal);
    auto found = std::find(names.begin(), names.end(), name);
    if (found == names.end()) throw std::runtime_error("Unknown parameter: " + name);
    auto index = size_t(found - names.begin());
    if (supplied[index]) throw std::runtime_error("Duplicate parameter: " + name);
    values[index] = finite_input(argument.substr(equal + 1));
    supplied[index] = true;
  }
  for (size_t i = 0; i < names.size(); ++i)
    if (!supplied[i]) throw std::runtime_error("Missing parameter: " + names[i]);
  return values;
}
'''


def literal(text):
    """Quote display names as UTF-8 C++ string literals, never identifiers."""
    return json.dumps(text, ensure_ascii=False)


def parameter_setup(lowering, start):
    parameters = lowering.get('parameters', [])
    names = ', '.join(literal(p['name']) for p in parameters)
    lines = [f'    const auto parameters = read_parameters(argc, argv, {start}, {{{names}}});']
    for i, parameter in enumerate(parameters):
        lines.append(f"    const double {parameter['cppName']} = parameters[{i}];")
    for assumption in lowering.get('assumptions', []):
        error = literal('Assumption failed (binary64 evaluation): ' + assumption['statement'])
        lines.append(f"    if (!({assumption['expression']})) throw std::runtime_error({error});")
    return '\n'.join(lines) + '\n'


CATCH = r'''  } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
}
'''


def emit_cpp(lowering):
    kind = lowering['kind']
    if kind == 'analytic_roots':
        branches = lowering.get('branches')
        if branches is None:  # Reports from the original closed-root exporter.
            branches = [{'condition': 'true', 'root': value} for value in lowering['roots']]
        body = ['    std::vector<double> roots;', '    bool all_reals = false;']
        # Evaluate all-real branches first: no irrelevant finite roots are needed.
        for branch in branches:
            if branch['root'] is None:
                body.append(f"    if ({branch['condition']}) all_reals = true;")
        body.append('    if (!all_reals) {')
        for branch in branches:
            if branch['root'] is not None:
                body.append(f"      if ({branch['condition']}) roots.push_back(checked({branch['root']}));")
        body.append('    }')
        return PRELUDE + '\nint main(int argc, char** argv) {\n  try {\n' + parameter_setup(lowering, 1) + '\n'.join(body) + r'''
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
    const char* set = all_reals ? "all_reals" : roots.empty() ? "empty" : "finite";
    std::cout << std::setprecision(17) << "{\"method\":\"analytic_roots\",\"kernelCertified\":false,\"solutionSet\":\""
              << set << "\",\"roots\":[";
    for (size_t i = 0; i < roots.size(); ++i) std::cout << (i ? "," : "") << roots[i];
    std::cout << "]}\n";
''' + CATCH
    if kind == 'bracketed_bisection':
        return PRELUDE + r'''
int main(int argc, char** argv) {
  try {
    if (argc < 3) throw std::runtime_error("Usage: equation LOWER UPPER [NAME=VALUE ...] (a continuous, sign-changing bracket)");
    double lo = finite_input(argv[1]), hi = finite_input(argv[2]);
    if (!(lo < hi)) throw std::runtime_error("LOWER must be less than UPPER");
''' + parameter_setup(lowering, 3) + '\n    const auto residual = [=](double x) { return ' + lowering['expression'] + r'''; };
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
''' + CATCH
    return ''
