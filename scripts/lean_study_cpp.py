"""Emit statistical C++ from a numerical planning specification elaborated in Lean.

The model produces an effect and per-independent-block variance, not a canned
sample size. All designs use the same two-sided Gaussian/Wald power inversion.
"""
from fractions import Fraction
import json
import math


def number(value):
    value = float(Fraction(value))
    if not math.isfinite(value):
        raise ValueError('Planning value is outside the finite numerical range.')
    return format(value, '.17g')


def emit_study_cpp(plan):
    stats, criterion = plan['statistics'], plan['criterion']
    alpha, target = float(Fraction(criterion['alpha'])), float(Fraction(criterion['targetPower']))
    if not 1e-12 <= alpha <= 1-1e-12 or not alpha < target <= 1-1e-12:
        raise ValueError('Numerical power evaluation requires 1e-12 <= alpha < targetPower <= 1-1e-12.')
    if float(Fraction(stats['variance'])) <= 0:
        raise ValueError('The positive planning variance is below numerical precision.')
    mode = plan['solveFor']
    code = r'''// Statistical design compiled from a Lean study. Standard-library C++17.
// The exact model arithmetic is Lean checked. Floating-point power is approximate.
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
double upper_tail(double z) { return 0.5*std::erfc(z/std::sqrt(2.0)); }
double critical_value(double a) {
  double lo=0,hi=40;
  for(int i=0;i<120;++i) {
    double mid=(lo+hi)/2;
    if(upper_tail(mid)>a/2) lo=mid; else hi=mid;
  }
  return (lo+hi)/2;
}
double normal_quantile(double probability) {
  double lo=-40,hi=40;
  for(int i=0;i<120;++i) {
    double mid=(lo+hi)/2;
    if(upper_tail(-mid)<probability) lo=mid; else hi=mid;
  }
  return (lo+hi)/2;
}
'''
    for name, value in [('variance', stats['variance']), ('effect', stats['effect']),
                        ('alpha', criterion['alpha']), ('target', criterion['targetPower'])]:
        code += f'constexpr double {name} = {number(value)};\n'
    code += r'''
double power_lambda(double lambda, double z) {
  double p=upper_tail(z+lambda)+upper_tail(z-lambda);
  if(!std::isfinite(p)) throw std::runtime_error("Non-finite power; rescale planning assumptions");
  return p;
}
double power_at(int64_t n, double z) {
  double lambda=std::abs(effect)/std::sqrt(variance)*std::sqrt(double(n));
  if(!std::isfinite(lambda)) throw std::runtime_error("Non-finite noncentrality; rescale planning assumptions");
  return power_lambda(lambda,z);
}
void counts(int64_t n) {
'''
    multiplier = sum(g['unitsPerN']*g['participantsPerUnit'] for g in plan['groups'])
    code += f'  std::cout << ",\\\"totalParticipants\\\":" << n*INT64_C({multiplier});\n'
    code += '  std::cout << ",\\\"groups\\\":[";\n'
    for i, group in enumerate(plan['groups']):
        prefix = (',' if i else '') + '{"label":' + json.dumps(group['label'], ensure_ascii=False) + ',"count":'
        code += f'  std::cout << {json.dumps(prefix, ensure_ascii=False)} << n*INT64_C({group["unitsPerN"]})'
        code += f' << ",\\\"participants\\\":" << n*INT64_C({group["unitsPerN"]*group["participantsPerUnit"]}) << "}}";\n'
    code += '  std::cout << "]";\n}\n'
    code += r'''
int main(int argc, char**) {
  try {
    if(argc!=1) throw std::runtime_error("Planning values are compiled from Lean source; edit the study and recompile");
    const double z=critical_value(alpha);
    std::cout << std::setprecision(17);
    std::cout << "{\"kernelCertified\":false,\"criticalValue\":" << z
              << ",\"varianceCoefficient\":" << variance << ",\"alpha\":" << alpha
              << ",\"targetPower\":" << target;
'''
    if mode == 'sample_size':
        code += r'''
    if(effect!=0) {
      double approximate=variance*std::pow(z+normal_quantile(target),2)/(effect*effect);
      if(std::isfinite(approximate)) std::cout << ",\"dominantTailApproximation\":" << approximate;
    }
    double previous=0;
'''
        code += f'    for(int64_t n={criterion["minN"]}; n<={criterion["maxN"]}; ++n) {{\n'
        code += r'''
      double p=power_at(n,z);
      if(std::abs(p-target)<=1e-12) {
        std::cout << ",\"status\":\"numerically_unresolved\",\"unresolvedN\":" << n << ",\"power\":" << p << "}\n"; return 0;
      }
      if(p>target) {
        std::cout << ",\"status\":\"found\",\"selectedN\":" << n << ",\"power\":" << p
                  << ",\"standardError\":" << std::sqrt(variance)/std::sqrt(double(n))
                  << ",\"noncentrality\":" << std::abs(effect)/std::sqrt(variance)*std::sqrt(double(n));
'''
        code += f'        std::cout << ",\\\"evaluated\\\":" << n-{criterion["minN"]}+1;\n'
        code += f'        if(n>{criterion["minN"]}) std::cout << ",\\\"predecessorN\\\":" << n-1 << ",\\\"predecessorPower\\\":" << previous;\n'
        code += r'''
        counts(n);
        std::cout << ",\"minimality\":\"all_smaller_admissible_counts_fail_in_binary64\"}\n"; return 0;
      }
      previous=p;
      if(effect==0) break;
    }
    std::cout << ",\"status\":\"not_found_in_domain\",\"powerAtLastEvaluated\":" << previous << "}\n";
'''
    elif mode == 'power':
        n = plan['fixedN']
        code += f'    const int64_t n={n};\n'
        code += r'''
    std::cout << ",\"status\":\"evaluated\",\"n\":" << n << ",\"power\":" << power_at(n,z)
              << ",\"standardError\":" << std::sqrt(variance)/std::sqrt(double(n))
              << ",\"noncentrality\":" << std::abs(effect)/std::sqrt(variance)*std::sqrt(double(n));
    counts(n); std::cout << "}\n";
'''
    else:
        code += f'    const int64_t n={plan["fixedN"]};\n'
        code += r'''
    double lo=0,hi=64;
    for(int i=0;i<120;++i) {
      double mid=(lo+hi)/2;
      if(mid==lo || mid==hi) break;
      if(power_lambda(mid,z)>=target) hi=mid; else lo=mid;
    }
    double scale=std::sqrt(variance)/std::sqrt(double(n));
    if(!std::isfinite(hi*scale)) throw std::runtime_error("Detectable effect overflowed");
    std::cout << ",\"status\":\"approximated\",\"n\":" << n
              << ",\"detectableMagnitude\":" << hi*scale
              << ",\"effectLower\":" << lo*scale << ",\"effectUpper\":" << hi*scale
              << ",\"powerAtLower\":" << power_lambda(lo,z) << ",\"powerAtUpper\":" << power_lambda(hi,z);
    counts(n); std::cout << "}\n";
'''
    return code + r'''  } catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
}
'''
